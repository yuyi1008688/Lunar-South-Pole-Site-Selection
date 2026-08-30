# -*- coding: utf-8 -*-
"""
L3_verify_energy.py —— L3 端到端验证：新 path_profile.csv 驱动（算法未改动的）
三态能量仿真，对比历史基准
================================================================
基准（历史值，竞赛期实测留痕）：全程最低 SoC = 215.2 Wh（必须 >200 热控底线）、
驻留充电 8 次、任务周期 8 个月昼（历史逐点 energy_curve.csv 未存档，故按
关键统计量对比）。

口径说明（重要，根因见报告）：
  - 本仓库交付剖面按提示词约定采用**双线性插值**取值（消除最近邻的台阶噪声）；
  - 归因实验表明"最近邻台阶剖面"是历史驻留/月昼次数的主要成因：
    最近邻变体 7 次/6 月昼（接近基准 8/8），双线性 3 次/5 月昼；
  - 无论哪种口径，三项工程验收判定（最低 SoC>200 / 驻留<15 / 月昼<10）
    全部保持——即"换了寻路实现与取值方式，下游工程结论不翻盘"。

本脚本做两件事：
  1. 对交付剖面（data/output/stage2/path_profile.csv，双线性）跑未改动的
     energy_simulation.py，与基准逐项对比并给 PASS/FAIL；
  2. 附最近邻敏感性变体（仅本验证用），量化取值方式的影响。

输出：results/decouple_verification/L3_energy_comparison.png + tests/L3_metrics.json
用法：python tests/L3_verify_energy.py
"""

import csv
import importlib.util
import json
import os
import shutil
import sys

import numpy as np

sys.stdout.reconfigure(encoding='utf-8')
_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_HERE, 'src'))
sys.path.insert(0, os.path.join(_HERE, 'src', 'utils'))
sys.path.insert(0, os.path.join(_HERE, 'src', 'ch06_path_planning'))

ESIM = os.path.join(_HERE, 'src', 'ch06_path_planning', 'energy_simulation.py')
FIG_DIR = os.path.join(_HERE, 'results', 'decouple_verification')

BASE = dict(min_soc=215.2, dwell=8, lunar=8)
TOL = dict(min_soc=5.0, dwell=1, lunar=1)   # 提示词 §6 预设容差


def run_esim(profile_csv, out_tag):
    """以独立命名空间运行未改动的 energy_simulation（仅替换输入输出路径）。"""
    src = open(ESIM, encoding='utf-8').read()
    src = src.replace('PROFILE_CSV = os.path.join(_base, "stage2", "path_profile.csv")',
                      f'PROFILE_CSV = r"{profile_csv}"')
    src = src.replace('OUT_DIR     = os.path.join(_base, "stage3_energy")',
                      f'OUT_DIR     = r"{os.path.join(_HERE, "data", "output", "_esim_" + out_tag)}"')
    ns = {'__name__': '__esim__'}
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(compile(src, ESIM, 'exec'), ns)
        ns['main']()
    out_dir = os.path.join(_HERE, 'data', 'output', '_esim_' + out_tag)
    curve = os.path.join(out_dir, 'energy_curve.csv')
    rows = []
    with open(curve, encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            rows.append(dict(s=float(r['s_m']), soc=float(r['soc_wh']),
                             state=r['state']))
    printed = buf.getvalue()
    stats = {}
    for line in printed.splitlines():
        # 只取首次出现（主指标行在验收行之前，验收行同样含关键词）
        if '最低 SoC' in line and 'min_soc' not in stats:
            stats['min_soc'] = float(line.split('Wh')[0].split()[-1])
        elif '驻留充电' in line and 'dwell' not in stats:
            stats['dwell'] = int(line.split('次')[0].split()[-1])
        elif '任务周期' in line and 'lunar' not in stats:
            stats['lunar'] = int(line.split('个月昼')[0].split()[-1].split()[-1])
    return stats, rows


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    metrics = {}
    failures = []
    print('── L3 端到端能量仿真验证（基准：SoC 215.2 Wh / 驻留 8 / 月昼 8） ──\n')

    variants = []
    # 1) 交付剖面（双线性）
    delivered = os.path.join(_HERE, 'data', 'output', 'stage2', 'path_profile.csv')
    if not os.path.exists(delivered):
        raise SystemExit('未找到交付剖面，请先运行 src/ch06_path_planning/stage2_purepython.py')
    stats, rows = run_esim(delivered, 'delivered')
    variants.append(('交付剖面（双线性）', stats, rows))

    # 2) 最近邻敏感性变体（仅验证用）：用最近邻重建剖面后重跑
    from raster_grid import RasterGrid
    from lcp_dijkstra import dijkstra_accumulation, trace_backlink
    from bspline_smooth import bspline_smooth, resample

    def rc_at(x, y):
        return int(round((46080.0 - y) / 240.0 - 0.5)), int(round((x + 46080.0) / 240.0 - 0.5))

    cost = RasterGrid.from_file(os.path.join(_HERE, 'data', 'rasters', 'cost_surface.tif')).arr
    G = {n: RasterGrid.from_file(os.path.join(_HERE, 'data', 'rasters', n + '.tif'))
         for n in ['DEM_fused', 'slope_deg', 'vrm_5x5', 'AVGVISIB_probability']}
    rS, cS = rc_at(44760.0, 10920.0)
    acc, bp = dijkstra_accumulation(cost + 0.001, (rS, cS))
    rE, cE = rc_at(40920.0, 10920.0)
    trace = trace_backlink(bp, (rE, cE), (rS, cS))[::-1]
    raw = [G['AVGVISIB_probability'].xy_at_rc(r, c) for r, c in trace]
    sm5 = bspline_smooth(raw, degree=5)
    nn_profile = os.path.join(_HERE, 'data', 'output', '_profile_nearest.csv')
    pts, _ = resample(sm5, 50.0)
    with open(nn_profile, 'w', newline='', encoding='utf-8-sig') as fo:
        w = csv.writer(fo)
        w.writerow(["id", "s_m", "x", "y", "elev", "slope", "vrm", "visib"])
        for i, (x, y, s_) in enumerate(pts):
            w.writerow([i, f"{s_:.1f}", f"{x:.1f}", f"{y:.1f}",
                        f"{G['DEM_fused'].value_at_xy(x, y, method='nearest'):.2f}",
                        f"{G['slope_deg'].value_at_xy(x, y, method='nearest'):.3f}",
                        f"{G['vrm_5x5'].value_at_xy(x, y, method='nearest'):.5f}",
                        f"{G['AVGVISIB_probability'].value_at_xy(x, y, method='nearest'):.5f}"])
    stats_nn, rows_nn = run_esim(nn_profile, 'nearest')
    variants.append(('最近邻敏感性变体', stats_nn, rows_nn))

    # 逐项对比
    for name, st, rows_ in variants:
        print(f'[{name}]')
        checks = []
        d_soc = st['min_soc'] - BASE['min_soc']
        ok1 = abs(d_soc) <= TOL['min_soc'] and st['min_soc'] > 200
        checks.append(('最低SoC', st['min_soc'], BASE['min_soc'],
                       f'Δ{d_soc:+.1f} Wh（±5 且 >200）', ok1))
        ok2 = abs(st['dwell'] - BASE['dwell']) <= TOL['dwell']
        checks.append(('驻留充电次数', st['dwell'], BASE['dwell'], '±1', ok2))
        ok3 = abs(st['lunar'] - BASE['lunar']) <= TOL['lunar']
        checks.append(('单程月昼数', st['lunar'], BASE['lunar'], '±1', ok3))
        ok_engineering = (st['min_soc'] > 200 and st['dwell'] < 15 and st['lunar'] < 10)
        checks.append(('工程验收三元组（>200/<15/<10）',
                       f"SoC {st['min_soc']} / 驻留 {st['dwell']} / 月昼 {st['lunar']}",
                       '全过', '结论不翻盘', ok_engineering))
        for cname, got, ref, tol, ok in checks:
            print(f"    [{'PASS' if ok else 'FAIL'}] {cname}: 实测 {got} / 基准 {ref}（{tol}）")
            if not ok and name.startswith('交付'):
                failures.append(f'{name}-{cname}')
        metrics[name] = dict(stats=st,
                             engineering_verdict=bool(ok_engineering),
                             min_soc_ok=bool(ok1), dwell_ok=bool(ok2),
                             lunar_ok=bool(ok3))
        print()

    # 对比图：两种口径 SoC 曲线 + 200Wh 红线
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = ['crimson', 'steelblue']
    for (name, st, rows_), color in zip(variants, colors):
        ss = [r['s'] / 1000 for r in rows_]
        vv = [r['soc'] for r in rows_]
        ax.plot(ss, vv, color=color, lw=1.8,
                label=f"{name}：min SoC {st['min_soc']:.1f} Wh，驻留 {st['dwell']}，{st['lunar']} 月昼")
    ax.axhline(200, color='k', ls='--', lw=1.5, label='E_min = 200 Wh（热控底线）')
    ax.axhline(BASE['min_soc'], color='gray', ls=':', lw=1.2,
               label='历史基准 min SoC = 215.2 Wh')
    ax.set_xlabel('里程 (km)')
    ax.set_ylabel('SoC (Wh)')
    ax.set_title('L3 端到端能量仿真：纯 Python 重构链 vs 历史基准口径')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig_path = os.path.join(FIG_DIR, 'L3_energy_comparison.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f'对比图已保存: {fig_path}')

    # 清理验证临时目录
    for tag in ('_esim_delivered', '_esim_nearest'):
        p = os.path.join(_HERE, 'data', 'output', tag)
        if os.path.isdir(p):
            shutil.rmtree(p, ignore_errors=True)

    metrics['failures'] = failures
    with open(os.path.join(_HERE, 'tests', 'L3_metrics.json'), 'w',
              encoding='utf-8') as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    print('\n' + '=' * 70)
    verdict = ('PASS' if not failures else
               '部分对齐（工程结论不翻盘；驻留/月昼绝对次数受取值口径影响，根因见报告）')
    print(f'L3 汇总: {verdict}')
    print('=' * 70)
    return 0 if not failures else 1


if __name__ == '__main__':
    sys.exit(main())
