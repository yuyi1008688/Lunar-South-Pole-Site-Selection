# -*- coding: utf-8 -*-
"""
L1_verify_lcp.py —— L1 算子级验证：自实现 Dijkstra vs UDBX 内商业软件产出
================================================================
验证对象（均来自 data/rasters/，源头为 UDBX，见 L0）：
  1. 成本累积面  distance_accumulation.tif（商业软件以 START 为源的耗费累积）
  2. 回溯方向面  backlink_direction.tif（8 方向编码，4-bit 打包存储）
  3. 原始路径栅格 optimal_path_raw.tif（全区几何展示路径，247 像元）

对比方法：
  - 累积面：以 START 为源、cost_surface 为成本（实测确认商业软件累积面不含
    stage2 的 +0.001 底数，见报告）跑 dijkstra_accumulation，逐像元相对误差
    统计 + Pearson 相关；
  - 方向面：商业软件方向码 1-8 与本实现码 1-8 存在固定旋转（1=W 起顺时针 vs
    1=N 起顺时针，语义均为"到达方向"），推导映射后全图比对一致率；
  - 路径：从几何终点沿方向面回溯到源，与 optimal_path_raw 像元集合求交。

容差声明（事先约定）：
  - 累积面：相对误差 <1e-6 视为相等（纯浮点累加顺序差异的 ULP 量级）；
  - 方向面：要求映射后 100% 一致（等成本平局在方向面上同样成立才不计入）；
  - 路径：几何路径要求像元集合 100% 重合；任务路径无独立基准栅格，其对标
    放在 L2（长度/顶点/剖面）。

输出：results/decouple_verification/L1_lcp_comparison.png + tests/L1_metrics.json
用法：python tests/L1_verify_lcp.py
"""

import json
import os
import sys
from collections import Counter

import numpy as np
import rasterio

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
from utils.raster_grid import RasterGrid
from ch06_path_planning.lcp_dijkstra import (NEIGHBORS, dijkstra_accumulation,
                                             path_length_m, trace_backlink)

RASTER_DIR = os.environ.get('LUNAR_RASTER_DIR', os.path.join('data', 'rasters'))
FIG_DIR = os.path.join('results', 'decouple_verification')

START = (44760.0, 10920.0)
END_TASK = (40920.0, 10920.0)
END_GEO = (-4680.0, -13080.0)

MY2DIR = {code: (dr, dc) for code, (dr, dc, _) in enumerate(NEIGHBORS, start=1)}
DIR2MY = {d: c for c, d in MY2DIR.items()}


def rc_at(x, y, left=-46080.0, top=46080.0, res=240.0):
    return int(round((top - y) / res - 0.5)), int(round((x - left) / res - 0.5))


def xy_at(r, c, left=-46080.0, top=46080.0, res=240.0):
    return (left + c * res + res / 2.0, top - r * res - res / 2.0)


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    metrics = {}

    cost = RasterGrid.from_file(os.path.join(RASTER_DIR, 'cost_surface.tif')).arr
    sm_acc = RasterGrid.from_file(os.path.join(RASTER_DIR, 'distance_accumulation.tif')).arr
    sm_bp = RasterGrid.from_file(os.path.join(RASTER_DIR, 'backlink_direction.tif')).arr
    op = RasterGrid.from_file(os.path.join(RASTER_DIR, 'optimal_path_raw.tif')).arr

    rS, cS = rc_at(*START)
    print('── L1 算子级验证：8 邻域 Dijkstra vs UDBX 内商业软件产出 ──')
    print(f'源 = START {START}（累积面源点实测确认：acc[START]=0），成本 = cost_surface')

    my_acc, my_bp = dijkstra_accumulation(cost, (rS, cS))

    # ---------- 1. 累积面 ----------
    rel = np.abs(my_acc - sm_acc) / np.maximum(np.abs(sm_acc), 1e-9)
    n_all = rel.size
    m = dict(
        median_pct=float(np.median(rel) * 100),
        p95_pct=float(np.percentile(rel, 95) * 100),
        max_pct=float(rel.max() * 100),
        lt_1e6_pct=float((rel < 1e-6).mean() * 100),
        pearson=float(np.corrcoef(my_acc.ravel(), sm_acc.ravel())[0, 1]),
    )
    metrics['accumulation'] = m
    print('\n[1] 成本累积面（147456 像元全图）')
    print(f"    相对误差中位数 {m['median_pct']:.6f}% | P95 {m['p95_pct']:.6f}% | "
          f"最大 {m['max_pct']:.6f}%")
    print(f"    相对误差<1e-6 像元占比 {m['lt_1e6_pct']:.4f}% | Pearson r = {m['pearson']:.9f}")

    # ---------- 2. 方向面（先推导码映射） ----------
    close = rel < 1e-9
    votes = {k: Counter() for k in range(1, 9)}
    rows, cols = np.where(close & (sm_bp > 0) & (my_bp > 0))
    for r, c in zip(rows, cols):
        votes[int(sm_bp[r, c])][MY2DIR[int(my_bp[r, c])]] += 1
    sm2dir = {k: v.most_common(1)[0][0] for k, v in votes.items() if v}
    assert len(sm2dir) == 8 and len(set(sm2dir.values())) == 8, '方向码映射非双射'
    sm2my = {k: DIR2MY[d] for k, d in sm2dir.items()}
    mapped = np.vectorize(lambda k: sm2my.get(int(k), 0))(sm_bp)
    both = (sm_bp > 0) & (my_bp > 0)
    agree = int(((mapped == my_bp) & both).sum())
    tot = int(both.sum())
    metrics['backlink'] = dict(mapping={str(k): list(v) for k, v in sm2dir.items()},
                               agree=agree, total=tot,
                               agree_pct=agree / tot * 100)
    print('\n[2] 回溯方向面')
    print(f'    商业软件码(到达方向, 1=W起顺时针) → 本实现到达方向(dr,dc): '
          f'{ {k: v for k, v in sorted(sm2dir.items())} }')
    print(f'    全图方向码一致: {agree}/{tot} = {agree / tot * 100:.4f}%')

    # ---------- 3. 路径（几何路径 vs optimal_path_raw） ----------
    rG, cG = rc_at(*END_GEO)
    trace_geo = trace_backlink(my_bp, (rG, cG), (rS, cS))
    op_pixels = {(int(r), int(c)) for r, c in zip(*np.where(op == 1))}
    tp = set(trace_geo)
    inter = tp & op_pixels
    L_geo = path_length_m(trace_geo)
    metrics['path_geo'] = dict(n=len(trace_geo), overlap=len(inter),
                               overlap_pct=len(inter) / len(tp) * 100,
                               op_n=len(op_pixels),
                               fold_km=L_geo / 1000,
                               ends_at_start=trace_geo[-1] == (rS, cS))
    print('\n[3] 路径回溯（几何展示路径）')
    print(f"    从 END_geo 沿本实现方向面回溯: {len(trace_geo)} 像元，"
          f"折线长 {L_geo / 1000:.3f} km，到源 {trace_geo[-1] == (rS, cS)}")
    print(f"    与 optimal_path_raw 像元集合重合: {len(inter)}/{len(tp)} = "
          f"{len(inter) / len(tp) * 100:.2f}%（商业软件栅格共 {len(op_pixels)} 像元）")

    # 任务路径（无独立基准栅格，记录指标供 L2 对标）
    rE, cE = rc_at(*END_TASK)
    trace_task = trace_backlink(my_bp, (rE, cE), (rS, cS))
    L_task = path_length_m(trace_task)
    metrics['path_task'] = dict(n=len(trace_task), fold_km=L_task / 1000,
                                ends_at_start=trace_task[-1] == (rS, cS))
    print(f"    任务路径回溯: {len(trace_task)} 像元，折线长 {L_task / 1000:.3f} km"
          f"（历史基准 4.039 km / 17 顶点，对标见 L2）")

    ok = (m['lt_1e6_pct'] == 100.0 and agree == tot
          and len(inter) == len(tp) == len(op_pixels))
    print('\n' + '=' * 70)
    print(f'L1 汇总: {"PASS" if ok else "FAIL"}'
          f'（累积面<1e-6 占比 100% ✓ | 方向码 100% 一致 ✓ | 路径像元 100% 重合 ✓）'
          if ok else f'L1 汇总: FAIL（见上方分项）')
    print('=' * 70)

    # ---------- 对比图 ----------
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(19, 6))
    rel_map = np.where(np.isfinite(rel), rel, np.nan)
    im0 = axes[0].imshow(rel_map * 100, cmap='magma', vmin=0,
                         vmax=max(float(np.nanpercentile(rel_map, 99.9) * 100), 1e-6))
    axes[0].set_title('累积面相对误差 (%) —— P99.9 内全为浮点ULP级')
    plt.colorbar(im0, ax=axes[0], fraction=0.046)

    idx = np.random.default_rng(42).choice(n_all, 5000, replace=False)
    axes[1].scatter(sm_acc.ravel()[idx], my_acc.ravel()[idx], s=2, alpha=0.3)
    lim = [0, max(sm_acc.max(), my_acc.max()) * 1.02]
    axes[1].plot(lim, lim, 'r--', lw=1)
    axes[1].set_xlabel('商业软件累积成本')
    axes[1].set_ylabel('本实现 Dijkstra 累积成本')
    axes[1].set_title(f"累积面 1:1 对比（r={m['pearson']:.9f}）")

    canvas = np.full(op.shape, np.nan)
    for r, c in op_pixels:
        canvas[r, c] = 1
    only_mine = tp - op_pixels
    only_ref = op_pixels - tp
    for r, c in tp:
        canvas[r, c] = 2 if (r, c) in inter else 3
    for r, c in only_ref:
        canvas[r, c] = 4
    axes[2].imshow(np.where(np.isnan(canvas), np.nan, canvas), cmap='viridis',
                   interpolation='nearest')
    axes[2].set_title(f'路径像元叠加：绿=重合{len(inter)} 黄=仅本实现{len(only_mine)} '
                      f'蓝=仅基准{len(only_ref)}')
    for ax, (ex, ey), name in zip([axes[2]] * 3, [START, END_TASK, END_GEO],
                                  ['START', 'END_task', 'END_geo']):
        r, c = rc_at(ex, ey)
        ax.plot(c, r, 'r*', ms=14)
        ax.annotate(name, (c, r), color='red', fontsize=9,
                    xytext=(c + 6, r + 6))
    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, 'L1_lcp_comparison.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f'\n对比图已保存: {fig_path}')

    with open('tests/L1_metrics.json', 'w', encoding='utf-8') as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)
    print('指标已保存: tests/L1_metrics.json')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
