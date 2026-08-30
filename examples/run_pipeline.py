# -*- coding: utf-8 -*-
"""
run_pipeline.py —— 全链路纯 Python 一键总控（阶段 DAG 编排器）
================================================================
Ch01–Ch06 主分析链已 100% 纯 Python（开源栈：标准库 + numpy + scipy +
rasterio），本脚本把它编排为**单命令端到端流水线**：

  python examples/run_pipeline.py --list              # 查看 DAG 与数据血缘
  python examples/run_pipeline.py --all               # 端到端全跑（断点续跑）
  python examples/run_pipeline.py --run ch05          # 只跑某章及其依赖
  python examples/run_pipeline.py --from ch05         # 从某章续跑到结束
  python examples/run_pipeline.py --all --force       # 忽略已有产物强制重算
  python examples/run_pipeline.py --all --dry-run     # 只打印将执行什么

特性：
  - 依赖自动满足：跑某阶段前检查输入产物，缺失则自动先跑上游；
  - 断点续跑：输出已存在且未加 --force 则跳过（up-to-date）；
  - 缺数据优雅跳过：如 Ch03 需要的 wang 冰点 xlsx 未按 data/README.md 布局
    准备时，该阶段标记 SKIPPED 并说明，不阻塞主链；
  - 新写模块（Ch02/Ch05）函数级调用（同进程、零额外 IO）；历史 __main__
    脚本 subprocess 调用，捕获退出码、失败即停并定位阶段；
  - 每阶段产物自动过铁基准断言（384×384@240m@±46080，见 src/utils/iron_grid.py）；
  - 运行结束写 data/output/manifest.json：各阶段输入/输出文件、字节数、
    sha256、min/max/mean、耗时、校验结果——论文可复现的产物血缘台账。

不依赖任何商业 GIS 软件，也不启动任何 MCP/外部服务。
"""

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time

sys.stdout.reconfigure(encoding='utf-8')

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join('data', 'rasters')          # 栅格库（udbx_extract 导出）
O = os.path.join('data', 'output')           # 产物根
UDBX_DEFAULT = os.path.join('..', 'A+ak9HzkqIZG', '2成果数据', '数据源', '输入.udbx')

sys.path.insert(0, os.path.join(REPO, 'src'))
sys.path.insert(0, os.path.join(REPO, 'src', 'utils'))


def _py(*rel):
    return os.path.join(REPO, *rel)


def _load_module(rel_path, name):
    spec = importlib.util.spec_from_file_location(name, _py(rel_path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# 阶段 DAG 定义：name → dict(deps, inputs, outputs, runner, desc)
# runner 两种形态：
#   ('func', 模块相对路径, 模块名, 可调用名, 参数字典)      —— 函数级调用
#   ('sub', 脚本相对路径, 环境变量字典, 相对工作目录)      —— 子进程调用
# ---------------------------------------------------------------------------

def _A(p):
    """runner 环境变量用的绝对路径（子进程 cwd 可能不同）。"""
    return os.path.abspath(os.path.join(REPO, p))


def build_stages():
    ch06_env = {'LUNAR_RASTER_DIR': _A(R), 'LUNAR_OUTPUT_DIR': _A(O),
                'MPLBACKEND': 'Agg'}
    return {
        'data': dict(
            deps=[], desc='UDBX → GeoTIFF 栅格解放（35 个，含逐栅格自检）',
            inputs=[os.environ.get('LUNAR_UDBX', UDBX_DEFAULT)],
            outputs=[os.path.join(R, n) for n in (
                'slope_deg.tif', 'AVGVISIB_probability.tif', 'combined_hazard.tif',
                'cost_surface.tif', 'DEM_fused.tif', 'vrm_5x5.tif',
                'distance_accumulation.tif', 'backlink_direction.tif',
                'optimal_path_raw.tif')],
            runner=('sub', 'src/utils/udbx_extract.py',
                    {'LUNAR_UDBX': _A(os.environ.get('LUNAR_UDBX', UDBX_DEFAULT)),
                     'LUNAR_RASTER_DIR': _A(R), 'MPLBACKEND': 'Agg'}, REPO),
        ),
        'ch01': dict(
            deps=['data'], desc='Ch01 数据底座：融合 DEM 山体阴影（4 方位角）',
            inputs=[os.path.join(R, 'DEM_fused.tif')],
            outputs=[os.path.join(O, 'ch01', f'hillshade_az{az:03d}_alt15.tif')
                     for az in (315, 45, 135, 225)],
            runner=('sub', 'src/ch01_data_foundation/run_hillshade.py',
                    {'LUNAR_DATA_DIR': _A(R), 'LUNAR_OUTPUT_DIR': _A(os.path.join(O, 'ch01')),
                     'MPLBACKEND': 'Agg'},
                    REPO),
        ),
        'ch02': dict(
            deps=['data'], desc='Ch02 光照分类：双 PSR / 五级分类 / F1 归一化（纯 Python）',
            inputs=[os.path.join(R, 'AVGVISIB_probability.tif')],
            outputs=[os.path.join(O, 'ch02', 'C2_illumination_class.tif'),
                     os.path.join(O, 'ch02', 'F1_illumination_norm.tif')],
            runner=('func', 'src/ch02_illumination/illumination.py', 'illumination',
                    'run', dict(input_tif=os.path.join(R, 'AVGVISIB_probability.tif'),
                                out_class=os.path.join(O, 'ch02', 'C2_illumination_class.tif'),
                                out_f1=os.path.join(O, 'ch02', 'F1_illumination_norm.tif'))),
        ),
        'ch03': dict(
            deps=['data'],
            desc='Ch03 水冰 KDE（可选）：需先把 wang 冰点 xlsx/PSR 掩膜按 data/ch03_staging/00_数据源/ 布局就位'
                 '（xlsx 在提交包 wang_ice_point/ 内；就位后本阶段自动纳入运行）',
            inputs=[os.path.join('data', 'ch03_staging', '00_数据源', '最新数据',
                                 'wang2025_ice_pixel_positions_spectra.xlsx'),
                    os.path.join('data', 'ch03_staging', '00_数据源', '最新数据', 'AVGVISIB_probability.tif'),
                    os.path.join('data', 'ch03_staging', '00_数据源', '最新数据', 'sPSR_mask.tif'),
                    os.path.join('data', 'ch03_staging', '00_数据源', '最新数据', 'subPSR_mask.tif'),
                    os.path.join('data', 'ch03_staging', '00_数据源', 'PSR_mask.tif')],
            outputs=[os.path.join(O, 'ch03', 'F2_wang_kde_final.tif')],
            runner=('sub', 'src/ch03_water_ice/wang_kde_to_f2.py',
                    {'LUNAR_PROJECT_ROOT': _A(os.path.join('data', 'ch03_staging')),
                     'MPLBACKEND': 'Agg'},
                    REPO),
        ),
        'ch04_fos': dict(
            deps=['data'], desc='Ch04 安全势场：FoS 25 组参数敏感性（C×φ 网格）',
            inputs=[os.path.join(R, 'slope_deg.tif')],
            outputs=[os.path.join(O, 'ch04_sensitivity', 'sensitivity_area_km2.csv')],
            runner=('sub', 'src/ch04_safety/fos_sensitivity.py',
                    {'LUNAR_SLOPE_TIF': _A(os.path.join(R, 'slope_deg.tif')),
                     'LUNAR_OUTPUT_DIR': _A(os.path.join(O, 'ch04_sensitivity')),
                     'MPLBACKEND': 'Agg'}, REPO),
        ),
        'ch04_dist': dict(
            deps=['data'], desc='Ch04 安全势场：F5 连续危险距离场（欧氏距离变换）',
            inputs=[os.path.join(R, 'combined_hazard.tif'), os.path.join(R, 'slope_deg.tif')],
            outputs=[os.path.join(O, 'ch04', 'distance_raw.tif')],
            runner=('sub', 'src/ch04_safety/euclidean_distance.py',
                    {'LUNAR_HAZARD_TIF': _A(os.path.join(R, 'combined_hazard.tif')),
                     'LUNAR_SLOPE_TIF': _A(os.path.join(R, 'slope_deg.tif')),
                     'MPLBACKEND': 'Agg'},
                    _A(os.path.join(O, 'ch04'))),
        ),
        'ch04_ecsa': dict(
            deps=['data'], desc='Ch04 安全势场：ECSA 独立性诊断（需对地可见概率原始 jp2）',
            inputs=[os.path.join('..', 'A+ak9HzkqIZG', '1原始数据', 'EARTH_look',
                                 'cartorder', 'browse_extras', 'avgvisib_65s_240m_earth.jp2'),
                    os.path.join(R, 'AVGVISIB_probability.tif')],
            outputs=[os.path.join(O, 'ch04', 'step0_ecsa_diagnostic_report.txt'),
                     os.path.join(O, 'ch04', 'step0_ecsa_result.json')],
            runner=('sub', 'src/ch04_safety/step0_ecsa_diagnostic.py',
                    {'LUNAR_LIGHT_TIF': _A(os.path.join(R, 'AVGVISIB_probability.tif')),
                     'MPLBACKEND': 'Agg',
                     'LUNAR_EARTH_TIF': _A(os.path.join(
                         '..', 'A+ak9HzkqIZG', '1原始数据', 'EARTH_look', 'cartorder',
                         'browse_extras', 'avgvisib_65s_240m_earth.jp2')),
                     'LUNAR_OUTPUT_DIR': _A(os.path.join(O, 'ch04'))}, REPO),
        ),
        'ch05': dict(
            deps=['ch02'], desc='Ch05 AHP-WLC 综合成图：加权→硬约束→Jenks 五级→Ⅰ级点（纯 Python）',
            inputs=[os.path.join(O, 'ch02', 'F1_illumination_norm.tif'),
                    os.path.join(R, 'F2_wang_kde_final_1.tif'),
                    os.path.join(R, 'F3_fos_safety.tif'),
                    os.path.join(R, 'F4_ecsa_sync.tif'),
                    os.path.join(R, 'F5_continuous_distance.tif'),
                    os.path.join(R, 'slope_deg.tif'),
                    os.path.join(R, 'AVGVISIB_probability.tif'),
                    os.path.join(R, 'hazard_distance.tif')],
            outputs=[os.path.join(O, 'ch05', 'suitability_final.tif'),
                     os.path.join(O, 'ch05', 'suitability_classes.tif'),
                     os.path.join(O, 'ch05', 'constraint_mask.tif'),
                     os.path.join(O, 'ch05', 'grade_I_points.csv')],
            runner=('func', 'src/ch05_ahp_site/ahp_wlc.py', 'ahp_wlc', 'run',
                    dict(f1_dir=os.path.join(O, 'ch02'), out_dir=os.path.join(O, 'ch05'))),
        ),
        'ch06_stage1': dict(
            deps=['ch05'], desc='Ch06 阶段一：参数定标（D_min/D_max、F2 阈值分布）',
            inputs=[os.path.join(R, 'ice_density_final.tif'),
                    os.path.join(R, 'slope_deg.tif'),
                    os.path.join(R, 'combined_hazard.tif'),
                    os.path.join(R, 'AVGVISIB_probability.tif')],
            outputs=[os.path.join(O, 'stage1_params.txt')],
            runner=('sub', 'src/ch06_path_planning/stage1_select_mining_target.py',
                    ch06_env, REPO),
        ),
        'ch06_validate': dict(
            deps=['ch05'], desc='Ch06 阶段一：终点 (40920,10920) 六项合规验证（控制台判定）',
            inputs=[os.path.join(R, 'slope_deg.tif'),
                    os.path.join(R, 'combined_hazard.tif'),
                    os.path.join(R, 'AVGVISIB_probability.tif'),
                    os.path.join(R, 'cost_surface.tif')],
            outputs=[],   # 判定型阶段，无文件产物
            runner=('sub', 'src/ch06_path_planning/stage1_validate_endpoint.py',
                    ch06_env, REPO),
        ),
        'ch06_stage2': dict(
            deps=['ch06_stage1', 'ch06_validate'],
            desc='Ch06 阶段二主链：四掩膜→candidate→Dijkstra LCP→B样条→50m 剖面',
            inputs=[os.path.join(R, 'slope_deg.tif'),
                    os.path.join(R, 'AVGVISIB_probability.tif'),
                    os.path.join(R, 'combined_hazard.tif'),
                    os.path.join(R, 'cost_surface.tif'),
                    os.path.join(R, 'DEM_fused.tif'),
                    os.path.join(R, 'vrm_5x5.tif')],
            outputs=[os.path.join(O, 'stage2', 'optimal_path_vertices.csv'),
                     os.path.join(O, 'stage2', 'optimal_path_wkt.txt'),
                     os.path.join(O, 'stage2', 'path_profile.csv')],
            runner=('sub', 'src/ch06_path_planning/stage2_purepython.py', ch06_env, REPO),
        ),
        'ch06_energy': dict(
            deps=['ch06_stage2'], desc='Ch06 阶段三：三态状态机能量递推（表 6-3 + 能量曲线）',
            inputs=[os.path.join(O, 'stage2', 'path_profile.csv')],
            outputs=[os.path.join(O, 'stage3_energy', 'energy_curve.csv'),
                     os.path.join(O, 'stage3_energy', 'energy_curve.png')],
            runner=('sub', 'src/ch06_path_planning/energy_simulation.py', ch06_env, REPO),
        ),
        'ch06_svg': dict(
            deps=['ch06_energy'], desc='Ch06 零依赖能量曲线 SVG（可选）',
            inputs=[os.path.join(O, 'stage3_energy', 'energy_curve.csv')],
            outputs=[os.path.join(O, 'stage3_energy', 'energy_curve.svg')],
            runner=('sub', 'src/ch06_path_planning/plot_energy_curve_svg.py', ch06_env, REPO),
        ),
        'ch08': dict(
            deps=['ch06_stage2'],
            desc='Ch08 数字孪生：浏览器直接打开 src/ch08_digital_twin/ThreeJS_scene/index.html（无 Python 依赖，不进自动链）',
            inputs=[], outputs=[], runner=None,
        ),
    }

ORDER = ['data', 'ch01', 'ch02', 'ch03', 'ch04_fos', 'ch04_dist', 'ch04_ecsa',
         'ch05', 'ch06_stage1', 'ch06_validate', 'ch06_stage2', 'ch06_energy',
         'ch06_svg', 'ch08']


def sha256_of(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(chunk), b''):
            h.update(blk)
    return h.hexdigest()


def tif_stats(path):
    import rasterio
    import numpy as np
    with rasterio.open(path) as src:
        a = src.read(1).astype(np.float64)
        nd = src.nodata
    valid = (a != nd) & np.isfinite(a) if nd is not None else np.isfinite(a)
    v = a[valid]
    return dict(min=float(v.min()), max=float(v.max()), mean=float(v.mean())) if v.size else {}


def file_entry(path):
    e = dict(path=path.replace('\\', '/'), bytes=os.path.getsize(path),
             sha256=sha256_of(path))
    if path.lower().endswith('.tif'):
        try:
            e.update(tif_stats(path))
        except Exception:
            pass
    return e


def outputs_up_to_date(st):
    return all(os.path.exists(p) for p in st['outputs'])


def inputs_present(st):
    missing = [p for p in st['inputs'] if not os.path.exists(p)]
    return (len(missing) == 0), missing


def topo_order(start=None, stages=None):
    """从 start（含）出发的依赖闭包拓扑序；start=None 时全图。"""
    stages = stages or build_stages()
    seen, out = set(), []

    def visit(n):
        if n in seen:
            return
        seen.add(n)
        for d in stages[n]['deps']:
            visit(d)
        out.append(n)
    if start is None:
        for n in ORDER:
            visit(n)
    else:
        visit(start)
    return [n for n in ORDER if n in seen]


def run_stage(name, st, force, manifest):
    t0 = time.perf_counter()
    rec = dict(stage=name, desc=st['desc'], status=None, inputs=[], outputs=[],
               duration_s=None)
    ok_inputs, missing = inputs_present(st)
    if not ok_inputs:
        rec['status'] = 'SKIPPED'
        rec['missing_inputs'] = missing
        print(f"  [SKIPPED] {name}: 缺输入 {len(missing)} 个（如 {os.path.basename(missing[0])}）"
              f" —— {st['desc']}")
        manifest[name] = rec
        return rec
    if not force and outputs_up_to_date(st) and st['outputs']:
        rec['status'] = 'UPTODATE'
        rec['outputs'] = [file_entry(p) for p in st['outputs']]
        print(f"  [uptodate] {name}: 产物已存在（--force 可重算）")
        manifest[name] = rec
        return rec

    print(f"  [RUN     ] {name}: {st['desc']}")
    kind = st['runner'][0]
    try:
        if kind == 'func':
            _, mod_path, mod_name, fn_name, kwargs = st['runner']
            mod = _load_module(mod_path, 'pipe_' + mod_name)
            getattr(mod, fn_name)(**kwargs)
        else:
            _, script, env, cwd = st['runner']
            if cwd:
                os.makedirs(cwd, exist_ok=True)
            e = dict(os.environ)
            e.update({k: v for k, v in env.items() if v is not None})
            r = subprocess.run([sys.executable, _py(script)], env=e, cwd=cwd or REPO)
            if r.returncode != 0:
                raise RuntimeError(f"子进程退出码 {r.returncode}")
        # 产物存在性 + 铁基准校验
        from iron_grid import assert_iron_grid, IronGridViolation
        import rasterio as _rio
        for p in st['outputs']:
            if not os.path.exists(p):
                raise RuntimeError(f"声明产物未生成: {p}")
            if p.lower().endswith('.tif'):
                with _rio.open(p) as src:
                    assert_iron_grid(f'{name}/{os.path.basename(p)}',
                                     src.read(1).astype('float64'),
                                     src.transform, src.crs,
                                     nodata=src.nodata if src.nodata is not None else -9999.0)
        rec['status'] = 'OK'
    except Exception as e:
        rec['status'] = 'FAILED'
        rec['error'] = f'{type(e).__name__}: {e}'
        manifest[name] = rec
        print(f"  [FAILED  ] {name}: {type(e).__name__}: {e}")
        raise
    rec['duration_s'] = round(time.perf_counter() - t0, 2)
    rec['inputs'] = [file_entry(p) for p in st['inputs']]
    rec['outputs'] = [file_entry(p) for p in st['outputs']]
    manifest[name] = rec
    print(f"  [OK      ] {name}: {len(st['outputs'])} 产物, {rec['duration_s']}s")
    return rec


def chain_status(names, force, stages):
    """预演各阶段状态（--dry-run / --list 用）。"""
    rows = []
    for n in names:
        st = stages[n]
        ok, missing = inputs_present(st)
        if not st['runner']:
            stat = 'HINT'
        elif not ok:
            stat = f'SKIPPED(缺{len(missing)}输入)'
        elif not force and outputs_up_to_date(st) and st['outputs']:
            stat = 'uptodate'
        else:
            stat = 'will-run'
        rows.append((n, stat))
    return rows


def main():
    ap = argparse.ArgumentParser(
        description='月球南极选址 · 全链路纯 Python 一键总控（无商业 GIS 依赖）')
    ap.add_argument('--list', action='store_true', help='列出阶段 DAG 与数据血缘')
    ap.add_argument('--all', action='store_true', help='端到端全跑（断点续跑）')
    ap.add_argument('--run', metavar='CH', help='只跑某阶段及其依赖（如 ch05）')
    ap.add_argument('--from', dest='frm', metavar='CH', help='从某阶段续跑到结束')
    ap.add_argument('--force', action='store_true', help='忽略已有产物强制重算')
    ap.add_argument('--dry-run', action='store_true', help='只打印将执行什么，不实际运行')
    args = ap.parse_args()

    stages = build_stages()

    if args.list or not (args.all or args.run or args.frm):
        print('=' * 78)
        print('月球南极选址 · 全链路纯 Python 阶段 DAG（铁基准 384×384@240m@±46080）')
        print('=' * 78)
        for n in ORDER:
            st = stages[n]
            ins = ', '.join(os.path.basename(p) for p in st['inputs']) or '—'
            outs = ', '.join(os.path.basename(p) for p in st['outputs']) or '—'
            print(f"\n[{n}] {st['desc']}")
            print(f"    依赖: {', '.join(st['deps']) or '—'}")
            print(f"    输入: {ins}")
            print(f"    输出: {outs}")
        print('\n用法: --all / --run CH / --from CH / --force / --dry-run')
        return 0

    if args.run:
        names = topo_order(args.run, stages)
    elif args.frm:
        all_names = topo_order(None, stages)
        names = all_names[all_names.index(args.frm):]
    else:
        names = topo_order(None, stages)

    plan = chain_status(names, args.force, stages)
    print('=' * 78)
    print('执行计划（依赖闭包，拓扑序）:')
    for n, s in plan:
        print(f'  [{s:<16}] {n}')
    print('=' * 78)
    if args.dry_run:
        print('--dry-run：未实际执行。')
        return 0

    manifest = {}
    if os.path.exists(os.path.join(O, 'manifest.json')) and not args.force:
        try:
            with open(os.path.join(O, 'manifest.json'), encoding='utf-8') as f:
                old = json.load(f)
            if isinstance(old, dict):
                manifest.update(old.get('stages', {}))
        except Exception:
            pass
    t0 = time.perf_counter()
    for n in names:
        st = stages[n]
        if not st['runner']:
            print(f"\n[HINT] {n}: {st['desc']}")
            manifest[n] = dict(stage=n, status='HINT', desc=st['desc'])
            continue
        print(f"\n── 阶段 {n} ──")
        rec = run_stage(n, st, args.force, manifest)
        if rec['status'] == 'FAILED':
            with open(os.path.join(O, 'manifest.json'), 'w', encoding='utf-8') as f:
                json.dump(dict(started=None, stages=manifest), f, ensure_ascii=False, indent=2)
            sys.exit(f"\n✗ 阶段 {n} 失败，流水线中止（详见上方定位信息与 manifest.json）")

    total = round(time.perf_counter() - t0, 1)
    out_manifest = dict(
        project='月球南极 Shackleton 环形山科研站选址（纯 Python 全链路）',
        generated_by='examples/run_pipeline.py --all',
        python=sys.version.split()[0],
        total_duration_s=total,
        stages=manifest)
    os.makedirs(O, exist_ok=True)
    with open(os.path.join(O, 'manifest.json'), 'w', encoding='utf-8') as f:
        json.dump(out_manifest, f, ensure_ascii=False, indent=2)
    print('\n' + '=' * 78)
    print(f'★ 流水线完成：{len(names)} 阶段 / 总耗时 {total}s')
    print(f'  产物血缘台账 -> {os.path.join(O, "manifest.json")}')
    print('  数字孪生：浏览器直接打开 src/ch08_digital_twin/ThreeJS_scene/index.html')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
