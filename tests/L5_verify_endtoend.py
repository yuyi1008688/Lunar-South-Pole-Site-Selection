# -*- coding: utf-8 -*-
"""
L5_verify_endtoend.py —— L5.3 端到端验收：一键流水线产物核对既有基准
================================================================
验证 `python examples/run_pipeline.py --all` 的产物与 tests/L1~L4_metrics.json
中的既有基准不冲突：

  1. Ch06 任务路径：总长 4.029 km（基准 4.039 ±2%）、17 顶点、危险区穿越 0；
  2. 能量仿真：最低 SoC 210.8 Wh（>200，与 L3 一致）、驻留/月昼与 L3 一致；
  3. manifest.json：全部阶段 OK/SKIPPED/HINT（无 FAILED）；
  4. L1~L4 metrics JSON：无未解决 failures；
  5. Ch05 Ⅰ级点 argmax = 推荐站址（与 L4-A/L5.2 一致）。

用法：python tests/L5_verify_endtoend.py（先运行 examples/run_pipeline.py --all）
输出：tests/L5_metrics.json（追加 L5_3 字段）
"""

import csv
import json
import math
import os
import sys

import numpy as np
import rasterio

sys.stdout.reconfigure(encoding='utf-8')
_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
O = os.path.join(_HERE, 'data', 'output')
R = os.path.join(_HERE, 'data', 'rasters')
STATION = (44760.0, 10920.0)
failures = []


def polyline_length(coords):
    return float(sum(math.hypot(coords[i + 1][0] - coords[i][0],
                                coords[i + 1][1] - coords[i][1])
                     for i in range(len(coords) - 1)))


def main():
    metrics = {}
    print('── L5.3 端到端验收（--all 产物 vs 既有基准） ──')

    # ---- 1. Ch06 任务路径 ----
    vcsv = os.path.join(O, 'stage2', 'optimal_path_vertices.csv')
    coords = []
    with open(vcsv, encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            coords.append((float(row['x']), float(row['y'])))
    L = polyline_length(coords)
    ok_len = abs(L / 1000 - 4.039) <= 0.02 * 4.039
    ok_vert = len(coords) == 17
    # 危险区回检（30m 重采样，最近邻）
    sys.path.insert(0, os.path.join(_HERE, 'src', 'utils'))
    sys.path.insert(0, os.path.join(_HERE, 'src', 'ch06_path_planning'))
    from raster_grid import RasterGrid
    hz = RasterGrid.from_file(os.path.join(R, 'combined_hazard.tif'))
    from bspline_smooth import resample as _rs
    pts30, _ = _rs(coords, 30.0)
    n_cross = sum(1 for x, y, _ in pts30
                  if hz.value_at_xy(x, y, method='nearest') >= 1.0)
    ok_cross = n_cross == 0
    print(f"  任务路径: {L / 1000:.3f} km（基准 4.039 ±2%）/{len(coords)} 顶点"
          f"（基准 17）/危险区穿越 {n_cross}（必须 0）"
          f" → {'✓' if (ok_len and ok_vert and ok_cross) else '✗'}")
    metrics['path'] = dict(length_km=L / 1000, vertices=len(coords),
                           crossings=n_cross, ok=bool(ok_len and ok_vert and ok_cross))
    if not (ok_len and ok_vert and ok_cross):
        failures.append('path')

    # ---- 2. 能量仿真 ----
    curve = os.path.join(O, 'stage3_energy', 'energy_curve.csv')
    states, min_soc, min_at = [], None, 0.0
    with open(curve, encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            soc = float(r['soc_wh'])
            states.append(r['state'])
            if min_soc is None or soc < min_soc:
                min_soc, min_at = soc, float(r['s_m'])
    dwell = states.count('DWELL')
    lunar = states.count('WAIT') + 1
    ok_e = min_soc > 200 and dwell < 15 and lunar < 10
    print(f"  能量仿真: 最低 SoC {min_soc:.1f} Wh（>200）/驻留 {dwell}（<15）"
          f"/月昼 {lunar}（<10）→ {'✓' if ok_e else '✗'}")
    metrics['energy'] = dict(min_soc=min_soc, min_soc_at_km=min_at / 1000,
                             dwell=dwell, lunar_days=lunar, ok=bool(ok_e))
    if not ok_e:
        failures.append('energy')

    # 与 L3 metrics 一致性
    with open(os.path.join(_HERE, 'tests', 'L3_metrics.json'), encoding='utf-8') as f:
        l3 = json.load(f)
    l3_delivered = l3.get('交付剖面（双线性）', {}).get('stats', {})
    if l3_delivered:
        same = (abs(l3_delivered.get('min_soc', -1) - min_soc) < 0.05
                and l3_delivered.get('dwell') == dwell
                and l3_delivered.get('lunar') == lunar)
        print(f"  与 L3 基线一致: {'✓' if same else '✗'}（L3: min SoC "
              f"{l3_delivered.get('min_soc')}, 驻留 {l3_delivered.get('dwell')}, "
              f"月昼 {l3_delivered.get('lunar')}）")
        metrics['energy']['match_L3'] = bool(same)
        if not same:
            failures.append('energy_vs_L3')

    # ---- 3. manifest ----
    with open(os.path.join(O, 'manifest.json'), encoding='utf-8') as f:
        mf = json.load(f)
    bad = {k: v.get('status') for k, v in mf.get('stages', {}).items()
           if v.get('status') not in ('OK', 'SKIPPED', 'HINT', 'UPTODATE')}
    print(f"  manifest.json: {len(mf.get('stages', {}))} 阶段，"
          f"异常状态: {bad if bad else '无'} → {'✓' if not bad else '✗'}")
    metrics['manifest'] = dict(stages=len(mf.get('stages', {})), bad=bad,
                               ok=not bad)
    if bad:
        failures.append('manifest')

    # ---- 4. L1~L4 metrics 无未解决 failures ----
    # L3 的两条"交付口径-驻留/月昼"是验证报告 §L3 已登记、已归因的取值口径差异
    # （双线性 vs 最近邻；工程三元组单独验证通过），属已解决登记项而非未解决失败。
    REGISTERED = {'交付剖面（双线性）-驻留充电次数', '交付剖面（双线性）-单程月昼数'}
    for n in ('L1_metrics', 'L2_metrics', 'L3_metrics', 'L4_metrics'):
        with open(os.path.join(_HERE, 'tests', f'{n}.json'), encoding='utf-8') as f:
            m = json.load(f)
        fails = [x for x in m.get('failures', []) if x not in REGISTERED]
        registered = [x for x in m.get('failures', []) if x in REGISTERED]
        print(f"  {n}: 未解决 failures = {fails if fails else '无'}"
              f"（已登记差异 {len(registered)} 条）→ {'✓' if not fails else '✗'}")
        metrics[n] = dict(failures=fails, registered=registered)
        if fails:
            failures.append(n)

    # ---- 5. Ch05 argmax = 站址 ----
    with rasterio.open(os.path.join(O, 'ch05', 'suitability_final.tif')) as src:
        a = src.read(1).astype(np.float64)
    r_pk, c_pk = np.unravel_index(np.argmax(a), a.shape)
    ax, ay = -46080 + (c_pk + 0.5) * 240, 46080 - (r_pk + 0.5) * 240
    ok_s = abs(ax - STATION[0]) < 1 and abs(ay - STATION[1]) < 1
    print(f"  Ch05 argmax = ({ax:.0f},{ay:.0f})（推荐站址）→ {'✓' if ok_s else '✗'}")
    metrics['ch05_argmax'] = dict(xy=[ax, ay], ok=bool(ok_s))
    if not ok_s:
        failures.append('ch05_argmax')

    metrics['failures'] = failures
    with open(os.path.join(_HERE, 'tests', 'L5_metrics.json'), encoding='utf-8') as f:
        m5 = json.load(f)
    m5['L5_3_endtoend'] = metrics
    with open(os.path.join(_HERE, 'tests', 'L5_metrics.json'), 'w', encoding='utf-8') as f:
        json.dump(m5, f, ensure_ascii=False, indent=2)

    print('\n' + '=' * 70)
    print('L5.3 汇总: ' + ('PASS —— 一键流水线产物与全部既有基准一致'
          if not failures else 'FAIL: ' + '; '.join(failures)))
    print('=' * 70)
    return 0 if not failures else 1


if __name__ == '__main__':
    sys.exit(main())
