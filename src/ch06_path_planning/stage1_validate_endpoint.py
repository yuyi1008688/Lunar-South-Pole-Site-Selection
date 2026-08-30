# -*- coding: utf-8 -*-
"""
stage1_validate_endpoint.py —— 阶段一收尾：终点六项合规验证（纯 Python 重构版）
================================================================
原版为商业 GIS 内嵌组件脚本，现已重构为纯 Python（rasterio 读标准 GeoTIFF），
判定逻辑与打印格式原样保留。

验证终点（降级三级 / PSR 边缘采样口径）六项：
  ① AVGVISIB < 0.001（PSR 内）      ② 位于 PSR 边界（8 邻域含非 PSR）
  ③ 坡度 < 10°                      ④ hazard < 1.0
  ⑤ 距站址 ≤ D_max                  ⑥ cost_surface < 1.0（LCP 可通行）
不满足时自动在 PSR 边界上搜索『满足全部约束且距站址最近』的替代终点。

实测参考（历史基准，终点 (40920,10920)）：
  visib=0.000000 ✓ 边界 ✓ 坡度2.33° ✓ hazard=0 ✓ 距离3,840≤5,840 ✓ cost=0.0101 ✓

运行环境：普通 Python 3.10+（numpy + rasterio）。
"""

import math
import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'utils'))
from raster_grid import RasterGrid

# ============================= 配置区 =============================
RASTER_DIR = os.environ.get("LUNAR_RASTER_DIR", os.path.join("data", "rasters"))
RASTERS = {"slope": "slope_deg.tif", "hazard": "combined_hazard.tif",
           "visib": "AVGVISIB_probability.tif", "cost": "cost_surface.tif"}
START = (44760.0, 10920.0)
END = (40920.0, 10920.0)
PSR_TH, SLOPE_MAX, HAZ_MAX = 0.001, 10.0, 1.0
DIST_MARGIN, DIST_CAP = 2000.0, 10000.0
# ==================================================================


def main():
    print("=" * 64)
    print("终点六项合规验证（纯 Python / rasterio）")
    print("=" * 64)

    grids = {k: RasterGrid.from_file(os.path.join(RASTER_DIR, v)).arr
             for k, v in RASTERS.items()}
    vis = grids["visib"]
    n_r, n_c = vis.shape
    res = 240.0
    xmin, ymax = -46080.0, 46080.0

    xs = xmin + (np.arange(n_c) + 0.5) * res
    ys = ymax - (np.arange(n_r) + 0.5) * res
    X, Y = np.meshgrid(xs, ys)
    dist = np.hypot(X - START[0], Y - START[1])

    psr = np.nan_to_num(vis, nan=9e9) < PSR_TH
    inter = psr.copy()
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == dc == 0:
                continue
            sh = np.zeros_like(psr)
            rs = slice(max(dr, 0), n_r + min(dr, 0))
            rd = slice(max(-dr, 0), n_r + min(-dr, 0))
            cs = slice(max(dc, 0), n_c + min(dc, 0))
            cd = slice(max(-dc, 0), n_c + min(-dc, 0))
            sh[rd, cd] = psr[rs, cs]
            inter &= sh
    bnd = psr & ~inter
    d_min = float(dist[bnd].min())
    d_max = min(d_min + DIST_MARGIN, DIST_CAP)

    er, ec = int((ymax - END[1]) / res), int((END[0] - xmin) / res)
    checks = [
        (f"AVGVISIB = {vis[er, ec]:.6f} < {PSR_TH}（PSR 内）",
         vis[er, ec] < PSR_TH),
        ("位于 PSR 边界（8 邻域含非 PSR）", bool(bnd[er, ec])),
        (f"坡度 = {grids['slope'][er, ec]:.2f}° < {SLOPE_MAX}°",
         grids['slope'][er, ec] < SLOPE_MAX),
        (f"hazard = {grids['hazard'][er, ec]:.3f} < {HAZ_MAX}",
         grids['hazard'][er, ec] < HAZ_MAX),
        (f"距站址 = {math.hypot(END[0] - START[0], END[1] - START[1]):,.1f} m "
         f"≤ D_max {d_max:,.1f} m",
         math.hypot(END[0] - START[0], END[1] - START[1]) <= d_max),
        (f"cost_surface = {grids['cost'][er, ec]:.4f} < 1.0（LCP 可通行）",
         grids['cost'][er, ec] < 1.0),
    ]
    print(f"\n终点 {END} 合规检查（D_min={d_min:,.0f}m, D_max={d_max:,.0f}m）")
    ok = True
    for msg, passed in checks:
        print(f"  [{'✓' if passed else '✗'}] {msg}")
        ok &= passed

    if ok:
        print("\n★ 全部通过：终点定稿，可进入 stage2 主分析链。")
    else:
        okm = bnd & (np.nan_to_num(grids['slope'], nan=99) < SLOPE_MAX) \
                  & (np.nan_to_num(grids['hazard'], nan=99) < HAZ_MAX) \
                  & (np.nan_to_num(grids['cost'], nan=1) < 1.0) & (dist <= d_max)
        if okm.any():
            k = int(np.argmin(np.where(okm, dist, np.inf)))
            r, c = np.unravel_index(k, psr.shape)
            print(f"\n替代终点：({X[r, c]:g}, {Y[r, c]:g})，距站址 {dist[r, c]:,.1f} m"
                  "\n→ 更新终点后重跑本脚本确认。")
        else:
            print("\nD_max 内无合规边界像元 → 放宽 D_max（如 10000）再试。")
    print("=" * 64)


if __name__ == "__main__":
    main()
