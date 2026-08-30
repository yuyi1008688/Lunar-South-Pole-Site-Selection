# -*- coding: utf-8 -*-
"""
stage1_select_mining_target.py —— 阶段一：参数定标（纯 Python 重构版）
================================================================
原版为商业 GIS 内嵌组件脚本（需专有运行环境），现已重构为
纯 Python（rasterio 读标准 GeoTIFF + numpy 向量化），统计/判定逻辑、
打印格式与对标值原样保留。

输出（对应《操作指南》第二节）：
  ① D_min：站址到最近 PSR 边界像元的距离
  ② 距离上限 D_max = min(D_min + 2km, 10km)   ← 「距离重分类」参数
  ③ F2 阈值-候选数分布表 + 建议值             ← 「F2 重分类」参数
  并写 stage1_params.txt 留档。

实测参考（历史基准）：D_min=3,840 m，D_max=5,840 m；
D_max 内无 F2>0 候选 → 采用降级三级（PSR 边缘采样）口径。

运行环境：普通 Python 3.10+（numpy + rasterio）。
前置：data/rasters/ 目录已由 src/utils/udbx_extract.py 导出（见 data/README.md）。
"""

import math
import os
import sys
from collections import deque

import numpy as np

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'utils'))
from raster_grid import RasterGrid

# ============================= 配置区 =============================
RASTER_DIR = os.environ.get("LUNAR_RASTER_DIR", os.path.join("data", "rasters"))
OUT_TXT = os.path.join(os.environ.get("LUNAR_OUTPUT_DIR", os.path.join("data", "output")),
                       "stage1_params.txt")

RASTERS = {                        # 栅格文件名（相对 RASTER_DIR）
    "f2":     "ice_density_final.tif",     # F2 水冰丰度因子（Wang KDE 最终交付版）
    "slope":  "slope_deg.tif",
    "hazard": "combined_hazard.tif",
    "visib":  "AVGVISIB_probability.tif",
}
STATION_X, STATION_Y = 44760.0, 10920.0
SLOPE_MAX, HAZ_MAX, PSR_TH = 10.0, 1.0, 0.001
DIST_MARGIN, DIST_CAP = 2000.0, 10000.0
MIN_AREA_KM2 = 0.1
SCAN_TH = [0.1, 0.2, 0.3, 0.5]
# ==================================================================


def label8(mask):
    """8 邻域连通分量标记（BFS），返回 (labels, sizes)。"""
    labels = np.zeros(mask.shape, dtype=np.int32)
    sizes = []
    rows, cols = mask.shape
    cur = 0
    for r0, c0 in np.argwhere(mask):
        if labels[r0, c0]:
            continue
        cur += 1
        n = 0
        q = deque([(r0, c0)])
        labels[r0, c0] = cur
        while q:
            r, c = q.popleft()
            n += 1
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    rr, cc = r + dr, c + dc
                    if 0 <= rr < rows and 0 <= cc < cols \
                       and mask[rr, cc] and not labels[rr, cc]:
                        labels[rr, cc] = cur
                        q.append((rr, cc))
        sizes.append(n)
    return labels, sizes


def main():
    print("=" * 64)
    print("阶段一：参数定标（纯 Python / rasterio）")
    print("=" * 64)

    grids = {k: RasterGrid.from_file(os.path.join(RASTER_DIR, v))
             for k, v in RASTERS.items()}
    f2 = grids["f2"]
    res = f2.res
    n_r, n_c = f2.arr.shape

    xs = f2.left + (np.arange(n_c) + 0.5) * res
    ys = f2.top - (np.arange(n_r) + 0.5) * res
    X, Y = np.meshgrid(xs, ys)
    dist = np.hypot(X - STATION_X, Y - STATION_Y)
    print(f"网格：{n_r}x{n_c} @ {res:g} m，站址 ({STATION_X:g},{STATION_Y:g})")

    # ---------- ① D_min ----------
    psr = np.nan_to_num(grids["visib"].arr, nan=9e9) < PSR_TH
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
    k = int(np.argmin(np.where(bnd, dist, np.inf)))
    br, bc = np.unravel_index(k, psr.shape)
    d_min = float(dist[br, bc])
    print("\n① D_min（站址→最近 PSR 边界像元）")
    print(f"  PSR 像元 {int(psr.sum())}（{psr.mean():.2%}），边界 {int(bnd.sum())} 个")
    print(f"  D_min = {d_min:,.1f} m，最近边界像元 ({X[br, bc]:g}, {Y[br, bc]:g})"
          f"（历史基准 3,840 m）")

    # ---------- ② D_max ----------
    d_max = min(d_min + DIST_MARGIN, DIST_CAP)
    print(f"\n② 距离上限 D_max = min(D_min+{DIST_MARGIN:g}, {DIST_CAP:g}) = "
          f"{d_max:,.1f} m（历史基准 5,840 m）")

    # ---------- ③ F2 阈值分布 ----------
    slope = grids["slope"].arr
    hazard = grids["hazard"].arr
    hard = (np.nan_to_num(slope, nan=99) < SLOPE_MAX) \
         & (np.nan_to_num(hazard, nan=99) < HAZ_MAX) & (dist <= d_max)
    f2v = np.nan_to_num(f2.arr, nan=0.0)
    pos = f2v[f2v > 0]
    print(f"\n③ F2 阈值分布（硬约束底盘 {int(hard.sum())} 像元）")
    ths = sorted(set(SCAN_TH + ([round(float(np.percentile(pos, p)), 4)
                                 for p in (50, 75, 90, 95)] if pos.size else [])))
    min_px = max(2, int(math.ceil(MIN_AREA_KM2 * 1e6 / (res * res))))
    best = None
    rows_out = []
    print(f"  {'阈值':>8} | {'候选':>5} | {'连通块':>5} | {'合格块':>5} | 评价")
    for th in ths:
        cand = hard & (f2v >= th)
        cnt = int(cand.sum())
        if cnt:
            _, sizes = label8(cand)
            nbig = sum(1 for s_ in sizes if s_ >= min_px)
        else:
            sizes, nbig = [], 0
        v = "✓ 合理" if (nbig and cnt <= 200) else ("候选偏多" if nbig else "✗ 无合格块")
        if v == "✓ 合理" and (best is None or th > best):
            best = th
        rows_out.append((th, cnt, len(sizes), nbig, v))
        print(f"  {th:8.4f} | {cnt:5d} | {len(sizes):5d} | {nbig:5d} | {v}")
    if best is not None:
        print(f"\n  建议：F2 阈值取 {best:g}")
    else:
        print("\n  建议：无合格连通块 → 按指南降级三级：改用 AVGVISIB<0.001（PSR 边缘"
              "采样），终点=距站址最近的合规 PSR 边界像元（历史实测 (40920,10920)）。")

    print("\n" + "=" * 64)
    print("★ 抄进配置的参数")
    print(f"  D_min = {d_min:,.1f} m ｜ D_max = {d_max:,.1f} m（距离重分类『≤此值→1』）")
    if best is not None:
        print(f"  F2 阈值 = {best:g}（F2 重分类『≥此值→1』）")
    else:
        print("  掩膜口径 = AVGVISIB < 0.001 → 1（降级三级，替代 f2_ok）")
    print("=" * 64)

    try:
        os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
        with open(OUT_TXT, "w", encoding="utf-8") as f:
            f.write(f"station=({STATION_X},{STATION_Y})\nD_min_m={d_min:.1f}\n"
                    f"D_max_m={d_max:.1f}\n")
            f.write("threshold,candidates,components,big,verdict\n")
            for r_ in rows_out:
                f.write(",".join(str(x) for x in r_) + "\n")
        print("参数已留档：", OUT_TXT)
    except OSError:
        pass


if __name__ == "__main__":
    main()
