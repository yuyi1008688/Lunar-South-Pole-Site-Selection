# -*- coding: utf-8 -*-
"""
stage2_purepython.py —— 阶段二主分析链（纯 Python 重构版）
四掩膜 → candidate → Dijkstra LCP → B样条平滑（5→3→2 降级）→ 50m 剖面
================================================================
原版主链脚本调用商业 GIS 内嵌组件完成栅格读取、栅格代数、
LCP 寻路与 B 样条平滑（原脚本已随解耦重构删除，历史版本见竞赛提交包）；
本版以纯开源栈等价重构，等价性经 L0-L3 四级验证（见 tests/）：

  | 原算子（商业组件）              | 纯 Python 等价实现                    |
  | GridReader + 8 变换运行时定标   | RasterGrid（rasterio 读标准 GeoTIFF， |
  |                                 | 方向确定，仅需轻量 sanity check）     |
  | expression_math_analyst 的 Con  | numpy 向量化布尔运算                  |
  | 缓冲区→栅格化→Con 补0（dist_ok）| 像元中心到站址平面距离 ≤ 5840 m        |
  | cost_path_line（LCP）           | lcp_dijkstra.dijkstra_accumulation    |
  |                                 | （8 邻域 Dijkstra + 回溯）            |
  | BSPLINE smoothDegree=5          | bspline_smooth（scipy splprep，       |
  |                                 | 5→3→2 自动降级 + 危险区回检）         |
  | 50m 剖面（get_value 最近邻）    | 50m 等距重采样 + 双线性插值取值        |

对标值（历史基准，竞赛期实测留痕）：
  dist_ok=1202 ✓ 终点处=1 ✓ | slope_ok/psr_ok/hazard_ok = 79545/30252/101577 ✓
  candidate=2，终点(40920,10920)处=1 ✓
  LCP（B样条5阶测试档）：总长 4.039 km（顶点17）、直线 3.840 km、
  绕路系数 1.052（<1.5 ✓）、危险区穿越 0 ✓
  剖面：82 点 @50m → path_profile.csv（能量仿真输入）

运行环境：普通 Python 3.10+（numpy + rasterio + scipy）。
前置：data/rasters/ 已由 src/utils/udbx_extract.py 导出。
"""

import csv
import math
import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding='utf-8')
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', 'utils'))
sys.path.insert(0, _HERE)
from raster_grid import RasterGrid, direction_sanity_check
from lcp_dijkstra import (dijkstra_accumulation, path_length_m, trace_backlink)
from bspline_smooth import (bspline_smooth, hazard_crossings, resample,
                            smooth_with_degradation)

# ============================= 配置区 =============================
RASTER_DIR = os.environ.get("LUNAR_RASTER_DIR", os.path.join("data", "rasters"))
OUT_DIR = os.path.join(os.environ.get("LUNAR_OUTPUT_DIR", os.path.join("data", "output")),
                       "stage2")

RASTERS = ["slope_deg", "AVGVISIB_probability", "combined_hazard",
           "cost_surface", "DEM_fused", "vrm_5x5"]

START = (44760.0, 10920.0)   # Ⅰ级站址（AHP 综合得分最高像元）
END = (40920.0, 10920.0)     # PSR 边缘采样点（降级三级口径，stage1 定稿）
FAR_ZERO = (-40000.0, 40000.0)
SLOPE_MAX, PSR_TH, HAZ_MAX = 10.0, 0.001, 1.0
DIST_BUFFER = 5840.0         # 站址缓冲半径（= D_min 3,840 + 2,000 余量）

SMOOTH_DEGREE = 5            # 先用最光滑的 5 阶，穿越则 5→3→2 降级
PROFILE_DS = 50.0            # 剖面采样步长

REF = {"slope_ok": 79545, "psr_ok": 30252, "hazard_ok": 101577,
       "dist_ok": 1202, "candidate": 2,
       "lcp_km": 4.039, "straight_km": 3.840, "ratio": 1.052,
       "profile_n": 82, "vertices": 17}
# ==================================================================


def check(name, got, ref, tol=0):
    ok = abs(got - ref) <= tol
    print(f"    ✔ 对标 {name}: 实测 {got} / 参考 {ref} {'✓' if ok else '← 检查！'}")
    return ok


def build_dist_ok(grid, station_xy, radius=DIST_BUFFER):
    """站址缓冲掩膜（替代原 GUI 工艺：缓冲区→矢量栅格化→Con(IsNull,0,1)）。

    像元中心到站址的平面距离 ≤ radius → 1。历史基准 1202 像元（含研究区
    东边界截断；缓冲圆心为站址 START 而非终点——与原 GUI 工艺一致，且
    恰好复现 1202/candidate=2 两个基准）。
    """
    xs = grid.left + (np.arange(grid.w) + 0.5) * grid.res
    ys = grid.top - (np.arange(grid.h) + 0.5) * grid.res
    X, Y = np.meshgrid(xs, ys)
    return (np.hypot(X - station_xy[0], Y - station_xy[1]) <= radius).astype(np.uint8)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 64)
    print("阶段二主分析链 · 纯 Python（四掩膜→candidate→LCP→B样条→剖面）")
    print("=" * 64)

    grids = {name: RasterGrid.from_file(os.path.join(RASTER_DIR, name + ".tif"))
             for name in RASTERS + ["ice_density_final"]}
    if not direction_sanity_check(grids):
        raise RuntimeError("方向 sanity check 未通过，终止（检查栅格数组方向）")

    # ---------- 动作1：dist_ok ----------
    print("\n[动作1] dist_ok（站址 5840m 缓冲，纯 Python 重建）")
    dist_ok = build_dist_ok(grids["AVGVISIB_probability"], START)
    check("dist_ok", int(dist_ok.sum()), REF["dist_ok"], tol=10)
    rE, cE = grids["AVGVISIB_probability"].rc_at_xy(*END)
    v = int(dist_ok[rE, cE])
    print(f"    终点处 dist_ok = {v} {'✓' if v == 1 else '✗'}")

    # ---------- 检查2：三掩膜 ----------
    print("\n[检查2] slope_ok / psr_ok / hazard_ok（numpy 栅格代数）")
    slope_ok = (grids["slope_deg"].arr < SLOPE_MAX).astype(np.uint8)
    psr_ok = (grids["AVGVISIB_probability"].arr < PSR_TH).astype(np.uint8)
    hazard_ok = (grids["combined_hazard"].arr < HAZ_MAX).astype(np.uint8)
    check("slope_ok", int(slope_ok.sum()), REF["slope_ok"])
    check("psr_ok", int(psr_ok.sum()), REF["psr_ok"])
    check("hazard_ok", int(hazard_ok.sum()), REF["hazard_ok"])

    # ---------- 动作3：candidate ----------
    print("\n[动作3] candidate（四掩膜求交）")
    candidate = slope_ok * psr_ok * hazard_ok * dist_ok
    check("candidate", int(candidate.sum()), REF["candidate"])
    v = int(candidate[rE, cE])
    print(f"    终点处 candidate = {v} {'✓ 里程碑达成' if v == 1 else '✗'}")

    # ---------- 动作5：LCP + 危险区回检 ----------
    print("\n[动作5] LCP（8 邻域 Dijkstra，源=站址）+ 危险区回检")
    cost_pos = grids["cost_surface"].arr + 0.001   # 避免 0 成本，与原流程一致
    rS, cS = grids["AVGVISIB_probability"].rc_at_xy(*START)
    acc, backlink = dijkstra_accumulation(cost_pos, (rS, cS))
    trace = trace_backlink(backlink, (rE, cE), (rS, cS))[::-1]  # 反转为起点→终点
    coords = [grids["AVGVISIB_probability"].xy_at_rc(r, c) for r, c in trace]
    print(f"    原始像素链 {len(coords)} 像元，折线长 {path_length_m(trace) / 1000:.3f} km")

    best, records = smooth_with_degradation(coords, grids["combined_hazard"],
                                            degrees=(5, 3, 2), haz_max=HAZ_MAX)
    final_rec = records[-1]
    hits, _, total = hazard_crossings(best, grids["combined_hazard"],
                                      haz_max=HAZ_MAX)
    coords = best

    # ---------- 路径产出 ----------
    vcsv = os.path.join(OUT_DIR, "optimal_path_vertices.csv")
    with open(vcsv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["id", "x", "y"])
        for i, (x, y) in enumerate(coords):
            w.writerow([i, f"{x:.2f}", f"{y:.2f}"])
    wkt = "LINESTRING (" + ", ".join(f"{x:.2f} {y:.2f}" for x, y in coords) + ")"
    with open(os.path.join(OUT_DIR, "optimal_path_wkt.txt"), "w",
              encoding="utf-8") as f:
        f.write(wkt)
    print(f"    路径顶点 -> {vcsv}（{len(coords)} 点）+ WKT")

    D = math.hypot(END[0] - START[0], END[1] - START[1])
    print("\n    ── 表6-2 指标 ──")
    print(f"    路径总长度   {total / 1000:.3f} km（对标 {REF['lcp_km']}）")
    print(f"    直线距离     {D / 1000:.3f} km")
    print(f"    绕路系数     {total / D:.3f}（对标 {REF['ratio']}，<1.5）")
    print(f"    危险区穿越   {len(hits)}（必须=0 {'✓' if not hits else '✗'}）")

    # ---------- 动作6：50m 剖面（双线性取值） ----------
    print("\n[动作6] 50m 剖面（等距重采样 + 双线性插值）")
    pts, _ = resample(coords, PROFILE_DS)
    csv_path = os.path.join(OUT_DIR, "path_profile.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.writer(fo)
        w.writerow(["id", "s_m", "x", "y", "elev", "slope", "vrm", "visib"])
        for i, (x, y, s_) in enumerate(pts):
            w.writerow([
                i, f"{s_:.1f}", f"{x:.1f}", f"{y:.1f}",
                f"{grids['DEM_fused'].value_at_xy(x, y, method='bilinear'):.2f}",
                f"{grids['slope_deg'].value_at_xy(x, y, method='bilinear'):.3f}",
                f"{grids['vrm_5x5'].value_at_xy(x, y, method='bilinear'):.5f}",
                f"{grids['AVGVISIB_probability'].value_at_xy(x, y, method='bilinear'):.5f}"])
    print(f"    {len(pts)} 点（对标 {REF['profile_n']}）-> {csv_path}")

    print("\n" + "=" * 64)
    print("★ 阶段二纯 Python 主链全部完成！")
    print(f"  成果：optimal_path_vertices.csv / optimal_path_wkt.txt / "
          f"path_profile.csv（{OUT_DIR}）")
    print("  下一步：energy_simulation.py（三态状态机能量递推，算法未改动）")
    print("=" * 64)


if __name__ == "__main__":
    main()
