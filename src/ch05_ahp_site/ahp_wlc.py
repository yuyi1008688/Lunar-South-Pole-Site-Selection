# -*- coding: utf-8 -*-
"""
ahp_wlc.py —— Ch05 AHP-WLC 综合选址成图（纯 Python 重构版）
================================================================
当年本章核心步（WLC 加权叠加 → 硬约束掩膜 → Jenks 五级分级 → Ⅰ级点导出）
在商业 GIS 栅格计算器/重分类工具中完成；本模块以 numpy 向量化等价复刻，
与存档产物逐像元对标（tests/L5_verify_ch02_ch05.py，L5.2）。

算法口径（全部经存档反解核实，算法零改动）：
  1. WLC 加权叠加（AHP 权重固定，CR=0.54%，不许改）：
       S = 0.331·F1 + 0.258·F2 + 0.191·F3 + 0.126·F4 + 0.094·F5
     NoData 传染修复：每层先 where(isfinite & ≠NoData, ·, 0)（等价当年
     Con(IsNull(f), 0, f)）。
  2. F2 有效口径（存档反解确证，见验证报告 §L5.2）：最终存档版 suitability
     在约束域内 F2 贡献**严格为 0**（max|sf−rest|=2.98e-8，纯 float32 量化
     噪声）——因为最终迭代的 F2 栅格在 PSR 外为 NoData，经 Con(IsNull,0,·)
     归零（与 Ch03"PSR 外强制赋 0"一致）。故本模块取
       F2_eff = F2 × (p < 0.001)   （PSR 掩膜口径）
     约束域内 F2_eff≡0，与存档逐像元一致。对照：unit5 早期迭代直接用 scipy
     F2（约束域边缘渗漏 ≤0.036），两版 sf 差 0.009，属迭代差异、已登记。
  3. 选址硬约束掩膜（注意与 Ch06 选采矿目标掩膜 10/0.001/1.0 区分）：
       M = (坡度 ≤ 20°) ∧ (光照概率 p ≥ 0.001) ∧ (距危险区 ≥ 240 m)
     第一性复现与存档 constraint_mask_v32 **80,749 像元全等**。光照下限
     0.001 为 subPSR 边界口径（GPA 模型文档"与 F1 下界保持一致"；0.20 是
     Ch04 的 hard_constraint_v32 口径，两套掩膜并存，勿混）。
     掩膜外 suitability = 0（与存档一致：存档无 NoData 像元）。
  4. Jenks 五级分级：断点取**存档分级栅格反推的全精度重分类表**
       [0.34634405, 0.47284603, 0.54648739, 0.57978046]（最高级至 0.65227）
     （GPA 模型文档记载的 0.332/0.467/0.542/0.572 对应 unit5 早期迭代，
      与最终存档分级栅格边界不同——迭代差异，已在验证报告登记。）
     分级时容差 tol=1e-7 用于吸收存档 float32 量化（非阈值更改）。
  5. Ⅰ级点导出：分级值 5（得分最高区间 [0.5798, 0.6523]）的全部像元，
     坐标公式按 P0 修复口径 x = tf.c + (col+0.5)·tf.a、y = tf.f + (row+0.5)·tf.e
     （严禁 tf.a/tf.c 混用——当年 P0 bug 曾让坐标跑到 ±17,671 km）。

历史基准（L5.2 对标）：suitability 与存档 max|diff|≈1e-7（存档 float32 量化
下限）；分级一致率 100%；argmax = 推荐站址 (44760, 10920)；约束掩膜 80,749。
注：历史文档记载的"Ⅰ级 5,214 点"为 2026-08-03 P0 修复时代的早期迭代计数，
其底层数据未归档，无法从最终存档复现（见验证报告 §L5.2 根因登记）。
"""

import csv
import os
import sys

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin

sys.stdout.reconfigure(encoding='utf-8')
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', 'utils'))
from iron_grid import assert_iron_grid

RASTER_DIR = os.environ.get('LUNAR_RASTER_DIR', os.path.join('data', 'rasters'))
F1_DIR = os.environ.get('LUNAR_CH02_DIR', RASTER_DIR)   # 管线模式下传 ch02 产物目录
OUT_DIR = os.path.join(os.environ.get('LUNAR_OUTPUT_DIR', os.path.join('data', 'output')),
                       'ch05')

CRS_PROJ4 = ('+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +x_0=0 +y_0=0 '
             '+R=1737400 +units=m +no_defs')
NODATA = -9999.0

# AHP 权重（F1..F5，CR=0.54%，固定不许改）
WEIGHTS = (0.331, 0.258, 0.191, 0.126, 0.094)
# 选址硬约束（GPA 模型文档口径：与 F1 下界一致的 subPSR 边界）
SLOPE_MAX_HARD = 20.0
VISIB_MIN_HARD = 0.001
HAZ_DIST_MIN = 240.0
# Jenks 重分类表（存档分级栅格反推的全精度断点；分级值 5 = Ⅰ级）
JENKS_BREAKS = (0.3463440537452698, 0.47284603118896484,
                0.5464873909950256, 0.5797804594039917)
JENKS_TOL = 1e-7   # 吸收存档 float32 量化（非阈值更改）

STATION = (44760.0, 10920.0)


def _read(name, directory=RASTER_DIR):
    with rasterio.open(os.path.join(directory, name + '.tif')) as src:
        arr = src.read(1).astype(np.float64)
        if src.nodata is not None:
            arr[arr == src.nodata] = np.nan
    return arr


def _nodata_zero(a):
    """NoData 传染修复（等价当年 Con(IsNull(f), 0, f)）。"""
    return np.where(np.isfinite(a), a, 0.0)


def hard_constraint_mask(slope, visib, hazard_distance):
    """选址硬约束：坡度≤20° ∧ 光照概率≥0.001 ∧ 距危险区≥240 m。

    第一性复现与存档 constraint_mask_v32 80,749 像元全等（L5.2 核对）。
    """
    return ((slope <= SLOPE_MAX_HARD) & (visib >= VISIB_MIN_HARD)
            & (hazard_distance >= HAZ_DIST_MIN))


def ch04_hard_constraint_v32(slope, visib, hazard_distance,
                             visib_min=0.20, n_ref=42344):
    """Ch04 口径的硬约束变体（光照下限 0.20，即 hard_constraint_v32）。

    可由 (slope≤20)∧(p≥0.20)∧(距险≥240) 第一性复现存档 42,344 像元全等
    （验证报告 §L5.2）。与选址掩膜（0.001 口径）并存、各司其职。
    """
    return ((slope <= SLOPE_MAX_HARD) & (visib >= visib_min)
            & (hazard_distance >= HAZ_DIST_MIN))


def f2_effective(f2, visib):
    """F2 有效口径：PSR 掩膜（存档反解确证，见模块 docstring §2）。"""
    return f2 * (visib < 0.001)


def wlc(f1, f2, f3, f4, f5):
    """WLC 加权叠加（NoData→0 后线性加权）。"""
    w = WEIGHTS
    return (w[0] * _nodata_zero(f1) + w[1] * _nodata_zero(f2)
            + w[2] * _nodata_zero(f3) + w[3] * _nodata_zero(f4)
            + w[4] * _nodata_zero(f5))


def jenks_classes(scores, tol=JENKS_TOL):
    """按存档重分类表分五级（1 最差 … 5=Ⅰ级最优）；tol 吸收 float32 量化。"""
    b1, b2, b3, b4 = JENKS_BREAKS
    return np.where(scores <= b1 + tol, 1.0,
           np.where(scores <= b2 + tol, 2.0,
           np.where(scores <= b3 + tol, 3.0,
           np.where(scores <= b4 + tol, 4.0, 5.0))))


def grade_i_points(scores, classes, transform):
    """Ⅰ级点导出（P0 修复坐标公式）：返回 [(x, y, score), ...]。"""
    rr, cc = np.where(classes == 5)
    xs = transform.c + (cc + 0.5) * transform.a
    ys = transform.f + (rr + 0.5) * transform.e
    return list(zip(xs.tolist(), ys.tolist(), scores[rr, cc].tolist()))


def run(f1_dir=F1_DIR, out_dir=OUT_DIR):
    """主入口：WLC → 硬约束 → 分级 → Ⅰ级点，铁基准校验后写盘。"""
    f1 = _read('F1_illumination_norm', f1_dir)
    f2 = _read('F2_wang_kde_final_1')
    f3 = _read('F3_fos_safety')
    f4 = _read('F4_ecsa_sync')
    f5 = _read('F5_continuous_distance')
    slope = _read('slope_deg')
    visib = _read('AVGVISIB_probability')
    hz_d = _read('hazard_distance')

    # 硬约束掩膜（第一性复现，与存档全等性在 L5.2 核对）
    mask = hard_constraint_mask(slope, visib, hz_d)
    print(f'[Ch05] 硬约束掩膜: {int(mask.sum())} 像元'
          f'（存档 constraint_mask_v32 基准 80,749）')

    # WLC（F2 取 PSR 掩膜口径）× 硬约束
    f2_eff = f2_effective(f2, visib)
    s_raw = wlc(f1, f2_eff, f3, f4, f5)
    suitability = s_raw * mask          # 掩膜外=0（与存档一致，无 NoData）
    classes = jenks_classes(suitability)

    # 站址核对（argmax 必须 = 推荐站址）
    r_pk, c_pk = np.unravel_index(np.argmax(suitability), suitability.shape)
    tx = from_origin(-46080, 46080, 240, 240)
    sx = tx.c + (c_pk + 0.5) * tx.a
    sy = tx.f + (r_pk + 0.5) * tx.e
    print(f'[Ch05] 适宜性最大值 {suitability.max():.4f} @ ({sx:.0f},{sy:.0f})'
          f'（推荐站址 (44760,10920)）{"✓" if abs(sx - STATION[0]) < 1 and abs(sy - STATION[1]) < 1 else "✗"}')

    # Ⅰ级点
    pts = grade_i_points(suitability, classes, tx)
    print(f'[Ch05] Ⅰ级点（分级值 5）: {len(pts)} 个'
          f'（历史文档 5,214 为 P0 时代早期迭代计数，见验证报告 §L5.2）')

    # 铁基准校验 + 写盘
    os.makedirs(out_dir, exist_ok=True)
    profile = dict(driver='GTiff', height=384, width=384, count=1, dtype='float32',
                   crs=CRS.from_proj4(CRS_PROJ4), transform=tx,
                   nodata=NODATA, compress='deflate')
    outs = {}
    for name, arr in [('suitability_final', suitability),
                      ('suitability_classes', classes),
                      ('constraint_mask', mask.astype(np.float64))]:
        assert_iron_grid(f'ch05/{name}', arr, tx, profile['crs'],
                         expect_full_valid=True)
        path = os.path.join(out_dir, name + '.tif')
        with rasterio.open(path, 'w', **profile) as dst:
            dst.write(arr.astype('float32'), 1)
        outs[name] = path

    csv_path = os.path.join(out_dir, 'grade_I_points.csv')
    with open(csv_path, 'w', newline='', encoding='utf-8-sig') as f:
        wcsv = csv.writer(f)
        wcsv.writerow(['x', 'y', 'score'])
        for x, y, s in pts:
            wcsv.writerow([f'{x:.2f}', f'{y:.2f}', f'{s:.6f}'])
    print(f'[Ch05] 已输出: {outs["suitability_final"]}\n'
          f'        {outs["suitability_classes"]}\n'
          f'        {outs["constraint_mask"]}\n        {csv_path}')
    return dict(suitability_tif=outs['suitability_final'],
                classes_tif=outs['suitability_classes'],
                mask_tif=outs['constraint_mask'],
                grade_csv=csv_path, n_grade_i=len(pts),
                argmax_xy=[sx, sy], mask_pixels=int(mask.sum()))


if __name__ == '__main__':
    run()
