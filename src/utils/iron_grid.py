# -*- coding: utf-8 -*-
"""
iron_grid.py —— 统一空间铁基准断言工具
================================================================
所有阶段产物写入前自动校验（任务要求 §2）：
  - 尺寸 384×384、像元 240 m、范围 ±46,080 m
  - CRS 为月球南极立体投影（proj4 含 R=1737400）
  - NoData 约定 -9999；统计有限值占比
  - 打印 min/max/mean 留痕

用法：
  from utils.iron_grid import assert_iron_grid
  info = assert_iron_grid('ch02/C2_illumination_class', arr, transform, crs)
  # 校验失败抛 IronGridViolation；通过返回统计 dict
"""

import numpy as np

GRID_N = 384
RES = 240.0
LIMIT = 46080.0
NODATA = -9999.0
CRS_KEY = '1737400'   # 月球半径（m），proj4/WKT 中必含


class IronGridViolation(AssertionError):
    """铁基准校验失败。"""


def assert_iron_grid(name, arr, transform, crs, nodata=NODATA, expect_full_valid=False):
    """校验 (arr, transform, crs) 满足铁基准，打印统计并返回摘要。

    参数
    ----
    name : 产物名（日志用）
    arr : (H, W) 数组（float）
    transform : rasterio Affine
    crs : rasterio CRS
    nodata : 本产物的 NoData 值（默认 -9999）
    expect_full_valid : True 时要求无 NoData 像元（如分级/掩膜产品）
    """
    problems = []
    arr = np.asarray(arr, dtype=np.float64)
    if arr.shape != (GRID_N, GRID_N):
        problems.append(f"尺寸 {arr.shape} ≠ (384, 384)")
    if abs(abs(transform.a) - RES) > 1e-9 or abs(abs(transform.e) - RES) > 1e-9:
        problems.append(f"像元尺寸 ({transform.a},{transform.e}) ≠ 240")
    b = (transform.c, transform.f + arr.shape[0] * transform.e,
         transform.c + arr.shape[1] * transform.a, transform.f)
    for got, want, label in [(b[0], -LIMIT, 'left'), (b[3], LIMIT, 'top'),
                             (b[2], LIMIT, 'right'), (b[1], -LIMIT, 'bottom')]:
        if abs(got - want) > 1e-6:
            problems.append(f"bounds.{label}={got} ≠ {want}")
    crs_text = str(crs) if crs is not None else ''
    if CRS_KEY not in crs_text:
        problems.append(f"CRS 不含月球半径 {CRS_KEY}: {crs_text[:60]}")

    valid = (arr != nodata) & np.isfinite(arr)
    n_valid = int(valid.sum())
    finite_frac = n_valid / arr.size
    if expect_full_valid and n_valid != arr.size:
        problems.append(f"要求全有效，但 NoData/非有限像元 {arr.size - n_valid} 个")
    if (arr[valid] < -1e9).any() if n_valid else False:
        problems.append("有效值中存在疑似 NoData 残留（< -1e9）")

    summary = dict(name=name, shape=list(arr.shape), valid=n_valid,
                   finite_frac=finite_frac, nodata=nodata)
    if n_valid:
        v = arr[valid]
        summary.update(min=float(v.min()), max=float(v.max()), mean=float(v.mean()))
        print(f"  [铁基准] {name}: 384×384@240m ✓ | 有效 {n_valid}/{arr.size}"
              f" ({finite_frac:.1%}) | min={v.min():.6g} max={v.max():.6g} "
              f"mean={v.mean():.6g}")
    else:
        print(f"  [铁基准] {name}: 384×384@240m ✓ | 有效 0（全 NoData）")
    if problems:
        msg = f"{name} 铁基准校验失败: " + '; '.join(problems)
        print(f"  [铁基准] ✗ {msg}")
        raise IronGridViolation(msg)
    return summary
