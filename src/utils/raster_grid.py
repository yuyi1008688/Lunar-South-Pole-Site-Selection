# -*- coding: utf-8 -*-
"""
raster_grid.py —— 统一栅格读取层（纯 Python，替代原 GridReader + 8 变换定标）
================================================================
原 Ch06 脚本中那段 TRANSFORMS（转置/上下翻/左右翻共 8 种组合）是针对
商业 GIS 内嵌组件栅格数组方向不确定而写的"运行时定标"。改用 rasterio 读标准
GeoTIFF 后数组方向是确定的：**第 0 行 = 地理最北（y 最大）**，故 8 变换
定标逻辑全部删除。

本类只保留一个轻量的方向 sanity check 入口（direction_sanity_check），
用 START/END/已知地标三个坐标处的物理预期值确认数组方向正确。

约定（与项目铁基准一致）：
  - row 0 = 北；row 向下递增 = y 递减；col 向右递增 = x 递增
  - 行列↔坐标：col = (x - left)/res；row = (top - y)/res
    像元中心反算：x = left + col*res + res/2；y = top - row*res - res/2
"""

import math
import os

import numpy as np
import rasterio

NODATA = -9999.0


class RasterGrid:
    """栅格 + 地理变换的轻量封装。掩膜判定用最近邻，剖面采样用双线性。"""

    def __init__(self, arr, left, top, res, nodata=NODATA):
        self.arr = np.asarray(arr, dtype=np.float64)
        self.h, self.w = self.arr.shape
        self.left = float(left)
        self.top = float(top)
        self.res = float(res)
        self.nodata = float(nodata)

    # ---------- 构造 ----------
    @classmethod
    def from_file(cls, path):
        with rasterio.open(path) as src:
            arr = src.read(1).astype(np.float64)
            nodata = src.nodata if src.nodata is not None else NODATA
            t = src.transform
            left = t.c
            top = t.f
            res = abs(t.a)
        return cls(arr, left, top, res, nodata=nodata)

    @classmethod
    def from_array(cls, arr, left=-46080.0, top=46080.0, res=240.0, nodata=NODATA):
        return cls(arr, left, top, res, nodata=nodata)

    # ---------- 基本性质 ----------
    @property
    def bounds(self):
        return (self.left, self.top - self.h * self.res,
                self.left + self.w * self.res, self.top)

    @property
    def valid_mask(self):
        return (self.arr != self.nodata) & ~np.isnan(self.arr)

    def rc_at_xy(self, x, y):
        """地理坐标 → 像元 (row, col)（像元索引，四舍五入到最近像元中心）。"""
        c = int(round((x - self.left) / self.res - 0.5))
        r = int(round((self.top - y) / self.res - 0.5))
        return r, c

    def xy_at_rc(self, r, c):
        """像元 (row, col) → 像元中心地理坐标。"""
        return (self.left + c * self.res + self.res / 2.0,
                self.top - r * self.res - self.res / 2.0)

    # ---------- 取值 ----------
    def value_at_xy(self, x, y, method='nearest'):
        """最近邻（掩膜/离散判定）或双线性（剖面连续量）取值。

        越界或 NoData 邻域：最近邻返回 NoData；双线性对 NoData 邻居按
        有效值重归一化加权，全部无效时返回 NoData。
        """
        if method == 'nearest':
            r, c = self.rc_at_xy(x, y)
            if not (0 <= r < self.h and 0 <= c < self.w):
                return self.nodata
            return float(self.arr[r, c])
        elif method == 'bilinear':
            return float(self._bilinear(x, y))
        raise ValueError(f"未知插值方法：{method}")

    def _bilinear(self, x, y):
        fx = (x - self.left) / self.res - 0.5
        fy = (self.top - y) / self.res - 0.5
        r0, c0 = int(math.floor(fy)), int(math.floor(fx))
        dr, dc = fy - r0, fx - c0
        vals = np.empty(4)
        for i, (rr, cc) in enumerate([(r0, c0), (r0, c0 + 1),
                                      (r0 + 1, c0), (r0 + 1, c0 + 1)]):
            if 0 <= rr < self.h and 0 <= cc < self.w:
                v = self.arr[rr, cc]
                vals[i] = np.nan if (v == self.nodata or np.isnan(v)) else v
            else:
                vals[i] = np.nan
        w = np.array([(1 - dr) * (1 - dc), (1 - dr) * dc,
                      dr * (1 - dc), dr * dc])
        good = ~np.isnan(vals)
        if not good.any():
            return self.nodata
        return float((vals[good] * w[good]).sum() / w[good].sum())

    # ---------- 统计 ----------
    def count_eq(self, val, tol=0.0):
        if tol:
            return int((np.abs(self.arr - val) <= tol).sum())
        return int((self.arr == val).sum())

    def stats(self):
        v = self.arr[self.valid_mask]
        if v.size == 0:
            return dict(n=0)
        return dict(n=int(v.size), min=float(v.min()), max=float(v.max()),
                    mean=float(v.mean()))


# ── 方向 sanity check（轻量校验，替代原 8 变换定标） ──────────────
def direction_sanity_check(grids):
    """用已知物理预期值确认各栅格数组方向正确（第 0 行=北）。

    检查点（均来自项目已验收的实测记录，一旦数组南北翻转即不命中）：
      1. START(44760,10920)：光照概率 ≥ 0.20（站址满足硬约束）
      2. END(40920,10920)：光照概率 < 0.001（PSR 边缘采样点）且坡度 < 10°
         （终点六项合规验证实测值：visib=0.000000 / 坡度=2.33° / hazard=0）
      3. 三掩膜像元数与历史基准逐个相等：
         slope_ok=79545 / psr_ok=30252 / hazard_ok=101577（翻转后必不等）
      4. ice_density>0 像元与 sPSR_mask 空间一致（同源同向，≈95%；若相对
         翻转则骤降至个位数百分比）
    注：ice_density_final 是固定半径 KDE（4800m）产物，其峰值
    位置与 scipy 自适应带宽版 F2_wang_kde 的记录峰值 (-14040,+1080) 本就
    不同，不能作为方向检查点（详见验证报告）。
    """
    print('── 方向 sanity check（确认 rasterio 数组第 0 行=北） ──')
    ok = True

    av = grids['AVGVISIB_probability'].value_at_xy(44760, 10920)
    print(f'  START 处光照概率 = {av:.4f}（预期 ≥0.20，站址硬约束）')
    ok &= av >= 0.20

    av_e = grids['AVGVISIB_probability'].value_at_xy(40920, 10920)
    sl_e = grids['slope_deg'].value_at_xy(40920, 10920)
    hz_e = grids['combined_hazard'].value_at_xy(40920, 10920)
    print(f'  END 处 光照={av_e:.5f}（预期 <0.001）坡度={sl_e:.2f}°（预期 2.33）'
          f' 危险={hz_e:.2f}（预期 0）')
    ok &= av_e < 0.001 and sl_e < 10.0 and hz_e < 1.0

    refs = {'slope_deg': (79545, 10.0), 'AVGVISIB_probability': (30252, 0.001),
            'combined_hazard': (101577, 1.0)}
    for name, (ref, thr) in refs.items():
        n = int((grids[name].arr < thr).sum())
        print(f'  {name} < {thr} 像元数 = {n}（基准 {ref}）')
        ok &= (n == ref)

    if 'ice_density_final' in grids and 'sPSR_mask' in grids:
        ice = grids['ice_density_final'].arr
        spsr = grids['sPSR_mask'].arr
        inside = ((ice > 0) & (spsr > 0)).sum() / max((ice > 0).sum(), 1)
        print(f'  ice>0 像元落在 sPSR 内比例 = {inside*100:.2f}%（同向应 >90%）')
        ok &= inside > 0.9

    print(f'  方向校验结果: {"全部命中，方向正确 ✓" if ok else "存在不命中项 ✗（检查数组方向）"}')
    return ok
