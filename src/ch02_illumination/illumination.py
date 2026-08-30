# -*- coding: utf-8 -*-
"""
illumination.py —— Ch02 光照时空分布分析与 PSR 界定（纯 Python 重构版）
================================================================
当年本章全部处理（AVGVISIB 概率解码、双 PSR 分类、五级光照分类、F1 归一化）
在商业 GIS 栅格计算器中完成，仅有 GPA 表达式留档；本模块以 numpy 向量化
等价复刻，逐像元对标存档产物（见 tests/L5_verify_ch02_ch05.py，L5.1）。

算法口径（与 GPA 模型节点 1~7 的 Con 嵌套表达式逐项对齐，算法零改动）：
  1. 概率解码：原始 int16 编码 = 概率×25000。先做量纲体检（max>1 → ÷25000；
     max≤1 说明数据源内存的已是概率，勿重复解码）。
  2. 双 PSR 分类（物理依据：水冰升华速率对温度的非线性，见 methodology §Ch2）：
       sPSR   : p < 1e-6      （严格永久阴影，水冰核心保存区）
       subPSR : 1e-6 ≤ p < 1e-3（次永久阴影，亚稳定）
       PSR    : p < 1e-3      （联合掩膜）
  3. 五级光照分类（整数 1–5，与 GPA Con 嵌套完全等价）：
       class 1: p ≤ 1e-6
       class 2: 1e-6 < p ≤ 1e-3
       class 3: 1e-3 < p ≤ 0.128
       class 4: 0.128 < p ≤ 0.264
       class 5: p > 0.264   （连续光照区，P75）
  4. F1 分段线性归一化（Ch05 五因子之一）：
       F1 = 0                    当 p < 0.001
       F1 = (p − 0.001)/0.263    当 0.001 ≤ p ≤ 0.264
       F1 = 1                    当 p > 0.264
     注意 0.20 是 Ch05 选址硬约束阈值，与 F1 归一化各司其职，不在本公式内。

历史基准（存档实测，L5.1 对标）：sPSR 19.5% / subPSR 1.0% / PSR 20.5% /
连续光照 25.0%；F1 与存档 max|diff|≈3e-8；五级分类一致率 100%。
"""

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
OUT_DIR = os.environ.get('LUNAR_OUTPUT_DIR', os.path.join('data', 'output'))
INPUT_TIF = os.path.join(RASTER_DIR, 'AVGVISIB_probability.tif')
OUT_CLASS = os.path.join(OUT_DIR, 'ch02', 'C2_illumination_class.tif')
OUT_F1 = os.path.join(OUT_DIR, 'ch02', 'F1_illumination_norm.tif')

CRS_PROJ4 = ('+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +x_0=0 +y_0=0 '
             '+R=1737400 +units=m +no_defs')
NODATA = -9999.0

# 阈值（GPA 模型变量口径，算法零改动）
TH_SPSR = 1e-6
TH_PSR = 1e-3
TH_CLASS3 = 0.128    # P25
TH_CLASS4 = 0.264    # P75，连续光照区下限
F1_LO, F1_HI = 0.001, 0.264


def decode_probability(raw):
    """AVGVISIB 量纲体检与概率解码：max>1 视为 int16×25000 编码并 ÷25000；
    max≤1 视为已是概率（不重复解码）。返回 (p, 说明)。"""
    raw = np.asarray(raw, dtype=np.float64)
    mx = float(np.nanmax(raw))
    if mx > 1.0:
        return raw / 25000.0, f'max={mx:.1f}>1，判定 int16×25000 编码 → ÷25000'
    return raw.copy(), f'max={mx:.4f}≤1，已是概率量纲 → 不再解码'


def psr_masks(p):
    """双 PSR 分类：返回 dict(spsr, subpsr, psr, continuous)（bool 数组）。"""
    return dict(
        spsr=p < TH_SPSR,
        subpsr=(p >= TH_SPSR) & (p < TH_PSR),
        psr=p < TH_PSR,
        continuous=p >= TH_CLASS4,
    )


def illumination_classes(p):
    """五级光照分类（1–5，float64；与 GPA Con 嵌套逐项等价）。"""
    return np.where(p <= TH_SPSR, 1.0,
           np.where(p <= TH_PSR, 2.0,
           np.where(p <= TH_CLASS3, 3.0,
           np.where(p <= TH_CLASS4, 4.0, 5.0))))


def f1_normalize(p):
    """F1 分段线性归一化（[0,1]；下界 0.001=subPSR 边界，上界 0.264=P75）。"""
    f1 = np.clip((p - F1_LO) / (F1_HI - F1_LO), 0.0, 1.0)
    return np.where(p < F1_LO, 0.0, f1)


def run(input_tif=INPUT_TIF, out_class=OUT_CLASS, out_f1=OUT_F1):
    """主入口：读概率栅格 → 双 PSR/五级分类/F1 → 铁基准校验 → 写盘。"""
    with rasterio.open(input_tif) as src:
        raw = src.read(1).astype(np.float64)
        transform, crs = src.transform, src.crs
    p, note = decode_probability(raw)
    print(f'[Ch02] 输入体检: {input_tif} → {note}')

    masks = psr_masks(p)
    cls = illumination_classes(p)
    f1 = f1_normalize(p)

    print('[Ch02] 覆盖率（历史基准: sPSR 19.5% / subPSR 1.0% / PSR 20.5% / 连续光照 25.0%）:')
    cover = {}
    for k, m in masks.items():
        cover[k] = float(m.mean())
        print(f'    {k:<12} {m.mean():.2%}')
    cls_cover = {int(k): float((cls == k).mean()) for k in range(1, 6)}
    print('[Ch02] 五级分类占比:', {k: f'{v:.2%}' for k, v in cls_cover.items()})

    os.makedirs(os.path.dirname(out_class), exist_ok=True)
    profile = dict(driver='GTiff', height=384, width=384, count=1, dtype='float32',
                   crs=CRS.from_proj4(CRS_PROJ4),
                   transform=transform if abs(transform.a - 240) < 1e-9
                   else from_origin(-46080, 46080, 240, 240),
                   nodata=NODATA, compress='deflate')

    assert_iron_grid('ch02/C2_illumination_class', cls, profile['transform'],
                     profile['crs'], expect_full_valid=True)
    with rasterio.open(out_class, 'w', **profile) as dst:
        dst.write(cls.astype('float32'), 1)
    assert_iron_grid('ch02/F1_illumination_norm', f1, profile['transform'],
                     profile['crs'], expect_full_valid=True)
    with rasterio.open(out_f1, 'w', **profile) as dst:
        dst.write(f1.astype('float32'), 1)
    print(f'[Ch02] 已输出: {out_class}\n        {out_f1}')
    return dict(cover=cover, class_cover=cls_cover,
                class_tif=out_class, f1_tif=out_f1)


if __name__ == '__main__':
    run()
