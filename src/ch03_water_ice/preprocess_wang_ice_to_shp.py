# -*- coding: utf-8 -*-
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch03 · 水冰分布：冰点 xlsx → Shapefile 预处理（使用 POINT_X/Y 精确投影坐标）
# 来源  : 竞赛提交包 3工程文件/Ch03_*/preprocess_wang_ice_to_shp.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
Wang 2025 冰点数据预处理 → 标准 Shapefile 点数据
=================================================
用途：把 Wang 2025 实测冰点（xlsx）转成标准 Shapefile 点数据，
     供核密度分析（GIS 软件 KDE 工具或 Python）使用。

处理逻辑：
  1. 读取 wang2025_ice_pixel_positions_spectra.xlsx
  2. 筛选研究区（Latitude <= -88.5，对应 ±46080 m 范围，得 5031 个点）
  3. 几何坐标用 POINT_X / POINT_Y（米，南极立体投影）
     —— 注意：不能用 Latitude/Longitude 列（精度不足，仅 19 个唯一值）
  4. 坐标系从 PSR_mask.tif 复用（南极立体投影，R=1737400）
  5. 输出 Shapefile 到 交付包/数据/wang_ice_points/

运行：python preprocess_wang_ice_to_shp.py
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import pyproj
from shapely.geometry import Point

# ---------------- 路径（按需修改） ----------------
WORKDIR = os.environ.get('LUNAR_DATA_DIR', os.path.dirname(os.path.abspath(__file__)))  # 可用环境变量覆盖
XLSX_PATH = os.path.join(WORKDIR, '最新数据', 'wang2025_ice_pixel_positions_spectra.xlsx')
PSR_PATH = os.path.join(WORKDIR, 'PSR_mask.tif')
OUT_DIR = os.path.join(WORKDIR, '交付包', '数据', 'wang_ice_points')

LIM = 46080.0  # 研究区半宽（米）


def main():
    print('>>> 读取冰点数据...')
    df = pd.read_excel(XLSX_PATH)
    print(f'  原始记录数：{len(df)}')

    # 筛选研究区（纬度 <= -88.5）
    roi = df[df['Latitude'] <= -88.5].copy()
    print(f'  研究区内冰点数（Latitude<=-88.5）：{len(roi)}')

    # 几何坐标用 POINT_X / POINT_Y（米）
    x = roi['POINT_X'].values
    y = roi['POINT_Y'].values
    valid = (np.abs(x) <= LIM) & (np.abs(y) <= LIM)
    roi = roi[valid].copy()
    x = roi['POINT_X'].values
    y = roi['POINT_Y'].values
    print(f'  坐标落在 ±46080m 内的点：{len(roi)}')
    print(f'  POINT_X 范围：{x.min():.0f} ~ {x.max():.0f} m')
    print(f'  POINT_Y 范围：{y.min():.0f} ~ {y.max():.0f} m')

    # 几何
    geometry = [Point(xi, yi) for xi, yi in zip(x, y)]

    # 属性字段（shapefile .dbf 字段名限 10 字符，Maximum temperature 缩写为 MaxTemp）
    attr = roi[['Latitude', 'Longitude', 'POINT_X', 'POINT_Y',
                'Maximum temperature']].rename(
        columns={'Maximum temperature': 'MaxTemp'})

    # 坐标系：从 PSR_mask.tif 复用（南极立体投影，R=1737400）
    with rasterio.open(PSR_PATH) as src:
        crs_wkt = src.crs.to_wkt()
    crs = pyproj.CRS.from_wkt(crs_wkt)
    print(f'  坐标系：{crs.to_proj4()}')

    gdf = gpd.GeoDataFrame(attr, geometry=geometry, crs=crs)

    # 输出 Shapefile
    os.makedirs(OUT_DIR, exist_ok=True)
    out_shp = os.path.join(OUT_DIR, 'wang_ice_points.shp')
    gdf.to_file(out_shp, driver='ESRI Shapefile', encoding='utf-8')
    print(f'  已输出 Shapefile：{out_shp}')
    print(f'  点数：{len(gdf)}，字段：{list(gdf.columns)}')

    print('\n>>> 预处理完成。交付包点数据已生成。')


if __name__ == '__main__':
    main()
