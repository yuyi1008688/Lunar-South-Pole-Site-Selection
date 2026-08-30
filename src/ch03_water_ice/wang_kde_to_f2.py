# -*- coding: utf-8 -*-
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch03 · 水冰分布（核心算法）：Wang 2025 实测冰点 → 高斯KDE → PSR分级掩膜 → P98归一化 → F2因子
# 来源  : 竞赛提交包 3工程文件/Ch03_*/wang_kde_to_f2.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
Wang 2025 实测冰点 KDE -> F2 水冰丰度因子
=================================================
将 Wang 2025 深度学习 M³ 冰识别点(5031个研究区内冰点)做核密度估计，
叠加 PSR 分级掩膜(sPSR/subPSR)，生成 F2 水冰丰度因子。

数据说明：
  - wang2025: 40,623 行冰点, Latitude/Longitude(-180~180), 研究区内 5031 点
  - AVGVISIB_probability.tif: 平均可见度(0=全年无光照), 384×384 与 PSR 同网格
  - sPSR_mask.tif / subPSR_mask.tif: 分级掩膜(与 avgvisib 阈值等价, 已验证)
  - PSR_mask = sPSR + subPSR (30252 = 28776 + 1476)

分级逻辑(Step 5):
  - sPSR(avgvisib<1e-6): 用完整 KDE 值
  - subPSR(avgvisib in [1e-6,0.001)): 用 KDE×0.8 降权
  - PSR 外部: 0
"""

import os
import numpy as np
import pandas as pd
import rasterio
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ---------------- 路径与常量 ----------------
WORKDIR = os.path.dirname(os.path.abspath(__file__))          # .../01_F2最终方案_WangKDE/脚本
ROOT = os.environ.get('LUNAR_PROJECT_ROOT',
                 os.path.abspath(os.path.join(WORKDIR, '..', '..')))  # 项目根目录（可用环境变量覆盖）
DATA_DIR = os.path.join(ROOT, '00_数据源', '最新数据')
WANG_PATH = os.path.join(DATA_DIR, 'wang2025_ice_pixel_positions_spectra.xlsx')
AVGVISIB_PATH = os.path.join(DATA_DIR, 'AVGVISIB_probability.tif')
SPSR_PATH = os.path.join(DATA_DIR, 'sPSR_mask.tif')
SUBPSR_PATH = os.path.join(DATA_DIR, 'subPSR_mask.tif')
PSR_PATH = os.path.join(ROOT, '00_数据源', 'PSR_mask.tif')
OUT_PATH = os.path.join(WORKDIR, '..', '结果', 'F2_wang_kde_final.tif')
PNG_PATH = os.path.join(WORKDIR, '..', '结果', 'F2_wang_kde_validation.png')

R = 1737400.0
LIM = 46080.0
GRID_N = 384
BW_METHOD = 0.11  # 使实际带宽≈0.11×std≈1600m（任务意图，已验证 scipy 实际带宽=bw×std）


def main():
    # ---------------- Step 1：读取 Wang 2025 数据 ----------------
    print('>>> Step 1 开始：读取Wang 2025数据...')
    df = pd.read_excel(WANG_PATH)
    print(f'  列名前3个: {list(df.columns[:3])}')
    df_roi = df[df['Latitude'] <= -88.5]
    print(f'  研究区内冰点数: {len(df_roi)}')

    # ---------------- Step 2：坐标（用精确投影坐标 POINT_X/Y） ----------------
    # 关键修正：Latitude/Longitude 列精度不足（5031点仅19个唯一值），
    # 而 POINT_X/POINT_Y 是 5031 个唯一精确投影坐标（南极立体，lat_ts=-90，
    # 已验证与 PSR 网格 rho 比值均值 1.004 基本一致），直接使用。
    print('>>> Step 2 开始：坐标处理（用精确投影坐标POINT_X/POINT_Y）...')
    x_ice = df_roi['POINT_X'].values
    y_ice = df_roi['POINT_Y'].values

    valid = (np.abs(x_ice) <= LIM) & (np.abs(y_ice) <= LIM)
    x_valid = x_ice[valid]
    y_valid = y_ice[valid]
    print(f'  有效点: {valid.sum()}（x范围 {x_valid.min():.0f}~{x_valid.max():.0f}m，'
          f'y范围 {y_valid.min():.0f}~{y_valid.max():.0f}m）')

    # ---------------- Step 3：KDE ----------------
    print('>>> Step 3 开始：核密度估计KDE...')
    points = np.vstack([x_valid, y_valid])
    kde = gaussian_kde(points, bw_method=BW_METHOD)

    x_grid = np.linspace(-LIM + 120, LIM - 120, GRID_N)
    # 关键：y 轴必须递减（北→南），与 rasterio 第 0 行=北极一致；
    # 若递增（南→北）会导致栅格南北翻转（与 validation_wang_kde.py 第 71 行写法一致）。
    y_grid = np.linspace(LIM - 120, -LIM + 120, GRID_N)
    xx, yy = np.meshgrid(x_grid, y_grid)
    grid_points = np.vstack([xx.ravel(), yy.ravel()])

    print('  正在计算KDE（可能需要1-2分钟）...')
    kde_values = kde(grid_points).reshape(GRID_N, GRID_N)
    print(f'  KDE值域: {kde_values.min():.6f} ~ {kde_values.max():.6f}')

    # ---------------- Step 4：归一化至[0,1] ----------------
    print('>>> Step 4 开始：归一化...')
    p98 = np.percentile(kde_values, 98)
    kde_norm = np.clip(kde_values / p98, 0, 1)
    print(f'  归一化后值域: {kde_norm.min():.4f} ~ {kde_norm.max():.4f}')

    # ---------------- Step 5：叠合 PSR 分级掩膜 ----------------
    print('>>> Step 5 开始：叠合PSR分级掩膜...')
    with rasterio.open(AVGVISIB_PATH) as src:
        avgvisib = src.read(1)
    with rasterio.open(PSR_PATH) as src:
        psr_mask = src.read(1)
        psr_transform = src.transform
        psr_crs = src.crs

    F2_wang = np.zeros((GRID_N, GRID_N), dtype=np.float32)

    # sPSR：完整 KDE 值
    spsr = avgvisib < 1e-6
    F2_wang[spsr] = kde_norm[spsr]

    # subPSR：KDE×0.8 降权
    subpsr = (avgvisib >= 1e-6) & (avgvisib < 0.001)
    F2_wang[subpsr] = kde_norm[subpsr] * 0.8

    print(f'  F2_wang非零像元: {(F2_wang > 0).sum()}，占比: {(F2_wang > 0).mean()*100:.2f}%')
    print(f'  F2_wang在sPSR内均值: {F2_wang[spsr].mean():.4f}')
    print(f'  F2_wang最大值: {F2_wang.max():.4f}')

    # ---------------- Step 6：输出 GeoTIFF ----------------
    print('>>> Step 6 开始：输出GeoTIFF...')
    with rasterio.open(OUT_PATH, 'w', driver='GTiff',
                       height=GRID_N, width=GRID_N, count=1,
                       dtype='float32', crs=psr_crs, transform=psr_transform,
                       nodata=-9999) as dst:
        dst.write(F2_wang, 1)
    print('  F2_wang_kde_final.tif 已输出')

    # ---------------- Step 7：验证图（1×3） ----------------
    print('>>> Step 7 开始：生成验证图...')
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    axes[0].scatter(x_valid, y_valid, s=1, alpha=0.3, c='blue')
    axes[0].set_title(f'Wang 2025实测冰点\n({valid.sum()}个，研究区内)')
    axes[0].set_xlim(-LIM, LIM)
    axes[0].set_ylim(-LIM, LIM)
    axes[0].set_aspect('equal')

    im2 = axes[1].imshow(kde_norm, origin='upper', cmap='Blues',
                         extent=[-LIM, LIM, -LIM, LIM])
    axes[1].set_title('KDE归一化（未掩膜）\n深蓝=冰点聚集处')
    plt.colorbar(im2, ax=axes[1])

    im3 = axes[2].imshow(F2_wang, origin='upper', cmap='Blues',
                         extent=[-LIM, LIM, -LIM, LIM])
    axes[2].set_title('F2_wang_kde（PSR掩膜后）\n期望：深蓝在PSR内冰点聚集区')
    plt.colorbar(im3, ax=axes[2])

    plt.tight_layout()
    plt.savefig(PNG_PATH, dpi=150, bbox_inches='tight')
    plt.close()
    print('  验证图已保存')

    # ---------------- Step 8：自检 ----------------
    print('>>> Step 8 开始：自检...')
    pct = (F2_wang > 0).mean() * 100
    max_row, max_col = np.unravel_index(F2_wang.argmax(), (GRID_N, GRID_N))
    print('======== Wang KDE F2自检 ========')
    print(f'非零占比：{pct:.2f}%（期望15-25%）')
    print(f'sPSR内均值：{F2_wang[spsr].mean():.4f}')
    print(f'subPSR内均值：{F2_wang[subpsr].mean():.4f}')
    print(f'F2最大值：{F2_wang.max():.4f}')
    print(f'F2最大值坐标：行={max_row}, 列={max_col}')
    print(f'F2最大值在sPSR内：{"[是]" if spsr[max_row, max_col] else "[否]"}')
    print('================================')


if __name__ == '__main__':
    main()
