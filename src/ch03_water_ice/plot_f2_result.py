# -*- coding: utf-8 -*-
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch03 · 水冰分布：F2 最终成果图绘制
# 来源  : 竞赛提交包 3工程文件/Ch03_*/plot_f2_result.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""F2 最终成果图（水冰丰度因子成图）"""
import os
import numpy as np
import rasterio
import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import uniform_filter

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

BASE = os.environ.get('LUNAR_PROJECT_ROOT', './')  # 复现时设置 LUNAR_PROJECT_ROOT 或改为本地路径
F2_PATH = BASE + '/01_F2最终方案_WangKDE/结果/F2_wang_kde_final.tif'
PSR_PATH = BASE + '/00_数据源/PSR_mask.tif'
SHP_PATH = BASE + '/03_交付物/交付包/数据/wang_ice_points/wang_ice_points.shp'
OUT_PATH = BASE + '/01_F2最终方案_WangKDE/结果/F2最终成果图.png'

LIM = 46080.0

# ---- 读取 F2 ----
with rasterio.open(F2_PATH) as src:
    f2 = src.read(1)
    tr = src.transform
    crs = src.crs
f2_ma = np.ma.masked_where(f2 <= 0, f2)  # PSR 外=0 → 背景

# ---- 读取 PSR 掩膜 ----
with rasterio.open(PSR_PATH) as src:
    psr = src.read(1)

# ---- 读取冰点 ----
g = gpd.read_file(SHP_PATH)
x_ice = g.geometry.x.values.astype(float)
y_ice = g.geometry.y.values.astype(float)

# ---- 统计 ----
nz_pct = (f2 > 0).mean() * 100
# 峰值检测：F2 经 P98 归一化后 clip 到 1.0，多个像元饱和，直接 argmax 会标在
# 第一个饱和像元而非真正核心。改用 5×5（约1.2km）邻域均值定位最密集冰点区。
f2_smooth = uniform_filter(f2, size=5)
r, c = np.unravel_index(np.argmax(f2_smooth), f2.shape)
peak_x = -LIM + tr.a * (c + 0.5)
peak_y = LIM + tr.e * (r + 0.5)  # tr.e 为负
peak_val = f2[r, c]

# ---- 自定义 sequential 色带（单 hue 浅→深，避开背景白）----
# magma 的裁剪版：低值淡、高值深紫红，感知均匀
cmap = plt.get_cmap('magma').copy()
cmap.set_bad('white')

fig, ax = plt.subplots(figsize=(10, 9))

im = ax.imshow(f2_ma, cmap=cmap, vmin=0, vmax=1,
               extent=[-LIM, LIM, -LIM, LIM], origin='upper',
               interpolation='nearest')

# PSR 边界轮廓
ax.contour(psr, levels=[0.5], colors='black', linewidths=0.5,
           extent=[-LIM, LIM, -LIM, LIM], origin='upper', alpha=0.55)

# 冰点叠加（淡蓝，展示集聚）
ax.scatter(x_ice, y_ice, s=0.6, c='#4da6ff', alpha=0.25, linewidths=0)

# 峰值标注
ax.plot(peak_x, peak_y, marker='x', color='white', markersize=12,
        markeredgewidth=2.5, markeredgecolor='white')
ax.annotate(f'峰值 F2={peak_val:.3f}\n({peak_x:.0f}, {peak_y:.0f} m)',
            xy=(peak_x, peak_y), xytext=(peak_x - 16000, peak_y - 12000),
            color='white', fontsize=10, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', fc='black', ec='white', alpha=0.7),
            arrowprops=dict(arrowstyle='->', color='white', lw=1.4))

cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.set_label('水冰丰度 F2（0=无水冰，1=最高）', fontsize=11)

ax.set_xlabel('投影坐标 X（米，东向为正）', fontsize=11)
ax.set_ylabel('投影坐标 Y（米，北向为正）', fontsize=11)
ax.set_title('F2 水冰丰度最终成果图（Wang 2025 KDE，PSR 掩膜后）\n'
             f'南极立体投影 | 384×384 @ 240m | 非零像元 {nz_pct:.1f}% | '
             f'冰点 {len(x_ice)} 个 | 峰值 ({peak_x:.0f}, {peak_y:.0f}) m',
             fontsize=12)
ax.set_aspect('equal')
ax.grid(False)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=150, bbox_inches='tight', facecolor='white')
print('已保存:', OUT_PATH)
print(f'非零像元占比: {nz_pct:.2f}%')
print(f'峰值: F2={peak_val:.4f} @ ({peak_x:.0f}, {peak_y:.0f}) m')
