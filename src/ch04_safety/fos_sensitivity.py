# ──────────────────────────────────────────────────────────────────
# 章节  : Ch04 · 安全势场：FoS 参数敏感性分析（C×φ 5×5=25组网格扫描，F3>0.8面积变化率）
# 来源  : 竞赛提交包 3工程文件/Ch04_*/fos_sensitivity.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
import rasterio
from rasterio.profiles import DefaultGTiffProfile
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import os
# 强制指定黑体作为默认字体
plt.rcParams['font.sans-serif'] = ['simhei']   # 指定黑体
plt.rcParams['axes.unicode_minus'] = False     # 解决负号显示为方块的问题

# ============================================================
# 配置
# ============================================================
# 输入数据路径（请改成你电脑上的实际路径）
SLOPE_PATH = os.environ.get('LUNAR_SLOPE_TIF', os.path.join('data', 'rasters', 'slope_deg.tif'))  # 坡度图（单位：度）
OUTPUT_DIR = os.environ.get('LUNAR_OUTPUT_DIR', os.path.join('data', 'output', 'ch04_sensitivity'))  # 输出目录

# FoS固定参数（不变）
RHO = 1500.0  # 月壤体积密度 (kg/m³)
G = 1.62  # 月球重力加速度 (m/s²)
H = 1.0  # 滑移面深度 (m)

# 参数网格
C_LIST = [0.5, 1.0, 1.5, 2.0, 3.0]  # 凝聚力 (kPa)
PHI_LIST = [30, 32, 35, 38, 40]  # 内摩擦角 (°)

BASE_C = 1.5  # 基准参数
BASE_PHI = 35

# ============================================================
# 读取坡度数据
# ============================================================
with rasterio.open(SLOPE_PATH) as src:
    slope_deg = src.read(1).astype(np.float64)
    meta = src.meta.copy()
    nodata = src.nodata

valid = (slope_deg != nodata) & ~np.isnan(slope_deg)
slope_rad = np.radians(slope_deg)

# 预计算 cos 和 sin 提升性能
cos_slope = np.cos(slope_rad)
sin_slope = np.sin(slope_rad)
cos_slope_sq = cos_slope ** 2

# 平地掩膜 (坡度 < 0.5°)
flat_mask = slope_deg < 0.5

# 固定参数组合
rho_g_H = RHO * G * H

# 存储结果
results = []

# ============================================================
# 主循环：遍历25组参数
# ============================================================
print("=" * 60)
print("开始 FoS 敏感性分析 (25组参数)")
print("=" * 60)

for C in C_LIST:
    for phi_deg in PHI_LIST:
        phi_rad = np.radians(phi_deg)
        tan_phi = np.tan(phi_rad)

        # FoS公式（向量化运算，速度极快）
        numerator = C * 1000 + rho_g_H * cos_slope_sq * tan_phi  # C转Pa
        denominator = rho_g_H * sin_slope * cos_slope

        # 防止分母为0
        denom_safe = np.where(denominator == 0, 1e-10, denominator)
        fos = numerator / denom_safe

        # 平地处理（坡度<0.5°赋3.0）
        fos = np.where(flat_mask, 3.0, fos)
        fos = np.where(~valid, np.nan, fos)

        # 截断归一化 → F3
        fos_clipped = np.clip(fos, 0.5, 3.0)
        F3 = (fos_clipped - 0.5) / 2.5
        F3 = np.where(~valid, np.nan, F3)

        # 统计 F3 > 0.8 的像元数
        high_area_pixels = np.nansum(F3 > 0.8)
        high_area_km2 = high_area_pixels * (240 * 240) / 1e6  # 转换为km²

        results.append({
            'C_kPa': C,
            'phi_deg': phi_deg,
            'F3_gt_0.8_pixels': int(high_area_pixels),
            'F3_gt_0.8_km2': high_area_km2
        })

# ============================================================
# 转换为DataFrame
# ============================================================
df = pd.DataFrame(results)
df_pivot = df.pivot(index='C_kPa', columns='phi_deg', values='F3_gt_0.8_km2')

print("\n✅ 敏感性分析完成！")
print(df_pivot.round(2))

# ============================================================
# 计算面积变化率（以基准参数 C=1.5, phi=35° 为参照）
# ============================================================
base_area = df_pivot.loc[BASE_C, BASE_PHI]
df_pivot_pct = (df_pivot - base_area) / base_area * 100

print("\n📊 面积变化率 (%) —— 相对于基准 (C=1.5, φ=35°)")
print(df_pivot_pct.round(2))

# ============================================================
# 保存CSV
# ============================================================
os.makedirs(OUTPUT_DIR, exist_ok=True)
df_pivot.to_csv(os.path.join(OUTPUT_DIR, 'sensitivity_area_km2.csv'))
df_pivot_pct.to_csv(os.path.join(OUTPUT_DIR, 'sensitivity_area_pct.csv'))

# ============================================================
# 绘制热力图
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# 图1：面积热力图 (km²)
sns.heatmap(df_pivot, annot=True, fmt='.1f', cmap='YlOrRd',
            xticklabels=PHI_LIST, yticklabels=C_LIST,
            ax=ax1, cbar_kws={'label': 'F3>0.8 面积 (km²)'})
ax1.set_xlabel('内摩擦角 φ (°)')
ax1.set_ylabel('凝聚力 C (kPa)')
ax1.set_title('FoS敏感性分析 — F3>0.8面积')

# 图2：变化率热力图 (%)
sns.heatmap(df_pivot_pct, annot=True, fmt='.1f', cmap='RdBu_r',
            xticklabels=PHI_LIST, yticklabels=C_LIST,
            ax=ax2, center=0, cbar_kws={'label': '变化率 (%)'})
ax2.set_xlabel('内摩擦角 φ (°)')
ax2.set_ylabel('凝聚力 C (kPa)')
ax2.set_title('FoS敏感性分析 — 相对基准面积变化率')

# 标注基准点
ax2.add_patch(plt.Rectangle((PHI_LIST.index(BASE_PHI), C_LIST.index(BASE_C)),
                            1, 1, fill=False, edgecolor='green', linewidth=3))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sensitivity_heatmap.png'), dpi=300, bbox_inches='tight')
plt.show()

print(f"\n✅ 结果已保存至: {OUTPUT_DIR}")
print("   文件: sensitivity_area_km2.csv, sensitivity_area_pct.csv, sensitivity_heatmap.png")