# -*- coding: utf-8 -*-
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch05 · AHP 选址：蒙特卡洛 1000 次权重敏感性分析（Dirichlet 随机权重，位移标准差/重叠率）
# 来源  : 竞赛提交包 3工程文件/Ch05_*/蒙特卡洛1000次敏感性分析.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
S8 蒙特卡洛1000次敏感性分析（修正版）
删除权重的排序硬约束，使用Dirichlet分布生成纯随机权重
输出：最优点散点分布、位移标准差、I级区重叠率估算
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import rasterio

# 图标中文字体
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


# ============================================================
# 1. 配置（请修改为你自己的文件路径）
# ============================================================
# 五个因子层路径（必须全部存在且空间对齐）
# ── 路径配置：复现时改为本地 F1~F5 因子层所在目录，或设环境变量 LUNAR_FACTOR_DIR ──
_factor_dir = os.environ.get("LUNAR_FACTOR_DIR", os.path.join("data", "rasters"))
F1_PATH = os.path.join(_factor_dir, "F1_illumination_norm.tif")
F2_PATH = os.path.join(_factor_dir, "F2_wang_kde_final_1.tif")
F3_PATH = os.path.join(_factor_dir, "F3_fos_safety.tif")
F4_PATH = os.path.join(_factor_dir, "F4_ecsa_sync.tif")
F5_PATH = os.path.join(_factor_dir, "F5_continuous_distance.tif")

# 掩膜（硬约束后的有效区域）—— 可选，如果提供则只在掩膜内找最优点
MASK_PATH = os.path.join(_factor_dir, "suitability_final.tif")  # 可选掩膜

# 输出目录
OUTPUT_DIR = os.environ.get("LUNAR_OUTPUT_DIR", os.path.join("data", "output", "ch05"))
os.makedirs(OUTPUT_DIR, exist_ok=True)  # 输出目录自建（原稿假设目录已存在）

# 蒙特卡洛迭代次数
N_ITER = 1000
RANDOM_SEED = 42

# ============================================================
# 2. 读取数据
# ============================================================
print("正在读取因子层...")
with rasterio.open(F1_PATH) as src:
    F1 = src.read(1).astype(np.float32)
    profile = src.profile
    transform = src.transform
    rows, cols = F1.shape

with rasterio.open(F2_PATH) as src:
    F2 = src.read(1).astype(np.float32)
with rasterio.open(F3_PATH) as src:
    F3 = src.read(1).astype(np.float32)
with rasterio.open(F4_PATH) as src:
    F4 = src.read(1).astype(np.float32)
with rasterio.open(F5_PATH) as src:
    F5 = src.read(1).astype(np.float32)

# 读取掩膜（如果有）
if os.path.exists(MASK_PATH):
    with rasterio.open(MASK_PATH) as src:
        mask = src.read(1).astype(np.uint8)
        # 强制用F1有效值作为掩膜，忽略外部掩膜
        valid_mask = (~np.isnan(F1)) & (F1 > 0)
else:
    valid_mask = (~np.isnan(F1)) & (F1 > 0)

print(f"有效像元数: {np.sum(valid_mask)}")

# 提取有效像元的坐标（行、列）
coords = np.array(np.where(valid_mask)).T  # (N, 2)
valid_F1 = F1[valid_mask]
valid_F2 = F2[valid_mask]
valid_F3 = F3[valid_mask]
valid_F4 = F4[valid_mask]
valid_F5 = F5[valid_mask]

# ============================================================
# 3. 蒙特卡洛主循环（无排序约束）
# ============================================================
print("开始蒙特卡洛采样（无排序约束）...")
np.random.seed(RANDOM_SEED)
best_points = []  # 存储每次迭代的最优点坐标 (x, y)

for i in range(N_ITER):
    # 使用Dirichlet分布生成随机权重，总和为1，无排序限制
    weights = np.random.dirichlet(np.ones(5))  # 均匀分布在单纯形上

    # 计算加权得分
    weighted = (weights[0] * valid_F1 +
                weights[1] * valid_F2 +
                weights[2] * valid_F3 +
                weights[3] * valid_F4 +
                weights[4] * valid_F5)

    # 找最大得分索引
    idx_max = np.argmax(weighted)
    row, col = coords[idx_max]
    # 转换为地理坐标（像素中心）
    # 注：坐标公式按方法论 §5.7 P0 修复记录采用 x = tf[2] + col×tf[0]、
    #     y = tf[5] + row×tf[4]（tf[0]=像宽, tf[2]=X原点）。原稿此处曾把
    #     tf[0]/tf[2] 混淆，导致坐标跑到百万级（P0 bug），已修正。
    x = transform[2] + (col + 0.5) * transform[0]
    y = transform[5] + (row + 0.5) * transform[4]
    best_points.append([x, y])

best_points = np.array(best_points)  # (N_ITER, 2)

# ============================================================
# 4. 统计结果
# ============================================================
x_mean, y_mean = np.mean(best_points, axis=0)
x_std, y_std = np.std(best_points, axis=0, ddof=1)
disp_std = np.sqrt(x_std ** 2 + y_std ** 2)

print("=" * 50)
print("蒙特卡洛结果（无排序约束）")
print("=" * 50)
print(f"迭代次数: {N_ITER}")
print(f"最优点X均值: {x_mean:.2f} m, 标准差: {x_std:.2f} m")
print(f"最优点Y均值: {y_mean:.2f} m, 标准差: {y_std:.2f} m")
print(f"综合位移标准差: {disp_std:.2f} m")
print(f"稳健性判定: {'稳健' if disp_std < 500 else '敏感'}")

# 估算I级区空间重叠率（以基准权重的最优点为中心，半径500m内的点占比）
base_weights = np.array([0.331, 0.258, 0.191, 0.126, 0.094])
weighted_base = (base_weights[0] * valid_F1 +
                 base_weights[1] * valid_F2 +
                 base_weights[2] * valid_F3 +
                 base_weights[3] * valid_F4 +
                 base_weights[4] * valid_F5)
idx_base = np.argmax(weighted_base)
row_base, col_base = coords[idx_base]
x_base = transform[2] + (col_base + 0.5) * transform[0]
y_base = transform[5] + (row_base + 0.5) * transform[4]

# 计算到基准点的距离
distances = np.sqrt((best_points[:, 0] - x_base) ** 2 + (best_points[:, 1] - y_base) ** 2)
overlap_500 = np.sum(distances <= 500) / N_ITER * 100
overlap_1000 = np.sum(distances <= 1000) / N_ITER * 100

print(f"基准最优点坐标: ({x_base:.1f}, {y_base:.1f})")
print(f"落在500m内的迭代次数占比: {overlap_500:.1f}%")
print(f"落在1000m内的迭代次数占比: {overlap_1000:.1f}%")
print("=" * 50)

# ============================================================
# 5. 可视化
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 子图1：散点图
ax = axes[0, 0]
ax.scatter(best_points[:, 0], best_points[:, 1], s=2, alpha=0.4, c='steelblue')
ax.scatter(x_base, y_base, c='red', s=80, marker='*', label='基准最优点', edgecolors='black')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_title(f'1000次蒙特卡洛最优点分布\n位移标准差={disp_std:.1f}m')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')

# 子图2：X/Y位移直方图
ax = axes[0, 1]
ax.hist(best_points[:, 0], bins=30, alpha=0.5, label='X位移', color='blue')
ax.hist(best_points[:, 1], bins=30, alpha=0.5, label='Y位移', color='orange')
ax.axvline(x_mean, color='blue', linestyle='--', label='X均值')
ax.axvline(y_mean, color='orange', linestyle='--', label='Y均值')
ax.axvline(x_mean - 500, color='gray', linestyle=':', alpha=0.7)
ax.axvline(x_mean + 500, color='gray', linestyle=':', alpha=0.7, label='±500m')
ax.set_xlabel('位移 (m)')
ax.set_ylabel('频数')
ax.set_title('最优点坐标位移分布')
ax.legend()

# 子图3：距离直方图（到基准点的距离）
ax = axes[1, 0]
ax.hist(distances, bins=30, color='green', alpha=0.7)
ax.axvline(500, color='red', linestyle='--', label='500m')
ax.axvline(1000, color='orange', linestyle='--', label='1000m')
ax.set_xlabel('到基准最优点距离 (m)')
ax.set_ylabel('频数')
ax.set_title(f'距离基准点分布 (500m内: {overlap_500:.1f}%)')
ax.legend()

# 子图4：重叠率估算
ax = axes[1, 1]
labels = ['基准最优点附近500m', '基准最优点附近1000m']
overlaps = [overlap_500, overlap_1000]
bars = ax.bar(labels, overlaps, color=['lightcoral', 'lightblue'])
ax.set_ylim(0, 100)
ax.set_ylabel('重叠率 (%)')
ax.set_title('I级区空间重叠率估算')
for bar, val in zip(bars, overlaps):
    ax.text(bar.get_x() + bar.get_width() / 2, val + 2, f'{val:.1f}%', ha='center', fontsize=12)

plt.tight_layout()
save_path = os.path.join(OUTPUT_DIR, 'monte_carlo_s8_fixed.png')
plt.savefig(save_path, dpi=150, bbox_inches='tight')
plt.show()
print(f"✅ 可视化保存至: {save_path}")

# 保存坐标数据
np.savetxt(os.path.join(OUTPUT_DIR, 'monte_carlo_points_fixed.csv'),
           best_points, delimiter=',', header='x,y', comments='')

print("✅ 蒙特卡洛分析完成！")