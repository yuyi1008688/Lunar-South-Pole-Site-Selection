# ──────────────────────────────────────────────────────────────────
# 章节  : Ch04 · 安全势场主执行（V3.2 终版）：ECSA诊断 / d_max锚定 / 曲率CDF拐点 / FoS→F3 / VRM / 25组敏感性
# 来源  : 竞赛提交包 3工程文件/Ch04_*/ch04_v32_execution.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
第4章 V3.2 综合执行脚本
========================
整合 Claude + Gemini 两位专家方案、经方案评审后的最终执行代码。
用途：月球南极科研站选址 · 安全性与通信覆盖分析

模块清单：
  Step 0: ECSA独立性诊断（必须第一个运行，决定F4方案走向）
  Step 1: d_max 物理锚定（数据驱动确定危险距离场饱和距离）
  Step 2: 曲率CDF拐点分析（确定形态学危险阈值）
  Step 3: FoS边坡安全系数计算 → F3因子层
  Step 4: VRM矢量崎岖度计算（5×5主产出 + 3×3备用）→ 第6章成本面
  Step 5: FoS敏感性分析（5×5全网格 = 25组参数组合）

依赖：numpy, rasterio, scipy, geopandas, matplotlib
输入：slope.tif, aspect.tif, DEM.tif, profile_curvature.tif,
      AVGVISIB.tif, AVGVISIB_EARTH.tif, robins_craters.shp
输出：F3_fos_safety.tif, vrm_5x5.tif, vrm_3x3.tif,
      fos_sensitivity_report, curvature_cdf_plot.png

作者：yuyi（方案设计与执行）
日期：2026-07-15
"""

import numpy as np
import rasterio
from scipy import stats
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # 无头模式
import geopandas as gpd
from itertools import product
import os
import sys
import time

# ============================================================
# 全局配置（与 config.yaml V3.2 一致）
# ============================================================
CONFIG = {
    # 文件路径（相对当前工作目录解析；复现时改为本地实际路径）
    'slope_path': 'data/中间数据/slope.tif',
    'aspect_path': 'data/中间数据/aspect.tif',
    'dem_path': 'data/中间数据/CH01_DEM_fused_v1.tif',
    'curvature_path': 'data/中间数据/profile_curvature.tif',
    'avgvisib_path': 'data/原始数据/AVGVISIB_65S_240M.tif',
    'avgvisib_earth_path': 'data/原始数据/AVGVISIB_EARTH_65S_240M.tif',
    'robins_path': 'data/原始数据/robins_craters.shp',
    
    # 输出路径
    'output_dir': 'data/最终产出/ch04/',
    
    # FoS参数（Mitchell et al. 1972, Apollo 14）
    'C_pa': 1500.0,           # 凝聚力 Pa
    'rho': 1500.0,            # 体积密度 kg/m³
    'phi_deg': 35.0,          # 内摩擦角
    'H_m': 1.0,               # 滑移面深度 m
    'g_moon': 1.62,           # 月球重力 m/s²
    'flat_threshold_deg': 0.5,  # 平地保护阈值
    'fos_range': (0.5, 3.0),  # 归一化截断区间
    
    # VRM参数
    'vrm_primary_window': 5,
    'vrm_secondary_window': 3,
    
    # 研究区边界
    'lat_min': -90.0,
    'lat_max': -88.5,
    
    # 分辨率
    'pixel_size': 240,
}


def read_raster(path):
    """读取栅格文件，返回 (data, meta, nodata)"""
    with rasterio.open(path) as src:
        data = src.read(1).astype(np.float64)
        meta = src.meta.copy()
        nodata = src.nodata
    if nodata is not None:
        data = np.where(data == nodata, np.nan, data)
    return data, meta, nodata


def save_raster(data, meta, output_path, nodata=-9999.0):
    """保存栅格文件"""
    data_save = np.where(np.isnan(data), nodata, data)
    meta.update(dtype=rasterio.float32, nodata=nodata)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with rasterio.open(output_path, 'w', **meta) as dst:
        dst.write(data_save.astype(np.float32), 1)
    print(f"  → 已保存: {output_path}")


# ============================================================
# Step 0: ECSA 独立性诊断（必须第一个运行！）
# ============================================================
def step0_ecsa_diagnostic():
    """
    ECSA独立性假设诊断
    决定F4因子采用哪条路径：
      路径A: CV<5% → F4移除
      路径B: CV≥5% 且 |r|<0.5 → ECSA乘积
      路径C: CV≥5% 且 |r|≥0.5 → 退回AVGVISIB_EARTH
    """
    print("\n" + "="*60)
    print("Step 0: ECSA 独立性诊断")
    print("="*60)
    
    light, _, _ = read_raster(CONFIG['avgvisib_path'])
    comm, _, _ = read_raster(CONFIG['avgvisib_earth_path'])
    
    # 展平并去除NaN
    light_flat = light.ravel()
    comm_flat = comm.ravel()
    valid = ~(np.isnan(light_flat) | np.isnan(comm_flat))
    light_v = light_flat[valid]
    comm_v = comm_flat[valid]
    
    # 统计量
    r_pearson, p_pearson = stats.pearsonr(light_v, comm_v)
    r_spearman, p_spearman = stats.spearmanr(light_v, comm_v)
    comm_cv = np.nanstd(comm_v) / np.nanmean(comm_v)
    light_cv = np.nanstd(light_v) / np.nanmean(light_v)
    
    print(f"\n  AVGVISIB (光照):")
    print(f"    均值 = {np.nanmean(light_v):.4f}")
    print(f"    标准差 = {np.nanstd(light_v):.4f}")
    print(f"    变异系数 CV = {light_cv:.4f}")
    print(f"\n  AVGVISIB_EARTH (通信):")
    print(f"    均值 = {np.nanmean(comm_v):.4f}")
    print(f"    标准差 = {np.nanstd(comm_v):.4f}")
    print(f"    变异系数 CV = {comm_cv:.4f}")
    print(f"\n  空间相关性:")
    print(f"    Pearson r  = {r_pearson:.4f} (p = {p_pearson:.4e})")
    print(f"    Spearman ρ = {r_spearman:.4f} (p = {p_spearman:.4e})")
    
    # 决策路径
    print(f"\n  {'─'*50}")
    if comm_cv < 0.05:
        print("  >>> 路径A: AVGVISIB_EARTH 空间变异极小 (CV < 5%)")
        print("  >>> F4因子无区分度，建议从AHP中移除")
        print("  >>> 权重0.12分配给F3(+0.06)和F5(+0.06)")
        decision = 'A'
    elif abs(r_pearson) >= 0.5:
        print(f"  >>> 路径C: 光照与通信显著相关 (r = {r_pearson:.3f})")
        print("  >>> 独立性假设不严格成立，退回AVGVISIB_EARTH单独作F4")
        decision = 'C'
    else:
        print(f"  >>> 路径B: 独立性假设基本成立 (CV = {comm_cv:.3f}, |r| = {abs(r_pearson):.3f})")
        print("  >>> ECSA乘积模型可用: P_sync = AVGVISIB × AVGVISIB_EARTH")
        decision = 'B'
    
    print(f"  {'─'*50}")
    
    return decision, {
        'comm_cv': comm_cv,
        'light_cv': light_cv,
        'r_pearson': r_pearson,
        'r_spearman': r_spearman,
    }


# ============================================================
# Step 1: d_max 物理锚定
# ============================================================
def step1_dmax_determination():
    """
    基于研究区内最大撞击坑CEB范围确定d_max
    d_max = max(D) × 0.75，上限5000m
    """
    print("\n" + "="*60)
    print("Step 1: d_max 物理锚定")
    print("="*60)
    
    craters = gpd.read_file(CONFIG['robins_path'])
    
    # 筛选研究区内
    lat_min = CONFIG['lat_min']
    lat_max = CONFIG['lat_max']
    
    # 尝试多种可能的字段名
    lat_col = None
    for col in ['LAT_CIRC_IMG', 'LAT', 'lat', 'latitude', 'LATITUDE']:
        if col in craters.columns:
            lat_col = col
            break
    
    if lat_col is None:
        print(f"  警告: 未找到纬度字段，可用字段: {list(craters.columns)}")
        print("  请手动指定字段名")
        return None
    
    diam_col = None
    for col in ['DIAM_CIRC_IMG', 'DIAM', 'diam', 'diameter', 'DIAMETER']:
        if col in craters.columns:
            diam_col = col
            break
    
    if diam_col is None:
        print(f"  警告: 未找到直径字段，可用字段: {list(craters.columns)}")
        return None
    
    study_craters = craters[
        (craters[lat_col] >= lat_min) & (craters[lat_col] <= lat_max)
    ].copy()
    
    # 直径换算为米（如果原始单位是km）
    study_craters['DIAM_M'] = study_craters[diam_col] * 1000
    
    max_diam_m = study_craters['DIAM_M'].max()
    max_idx = study_craters['DIAM_M'].idxmax()
    
    # 尝试获取坑名
    crater_name = "未知"
    for col in ['CRATER_NAME', 'NAME', 'name', 'crater_name']:
        if col in study_craters.columns:
            crater_name = study_craters.loc[max_idx, col]
            break
    
    d_max_physical = max_diam_m * 0.75
    d_max_capped = min(d_max_physical, 5000.0)
    
    # 向上取整到像元整数倍
    pixel = CONFIG['pixel_size']
    d_max_rounded = np.ceil(d_max_capped / pixel) * pixel
    
    print(f"\n  研究区内最大撞击坑: {crater_name}")
    print(f"  直径: {max_diam_m/1000:.1f} km")
    print(f"  CEB外边界半径 (0.75×D): {d_max_physical:.0f} m")
    if d_max_physical > 5000:
        print(f"  超过上限5000m，截断为5000m")
    print(f"  取整到像元边界: {d_max_rounded:.0f} m ({d_max_rounded/pixel:.0f}个像元)")
    print(f"\n  >>> 建议 d_max = {d_max_rounded:.0f} m")
    
    # Shackleton参考
    shackleton_ceb = 21000 * 0.75
    print(f"\n  参考: Shackleton (~21km) CEB半径 = {shackleton_ceb:.0f} m")
    
    return d_max_rounded


# ============================================================
# Step 2: 曲率CDF拐点分析
# ============================================================
def step2_curvature_threshold():
    """
    融合方案：CDF二阶导拐点法 + 高坡度区核密度双峰验证
    """
    print("\n" + "="*60)
    print("Step 2: 曲率CDF拐点分析")
    print("="*60)
    
    curvature, _, _ = read_raster(CONFIG['curvature_path'])
    slope, _, _ = read_raster(CONFIG['slope_path'])
    
    curv_flat = curvature.ravel()
    slope_flat = slope.ravel()
    valid = ~(np.isnan(curv_flat) | np.isnan(slope_flat))
    curv_v = curv_flat[valid]
    slope_v = slope_flat[valid]
    
    # ---- 方法1: 全区曲率CDF二阶导拐点 ----
    curv_sorted = np.sort(curv_v)
    n = len(curv_sorted)
    cdf = np.arange(1, n+1) / n
    
    # 平滑CDF曲线
    from scipy.ndimage import gaussian_filter1d
    cdf_smooth = gaussian_filter1d(cdf, sigma=max(1, n//200))
    
    # 二阶导数
    d1 = np.gradient(cdf_smooth, curv_sorted)
    d2 = np.gradient(d1, curv_sorted)
    
    # 找二阶导过零点（CDF拐点）
    # 在高百分位区间（P90-P99）搜索
    p90_idx = np.searchsorted(curv_sorted, np.percentile(curv_v, 90))
    p99_idx = np.searchsorted(curv_sorted, np.percentile(curv_v, 99))
    
    inflection_candidates = []
    for i in range(p90_idx, min(p99_idx, len(d2)-1)):
        if d2[i] * d2[i+1] < 0:  # 二阶导变号
            inflection_candidates.append(curv_sorted[i])
    
    if inflection_candidates:
        threshold_method1 = inflection_candidates[-1]  # 取最后一个拐点
        print(f"\n  方法1 (CDF拐点法):")
        print(f"    在P90-P99区间找到 {len(inflection_candidates)} 个拐点")
        print(f"    最终阈值 = {threshold_method1:.6f}")
    else:
        threshold_method1 = np.percentile(curv_v, 95)
        print(f"\n  方法1 (CDF拐点法):")
        print(f"    未找到明显拐点，退回P95")
        print(f"    阈值 = {threshold_method1:.6f}")
    
    # ---- 方法2: 高坡度区曲率核密度双峰验证 ----
    high_slope_mask = slope_v > 15.0
    if np.sum(high_slope_mask) > 100:
        curv_high_slope = curv_v[high_slope_mask]
        
        # 核密度估计
        kde = stats.gaussian_kde(curv_high_slope)
        curv_range = np.linspace(np.percentile(curv_high_slope, 80),
                                  np.percentile(curv_high_slope, 99.5), 500)
        density = kde(curv_range)
        
        # 找密度的局部最小值（双峰分布的谷）
        from scipy.signal import argrelextrema
        minima_idx = argrelextrema(density, np.less, order=10)[0]
        
        if len(minima_idx) > 0:
            # 取最靠近P95的谷值
            p95_curv = np.percentile(curv_high_slope, 95)
            best_minima = minima_idx[np.argmin(np.abs(curv_range[minima_idx] - p95_curv))]
            threshold_method2 = curv_range[best_minima]
            print(f"\n  方法2 (高坡度区双峰验证):")
            print(f"    高坡度(>15°)像元数: {np.sum(high_slope_mask)}")
            print(f"    双峰谷值位置 = {threshold_method2:.6f}")
        else:
            threshold_method2 = None
            print(f"\n  方法2 (高坡度区双峰验证):")
            print(f"    未检测到明显双峰结构")
    else:
        threshold_method2 = None
        print(f"\n  方法2 (高坡度区双峰验证):")
        print(f"    高坡度(>15°)像元不足100个，跳过")
    
    # ---- 交叉验证 ----
    print(f"\n  {'─'*50}")
    if threshold_method2 is not None:
        diff_pct = abs(threshold_method1 - threshold_method2) / threshold_method1 * 100
        print(f"  交叉验证: 两种方法差异 = {diff_pct:.1f}%")
        if diff_pct < 10:
            final_threshold = (threshold_method1 + threshold_method2) / 2
            print(f"  差异<10%，取平均值作为最终阈值: {final_threshold:.6f}")
        else:
            final_threshold = threshold_method1
            print(f"  差异>10%，以CDF拐点法为准: {final_threshold:.6f}")
    else:
        final_threshold = threshold_method1
        print(f"  仅CDF拐点法可用，最终阈值: {final_threshold:.6f}")
    
    print(f"  {'─'*50}")
    
    # ---- 可视化 ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # 左图: CDF + 二阶导
    ax1 = axes[0]
    ax1.plot(curv_sorted[p90_idx:p99_idx], cdf_smooth[p90_idx:p99_idx], 'b-', linewidth=2, label='CDF (smoothed)')
    ax1_twin = ax1.twinx()
    ax1_twin.plot(curv_sorted[p90_idx:p99_idx], d2[p90_idx:p99_idx], 'r--', linewidth=1, alpha=0.7, label='d²CDF/dx²')
    ax1.axvline(x=final_threshold, color='green', linestyle='-', linewidth=2, label=f'Threshold = {final_threshold:.4f}')
    ax1.set_xlabel('Profile Curvature')
    ax1.set_ylabel('CDF', color='b')
    ax1_twin.set_ylabel('Second Derivative', color='r')
    ax1.set_title('CDF Inflection Method')
    ax1.legend(loc='upper left')
    ax1_twin.legend(loc='upper right')
    
    # 右图: 高坡度区核密度
    ax2 = axes[1]
    if threshold_method2 is not None:
        ax2.plot(curv_range, density, 'b-', linewidth=2)
        ax2.axvline(x=threshold_method2, color='orange', linestyle='--', linewidth=2, 
                     label=f'Bimodal valley = {threshold_method2:.4f}')
        ax2.axvline(x=final_threshold, color='green', linestyle='-', linewidth=2, 
                     label=f'Final threshold = {final_threshold:.4f}')
        ax2.legend()
    else:
        ax2.plot(curv_range, density, 'b-', linewidth=2)
        ax2.axvline(x=final_threshold, color='green', linestyle='-', linewidth=2, 
                     label=f'Threshold (CDF only) = {final_threshold:.4f}')
        ax2.legend()
    ax2.set_xlabel('Profile Curvature (slope > 15°)')
    ax2.set_ylabel('Kernel Density')
    ax2.set_title('High-Slope Curvature Distribution')
    
    plt.tight_layout()
    plot_path = os.path.join(CONFIG['output_dir'], 'curvature_cdf_inflection.png')
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  可视化已保存: {plot_path}")
    
    return final_threshold


# ============================================================
# Step 3: FoS 边坡安全系数计算 → F3
# ============================================================
def step3_fos_calculation():
    """
    无限边坡安全系数模型
    FoS = (C + ρ·g·H·cos²θ·tanφ) / (ρ·g·H·sinθ·cosθ)
    输出: F3_fos_safety.tif (归一化到[0,1])
    """
    print("\n" + "="*60)
    print("Step 3: FoS 边坡安全系数计算 → F3")
    print("="*60)
    
    slope, meta, nodata = read_raster(CONFIG['slope_path'])
    slope_rad = np.deg2rad(slope)
    
    # 物理参数
    C = CONFIG['C_pa']
    rho = CONFIG['rho']
    g = CONFIG['g_moon']
    H = CONFIG['H_m']
    phi = np.deg2rad(CONFIG['phi_deg'])
    flat_thresh = np.deg2rad(CONFIG['flat_threshold_deg'])
    fos_min, fos_max = CONFIG['fos_range']
    
    # 防除零保护
    flat_mask = slope < flat_thresh
    invalid_mask = np.isnan(slope)
    
    sin_cos = np.sin(slope_rad) * np.cos(slope_rad)
    sin_cos_safe = np.where(flat_mask | invalid_mask, 1.0, sin_cos)
    
    # FoS计算
    numerator = C + rho * g * H * (np.cos(slope_rad)**2) * np.tan(phi)
    denominator = rho * g * H * sin_cos_safe
    
    fos = numerator / denominator
    fos = np.where(flat_mask, fos_max, fos)  # 平地赋最大值
    
    # 统计
    valid_fos = fos[~invalid_mask]
    print(f"\n  FoS 统计:")
    print(f"    最小值: {np.nanmin(valid_fos):.3f}")
    print(f"    最大值: {np.nanmax(valid_fos):.3f}")
    print(f"    均值:   {np.nanmean(valid_fos):.3f}")
    print(f"    FoS < 1.0 (不稳定) 像元比例: {np.nansum(valid_fos < 1.0) / len(valid_fos) * 100:.2f}%")
    print(f"    FoS > 3.0 (极安全) 像元比例: {np.nansum(valid_fos > fos_max) / len(valid_fos) * 100:.2f}%")
    
    # 归一化: 截断[fos_min, fos_max] → [0, 1]
    fos_clipped = np.clip(fos, fos_min, fos_max)
    f3 = (fos_clipped - fos_min) / (fos_max - fos_min)
    f3[invalid_mask] = np.nan
    
    # 保存
    output_path = os.path.join(CONFIG['output_dir'], 'F3_fos_safety.tif')
    save_raster(f3, meta, output_path)
    
    print(f"\n  F3因子层已生成: {output_path}")
    print(f"  F3 > 0.8 区域面积: {np.nansum(f3 > 0.8) * (CONFIG['pixel_size']**2) / 1e6:.2f} km²")
    
    return f3


# ============================================================
# Step 4: VRM 矢量崎岖度计算
# ============================================================
def step4_vrm_calculation():
    """
    VRM = 1 - |Σn_i| / N
    其中 n_i 为邻域内各像元的三维地形法向量
    主产出: 5×5窗口（1200m）
    备用: 3×3窗口（720m）
    """
    print("\n" + "="*60)
    print("Step 4: VRM 矢量崎岖度计算")
    print("="*60)
    
    dem, meta, nodata = read_raster(CONFIG['dem_path'])
    slope, _, _ = read_raster(CONFIG['slope_path'])
    aspect, _, _ = read_raster(CONFIG['aspect_path'])
    
    slope_rad = np.deg2rad(slope)
    aspect_rad = np.deg2rad(aspect)
    res = CONFIG['pixel_size']
    
    # 三维方向余弦
    x = np.sin(slope_rad) * np.sin(aspect_rad)
    y = np.sin(slope_rad) * np.cos(aspect_rad)
    z = np.cos(slope_rad)
    
    invalid_mask = np.isnan(dem)
    
    for window, label in [(5, '5x5'), (3, '3x3')]:
        print(f"\n  计算 VRM {label} (窗口 {window}×{window}, "
              f"覆盖 {window*res:.0f}m×{window*res:.0f}m)...")
        
        n = window ** 2
        sum_x = uniform_filter(x, size=window, mode='nearest') * n
        sum_y = uniform_filter(y, size=window, mode='nearest') * n
        sum_z = uniform_filter(z, size=window, mode='nearest') * n
        
        R = np.sqrt(sum_x**2 + sum_y**2 + sum_z**2)
        vrm = 1.0 - (R / n)
        vrm = np.clip(vrm, 0, 1)
        vrm[invalid_mask] = np.nan
        
        output_path = os.path.join(CONFIG['output_dir'], f'vrm_{label}.tif')
        save_raster(vrm, meta, output_path)
        
        valid_vrm = vrm[~invalid_mask]
        print(f"    VRM 值域: {np.nanmin(valid_vrm):.4f} ~ {np.nanmax(valid_vrm):.4f}")
        print(f"    VRM 均值: {np.nanmean(valid_vrm):.4f}")
        print(f"    VRM > 0.1 (崎岖) 比例: {np.nansum(valid_vrm > 0.1) / len(valid_vrm) * 100:.1f}%")
    
    print(f"\n  VRM 主产出 (5×5) 用于第6章路径成本面")


# ============================================================
# Step 5: FoS 敏感性分析
# ============================================================
def step5_sensitivity_analysis():
    """
    5×5参数网格扫描
    C ∈ [0.5, 1.0, 1.5, 2.0, 3.0] kPa
    φ ∈ [30°, 32°, 35°, 38°, 40°]
    关键指标: F3>0.8区域面积变化率
    """
    print("\n" + "="*60)
    print("Step 5: FoS 敏感性分析 (5×5 = 25组)")
    print("="*60)
    
    slope, meta, nodata = read_raster(CONFIG['slope_path'])
    slope_rad = np.deg2rad(slope)
    
    rho = CONFIG['rho']
    g = CONFIG['g_moon']
    H = CONFIG['H_m']
    flat_thresh = np.deg2rad(CONFIG['flat_threshold_deg'])
    fos_min, fos_max = CONFIG['fos_range']
    pixel = CONFIG['pixel_size']
    
    flat_mask = slope < flat_thresh
    invalid_mask = np.isnan(slope)
    
    sin_cos = np.sin(slope_rad) * np.cos(slope_rad)
    sin_cos_safe = np.where(flat_mask | invalid_mask, 1.0, sin_cos)
    
    C_values = [500, 1000, 1500, 2000, 3000]  # Pa
    phi_values = [30, 32, 35, 38, 40]          # 度
    
    # 基准方案
    C_base, phi_base = 1500.0, np.deg2rad(35)
    num_base = C_base + rho * g * H * (np.cos(slope_rad)**2) * np.tan(phi_base)
    den_base = rho * g * H * sin_cos_safe
    fos_base = np.where(flat_mask, fos_max, np.clip(num_base / den_base, fos_min, fos_max))
    f3_base = (fos_base - fos_min) / (fos_max - fos_min)
    f3_base[invalid_mask] = np.nan
    base_area = np.nansum(f3_base > 0.8) * (pixel**2) / 1e6
    
    print(f"\n  基准方案 (C=1.5kPa, φ=35°): F3>0.8面积 = {base_area:.2f} km²")
    print(f"\n  {'C(kPa)':<8} {'φ(°)':<6} {'F3>0.8面积':<14} {'变化率':<10} {'最大差值'}")
    print(f"  {'─'*55}")
    
    results = []
    for C_val, phi_deg in product(C_values, phi_values):
        phi_rad = np.deg2rad(phi_deg)
        num = C_val + rho * g * H * (np.cos(slope_rad)**2) * np.tan(phi_rad)
        fos = np.where(flat_mask, fos_max, np.clip(num / den_base, fos_min, fos_max))
        f3 = (fos - fos_min) / (fos_max - fos_min)
        f3[invalid_mask] = np.nan
        
        area = np.nansum(f3 > 0.8) * (pixel**2) / 1e6
        change_pct = (area - base_area) / base_area * 100 if base_area > 0 else 0
        max_diff = np.nanmax(np.abs(f3 - f3_base))
        
        marker = " ← 基准" if (C_val == 1500 and phi_deg == 35) else ""
        print(f"  {C_val/1000:<8.1f} {phi_deg:<6} {area:<14.2f} {change_pct:>+8.1f}%   {max_diff:.3f}{marker}")
        
        results.append({
            'C': C_val, 'phi': phi_deg,
            'area': area, 'change_pct': change_pct,
            'max_diff': max_diff
        })
    
    # 结论
    max_change = max(abs(r['change_pct']) for r in results)
    print(f"\n  {'─'*55}")
    print(f"  参数扰动导致 F3>0.8 区域面积最大变化: {max_change:.1f}%")
    
    if max_change < 20:
        print("  结论稳健: 参数选择对F3空间分布影响有限")
    else:
        print("  结论对参数敏感: 需要在报告中重点讨论")
    
    # 保存报告图
    fig, ax = plt.subplots(figsize=(10, 6))
    C_arr = [r['C']/1000 for r in results]
    changes = [r['change_pct'] for r in results]
    phis = [r['phi'] for r in results]
    
    for phi in phi_values:
        subset = [(r['C']/1000, r['change_pct']) for r in results if r['phi'] == phi]
        cs, chs = zip(*subset)
        ax.plot(cs, chs, 'o-', label=f'φ={phi}°', linewidth=1.5, markersize=5)
    
    ax.axhline(y=0, color='green', linestyle='--', alpha=0.5, label='基准线')
    ax.axhline(y=20, color='red', linestyle=':', alpha=0.5, label='±20%阈值')
    ax.axhline(y=-20, color='red', linestyle=':', alpha=0.5)
    ax.set_xlabel('Cohesion C (kPa)')
    ax.set_ylabel('F3>0.8 Area Change (%)')
    ax.set_title('FoS Sensitivity Analysis: F3>0.8 Area Variation')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plot_path = os.path.join(CONFIG['output_dir'], 'fos_sensitivity_report.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  敏感性分析图已保存: {plot_path}")
    
    return results


# ============================================================
# 主流程
# ============================================================
def main():
    print("="*60)
    print("第4章 V3.2 综合执行脚本")
    print("LESF: 月面工程适宜性势场")
    print("="*60)
    print(f"开始时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 检查输入文件
    print("\n检查输入文件...")
    missing = []
    for key in ['slope_path', 'aspect_path', 'dem_path', 'curvature_path', 
                'avgvisib_path', 'avgvisib_earth_path', 'robins_path']:
        path = CONFIG[key]
        if not os.path.exists(path):
            missing.append(f"  {key}: {path}")
    
    if missing:
        print("  以下文件未找到:")
        for m in missing:
            print(m)
        print("\n  请修改脚本顶部 CONFIG 中的文件路径后重新运行。")
        print("  如果某些文件尚未生成（如 aspect.tif），可以先运行其他步骤。")
    
    # Step 0: ECSA 诊断（最关键，第一个运行）
    try:
        ecsa_decision, ecsa_stats = step0_ecsa_diagnostic()
    except Exception as e:
        print(f"\n  Step 0 失败: {e}")
        print("  请检查 AVGVISIB 和 AVGVISIB_EARTH 文件路径")
        ecsa_decision = 'UNKNOWN'
        ecsa_stats = {}
    
    # Step 1: d_max
    try:
        d_max = step1_dmax_determination()
    except Exception as e:
        print(f"\n  Step 1 失败: {e}")
        d_max = None
    
    # Step 2: 曲率阈值
    try:
        curv_threshold = step2_curvature_threshold()
    except Exception as e:
        print(f"\n  Step 2 失败: {e}")
        curv_threshold = None
    
    # Step 3: FoS → F3
    try:
        f3 = step3_fos_calculation()
    except Exception as e:
        print(f"\n  Step 3 失败: {e}")
    
    # Step 4: VRM
    try:
        step4_vrm_calculation()
    except Exception as e:
        print(f"\n  Step 4 失败: {e}")
    
    # Step 5: 敏感性分析
    try:
        sens_results = step5_sensitivity_analysis()
    except Exception as e:
        print(f"\n  Step 5 失败: {e}")
    
    # 汇总
    print("\n" + "="*60)
    print("执行汇总")
    print("="*60)
    print(f"  ECSA诊断结果: 路径{ecsa_decision}")
    if ecsa_stats:
        print(f"    AVGVISIB_EARTH CV = {ecsa_stats.get('comm_cv', 'N/A')}")
        print(f"    Pearson r = {ecsa_stats.get('r_pearson', 'N/A')}")
    if d_max:
        print(f"  d_max = {d_max:.0f} m")
    if curv_threshold:
        print(f"  曲率阈值 = {curv_threshold:.6f}")
    print(f"\n  后续手动步骤（GIS 栅格代数）:")
    print(f"    1. 连续危险距离场: 欧氏距离变换 → F5_continuous_distance.tif")
    print(f"    2. 硬性约束掩膜更新: distance < 240m → 排除")
    if ecsa_decision == 'B':
        print(f"    3. F4计算: P_sync = [AVGVISIB] × [AVGVISIB_EARTH]")
    elif ecsa_decision == 'A':
        print(f"    3. F4移除: 权重重新分配 F3=0.24, F5=0.18")
    elif ecsa_decision == 'C':
        print(f"    3. F4退回: 直接用AVGVISIB_EARTH作为F4")
    
    print(f"\n完成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == '__main__':
    main()
