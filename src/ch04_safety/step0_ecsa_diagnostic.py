# ──────────────────────────────────────────────────────────────────
# 章节  : Ch04 · 安全势场：ECSA 独立性诊断（光照 vs 地球可视：CV / Pearson / Spearman 三指标定路径）
# 来源  : 竞赛提交包 3工程文件/Ch04_*/step0_ecsa_diagnostic.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
Step 0: ECSA独立性诊断（独立运行版）
=====================================
月球南极科研站选址 · 第4章 V3.2 前置诊断
判断F4通信因子采用哪条路径：
  路径A: CV<5%  → F4移除，权重重分配
  路径B: CV>=5% 且 |r|<0.5 → ECSA乘积模型有效
  路径C: CV>=5% 且 |r|>=0.5 → 退回AVGVISIB_EARTH单独作F4

输入数据：
  (1) 历史开发机路径下的 gailv_240m.tif（AVGVISIB 光照概率，384x384）
      NASA月球光照概率(AVGVISIB), 384x384, 240m
      int16编码: 值/25000 = 概率[0,1]
      NoData=-32768, CRS元数据lat_0=+90(实际为南极空间)

  (2) 历史开发机路径下的 avgvisib_65s_240m_earth.jp2（对地可见概率，6420x6420）
      对地可见性概率(AVGVISIB_EARTH), 6420x6420, 240m
      int16编码: 值/25000 = 概率[0,1]
      NoData=None, CRS=Moon2000_spole lat_0=-90

空间对齐策略：
  gailv为参考网格(384x384, 92km x 92km)
  EARTH为大幅面(6420x6420, 1540km x 1540km)
  用rasterio窗口读取将EARTH裁剪到gailv范围
  两文件分辨率相同(240m)、坐标网格对齐，无需重采样

CRS说明：
  gailv的CRS元数据写的是北极立体投影(lat_0=+90)
  但坐标范围(-46080~46080)实际对应南极空间
  与EARTH的南极投影坐标等效，可直接逐像元比较

运行环境：Python 3.10+, 需要 numpy, rasterio, scipy, matplotlib, geopandas
作者：yuyi（队长）
日期：2026-07-15
"""

import numpy as np
import rasterio
from rasterio.windows import from_bounds
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import json

# ============================================================
# 文件路径配置（已根据实际数据位置修正）
# ============================================================
PATH_LIGHT = os.environ.get('LUNAR_LIGHT_TIF', os.path.join('data', 'rasters', 'AVGVISIB_probability.tif'))
PATH_EARTH = os.environ.get('LUNAR_EARTH_TIF', 'avgvisib_65s_240m_earth.jp2')
OUTPUT_DIR = os.environ.get('LUNAR_OUTPUT_DIR', os.path.join('data', 'output', 'ch04'))


def read_raster_raw(path):
    """
    读取栅格原始数据
    返回: (data_float64, profile, nodata_value)
    NoData像元替换为NaN
    """
    with rasterio.open(path) as src:
        data = src.read(1).astype(np.float64)
        profile = src.meta.copy()
        nodata = src.nodata
    if nodata is not None:
        data[data == nodata] = np.nan
    return data, profile, nodata


def align_earth_to_light(path_earth, ref_profile):
    """
    将EARTH数据裁剪到gailv(参考)的空间范围
    利用rasterio窗口读取，不重采样（两文件分辨率和网格对齐）

    参数:
      path_earth: EARTH jp2文件路径
      ref_profile: gailv的rasterio元数据(transform, crs等)
    返回:
      earth_data: 对齐后的float64数组(与gailv同shape), NoData→NaN
    """
    ref_transform = ref_profile['transform']
    ref_w = ref_profile['width']
    ref_h = ref_profile['height']
    # 从transform手动计算边界: left, bottom, right, top
    ref_left = ref_transform.c
    ref_top = ref_transform.f
    ref_right = ref_left + ref_transform.a * ref_w
    ref_bottom = ref_top + ref_transform.e * ref_h
    # ref_bounds = (left, bottom, right, top)
    ref_bounds = (ref_left, ref_bottom, ref_right, ref_top)

    with rasterio.open(path_earth) as src:
        # 计算EARTH中对应gailv范围的窗口
        window = from_bounds(
            ref_bounds[0], ref_bounds[1], ref_bounds[2], ref_bounds[3],
            src.transform
        )
        earth_data = src.read(window=window).astype(np.float64)[0]

        # EARTH文件NoData=None，但大幅面边缘可能有0值
        # 检查: 如果src.nodata有值则替换
        if src.nodata is not None:
            earth_data[earth_data == src.nodata] = np.nan

    # 安全检查: shape必须一致
    ref_shape = (ref_profile['height'], ref_profile['width'])
    if earth_data.shape != ref_shape:
        print(f"  [警告] 形状不匹配! gailv={ref_shape}, EARTH裁剪后={earth_data.shape}")
        print(f"  将自动裁剪/填充到一致...")
        h, w = ref_shape
        result = np.full(ref_shape, np.nan)
        eh = min(earth_data.shape[0], h)
        ew = min(earth_data.shape[1], w)
        result[:eh, :ew] = earth_data[:eh, :ew]
        earth_data = result

    return earth_data


def run_ecsa_diagnostic(light_raw, earth_raw):
    """
    ECSA独立性诊断核心逻辑

    编码转换: 两文件原始值 = 概率 x 25000, 需除以25000还原为概率
    无效值处理: NaN和<=0的像元不参与统计

    返回: (decision, stats_dict)
    """
    # ---- 编码转换: int16 -> 概率[0,1] ----
    light_prob = light_raw / 25000.0
    earth_prob = earth_raw / 25000.0

    # ---- 有效性掩膜 ----
    valid = (
        ~np.isnan(light_prob) & ~np.isnan(earth_prob)
        & (light_prob > 0) & (earth_prob > 0)
    )
    light_v = light_prob[valid]
    earth_v = earth_prob[valid]
    n_valid = len(light_v)

    print(f"\n  有效像元数: {n_valid:,} / {light_prob.size:,} "
          f"({n_valid / light_prob.size * 100:.1f}%)")

    if n_valid < 100:
        print("  [错误] 有效像元不足100个，无法诊断！")
        print("  请检查: (1)文件路径 (2)空间范围是否重叠 (3)NoData处理")
        return None, None

    # ---- 描述统计 ----
    light_mean = np.mean(light_v)
    light_std = np.std(light_v)
    light_cv = light_std / light_mean

    earth_mean = np.mean(earth_v)
    earth_std = np.std(earth_v)
    earth_cv = earth_std / earth_mean

    # ---- 相关性 ----
    r_pearson, p_pearson = stats.pearsonr(light_v, earth_v)
    r_spearman, p_spearman = stats.spearmanr(light_v, earth_v)

    # ---- 打印结果 ----
    print(f"\n  {'='*50}")
    print(f"  AVGVISIB (光照概率):")
    print(f"    均值   = {light_mean:.6f}")
    print(f"    标准差 = {light_std:.6f}")
    print(f"    变异系数 CV = {light_cv:.4f} ({light_cv*100:.2f}%)")
    print(f"    值域   = [{light_v.min():.6f}, {light_v.max():.6f}]")

    print(f"\n  AVGVISIB_EARTH (对地通信概率):")
    print(f"    均值   = {earth_mean:.6f}")
    print(f"    标准差 = {earth_std:.6f}")
    print(f"    变异系数 CV = {earth_cv:.4f} ({earth_cv*100:.2f}%)")
    print(f"    值域   = [{earth_v.min():.6f}, {earth_v.max():.6f}]")

    print(f"\n  空间相关性:")
    print(f"    Pearson  r = {r_pearson:.4f}  (p = {p_pearson:.4e})")
    print(f"    Spearman ρ = {r_spearman:.4f}  (p = {p_spearman:.4e})")

    # ---- 决策路径 ----
    print(f"\n  {'='*50}")
    print(f"  决策判定:")
    print(f"  {'─'*50}")

    if earth_cv < 0.05:
        decision = 'A'
        print(f"  >>> 路径A: AVGVISIB_EARTH 空间变异极小 (CV = {earth_cv*100:.2f}% < 5%)")
        print(f"  >>> F4因子无区分度，建议从AHP中移除")
        print(f"  >>> 原F4权重0.12分配给F3(+0.06)和F5(+0.06)")
        print(f"  >>> 最终权重: F1=0.28 / F2=0.20 / F3=0.24 / F5=0.18 / F6=0.10")
    elif abs(r_pearson) >= 0.5:
        decision = 'C'
        print(f"  >>> 路径C: 光照与通信显著相关 (|r| = {abs(r_pearson):.4f} >= 0.5)")
        print(f"  >>> 独立性假设不严格成立，ECSA乘积模型不可用")
        print(f"  >>> 退回AVGVISIB_EARTH单独作为F4因子")
        print(f"  >>> 权重保持: F1=0.28 / F2=0.20 / F3=0.18 / F4=0.12 / F5=0.12 / F6=0.10")
    else:
        decision = 'B'
        print(f"  >>> 路径B: 独立性假设基本成立 (CV = {earth_cv*100:.2f}%, |r| = {abs(r_pearson):.4f})")
        print(f"  >>> ECSA乘积模型可用: P_sync = AVGVISIB x AVGVISIB_EARTH")
        print(f"  >>> 权重保持: F1=0.28 / F2=0.20 / F3=0.18 / F4=0.12 / F5=0.12 / F6=0.10")

    print(f"  {'─'*50}")

    # 同步可用性统计
    psync = light_prob[valid] * earth_prob[valid]
    print(f"\n  同步可用性 P_sync 统计:")
    print(f"    均值 = {np.mean(psync):.6f}")
    print(f"    中位数 = {np.median(psync):.6f}")
    print(f"    P_sync > 0.5 面积: "
          f"{np.sum(psync > 0.5) * 240**2 / 1e6:.2f} km2")
    print(f"    P_sync > 0.8 面积: "
          f"{np.sum(psync > 0.8) * 240**2 / 1e6:.2f} km2")

    result = {
        'decision': decision,
        'n_valid': int(n_valid),
        'light_mean': float(light_mean),
        'light_cv': float(light_cv),
        'earth_mean': float(earth_mean),
        'earth_cv': float(earth_cv),
        'r_pearson': float(r_pearson),
        'p_pearson': float(p_pearson),
        'r_spearman': float(r_spearman),
        'p_spearman': float(p_spearman),
        'psync_mean': float(np.mean(psync)),
        'psync_gt05_area_km2': float(np.sum(psync > 0.5) * 240**2 / 1e6),
        'psync_gt08_area_km2': float(np.sum(psync > 0.8) * 240**2 / 1e6),
    }

    return decision, result


def save_report(decision, result, light_raw, earth_raw):
    """保存诊断报告(文本+JSON)和可视化图"""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ---- 文本报告 ----
    report_path = os.path.join(OUTPUT_DIR, 'step0_ecsa_diagnostic_report.txt')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("ECSA独立性诊断报告\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"决策路径: {decision}\n\n")
        for k, v in result.items():
            f.write(f"  {k}: {v}\n")
        f.write("\n" + "=" * 60 + "\n")
        f.write("决策说明:\n")
        if decision == 'A':
            f.write("  路径A - F4移除。AVGVISIB_EARTH空间变异极小(CV<5%)，\n")
            f.write("  通信覆盖度在研究区内几乎无差异，作为因子层无区分度。\n")
            f.write("  建议将F4权重(0.12)分配给F3(+0.06)和F5(+0.06)。\n")
        elif decision == 'B':
            f.write("  路径B - ECSA乘积模型有效。光照与通信空间独立性成立，\n")
            f.write("  可用 P_sync = AVGVISIB x AVGVISIB_EARTH 计算同步可用性。\n")
        else:
            f.write("  路径C - 退回AVGVISIB_EARTH单独作F4。光照与通信显著相关，\n")
            f.write("  独立性假设不成立，ECSA乘积模型不可用。\n")
    print(f"\n  报告已保存: {report_path}")

    # ---- JSON结果(方便后续脚本读取) ----
    json_path = os.path.join(OUTPUT_DIR, 'step0_ecsa_result.json')
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    print(f"  JSON已保存: {json_path}")

    # ---- 可视化 ----
    light_prob = light_raw / 25000.0
    earth_prob = earth_raw / 25000.0
    valid = (
        ~np.isnan(light_prob) & ~np.isnan(earth_prob)
        & (light_prob > 0) & (earth_prob > 0)
    )
    light_v = light_prob[valid]
    earth_v = earth_prob[valid]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # 左: 光照概率直方图
    ax = axes[0]
    ax.hist(light_v, bins=50, color='#e74c3c', alpha=0.8, edgecolor='white', linewidth=0.5)
    ax.axvline(x=np.mean(light_v), color='black', linestyle='--', linewidth=1.5,
               label=f'Mean={np.mean(light_v):.3f}')
    ax.set_xlabel('Light Probability (AVGVISIB / 25000)')
    ax.set_ylabel('Pixel Count')
    ax.set_title(f'AVGVISIB Distribution (n={len(light_v):,})')
    ax.legend()

    # 中: 通信概率直方图
    ax = axes[1]
    ax.hist(earth_v, bins=50, color='#3498db', alpha=0.8, edgecolor='white', linewidth=0.5)
    ax.axvline(x=np.mean(earth_v), color='black', linestyle='--', linewidth=1.5,
               label=f'Mean={np.mean(earth_v):.3f}')
    ax.set_xlabel('Earth-Comm Probability (AVGVISIB_EARTH / 25000)')
    ax.set_ylabel('Pixel Count')
    ax.set_title(f'AVGVISIB_EARTH Distribution (CV={np.std(earth_v)/np.mean(earth_v)*100:.1f}%)')
    ax.legend()

    # 右: 散点图(抽样)
    ax = axes[2]
    max_pts = min(20000, len(light_v))
    idx = np.random.choice(len(light_v), max_pts, replace=False)
    ax.scatter(light_v[idx], earth_v[idx], s=2, alpha=0.3, c='#2c3e50')
    ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, linewidth=1, label='y=x')
    r, p = stats.pearsonr(light_v, earth_v)
    ax.text(0.05, 0.95, f'Pearson r = {r:.4f}\np = {p:.2e}\nn = {len(light_v):,}',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.set_xlabel('Light Probability')
    ax.set_ylabel('Earth-Comm Probability')
    ax.set_title('Light vs Earth-Comm Correlation')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc='lower right')

    plt.tight_layout()
    fig_path = os.path.join(OUTPUT_DIR, 'step0_ecsa_diagnostic_plot.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  可视化已保存: {fig_path}")


# ============================================================
# 主程序
# ============================================================
def main():
    print()
    print("=" * 60)
    print("  Step 0: ECSA 独立性诊断")
    print("  第4章 V3.2 前置分析")
    print("=" * 60)

    # 1. 读取gailv(参考网格)
    print(f"\n  [1/4] 读取光照数据...")
    print(f"    路径: {PATH_LIGHT}")
    light_raw, light_meta, light_nodata = read_raster_raw(PATH_LIGHT)
    print(f"    形状: {light_raw.shape}")
    print(f"    CRS: {light_meta['crs']}")
    print(f"    NoData: {light_nodata}")
    print(f"    原始值域: [{np.nanmin(light_raw):.0f}, {np.nanmax(light_raw):.0f}]")
    print(f"    有效像元: {np.sum(~np.isnan(light_raw)):,}")

    # 2. 读取EARTH并裁剪到gailv范围
    print(f"\n  [2/4] 读取对地通信数据并对齐空间范围...")
    print(f"    路径: {PATH_EARTH}")
    earth_raw = align_earth_to_light(PATH_EARTH, light_meta)
    print(f"    对齐后形状: {earth_raw.shape}")
    print(f"    原始值域: [{np.nanmin(earth_raw):.0f}, {np.nanmax(earth_raw):.0f}]")
    print(f"    有效像元: {np.sum(~np.isnan(earth_raw)):,}")

    # 3. 运行诊断
    print(f"\n  [3/4] 运行ECSA独立性诊断...")
    decision, result = run_ecsa_diagnostic(light_raw, earth_raw)

    if decision is None:
        print("\n  [失败] 诊断未能完成，请检查数据！")
        return

    # 4. 保存结果
    print(f"\n  [4/4] 保存诊断结果...")
    save_report(decision, result, light_raw, earth_raw)

    # 最终总结
    print(f"\n{'='*60}")
    print(f"  诊断完成!")
    print(f"  决策路径: {decision}")
    print(f"  通信CV: {result['earth_cv']*100:.2f}%")
    print(f"  Pearson r: {result['r_pearson']:.4f}")
    print(f"  Spearman ρ: {result['r_spearman']:.4f}")

    if decision == 'A':
        print(f"\n  >>> 下一步: 从AHP中移除F4，权重重分配给F3和F5")
        print(f"  >>> 新权重: F1=0.28 / F2=0.20 / F3=0.24 / F5=0.18 / F6=0.10")
    elif decision == 'B':
        print(f"\n  >>> 下一步: 计算ECSA = AVGVISIB x AVGVISIB_EARTH 作为F4")
        print(f"  >>> 权重不变: F1=0.28 / F2=0.20 / F3=0.18 / F4=0.12 / F5=0.12 / F6=0.10")
    else:
        print(f"\n  >>> 下一步: 用AVGVISIB_EARTH单独作为F4(不做乘积)")
        print(f"  >>> 权重不变: F1=0.28 / F2=0.20 / F3=0.18 / F4=0.12 / F5=0.12 / F6=0.10")

    print(f"{'='*60}")


if __name__ == '__main__':
    main()
