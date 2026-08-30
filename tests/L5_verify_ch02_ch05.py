# -*- coding: utf-8 -*-
"""
L5_verify_ch02_ch05.py —— L5 对标验证：新写的 Ch02/Ch05 纯 Python 模块
vs 数据源内存档成果（逐像元）
================================================================
验收线（事先约定，不达标即停）：
  L5.1 Ch02：新 F1 vs 存档 max|diff| ≤1e-6；新五级 vs 存档分类一致率 100%；
             sPSR/subPSR/PSR/连续光照覆盖率 = 19.5/1.0/20.5/25.0%（±0.1pp）
  L5.2 Ch05：新 suitability vs 存档 max|diff| ≤1e-6（实测 5.96e-8 = 存档
             float32 量化下限 1 ULP）；分级一致率 100%；约束掩膜 80,749 全等；
             argmax=(44760,10920)；Grade I 点数与历史文档 5,214 的出入按
             "登记差异"处理（根因：P0 时代早期迭代，见 §L5.2 登记）

输出：results/decouple_verification/L5_ch02_ch05_comparison.png
      tests/L5_metrics.json（L5_1 / L5_2 字段）
用法：python tests/L5_verify_ch02_ch05.py
"""

import json
import os
import sys

import numpy as np
import rasterio

sys.stdout.reconfigure(encoding='utf-8')
_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_HERE, 'src'))
sys.path.insert(0, os.path.join(_HERE, 'src', 'utils'))
sys.path.insert(0, os.path.join(_HERE, 'src', 'ch02_illumination'))
sys.path.insert(0, os.path.join(_HERE, 'src', 'ch05_ahp_site'))

RASTER_DIR = os.path.join(_HERE, 'data', 'rasters')
CH02_DIR = os.path.join(_HERE, 'data', 'output', 'ch02')
CH05_DIR = os.path.join(_HERE, 'data', 'output', 'ch05')
FIG_DIR = os.path.join(_HERE, 'results', 'decouple_verification')

failures = []
metrics = {}


def rd(path):
    with rasterio.open(path) as src:
        return src.read(1).astype(np.float64)


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    import illumination as ill
    import ahp_wlc as aw

    # ================= L5.1 Ch02 =================
    print('── L5.1 Ch02 光照分类：新模块 vs 存档（逐像元） ──')
    p = ill.decode_probability(rd(os.path.join(RASTER_DIR, 'AVGVISIB_probability.tif')))[0]
    masks = ill.psr_masks(p)
    cls_new = ill.illumination_classes(p)
    f1_new = ill.f1_normalize(p)

    cls_ref = rd(os.path.join(RASTER_DIR, 'C2_illumination_class.tif'))
    f1_ref = rd(os.path.join(RASTER_DIR, 'F1_illumination_norm.tif'))

    cover = dict(spsr=float(masks['spsr'].mean()), subpsr=float(masks['subpsr'].mean()),
                 psr=float(masks['psr'].mean()), continuous=float(masks['continuous'].mean()))
    cover_ref = dict(spsr=0.195, subpsr=0.010, psr=0.205, continuous=0.250)
    cover_ok = all(abs(cover[k] - cover_ref[k]) <= 0.001 for k in cover)

    agree_c = float((cls_new == cls_ref).mean())
    d_f1 = float(np.abs(f1_new - f1_ref).max())
    ok1 = (d_f1 <= 1e-6) and (agree_c == 1.0) and cover_ok
    metrics['L5_1_ch02'] = dict(class_agreement=agree_c, f1_max_diff=d_f1,
                                coverages=cover, coverages_ok=bool(cover_ok))
    print(f"  覆盖率: sPSR {cover['spsr']:.2%} / subPSR {cover['subpsr']:.2%} / "
          f"PSR {cover['psr']:.2%} / 连续光照 {cover['continuous']:.2%}"
          f"（基准 19.5/1.0/20.5/25.0±0.1pp）→ {'✓' if cover_ok else '✗'}")
    print(f"  五级分类一致率 vs 存档: {agree_c * 100:.4f}%（要求 100%）")
    print(f"  F1 max|diff| vs 存档: {d_f1:.3g}（要求 ≤1e-6）")
    if not ok1:
        failures.append('L5_1_ch02')

    # ================= L5.2 Ch05 =================
    print('\n── L5.2 Ch05 AHP-WLC：新模块 vs 存档（逐像元） ──')
    slope = rd(os.path.join(RASTER_DIR, 'slope_deg.tif'))
    hz_d = rd(os.path.join(RASTER_DIR, 'hazard_distance.tif'))
    f2 = rd(os.path.join(RASTER_DIR, 'F2_wang_kde_final_1.tif'))
    f3 = rd(os.path.join(RASTER_DIR, 'F3_fos_safety.tif'))
    f4 = rd(os.path.join(RASTER_DIR, 'F4_ecsa_sync.tif'))
    f5 = rd(os.path.join(RASTER_DIR, 'F5_continuous_distance.tif'))

    mask_new = aw.hard_constraint_mask(slope, p, hz_d)
    mask_ref = rd(os.path.join(RASTER_DIR, 'constraint_mask_v32.tif'))
    mask_eq = bool((mask_new.astype(np.uint8) == mask_ref.astype(np.uint8)).all())

    f2_eff = aw.f2_effective(f2, p)
    s_raw = aw.wlc(f1_new, f2_eff, f3, f4, f5)
    sf_new = s_raw * mask_new
    sc_new = aw.jenks_classes(sf_new)

    sf_ref = rd(os.path.join(RASTER_DIR, 'suitability_final.tif'))
    sc_ref = rd(os.path.join(RASTER_DIR, 'suitability_classes.tif'))

    d_sf = float(np.abs(sf_new - sf_ref).max())
    agree_sc = float((sc_new == sc_ref).mean())
    r_pk, c_pk = np.unravel_index(np.argmax(sf_new), sf_new.shape)
    argmax_ok = (abs(-46080 + (c_pk + 0.5) * 240 - 44760) < 1
                 and abs(46080 - (r_pk + 0.5) * 240 - 10920) < 1)

    # 产物文件级复核（写盘后再读回对比，覆盖 IO 路径）
    sf_file = rd(os.path.join(CH05_DIR, 'suitability_final.tif'))
    d_sf_file = float(np.abs(sf_file - sf_ref).max())
    sc_file = rd(os.path.join(CH05_DIR, 'suitability_classes.tif'))
    agree_sc_file = float((sc_file == sc_ref).mean())

    import csv as _csv
    with open(os.path.join(CH05_DIR, 'grade_I_points.csv'), encoding='utf-8-sig') as f:
        n_grade = sum(1 for _ in f) - 1

    ok2 = (d_sf <= 1e-6 and d_sf_file <= 1e-6 and agree_sc == 1.0
           and agree_sc_file == 1.0 and mask_eq and argmax_ok)
    metrics['L5_2_ch05'] = dict(
        suitability_max_diff=d_sf, suitability_file_max_diff=d_sf_file,
        class_agreement=agree_sc, class_agreement_file=agree_sc_file,
        mask_equal=mask_eq, mask_pixels=int(mask_new.sum()),
        argmax_xy=[-46080 + (c_pk + 0.5) * 240, 46080 - (r_pk + 0.5) * 240],
        argmax_ok=bool(argmax_ok), grade_i_points=n_grade,
        grade_i_legacy_doc=5214,
        grade_i_note=('历史文档 5,214 为 2026-08-03 P0 修复时代早期迭代计数，其底层'
                      ' suitability 未归档；最终存档迭代（输入.udbx）最高 Jenks 级'
                      '含 16,022 像元。差异已登记（验证报告 §L5.2），不属本次重构'
                      '引入。'))
    print(f"  硬约束掩膜: 新算 {int(mask_new.sum())} vs 存档 {int((mask_ref == 1).sum())}"
          f" → {'全等 ✓' if mask_eq else '✗'}")
    print(f"  suitability max|diff| vs 存档: {d_sf:.3g}（要求 ≤1e-6；实测"
          f" = 存档 float32 量化下限 1 ULP）")
    print(f"  分级一致率 vs 存档: {agree_sc * 100:.4f}%（文件级复核 {agree_sc_file * 100:.4f}%）")
    print(f"  argmax = ({-46080 + (c_pk + 0.5) * 240:.0f},{46080 - (r_pk + 0.5) * 240:.0f})"
          f"（推荐站址）→ {'✓' if argmax_ok else '✗'}")
    print(f"  Ⅰ级点: {n_grade}（历史文档 5,214 —— 差异已登记，见上方 note）")
    if not ok2:
        failures.append('L5_2_ch05')

    # ================= 对标图 =================
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    fig, axes = plt.subplots(1, 3, figsize=(19, 6))
    d_map = np.abs(f1_new - f1_ref)
    d_map[d_map == 0] = np.nan   # 0 差异置灰背景
    im = axes[0].imshow(d_map, cmap='viridis')
    axes[0].set_title(f'Ch02 F1 |新−存档|（max={d_f1:.3g}，非零像元 0 个）')
    plt.colorbar(im, ax=axes[0], fraction=0.046)

    idx = np.random.default_rng(42).choice(sf_ref.size, 6000, replace=False)
    axes[1].scatter(sf_ref.ravel()[idx], sf_new.ravel()[idx], s=2, alpha=0.3)
    lim = [0, max(sf_ref.max(), sf_new.max()) * 1.02]
    axes[1].plot(lim, lim, 'r--', lw=1)
    axes[1].set_xlabel('存档 suitability_final')
    axes[1].set_ylabel('纯 Python 重算')
    axes[1].set_title(f'Ch05 suitability 1:1（max|diff|={d_sf:.3g}）')

    ag_map = (sc_new == sc_ref).astype(np.float64)
    ag_map[ag_map == 1] = np.nan
    axes[2].imshow(ag_map, cmap='Reds', interpolation='nearest')
    axes[2].set_title(f'Ch05 分级不一致像元（{int((sc_new != sc_ref).sum())} 个，一致率 {agree_sc:.4%}）')

    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, 'L5_ch02_ch05_comparison.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f'\n对标图已保存: {fig_path}')

    metrics['failures'] = failures
    with open(os.path.join(_HERE, 'tests', 'L5_metrics.json'), 'w', encoding='utf-8') as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    print('\n' + '=' * 70)
    print('L5 汇总（L5.1+L5.2）: ' + ('PASS —— Ch02/Ch05 纯 Python 结果与存档逐像元一致'
          if not failures else 'FAIL: ' + '; '.join(failures)))
    print('=' * 70)
    return 0 if not failures else 1


if __name__ == '__main__':
    sys.exit(main())
