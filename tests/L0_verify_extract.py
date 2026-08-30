# -*- coding: utf-8 -*-
"""
L0_verify_extract.py —— L0 数据层验证：UDBX → GeoTIFF 提取正确性
================================================================
验证内容：
  1. 每个导出 TIFF 的网格规格：384×384、分辨率 240m、范围 ±46080、CRS 为月球南极立体投影；
  2. 数值级回环：从 UDBX 重新独立解码同一栅格（走一遍完整解码管线），与导出 TIFF
     逐像元比较——完全相等像元占比必须 = 100%（提取只是"解压搬运+NoData 归一"，
     未做任何重采样，浮点允许用 == 精确比较）；
  3. 打印每个栅格的 min/max/mean/有效像元数留档。

容差说明：不设任何数值容差——搬运过程不改变数值，必须 100% 一致。
NoData 处理：UDBX 内部填充值（-32768/-9999/255/65535 等）按注册表 SmNovalue 识别，
导出统一映射为 -9999；回环比较时同样先做该映射再比较。

用法：python tests/L0_verify_extract.py
环境变量：LUNAR_UDBX（默认 输入.udbx）、LUNAR_RASTER_DIR（默认 data/rasters）
"""

import os
import sys

import numpy as np
import rasterio

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
from utils.udbx_extract import (DEFAULT_EXPORT, GRID_N, RES, LEFT, TOP,
                               extract_raster, list_rasters)

UDBX = os.environ.get('LUNAR_UDBX', os.path.join(
    '..', 'A+ak9HzkqIZG', '2成果数据', '数据源', '输入.udbx'))
RASTER_DIR = os.environ.get('LUNAR_RASTER_DIR', os.path.join('data', 'rasters'))


def check_grid(src):
    """网格规格检查：尺寸/分辨率/范围/CRS。"""
    issues = []
    if src.width != GRID_N or src.height != GRID_N:
        issues.append(f"尺寸 {src.width}x{src.height} ≠ {GRID_N}x{GRID_N}")
    rx, ry = src.res
    if abs(rx - RES) > 1e-6 or abs(ry - RES) > 1e-6:
        issues.append(f"分辨率 ({rx},{ry}) ≠ 240")
    b = src.bounds
    for got, want in [(b.left, LEFT), (b.top, TOP), (b.right, -LEFT), (b.bottom, -TOP)]:
        if abs(got - want) > 1e-6:
            issues.append(f"bounds 分量 {got} ≠ {want}")
            break
    if '1737400' not in (src.crs.to_proj4() if src.crs else ''):
        issues.append(f"CRS 异常: {src.crs}")
    return issues


def main():
    # 只验证 UDBX 中实际存在的数据集
    present = {i['name'] for i in list_rasters(UDBX)}
    targets = [(ds, tif) for ds, tif in DEFAULT_EXPORT.items() if ds in present]
    tif_path = lambda t: os.path.join(RASTER_DIR, t + '.tif')
    targets = [(ds, t) for ds, t in targets if os.path.exists(tif_path(t))]

    print('=' * 78)
    print('L0 数据层验证：UDBX → GeoTIFF 提取回环（基准：铁基准 384×384@240m@±46080）')
    print('=' * 78)
    print(f'UDBX: {UDBX}')
    print(f'TIFF 目录: {RASTER_DIR}（待验证 {len(targets)} 个）\n')

    n_pass, n_fail = 0, 0
    rows = []
    for ds_name, tif_name in targets:
        path = tif_path(tif_name)
        with rasterio.open(path) as src:
            tif_arr = src.read(1)
            issues = check_grid(src)
        # 从 UDBX 独立重新解码
        udbx_arr, meta = extract_raster(UDBX, ds_name)
        nv = meta['novalue']
        if nv is not None:
            nod = np.isclose(udbx_arr, float(nv)) | np.isnan(udbx_arr)
        else:
            nod = np.isnan(udbx_arr)
        udbx_norm = np.where(nod, -9999.0, udbx_arr).astype(np.float32)

        equal = (tif_arr == udbx_norm)
        pct = equal.mean() * 100
        valid = ~nod
        vmin = float(udbx_arr[valid].min()) if valid.any() else float('nan')
        vmax = float(udbx_arr[valid].max()) if valid.any() else float('nan')
        vmean = float(udbx_arr[valid].mean()) if valid.any() else float('nan')

        ok = (not issues) and pct == 100.0
        n_pass += ok
        n_fail += not ok
        tag = 'PASS' if ok else 'FAIL'
        issue_str = ('; '.join(issues)) if issues else '-'
        print(f"[{tag}] {tif_name:32s} 一致像元 {pct:6.2f}%  "
              f"有效 {int(valid.sum()):6d}  min={vmin:12.4f} max={vmax:12.4f} mean={vmean:10.4f}"
              + (f"  问题: {issue_str}" if issues else ''))
        rows.append((tag, tif_name, pct, int(valid.sum())))

    print('\n' + '=' * 78)
    print(f'L0 汇总: {n_pass} PASS / {n_fail} FAIL（共 {len(targets)} 个栅格）')
    print('容差声明: 数值级回环要求 100% 逐像元相等（无重采样、无坐标变换，仅解压搬运）')
    print('=' * 78)
    return 0 if n_fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
