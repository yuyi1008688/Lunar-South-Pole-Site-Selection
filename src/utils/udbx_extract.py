# -*- coding: utf-8 -*-
"""
udbx_extract.py —— UDBX（SQLite 格式）→ GeoTIFF 纯 Python 数据导出器
================================================================
UDBX 本质是 SQLite 数据库（可直接 sqlite3.connect，无需任何 GIS 软件）。
本模块把 UDBX 内的栅格数据集"解放"为标准 GeoTIFF，供后续纯 Python 分析链使用。

UDBX 栅格存储结构（逆向确认，详见 tests/解耦与精度验证报告.md 附录）：
  SmImgRegister   每个栅格数据集一行：SmDatasetName（数据集名）、SmTableName（真实分块表名，
                  可能被截断到 ~20 字符，必须用它查表）、SmWidth/SmHeight、
                  SmGeoLeft/Top/Right/Bottom（地理范围）、SmeBlockSize（分块边长，256 或 128）
  SmBandRegister  每个波段一行：SmPixelFormat（像素格式）、SmNovalue（NoData 值）、
                  SmMaxBlockSize（块最大字节数 = 块像素数×每像素字节 + 4）
  <SmTableName>   分块数据表：SmRow/SmColumn（块网格坐标，从左上角 (0,0) 起）、
                  SmBand（块字节流，zlib 压缩或原始）

像素格式 → 每像素字节数（以 SmMaxBlockSize 实测为准，枚举值仅作参考）：
  3200 → 4 字节 float32；6400 → 8 字节 float64；8 → 1 字节 uint8；
  160 → 2 字节 uint16；4 → 0.5 字节（UBIT4 半字节打包，高半字节在前）；
  1 → 0.125 字节（UBIT1 位打包，每字节 8 像素，MSB 在前）

铁基准：384×384 @ 240m @ ±46080m，Moon South Polar Stereographic，NoData 统一 -9999。
凡是导出结果出现 383 尺寸或 ±46000 范围都视为错误并报警（历史教训：范围不整除时
GIS 软件会偷偷修改分辨率，商业 GIS 软件自带导出功能产出的 tif 就是错的，UDBX 内部才是
source of truth）。

用法：
  python src/utils/udbx_extract.py --udbx <路径.udbx> [--out data/rasters] [--only 名称1,名称2]

环境变量：LUNAR_UDBX（默认 udbx 路径）、LUNAR_RASTER_DIR（默认输出目录）。
"""

import argparse
import os
import sqlite3
import sys
import zlib

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin

sys.stdout.reconfigure(encoding='utf-8')

# ── 铁基准常量（与项目全局约定一致） ─────────────────────────────
GRID_N = 384
RES = 240.0
LEFT, TOP = -46080.0, 46080.0
CRS_PROJ4 = ('+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +x_0=0 +y_0=0 '
             '+R=1737400 +units=m +no_defs')
NODATA_OUT = -9999.0

# 数据集名 → 输出 GeoTIFF 名（Ch06 必需 9 个 + 建议一并提取的全流程因子层）
DEFAULT_EXPORT = {
    # Ch06 重写必需
    'Slope':                 'slope_deg',
    'AVGVISIB_probability':  'AVGVISIB_probability',
    'combined_hazard':       'combined_hazard',
    'cost_surface':          'cost_surface',
    'DEM_fused':             'DEM_fused',
    'vrm_5x5_float32':       'vrm_5x5',          # 优先 float32 版（vrm_5x5 原表是 float64）
    'distance_accumulation': 'distance_accumulation',
    'backlink_direction':    'backlink_direction',
    'optimal_path_raw':      'optimal_path_raw',
    # 全流程因子层（论文复现/出图用）
    'F1_illumination_norm':      'F1_illumination_norm',
    'F3_fos_safety':             'F3_fos_safety',
    'F4_ecsa_sync':              'F4_ecsa_sync',
    'F5_continuous_distance':    'F5_continuous_distance',
    'suitability_final':         'suitability_final',
    'suitability_classes':       'suitability_classes',
    'constraint_mask_v32':       'constraint_mask_v32',
    'hard_constraint_v32':       'hard_constraint_v32',
    'sPSR_mask':                 'sPSR_mask',
    'subPSR_mask':               'subPSR_mask',
    'C2_PSR_mask':               'C2_PSR_mask',
    'C2_illumination_class':     'C2_illumination_class',
    'C2_continuous_light_mask':  'C2_continuous_light_mask',
    'FoS_raw':                   'FoS_raw',
    'FoS_fixed':                 'FoS_fixed',
    'ice_density_final':         'ice_density_final',
    'Aspect':                    'aspect',
    'result_profileCurvature':   'profile_curvature',
    'result_planCurvature':      'plan_curvature',
    'result_averageCurvature':   'average_curvature',
    'hazard_distance':           'hazard_distance',
    'C_slope':                   'C_slope',
    'C_hazard':                  'C_hazard',
    'C_vrm':                     'C_vrm',
    'LOLA_240m_dem':             'LOLA_240m_dem',
    'CE2':                       'CE2',
}

# UDBX PixelFormat 枚举 → 每像素字节数（兜底表；实际以 SmMaxBlockSize 实测为准）
PIXEL_FORMAT_BYTES = {1: 1/8, 4: 1/2, 8: 1.0, 16: 2.0, 32: 4.0,
                      160: 2.0, 3200: 4.0, 6400: 8.0}


def list_rasters(udbx_path):
    """枚举 SmImgRegister 全部栅格数据集，返回信息字典列表。"""
    conn = sqlite3.connect(udbx_path)
    cur = conn.cursor()
    cur.execute(
        "SELECT SmDatasetID, SmDatasetName, SmTableName, SmWidth, SmHeight, "
        "SmeBlockSize, SmGeoLeft, SmGeoTop, SmGeoRight, SmGeoBottom "
        "FROM SmImgRegister ORDER BY SmDatasetID")
    imgs = []
    for did, name, tbl, w, h, bs, l, t, r, b in cur.fetchall():
        imgs.append(dict(dataset_id=did, name=name, table=tbl, width=w, height=h,
                         block_size=bs, left=l, top=t, right=r, bottom=b))
    conn.close()
    return imgs


def _pixel_bytes(cur, dataset_id, block_size):
    """从 SmBandRegister 推断每像素字节数（MaxBlockSize=块像素×字节+4）。"""
    cur.execute("SELECT SmPixelFormat, SmNovalue, SmMaxBlockSize "
                "FROM SmBandRegister WHERE SmDatasetID=?", (dataset_id,))
    row = cur.fetchone()
    if row is None:
        raise RuntimeError(f"SmBandRegister 中找不到 DatasetID={dataset_id}")
    pixfmt, novalue, maxblk = row
    n_block_px = block_size * block_size
    if maxblk:
        bpp = (maxblk - 4) / n_block_px
        if bpp <= 0 or bpp > 8:
            bpp = PIXEL_FORMAT_BYTES.get(pixfmt)
    else:
        bpp = PIXEL_FORMAT_BYTES.get(pixfmt)
    if bpp is None:
        raise RuntimeError(f"未知像素格式 SmPixelFormat={pixfmt}")
    return float(bpp), pixfmt, novalue


def decode_blob(blob, bpp, n_pixels):
    """块字节流 → 一维像素值数组（float64）。zlib 压缩自动检测+兜底。"""
    if blob[:1] == b'\x78':
        try:
            blob = zlib.decompress(blob)
        except zlib.error:
            pass  # 魔数巧合，按原始字节处理
    if bpp == 4.0:
        vals = np.frombuffer(blob, dtype='<f4').astype(np.float64)
    elif bpp == 8.0:
        vals = np.frombuffer(blob, dtype='<f8').astype(np.float64)
    elif bpp == 1.0:
        vals = np.frombuffer(blob, dtype=np.uint8).astype(np.float64)
    elif bpp == 2.0:
        vals = np.frombuffer(blob, dtype='<u2').astype(np.float64)
    elif bpp == 0.5:   # UBIT4：每字节 2 个半字节，高半字节在前
        b = np.frombuffer(blob, dtype=np.uint8).astype(np.uint8)
        vals = np.empty(b.size * 2, dtype=np.float64)
        vals[0::2] = (b >> 4) & 0x0F
        vals[1::2] = b & 0x0F
    elif bpp == 0.125:  # UBIT1：每字节 8 位，MSB 在前
        b = np.frombuffer(blob, dtype=np.uint8)
        vals = np.unpackbits(b).astype(np.float64)
    else:
        raise RuntimeError(f"不支持的每像素字节数：{bpp}")
    return vals[:n_pixels]


def extract_raster(udbx_path, dataset_name):
    """提取单个栅格为 (arr, meta)。arr 形状 (H, W)，第 0 行=地理最北。"""
    conn = sqlite3.connect(udbx_path)
    cur = conn.cursor()
    cur.execute("SELECT SmDatasetID, SmTableName, SmWidth, SmHeight, SmeBlockSize, "
                "SmGeoLeft, SmGeoTop, SmGeoRight, SmGeoBottom "
                "FROM SmImgRegister WHERE SmDatasetName=?", (dataset_name,))
    row = cur.fetchone()
    if row is None:
        conn.close()
        raise RuntimeError(f"UDBX 中不存在数据集：{dataset_name}")
    did, tbl, w, h, bs, l, t, r_, b_ = row
    bs = bs or 256
    bpp, pixfmt, novalue = _pixel_bytes(cur, did, bs)

    arr = np.full((h, w), np.nan, dtype=np.float64)
    cur.execute(f"SELECT SmRow, SmColumn, SmBand FROM [{tbl}]")
    n_block_px = bs * bs
    for br, bc, blob in cur.fetchall():
        vals = decode_blob(blob, bpp, n_block_px)
        if vals.size < n_block_px:  # 理论上块都是满尺寸；不足则跳过
            continue
        block = vals.reshape(bs, bs)
        r0, c0 = br * bs, bc * bs
        r1, c1 = min(r0 + bs, h), min(c0 + bs, w)
        if r0 >= h or c0 >= w:
            continue
        arr[r0:r1, c0:c1] = block[:r1 - r0, :c1 - c0]
    conn.close()

    meta = dict(width=w, height=h, block_size=bs, pixfmt=pixfmt,
                novalue=novalue, left=l, top=t, right=r_, bottom=b_)
    return arr, meta


def export_raster(udbx_path, dataset_name, out_tif, check=True):
    """提取并导出单个 GeoTIFF（float32, NoData=-9999, deflate 压缩），返回自检摘要。"""
    arr, meta = extract_raster(udbx_path, dataset_name)
    h, w = arr.shape

    # NoData 归一：UDBX 内部填充值（-32768 / -9999 / 255 / 65535 等）→ -9999
    nv = meta['novalue']
    if nv is not None:
        with np.errstate(invalid='ignore'):
            nodata_mask = np.isclose(arr, float(nv)) | np.isnan(arr)
    else:
        nodata_mask = np.isnan(arr)
    arr_out = np.where(nodata_mask, NODATA_OUT, arr).astype(np.float32)

    # 地理配准：优先 UDBX 注册范围，异常时回退铁基准
    left, top = meta['left'], meta['top']
    if left is None or abs(left - LEFT) > 1 or abs(top - TOP) > 1:
        print(f"  [警告] {dataset_name} 地理范围 ({left},{top}) 与铁基准不符，按铁基准输出")
        left, top = LEFT, TOP
    transform = from_origin(left, top, RES, RES)
    crs = CRS.from_proj4(CRS_PROJ4)

    profile = dict(driver='GTiff', height=h, width=w, count=1, dtype='float32',
                   crs=crs, transform=transform, nodata=NODATA_OUT, compress='deflate')
    with rasterio.open(out_tif, 'w', **profile) as dst:
        dst.write(arr_out, 1)

    summary = dict(dataset=dataset_name, tif=os.path.basename(out_tif),
                   width=w, height=h, nodata_src=nv,
                   valid=int((~nodata_mask).sum()), total=int(arr.size))
    if check:
        ok_size = (w == GRID_N and h == GRID_N)
        summary['size_ok'] = ok_size
        if not ok_size:
            print(f"  [报警] {dataset_name}: 尺寸 {w}x{h} ≠ 铁基准 {GRID_N}x{GRID_N}")
        valid = arr_out[~nodata_mask]
        if valid.size:
            summary.update(min=float(valid.min()), max=float(valid.max()),
                           mean=float(valid.mean()))
    return summary


def export_batch(udbx_path, out_dir, only=None):
    """批量导出 DEFAULT_EXPORT 清单（或 --only 指定子集），打印逐栅格自检。"""
    os.makedirs(out_dir, exist_ok=True)
    imgs = {i['name']: i for i in list_rasters(udbx_path)}
    names = only if only else list(DEFAULT_EXPORT)
    results = []
    for name in names:
        tif_name = DEFAULT_EXPORT.get(name, name)
        if name not in imgs:
            print(f"[跳过] UDBX 中无数据集：{name}")
            continue
        info = imgs[name]
        if info['width'] == 192 or str(name).endswith('Tier1'):
            print(f"[跳过] 金字塔概化层：{name} ({info['width']}x{info['height']})")
            continue
        if info['width'] == 383:
            print(f"[跳过] 旧研究区错网格(383/±46000)：{name}")
            continue
        out_tif = os.path.join(out_dir, tif_name + '.tif')
        s = export_raster(udbx_path, name, out_tif)
        print(f"[导出] {s['dataset']:24s} -> {s['tif']:32s} "
              f"{s['width']}x{s['height']} 有效 {s['valid']}/{s['total']} "
              f"min={s.get('min', float('nan')):.4f} max={s.get('max', float('nan')):.4f} "
              f"mean={s.get('mean', float('nan')):.4f}")
        results.append(s)
    return results


def main():
    ap = argparse.ArgumentParser(description='UDBX → GeoTIFF 纯 Python 导出器')
    ap.add_argument('--udbx', default=os.environ.get(
        'LUNAR_UDBX', r'A+ak9HzkqIZG/2成果数据/数据源/输入.udbx'))
    ap.add_argument('--out', default=os.environ.get('LUNAR_RASTER_DIR', 'data/rasters'))
    ap.add_argument('--only', help='逗号分隔的数据集名列表（默认导出全部清单）')
    ap.add_argument('--list', action='store_true', help='仅列出 UDBX 内的栅格数据集')
    args = ap.parse_args()

    if args.list:
        for i in list_rasters(args.udbx):
            print(f"{i['dataset_id']:4d} {i['name']:30s} table={i['table']:24s} "
                  f"{i['width']}x{i['height']} block={i['block_size']}")
        return

    only = args.only.split(',') if args.only else None
    export_batch(args.udbx, args.out, only)


if __name__ == '__main__':
    main()
