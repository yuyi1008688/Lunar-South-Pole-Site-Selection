# ──────────────────────────────────────────────────────────────────
# 章节  : Utils · 跨章节工具：铁基准对齐（任意 tif 统一重采样至 384×384@240m@±46080）
# 来源  : 竞赛提交包 3工程文件/Utils_*/align_to_standard_grid.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
铁基准对齐脚本
==============
将所有非标准栅格统一转换为 384×384 @ 240m @ ±46080（Moon South Polar Stereographic）

用法：
  python align_to_standard_grid.py

扫描范围（用正斜杠避免转义）：
  - 产出栅格目录（现默认 data/rasters）
  - 产出栅格目录（现默认 data/rasters，可用 LUNAR_SCAN_DIRS 覆盖）

行为：
  - 已经合规的文件（384×384 @ 240m @ ±46080）→ 跳过
  - 不合规的文件 → 重采样到标准网格，覆盖原文件
  - 输出对齐报告

提交前跑一次即可。
"""

import os
import sys
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_bounds

# ============================================================
# 铁基准参数
# ============================================================
STANDARD_SHAPE = (384, 384)       # (height, width)
STANDARD_RES = 240.0              # 米/像素
STANDARD_BOUNDS = (-46080, -46080, 46080, 46080)  # (left, bottom, right, top)
STANDARD_CRS = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 "
    "+x_0=0 +y_0=0 +R=1737400 +units=m"
)

STANDARD_TRANSFORM = from_bounds(
    *STANDARD_BOUNDS, *STANDARD_SHAPE[::-1]  # width, height
)

# ============================================================
# 扫描目录
# ============================================================
# 扫描目录：默认原开发机路径；复现时设环境变量 LUNAR_SCAN_DIRS（os.pathsep 分隔多目录）
SCAN_DIRS = os.environ.get("LUNAR_SCAN_DIRS", os.pathsep.join([
    os.path.join("data", "rasters"),
])).split(os.pathsep)

# 跳过的文件模式（原始全幅数据不需要裁剪）
SKIP_PATTERNS = [
    "avgvisib_65s_240m.jp2",       # 全幅原始数据
    "AVGVISIB_study_5m",           # 5m中间文件（超高分辨率，不强制对齐）
]


def is_standard(filepath):
    """检查文件是否已经符合铁基准"""
    try:
        with rasterio.open(filepath) as ds:
            if ds.width != 384 or ds.height != 384:
                return False
            if abs(ds.res[0] - 240.0) > 1.0:
                return False
            ox, oy = ds.transform.c, ds.transform.f
            if abs(ox - (-46080)) > 1.0 or abs(oy - 46080) > 1.0:
                return False
            return True
    except Exception:
        return False


def should_skip(filepath):
    """判断是否应该跳过（原始全幅数据等）"""
    basename = os.path.basename(filepath)
    for pattern in SKIP_PATTERNS:
        if pattern in basename:
            return True
    # 跳过5m超高分辨率文件（18400×18400那种）
    try:
        with rasterio.open(filepath) as ds:
            if ds.res[0] < 50:  # 分辨率高于50m的不处理
                return True
    except Exception:
        pass
    return False


def align_file(filepath):
    """将单个文件对齐到铁基准网格"""
    with rasterio.open(filepath) as src:
        src_crs = src.crs
        src_nodata = src.nodata
        dtype = src.dtypes[0]

        # 逐波段重采样到标准网格
        dst_bands = []
        for band_idx in range(1, src.count + 1):
            src_band = src.read(band_idx)
            dst_band = np.zeros((384, 384), dtype=dtype)

            reproject(
                source=src_band,
                destination=dst_band,
                src_transform=src.transform,
                src_crs=src_crs,
                dst_transform=STANDARD_TRANSFORM,
                dst_crs=STANDARD_CRS,
                resampling=Resampling.bilinear,
                src_nodata=src_nodata,
                dst_nodata=src_nodata,
            )
            dst_bands.append(dst_band)

        # 合并多波段
        if len(dst_bands) == 1:
            dst_array = dst_bands[0][np.newaxis, :, :]
        else:
            dst_array = np.stack(dst_bands, axis=0)

        band_count = src.count

    # 写回原文件
    profile = {
        'driver': 'GTiff',
        'height': 384,
        'width': 384,
        'count': band_count,
        'dtype': dtype,
        'crs': STANDARD_CRS,
        'transform': STANDARD_TRANSFORM,
        'nodata': src_nodata,
    }
    with rasterio.open(filepath, 'w', **profile) as dst:
        dst.write(dst_array)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    print("=" * 70)
    print("铁基准对齐工具 -- 384x384 @ 240m @ +-46080")
    print("=" * 70)

    aligned = 0
    skipped_ok = 0
    skipped_other = 0
    errors = 0

    for scan_dir in SCAN_DIRS:
        if not os.path.exists(scan_dir):
            print(f"\n[跳过目录] {scan_dir} 不存在")
            continue

        print(f"\n--- 扫描: {scan_dir} ---")

        for root, dirs, files in os.walk(scan_dir):
            for fname in sorted(files):
                if not fname.lower().endswith('.tif'):
                    continue

                filepath = os.path.join(root, fname)

                # 跳过不需要处理的文件
                if should_skip(filepath):
                    skipped_other += 1
                    continue

                # 检查是否已经合规
                if is_standard(filepath):
                    skipped_ok += 1
                    rel = os.path.relpath(filepath, scan_dir)
                    print(f"  [OK] {rel}")
                    continue

                # 需要对齐
                try:
                    with rasterio.open(filepath) as ds:
                        old_shape = f"{ds.width}x{ds.height}"
                        old_res = f"{ds.res[0]:.1f}m"
                        old_origin = f"({ds.transform.c:.0f}, {ds.transform.f:.0f})"

                    align_file(filepath)
                    aligned += 1
                    rel = os.path.relpath(filepath, scan_dir)
                    print(f"  [DONE] {rel}")
                    print(f"           {old_shape}@{old_res}@{old_origin} -> 384x384@240m@+-46080")

                except Exception as e:
                    errors += 1
                    rel = os.path.relpath(filepath, scan_dir)
                    print(f"  [FAIL] {rel}: {e}")

    print("\n" + "=" * 70)
    print("对齐报告")
    print("=" * 70)
    print(f"  已对齐:  {aligned} 个文件")
    print(f"  已合规:  {skipped_ok} 个文件(无需处理)")
    print(f"  已跳过:  {skipped_other} 个文件(原始数据/超高分辨率)")
    print(f"  失败:    {errors} 个文件")
    print("=" * 70)

    if errors > 0:
        sys.exit(1)
    else:
        print("\n全部完成! 所有产出栅格已对齐铁基准。")


if __name__ == "__main__":
    main()
