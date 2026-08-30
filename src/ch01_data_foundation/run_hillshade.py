# ──────────────────────────────────────────────────────────────────
# 章节  : Ch01 · 数据底座：融合DEM山体阴影生成（4方位角，高度角15°，Z因子3）
# 来源  : 竞赛提交包 3工程文件/Ch01_*/_run_hillshade.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
为月球南极 DEM 计算 Hillshade（山体阴影），4 个方位角，高度角 15°，配 Z 因子 3。
输出 4 张 GeoTIFF + 1 张 4 合 1 对比 PNG。
"""
import os
import numpy as np
import rasterio
from rasterio.transform import Affine

# ── 路径配置：复现时改为本地实际路径，或设置环境变量 LUNAR_DATA_DIR ──
_DATA_DIR = os.environ.get("LUNAR_DATA_DIR", "./data")  # 复现时设置 LUNAR_DATA_DIR 或改为本地数据目录
DEM_PATH = os.path.join(_DATA_DIR, "DEM_fused.tif")
OUT_DIR = os.environ.get("LUNAR_OUTPUT_DIR", os.path.join("data", "output", "ch01"))  # 输出与输入分离

ALTITUDE_DEG = 15.0   # 太阳高度角
Z_FACTOR = 3.0        # 高程夸张系数（与三维展示场景的夸张设置一致）
AZIMUTHS = [315, 45, 135, 225]  # 右上 / 左上 / 左下 / 右下


def hillshade(dem, dx, dy, altitude_deg, azimuth_deg, z_factor):
    """标准 hillshade：255*(cosZ*cosS + sinZ*sinS*cos(A-aspect))"""
    # 梯度（numpy 数组：row=南, col=东；故 dy 物理上是"南向"，dx 是"东向"）
    grad_y, grad_x = np.gradient(dem, dy, dx)

    slope = np.arctan(z_factor * np.sqrt(grad_x ** 2 + grad_y ** 2))

    # aspect（从北顺时针，单位弧度，[0, 2π)）
    aspect = np.arctan2(-grad_y, grad_x)  # 北向分量 = -grad_y（南为正）
    aspect = np.where(aspect < 0, aspect + 2 * np.pi, aspect)

    zenith = np.radians(90.0 - altitude_deg)
    azim = np.radians(azimuth_deg)

    hs = 255.0 * (
        np.cos(zenith) * np.cos(slope)
        + np.sin(zenith) * np.sin(slope) * np.cos(azim - aspect)
    )
    return np.clip(hs, 0, 255).astype(np.uint8)


with rasterio.open(DEM_PATH) as src:
    dem = src.read(1).astype(np.float32)
    nodata = src.nodata
    transform = src.transform
    profile = src.profile.copy()

# 像素物理大小（米）
pixel_x = abs(transform.a)
pixel_y = abs(transform.e)

# nodata 掩膜
if nodata is not None:
    mask = (dem == nodata) | np.isnan(dem)
else:
    mask = np.isnan(dem)

os.makedirs(OUT_DIR, exist_ok=True)  # 输出目录自建
print(f"DEM 形状: {dem.shape}, 像素大小: {pixel_x} x {pixel_y} m, nodata={nodata}")

# 生成 4 个方位角的 hillshade
out_paths = {}
for az in AZIMUTHS:
    hs = hillshade(dem, pixel_x, pixel_y, ALTITUDE_DEG, az, Z_FACTOR)
    hs[mask] = 0  # nodata 区填 0（黑色）

    out_path = f"{OUT_DIR}\\hillshade_az{az:03d}_alt{int(ALTITUDE_DEG):02d}.tif"
    out_profile = profile.copy()
    out_profile.update(dtype="uint8", count=1, nodata=0, compress="lzw")

    with rasterio.open(out_path, "w", **out_profile) as dst:
        dst.write(hs, 1)
    out_paths[az] = out_path
    print(f"[OK] {out_path}")

# 生成 4 合 1 对比图（PNG 方便用户挑）
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()
    for ax, az in zip(axes, AZIMUTHS):
        with rasterio.open(out_paths[az]) as src:
            arr = src.read(1)
        im = ax.imshow(arr, cmap="gray", vmin=0, vmax=255)
        ax.set_title(f"Azimuth = {az}° (Height = {ALTITUDE_DEG}°, Z = {Z_FACTOR})",
                     fontsize=13, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Hillshade (0=dark, 255=light)")

    fig.suptitle("Lunar South Pole Hillshade — 4 Azimuths",
                 fontsize=16, fontweight="bold")
    plt.tight_layout()
    preview_path = f"{OUT_DIR}\\hillshade_4az_preview.png"
    plt.savefig(preview_path, dpi=120, bbox_inches="tight")
    print(f"[OK] 对比图: {preview_path}")
except Exception as e:
    print(f"[WARN] 对比图生成失败: {e}")

print("\n全部完成。任一 hillshade_az***.tif 可作为影像图层叠加到地形底图上\n（建议透明度 50% + 灰度色表），用于制图与地形判读。")