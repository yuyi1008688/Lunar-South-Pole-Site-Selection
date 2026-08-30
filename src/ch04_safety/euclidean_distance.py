# ──────────────────────────────────────────────────────────────────
# 章节  : Ch04 · 安全势场：F5 连续危险距离场（combined_hazard 的欧氏距离变换，240m像元换算米）
# 来源  : 竞赛提交包 3工程文件/Ch04_*/欧式距离.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
import os
import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt

# 1. 读取 combined_hazard 的数据（只取数值，不取元数据）
HAZARD_TIF = os.environ.get("LUNAR_HAZARD_TIF", "combined_hazard.tif")  # 复现时改为本地路径
with rasterio.open(HAZARD_TIF) as src:
    data = src.read(1).astype(np.float32)
    # 注意：不从这里取 meta

# 2. 以 unit4.Slope 为模板（读取正确的范围和坐标系）
SLOPE_TIF = os.environ.get("LUNAR_SLOPE_TIF", "slope.tif")  # 铁基准模板栅格
with rasterio.open(SLOPE_TIF) as template:
    meta = template.meta.copy()  # 复制正确的投影、范围、分辨率
    template_transform = template.transform
    template_crs = template.crs
    template_bounds = template.bounds

# 3. 强制二值化
source = (data > 0.5).astype(np.uint8)

# 4. 计算欧氏距离（像素单位）
dist_pixels = distance_transform_edt(source == 0)
dist_meters = dist_pixels * 240.0

# 5. 更新元数据：保持模板的投影、范围、分辨率，只改数据类型
meta.update({
    'dtype': 'float32',
    'nodata': -9999,
    'crs': template_crs,
    'transform': template_transform,
})

# 6. 保存
with rasterio.open("distance_raw.tif", 'w', **meta) as dst:
    dst.write(dist_meters.astype(np.float32), 1)

print(f"✅ distance_raw.tif 生成完成！")
print(f"   坐标系: {meta['crs']}")
print(f"   范围: {template_bounds}")
print(f"   行列数: {meta['width']} x {meta['height']}")