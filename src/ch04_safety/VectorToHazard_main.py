# -*- coding: utf-8 -*-
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch04 · 安全势场：矢量撞击坑 → CEB缓冲区(0.75×D) → 危险区栅格
# 来源  : 竞赛提交包 3工程文件/Ch04_*/VectorToHazard_main.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
矢量撞击坑 → 缓冲区 → 栅格工具（可作为宿主工作流的外部计算工具）
输入：矢量数据（Shapefile 等矢量数据源）+ DEM 参数
输出：crater_hazard.tif（二值掩膜，1=危险区）
"""

import os
import sys
import traceback
import numpy as np

# ===== 定义外部调用输入输出参数 =====
# input, vector_path, string, 矢量数据路径(含撞击坑),
# input, diam_field, string, 直径字段名,
# input, dem_xmin, double, DEM左边界(-46080),
# input, dem_xmax, double, DEM右边界(46080),
# input, dem_ymin, double, DEM下边界(-46080),
# input, dem_ymax, double, DEM上边界(46080),
# input, dem_cols, int, DEM列数(384),
# input, dem_rows, int, DEM行数(384),
# input, output_dir, string, 输出目录,
# output, crater_hazard, string, 输出crater_hazard路径,


class VectorToHazard(object):
    """矢量撞击坑 → 缓冲区 → 栅格（用ogr读取矢量）"""

    def execute(self, keyargs):
        try:
            # ===== 1. 获取参数 =====
            vector_path = keyargs.get("vector_path", "")
            diam_field = keyargs.get("diam_field", "DIAMETER_KM")
            xmin = float(keyargs.get("dem_xmin", -46080))
            xmax = float(keyargs.get("dem_xmax", 46080))
            ymin = float(keyargs.get("dem_ymin", -46080))
            ymax = float(keyargs.get("dem_ymax", 46080))
            cols = int(keyargs.get("dem_cols", 384))
            rows = int(keyargs.get("dem_rows", 384))
            output_dir = keyargs.get("output_dir", "")

            if not vector_path:
                return {"status": "ERROR: 缺少 vector_path 参数"}
            if not os.path.exists(vector_path):
                return {"status": f"ERROR: 矢量文件不存在: {vector_path}"}
            if not output_dir:
                output_dir = os.path.dirname(vector_path)
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)

            res = (xmax - xmin) / cols
            buf_scale = 0.75

            # ===== 2. 用ogr读取矢量数据 =====
            try:
                from osgeo import ogr, osr, gdal
            except ImportError:
                return {"status": "ERROR: GDAL未安装，无法读取矢量"}

            # 打开矢量数据
            ds = ogr.Open(vector_path, 0)  # 0=只读
            if not ds:
                return {"status": f"ERROR: 无法打开矢量: {vector_path}"}

            lyr = ds.GetLayer(0)
            if not lyr:
                return {"status": "ERROR: 无法获取图层"}

            # 获取空间参考
            srs = lyr.GetSpatialRef()
            if srs is None:
                return {"status": "ERROR: 矢量数据无空间参考"}

            # 创建目标空间参考（月球南极投影 EPSG:104903）
            target_srs = osr.SpatialReference()
            target_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            target_srs.ImportFromEPSG(104903)

            # 创建坐标转换（如果矢量坐标不是月球投影）
            transform = None
            if not srs.IsSame(target_srs):
                transform = osr.CoordinateTransformation(srs, target_srs)

            # ===== 3. 提取撞击坑（坐标+直径） =====
            craters = []
            for feature in lyr:
                geom = feature.GetGeometryRef()
                if geom is None:
                    continue

                # 获取直径
                diam = feature.GetField(diam_field)
                if diam is None or diam <= 0:
                    continue

                # 获取中心点坐标
                if geom.GetGeometryType() == ogr.wkbPoint:
                    x, y = geom.GetX(), geom.GetY()
                elif geom.GetGeometryType() in [ogr.wkbPolygon, ogr.wkbMultiPolygon]:
                    centroid = geom.Centroid()
                    x, y = centroid.GetX(), centroid.GetY()
                else:
                    continue

                # 坐标转换
                if transform:
                    pt = ogr.Geometry(ogr.wkbPoint)
                    pt.AddPoint(x, y)
                    pt.Transform(transform)
                    x, y = pt.GetX(), pt.GetY()

                radius_m = diam * 1000 * buf_scale  # 假设直径单位是km，转为米
                craters.append({'x': x, 'y': y, 'radius': radius_m})

            ds = None

            if not craters:
                return {"status": "ERROR: 未提取到有效撞击坑"}

            print(f"✅ 提取到 {len(craters)} 个撞击坑")

            # ===== 4. 生成栅格 =====
            raster = np.zeros((rows, cols), dtype=np.uint8)

            for c in craters:
                col = int((c['x'] - xmin) / res)
                row = int((ymax - c['y']) / res)
                radius_px = int(c['radius'] / res)

                if radius_px <= 0:
                    continue

                # 绘制圆形缓冲区
                for dr in range(-radius_px, radius_px + 1):
                    for dc in range(-radius_px, radius_px + 1):
                        if dr * dr + dc * dc <= radius_px * radius_px:
                            r = row + dr
                            c_idx = col + dc
                            if 0 <= r < rows and 0 <= c_idx < cols:
                                raster[r, c_idx] = 1

            # ===== 5. 保存TIFF =====
            output_path = os.path.join(output_dir, 'crater_hazard.tif')

            driver = gdal.GetDriverByName('GTiff')
            ds = driver.Create(output_path, cols, rows, 1, gdal.GDT_Byte)
            ds.SetGeoTransform([xmin, res, 0, ymax, 0, -res])
            ds.SetProjection(target_srs.ExportToWkt())
            ds.GetRasterBand(1).WriteArray(raster)
            ds.GetRasterBand(1).SetNoDataValue(0)
            ds.FlushCache()
            ds = None

            hazard_pixels = int(np.sum(raster))
            hazard_area = hazard_pixels * res * res / 1e6

            print(f"✅ 危险区像元数: {hazard_pixels}")
            print(f"✅ 危险区面积: {hazard_area:.2f} km²")

            msg = f"OK: {len(craters)}个坑 | 危险区 {hazard_pixels} 像元 ({hazard_area:.2f} km²)"
            return {"status": msg, "crater_hazard": output_path}

        except Exception as e:
            tb = traceback.format_exc()
            print(f"ERROR: {e}\n{tb}")
            return {"status": f"ERROR: {str(e)}"}


if __name__ == "__main__":
    # 独立运行示例（路径改为本地实际值；被宿主工作流调用时由参数传入）
    tool = VectorToHazard()
    result = tool.execute({
        "vector_path": r"./craters_cleaned.shp",
        "diam_field": "DIAMETER_KM",
        "dem_xmin": "-46080",
        "dem_xmax": "46080",
        "dem_ymin": "-46080",
        "dem_ymax": "46080",
        "dem_cols": "384",
        "dem_rows": "384",
        "output_dir": r"./"
    })
    print(result)