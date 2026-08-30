# -*- coding: utf-8 -*-
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch08 · 数字孪生：PyVista 离屏渲染（地形+PSR+选址点+路径+水冰体元，6截图+360°GIF）
# 来源  : 竞赛提交包 3工程文件/Ch08_*/Ch08_pyvista_render.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
Ch08 数字孪生可视化 - pyvista 完整渲染脚本
地形 + PSR + 选址点 + 路径 + 体元栅格（水冰/适宜性三维体素）
输出：6张高清截图 + 360度旋转GIF
"""
import numpy as np
import rasterio
import geopandas as gpd
import pyvista as pv
import csv
import os

# ========== 配置 ==========
BASE = os.environ.get("LUNAR_REVIEW_DIR", "./data")  # 复现时设置 LUNAR_REVIEW_DIR 或改为本地审阅目录
PATHS = {
    "dem":  os.path.join(BASE, r"01_Ch01_数据底座\输出数据\DEM_fused.tif"),
    "psr":  os.path.join(BASE, r"02_Ch02_光照分类\输出数据\PSR_mask.tif"),
    "ice":  os.path.join(BASE, r"03_Ch03_水冰分布\输出数据\ice_density_final.tif"),
    "suit": os.path.join(BASE, r"05_Ch05_AHP选址\输出数据\suitability_classes_correct.tif"),
    "csv":  os.path.join(BASE, r"05_Ch05_AHP选址\输出数据\grade_I_coordinates.csv"),
    "path": os.path.join(BASE, r"06_Ch06_路径规划\输出数据\optimal_path_smooth_1.shp"),
}
OUT_DIR = os.path.join(BASE, "08_Ch08_截图")
os.makedirs(OUT_DIR, exist_ok=True)

ELEV_SCALE = 3      # 高程夸张系数
ICE_DEPTH_SCALE = 800  # 体素高度缩放（ice_density * 800 = 体素高度m）
VOXEL_LAYERS = 8    # 体素层数（地下层数）

# ========== 数据加载 ==========
print("[1/7] 加载数据...")

with rasterio.open(PATHS["dem"]) as ds:
    dem = ds.read(1).astype(float)
    dem = np.where(dem == ds.nodata, np.nan, dem)
    bounds = ds.bounds
    res = ds.res[0]
    ny, nx = dem.shape

with rasterio.open(PATHS["psr"]) as ds:
    psr = ds.read(1).astype(float)

with rasterio.open(PATHS["ice"]) as ds:
    ice = ds.read(1).astype(float)
    ice = np.where(ice == ds.nodata, np.nan, ice)

with rasterio.open(PATHS["suit"]) as ds:
    suit = ds.read(1).astype(float)

with open(PATHS["csv"], "r") as f:
    reader = csv.DictReader(f)
    pt_data = list(reader)
    pt_x = np.array([float(r["x"]) for r in pt_data])
    pt_y = np.array([float(r["y"]) for r in pt_data])
    pt_score = np.array([float(r["score"]) for r in pt_data])

gdf = gpd.read_file(PATHS["path"])
path_coords = np.array(list(gdf.geometry.iloc[0].coords))

print(f"  DEM: {dem.shape}, [{np.nanmin(dem):.0f}, {np.nanmax(dem):.0f}]m")
print(f"  PSR: {(psr==1).sum()} pixels")
print(f"  Ice: {np.sum(~np.isnan(ice))} valid pixels")
print(f"  Points: {len(pt_x)}")
print(f"  Path: {len(path_coords)} vertices")

# ========== 坐标网格 ==========
print("[2/7] 构建坐标网格...")
x = np.linspace(bounds.left, bounds.right, nx)
y = np.linspace(bounds.top, bounds.bottom, ny)
X, Y = np.meshgrid(x, y)
Z = dem * ELEV_SCALE

# ========== 场景构建 ==========
print("[3/7] 构建场景...")
plotter = pv.Plotter(window_size=(1920, 1080), off_screen=True)
plotter.set_background("black")

# --- 地形网格 ---
grid = pv.StructuredGrid(X, Y, np.nan_to_num(Z))
grid.point_data["elevation"] = np.nan_to_num(Z).flatten(order="F")
plotter.add_mesh(grid, scalars="elevation", cmap="terrain",
                 lighting=True, smooth_shading=True,
                 clim=[np.nanmin(Z), np.nanmax(Z)], opacity=1.0,
                 name="terrain")

# --- PSR掩膜（青色半透明点云） ---
psr_mask_flat = (psr == 1).flatten(order="F")
psr_z = np.nan_to_num(Z) + 3
psr_grid = pv.StructuredGrid(X, Y, psr_z)
psr_pts = psr_grid.points[psr_mask_flat]
if len(psr_pts) > 0:
    plotter.add_mesh(pv.PolyData(psr_pts), color="cyan", opacity=0.35,
                     point_size=3, render_points_as_spheres=True, name="psr")

# --- 体元栅格：水冰三维体素（地下） ---
print("[4/7] 构建水冰体元栅格...")
ice_valid_mask = ~np.isnan(ice)
if np.any(ice_valid_mask):
    # 为每个有冰像元创建地下体素柱
    voxel_centers = []
    voxel_vals = []
    ice_max = np.nanmax(ice)

    for iy in range(ny):
        for ix in range(nx):
            if ice_valid_mask[iy, ix]:
                val = ice[iy, ix]
                # 体素高度 = ice_density / max * ICE_DEPTH_SCALE
                col_height = (val / ice_max) * ICE_DEPTH_SCALE
                n_layers = max(1, int(VOXEL_LAYERS * (val / ice_max)))
                surf_z = Z[iy, ix]
                for il in range(n_layers):
                    vz = surf_z - (il + 0.5) * (col_height / n_layers)
                    voxel_centers.append([X[iy, ix], Y[iy, ix], vz])
                    voxel_vals.append(val)

    if voxel_centers:
        voxel_cloud = pv.PolyData(np.array(voxel_centers))
        voxel_cloud["ice_density"] = np.array(voxel_vals)
        plotter.add_mesh(voxel_cloud, scalars="ice_density",
                         cmap="Blues", opacity=0.6,
                         point_size=8, render_points_as_spheres=True,
                         clim=[0, ice_max], name="ice_voxel")

# --- 体元栅格：适宜性三维体素（地上） ---
print("[5/7] 构建适宜性体元栅格...")
suit_colors = {5: "#FFD700", 4: "#9ACD32", 3: "#FFFF00", 2: "#FF8C00", 1: "#FF0000"}
suit_heights = {5: 1500, 4: 1200, 3: 900, 2: 600, 1: 300}

for cls_val in [5, 4, 3, 2, 1]:
    cls_mask = (suit == cls_val)
    if not np.any(cls_mask):
        continue
    cls_flat = cls_mask.flatten(order="F")
    cls_z = np.nan_to_num(Z) + suit_heights[cls_val] / 2
    cls_grid = pv.StructuredGrid(X, Y, cls_z)
    cls_pts = cls_grid.points[cls_flat]
    if len(cls_pts) > 0:
        # 用GeometricObjects构建体素柱
        for pt in cls_pts[::20]:  # 降采样，每20个取1个
            cube = pv.Cube(center=pt, x_length=res*0.8, y_length=res*0.8,
                          z_length=suit_heights[cls_val])
            plotter.add_mesh(cube, color=suit_colors[cls_val], opacity=0.15,
                             name=f"suit_voxel_{cls_val}_{int(pt[0])}_{int(pt[1])}")

# --- I级选址点 ---
print("[6/7] 添加选址点和路径...")
pt_z = np.zeros(len(pt_x))
for i in range(len(pt_x)):
    col = max(0, min(int((pt_x[i] - bounds.left) / res), nx-1))
    row = max(0, min(int((bounds.top - pt_y[i]) / res), ny-1))
    pt_z[i] = Z[row, col] + 20

pt_cloud = pv.PolyData(np.column_stack([pt_x, pt_y, pt_z]))
pt_cloud["score"] = pt_score
plotter.add_mesh(pt_cloud, scalars="score", cmap="YlOrRd",
                 point_size=4, render_points_as_spheres=True,
                 clim=[0, 0.65], name="grade_I_points")

# --- 最优路径线 ---
path_z = np.zeros(len(path_coords))
for i, (px, py) in enumerate(path_coords):
    col = max(0, min(int((px - bounds.left) / res), nx-1))
    row = max(0, min(int((bounds.top - py) / res), ny-1))
    path_z[i] = Z[row, col] + 25
path_3d = np.column_stack([path_coords[:, 0], path_coords[:, 1], path_z])
path_line = pv.Spline(path_3d, len(path_3d) * 3)
plotter.add_mesh(path_line, color="#FF5000", line_width=4, name="optimal_path")

# --- 光照 ---
light = pv.Light(position=(30000, 30000, 15000), focal_point=(0, 0, -2000),
                 intensity=0.8, positional=False)
plotter.add_light(light)
plotter.enable_anti_aliasing("ssaa")

# ========== 截图 ==========
print("[7/7] 截图...")

def screenshot(name, cam_pos):
    plotter.camera.position = cam_pos
    plotter.camera.focal_point = (0, 0, -2000)
    plotter.render()
    path = os.path.join(OUT_DIR, name)
    plotter.screenshot(path)
    print(f"  已保存: {name}")

screenshot("Ch08_B01_远景俯瞰.png", (0, 0, 90000))
screenshot("Ch08_B02_东南斜视.png", (60000, 60000, 35000))
screenshot("Ch08_B03_坑壁特写.png", (15000, 15000, 5000))
screenshot("Ch08_B04_选址点群.png", (10000, 10000, 15000))
screenshot("Ch08_B05_路径特写.png", (20000, -5000, 8000))
screenshot("Ch08_B06_体元栅格特写.png", (5000, 5000, 5000))

# ========== 360度旋转GIF ==========
print("生成360度旋转动画...")
gif_path = os.path.join(OUT_DIR, "Ch08_B_360旋转.gif")
plotter.open_gif(gif_path, fps=30)
for theta in np.linspace(0, 2*np.pi, 120):
    r = 70000
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    plotter.camera.position = (x, y, 40000)
    plotter.camera.focal_point = (0, 0, -2000)
    plotter.render()
    plotter.write_frame()
plotter.close()

print("\n=== pyvista 渲染完成 ===")
print(f"截图目录: {OUT_DIR}")
