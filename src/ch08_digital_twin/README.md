# Ch08 数字孪生

本项目语境下的"数字孪生"指**任务规划阶段（pre-mission）的静态多源数据融合与交互式决策仿真系统**（区别于工业界实时监测反馈型）。核心原则：只做展示与交互，不引入新的分析逻辑——每个场景元素都能追溯到第 2–6 章的具体产出文件。

## 一、ThreeJS_scene/（浏览器端，推荐入口）★

**无需任何环境**：Chrome / Edge 直接双击 `ThreeJS_scene/index.html`。

- 场景：融合 DEM 三维地形（高程夸张 ×3）、光照分类纹理叠加、PSR 半透明覆盖、地下水分冰体元（按丰度着色）、5,214 个Ⅰ级选址点（得分>0.55 金色 / >0.40 橙 / 其余红）、最优巡视路径（橙色）、地球远景 + 星空。
- 交互：鼠标拖拽旋转、滚轮缩放（5–400 km）；「▶ 开始飞行漫游」播放 8 阶段相机动画（火箭起飞 → 飞向月球 → 进入轨道 → 俯瞰 Shackleton → 发现水冰 → 选址点群 → 巡视路径 → 拉远收尾）。
- 若浏览器限制本地 fetch，在本目录执行 `python -m http.server 8000` 后访问 `http://localhost:8000`。

数据文件（384×384 铁基准，与 `index.html` 同目录）：

| 文件 | 内容 | 来源章节 |
| --- | --- | --- |
| dem.json | 中心区高程网格 | Ch1 DEM_fused |
| dem_texture.png / illumination_texture.png / psr_texture.png | 灰度/光照分类/PSR 纹理 | Ch1–Ch2 |
| ice_voxels.json | 地下水冰体元（L1 表层，丰度驱动） | Ch3 |
| points.json | Ⅰ级选址点（坐标+得分） | Ch5 |
| path.json | 最优巡视路径顶点（264 点） | Ch6 |

> 提交包中另有升级版大文件（`dem_outer_ring.json` 11 MB、`dem_outer_ring_1km.json` 2.7 MB、`lunar_south_dem_merged.tif` 6 MB 等，用于外围过渡地形与球面曲率方案），未被当前 `index.html` 引用，为控制仓库体积未入库；需要时可从竞赛提交包 `3工程文件/Ch08_数字孪生/ThreeJS_scene/` 获取，对应视觉需求见提交包内《Ch08_C类_ThreeJS视觉调整需求.md》。

## 二、Python 端（可选，需 pyvista 环境）

| 脚本 | 功能 |
| --- | --- |
| `Ch08_pyvista_render.py` | PyVista 离屏完整渲染：地形 + PSR 点云 + 选址点 + 路径 + 水冰体元，输出 6 张高清截图 + 360° 旋转 GIF（输入路径经 `LUNAR_REVIEW_DIR` 环境变量配置） |
| `voxel_viz_platform.py` | PyQt5 + PyVista 交互式水冰体元可视化平台（阈值/配色/图层可调；命令行 `--test` 无界面自检） |
| `VoxelVizLauncher.py` | 独立启动器：一键启动体元平台，`VOXEL_PYTHON` 环境变量可覆盖解释器 |

依赖：pyvista、vtk、PyQt5、numpy、rasterio、geopandas、matplotlib。
