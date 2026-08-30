# 数据说明（Data）

> 月球原始数据体积大（LOLA DEM 全幅 GB 级），**不随仓库分发**。本页给出 8 项数据的来源、下载入口与预处理规范，全部数据可免费获取。

## 一、数据来源总览

| # | 数据名称 | 来源 | 获取方式 | 分辨率 | 主要用途 |
| --- | --- | --- | --- | --- | --- |
| 1 | LOLA DEM | NASA PDS | 官网下载 | 5 m | 地形分析主底图 |
| 2 | CE-2 GRAS DEM | 嫦娥二号数据门户 | 官网下载 | 20 m | 地形融合（30% 权重） |
| 3 | 光照概率 AVGVISIB | NASA PDS | 官网下载 | 240 m | 光照分类基础 |
| 4 | 地球可视概率 AVGVISIB_EARTH | NASA PDS | 官网下载 | 240 m | 通信覆盖因子 |
| 5 | Diviner 温度+冰深 | NASA PDS | 官网下载 | ~1 km | 冰稳定深度（体元分层） |
| 6 | LEND 超热中子 | NASA PDS | 官网下载 | ~10 km | 水冰辅助验证 |
| 7 | Wang 2025 冰点（M³ 冰识别） | 论文补充材料 | 论文补充下载 | ~140 m（点） | Ch03 F2 主证据层 |
| 8 | Robbins 撞击坑数据库 | NASA PDS | 官网下载 | — | F5 避障因子 |

**下载入口**：

- NASA PDS / ODE 行星数据系统：https://ode.rsl.wustl.edu/moon/ （PDS3/PDS4 归档格式）
- 嫦娥二号数据门户（国家航天局）：https://moon.bao.ac.cn/
- Wang et al. (2025, Icarus)：论文补充材料（`wang2025_ice_pixel_positions_spectra.xlsx`）
- LPNS 氢丰度（配套交叉验证）：随 Wang 数据一同使用的 `LPNS_hydrogen_south_of_85S_subset.csv`（Lunar Prospector 中子谱仪）

## 二、各数据详情与 PDS 编号

1. **LOLA DEM**：`ldem_875s_5m_float.img + .lbl`，PDS3 ID `LRO-L-LOLA-4-GDR-V1.0`，float32 IMG + LBL 标签，全球 87.5°S 区域。
2. **光照概率 AVVGISIB**：`avgvisib_65s_240m.img`，同上 PDS3 ID，int16 编码（概率×25000，缩放因子 0.00004），参考文献 Mazarico et al. 2011。
3. **地球可视概率**：`avgvisib_65s_240m_earth.img`，记录地球方向可视概率，编码同上。
4. **Diviner**：`dlre_prp_south.tab + .lbl + .xml`，PDS3 ID `LRO-L-DLRE-3-RDR-V1.0`，含年均温/最高温/冰稳定深度等字段。⚠️ 温度字段已废弃（IDW 插值在 PSR 内严重平滑），仅保留冰稳定深度字段。
5. **LEND**：`lend_rdr_alds_20111211.dat`，PDS3 ID `LRO-L-LEND-4/5-RDR-V1.0`，观测期 2011-12-11 至 2016-09-14，原始 ~10 km。保持 10 km 不降尺度。
6. **Robbins 撞击坑目录**：`lunar_crater_database_robbins_2018_bundle/`（PDS4 bundle），含直径/深度/坐标/形态参数。
7. **Wang 2025 冰点**：原始 40,623 行，研究区内 5,031 个唯一精确投影坐标；含 `POINT_X/POINT_Y`（南极立体 lat_ts=-90）与 M³ 光谱特征。**务必使用 POINT_X/POINT_Y 列**（Latitude/Longitude 精度不足）。
8. **CE-2 DEM**：`CE2_GRAS_DEM_20m_N001_87S000W_A.tif`，直接 GeoTIFF 无需格式转换。

## 三、预处理规范（所有栅格统一目标）

| 参数 | 规格 |
| --- | --- |
| 投影 | Moon South Polar Stereographic（南极立体，中央经线 0°，标准纬线 -90°） |
| 椭球 | IAU_2000_Moon，R = 1,737,400 m |
| 分辨率 | 240 m × 240 m |
| 网格范围 | ±46,080 m（384 × 384 像元），对应约 88.5°S–90°S |
| 坐标基准 | 原点 (0, 0) |

处理声明：所有预处理严格遵循数据官方发布方的专业处理方式——PDS3/PDS4 数据按 PDS 官网规范做格式转换/投影定义/缩放因子还原；Wang 冰点按论文配套方法处理；工具链为 QGIS（格式转换）→ SuperMap iDesktopX（空间分析）→ Python（LEND 二进制解析、撞击坑提取）。

关键步骤摘要：

1. **LOLA**：QGIS 读 PDS3 IMG（自动读 LBL 投影）→ 裁剪 ±46,080 m → 双线性重采样 240 m → GeoTIFF。
2. **CE-2**：确认投影 → 裁剪 → 重采样至与 LOLA 相同网格。
3. **融合 DEM**：栅格计算器 `fused = 0.7×LOLA + 0.3×CE2`，缺失像元互填。
4. **AVGVISIB / EARTH**：`存储值 × 0.00004` 还原概率 [0,1] → 裁剪 → GeoTIFF。
5. **Diviner**：TAB → 点图层 → 投影转换 → 裁剪（仅保留 ICE_DEPTH 字段）。
6. **LEND**：Python/GDAL 解析二进制 → 保持 10 km 导出。
7. **Wang 冰点**：读 POINT_X/POINT_Y → KDE（bw=0.11）→ sPSR/subPSR 掩膜 → P98 归一化。
8. **Robbins**：筛选 88.5°S–90°S 且直径 ≥0.5 km（223 坑）→ 投影转换 → CEB 缓冲半径 = 0.75×D。

## 四、栅格解放：一行命令导出标准 GeoTIFF

竞赛期的中间成果栅格存放在工程化数据源文件（SQLite 格式 `.udbx`）中。本仓库提供
纯 Python 导出器（`src/utils/udbx_extract.py`，仅需标准库 + numpy + rasterio），
从数据源内部直接提取并导出为带正确 CRS/transform/NoData 的 GeoTIFF——
**不需要安装任何 GIS 软件**：

```bash
# 全量导出（35 个栅格，含逐栅格自检：尺寸 384×384 / 范围 ±46080 / 值域统计）
python src/utils/udbx_extract.py --udbx <数据源文件路径> --out data/rasters

# 仅导出 Ch06 必需的 9 个
python src/utils/udbx_extract.py --udbx <数据源文件路径> --out data/rasters     --only Slope,AVGVISIB_probability,combined_hazard,cost_surface,DEM_fused,vrm_5x5_float32,distance_accumulation,backlink_direction,optimal_path_raw

# 查看数据源内全部栅格清单
python src/utils/udbx_extract.py --udbx <数据源文件路径> --list
```

导出后 `data/rasters/` 即为全部分析脚本的默认输入目录（可用 `LUNAR_RASTER_DIR`
覆盖）。原始数据源文件不入库；每个栅格导出后立即自检，出现 383 尺寸或
±46000 范围会明确报警（历史教训：范围/分辨率不整除时 GIS 软件会偷改分辨率）。

## 五、原始输入文件的组织建议

Wang 2025 冰点 xlsx、LPNS csv、撞击坑 shp 等**非数据源栅格**的原始文件
（获取方式见本文档第一、二节）建议按如下结构组织（与脚本环境变量对应）：

```
your_data_root/
├── wang2025_ice_pixel_positions_spectra.xlsx   # Ch03 冰点主数据
├── LPNS_hydrogen_south_of_85S_subset.csv      # Ch03 交叉验证
├── robins_craters.shp                          # Ch04 撞击坑（或自行从 PDS4 bundle 提取）
└── avgvisib_65s_240m_earth.jp2                 # Ch04 对地可见概率（原始 PDS3）
```

## 五、未纳入数据说明

以下数据曾下载但未用于分析：Unified_Geologic_Map（未参与任何流程）、M³ 水冰 ROI（Lemelin et al. 2021，旧克里金方案输入，LOOCV R²=0.229 已废弃）、lro_lend 2022 版（使用 2011 版）。
