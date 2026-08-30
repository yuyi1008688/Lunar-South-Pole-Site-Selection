# 快速上手（Quickstart）

本项目**不依赖任何商业 GIS 软件**：主分析链（含 LCP 寻路、B 样条平滑、能量仿真与四级精度验证）全部由开源 Python 栈实现。按投入程度分三级体验，从 5 分钟到完整复现。

## 级别 ①：5 分钟看懂结果（零依赖）

1. **浏览数字孪生**：用 Chrome / Edge 打开 [`src/ch08_digital_twin/ThreeJS_scene/index.html`](../src/ch08_digital_twin/ThreeJS_scene/index.html)
   - 鼠标拖拽旋转、滚轮缩放；
   - 点右上角「▶ 开始飞行漫游」看 8 阶段动画（约 43 秒）；
   - 图例：青色=PSR，蓝色体素=地下水冰，红/金点=Ⅰ级选址点，橙色线=最优巡视路径。
2. **看成果图**：[docs/results.md](results.md)（10 张成果地图 + 9 张统计验证图 + 关键数字速查表）。
3. **看解耦验证**：[tests/解耦与精度验证报告.md](../tests/解耦与精度验证报告.md)（L0–L4 五级验证：数据提取 100%、Dijkstra 累积面 ULP 级一致、路径 100% 重合、能量结论不翻盘、参数逐项对账+FoS 敏感性 25/25 复现）。

> 如果浏览器直接打开 index.html 出现数据加载失败（个别浏览器限制本地 fetch），在 `ThreeJS_scene/` 目录下执行 `python -m http.server 8000`，然后访问 `http://localhost:8000`。

## 级别 ②：跑通纯 Python 分析链

### 环境要求

- Python 3.10+（建议 3.10，匹配 GDAL / tensorflow 预编译版本）
- Windows / Linux / macOS 均可（分析链无平台专属依赖）
- 建议 8 核 CPU / 16 GB 内存（蒙特卡洛与 KDE 计算密集）

### 安装依赖

```bash
pip install -r requirements.txt
```

本项目不依赖任何商业 GIS 软件；GDAL（osgeo）若 pip 安装失败，可从 conda-forge 获取预编译 wheel；tensorflow 装 CPU 版即可（Ch03 预训练模型已随附，无需重训）。

### 第一步：导出栅格数据（一次性）

原始数据获取见 [data/README.md](../data/README.md)。拿到数据源文件后，用纯 Python 导出器解放为标准 GeoTIFF：

```bash
python src/utils/udbx_extract.py --udbx <数据源文件路径> --out data/rasters
# 自带逐栅格自检：尺寸必须 384×384、范围 ±46080，否则报警
```

### 第二步：按章节运行

```bash
python examples/run_pipeline.py --list          # 查看阶段 DAG 与数据血缘
python examples/run_pipeline.py --all           # 端到端全跑（断点续跑，约 20 秒，推荐）
python examples/run_pipeline.py --run ch05      # 只跑某章及其依赖

# 或单独运行（详见各脚本头部注释）：
#   Ch01 山体阴影 | Ch03 冰点→KDE→验证→成图 | Ch04 ECSA→FoS→VRM→距离场
#   Ch05 蒙特卡洛 | Ch06 定标→验证→LCP→平滑→剖面→能量仿真
```

Ch06 主链的一键顺序：

```bash
python src/ch06_path_planning/stage1_select_mining_target.py   # ① D_min/D_max 定标
python src/ch06_path_planning/stage1_validate_endpoint.py      # ② 终点六项合规验证
python src/ch06_path_planning/stage2_purepython.py             # ③ 四掩膜→LCP→B样条→剖面
python src/ch06_path_planning/energy_simulation.py             # ④ 三态能量递推
```

### 第三步：四级精度验证

```bash
python tests/L0_verify_extract.py        # 数据层：提取回环 100%
python tests/L1_verify_lcp.py            # 算子层：Dijkstra vs 基准累积面/方向面/路径
python tests/L2_verify_path_metrics.py   # 路径层：长度/绕路系数/穿越/剖面点数
python tests/L3_verify_energy.py         # 端到端：能量仿真结论不翻盘
python tests/L4_verify_params_sensitivity.py  # 参数层：权重/F1/Jenks/ECSA/FoS 对账 + 敏感性复现
python tests/L5_verify_ch02_ch05.py           # 全链对标：Ch02/Ch05 逐像元 vs 存档
python tests/L5_verify_endtoend.py            # 一键流水线产物 vs 全部既有基准
```

### 路径配置方式

各脚本以环境变量配置路径（默认值均为相对路径，开箱即用）：

```bash
export LUNAR_RASTER_DIR=data/rasters     # 栅格目录（udbx_extract 输出）
export LUNAR_OUTPUT_DIR=data/output      # 分析产物输出目录
export LUNAR_DATA_DIR=/path/to/data      # Ch01/Ch03 数据目录（可选）
export LUNAR_FACTOR_DIR=/path/to/F1_F5   # Ch05 五因子层目录（可选）
```

## 级别 ③：竞赛期工程化封装留档（可选）

`gpa_model/` 保留了竞赛期把全流程封装为**处理自动化建模工作流**（50 步）的历史版本，用于展示"数据进—模型出"的工程化思路。它需要相应商业 GIS 平台才能打开（软件本体不随仓库分发），**主分析链不依赖它**——在纯 Python 环境下完成 Ch01–Ch06 全部分析请走级别 ②。

## 常见问题

**Q：脚本提示找不到数据文件？**
A：栅格数据不随仓库分发（体积大），请先按 [data/README.md](../data/README.md) 获取原始数据并运行 `src/utils/udbx_extract.py` 导出到 `data/rasters/`。

**Q：Three.js 场景数据加载失败？**
A：用本地 HTTP 服务（见级别 ① 的提示），或将目录托管到任意静态服务器。

**Q：KDE 结果与论文数字有偏差？**
A：确认使用 `POINT_X/POINT_Y` 列（不是 Latitude/Longitude）且 `bw_method=0.11`；scipy 的实际带宽 = bw×std，与冰点分布的标准差直接相关。

**Q：蒙特卡洛"位移标准差<500m"能复现吗？**
A：不能——那是历史文档的过度声明（详见验证报告 §L4.2）。可复现的稳健性结论是：基准最优点在三种掩膜口径下均为推荐站址 (44760,10920)；系统性小扰动 40 组中 32 组原地不动；±10% 权重扰动下适宜性空间格局 Spearman≈1.000。存档的蒙特卡洛 CSV 是 P0 坐标 bug 修复前产物，其坐标为百万级、统计量不可用。

**Q：L3 能量仿真的驻留次数与历史基准 8 次不同？**
A：历史剖面用最近邻取值（台阶光照剖面），本仓库交付剖面按双线性插值（消除台阶噪声），二者驱动的 DWELL/WAIT 切换次数不同；工程验收结论（SoC>200 / 驻留<15 / 月昼<10）在两种口径下均不翻盘。归因实验见验证报告 §L3。
