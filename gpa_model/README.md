# SuperMap GPA 地理处理自动化模型

> **格式说明**：本目录下的 `GPA.smwu` 和 `GPA模型.xml` 为 **SuperMap 专有格式**工作流文件，需 **SuperMap iDesktopX 2026 (V12.0.0+)** 及以上版本打开运行，其他 GIS 软件无法直接读取。这些文件记录的是本项目团队设计的工作流节点连接与参数配置（原创设计），不包含 SuperMap 软件本身的任何代码或许可。
>
> 竞赛官方加分项「SuperMap 地理处理建模」的实体承载：串联第 1～6 章的**数据进—模型出**端到端工作流，50 个处理步骤一键运行。

## 文件清单

| 文件 | 大小 | 说明 |
| --- | --- | --- |
| `GPA模型.xml` | 288 KB | 50 步模型（iDesktopX「处理自动化 → 导入模型」直接导入） |
| `GPA.smwu` | 696 KB | GPA 总工作空间（打开后需重连数据源到本地 udbx 目录，模型原引用开发机 `e:/chaotu/...` 绝对路径） |

> XML 由 GPA 12.0.1（iDesktopX 试用版）导出；Ch6 重建手册基于 iDesktopX 12.1 编写，两者向下兼容。

## 设计目标

- **可复现**：所有中间产物由模型节点自动生成，参数显式写入模型，杜绝手工操作的随机性；
- **可审查**：每个章节对应一组可追溯的模型节点，计算过程透明；
- **铁基准自洽**：全流程统一在 384×384 栅格、240 m 分辨率、±46,080 m、月球南极立体投影下运行，跨章栅格可直接代数运算。

## 架构

「多数据源 + 栅格代数串联 + 专项分析工具 + Python 工具」四层结构；第 4 章（unit4）、第 6 章（M8）由独立子模型产出，主干模型调用其 UDBX 结果。

**数据源组织（7 个 UDBX）**：

| 数据源 | 职责 | 关键数据集 |
| --- | --- | --- |
| 输入.udbx | 原始数据底座 | LOLA_240m_dem、CE2（融合输入）、AVGVISIB 等 |
| 第二章.udbx | 光照产物 | AVGVISIB_probability、五级分类、PSR 掩膜 |
| 第三章.udbx | 水冰产物 | F2_wang_kde_final |
| unit4.udbx | 安全势场产物 | FoS_raw、VRM、combined_hazard、distance_raw |
| 第五章.udbx | 选址产物 | F1~F5、suitability 系列、约束掩膜 |
| 输出.udbx | 最终产出集 | DEM_fused、F3/F4/F5、optimal_path 等 |
| 第六章.udbx | 路径规划 | cost_surface 等 |

**节点工具族（14 类，共 50 步）**：栅格代数 rastermathanalyst（28）｜邻域统计 neighbourhoodstatistics（6）｜坡度/坡向/曲率（各 2）｜栅格矢量化 rastertovector（2）｜重分类 reclass（1）｜自定义 Python 工具「Ⅰ级区统计」（1）｜加权总和 rasterweightsum（1）｜直线距离 straightdistance（1）｜核密度 kerneldensity（1）｜重设坐标系 setprojection（1）｜耗费路径 costpath（1）｜耗费距离 costdistance（1）。

## 50 步节点清单

| # | 节点 | 工具 | 所属环节 |
| --- | --- | --- | --- |
| 1 | 五级光照分类 | rastermathanalyst | Ch2 |
| 2–5 | 掩膜 1–4（sPSR/subPSR/PSR/连续光照） | rastermathanalyst | Ch2 |
| 6 | F1 归一化 | rastermathanalyst | Ch5 |
| 7 | 硬约束掩膜 | rastermathanalyst | Ch5 |
| 8 | WLC 加权叠加 | rastermathanalyst | Ch5 |
| 9 | 应用硬约束 | rastermathanalyst | Ch5 |
| 10 | Jenks 五级分类 | reclass | Ch5 |
| 11/24 | 坡度分析（×2） | calculateslope | Ch1/Ch4 |
| 12 | Ⅰ级区统计 | 自定义 Python 工具 | Ch5 |
| 13 | 加权总和（DEM 融合 0.7/0.3） | rasterweightsum | Ch1 |
| 14/29 | 坡向分析（×2） | calculateaspect | Ch1/Ch4 |
| 15/18 | DEM 曲率计算（×2） | calculatecurvature | Ch1/Ch4 |
| 16 | A3 FoS_raw 前置计算 | rastermathanalyst | Ch4 |
| 17/27 | VRM 5×5 / 3×3 矢量崎岖度 | rastermathanalyst | Ch4 |
| 19 | 合并危险区 | rastermathanalyst | Ch4 |
| 20–22, 31, 33, 38–39 | 邻域统计（nx/ny/nz × 3×3/5×5 求和，共 6 步） | neighbourstatistics | Ch4 |
| 23 | 转 NoData | rastermathanalyst | Ch4 |
| 25–26, 33 | nx / ny / nz 分量 | rastermathanalyst | Ch4 |
| 28 | 曲率危险 | rastermathanalyst | Ch4 |
| 30 | 硬约束掩膜_1 | rastermathanalyst | Ch4 |
| 32 | 计算 F4 通信同步覆盖度 | rastermathanalyst | Ch4 |
| 34 | 生成距离栅格（直线） | straightdistance | Ch4 |
| 35 | F5 归一化 | rastermathanalyst | Ch4 |
| 36 | 平地处理（FoS 封顶 3.0） | rastermathanalyst | Ch4 |
| 37 | F3 归一化 | rastermathanalyst | Ch4 |
| 40 | PSR 掩膜 | rastermathanalyst | Ch3 |
| 41 | 核密度分析（KDE，半径 4,800 m） | kerneldensity | Ch3 |
| 42 | 重设坐标系 | setprojection | Ch3 |
| 43 | 归一化（P98） | rastermathanalyst | Ch3 |
| 44 | 栅格矢量化（终点提取） | rastertovector | Ch6 |
| 45 | 栅格路径矢量化兼 B 样条路径平滑 | rastertovector | Ch6 |
| 46 | 坡度成本 | rastermathanalyst | Ch6 |
| 47 | 地形崎岖度成本 | rastermathanalyst | Ch6 |
| 48 | 计算最短路径（数据集，LCP） | costpath | Ch6 |
| 49 | 生成距离栅格（自定义耗费） | costdistance | Ch6 |
| 50 | 综合成本面合成 | rastermathanalyst | Ch6 |

## 关键公式台账（栅格代数实际执行版）

```text
五级光照分类   Con(AVGVISIB<=1e-6,1, Con(<=0.001,2, Con(<=0.128,3, Con(<=0.264,4,5))))
DEM 融合       0.7×LOLA + 0.3×CE2
FoS_raw        (1500+1500·1.62·1.0·cos²(slope)·tan35°)/(1500·1.62·1.0·sin(slope)·cos(slope))
F3             Con(FoS_fixed<0.5,0, Con(FoS_fixed>3,1,(FoS_fixed-0.5)/2.5))
F4             AVGVISIB × AVGVISIB_EARTH
F5             Con(distance_raw>5000,1.0, distance_raw/5000.0)
曲率危险       (profileCurvature>0.001 AND slope>20) OR (profileCurvature>0.002)
硬约束         Con(slope<=20 ∧ AVGVISIB>=0.001 ∧ distance>=240,1,0)
WLC            0.331·F1+0.258·F2+0.191·F3+0.126·F4+0.094·F5
应用硬约束     suitability_final = suitability_raw × constraint_mask
Jenks 断点     约 0.332 / 0.467 / 0.542 / 0.572 / 0.652（最高级即Ⅰ级）
```

## 工具链选择逻辑

**优先 GPA 原生工具**（坡度/坡向/曲率、栅格计算器、邻域统计、重分类、LCP、B 样条平滑、栅格转点）；**Python 启动器模式兜底**原生工具跑不通或无法一步到位的环节：

- 撞击坑 CSV 解析 → 缓冲 → 转栅格（`src/ch04_safety/VectorToHazard*.py`）；
- 欧氏距离变换（`src/ch04_safety/euclidean_distance.py`）；
- Ⅰ级区统计（SuperMapPythonProcessFactory 自定义工具）；
- 剖面采样与能量递推——**能量递推是带状态机的时序循环，DAG 原理上无法表达**，故封装边界止于 GIS 分析链，能量仿真由外部脚本衔接（`src/ch06_path_planning/energy_simulation.py`），这一分工本身即方法论。

## 部署步骤

1. iDesktopX 2026（V12.0.0+）→「处理自动化」→「导入模型」→ 选择 `GPA模型.xml`；
2. 模型属性中把各数据源节点重连到本地 udbx（默认引用 `e:/chaotu/...` 开发机路径）；
3. 运行前关闭 Python 窗口与相关数据标签页，解除 UDBX 文件锁；
4. 运行后可按下列对标值核对：
   - slope_ok / psr_ok / hazard_ok = 79,545 / 30,252 / 101,577；
   - candidate = 2 像元，终点 (40,920, 10,920) 处 = 1；
   - LCP 总长 4.039 km（B 样条系数 5），危险区穿越 0。

## 关键参数（模型变量化）

| 变量 | 值 |
| --- | --- |
| var_spsr_thresh / var_subpsr_thresh | 1e-6 / 0.001 |
| var_f1_min / var_f1_max | 0.001 / 0.264 |
| var_slope_max_hard | 20° |
| var_hazard_buffer | 240 m |
| var_fos_C / var_fos_phi | 1,500 Pa / 35° |
| var_w1~var_w5 | 0.331 / 0.258 / 0.191 / 0.126 / 0.094 |
| var_cw_slope/vrm/hazard | 0.45 / 0.25 / 0.30 |
| var_analysis_res | 240 m |
| KDE 搜索半径 | 4,800 m |
| B 样条系数 | 5 |
