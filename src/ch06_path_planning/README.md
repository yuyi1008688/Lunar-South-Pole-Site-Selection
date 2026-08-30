# Ch06 路径规划与能量仿真（纯 Python 实现）

> 本章主分析链已**完全脱离商业 GIS 软件**，仅用开源栈（标准库 + numpy + scipy + rasterio）实现，等价性经 L0–L4 五级逐像元/逐指标/逐参数验证（见 `tests/解耦与精度验证报告.md`）。竞赛期的工程化封装（处理自动化建模工作流）仅在 `gpa_model/` 作为历史留档，主分析链不依赖它。

## 文件清单与执行顺序

```bash
# 前置：按 data/README.md 导出栅格到 data/rasters/

python stage1_select_mining_target.py     # ① 参数定标：D_min/D_max、F2 阈值分布
python stage1_validate_endpoint.py        # ② 终点 (40920,10920) 六项合规验证
python stage2_purepython.py               # ③ 主分析链：四掩膜→candidate→LCP→B样条→50m剖面
python energy_simulation.py               # ④ 三态状态机能量递推（表6-3 + 能量曲线）
python plot_energy_curve_svg.py           # ⑤ 零依赖能量曲线 SVG（可选）
```

## 纯 Python 等价实现对照

| 原算子（商业 GIS 内嵌组件） | 本实现 | 文件 |
| --- | --- | --- |
| 栅格读取 + 8 变换运行时定标 | `RasterGrid`（rasterio 读标准 GeoTIFF，方向确定，仅轻量 sanity check） | `../utils/raster_grid.py` |
| 栅格代数 Con 表达式（掩膜/求交） | numpy 向量化布尔运算 | `stage2_purepython.py` |
| 站址缓冲→栅格化→Con 补 0（dist_ok） | 像元中心到站址平面距离 ≤ 5840 m | `stage2_purepython.py::build_dist_ok` |
| LCP（两点最小耗费路径） | 自实现 8 邻域 Dijkstra（累积面 + 回溯方向面 + 路径回溯） | `lcp_dijkstra.py` |
| BSPLINE smoothDegree=5（5→3→2 降级） | scipy `splprep/splev`（弦长参数化 + 危险区回检自动降级） | `bspline_smooth.py` |
| 50m 剖面（最近邻取值） | 50m 等距重采样 + **双线性**插值取值 | `stage2_purepython.py` |

LCP 的数学定义与商业 GIS 耗费距离同构：边成本 = (两端成本均值) × 移动距离（正交 240 m、对角 240√2 m）。

## 对标值（历史基准 vs 纯 Python 重现，均 PASS）

| 指标 | 历史基准 | 纯 Python 重现 |
| --- | --- | --- |
| dist_ok 像元数 | 1,202 | 1,202 |
| slope_ok / psr_ok / hazard_ok | 79,545 / 30,252 / 101,577 | 完全一致 |
| candidate | 2（终点处=1） | 完全一致 |
| 成本累积面（147,456 像元） | 商业软件产出 | 相对误差 <1e-6 占比 100%（ULP 级） |
| 回溯方向面 | 8 方向编码 | 映射后 100% 一致 |
| 几何路径像元（optimal_path_raw） | 247 | 247（100% 重合） |
| 任务路径（5阶测试档） | 4.039 km / 17 顶点 / 穿越 0 / 剖面 82 点 | 4.029 km（-0.26%）/ 17 顶点 / 穿越 0 / 剖面 82 点 |
| 几何路径（2阶交付档） | 65.295 km（path.json） | 65.558 km（+0.40%）/ 穿越 0 |
| 能量仿真最低 SoC | 215.2 Wh（>200 ✓） | 210.8 Wh（>200 ✓，Δ-4.4） |
| 驻留充电 / 月昼 | 8 次 / 8 个月昼（<15/<10 ✓） | 3 次 / 5 个月昼（<15/<10 ✓，口径差异见报告 §L3） |

## 关键口径说明

1. **双路径口径**：任务采样路径 START(44760,10920)→END(40920,10920)（4 km，喂能量仿真）；全区几何展示路径 START→(-4680,-13080)（65 km 级，仅三维展示）。
2. **B样条双档**：5 阶测试档（对齐历史 17 顶点/4.039 km）+ 2 阶交付档（对齐 path.json 的保守安全档）；5→3→2 自动降级机制保留（实测几何路径上降级链真实触发：5 阶穿 2 处 → 3 阶穿 1 处 → 2 阶 0 穿越）。
3. **dist_ok 缓冲圆心**为站址 START（非终点）——与原 GUI 工艺一致，且恰好复现 1,202 像元 / candidate=2 两个基准。
4. **剖面取值**用双线性（提示词口径）；最近邻为历史口径，其对能量仿真驻留/月昼计数的影响见 `tests/解耦与精度验证报告.md` §L3 归因实验。
