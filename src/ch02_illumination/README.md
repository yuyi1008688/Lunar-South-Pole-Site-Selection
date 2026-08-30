# Ch02 光照分类（纯 Python 模块）

本章核心处理已补写为**纯 Python 模块** `illumination.py`（numpy 向量化），
与存档产物逐像元对标一致（L5.1：五级分类一致率 100%、F1 max|diff|=2.97e-8、
覆盖率 sPSR 19.51% / subPSR 1.00% / PSR 20.52% / 连续光照 25.05%）：

```bash
python src/ch02_illumination/illumination.py
# 输入 data/rasters/AVGVISIB_probability.tif（量纲自动体检：max>1 才 ÷25000）
# 输出 data/output/ch02/C2_illumination_class.tif + F1_illumination_norm.tif
```

算法口径（与当年栅格代数 Con 嵌套表达式逐项等价，算法零改动）：

- 双 PSR 分类：sPSR（p<1e-6）/ subPSR（[1e-6,0.001)）/ PSR 联合掩膜（<0.001）；
- 五级光照分类：class1 p≤1e-6 → class5 p>0.264（断点 1e-6/0.001/0.128/0.264）；
- F1 分段线性归一化：p<0.001 → 0；[0.001,0.264] 线性；>0.264 → 1
  （0.20 是 Ch05 选址硬约束阈值，与 F1 公式各司其职）。

历史表达式（当年栅格代数口径，留档对照；亦可导入 `gpa_model/GPA模型.xml` 节点 1 查看）：

```text
C2_illumination_class = Con(AVGVISIB<=1e-6, 1,
                       Con(AVGVISIB<=0.001, 2,
                       Con(AVGVISIB<=0.128, 3,
                       Con(AVGVISIB<=0.264, 4, 5))))
```

方法与阈值依据详见 [docs/methodology.md 第 2 章](../../docs/methodology.md)。
