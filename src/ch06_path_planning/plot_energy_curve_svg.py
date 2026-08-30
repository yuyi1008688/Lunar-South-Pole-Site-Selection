# -*- coding: utf-8 -*-
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch06 · 路径规划：零依赖能量曲线 SVG 出图（三面板：SoC+200Wh红线 / 光照 / 坡度）
# 来源  : 竞赛提交包 3工程文件/Ch06_*/plot_energy_curve_svg.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
plot_energy_curve_svg.py —— 零依赖出图：能量曲线（含 200Wh 红线）→ SVG

不需要 matplotlib / 任何第三方库，纯 Python 标准库，
任意 Python 3.10+ 环境直接运行（纯标准库，零第三方依赖）。

输入：data/output/stage3_energy/energy_curve.csv
输出：data/output/stage3_energy/energy_curve.svg
      · 浏览器双击打开即可查看/打印
      · Word/WPS 可直接插入 SVG；需要 PNG 时浏览器截图或

图样式对照论文图10：三面板 = SoC曲线+200Wh红线+驻留/等待标记、
沿程光照率、沿程坡度。
"""

import csv
import os

_base   = os.environ.get("LUNAR_OUTPUT_DIR", os.path.join("data", "output"))
CSV_IN  = os.path.join(_base, "stage3_energy", "energy_curve.csv")
SVG_OUT = os.path.join(_base, "stage3_energy", "energy_curve.svg")

E_MIN, SOC_DWELL, E_CAP = 200.0, 600.0, 1040.0

# ---------------- 读数据 ----------------
rows = []
with open(CSV_IN, newline="", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        rows.append((float(r["s_m"]) / 1000.0, float(r["soc_wh"]),
                     r["state"], float(r["visib"]), float(r["slope_deg"])))
if not rows:
    raise SystemExit(f"读不到数据：{CSV_IN}")

ss   = [r[0] for r in rows]
soc  = [r[1] for r in rows]
vis  = [r[3] for r in rows]
slp  = [r[4] for r in rows]
smax = max(ss) or 1.0
vmax = 1.0
pmax = max(max(slp) * 1.15, 1.0)
min_soc  = min(soc)
dwell_xs = [r[0] for r in rows if r[2] == "DWELL"]
wait_xs  = [r[0] for r in rows if r[2] == "WAIT"]

# ---------------- SVG 布局 ----------------
W, H = 1100, 780
ML, MR, MT = 70, 30, 46            # 边距
PW = W - ML - MR                   # 绘图区宽
# 三个面板 (y顶, 高)
P1 = (MT, 400)
P2 = (MT + 400 + 46, 120)
P3 = (MT + 400 + 46 + 120 + 46, 120)


def X(s):
    return ML + s / smax * PW


def Y1(v):
    top, h = P1
    return top + h - max(0.0, min(v, E_CAP * 1.05)) / (E_CAP * 1.05) * h


def Y2(v):
    top, h = P2
    return top + h - max(0.0, min(v, vmax)) / vmax * h


def Y3(v):
    top, h = P3
    return top + h - max(0.0, min(v, pmax)) / pmax * h


def polyline(pts, color, width=2, fill="none", opacity=1.0):
    p = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    return (f'<polyline points="{p}" fill="{fill}" stroke="{color}" '
            f'stroke-width="{width}" opacity="{opacity}" '
            f'stroke-linejoin="round"/>')


def area(pts, ybase, color, opacity=0.65):
    p = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    x0, xn = pts[0][0], pts[-1][0]
    return (f'<polygon points="{x0:.1f},{ybase:.1f} {p} {xn:.1f},{ybase:.1f}" '
            f'fill="{color}" opacity="{opacity}" stroke="none"/>')


def text(x, y, s, size=13, color="#333", anchor="start", bold=False, rot=None):
    w = ' font-weight="bold"' if bold else ""
    tr = f' transform="rotate({rot} {x} {y})"' if rot is not None else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{color}" '
            f'text-anchor="{anchor}" font-family="Microsoft YaHei, SimHei, '
            f'sans-serif"{w}{tr}>{s}</text>')


def hline(y, x0, x1, color, width=1, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x1:.1f}" y2="{y:.1f}" '
            f'stroke="{color}" stroke-width="{width}"{d}/>')


def vband(x, top, h, color, opacity=0.4, w=6):
    return (f'<rect x="{x - w/2:.1f}" y="{top}" width="{w}" height="{h}" '
            f'fill="{color}" opacity="{opacity}"/>')


svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'viewBox="0 0 {W} {H}" style="background:#fff">',
       f'<rect width="{W}" height="{H}" fill="white"/>']

# ===== 标题 =====
svg.append(text(W / 2, 26,
                f"巡视器能量仿真曲线（最低 SoC {min_soc:.1f} Wh &gt; 200，"
                f"驻留 {len(dwell_xs)} 次，跨月昼 {len(wait_xs)} 次）",
                size=17, anchor="middle", bold=True))

# ===== 面板1：SoC =====
top, h = P1
svg.append(f'<rect x="{ML}" y="{top}" width="{PW}" height="{h}" '
           f'fill="none" stroke="#999"/>')
# 网格+刻度
for v in (0, 200, 400, 600, 800, 1000):
    y = Y1(v)
    svg.append(hline(y, ML, ML + PW, "#e5e5e5", 1))
    svg.append(text(ML - 8, y + 4, str(v), size=11, anchor="end"))
# 驻留/等待竖带（画在曲线下层）
for x in dwell_xs:
    svg.append(vband(X(x), top, h, "orange", 0.45))
for x in wait_xs:
    svg.append(vband(X(x), top, h, "purple", 0.35))
# 600 灰虚线 / 200 红虚线
svg.append(hline(Y1(SOC_DWELL), ML, ML + PW, "#888", 1, "3,4"))
svg.append(hline(Y1(E_MIN), ML, ML + PW, "red", 2.2, "8,5"))
# SoC 曲线
svg.append(polyline([(X(s), Y1(v)) for s, v in zip(ss, soc)], "#1f6fb2", 2.2))
svg.append(text(ML - 46, top + h / 2, "SoC (Wh)", size=13,
                anchor="middle", rot=-90))
# 图例
lx, ly = ML + 14, top + 16
svg.append(f'<rect x="{lx-8}" y="{ly-14}" width="270" height="96" '
           f'fill="white" opacity="0.85" stroke="#ccc"/>')
svg.append(hline(ly, lx, lx + 28, "#1f6fb2", 3))
svg.append(text(lx + 34, ly + 4, "电池电量 SoC", size=12))
svg.append(hline(ly + 20, lx, lx + 28, "red", 2.2, "8,5"))
svg.append(text(lx + 34, ly + 24, "E_min = 200 Wh（热控生存底线）", size=12))
svg.append(vband(lx + 14, ly + 32, 14, "orange", 0.5, 10))
svg.append(text(lx + 34, ly + 44, f"驻留充电 ×{len(dwell_xs)}", size=12))
svg.append(vband(lx + 14, ly + 52, 14, "purple", 0.4, 10))
svg.append(text(lx + 34, ly + 64, f"跨月昼等待 ×{len(wait_xs)}", size=12))

# ===== 面板2：光照率 =====
top, h = P2
svg.append(f'<rect x="{ML}" y="{top}" width="{PW}" height="{h}" '
           f'fill="none" stroke="#999"/>')
svg.append(area([(X(s), Y2(v)) for s, v in zip(ss, vis)], top + h, "#f4b942"))
svg.append(hline(Y2(0.3), ML, ML + PW, "#b26f1f", 1.2, "5,4"))
for v in (0, 0.5, 1.0):
    svg.append(text(ML - 8, Y2(v) + 4, f"{v:g}", size=11, anchor="end"))
svg.append(text(ML - 46, top + h / 2, "光照率", size=13,
                anchor="middle", rot=-90))
svg.append(text(ML + PW - 6, Y2(0.3) - 5, "0.3 行驶/充电阈值",
                size=11, color="#b26f1f", anchor="end"))

# ===== 面板3：坡度 =====
top, h = P3
svg.append(f'<rect x="{ML}" y="{top}" width="{PW}" height="{h}" '
           f'fill="none" stroke="#999"/>')
svg.append(area([(X(s), Y3(v)) for s, v in zip(ss, slp)], top + h, "#7aa374"))
for v in range(0, int(pmax) + 1, 2):
    svg.append(text(ML - 8, Y3(v) + 4, str(v), size=11, anchor="end"))
svg.append(text(ML - 46, top + h / 2, "坡度(°)", size=13,
                anchor="middle", rot=-90))

# X 轴刻度（公里）
step = 0.5 if smax <= 6 else 1.0
k = 0.0
while k <= smax + 1e-9:
    x = X(k)
    svg.append(hline(top + h, x, x, "#999"))
    svg.append(f'<line x1="{x:.1f}" y1="{top+h}" x2="{x:.1f}" '
               f'y2="{top+h+5}" stroke="#999"/>')
    svg.append(text(x, top + h + 20, f"{k:g}", size=11, anchor="middle"))
    k += step
svg.append(text(ML + PW / 2, top + h + 40, "沿程距离 (km)", size=13,
                anchor="middle"))

svg.append("</svg>")

os.makedirs(os.path.dirname(SVG_OUT), exist_ok=True)
with open(SVG_OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(svg))

print(f"能量曲线图（SVG）已生成 -> {SVG_OUT}")
print(f"  最低 SoC {min_soc:.1f} Wh | 驻留 {len(dwell_xs)} | 跨月昼 {len(wait_xs)}")
print("  浏览器双击打开查看；Word 可直接插入 SVG；")
print("  需要 PNG：浏览器打开 SVG 后截图或另存即可")
