# -*- coding: utf-8 -*-
"""
bspline_smooth.py —— B 样条路径平滑 + 5→3→2 自动降级 + 危险区回检（纯 scipy 实现）
================================================================
栅格 LCP 路径只能沿 8 方向像元中心走线，天然是横平竖直 + 45° 斜角的锯齿折线，
必须用 B 样条拟合光滑曲线（替代商业 GIS 的 BSPLINE smooth_method）。

与商业 GIS 行为的对齐要点（均经 L2 验证，详见 tests/解耦与精度验证报告.md）：
  1. 阶数语义：k = 平滑阶数（5 最光滑、最易"抄近道"蹭进危险区；2 最贴栅格
     原路径、安全冗余最高）。scipy 的 k 与商业软件 smoothDegree 数值语义
     不保证等价，以结果指标反推标定，不直接套数值；
  2. 参数化：对输入像素链按累计弧长（chord length）参数化后再 splprep，
     避免均匀参数化在折线密集/稀疏段产生震荡；
  3. 输出顶点数：与商业软件"输入 n 点 → 输出 n 点"的采样约定一致，
     n_out 默认 = len(coords)（任务路径 17 像元 → 17 顶点，与历史基准一致）；
  4. 自动降级：先用 5 阶，把平滑结果按 30m 等间距重采样逐点查危险区，
     若存在 hazard ≥ 1.0 的穿越点则降 3 阶重试，再穿降到 2 阶；
     2 阶仍穿越则如实接受并报告（"0 穿越硬约束下尽量光滑"是判据，
     不是越光滑越好）。

平滑因子 s：s=0 为插值样条（过全部像素链顶点）。经标定，s=0 已同时满足
"0 穿越 + 长度/顶点最接近历史基准"（见 L2 标定记录），故默认 s=0；
如需更松的近似可调大 s（scan_s() 提供扫描工具）。
"""

import math

import numpy as np
from scipy.interpolate import splev, splprep


def chord_parameter(coords):
    """累计弧长参数化：返回各点在 [0,1] 上的参数 u（总弧长归一化）。"""
    pts = np.asarray(coords, dtype=np.float64)
    seg = np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    total = cum[-1]
    if total <= 0:
        return np.linspace(0.0, 1.0, len(coords))
    return cum / total


def bspline_smooth(coords, degree=5, n_out=None, s_param=0.0):
    """B 样条平滑一条折线。

    参数
    ----
    coords : [(x, y), ...]  像元中心坐标序列（起点→终点）
    degree : 样条阶数 k（5/3/2）
    n_out  : 输出顶点数，默认 = len(coords)（对齐商业软件采样约定）
    s_param: scipy splprep 平滑因子（0 = 插值样条）

    返回
    ----
    smoothed : [(x, y), ...]  平滑后顶点序列
    """
    pts = np.asarray(coords, dtype=np.float64)
    n = len(pts)
    if n < degree + 1:
        # 点数不足以支撑该阶数：直接降为可行最高阶
        degree = max(1, n - 1)
    if n_out is None:
        n_out = n
    u = chord_parameter(coords)
    # splprep 的 u 必须严格递增；若存在重复点导致的零段，给 u 加微小增量
    eps = 1e-12
    for i in range(1, len(u)):
        if u[i] <= u[i - 1]:
            u[i] = u[i - 1] + eps
    tck, _ = splprep([pts[:, 0], pts[:, 1]], u=u, s=s_param, k=degree)
    unew = np.linspace(0.0, 1.0, n_out)
    xn, yn = splev(unew, tck)
    return list(zip(xn.tolist(), yn.tolist()))


def resample(coords, step):
    """折线按 step 等间距重采样（保留原实现：累计长度 + 线性插值）。

    返回 ([(x, y, s_m), ...], total_length_m)；最后一点强制为折线终点。
    """
    a = list(coords)
    seg = [math.hypot(a[i + 1][0] - a[i][0], a[i + 1][1] - a[i][1])
           for i in range(len(a) - 1)]
    cum = [0.0]
    for s_ in seg:
        cum.append(cum[-1] + s_)
    total = cum[-1]
    st = [i * step for i in range(int(total / step) + 1)] + [total]
    out, j = [], 0
    for s_ in st:
        while j < len(seg) - 1 and cum[j + 1] < s_:
            j += 1
        t = 0.0 if seg[j] == 0 else (s_ - cum[j]) / seg[j]
        out.append((a[j][0] + t * (a[j + 1][0] - a[j][0]),
                    a[j][1] + t * (a[j + 1][1] - a[j][1]), s_))
    return out, total


def hazard_crossings(coords_xy, hazard_grid, haz_max=1.0, check_step=30.0):
    """平滑路径的危险区回检：30m 等间距逐点查 hazard（最近邻，判定用）。

    返回 (穿越点列表 [(x, y, hazard), ...], 检查点数, 路径总长 m)。
    """
    pts, total = resample(coords_xy, check_step)
    hits = [(x, y, hazard_grid.value_at_xy(x, y, method='nearest'))
            for x, y, _ in pts
            if hazard_grid.value_at_xy(x, y, method='nearest') >= haz_max]
    return hits, len(pts), total


def smooth_with_degradation(coords, hazard_grid, degrees=(5, 3, 2),
                            haz_max=1.0, check_step=30.0, n_out=None, s_param=0.0,
                            verbose=True):
    """5→3→2 自动降级平滑（判据：0 穿越硬约束下尽量光滑）。

    返回 (best_coords, record_list)；record 每项含 degree/vertices/length_km/
    crossings。2 阶仍穿越则接受并如实记录（record 带 accepted=False）。
    """
    records = []
    best = None
    for degree in degrees:
        smoothed = bspline_smooth(coords, degree=degree, n_out=n_out,
                                  s_param=s_param)
        hits, n_chk, total = hazard_crossings(smoothed, hazard_grid,
                                              haz_max=haz_max,
                                              check_step=check_step)
        rec = dict(degree=degree, vertices=len(smoothed),
                   length_km=total / 1000.0, crossings=len(hits),
                   check_points=n_chk, accepted=len(hits) == 0)
        records.append(rec)
        if verbose:
            print(f"    平滑{degree}: 顶点 {len(smoothed)}, "
                  f"长度 {total / 1000:.3f} km, 穿越 {len(hits)}")
        best = smoothed
        if len(hits) == 0 or degree <= min(degrees):
            break
    return best, records


def scan_s(coords, hazard_grid, s_values=(0.0, 1.0, 10.0, 100.0), degree=5,
           haz_max=1.0, check_step=30.0, n_out=None):
    """平滑因子 s 扫描标定工具：打印各 s 下的顶点数/长度/穿越数。

    返回满足 0 穿越且长度最接近 |coords| 像元折线长的 s 值。
    """
    seg = np.hypot(*(np.diff(np.asarray(coords, dtype=float), axis=0).T))
    fold_len = float(seg.sum())
    best_s, best_gap = None, None
    print(f'  s 扫描（degree={degree}, 像元折线长基准 {fold_len / 1000:.3f} km）:')
    for s_ in s_values:
        try:
            smoothed = bspline_smooth(coords, degree=degree, n_out=n_out,
                                      s_param=s_)
        except Exception as e:
            print(f'    s={s_:<8g} 拟合失败: {e}')
            continue
        hits, _, total = hazard_crossings(smoothed, hazard_grid,
                                          haz_max=haz_max,
                                          check_step=check_step)
        gap = abs(total - fold_len)
        print(f'    s={s_:<8g} 顶点 {len(smoothed):4d} 长度 {total / 1000:8.3f} km '
              f'穿越 {len(hits)}')
        if not hits and (best_gap is None or gap < best_gap):
            best_s, best_gap = s_, gap
    return best_s
