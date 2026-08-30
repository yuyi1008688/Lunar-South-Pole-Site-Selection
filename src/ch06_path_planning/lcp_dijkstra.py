# -*- coding: utf-8 -*-
"""
lcp_dijkstra.py —— 8 邻域 Dijkstra 最低成本路径（纯 Python，替代商业 GIS 的
cost_distance / cost_path 算子）
================================================================
数学定义（与 ArcGIS / 商业 GIS 的耗费距离一致）：

  成本栅格 cost 给出"穿过该像元的单位距离成本"。从像元 u 移动到相邻像元 v 的
  边成本 = (cost[u] + cost[v]) / 2 × dist(u,v)，
  其中正交邻居 dist = 240 m，对角邻居 dist = 240×√2 m（"平均成本×移动距离"
  的标准耗费距离定义）。

  成本累积面 accumulation：从源像元出发做 Dijkstra，得到每个像元到达源的最小
  累积成本；回溯方向面 backlink：每个像元记录"通往源的路径上一步走哪个邻居"
  （8 方向编码）；最优路径：从任意目标像元沿 backlink 一路回溯到源。

实现要点：
  1. heapq 实现的 Dijkstra；384×384=147456 节点，秒级完成；不用 networkx；
  2. NoData/NaN/负成本像元视为不可通行（inf）；
  3. 平局处理：多个邻居累积成本相等时，固定扫描顺序（正北起顺时针）取第一个，
     保证结果确定性；等成本路径不唯一属正常现象（验证报告单独量化）；
  4. 支持任意 source → 一次实现复用于任务路径与全区几何路径。

方向编码（本模块约定，NEIGHBORS 顺序，与商业 GIS 的 backlink 同语义）：
  1=北(N) 2=东北(NE) 3=东(E) 4=东南(SE) 5=南(S) 6=西南(SW) 7=西(W) 8=西北(NW)
  方向语义 = **到达方向**："路径从前驱像元进入当前像元"的方向，即前驱位于
  当前像元的反方向（回溯到源时沿 -方向码 前进）。
"""

import heapq
import math

import numpy as np

# 邻域扫描顺序：正北起顺时针（dr, dc, 距离倍数）
#   1=北(dr=-1,dc=0) 2=东北(-1,+1) 3=东(0,+1) 4=东南(+1,+1)
#   5=南(+1,0) 6=西南(+1,-1) 7=西(0,-1) 8=西北(-1,-1)
NEIGHBORS = [(-1, 0, 1.0), (-1, 1, math.sqrt(2)), (0, 1, 1.0), (1, 1, math.sqrt(2)),
             (1, 0, 1.0), (1, -1, math.sqrt(2)), (0, -1, 1.0), (-1, -1, math.sqrt(2))]

NODATA = -9999.0


def _passable(cost):
    """可通行掩膜：非 NoData、非 NaN、成本有限且 >= 0。"""
    return np.isfinite(cost) & (cost != NODATA) & (cost >= 0)


def dijkstra_accumulation(cost, source_rc, cell_size=240.0):
    """以 source_rc 为源做 8 邻域 Dijkstra。

    参数
    ----
    cost : (H, W) array  单位距离成本栅格（NoData 用 -9999 或 NaN 标记）
    source_rc : (row, col) 源像元
    cell_size : 像元边长（米），正交步长；对角步长 = cell_size×√2

    返回
    ----
    accumulation : (H, W) float64  累积成本（cost×米；源=0，不可达=inf）
    backlink     : (H, W) int8     回溯方向码（1-8，语义=指向回溯前驱；
                                   源=0，不可达/未更新=0）
    """
    H, W = cost.shape
    passable = _passable(cost)
    sr, sc = source_rc
    if not (0 <= sr < H and 0 <= sc < W) or not passable[sr, sc]:
        raise ValueError(f"源像元不可通行：{source_rc}")

    inf = np.inf
    accumulation = np.full((H, W), inf, dtype=np.float64)
    backlink = np.zeros((H, W), dtype=np.int8)
    accumulation[sr, sc] = 0.0

    # 堆元素：(累积成本, row, col)；惰性删除（弹出时比对 accumulation）
    heap = [(0.0, sr, sc)]
    push, pop = heapq.heappush, heapq.heappop
    c_flat = cost  # 局部引用提速
    while heap:
        d, r, c = pop(heap)
        if d > accumulation[r, c]:
            continue  # 过期条目
        cu = c_flat[r, c]
        for code, (dr, dc, mult) in enumerate(NEIGHBORS, start=1):
            nr, nc = r + dr, c + dc
            if not (0 <= nr < H and 0 <= nc < W):
                continue
            if not passable[nr, nc]:
                continue
            nd = d + (cu + c_flat[nr, nc]) * 0.5 * cell_size * mult
            if nd < accumulation[nr, nc]:
                accumulation[nr, nc] = nd
                backlink[nr, nc] = code
                push(heap, (nd, nr, nc))
    return accumulation, backlink


def trace_backlink(backlink, target_rc, source_rc):
    """沿回溯方向面从 target 走到 source，返回 (row, col) 有序列表。

    方向码语义 = 到达方向（前驱→当前），故回溯时每步沿 **反方向** 前进：
    前驱 = 当前 - (dr, dc)。带步数上限保护（防死循环）。
    """
    H, W = backlink.shape
    tr, tc = target_rc
    sr, sc = source_rc
    path = [(tr, tc)]
    r, c = tr, tc
    max_steps = H * W  # 防御上限
    for _ in range(max_steps):
        if (r, c) == (sr, sc):
            return path
        code = int(backlink[r, c])
        if code == 0:
            break  # 到源（码=0）或断链
        dr, dc, _ = NEIGHBORS[code - 1]
        r, c = r - dr, c - dc  # 前驱在到达方向的反方向
        if not (0 <= r < H and 0 <= c < W):
            break
        path.append((r, c))
    return path


def path_cells_to_coords(path_rc, left=-46080.0, top=46080.0, res=240.0):
    """路径像元序列 → 像元中心地理坐标 [(x, y), ...]。"""
    return [(left + c * res + res / 2.0, top - r * res - res / 2.0) for r, c in path_rc]


def path_length_m(path_rc, cell_size=240.0):
    """像元序列的折线长度（米）：相邻像元中心距（正交 cell_size，对角 ×√2）。"""
    total = 0.0
    for (r1, c1), (r2, c2) in zip(path_rc[:-1], path_rc[1:]):
        dr, dc = r2 - r1, c2 - c1
        steps = {(-1, 0), (1, 0), (0, -1), (0, 1)}
        mult = 1.0 if (dr, dc) in steps else math.sqrt(2)
        total += cell_size * mult
    return total
