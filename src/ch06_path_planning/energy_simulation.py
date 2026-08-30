# -*- coding: utf-8 -*-
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch06 · 路径规划 阶段三：巡视器三态状态机能量递推（DRIVE/DWELL/WAIT，SoC≥200Wh）
# 来源  : 竞赛提交包 3工程文件/Ch06_*/energy_simulation.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
energy_simulation.py (v2) —— 第六章 阶段三：巡视器三态状态机能量递推

输入：path_profile.csv（阶段二产出，82 点 @50m）
输出：表 6-3 指标 + energy_curve.csv + energy_curve.png（200Wh 红线）

三态状态机（v2，物理修正版）：
  DRIVE：SoC 足够 → 行驶一步（Δs=50m，Δt=0.5h）
         · 亮区（visib≥0.3）：行驶阈值 SoC>600（保留工程储备）
         · 暗区（visib<0.3）：允许电池推进至热控底线 E_min=200
           —— 与论文 PSR 边缘论证同源（840Wh 可用电量 ÷ 220W 净功耗）
  DWELL：亮区 SoC≤600 → 驻留充电至 900 恢复（计 1 次驻留）
  WAIT ：暗区电量不足以再走一步 → 原地跨月昼等待（计 1 个月昼），
         利用下一月昼的累计光照满充至 E_cap。
         物理依据：AVGVISIB 为整太阴月光照时间占比，
         驻留一个月昼可积累 P_solar × visib × 708h；
         即使 visib=0.01 也有 ~991Wh > 840Wh 可用容量，
         仅真 PSR（<0.001）无法补能（终点即边界，到达即止）。

能量递推（下界=E_min，SoC=0 问题的根本修复）：
  E_{k+1} = max(E_min, min(E_cap, E_k + ΔE_charge − ΔE_discharge))
  充电：ΔE_charge = P_solar · visib · Δt（行驶中同步充电）
  耗电：ΔE_discharge = [P_base + P_drive·(1 + sinθ·2.5 + k_vrm·VRM)]·Δt
        下坡（θ<0）：行驶功耗降至最低安全代偿 30W

运行：普通 Python 3.10+ 直接运行（纯标准库，无第三方依赖）；
      环境无 matplotlib 时自动只出 CSV，用 plot_energy_curve_svg.py 零依赖出图
"""

import os
import csv
import math

# ============================= 配置区 =============================

_base        = os.environ.get("LUNAR_OUTPUT_DIR", os.path.join("data", "output"))
PROFILE_CSV = os.path.join(_base, "stage2", "path_profile.csv")
OUT_DIR     = os.path.join(_base, "stage3_energy")

E_MIN     = 200.0     # Wh 热控生存底线（递推下界）
E_CAP     = 1040.0    # Wh SoC 上限 = E_MIN + E_usable(840)
E_START   = 1040.0    # Wh 出发满电
SOC_DWELL = 600.0     # Wh 亮区驻留充电阈值
SOC_RESUME= 900.0     # Wh 驻留充至此恢复
P_BASE    = 120.0     # W  基础功耗
P_DRIVE   = 100.0     # W  平地行驶功耗
P_DOWNHILL= 30.0      # W  下坡最低安全代偿
P_STANDBY = 30.0      # W  驻留基础功耗
P_SOLAR   = 140.0     # W  太阳能板等效输出
V_KMH     = 0.1       # km/h
K_VRM     = 2.5       # VRM 阻力系数
VISIB_CHG = 0.3       # 亮区判定（可即时充电）
LUNAR_LIT_H = 708.0   # 一个太阴月小时数（29.5天），WAIT 补能积分用

REQ_MIN_SOC   = 200.0
REQ_MAX_DWELL = 15
REQ_MAX_LUNAR = 10

# ==================================================================


def load_profile(path):
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            rows.append(dict(
                s=float(r["s_m"]), x=float(r["x"]), y=float(r["y"]),
                elev=float(r["elev"]), slope=float(r["slope"]),
                vrm=float(r["vrm"]), visib=float(r["visib"])))
    if len(rows) < 2:
        raise RuntimeError("剖面表少于 2 点")
    return rows


def drive_power(theta_deg, vrm):
    if theta_deg < 0:
        return P_BASE + P_DOWNHILL
    return P_BASE + P_DRIVE * (1.0 + math.sin(math.radians(theta_deg)) * 2.5
                               + K_VRM * vrm)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    prof = load_profile(PROFILE_CSV)
    n = len(prof)
    print("=" * 64)
    print("阶段三：巡视器能量仿真（三态状态机递推 v2）")
    print("=" * 64)
    dark = sum(1 for p in prof if p["visib"] < VISIB_CHG)
    print(f"剖面：{n} 点，总长 {prof[-1]['s']/1000:.3f} km，"
          f"平均光照 {sum(p['visib'] for p in prof)/n:.4f}，"
          f"暗区(<{VISIB_CHG:g})占 {dark}/{n} 点")

    soc = E_START
    total_dis = total_chg = 0.0
    min_soc, min_soc_at = soc, 0.0
    dwell_cnt, lunar_days = 0, 1
    log = [(0.0, soc, "START", prof[0]["visib"], prof[0]["slope"])]

    i = 1
    guard = 0
    while i < n:
        guard += 1
        if guard > 20000:
            raise RuntimeError("迭代超限，检查参数")
        p0, p1 = prof[i - 1], prof[i]
        ds = p1["s"] - p0["s"]
        dt = (ds / 1000.0) / V_KMH
        theta = math.degrees(math.atan2(p1["elev"] - p0["elev"], max(ds, 1e-6)))
        visib0 = p0["visib"]          # 当前位置光照（决定能否充电）
        visib1 = p1["visib"]          # 目标点光照（行驶中充电按沿途均值近似）
        v_avg = (visib0 + visib1) / 2.0

        p_use = drive_power(theta, p1["vrm"])
        e_dis = p_use * dt
        e_chg = P_SOLAR * v_avg * dt
        lit_here = visib0 >= VISIB_CHG

        # ---- 能否走这一步？ ----
        floor = SOC_DWELL if lit_here else E_MIN
        if soc + e_chg - e_dis >= floor or soc >= E_CAP * 0.999 and not lit_here:
            can_drive = (soc + e_chg - e_dis) >= (SOC_DWELL if lit_here else E_MIN)
        can_drive = (soc + e_chg - e_dis) >= (SOC_DWELL if lit_here else E_MIN)

        if can_drive:
            total_dis += e_dis
            total_chg += e_chg
            soc = max(E_MIN, min(E_CAP, soc + e_chg - e_dis))
            if soc < min_soc:
                min_soc, min_soc_at = soc, p1["s"]
            log.append((p1["s"], soc, "DRIVE", visib1, p1["slope"]))
            i += 1
            continue

        if lit_here:
            # DWELL：亮区驻留充电 600→900
            net = P_SOLAR * visib0 - P_STANDBY
            if net > 0:
                need = SOC_RESUME - soc
                t_chg = need / net
                total_chg += P_SOLAR * visib0 * t_chg
                total_dis += P_STANDBY * t_chg
                soc = SOC_RESUME
                dwell_cnt += 1
                log.append((p0["s"], soc, "DWELL", visib0, p0["slope"]))
                continue
            # 亮区但净充电≤0（不应发生），退化为 WAIT

        # WAIT：暗区（或无法正充）跨月昼补能
        e_lunar = P_SOLAR * visib0 * LUNAR_LIT_H     # 下一月昼可积累能量
        if e_lunar < (e_dis - e_chg) and visib0 < 0.001:
            # 真 PSR 内无法补能——按设计不应发生（终点=边界即止）
            print(f"  [警告] s={p0['s']:.0f}m visib={visib0:.4g} 处无法补能，"
                  f"任务在此终止条件下不可行（检查路径是否深入 PSR）")
            break
        lunar_days += 1
        soc = min(E_CAP, max(soc, E_MIN) + min(e_lunar, E_CAP))
        soc = min(E_CAP, soc)
        total_chg += min(e_lunar, E_CAP - E_MIN)
        log.append((p0["s"], soc, "WAIT", visib0, p0["slope"]))

    # ---------------- 输出 ----------------
    net_e = total_chg - total_dis
    print("\n── 表 6-3 能量收支指标 ──")
    print(f"  单程总耗能      {total_dis:,.1f} Wh")
    print(f"  单程太阳能补给  {total_chg:,.1f} Wh")
    print(f"  净能量差        {net_e:,.1f} Wh")
    print(f"  最低 SoC        {min_soc:.1f} Wh（里程 {min_soc_at/1000:.2f} km 处）")
    print(f"  驻留充电        {dwell_cnt} 次")
    print(f"  跨月昼等待      {lunar_days - 1} 次 → 任务周期 {lunar_days} 个月昼（单程）")

    print("\n── 验收判定 ──")
    ok1 = min_soc > REQ_MIN_SOC
    ok1b = min_soc >= REQ_MIN_SOC
    tag1 = "✓" if ok1 else ("△" if ok1b else "✗")
    print(f"  [{tag1}] 全程最低 SoC > {REQ_MIN_SOC:g} Wh（实测 {min_soc:.1f}）"
          + ("" if ok1 else "（=下界即触底保护，答辩需说明暗区推进段贴底线运行）"
             if ok1b else ""))
    ok2 = dwell_cnt < REQ_MAX_DWELL
    print(f"  [{'✓' if ok2 else '✗'}] 驻留充电 < {REQ_MAX_DWELL} 次（实测 {dwell_cnt}）")
    ok3 = lunar_days < REQ_MAX_LUNAR
    print(f"  [{'✓' if ok3 else '✗'}] 单程月昼 < {REQ_MAX_LUNAR} 个（实测 {lunar_days}）")

    curve_csv = os.path.join(OUT_DIR, "energy_curve.csv")
    with open(curve_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["s_m", "soc_wh", "state", "visib", "slope_deg"])
        for row in log:
            w.writerow([f"{row[0]:.1f}", f"{row[1]:.1f}", row[2],
                        f"{row[3]:.4f}", f"{row[4]:.2f}"])
    print(f"\n曲线数据 -> {curve_csv}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        for fam in ["SimHei", "Microsoft YaHei", "Noto Sans CJK SC"]:
            try:
                import matplotlib.font_manager as fm
                if any(fam.lower() in f.name.lower() for f in fm.fontManager.ttflist):
                    plt.rcParams["font.sans-serif"] = [fam]
                    break
            except Exception:
                pass
        plt.rcParams["axes.unicode_minus"] = False
        ss = [l[0] / 1000 for l in log]
        es = [l[1] for l in log]
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.plot(ss, es, lw=1.6, color="#1f6fb2", label="SoC")
        ax.axhline(200, color="red", lw=1.5, ls="--", label="E_min = 200 Wh")
        ax.axhline(SOC_DWELL, color="#888", lw=0.8, ls=":", label="600 Wh")
        first_dw = first_wt = True
        for l in log:
            if l[2] == "DWELL":
                ax.axvline(l[0] / 1000, color="orange", alpha=0.4, lw=3,
                           label="DWELL" if first_dw else None)
                first_dw = False
            elif l[2] == "WAIT":
                ax.axvline(l[0] / 1000, color="purple", alpha=0.35, lw=3,
                           label="WAIT (next lunar day)" if first_wt else None)
                first_wt = False
        ax.set_xlabel("distance (km)")
        ax.set_ylabel("SoC (Wh)")
        ax.set_title(f"Rover energy simulation  (min SoC {min_soc:.0f} Wh, "
                     f"dwell x{dwell_cnt}, {lunar_days} lunar days)")
        ax.set_ylim(0, E_CAP * 1.05)
        ax.legend(loc="lower left", fontsize=9)
        ax.grid(alpha=0.3)
        png = os.path.join(OUT_DIR, "energy_curve.png")
        fig.savefig(png, dpi=150, bbox_inches="tight")
        print(f"能量曲线图 -> {png}（含 200 Wh 红线）")
    except ImportError:
        print("[提示] 无 matplotlib，未出图；energy_curve.csv 可用 Excel 画。")

    print("=" * 64)


if __name__ == "__main__":
    main()
