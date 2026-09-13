# -*- coding: utf-8 -*-
"""
2026 CUMCM A题(药材热风干燥) 论文图表生成脚本
输出: SVG 矢量图 -> C:/Users/Admin/Desktop/CUMCM/figures/
风格: 工程热物理学术风, 低饱和配色, 图题位于图下方, 数值四位小数
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

BASE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(BASE, "支撑材料")
OUT = os.path.join(os.path.dirname(BASE), "figures")
os.makedirs(OUT, exist_ok=True)

# ---------------- 全局样式 ----------------
plt.rcParams.update({
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "svg.fonttype": "path",          # 文字转曲, 保证跨机器一致
    "pdf.fonttype": 3,               # 文字转曲线, PDF 文本提取无乱码
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#555555",
    "axes.linewidth": 0.9,
    "axes.grid": True,
    "grid.color": "#DDDDDD",
    "grid.linewidth": 0.6,
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "font.size": 10.5,
})
# 低饱和配色
C_BLUE, C_TEAL, C_ORANGE, C_RED, C_PURPLE, C_GRAY = \
    "#54789E", "#6FA191", "#D9A05B", "#C07878", "#8E84A8", "#7F7F7F"
PALETTE = [C_BLUE, C_TEAL, C_ORANGE, C_RED, C_PURPLE, C_GRAY]
CMAP_HEAT = "YlGnBu"   # 低饱和、彩色打印清晰

RADII = [0, 0.5, 1, 1.5, 2]


def style_ax(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_axisbelow(True)


def save_fig(fig, name, title):
    """图题置于图下方, 保存 SVG"""
    fig.text(0.5, 0.015, title, ha="center", va="bottom", fontsize=12)
    fig.savefig(os.path.join(OUT, name), format="svg",
                bbox_inches="tight", pad_inches=0.15)
    if name.startswith(("fig07", "fig08", "fig09", "fig10", "fig11")):
        fig.savefig(os.path.join(OUT, name.replace(".svg", ".pdf")),
                    format="pdf", bbox_inches="tight", pad_inches=0.15)
    if os.environ.get("FIG_PREVIEW"):
        prev = os.environ["FIG_PREVIEW"]
        os.makedirs(prev, exist_ok=True)
        fig.savefig(os.path.join(prev, name.replace(".svg", ".png")),
                    format="png", dpi=110, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print("saved:", name)


def read_result(fname, sheet=0):
    df = pd.read_excel(os.path.join(BASE, fname), sheet_name=sheet)
    t = df.iloc[:, 0].to_numpy(dtype=float)
    r = np.array([float(c) for c in df.columns[1:22]])
    data = df.iloc[:, 1:22].to_numpy(dtype=float)
    extra = df.iloc[:, 22].to_numpy(dtype=float) if df.shape[1] > 22 else None
    return t, r, data, extra


# ================= 图1 问题一 温度径向分布 =================
t1, r1, T1, _ = read_result("result1.xlsx", sheet="温度")
times_q1 = [100, 300, 600, 900, 1200, 1500, 1800]
t1_map = {int(v): i for i, v in enumerate(t1)}
fig, ax = plt.subplots(figsize=(7.0, 4.4))
for i, tt in enumerate(times_q1):
    row = t1_map[tt]
    ax.plot(r1, T1[row], color=PALETTE[i % 6], lw=1.7,
            marker="o", ms=3.2, label=f"t = {tt} s")
ax.set_xlabel("到药材中心的距离 r / cm")
ax.set_ylabel("温度 T / °C")
ax.legend(ncol=4, frameon=False, fontsize=9, loc="lower right")
style_ax(ax)
fig.subplots_adjust(left=0.11, right=0.97, top=0.96, bottom=0.24)
save_fig(fig, "fig01_q1_temperature_radial_profiles.svg",
         "图 1  问题一预热平衡阶段药材温度的径向分布")

# ================= 图2 问题一 水分浓度径向分布 =================
_, _, M1, _ = read_result("result1.xlsx", sheet="水分浓度")
fig, ax = plt.subplots(figsize=(7.0, 4.4))
for i, tt in enumerate(times_q1):
    row = t1_map[tt]
    ax.plot(r1, M1[row], color=PALETTE[i % 6], lw=1.7,
            marker="s", ms=3.0, label=f"t = {tt} s")
ax.set_xlabel("到药材中心的距离 r / cm")
ax.set_ylabel("水分浓度 C / (kg/kg)")
ax.legend(ncol=4, frameon=False, fontsize=9, loc="lower left")
style_ax(ax)
fig.subplots_adjust(left=0.11, right=0.97, top=0.96, bottom=0.24)
save_fig(fig, "fig02_q1_moisture_radial_profiles.svg",
         "图 2  问题一预热平衡阶段药材水分浓度的径向分布")

# ================= 附件1 烘房条件 =================
a1 = pd.read_excel(os.path.join(SUPP, "附件1.xlsx"), header=0)
a1.columns = ["t", "Ta", "Ca"]
ta_h = a1["t"].to_numpy(dtype=float) / 3600.0
Ta = a1["Ta"].to_numpy(dtype=float)
Ca = a1["Ca"].to_numpy(dtype=float)

# ================= 图3 问题二 温度随时间演化 =================
t2, r2, T2, _ = read_result("result2.xlsx", sheet="温度")
t2_h = t2 / 3600.0
step = 60  # 每 60 s 取点
fig, ax = plt.subplots(figsize=(7.0, 4.4))
for j, rr in enumerate(RADII):
    col = int(round(rr / 0.1))
    ax.plot(t2_h[::step], T2[::step, col], color=PALETTE[j], lw=1.7,
            label=f"r = {rr:g} cm")
ax.plot(ta_h, Ta, color="#444444", lw=1.4, ls="--", label="烘房温度 $T_a$")
ax.set_xlabel("时间 t / h")
ax.set_ylabel("温度 T / °C")
ax.set_xlim(0, 3)
ax.legend(ncol=3, frameon=False, fontsize=9)
style_ax(ax)
fig.subplots_adjust(left=0.11, right=0.97, top=0.96, bottom=0.24)
save_fig(fig, "fig03_q2_temperature_time_evolution.svg",
         "图 3  问题二药材温度随时间的演化（前 3 h）")

# ================= 图4 问题二 水分浓度随时间演化 =================
_, _, M2, _ = read_result("result2.xlsx", sheet="水分浓度")
fig, ax = plt.subplots(figsize=(7.0, 4.4))
for j, rr in enumerate(RADII):
    col = int(round(rr / 0.1))
    ax.plot(t2_h[::step], M2[::step, col], color=PALETTE[j], lw=1.7,
            label=f"r = {rr:g} cm")
ax.plot(ta_h, Ca, color="#444444", lw=1.4, ls="--", label="烘房水分浓度 $C_a$")
ax.set_xlabel("时间 t / h")
ax.set_ylabel("水分浓度 C / (kg/kg)")
ax.set_xlim(0, 3)
ax.legend(ncol=3, frameon=False, fontsize=9)
style_ax(ax)
fig.subplots_adjust(left=0.11, right=0.97, top=0.96, bottom=0.24)
save_fig(fig, "fig04_q2_moisture_time_evolution.svg",
         "图 4  问题二药材水分浓度随时间的演化（前 3 h）")

# ================= 图5 问题三 全程水分浓度演化 =================
t3, r3, M3, _ = read_result("result3.xlsx", sheet=0)
t3_h = t3 / 3600.0
TSTAR3 = 57.1392
fig, ax = plt.subplots(figsize=(7.2, 4.5))
for j, rr in enumerate(RADII):
    col = int(round(rr / 0.1))
    ax.plot(t3_h, M3[:, col], color=PALETTE[j], lw=1.6,
            label=f"r = {rr:g} cm")
ax.axhline(0.15, color=C_RED, lw=1.2, ls="--", label="烘干判据 0.1500 kg/kg")
ax.axvline(TSTAR3, color=C_GRAY, lw=1.2, ls=":")
ax.annotate(f"t* = {TSTAR3:.4f} h", xy=(TSTAR3, 0.15), xytext=(TSTAR3 - 17, 0.62),
            fontsize=10, color="#444444",
            arrowprops=dict(arrowstyle="->", color="#444444", lw=0.9))
ax.set_xlabel("时间 t / h")
ax.set_ylabel("水分浓度 C / (kg/kg)")
ax.set_xlim(0, 60)
ax.legend(ncol=3, frameon=False, fontsize=9, loc="upper right")
style_ax(ax)
fig.subplots_adjust(left=0.11, right=0.97, top=0.96, bottom=0.24)
save_fig(fig, "fig05_q3_moisture_time_evolution.svg",
         "图 5  问题三药材水分浓度随时间的全程演化")

# ================= 图6 问题三 水分浓度时空分布热图 =================
tstep = 30  # 每 30 行(1800 s)取样
tt = t3_h[::tstep]
rr = r3
ZZ = M3[::tstep, :].T  # (r, t)
from matplotlib.colors import PowerNorm
fig, ax = plt.subplots(figsize=(7.2, 4.5))
mesh = ax.pcolormesh(tt, rr, ZZ, cmap=CMAP_HEAT, shading="auto",
                     norm=PowerNorm(gamma=0.5, vmin=float(np.nanmin(ZZ)),
                                    vmax=float(np.nanmax(ZZ))))
cb = fig.colorbar(mesh, ax=ax, pad=0.02)
cb.set_label("水分浓度 C / (kg/kg)")
cb.outline.set_visible(False)
ax.set_xlabel("时间 t / h")
ax.set_ylabel("到药材中心的距离 r / cm")
style_ax(ax)
ax.grid(False)
fig.subplots_adjust(left=0.11, right=0.94, top=0.96, bottom=0.24)
save_fig(fig, "fig06_q3_moisture_spatiotemporal_heatmap.svg",
         "图 6  问题三药材水分浓度的时空分布")

# ================= 图7 问题四 半径收缩曲线 =================
a2 = pd.read_excel(os.path.join(SUPP, "附件2.xlsx"), header=0)
a2.columns = ["t", "R"]
tr_h = a2["t"].to_numpy(dtype=float) / 3600.0
Rcm = a2["R"].to_numpy(dtype=float)
fig, ax = plt.subplots(figsize=(7.0, 4.4))
tg = np.linspace(tr_h[0], tr_h[-1], 1200)
Rp = PchipInterpolator(tr_h, Rcm)(tg)
ax.plot(tg, Rp, color=C_BLUE, lw=1.8, label="Pchip 插值 $R(t)$")
ax.plot(tr_h, Rcm, ls="none", marker="o", ms=2.6, mfc="white",
        mec=C_BLUE, mew=0.8, label="附件 2 实测半径")
ax.annotate(f"R(0) = {Rcm[0]:.4f} cm", xy=(0, Rcm[0]), xytext=(6, 1.90),
            fontsize=10, color="#444444",
            arrowprops=dict(arrowstyle="->", color="#444444", lw=0.9))
ax.annotate(f"R(72 h) = {Rcm[-1]:.4f} cm", xy=(tr_h[-1], Rcm[-1]),
            xytext=(48, 1.38), fontsize=10, color="#444444",
            arrowprops=dict(arrowstyle="->", color="#444444", lw=0.9))
ax.set_xlabel("时间 t / h")
ax.set_ylabel("药材半径 R / cm")
ax.set_xlim(0, 72)
ax.legend(frameon=False, fontsize=9, loc="upper right")
style_ax(ax)
fig.subplots_adjust(left=0.11, right=0.97, top=0.96, bottom=0.24)
save_fig(fig, "fig07_q4_radius_shrinkage_curve.svg",
         "图 7  问题四药材半径随时间的收缩曲线")

# ================= 图8 问题四 含收缩水分浓度演化 =================
t4, r4, M4, surf4 = read_result("result4.xlsx", sheet=0)
t4_h = t4 / 3600.0
TSTAR4 = 50.8049
fig, ax = plt.subplots(figsize=(7.2, 4.5))
for j, rr in enumerate([0, 0.5, 1]):
    col = int(round(rr / 0.1))
    ax.plot(t4_h, M4[:, col], color=PALETTE[j], lw=1.6,
            label=f"r = {rr:g} cm")
ax.plot(t4_h, surf4, color=C_PURPLE, lw=1.6, label="药材表面")
ax.axhline(0.15, color=C_RED, lw=1.2, ls="--", label="烘干判据 0.1500 kg/kg")
ax.axvline(TSTAR4, color=C_GRAY, lw=1.2, ls=":")
ax.annotate(f"t* = {TSTAR4:.4f} h", xy=(TSTAR4, 0.15), xytext=(TSTAR4 - 16, 0.62),
            fontsize=10, color="#444444",
            arrowprops=dict(arrowstyle="->", color="#444444", lw=0.9))
ax.set_xlabel("时间 t / h")
ax.set_ylabel("水分浓度 C / (kg/kg)")
ax.set_xlim(0, 54)
ax.legend(ncol=3, frameon=False, fontsize=9, loc="upper right")
style_ax(ax)
fig.subplots_adjust(left=0.11, right=0.97, top=0.96, bottom=0.24)
save_fig(fig, "fig08_q4_moisture_time_evolution.svg",
         "图 8  问题四含收缩药材水分浓度随时间的全程演化")

# ================= 图9 物性×几何 2×2 消融对照 =================
groups = ["附录 3 物性", "附录 4 物性"]
fixed = [57.0884, 129.0164]
shrink = [25.1162, 50.8049]
x = np.arange(2)
w = 0.32
fig, ax = plt.subplots(figsize=(6.6, 4.4))
b1 = ax.bar(x - w / 2, fixed, w, color=C_BLUE, label="固定半径 R = 2 cm")
b2 = ax.bar(x + w / 2, shrink, w, color=C_ORANGE, label="收缩 R(t)（附件 2）")
for bars in (b1, b2):
    for b in bars:
        ax.annotate(f"{b.get_height():.4f}",
                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=9.5, color="#333333")
ax.set_xticks(x)
ax.set_xticklabels(groups)
ax.set_ylabel("烘干时长 t* / h")
ax.set_ylim(0, 145)
ax.legend(frameon=False, fontsize=9.5)
style_ax(ax)
fig.subplots_adjust(left=0.12, right=0.97, top=0.96, bottom=0.24)
save_fig(fig, "fig09_ablation_property_geometry.svg",
         "图 9  物性 × 几何 2×2 消融对照的烘干时长")

# ================= 图10 参数灵敏度龙卷风图 =================
params = ["$h_m$ −50%", "$T_a$ −2 °C", "$T_a$ +2 °C", "$h_m$ +50%",
          "$C_a$ +20%", "$C_a$ −20%", "$h$ −50%", "$h$ +50%"]
dq3 = [7.5577, 3.8598, -3.5296, -2.0362, 0.3607, -0.3066, 0.0800, -0.0261]
dq4 = [3.3384, 3.4280, -3.1393, -0.7469, 0.3269, -0.2430, 0.0748, -0.0237]
y = np.arange(len(params))
h = 0.36
fig, ax = plt.subplots(figsize=(7.4, 4.8))
b1 = ax.barh(y + h / 2, dq3, h, color=C_BLUE, label="问题三（基准 56.9672 h）")
b2 = ax.barh(y - h / 2, dq4, h, color=C_ORANGE, label="问题四（基准 50.7713 h）")
for bars in (b1, b2):
    for b in bars:
        v = b.get_width()
        off = 0.12 if v >= 0 else -0.12
        ax.annotate(f"{v:+.4f}", xy=(v, b.get_y() + b.get_height() / 2),
                    xytext=(off, 0), textcoords="offset points",
                    ha="left" if v >= 0 else "right", va="center",
                    fontsize=8.6, color="#333333")
ax.axvline(0, color="#555555", lw=0.9)
ax.set_yticks(y)
ax.set_yticklabels(params)
ax.set_xlabel("烘干时长变化量 Δt* / h")
ax.set_xlim(-4.6, 8.9)
ax.legend(frameon=False, fontsize=9, loc="upper right")
style_ax(ax)
fig.subplots_adjust(left=0.15, right=0.95, top=0.96, bottom=0.24)
save_fig(fig, "fig10_sensitivity_tstar_tornado.svg",
         "图 10  烘干时长 t* 的参数灵敏度（问题三与问题四）")

# ================= 图11 网格收敛 =================
dr_q3 = [0.5, 0.25, 0.125, 0.0625, 0.03125]
ts_q3 = [56.7064, 56.9672, 57.0874, 57.1392, 57.1597]
nz_q4 = [20, 40, 80, 160, 320]
ts_q4 = [50.5416, 50.6928, 50.7713, 50.8049, 50.8170]
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
axL, axR = axes
axL.plot(dr_q3, ts_q3, color=C_BLUE, lw=1.7, marker="o", ms=4.5)
for xx, yy in zip(dr_q3, ts_q3):
    axL.annotate(f"{yy:.4f}", xy=(xx, yy), xytext=(0, 7),
                 textcoords="offset points", ha="center",
                 fontsize=8.4, color="#333333")
axL.set_xscale("log", base=2)
axL.set_xticks(dr_q3)
axL.set_xticklabels([f"{v:g}" for v in dr_q3])
axL.invert_xaxis()  # 细化程度向右递增, 与右图一致
axL.set_xlabel("空间步长 Δr / mm")
axL.set_ylabel("烘干时长 t* / h")
axL.set_title("问题三（固定半径）", fontsize=10.5, color="#333333")
axL.set_ylim(56.4, 57.45)
axR.plot(nz_q4, ts_q4, color=C_TEAL, lw=1.7, marker="s", ms=4.5)
for xx, yy in zip(nz_q4, ts_q4):
    axR.annotate(f"{yy:.4f}", xy=(xx, yy), xytext=(0, 7),
                 textcoords="offset points", ha="center",
                 fontsize=8.4, color="#333333")
axR.set_xscale("log", base=2)
axR.set_xticks(nz_q4)
axR.set_xticklabels([f"1/{v}" for v in nz_q4])
axR.set_xlabel("随体坐标步长 Δζ")
axR.set_ylabel("烘干时长 t* / h")
axR.set_title("问题四（含收缩）", fontsize=10.5, color="#333333")
axR.set_ylim(50.3, 51.05)
for a in axes:
    style_ax(a)
fig.subplots_adjust(left=0.09, right=0.98, top=0.90, bottom=0.26, wspace=0.28)
save_fig(fig, "fig11_grid_convergence_tstar.svg",
         "图 11  烘干时长 t* 随空间离散加密的收敛情况")

print("ALL DONE")
