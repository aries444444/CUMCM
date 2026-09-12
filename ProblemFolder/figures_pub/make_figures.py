# -*- coding: utf-8 -*-
"""
CUMCM 2026 A 题（中药材热风干燥）论文配套科研级静态图
依据 alterlab-scientific-viz 出版规范：Okabe-Ito 色盲友好调色板、
无 chart junk、矢量 PDF + 300dpi PNG 双格式导出。
数据：result1-4.xlsx（四位小数）、附件1/2.xlsx、论文表 7。
运行：python make_figures.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------- 路径 ----------------
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)  # ProblemFolder
ATT = os.path.join(ROOT, 'OriginalMaterial', 'A题', '附件')
OUT = BASE  # 图片输出到本脚本同目录

# ---------------- 出版级样式（参照 alterlab-scientific-viz publication 预设） ----------------
OKABE_ITO = ['#0072B2', '#E69F00', '#009E73', '#D55E00',
             '#56B4E9', '#CC79A7', '#F0E442', '#000000']
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'Arial'],
    'axes.unicode_minus': False,
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'xtick.labelsize': 8.5,
    'ytick.labelsize': 8.5,
    'legend.fontsize': 8.5,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'lines.linewidth': 1.4,
    'axes.prop_cycle': mpl.cycler(color=OKABE_ITO),
    'figure.dpi': 110,
    'savefig.dpi': 300,
})


def despine(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def save(fig, name):
    png = os.path.join(OUT, name + '.png')
    pdf = os.path.join(OUT, name + '.pdf')
    fig.savefig(png, dpi=300, bbox_inches='tight', pad_inches=0.08, facecolor='white')
    fig.savefig(pdf, bbox_inches='tight', pad_inches=0.08, facecolor='white')
    plt.close(fig)
    print('saved', name)


def read_result(fname):
    """返回 (t_s 数组, dist_cm 数组, 数据矩阵[t, r], 原始df)"""
    df = pd.read_excel(os.path.join(ROOT, fname))
    t = df.iloc[:, 0].to_numpy(dtype=float)
    dists, cols = [], []
    for c in df.columns[1:]:
        try:
            dists.append(float(c))
            cols.append(c)
        except (TypeError, ValueError):
            pass  # 非数值列（如“表面”）单独处理
    M = df[cols].to_numpy(dtype=float)
    return t, np.array(dists), M, df


def panel_label(ax, s):
    ax.text(-0.14, 1.04, s, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top')


# ============ 图1 附件1：烘房温度与水分浓度边界 ============
def fig01_boundary():
    df = pd.read_excel(os.path.join(ATT, '附件1.xlsx'))
    t_h = df.iloc[:, 0].to_numpy(float) / 3600.0
    Ta = df.iloc[:, 1].to_numpy(float)
    Ca = df.iloc[:, 2].to_numpy(float)

    fig, ax1 = plt.subplots(figsize=(5.6, 3.2))
    ax1.plot(t_h, Ta, color=OKABE_ITO[3], label='烘房温度 $T_a$')
    ax1.set_xlabel('时间 $t$ (h)')
    ax1.set_ylabel('温度 (°C)', color=OKABE_ITO[3])
    ax1.tick_params(axis='y', labelcolor=OKABE_ITO[3])
    ax1.set_xlim(0, 4)

    ax2 = ax1.twinx()
    ax2.plot(t_h, Ca, color=OKABE_ITO[0], label='烘房水分浓度 $C_a$')
    ax2.set_ylabel('水分浓度 (kg/kg)', color=OKABE_ITO[0])
    ax2.tick_params(axis='y', labelcolor=OKABE_ITO[0])

    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [l.get_label() for l in lines], frameon=False, loc='center right')
    despine(ax1)
    ax2.spines['top'].set_visible(False)
    save(fig, 'fig01_附件1_烘房边界条件')


# ============ 图2 附件2：半径收缩 ============
def fig02_shrinkage():
    df = pd.read_excel(os.path.join(ATT, '附件2.xlsx'))
    t_h = df.iloc[:, 0].to_numpy(float) / 3600.0
    R = df.iloc[:, 1].to_numpy(float)

    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    ax.plot(t_h, R, color=OKABE_ITO[2])
    ax.set_xlabel('时间 $t$ (h)')
    ax.set_ylabel('药材半径 $R(t)$ (cm)')
    ax.set_xlim(0, t_h.max())
    ax.set_ylim(1.0, 2.05)
    pct = (R[0] - R[-1]) / R[0] * 100
    ax.annotate(f'累计收缩约 {pct:.0f}%\n$R$: {R[0]:.3f} → {R[-1]:.3f} cm',
                xy=(t_h[-1], R[-1]), xytext=(0.55, 0.25), textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=0.8, color='0.3'), fontsize=9)
    despine(ax)
    save(fig, 'fig02_附件2_半径收缩曲线')


# ============ 图3 Q1：表面/中心温度与水分演化 ============
def fig03_q1_evolution():
    t, d, T, _ = read_result('result1.xlsx')
    xl = pd.ExcelFile(os.path.join(ROOT, 'result1.xlsx'))
    dfC = xl.parse('水分浓度')
    C = dfC[[c for c in dfC.columns[1:] if _is_num(c)]].to_numpy(float)
    t_min = t / 60.0

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.1))
    i_c, i_s = 0, -1
    a1.plot(t_min, T[:, i_c], color=OKABE_ITO[0], label='中心 $r=0$')
    a1.plot(t_min, T[:, i_s], color=OKABE_ITO[3], label='表面 $r=2$ cm')
    a1.set_xlabel('时间 $t$ (min)')
    a1.set_ylabel('温度 (°C)')
    a1.annotate(f'表面 30 min: {T[-1, i_s]:.2f} °C', xy=(t_min[-1], T[-1, i_s]),
                xytext=(0.30, 0.55), textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=0.8, color='0.3'), fontsize=8.5)
    a1.legend(frameon=False, loc='lower right')
    despine(a1)

    a2.plot(t_min, C[:, i_c], color=OKABE_ITO[0], label='中心 $r=0$')
    a2.plot(t_min, C[:, i_s], color=OKABE_ITO[1], label='表面 $r=2$ cm')
    a2.set_xlabel('时间 $t$ (min)')
    a2.set_ylabel('水分浓度 (kg/kg)')
    a2.legend(frameon=False, loc='center right')
    despine(a2)

    panel_label(a1, '(a)')
    panel_label(a2, '(b)')
    fig.tight_layout(w_pad=2.5)
    save(fig, 'fig03_问题1_温度水分演化')


def _is_num(c):
    try:
        float(c)
        return True
    except (TypeError, ValueError):
        return False


# ============ 通用时空热图 ============
def spacetime(name, t_s, d, M, t_unit_div, t_label, cbar_label, cmap, t_max=None, gamma=None):
    t = t_s / t_unit_div
    if t_max is not None:
        sel = t <= t_max
        t, M = t[sel], M[sel]
    norm = mpl.colors.PowerNorm(gamma=gamma, vmin=np.nanmin(M), vmax=np.nanmax(M)) if gamma else None
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    pcm = ax.pcolormesh(t, d, M.T, cmap=cmap, shading='auto', norm=norm)
    cb = fig.colorbar(pcm, ax=ax, pad=0.02)
    cb.set_label(cbar_label)
    cb.outline.set_linewidth(0.8)
    ax.set_xlabel(t_label)
    ax.set_ylabel('距中心距离 $r$ (cm)')
    ax.set_xlim(t.min(), t.max())
    ax.set_ylim(d.min(), d.max())
    save(fig, name)


def fig04_q1_spacetime():
    t, d, T, _ = read_result('result1.xlsx')
    spacetime('fig04_问题1_温度场时空分布', t, d, T, 60.0,
              '时间 $t$ (min)', '温度 (°C)', 'viridis')


def fig05_q2_evolution():
    t, d, T, _ = read_result('result2.xlsx')
    t_h = t / 3600.0
    picks = [0.0, 0.5, 1.0, 1.5, 2.0]
    fig, ax = plt.subplots(figsize=(5.8, 3.3))
    for k, p in enumerate(picks):
        j = int(np.argmin(np.abs(d - p)))
        lab = '中心 $r=0$' if p == 0 else (f'表面 $r={p:g}$ cm' if p == 2.0 else f'$r={p:g}$ cm')
        ax.plot(t_h, T[:, j], color=OKABE_ITO[k], label=lab)
    ax.set_xlabel('时间 $t$ (h)')
    ax.set_ylabel('温度 (°C)')
    ax.set_xlim(0, 3)
    ax.legend(frameon=False, loc='lower right')
    despine(ax)
    save(fig, 'fig05_问题2_各位置温度演化')


def fig06_q2_spacetime():
    t, d, T, _ = read_result('result2.xlsx')
    spacetime('fig06_问题2_温度场时空分布', t, d, T, 3600.0,
              '时间 $t$ (h)', '温度 (°C)', 'viridis')


# ============ Q3/Q4 含水率曲线（含判据线） ============
def moisture_curve(name, fname, tstar_h, surface_mode='2.0'):
    t, d, C, df = read_result(fname)
    t_h = t / 3600.0
    i_c = 0
    if surface_mode == 'col':
        # 取最后一个数值列作为表面
        i_s = -1
    else:
        i_s = int(np.argmin(np.abs(d - float(surface_mode))))

    fig, ax = plt.subplots(figsize=(5.8, 3.3))
    ax.plot(t_h, C[:, i_c], color=OKABE_ITO[0], label='中心 $r=0$')
    ax.plot(t_h, C[:, i_s], color=OKABE_ITO[1], label='表面')
    ax.axhline(0.15, color=OKABE_ITO[3], ls='--', lw=1.0)
    ax.text(t_h.max() * 0.02, 0.19, '烘干判据 $C=0.15$ kg/kg',
            color=OKABE_ITO[3], fontsize=8.5)
    ax.axvline(tstar_h, color='0.35', ls=':', lw=1.0)
    ax.annotate(f'$t^*={tstar_h:.2f}$ h', xy=(tstar_h, 0.15),
                xytext=(0.62, 0.30), textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=0.8, color='0.3'), fontsize=9)
    ax.set_xlabel('时间 $t$ (h)')
    ax.set_ylabel('水分浓度 $C$ (kg/kg)')
    ax.set_xlim(0, t_h.max() * 1.02)
    ax.set_ylim(0, None)
    ax.legend(frameon=False, loc='upper right')
    despine(ax)
    save(fig, name)


def fig07_q3():
    t, d, C, _ = read_result('result3.xlsx')
    tstar_h = t[-1] / 3600.0
    moisture_curve('fig07_问题3_含水率演化与烘干判据', 'result3.xlsx', tstar_h, '2.0')


def fig08_q3_spacetime():
    t, d, C, _ = read_result('result3.xlsx')
    spacetime('fig08_问题3_含水率场时空分布', t, d, C, 3600.0,
              '时间 $t$ (h)', '水分浓度 (kg/kg)', 'YlGnBu_r', gamma=0.45)


def fig09_q4():
    t, d, C, _ = read_result('result4.xlsx')
    tstar_h = t[-1] / 3600.0
    moisture_curve('fig09_问题4_含收缩含水率演化', 'result4.xlsx', tstar_h, 'col')


def fig10_q4_spacetime():
    t, d, C, _ = read_result('result4.xlsx')
    spacetime('fig10_问题4_含收缩含水率场时空分布', t, d, C, 3600.0,
              '时间 $t$ (h)', '水分浓度 (kg/kg)', 'YlGnBu_r', gamma=0.45)


# ============ 图11 Q3 vs Q4 中心含水率对比 ============
def fig11_compare():
    t3, _, C3, _ = read_result('result3.xlsx')
    t4, _, C4, _ = read_result('result4.xlsx')
    t3h, t4h = t3 / 3600.0, t4 / 3600.0

    fig, ax = plt.subplots(figsize=(5.8, 3.3))
    ax.plot(t3h, C3[:, 0], color=OKABE_ITO[0], label='问题3：不考虑收缩')
    ax.plot(t4h, C4[:, 0], color=OKABE_ITO[3], label='问题4：考虑收缩')
    ax.axhline(0.15, color='0.35', ls='--', lw=1.0)
    ax.text(1.0, 0.20, '烘干判据 $C=0.15$ kg/kg', color='0.35', fontsize=8.5)
    for tt, cc, col in ((t3h[-1], C3[-1, 0], OKABE_ITO[0]), (t4h[-1], C4[-1, 0], OKABE_ITO[3])):
        ax.plot([tt], [cc], 'o', ms=5, color=col)
        ax.annotate(f'{tt:.2f} h', xy=(tt, cc), xytext=(tt - 8.5, cc + 0.28),
                    fontsize=8.5, color=col)
    shrink = (t3h[-1] - t4h[-1]) / t3h[-1] * 100
    ax.text(0.34, 0.60, f'收缩使烘干时长缩短约 {shrink:.1f}%',
            transform=ax.transAxes, fontsize=9)
    ax.set_xlabel('时间 $t$ (h)')
    ax.set_ylabel('中心水分浓度 $C(0,t)$ (kg/kg)')
    ax.set_xlim(0, t3h.max() * 1.02)
    ax.set_ylim(0, None)
    ax.legend(frameon=False, loc='upper right')
    despine(ax)
    save(fig, 'fig11_问题3与4_中心含水率对比')


# ============ 图12 参数灵敏度（表7） ============
def fig12_sensitivity():
    labels = [r'$T_a+2$°C', r'$T_a-2$°C', r'$C_a+20\%$', r'$C_a-20\%$',
              r'$h+50\%$', r'$h-50\%$', r'$h_m+50\%$', r'$h_m-50\%$']
    q3 = np.array([-3.44, 3.77, 0.34, -0.29, -0.03, 0.08, -2.04, 7.56])
    q4 = np.array([-3.04, 3.32, 0.32, -0.24, -0.02, 0.08, -0.75, 3.34])
    order = np.argsort(np.abs(q3))[::-1]

    y = np.arange(len(labels))
    h = 0.38
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.barh(y - h / 2, q3[order], height=h, color=OKABE_ITO[0], label='问题3')
    ax.barh(y + h / 2, q4[order], height=h, color=OKABE_ITO[1], label='问题4')
    ax.set_yticks(y)
    ax.set_yticklabels([labels[i] for i in order])
    ax.invert_yaxis()
    ax.axvline(0, color='0.2', lw=0.8)
    for yi, i in enumerate(order):
        for v, off in ((q3[i], -h / 2), (q4[i], h / 2)):
            ax.text(v + (0.15 if v >= 0 else -0.15), yi + off, f'{v:+.2f}',
                    va='center', ha='left' if v >= 0 else 'right', fontsize=7.5, color='0.25')
    ax.set_xlabel('烘干时长变化 $\\Delta t^*$ (h)')
    ax.set_xlim(-4.6, 9.2)
    ax.legend(frameon=False, loc='lower right')
    despine(ax)
    fig.tight_layout()
    save(fig, 'fig12_烘干时长参数灵敏度')


if __name__ == '__main__':
    fig01_boundary()
    fig02_shrinkage()
    fig03_q1_evolution()
    fig04_q1_spacetime()
    fig05_q2_evolution()
    fig06_q2_spacetime()
    fig07_q3()
    fig08_q3_spacetime()
    fig09_q4()
    fig10_q4_spacetime()
    fig11_compare()
    fig12_sensitivity()
    print('ALL DONE')
