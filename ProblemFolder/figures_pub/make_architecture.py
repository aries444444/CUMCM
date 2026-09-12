# -*- coding: utf-8 -*-
"""
CUMCM 2026 A 题 论文模型总体架构图（Nature 风格）
设计依据：scientific-illustration-guide（图形摘要构图原则、色盲友好、
矢量导出、期刊合规）+ alterlab-scientific-viz（Nature 出版样式）。
运行：python make_architecture.py
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = os.path.dirname(os.path.abspath(__file__))

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'Arial'],
    'axes.unicode_minus': False,
})

# Nature 风格配色：低饱和底色 + 深色描边 + 少量强调色
C_INPUT = ('#E8F0E9', '#2E7D46')   # 输入：浅绿底/深绿边
C_MODEL = ('#E7EEF7', '#1F4E9C')   # 模型：浅蓝底/深蓝边
C_NUM = ('#F3ECF5', '#6A3D9A')     # 求解：浅紫底/深紫边
C_OUT = ('#FBEFE3', '#C05A12')     # 输出：浅橙底/深橙边
C_CHECK = ('#FBE7E7', '#B2182B')   # 检验：浅红底/深红边
TXT = '#1A1A1A'
SUB = '#444444'
ARROW = '#666666'


def box(ax, x, y, w, h, title, sub=None, colors=C_MODEL, title_fs=9,
        sub_fs=7.8, lw=1.1, radius=0.8, title_dy=None):
    fc, ec = colors
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0,rounding_size={radius}",
                       facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2)
    ax.add_patch(p)
    cx = x + w / 2
    if sub is None:
        ax.text(cx, y + h / 2, title, ha='center', va='center',
                fontsize=title_fs, fontweight='bold', color=ec, zorder=3)
    else:
        dy = title_dy if title_dy is not None else h * 0.30
        ax.text(cx, y + h / 2 + dy, title, ha='center', va='center',
                fontsize=title_fs, fontweight='bold', color=ec, zorder=3)
        ax.text(cx, y + h / 2 - h * 0.22, sub, ha='center', va='center',
                fontsize=sub_fs, color=SUB, zorder=3, linespacing=1.5)


def arrow(ax, x1, y1, x2, y2, color=ARROW, lw=1.2, style='-|>', rad=0.0,
          ls='-', ms=11):
    a = FancyArrowPatch((x1, y1), (x2, y2),
                        arrowstyle=style, mutation_scale=ms,
                        color=color, lw=lw, linestyle=ls,
                        connectionstyle=f"arc3,rad={rad}", zorder=1,
                        shrinkA=1, shrinkB=1)
    ax.add_patch(a)


def band_label(ax, y, text, color):
    ax.text(0.5, y, text, ha='left', va='center', fontsize=9.5,
            fontweight='bold', color=color, zorder=4,
            bbox=dict(facecolor='white', edgecolor='none', pad=1.6))


fig, ax = plt.subplots(figsize=(9.8, 7.8))
ax.set_xlim(0, 100)
ax.set_ylim(-13.5, 74)
ax.axis('off')

# ============ 第一层：数据与参数输入 ============
band_label(ax, 70.5, '数据与参数输入', C_INPUT[1])
box(ax, 6, 58, 26, 9.5, '附件 1：烘房边界条件',
    '烘房温度 $T_a(t)$、水分浓度 $C_a(t)$\n0–14400 s 线性插值，其后恒温段',
    C_INPUT)
box(ax, 37, 58, 26, 9.5, '附件 2：收缩半径 $R(t)$',
    '145 个测点（0–72 h）\nPchip 保形插值 + 解析导数',
    C_INPUT)
box(ax, 68, 58, 26, 9.5, '附录 2–4：物性经验式',
    '$\\rho,\\ c_p,\\ k,\\ D$ 随 $C,T$ 变化\n$h=25$，$h_m=8\\times10^{-7}$',
    C_INPUT)

# ============ 第二层：热-质耦合模型 ============
band_label(ax, 54.0, '径向一维轴对称热–质双向耦合模型', C_MODEL[1])
# 大容器（虚线框）
p = FancyBboxPatch((5, 29.5), 76, 22.0,
                   boxstyle="round,pad=0,rounding_size=1.0",
                   facecolor='#F6F9FD', edgecolor=C_MODEL[1],
                   linewidth=1.4, linestyle=(0, (4, 2)), zorder=1)
ax.add_patch(p)

box(ax, 8, 34.5, 30, 14.5, '热传导方程（温度场 $T$）',
    '$\\rho(C)c_p(C)\\dfrac{\\partial T}{\\partial t}'
    '=\\dfrac{1}{r}\\dfrac{\\partial}{\\partial r}'
    '\\left(r\\,k(C)\\dfrac{\\partial T}{\\partial r}\\right)$',
    C_MODEL, title_fs=9.5, sub_fs=9, title_dy=4.6)
box(ax, 48, 34.5, 30, 14.5, '水分扩散方程（水分场 $C$）',
    '$\\dfrac{\\partial C}{\\partial t}'
    '=\\dfrac{1}{r}\\dfrac{\\partial}{\\partial r}'
    '\\left(r\\,D(C,T)\\dfrac{\\partial C}{\\partial r}\\right)$',
    C_MODEL, title_fs=9.5, sub_fs=9, title_dy=4.6)
# 耦合双向箭头（间隙 38–48）
arrow(ax, 38.8, 44.5, 47.2, 44.5, color='#B2182B', lw=1.4, rad=-0.30)
arrow(ax, 47.2, 37.5, 38.8, 37.5, color='#1F4E9C', lw=1.4, rad=-0.30)
ax.text(43, 46.3, '$T\\to D$', ha='center', fontsize=7.5, color='#B2182B',
        fontweight='bold')
ax.text(43, 34.6, '$C\\to \\rho,c_p,k$', ha='center', fontsize=7.5,
        color='#1F4E9C', fontweight='bold')
# 耦合与边界说明（容器底部两行）
ax.text(43, 32.5,
        '耦合：$T$ 经 Arrhenius 因子 $\\mathrm{e}^{-3850/T}$ 影响 $D$；'
        '$C$ 经 $\\rho,c_p,k$ 影响温度场',
        ha='center', fontsize=7.4, color=SUB)
ax.text(43, 30.4,
        '对称边界 $\\partial_r\\,\\cdot\\,(0,t)=0$　｜　'
        'Robin 对流边界：热 $h(T_R-T_a)$，质 $h_m(C_R-C_a)$',
        ha='center', fontsize=7.4, color=SUB)
# 动边界模块（问题4，容器右侧）
box(ax, 84, 34.5, 13.5, 14.5, '动边界\n（问题 4）',
    '随体坐标\n$\\zeta=r/R(t)$\n$D/R^2$ 增大\n路径缩短',
    C_MODEL, title_fs=8.6, sub_fs=7.2, title_dy=4.0)

# ============ 第三层：数值求解框架 ============
band_label(ax, 25.0, '数值求解框架（四问题共用）', C_NUM[1])
box(ax, 6, 14.5, 26, 8.5, '守恒型有限体积离散',
    '控制体积分 + 界面通量\n逐层严格守恒',
    C_NUM)
box(ax, 37, 14.5, 26, 8.5, 'Crank–Nicolson 时间推进',
    '无条件稳定、二阶精度\n三对角阵 + Thomas 算法',
    C_NUM)
box(ax, 68, 14.5, 26, 8.5, '交错 Picard 迭代（两轮）',
    '先 $T$ 后 $C$ 解耦耦合\n残差 < 单步变化 1%',
    C_NUM)

# ============ 第四层：问题分支与结果输出 ============
band_label(ax, 10.0, '问题分支与结果输出', C_OUT[1])
box(ax, 3, 0.5, 21, 8, '问题 1', '预热平衡 0–1800 s\n表 1/2 · result1.xlsx', C_OUT)
box(ax, 27.5, 0.5, 21, 8, '问题 2', '全过程前 3 h\n表 3/4 · result2.xlsx', C_OUT)
box(ax, 52, 0.5, 21, 8, '问题 3', '判据 $\\max_r C\\leq0.15$\n$t^*=57.09$ h · 表 5', C_OUT)
box(ax, 76.5, 0.5, 21, 8, '问题 4', '含收缩动边界\n$t^*=50.77$ h · 表 6', C_OUT)

# ============ 第五层：模型检验（底部横条 + 右侧反馈） ============
band_label(ax, -4.0, '模型检验', C_CHECK[1])
box(ax, 6, -12.5, 88, 7, None, None, C_CHECK)
ax.text(50, -9.0,
        '空间 / 时间收敛性检验（加密网格，差值 $O(10^{-5})$）　｜　'
        '离散守恒恒等式　｜　参数灵敏度分析（表 7，最不利组合增幅 < 20%）',
        ha='center', va='center', fontsize=7.8, color=SUB, zorder=3)
# 反馈虚线箭头：检验 → 模型层（沿右缘）
arrow(ax, 94.5, -5.5, 94.5, 33.0, color=C_CHECK[1], lw=1.3, ls=(0, (4, 2)))
ax.text(96.0, 14, '反馈修正', rotation=90, fontsize=7.5, color=C_CHECK[1],
        ha='center', va='center')

# ============ 层间主箭头 ============
for x0 in (19, 50, 81):
    arrow(ax, x0, 57.8, x0, 52.0, lw=1.4)
for x0 in (19, 50, 81):
    arrow(ax, x0, 29.3, x0, 23.4, lw=1.4)
for x0 in (13.5, 38, 62.5, 87):
    arrow(ax, x0, 14.3, x0, 8.8, lw=1.4)

fig.savefig(os.path.join(BASE, 'fig13_模型总体架构图.png'),
            dpi=300, bbox_inches='tight', pad_inches=0.1, facecolor='white')
fig.savefig(os.path.join(BASE, 'fig13_模型总体架构图.pdf'),
            bbox_inches='tight', pad_inches=0.1, facecolor='white')
plt.close(fig)
print('saved fig13')
