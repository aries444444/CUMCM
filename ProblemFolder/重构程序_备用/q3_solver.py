# -*- coding: utf-8 -*-
"""2026 CUMCM A 题 问题3 求解器
整个烘干过程模型，经验公式统一采用附录3（物性随 C、T 变），
热质双向耦合，交错 Picard + Crank--Nicolson；公共数值内核见 common.py。
问题3：Δr=0.0625 mm (N=321), Δt=2.5 s，从初态直至 max C <= 0.15 -> 烘干时长 t*、表5
       + result3.xlsx（每 60 s、每 0.1 cm）。
  步长依据：0.125 mm 生产离散的表 5 对齐行 Richardson 估计（按 t* 链实测阶 p=1.22）约 7.8e-5，
  高于四位小数对应的 0.5e-4 阈值；0.0625 mm 的估计（按 p=1.33）约 2.7e-5，达标，故生产离散取 0.0625 mm。
  收敛阶低于整阶二阶的原因：表 5 最大差位于 42 h 行、1.5 cm 列（时间-浓度剖面最陡处），且末行
  由 t* 时间插值取值，行收敛阶 1.1--1.3（V3 逐格定位与 t* 链给出），按实测阶做 Richardson 估计。
烘房条件：0--14400 s 用附件1 插值，之后恒温段取附件1 末值 (50.165 °C, 0.04986)。
事件时间：g(t)=max C(r,t)-0.15 在相邻时间层间线性插值求根（见 common.integrate_until_dry）。
验证：V1 离散质量守恒（机器精度，末行插值后为 1e-11 量级）、V2 时间收敛、
      V3 空间收敛（t* 五点链 + 表5 逐格差定位）、V4 渐近。
"""
import os
import numpy as np
import openpyxl

from common import (ROOT, crank_radial, integrate_until_dry, kernel,
                    load_annex, picard_step, setup_console, thomas_band)

setup_console()

R, dr, N = 0.02, 6.25e-5, 321          # 半径 2 cm，主离散 0.0625 mm，321 节点
r = np.arange(N) * dr
T0, C0 = 28.0, 2.55
h, hm = 25.0, 8e-7
TEND, CEND = 50.165, 0.04986          # 附件1 末值（恒温干燥段）

# 附录3 经验式（问题2、问题3 统一采用）
rho23 = lambda C: 650.0 + 128.0 * C
cp23 = lambda C: 1450.0 + 2736.0 * C / (C + 1.0)
k23 = lambda C: 0.21 + 0.38 * C / (C + 1.0)
D23 = lambda C, T: 2.4e-3 * np.exp(-0.45 / C) * np.exp(-3850.0 / T)   # T 用开尔文

d = load_annex('附件1.xlsx')
T_air = lambda t: np.interp(t, d[:, 0], d[:, 1]) if t <= 14400 else TEND
C_air = lambda t: np.interp(t, d[:, 0], d[:, 2]) if t <= 14400 else CEND

def crank_step(u, K, Cp, beta, bc_prev, bc_next, dt):
    """守恒型有限体积 + Crank--Nicolson 一步（固定半径；见 common.crank_radial）。"""
    return crank_radial(u, K, Cp, beta, bc_prev, bc_next, dt, N, dr, thomas_band)

def step(T, C, t_prev, t_cur, dt):
    """一个时间步：固定半径 + 附录3 物性的交错 Picard 耦合步。"""
    crk = kernel('radial', dt, N, dr, thomas_band)
    return picard_step(crk, T, C, dt, t_prev, t_cur,
                       rho23, cp23, k23, D23, T_air, C_air, h, hm)

def run_drying(drg, Ng, dtt, sample6h=True, rows_every=None, ri10=None):
    """全程求解直至 max C <= 0.15，t* 由判据事件求根给出。
    返回 (tstar, tab5, Cstar, CR_hist, iters, rels, rows60)。
    tab5: 每 6 h 一行 + 烘干结束时刻末行（0/0.5/1/1.5/2 cm，末行为 t* 处插值状态）。
    rows60: rows_every 秒记一行 C[ri10]（含 t* 处插值末行），rows_every=None 时不记。"""
    global N, dr, r
    N, dr, r = Ng, drg, np.arange(Ng) * drg
    T, C = np.full(N, T0), np.full(N, C0)
    ri5 = [0, Ng // 4, Ng // 2, 3 * Ng // 4, Ng - 1]
    tab5, rows60 = [], []
    CR_hist = [hm * (C0 - C_air(0.0))]        # F_n = h_m(C^n_R - C_air(t_n))，C_R 即表面节点
    iters, rels = [], []

    def one_step(T, C, t_prev, t_cur):
        T, C, it, rel = step(T, C, t_prev, t_cur, dtt)
        iters.append(it); rels.append(rel)
        return T, C

    def record(n, t, T, C):
        CR_hist.append(hm * (C[-1] - C_air(t)))
        if sample6h and n % int(6 * 3600 / dtt) == 0:
            tab5.append(C[ri5].copy())
        if rows_every is not None and n % int(rows_every / dtt) == 0:
            rows60.append((t, C[ri10].copy()))

    tstar, _, Cstar, _ = integrate_until_dry(one_step, T, C, dtt, record=record)
    if sample6h:
        tab5.append(Cstar[ri5].copy())     # 末行：烘干结束时刻 t*
    if rows_every is not None:
        rows60.append((tstar, Cstar[ri10].copy()))
    return tstar, np.array(tab5), Cstar, np.array(CR_hist), np.array(iters), np.array(rels), rows60

if __name__ == "__main__":
    # ================= 问题3：Δr=0.0625 mm, Δt=2.5 s，直至 max C <= 0.15 =================
    dt3 = 2.5
    ri10 = np.arange(0, N, int(round(0.001 / dr)))      # 0:0.1:2 cm，0.0625 mm 下每 16 节点
    tstar, tab5, C3f, CR_hist, iters, rels, rows60 = run_drying(dr, N, dt3, rows_every=60, ri10=ri10)
    print('问题3 烘干时长: t* = %.3f s = %.4f h (max C = %.6f)'
          % (tstar, tstar / 3600, C3f.max()))
    print('Picard 统计: 平均 %.2f 轮，最大 %d 轮，退出残差 中位 %.2e / 最大 %.2e'
          % (iters.mean(), iters.max(), np.median(rels), rels.max()))
    tab5r = np.round(tab5, 4)
    print('表5 水分浓度 (kg/kg)，行=6,12,...h 直至 t*，列=0/0.5/1/1.5/2 cm:')
    print(tab5r)

    # result3.xlsx：每 60 s、每 0.1 cm（与上同一次求解），末行补 t* 处插值状态
    wb = openpyxl.Workbook()
    ws = wb.active; ws.title = 'Sheet1'
    ws.cell(1, 1, '时间\\到药材中心的距离')
    for j, i in enumerate(ri10):
        ws.cell(1, 2 + j, round(r[i] * 100, 1))
    for k, (tt, row) in enumerate(rows60):
        ws.cell(2 + k, 1, round(tt, 1))
        for j in range(len(ri10)):
            ws.cell(2 + k, 2 + j, round(float(row[j]), 4))
    wb.save(os.path.join(ROOT, 'result3.xlsx'))
    print('已写出 result3.xlsx（%d 行）' % len(rows60))

    # ================= 验证 =================
    # V1 离散质量守恒（全程：ΣvΔC vs 表面通量积分；末行为 t* 处插值状态，
    # 插值使闭合差由机器精度升至 1e-11 量级，相对总失水量仍可忽略）
    v = np.zeros(N); v[0] = dr ** 2 / 8
    v[1:N - 1] = np.arange(1, N - 1) * dr ** 2
    v[N - 1] = ((N - 1) ** 2 - (N - 1.5) ** 2) * dr ** 2 / 2
    lhs = (C3f * v).sum() - C0 * v.sum()
    rhs = -R * np.trapezoid(CR_hist, np.arange(0, len(CR_hist)) * dt3)
    print('验证V1 离散守恒(0.0625mm): ΣvΔC = %.6e，-R∫h_m(C_R-C_air)dt = %.6e，差 %.2e'
          % (lhs, rhs, lhs - rhs))

    # V2 时间收敛：Δt=2.5 vs 1.25（前 6 h 表5 首行，主离散）
    def row6h(dt):
        global N, dr, r
        N, dr, r = 321, 6.25e-5, np.arange(321) * 6.25e-5
        T, C = np.full(N, T0), np.full(N, C0)
        for n in range(1, int(21600 / dt) + 1):
            T, C, it, rel = step(T, C, (n - 1) * dt, n * dt, dt)
        ri5 = [0, 80, 160, 240, 320]
        return C[ri5].copy()
    r25 = row6h(2.5); r125 = row6h(1.25)
    print('验证V2 时间收敛(6h): dt=2.5 vs 1.25 最大差 %.2e kg/kg' % np.max(np.abs(r25 - r125)))

    # V3 空间收敛：t* 五点链（0.5/0.25/0.125/0.0625/0.03125 mm，同 Δt=2.5 s）
    # + 表5 对齐行逐格差（定位最大差所在行/列）
    ts = []; tabs = []
    for drg, Ng in ((5e-4, 41), (2.5e-4, 81), (1.25e-4, 161), (6.25e-5, 321), (3.125e-5, 641)):
        t_, tab_, _, _, it_, rl_, _ = run_drying(drg, Ng, dt3)
        ts.append(t_); tabs.append(tab_)
        print('V3 链 Δr=%.5f mm: t* = %.4f h (%d s), 平均 %.2f 轮'
              % (drg * 1e3, t_ / 3600, round(t_), it_.mean()))
    th = [t / 3600.0 for t in ts]
    dd1, dd2, dd3, dd4 = th[1] - th[0], th[2] - th[1], th[3] - th[2], th[4] - th[3]
    p12, p23, p34 = np.log2(dd1 / dd2), np.log2(dd2 / dd3), np.log2(dd3 / dd4)
    ext = th[4] + dd4 / (2 ** p34 - 1.0)
    print('V3 实测收敛阶: p12=%.3f (0.5->0.25->0.125), p23=%.3f (0.25->0.125->0.0625),'
          ' p34=%.3f (0.125->0.0625->0.03125)' % (p12, p23, p34))
    print('V3 按实测阶外推: t* = %.4f h；最细两档差 = %.4f h；外推-最细 = %.4f h'
          % (ext, dd4, ext - th[4]))
    for lab, a, b in (('0.125 vs 0.0625 mm', tabs[2][:-1], tabs[3][:-1]),
                      ('0.0625 vs 0.03125 mm', tabs[3][:-1], tabs[4][:-1])):
        diff = np.abs(a - b)
        i, j = np.unravel_index(np.argmax(diff), diff.shape)
        print('V3 表5 逐格差 %s: 最大 %.2e 位于 第 %d 行 (%.0f h) 第 %d 列 (%.1f cm)，'
              '整表最大差 %.2e'
              % (lab, diff.max(), i + 1, (i + 1) * 6.0, j, j * 0.5, diff.max()))
    print('V3 表5 末行差(0.0625 vs 0.03125): %.2e kg/kg（t* 差 %.4f h）'
          % (np.max(np.abs(tabs[3][-1] - tabs[4][-1])), th[4] - th[3]))

    # V4 渐近：常数边界 (50.165, 0.04986)，长时程
    N, dr, r = 21, 1e-3, np.arange(21) * 1e-3
    T, C = np.full(N, T0), np.full(N, C0)
    for n in range(1, int(3e5 / dt3) + 1):
        T, C, it, rel = step(T, C, (n - 1) * dt3, n * dt3, dt3)
        if n * dt3 in (100000, 200000, 300000):
            print('验证V4 渐近: t=%ds 中心 T=%.4f (目标 %.4f)，表面 C=%.5f (目标 %.5f)'
                  % (n * dt3, T[0], TEND, C[-1], CEND))
