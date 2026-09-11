# -*- coding: utf-8 -*-
"""2026 CUMCM A 题 问题3 求解器
整个烘干过程模型，经验公式统一采用附录3（物性随 C、T 变），
热质双向耦合，交错迭代 + Crank--Nicolson + Thomas。
问题3：Δr=0.25 mm, Δt=2.5 s，从初态直至 max C <= 0.15 -> 烘干时长 t*、表5
       + result3.xlsx（每 60 s、每 0.1 cm）。
烘房条件：0--14400 s 用附件1 插值，之后恒温段取附件1 末值 (50.165 °C, 0.04986)。
验证：V1 离散质量守恒（机器精度）、V2 时间收敛、V3 网格收敛（0.5 mm 对照）、V4 渐近。
"""
import sys, os
import numpy as np
import openpyxl

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

ROOT = os.path.dirname(os.path.abspath(__file__))
ANNEX1 = os.path.join(ROOT, 'OriginalMaterial', 'A题', '附件', '附件1.xlsx')

R, dr, N = 0.02, 2.5e-4, 81            # 半径 2 cm，主网格 0.25 mm，81 节点
r = np.arange(N) * dr
T0, C0 = 28.0, 2.55
h, hm = 25.0, 8e-7
TEND, CEND = 50.165, 0.04986          # 附件1 末值（恒温干燥段）
CTARGET = 0.15                        # 问题3 烘干判据

# 附录3 经验式（问题2、问题3 统一采用）
rho23 = lambda C: 650.0 + 128.0 * C
cp23 = lambda C: 1450.0 + 2736.0 * C / (C + 1.0)
k23 = lambda C: 0.21 + 0.38 * C / (C + 1.0)
D23 = lambda C, T: 2.4e-3 * np.exp(-0.45 / C) * np.exp(-3850.0 / T)   # T 用开尔文

d = np.array([row for row in
              openpyxl.load_workbook(ANNEX1, data_only=True)['Sheet1'].iter_rows(values_only=True)][1:],
             float)
T_air = lambda t: np.interp(t, d[:, 0], d[:, 1]) if t <= 14400 else TEND
C_air = lambda t: np.interp(t, d[:, 0], d[:, 2]) if t <= 14400 else CEND

def thomas(a, b, c, rhs):
    n = len(rhs); cp = np.zeros(n); dp = np.zeros(n)
    cp[0] = c[0] / b[0]; dp[0] = rhs[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]; cp[i] = c[i] / m; dp[i] = (rhs[i] - a[i] * dp[i - 1]) / m
    x = np.zeros(n); x[-1] = dp[-1]
    for i in range(n - 2, -1, -1): x[i] = dp[i] - cp[i] * x[i + 1]
    return x

def crank_step(u, K, Cp, beta, bc_prev, bc_next, dt):
    """Cp_i (u^{n+1}_i - u^n_i)/dt = Lu + s；K: 每节点扩散系数；Cp: 每节点容量系数。
    热方程: K=k(C), Cp=rho(C)c_p(C), beta=h/k_s；质方程: K=D, Cp=1, beta=h_m/D_s。"""
    Km = (K[:-1] + K[1:]) / 2
    lo = np.zeros(N); di = np.zeros(N); up = np.zeros(N); s = np.zeros(N)
    di[0] = -4 * Km[0] / dr ** 2; up[0] = 4 * Km[0] / dr ** 2
    for i in range(1, N - 1):
        lo[i] = Km[i - 1] * (i - 0.5) / (i * dr ** 2)
        di[i] = -(Km[i] * (i + 0.5) + Km[i - 1] * (i - 0.5)) / (i * dr ** 2)
        up[i] = Km[i] * (i + 0.5) / (i * dr ** 2)
    # 表面控制体 [r_{N-3/2}, R]：表面节点 r_{N-1}=(N-1)dr=R 位于真实表面，
    # Robin 条件直接取节点值：表面通量 -R·beta·K_s·(u_{N-1}-bc)
    # （热: -R·h(T_R-T_air)；质: -R·h_m(C_R-C_air)），与内部界面通量逐时间层严格守恒。
    vhat = ((N - 1) ** 2 - (N - 1.5) ** 2) * dr ** 2 / 2
    Bc = (N - 1) * dr * beta * K[-1]
    lo[N - 1] = Km[N - 2] * (N - 1.5) / vhat
    di[N - 1] = -(lo[N - 1] + Bc / vhat)
    up[N - 1] = 0.0
    def src(bc):
        s = np.zeros(N); s[N - 1] = Bc * bc / vhat; return s
    # Crank--Nicolson: Cp_i (u^{n+1}-u^n)/dt = 1/2(Lu^{n+1}+Lu^n) + 1/2(s^{n+1}+s^n)
    lhs_lo = -0.5 * lo
    lhs_di = Cp / dt - 0.5 * di
    lhs_up = -0.5 * up
    rhs = Cp / dt * u + 0.5 * (lo * np.roll(u, 1) + di * u + up * np.roll(u, -1)) \
        + 0.5 * (src(bc_next) + src(bc_prev))
    rhs[0] = Cp[0] / dt * u[0] + 0.5 * (di[0] * u[0] + up[0] * u[1]) \
        + 0.5 * (src(bc_next)[0] + src(bc_prev)[0])
    rhs[N - 1] = Cp[N - 1] / dt * u[N - 1] \
        + 0.5 * (lo[N - 1] * u[N - 2] + di[N - 1] * u[N - 1]) \
        + 0.5 * (src(bc_next)[N - 1] + src(bc_prev)[N - 1])
    return thomas(lhs_lo, lhs_di, lhs_up, rhs)

def step(T, C, t_prev, t_cur, dt):
    """一个时间步：交错耦合，两轮 Picard（每轮都从 u^n 出发，重做同一步）。"""
    Tn, Cn = T.copy(), C.copy()
    for _ in range(2):
        Kk = k23(C); Rc = rho23(C) * cp23(C)
        T = crank_step(Tn, Kk, Rc, h / Kk[-1], T_air(t_prev), T_air(t_cur), dt)
        Tk = T + 273.15
        Dc = D23(C, Tk)
        C = crank_step(Cn, Dc, np.ones(N), hm / Dc[-1], C_air(t_prev), C_air(t_cur), dt)
    return T, C

def run_drying(drg, Ng, dtt, sample6h=True):
    """全程求解直至 max C <= 0.15。返回 (tstar, tab5, Cfinal, CR_hist)。
    tab5: 每 6 h 一行 + 烘干结束时刻末行（0/0.5/1/1.5/2 cm）。"""
    global N, dr, r
    N, dr, r = Ng, drg, np.arange(Ng) * drg
    T, C = np.full(N, T0), np.full(N, C0)
    tstar, nstar = None, None
    ri5 = [0, Ng // 4, Ng // 2, 3 * Ng // 4, Ng - 1]
    tab5 = []
    CR_hist = [hm * (C0 - C_air(0.0))]
    nmax = int(300 * 3600 / dtt)       # 最多 300 h
    for n in range(1, nmax + 1):
        T, C = step(T, C, (n - 1) * dtt, n * dtt, dtt)
        CR_hist.append(hm * (C[-1] - C_air(n * dtt)))
        if sample6h and n % int(6 * 3600 / dtt) == 0:
            tab5.append(C[ri5].copy())
        if C.max() <= CTARGET:
            tstar, nstar = n * dtt, n
            break
    if tstar is None:
        raise RuntimeError('300 h 内未达烘干判据')
    if sample6h:
        tab5.append(C[ri5].copy())     # 末行：烘干结束时刻 t*
    return tstar, np.array(tab5), C.copy(), np.array(CR_hist)

# ================= 问题3：Δr=0.25 mm, Δt=2.5 s，直至 max C <= 0.15 =================
dt3 = 2.5
tstar, tab5, C3f, CR_hist = run_drying(2.5e-4, 81, dt3)
nstar = int(tstar / dt3)
print('问题3 烘干时长: t* = %d s = %.4f h (max C = %.6f)' % (tstar, tstar / 3600, C3f.max()))
tab5r = np.round(tab5, 4)
print('表5 水分浓度 (kg/kg)，行=6,12,...h 直至 t*，列=0/0.5/1/1.5/2 cm:')
print(tab5r)

# result3.xlsx：每 60 s、每 0.1 cm（0.25 mm 网格上每 4 节点一列）
ri10 = np.arange(0, N, 4)
T, C = np.full(N, T0), np.full(N, C0)
wb = openpyxl.Workbook()
ws = wb.active; ws.title = 'Sheet1'
ws.cell(1, 1, '时间\\到药材中心的距离')
for j, i in enumerate(ri10):
    ws.cell(1, 2 + j, round(r[i] * 100, 1))
k = 0
nmax = int(300 * 3600 / dt3)
for n in range(1, nmax + 1):
    T, C = step(T, C, (n - 1) * dt3, n * dt3, dt3)
    if n % 24 == 0:                    # 每 60 s 记一行
        k += 1
        ws.cell(1 + k, 1, n * dt3)
        for j in range(len(ri10)):
            ws.cell(1 + k, 2 + j, round(float(C[ri10[j]]), 4))
    if C.max() <= CTARGET:
        break
if n % 24:                             # 末行补烘干结束时刻
    k += 1
    ws.cell(1 + k, 1, n * dt3)
    for j in range(len(ri10)):
        ws.cell(1 + k, 2 + j, round(float(C[ri10[j]]), 4))
wb.save(os.path.join(ROOT, 'result3.xlsx'))
print('已写出 result3.xlsx（%d 行）' % (k + 1))

# ================= 验证 =================
# V1 离散质量守恒（全程：ΣvΔC vs 表面通量积分，v 为控制体体积、C_R 即表面节点值）
v = np.zeros(N); v[0] = dr ** 2 / 8
v[1:N - 1] = np.arange(1, N - 1) * dr ** 2
v[N - 1] = ((N - 1) ** 2 - (N - 1.5) ** 2) * dr ** 2 / 2
lhs = (C3f * v).sum() - C0 * v.sum()
rhs = -R * np.trapezoid(CR_hist, np.arange(0, len(CR_hist)) * dt3)
print('验证V1 离散守恒(0.25mm): ΣvΔC = %.6e，-R∫h_m(C_R-C_air)dt = %.6e，差 %.2e'
      % (lhs, rhs, lhs - rhs))

# V2 时间收敛：Δt=2.5 vs 1.25（前 6 h 表5 首行，主网格）
def row6h(dt):
    global N, dr, r
    N, dr, r = 81, 2.5e-4, np.arange(81) * 2.5e-4
    T, C = np.full(N, T0), np.full(N, C0)
    for n in range(1, int(21600 / dt) + 1):
        T, C = step(T, C, (n - 1) * dt, n * dt, dt)
    ri5 = [0, 20, 40, 60, 80]
    return C[ri5].copy()
r25 = row6h(2.5); r125 = row6h(1.25)
print('验证V2 时间收敛(6h): dt=2.5 vs 1.25 最大差 %.2e kg/kg' % np.max(np.abs(r25 - r125)))

# V3 网格收敛：Δr=0.5 mm 全程对照（t* 与 表5）
tstar_c, tab5_c, _, _ = run_drying(5e-4, 41, 2.5)
print('验证V3 网格收敛: 0.5mm t* = %.4f h vs 0.25mm t* = %.4f h (差 %.3f h)'
      % (tstar_c / 3600, tstar / 3600, (tstar_c - tstar) / 3600))
print('              表5 最大差 0.5mm vs 0.25mm: %.2e kg/kg'
      % np.max(np.abs(tab5_c - tab5)))

# V4 渐近：常数边界 (50.165, 0.04986)，长时程
N, dr, r = 21, 1e-3, np.arange(21) * 1e-3
T, C = np.full(N, T0), np.full(N, C0)
for n in range(1, int(3e5 / dt3) + 1):
    T, C = step(T, C, (n - 1) * dt3, n * dt3, dt3)
    if n * dt3 in (100000, 200000, 300000):
        print('验证V4 渐近: t=%ds 中心 T=%.4f (目标 %.4f)，表面 C=%.5f (目标 %.5f)'
              % (n * dt3, T[0], TEND, C[-1], CEND))
