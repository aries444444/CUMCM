# -*- coding: utf-8 -*-
"""2026 CUMCM A 题 问题4 求解器（含收缩）
附录4 经验式（物性随 C、T 变），热质双向耦合，交错迭代 + Crank--Nicolson + Thomas。
收缩：附件2 半径 R(t)（Pchip 插值），拉格朗日坐标 zeta = r/R(t) 随物料变形，
      菲克/傅里叶定律相对固相骨架，zeta 系中无对流项：
      dC/dt = 1/(R^2 zeta) d(zeta D dC/dzeta)/dzeta，
      rho c_p dT/dt = 1/(R^2 zeta) d(zeta k dT/dzeta)/dzeta（d/dt 为固定 zeta 的随体导数）。
输出：表6（每 6 h x 0/0.5/1/1.5/2 cm + 表面，域外记 '-'）+ result4.xlsx（每 60 s，22 列）。
验证：V1 动边界质量守恒（含压缩项 2R'/R ∫C r dr）、V2 时间收敛、V3 网格收敛、V4 渐近。
"""
import sys, os
import numpy as np
import openpyxl
from scipy.linalg import solve_banded
from scipy.interpolate import PchipInterpolator

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

ROOT = os.path.dirname(os.path.abspath(__file__))
ANNEX1 = os.path.join(ROOT, 'OriginalMaterial', 'A题', '附件', '附件1.xlsx')
ANNEX2 = os.path.join(ROOT, 'OriginalMaterial', 'A题', '附件', '附件2.xlsx')

R0, T0, C0 = 0.02, 28.0, 2.55
h, hm = 25.0, 8e-7
TEND, CEND = 50.165, 0.04986        # 附件1 末值（恒温干燥段）
CTARGET = 0.15                      # 问题4 烘干判据
DZ, NZ = 1.0 / 40.0, 41             # zeta 网格 41 节点
zeta = np.arange(NZ) * DZ

# 附录4 经验式（C<0.02 时取 0.02，仅为数值保险：物理解 C >= C_air=0.04986）
rho4 = lambda C: 760.0 + 90.0 * C
cp4 = lambda C: 1850.0 + 2150.0 * C / (C + 1.0)
k4 = lambda C: 0.12 + 0.20 * C / (C + 1.0)
D4 = lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 0.02)) * np.exp(-3850.0 / T)  # T 开尔文

d1 = np.array([row for row in
               openpyxl.load_workbook(ANNEX1, data_only=True)['Sheet1'].iter_rows(values_only=True)][1:],
              float)
T_air = lambda t: np.interp(t, d1[:, 0], d1[:, 1]) if t <= 14400 else TEND
C_air = lambda t: np.interp(t, d1[:, 0], d1[:, 2]) if t <= 14400 else CEND

d2 = np.array([row for row in
               openpyxl.load_workbook(ANNEX2, data_only=True)['Sheet1'].iter_rows(values_only=True)][1:],
              float)
T2END = d2[-1, 0]; REND = d2[-1, 1] / 100.0    # 259200 s 后半径保持 1.198 cm（附件2 单位为 cm）
_Rp = PchipInterpolator(d2[:, 0], d2[:, 1] / 100.0, extrapolate=False)
RofT = lambda t: _Rp(t) if t <= T2END else REND
_Rpd = _Rp.derivative()
RpofT = lambda t: float(_Rpd(t)) if t <= T2END else 0.0

def thomas(a, b, c, rhs):
    """三对角求解:与追赶法数学等价,调用 LAPACK 带状求解(提速,结果不变)。"""
    n = len(rhs)
    ab = np.zeros((3, n))
    ab[0, 1:] = c[:-1]
    ab[1, :] = b
    ab[2, :-1] = a[1:]
    return solve_banded((1, 1), ab, rhs)

def crank_step(u, K, Cp, R, beta, bc_prev, bc_next, dt):
    """zeta 系 Crank--Nicolson 一步（物料随体坐标，无对流项；扩散系数含 1/R^2）。
    beta: 表面 Robin 系数（热 h*R/k_s；质 h_m*R/D_s）。"""
    Km = (K[:-1] + K[1:]) / 2
    R2 = R * R
    lo = np.zeros(NZ); di = np.zeros(NZ); up = np.zeros(NZ); s = np.zeros(NZ)
    di[0] = -4 * Km[0] / (DZ ** 2 * R2); up[0] = 4 * Km[0] / (DZ ** 2 * R2)
    for i in range(1, NZ - 1):
        lo[i] = Km[i - 1] * (i - 0.5) / (i * DZ ** 2 * R2)
        di[i] = -(Km[i] * (i + 0.5) + Km[i - 1] * (i - 0.5)) / (i * DZ ** 2 * R2)
        up[i] = Km[i] * (i + 0.5) / (i * DZ ** 2 * R2)
    # 表面控制体 [zeta_{NZ-3/2}, 1]：表面节点位于真实表面 zeta=1（ζ 网格），
    # Robin 通量直接取 -R·beta·K_s·(u_{NZ-1}-bc)，逐时间层严格守恒。
    # 界面通量 = zeta_{NZ-3/2}·D_m·(u_{NZ-2}-u)/DZ = (NZ-1.5)·D_m·(u_{NZ-2}-u)（DZ 抵消）。
    vhat = (1.0 - ((NZ - 1.5) * DZ) ** 2) / 2.0
    Bc = beta * K[-1]                              # beta 已含 R：热 h·R/k_s，质 h_m·R/D_s
    lo[NZ - 1] = Km[NZ - 2] * (NZ - 1.5) / (vhat * R2)
    di[NZ - 1] = -(lo[NZ - 1] + Bc / (vhat * R2))
    up[NZ - 1] = 0.0
    def src(bc):
        s = np.zeros(NZ); s[NZ - 1] = Bc * bc / (vhat * R2); return s
    lhs_lo = -0.5 * lo
    lhs_di = Cp / dt - 0.5 * di
    lhs_up = -0.5 * up
    rhs = Cp / dt * u + 0.5 * (lo * np.roll(u, 1) + di * u + up * np.roll(u, -1)) \
        + 0.5 * (src(bc_next) + src(bc_prev))
    rhs[0] = Cp[0] / dt * u[0] + 0.5 * (di[0] * u[0] + up[0] * u[1]) \
        + 0.5 * (src(bc_next)[0] + src(bc_prev)[0])
    rhs[NZ - 1] = Cp[NZ - 1] / dt * u[NZ - 1] \
        + 0.5 * (lo[NZ - 1] * u[NZ - 2] + di[NZ - 1] * u[NZ - 1]) \
        + 0.5 * (src(bc_next)[NZ - 1] + src(bc_prev)[NZ - 1])
    return thomas(lhs_lo, lhs_di, lhs_up, rhs)

def step(T, C, t_prev, t_cur, dt):
    """一个时间步：R、R' 取半步中点，交错耦合，两轮 Picard（每轮从 u^n 出发）。"""
    tm = 0.5 * (t_prev + t_cur)
    R = RofT(tm)
    Tn, Cn = T.copy(), C.copy()
    for _ in range(2):
        Kk = k4(C); Rc = rho4(C) * cp4(C)
        T = crank_step(Tn, Kk, Rc, R, h * R / Kk[-1], T_air(t_prev), T_air(t_cur), dt)
        Tk = T + 273.15
        Dc = D4(C, Tk)
        C = crank_step(Cn, Dc, np.ones(NZ), R, hm * R / Dc[-1], C_air(t_prev), C_air(t_cur), dt)
    return T, C

def water_int(C, R):
    """W = Σ v_i C_i（随体域离散控制体体积积分，v_i 与 crank_step 一致）。"""
    v = np.zeros(NZ); v[0] = DZ ** 2 / 8
    v[1:NZ - 1] = zeta[1:NZ - 1] * DZ
    v[NZ - 1] = (1.0 - ((NZ - 1.5) * DZ) ** 2) / 2.0
    return R * R * (C * v).sum()

if __name__ == "__main__":
    # ================= 问题4：dt=5 s，直到 max C <= 0.15 =================
    DT = 5.0
    T, C = np.full(NZ, T0), np.full(NZ, C0)
    dist = np.array([0.0, 0.005, 0.01, 0.015, 0.02])     # 表6 固定距离列（m）
    xlsx_d = np.arange(0.0, 0.02001, 0.001)             # result4 0:0.1:2 cm（m）
    tab6, rows4 = [], []
    Whist, Cshist, Rhist, Rphist = [], [], [], []
    tstar, nstar = None, None
    nmax = int(300 * 3600 / DT)
    for n in range(1, nmax + 1):
        tn = n * DT
        T, C = step(T, C, tn - DT, tn, DT)
        R = RofT(tn); Rp = RpofT(tn)
        W = water_int(C, R)
        Whist.append(W); Cshist.append(C[-1]); Rhist.append(R); Rphist.append(Rp)
        if n % 4320 == 0:                       # 每 6 h 记表6 一行
            row = []
            for dj in dist:
                row.append(float(np.interp(dj / R, zeta, C)) if dj <= R else np.nan)
            row.append(C[-1])                   # 表面列
            tab6.append(row)
        if n % 12 == 0:                         # 每 60 s 记 result4 一行
            row = []
            for dj in xlsx_d:
                row.append(float(np.interp(dj / R, zeta, C)) if dj <= R else None)
            row.append(C[-1])
            rows4.append((tn, row))
        if C.max() <= CTARGET:
            tstar, nstar = tn, n
            break
    if tstar is None:
        raise RuntimeError('300 h 内未达烘干判据')
    R = RofT(tstar)
    row = [float(np.interp(dj / R, zeta, C)) if dj <= R else np.nan for dj in dist] + [C[-1]]
    tab6.append(row)
    if nstar % 12:
        row = [float(np.interp(dj / R, zeta, C)) if dj <= R else None for dj in xlsx_d] + [C[-1]]
        rows4.append((tstar, row))
    print('问题4 烘干时长: t* = %d s = %.4f h (max C = %.6f, 此时 R = %.4f cm)'
          % (tstar, tstar / 3600, C.max(), R * 100))
    np.set_printoptions(precision=4, suppress=True)
    print('表6 水分浓度 (kg/kg)，行=6,12,...h 直至 t*，列=0/0.5/1/1.5/2 cm + 表面:')
    for k, row in enumerate(tab6):
        rr = ' '.join('  --' if np.isnan(v) else '%6.4f' % v for v in row)
        trow = (k + 1) * 6.0 if k + 1 < len(tab6) else tstar / 3600.0
        print('%7.1fh |%s' % (trow, rr))

    wb = openpyxl.Workbook()
    ws = wb.active; ws.title = 'Sheet1'
    ws.cell(1, 1, '时间\\到药材中心的距离')
    for j in range(len(xlsx_d)):
        ws.cell(1, 2 + j, round(xlsx_d[j] * 100, 1))
    ws.cell(1, 2 + len(xlsx_d), '药材表面')
    for k, (tt, row) in enumerate(rows4):
        ws.cell(2 + k, 1, tt)
        for j, v in enumerate(row):
            if v is not None:
                ws.cell(2 + k, 2 + j, round(v, 4))
    wb.save(os.path.join(ROOT, 'result4.xlsx'))
    print('已写出 result4.xlsx（%d 行）' % len(rows4))

    # ================= 验证 =================
    # V1 动边界质量守恒: d/dt ∫C r dr = 2R'/R·∫C r dr - R·h_m(C_s-C_air)
    W0 = water_int(np.full(NZ, C0), R0)
    lhs = Whist[-1] - W0
    tt = np.arange(1, nstar + 1) * DT
    f1 = [2 * rp / r * w for rp, r, w in zip(Rphist, Rhist, Whist)]
    f2 = [-r * hm * (cs - C_air(t)) for r, cs, t in zip(Rhist, Cshist, tt)]
    rhs = np.trapezoid([a + b for a, b in zip(f1, f2)], tt)
    print('验证V1 离散守恒: ΣvΔC = %.6e，∫[2R\'/R·W - R h_m(C_R-C_air)]dt = %.6e，差 %.2e'
          % (lhs, rhs, lhs - rhs))

    # V2 时间收敛: 前 6 h，dt=5 vs 2.5（同网格）
    def run6h(dt):
        T, C = np.full(NZ, T0), np.full(NZ, C0)
        for n in range(1, int(21600 / dt) + 1):
            T, C = step(T, C, (n - 1) * dt, n * dt, dt)
        R = RofT(21600)
        row = [float(np.interp(dj / R, zeta, C)) if dj <= R else np.nan for dj in dist] + [C[-1]]
        return np.array(row)
    r5 = run6h(5.0); r25 = run6h(2.5)
    print('验证V2 时间收敛(6h 行): dt=5 vs 2.5 最大差 %.2e' % np.nanmax(np.abs(r5 - r25)))

    # V3 网格收敛: dZ=1/20（NZ=21）vs 1/40，前 6 h
    _NZ_save, _DZ_save, _zeta_save = NZ, DZ, zeta
    def run6h_coarse():
        global NZ, DZ, zeta
        NZ, DZ = 21, 1.0 / 20.0
        zeta = np.arange(NZ) * DZ
        T, C = np.full(NZ, T0), np.full(NZ, C0)
        for n in range(1, int(21600 / 5.0) + 1):
            T, C = step(T, C, (n - 1) * 5.0, n * 5.0, 5.0)
        R = RofT(21600)
        row = [float(np.interp(dj / R, zeta, C)) if dj <= R else np.nan for dj in dist] + [C[-1]]
        NZ, DZ = _NZ_save, _DZ_save
        zeta = _zeta_save
        return np.array(row)
    rc = run6h_coarse()
    print('验证V3 网格收敛(6h 行): dZ=1/20 vs 1/40 最大差 %.2e' % np.nanmax(np.abs(rc - r5)))

    # V4 渐近: 常数边界，长时程（R 仍按附件2 收缩并保持 1.198 cm）
    T, C = np.full(NZ, T0), np.full(NZ, C0)
    for n in range(1, int(2.5e5 / DT) + 1):
        T, C = step(T, C, (n - 1) * DT, n * DT, DT)
        if n * DT in (50000, 150000, 250000):
            print('验证V4 渐近: t=%ds 中心 T=%.4f (目标 %.4f)，表面 C=%.5f (目标 %.5f)，R=%.4f cm'
                  % (n * DT, T[0], TEND, C[-1], CEND, RofT(n * DT) * 100))
