# -*- coding: utf-8 -*-
"""q4 短时演化探针:打印前若干步的表面/中心值,检查表面行系数。"""
import numpy as np
import openpyxl
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.interpolate import PchipInterpolator

ROOT = os.path.dirname(os.path.abspath(__file__))
ANNEX1 = os.path.join(ROOT, 'OriginalMaterial', 'A题', '附件', '附件1.xlsx')
ANNEX2 = os.path.join(ROOT, 'OriginalMaterial', 'A题', '附件', '附件2.xlsx')

R0, T0, C0 = 0.02, 28.0, 2.55
h, hm = 25.0, 8e-7
TEND, CEND = 50.165, 0.04986
DZ, NZ = 1.0 / 40.0, 41
zeta = np.arange(NZ) * DZ

rho4 = lambda C: 760.0 + 90.0 * C
cp4 = lambda C: 1850.0 + 2150.0 * C / (C + 1.0)
k4 = lambda C: 0.12 + 0.20 * C / (C + 1.0)
D4 = lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 0.02)) * np.exp(-3850.0 / T)

d1 = np.array([row for row in
               openpyxl.load_workbook(ANNEX1, data_only=True)['Sheet1'].iter_rows(values_only=True)][1:],
              float)
T_air = lambda t: np.interp(t, d1[:, 0], d1[:, 1]) if t <= 14400 else TEND
C_air = lambda t: np.interp(t, d1[:, 0], d1[:, 2]) if t <= 14400 else CEND

d2 = np.array([row for row in
               openpyxl.load_workbook(ANNEX2, data_only=True)['Sheet1'].iter_rows(values_only=True)][1:],
              float)
T2END = d2[-1, 0]; REND = d2[-1, 1] / 100.0
_Rp = PchipInterpolator(d2[:, 0], d2[:, 1] / 100.0, extrapolate=False)
RofT = lambda t: _Rp(t) if t <= T2END else REND
_Rpd = _Rp.derivative()
RpofT = lambda t: float(_Rpd(t)) if t <= T2END else 0.0

def thomas(a, b, c, rhs):
    n = len(rhs); cp = np.zeros(n); dp = np.zeros(n)
    cp[0] = c[0] / b[0]; dp[0] = rhs[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]; cp[i] = c[i] / m; dp[i] = (rhs[i] - a[i] * dp[i - 1]) / m
    x = np.zeros(n); x[-1] = dp[-1]
    for i in range(n - 2, -1, -1): x[i] = dp[i] - cp[i] * x[i + 1]
    return x

def crank_step(u, K, Cp, R, beta, bc_prev, bc_next, dt):
    Km = (K[:-1] + K[1:]) / 2
    R2 = R * R
    lo = np.zeros(NZ); di = np.zeros(NZ); up = np.zeros(NZ); s = np.zeros(NZ)
    di[0] = -4 * Km[0] / (DZ ** 2 * R2); up[0] = 4 * Km[0] / (DZ ** 2 * R2)
    for i in range(1, NZ - 1):
        lo[i] = Km[i - 1] * (i - 0.5) / (i * DZ ** 2 * R2)
        di[i] = -(Km[i] * (i + 0.5) + Km[i - 1] * (i - 0.5)) / (i * DZ ** 2 * R2)
        up[i] = Km[i] * (i + 0.5) / (i * DZ ** 2 * R2)
    vhat = (1.0 - ((NZ - 1.5) * DZ) ** 2) / 2.0
    Bc = beta * K[-1]
    lo[NZ - 1] = Km[NZ - 2] * (NZ - 1.5) / (DZ * vhat * R2)
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

T, C = np.full(NZ, T0), np.full(NZ, C0)
print('初始: D4(C0,300K)=%.3e  R(0)=%.4f  C_air(0)=%.5f' % (D4(C0, 28+273.15), RofT(0), C_air(0)))
for n in range(1, 13):
    T, C = step(T, C, (n - 1) * 5.0, n * 5.0, 5.0)
    if n in (1, 2, 4, 8, 12):
        print('t=%3ds 中心T=%.4f 表面T=%.4f | 中心C=%.5f 表面C=%.5f'
              % (n * 5, T[0], T[-1], C[0], C[-1]))
