# -*- coding: utf-8 -*-
"""诊断：纯热传导（常数 k，跳变 BC）下当前表面行格式的收敛阶与误差分布。
对照 q1_solver.crank_step 的行结构。"""
import numpy as np
from scipy.special import j0, j1
from scipy.optimize import brentq

R, k, h, rho, cp = 0.02, 0.36, 25.0, 820.0, 2600.0
alpha = k / (rho * cp)
T0, Ta = 28.0, 41.513
Bi = h * R / k
f = lambda b: b * j1(b) - Bi * j0(b)
bs = [brentq(f, lo, hi) for lo, hi in
      zip(np.arange(0.05, 400, 0.5), np.arange(0.55, 400.5, 0.5)) if f(lo) * f(hi) < 0][:40]
bs = np.array(bs)
An = 2 * j1(bs) / (bs * (j0(bs) ** 2 + j1(bs) ** 2))

def T_exact(rr, t):
    return Ta + (T0 - Ta) * (An * np.exp(-alpha * bs ** 2 * t / R ** 2) * j0(bs * rr / R)).sum()

def thomas(a, b, c, rhs):
    n = len(rhs); cp = np.zeros(n); dp = np.zeros(n)
    cp[0] = c[0] / b[0]; dp[0] = rhs[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]; cp[i] = c[i] / m; dp[i] = (rhs[i] - a[i] * dp[i - 1]) / m
    x = np.zeros(n); x[-1] = dp[-1]
    for i in range(n - 2, -1, -1): x[i] = dp[i] - cp[i] * x[i + 1]
    return x

def crank_row(u, N, dr, dt, bc_prev, bc_next, variant='node'):
    """variant: 'node' = 当前格式（表面通量直接取节点值）;
    'ghost' = 鬼点梯度格式; 'cap' = 节点格式+表面控制体容量二次修正."""
    Km = alpha  # 常数 k
    lo = np.zeros(N); di = np.zeros(N); up = np.zeros(N)
    di[0] = -4 * Km / dr ** 2; up[0] = 4 * Km / dr ** 2
    for i in range(1, N - 1):
        lo[i] = Km * (i - 0.5) / (i * dr ** 2)
        di[i] = -(2 * Km) / dr ** 2
        up[i] = Km * (i + 0.5) / (i * dr ** 2)
    vhat = (N * N - (N - 1.5) ** 2) * dr ** 2 / 2
    Bc = N * dr * (h / k) * Km                      # = R·h/(rho cp)
    lo[N - 1] = Km * (N - 1.5) / vhat
    di[N - 1] = -(lo[N - 1] + Bc / vhat)
    up[N - 1] = 0.0
    if variant == 'ghost':
        kap = dr * h / k
        Bc = N * dr * (h / k) * Km * 2 / (1 + kap)  # 鬼点梯度: -2R·h(u-bc)/(1+kappa)
        di[N - 1] = -(lo[N - 1] + Bc / vhat)
    def src(bc):
        s = np.zeros(N); s[N - 1] = Bc * bc / vhat; return s
    lhs_lo = -0.5 * dt * lo
    lhs_di = 1 - 0.5 * dt * di
    lhs_up = -0.5 * dt * up
    rhs = u + 0.5 * dt * (lo * np.roll(u, 1) + di * u + up * np.roll(u, -1)) \
        + 0.5 * dt * (src(bc_next) + src(bc_prev))
    rhs[0] = u[0] + 0.5 * dt * (di[0] * u[0] + up[0] * u[1])
    rhs[N - 1] = u[N - 1] + 0.5 * dt * (lo[N - 1] * u[N - 2] + di[N - 1] * u[N - 1]) \
        + 0.5 * dt * (src(bc_next)[N - 1] + src(bc_prev)[N - 1])
    return thomas(lhs_lo, lhs_di, lhs_up, rhs)

def run(dr, dt, variant='node'):
    N = int(R / dr) + 1
    u = np.full(N, T0)
    for n in range(1, int(1800 / dt) + 1):
        u = crank_row(u, N, dr, dt, Ta, Ta, variant)
    return u

def run_dirichlet(dr, dt):
    """表面取解析精确值 u_R(t)=T_exact(R,t)，其余行与 'node' 相同——隔离 Robin 行。"""
    N = int(R / dr) + 1
    u = np.full(N, T0)
    Km = alpha
    lo = np.zeros(N); di = np.zeros(N); up = np.zeros(N)
    di[0] = -4 * Km / dr ** 2; up[0] = 4 * Km / dr ** 2
    for i in range(1, N - 1):
        lo[i] = Km * (i - 0.5) / (i * dr ** 2)
        di[i] = -2 * Km / dr ** 2
        up[i] = Km * (i + 0.5) / (i * dr ** 2)
    up_save = up[N - 2]                     # 行 N-2 对外界面通量系数
    up[N - 2] = 0.0
    lhs_lo = -0.5 * dt * lo[:-1]; lhs_di = 1 - 0.5 * dt * di[:-1]; lhs_up = -0.5 * dt * up[:-1]
    for n in range(1, int(1800 / dt) + 1):
        t = n * dt
        uR_prev, uR_next = T_exact(R, t - dt), T_exact(R, t)
        rhs = u[:-1] + 0.5 * dt * (lo[:-1] * np.roll(u[:-1], 1) + di[:-1] * u[:-1]
                                   + up[:-1] * np.roll(u[:-1], -1))
        rhs[0] = u[0] + 0.5 * dt * (di[0] * u[0] + up[0] * u[1])
        rhs[N - 2] = u[N - 2] + 0.5 * dt * (lo[N - 2] * u[N - 3] + di[N - 2] * u[N - 2]) \
            + 0.5 * dt * up_save * (uR_prev + uR_next)
        u[:-1] = thomas(lhs_lo, lhs_di, lhs_up, rhs)
        u[-1] = uR_next
    return u

for variant in ('node', 'ghost'):
    print('== variant:', variant)
    prev = None
    for dr, dt in ((1e-3, 1.0), (5e-4, 0.5), (2.5e-4, 0.25)):
        u = run(dr, dt, variant)
        e = [u[0] - T_exact(0, 1800), u[-1] - T_exact(R, 1800)]
        print('dr=%g dt=%g 中心err=%+.4f 表面err=%+.4f' % (dr, dt, e[0], e[1]))
        if prev is not None:
            print('   误差比(中心/表面): %.2f / %.2f' % ((prev[0] / e[0]) if e[0] else 0,
                  (prev[1] / e[1]) if e[1] else 0))
        prev = e

print('== variant: dirichlet(解析表面值)')
prev = None
for dr, dt in ((1e-3, 1.0), (5e-4, 0.5), (2.5e-4, 0.25)):
    u = run_dirichlet(dr, dt)
    e = [u[0] - T_exact(0, 1800), u[-1] - T_exact(R, 1800)]
    print('dr=%g dt=%g 中心err=%+.4f 表面err=%+.4f' % (dr, dt, e[0], e[1]))
    if prev is not None:
        print('   误差比(中心/表面): %.2f / %.2f' % ((prev[0] / e[0]) if e[0] else 0,
              (prev[1] / e[1]) if e[1] else 0))
    prev = e
