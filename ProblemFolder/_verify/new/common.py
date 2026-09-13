# -*- coding: utf-8 -*-
"""2026 CUMCM A 题 公共数值内核（四个求解器共用，避免重复实现）

控制方程（u 为温度 T 或水分浓度 C）：
    rho*cp*dT/dt = (1/r) d(r*k*dT/dr)/dr ,    dC/dt = (1/r) d(r*D*dC/dr)/dr
表面 r = R 取第三类（Robin）边界：
    -k*dT/dr|_R = h*(T_R - T_air) ,   -D*dC/dr|_R = h_m*(C_R - C_air)
中心 r = 0 由轴对称性 + 洛必达法则化为 4*(u_1 - u_0)/dr^2。
问题 4 含收缩时改用随体坐标 zeta = r/R(t)，方程形式不变且无对流项。

集中实现（各求解器只保留自己的物性经验式、边界数据与输出格式）：
    setup_console    控制台 UTF-8 输出
    annex / load_annex  赛题附件定位与读入
    thomas           纯 Python 追赶法（与教材逐式对照）
    thomas_band      LAPACK 带状求解（与追赶法数学等价，节点多时提速）
    crank_radial     固定半径：守恒型有限体积 + Crank--Nicolson 单步
    crank_zeta       随体坐标（动边界）：同一离散形式的单步
    picard_step      交错 Picard 耦合步（热--质双向弱耦合）
    integrate_until_dry  烘干时长：判据事件求根

两套离散内核共用同一 Robin 通量闭合方式：表面节点落在真实表面，表面通量与
内部界面通量在每一时间层逐项对消，故离散守恒恒等式严格成立（见 verify_q*.py）。
"""
import os
import sys
from functools import partial

import numpy as np
import openpyxl
from scipy.linalg import solve_banded

ROOT = os.path.dirname(os.path.abspath(__file__))
# 交错 Picard 迭代：T、C 无穷范数相对变化小于 1e-8 或满 20 轮停止。
# 依据：时间步误差量级实测为 O(1e-5)（Δt 减半差值），迭代残差须再低 2--3 个量级。
PICARD_TOL, PICARD_MAXIT = 1e-8, 20
CTARGET = 0.15          # 烘干判据：各处水分浓度 <= 0.15 kg/kg
TMAX = 300 * 3600       # 最长积分时长（s），超时视为未达判据


def setup_console():
    """控制台按 UTF-8 输出（Windows 终端下中文与符号不乱码）。"""
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass


def annex(name):
    """定位赛题附件数据：优先脚本同目录（支撑材料布局），其次仓库目录结构。"""
    for p in (os.path.join(ROOT, name),
              os.path.join(ROOT, 'OriginalMaterial', 'A题', '附件', name)):
        if os.path.exists(p):
            return p
    raise FileNotFoundError('未找到 %s：请将赛题附件文件放在脚本同目录下' % name)


def load_annex(name):
    """读入附件首个工作表，跳过表头，返回 float 数组。"""
    ws = openpyxl.load_workbook(annex(name), data_only=True)['Sheet1']
    return np.array([row for row in ws.iter_rows(values_only=True)][1:], float)


def thomas(a, b, c, rhs):
    """三对角方程组（下/主/上对角线 a,b,c）：纯 Python 追赶法。"""
    n = len(rhs); cp = np.zeros(n); dp = np.zeros(n)
    cp[0] = c[0] / b[0]; dp[0] = rhs[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]; cp[i] = c[i] / m; dp[i] = (rhs[i] - a[i] * dp[i - 1]) / m
    x = np.zeros(n); x[-1] = dp[-1]
    for i in range(n - 2, -1, -1): x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def thomas_band(a, b, c, rhs):
    """三对角求解：与追赶法数学等价，调用 LAPACK 带状求解（节点多时提速，结果不变）。"""
    n = len(rhs)
    ab = np.zeros((3, n))
    ab[0, 1:] = c[:-1]
    ab[1, :] = b
    ab[2, :-1] = a[1:]
    return solve_banded((1, 1), ab, rhs)


def crank_radial(u, K, Cp, beta, bc_prev, bc_next, dt, N, dr, solver=thomas):
    """固定半径：Cp_i (u^{n+1}_i - u^n_i)/dt = Lu + s 的 Crank--Nicolson 单步。

    K：逐节点扩散系数（热 alpha 或 k(C)，质 D）；Cp：逐节点容量系数（热 rho*cp，质 1）；
    beta：表面 Robin 系数 h/k_s 或 h_m/D_s（不含半径，表面项按 R=(N-1)dr 折算）。
    """
    Km = (K[:-1] + K[1:]) / 2
    lo = np.zeros(N); di = np.zeros(N); up = np.zeros(N)
    di[0] = -4 * Km[0] / dr ** 2; up[0] = 4 * Km[0] / dr ** 2
    for i in range(1, N - 1):
        lo[i] = Km[i - 1] * (i - 0.5) / (i * dr ** 2)
        di[i] = -(Km[i] * (i + 0.5) + Km[i - 1] * (i - 0.5)) / (i * dr ** 2)
        up[i] = Km[i] * (i + 0.5) / (i * dr ** 2)
    # 表面控制体 [r_{N-3/2}, R]：表面节点 r_{N-1}=(N-1)dr=R 位于真实表面，
    # Robin 条件直接取节点值：表面通量 -R*beta*K_s*(u_{N-1}-bc)
    # （热: -R*h(T_R-T_air)；质: -R*h_m(C_R-C_air)），与内部界面通量逐时间层严格守恒。
    vhat = ((N - 1) ** 2 - (N - 1.5) ** 2) * dr ** 2 / 2   # 控制体体积 ∫ r dr
    Bc = (N - 1) * dr * beta * K[-1]                      # = R*beta*K_s
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
    return solver(lhs_lo, lhs_di, lhs_up, rhs)


def crank_zeta(u, K, Cp, beta, bc_prev, bc_next, dt, NZ, DZ, R, solver=thomas_band):
    """随体坐标 zeta = r/R(t) 下的同格式单步（无对流项，扩散项含 1/R^2）。

    beta：表面 Robin 系数，已含半径因子（热 h*R/k_s，质 h_m*R/D_s）。
    """
    Km = (K[:-1] + K[1:]) / 2
    R2 = R * R
    lo = np.zeros(NZ); di = np.zeros(NZ); up = np.zeros(NZ)
    di[0] = -4 * Km[0] / (DZ ** 2 * R2); up[0] = 4 * Km[0] / (DZ ** 2 * R2)
    for i in range(1, NZ - 1):
        lo[i] = Km[i - 1] * (i - 0.5) / (i * DZ ** 2 * R2)
        di[i] = -(Km[i] * (i + 0.5) + Km[i - 1] * (i - 0.5)) / (i * DZ ** 2 * R2)
        up[i] = Km[i] * (i + 0.5) / (i * DZ ** 2 * R2)
    # 表面控制体 [zeta_{NZ-3/2}, 1]：表面节点位于真实表面 zeta=1，
    # Robin 通量直接取 -R*beta*K_s*(u_{NZ-1}-bc)，逐时间层严格守恒。
    # 界面通量 = zeta_{NZ-3/2}*D_m*(u_{NZ-2}-u)/DZ = (NZ-1.5)*D_m*(u_{NZ-2}-u)（DZ 抵消）。
    vhat = (1.0 - ((NZ - 1.5) * DZ) ** 2) / 2.0
    Bc = beta * K[-1]                          # beta 已含 R：热 h*R/k_s，质 h_m*R/D_s
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
    return solver(lhs_lo, lhs_di, lhs_up, rhs)


def kernel(kind, dt, n, h, solver=None):
    """按几何绑定离散内核，返回可调用的 crank(u, K, Cp, beta, bc_prev, bc_next)。

    kind='radial'：固定半径，步长 h（m），beta 不含半径因子；
    kind='zeta'  ：随体坐标，步长 h=1/(n-1)，beta 已含当前半径 R。
    """
    if kind == 'radial':
        return partial(crank_radial, dt=dt, N=n, dr=h, solver=solver or thomas)
    return partial(crank_zeta, dt=dt, NZ=n, DZ=h, solver=solver or thomas_band)


def picard_step(crank, T, C, dt, t_prev, t_cur, rhoF, cpF, kF, DF,
                T_air, C_air, heff, hmeff, tol=PICARD_TOL, maxit=PICARD_MAXIT):
    """交错 Picard 耦合步：同一时间步内先由当前水分场求温度场，再用新温度场求水分场。

    crank：已绑定几何与 dt 的离散内核（见 kernel()）；
    rhoF/cpF/kF/DF：物性经验式；heff/hmeff：表面 Robin 系数（热 h、质 h_m，
    已含动边界的 R 因子）；温度对 D 的影响经 Arrhenius 因子 e^{-3850/T}（T 取开尔文）。
    每轮均从 u^n 出发；T、C 的相对变化（无穷范数，分母取 max(||u||inf,1)）小于 tol
    或满 maxit 轮停止。返回 (T, C, 轮数, 退出残差)。
    """
    Tn, Cn = T.copy(), C.copy()
    rel = 0.0
    for it in range(1, maxit + 1):
        Kk = kF(C); Rc = rhoF(C) * cpF(C)
        T = crank(Tn, Kk, Rc, heff / Kk[-1], T_air(t_prev), T_air(t_cur))
        Dc = DF(C, T + 273.15)
        C = crank(Cn, Dc, np.ones(len(C)), hmeff / Dc[-1], C_air(t_prev), C_air(t_cur))
        if it > 1:
            rel = max(np.max(np.abs(T - To)) / max(np.max(np.abs(T)), 1.0),
                      np.max(np.abs(C - Co)) / max(np.max(np.abs(C)), 1.0))
            if rel < tol:
                return T, C, it, rel
        To, Co = T.copy(), C.copy()
    return T, C, maxit, rel


def integrate_until_dry(step, T, C, dt, ct=CTARGET, record=None, nmax=None):
    """推进至烘干判据 max C <= ct 成立，返回 (t*, T*, C*, n*)。

    step(T, C, t_prev, t_cur)：已绑定几何与 dt 的一步推进，返回 (T, C)；
    record(n, t, T, C)：可选，每个时间层推进后回调（记录输出用）。
    事件求根：g(t) = max C - ct 在相邻时间层间线性插值求根，末行取 t* 处插值状态。
    依据：判据由中心点控制、C(t) 平滑，线性求根误差 O(dt^2)，远小于空间离散不确定度。
    """
    if nmax is None:
        nmax = int(TMAX / dt)
    for n in range(1, nmax + 1):
        Tprev, Cprev = T.copy(), C.copy()
        T, C = step(T, C, (n - 1) * dt, n * dt)
        if record is not None:
            record(n, n * dt, T, C)
        if C.max() <= ct:
            gp = Cprev.max() - ct
            gn = C.max() - ct
            frac = gp / (gp - gn)
            return (n - 1) * dt + dt * frac, \
                Tprev + frac * (T - Tprev), Cprev + frac * (C - Cprev), n
    raise RuntimeError('%.0f h 内未达烘干判据' % (TMAX / 3600.0))
