# -*- coding: utf-8 -*-
"""2026 CUMCM A 题 问题1 求解器
圆柱轴对称热传导 + 水分扩散（非线性 D(C)），守恒型有限体积离散 + Crank--Nicolson + Thomas。
输出：表1/表2（控制台）、result1.xlsx（两个 sheet，1800 行 x 21 列，保留 4 位小数）、
四重验证（Bessel 解析锚点 / 离散收敛 / 质量守恒 / 渐近行为）。
"""
import sys, os
import numpy as np
import openpyxl
from scipy.special import j0, j1
from scipy.optimize import brentq

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

ROOT = os.path.dirname(os.path.abspath(__file__))

def _annex(name):
    """定位赛题附件数据:优先脚本同目录(支撑材料布局),其次仓库目录结构。"""
    for p in (os.path.join(ROOT, name),
              os.path.join(ROOT, 'OriginalMaterial', 'A题', '附件', name)):
        if os.path.exists(p):
            return p
    raise FileNotFoundError('未找到 %s:请将赛题附件文件放在脚本同目录下' % name)

ANNEX1 = _annex('附件1.xlsx')
ANNEX2 = _annex('附件2.xlsx')
OUTXLSX = os.path.join(ROOT, 'result1.xlsx')

# ---------- 参数（附录2） ----------
R, dr, N = 0.02, 3.125e-5, 641     # 半径 2 cm，主离散空间步长 0.03125 mm，641 个节点
                                    # （0.125 mm 下 (100s,表面) 误差约 4.7e-4，加密至 0.03125 mm 后约 3e-5，四位小数可靠）
r = np.arange(N) * dr              # 节点半径 r_i = i*dr（r_{N-1} = R 为表面节点）
rho, cp, k, h, hm = 820.0, 2600.0, 0.36, 25.0, 8e-7
alpha = k / (rho * cp)             # 热扩散率
T0, C0 = 28.0, 2.55                # 初值
DofC = lambda C: 7e-9 * np.exp(-0.89 / C)   # 附录2：指数为 -0.89/C

# ---------- 附件1：烘房环境（线性插值） ----------
d = np.array([row for row in
              openpyxl.load_workbook(ANNEX1, data_only=True)['Sheet1'].iter_rows(values_only=True)][1:],
             float)
T_air = lambda t: np.interp(t, d[:, 0], d[:, 1])
C_air = lambda t: np.interp(t, d[:, 0], d[:, 2])

def thomas(a, b, c, rhs):
    n = len(rhs); cp = np.zeros(n); dp = np.zeros(n)
    cp[0] = c[0] / b[0]; dp[0] = rhs[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]; cp[i] = c[i] / m; dp[i] = (rhs[i] - a[i] * dp[i - 1]) / m
    x = np.zeros(n); x[-1] = dp[-1]
    for i in range(n - 2, -1, -1): x[i] = dp[i] - cp[i] * x[i + 1]
    return x

def crank_step(u, K, beta, bc_prev, bc_next, dt):
    """守恒型空间离散 + Crank--Nicolson 一步。
    K: 每节点扩散系数数组（热: alpha；质: D(C_i)）；beta: 表面 Robin 系数（h/k 或 h_m/D_s）。
    离散算子: (1/r) d(rK du/dr)/dr，r=0 用对称性+洛必达（4(u1-u0)/dr^2），
    表面控制体 [r_{N-3/2}, R] 上以 Robin 边界条件重构真实表面值 u_R，通量严格守恒。"""
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
    vhat = ((N - 1) ** 2 - (N - 1.5) ** 2) * dr ** 2 / 2  # 控制体体积 ∫_{r_{N-3/2}}^R r dr
    Bc = (N - 1) * dr * beta * K[-1]                      # = R·beta·K_s（热: R·h；质: R·h_m）
    lo[N - 1] = Km[N - 2] * (N - 1.5) / vhat
    di[N - 1] = -(lo[N - 1] + Bc / vhat)
    up[N - 1] = 0.0
    def src(bc):
        s = np.zeros(N); s[N - 1] = Bc * bc / vhat; return s
    lhs_lo = -0.5 * dt * lo; lhs_di = 1 - 0.5 * dt * di; lhs_up = -0.5 * dt * up
    rhs = u + 0.5 * dt * (lo * np.roll(u, 1) + di * u + up * np.roll(u, -1)) \
        + 0.5 * dt * (src(bc_next) + src(bc_prev))
    rhs[0] = u[0] + 0.5 * dt * (di[0] * u[0] + up[0] * u[1]) \
        + 0.5 * dt * (src(bc_next)[0] + src(bc_prev)[0])
    rhs[N - 1] = u[N - 1] + 0.5 * dt * (lo[N - 1] * u[N - 2] + di[N - 1] * u[N - 1]) \
        + 0.5 * dt * (src(bc_next)[N - 1] + src(bc_prev)[N - 1])
    return thomas(lhs_lo, lhs_di, lhs_up, rhs)

def solve(nt, dt, Tair, Cair, Kheat=None, DofC=DofC, record=None):
    """求解 0..nt*dt 秒。返回 (T(t_ev, r_ev), C(t_ev, r_ev), T_full, C_full)。
    record: dict(t -> 半径位置列表)。"""
    Kheat = np.full(N, alpha) if Kheat is None else Kheat
    T = np.full(N, T0); C = np.full(N, C0)
    rec = {tt: [] for tt in record}
    for n in range(1, nt + 1):
        t_prev, t_cur = (n - 1) * dt, n * dt
        T = crank_step(T, Kheat, h / k, Tair(t_prev), Tair(t_cur), dt)
        Dc = DofC(C)
        C = crank_step(C, Dc, hm / Dc[-1], Cair(t_prev), Cair(t_cur), dt)
        if n * dt in record:
            rec[n * dt].append((T.copy(), C.copy()))
    Tout = np.array([[rec[tt][0][0][i] for i in record[tt]] for tt in record], float)
    Cout = np.array([[rec[tt][0][1][i] for i in record[tt]] for tt in record], float)
    return Tout, Cout, T, C

if __name__ == "__main__":
    # ---------- 主计算：附件1 时变边界，0--1800 s，Δr=0.03125 mm，Δt=0.125 s ----------
    # 表2 表面处存在陡峭的浓度边界层，0.125 mm 离散下 (100s,表面) 误差约 4.7e-4，
    # 故主离散取 0.03125 mm（误差约 3e-5，四位小数可靠）；result1.xlsx 仍按模板输出每 1 s、每 0.1 cm。
    t_ev = [100, 300, 600, 900, 1200, 1500, 1800]
    ri = [0, 160, 320, 480, 640]     # 0, 0.5, 1, 1.5, 2 cm
    record = {tt: ri for tt in t_ev}
    Ttab, Ctab, T, C = solve(14400, 0.125, T_air, C_air, record=record)
    np.set_printoptions(precision=4, suppress=True)
    print('表1 温度 (°C)，行=100..1800 s，列=0/0.5/1/1.5/2 cm:')
    print(np.round(Ttab, 4))
    print('表2 水分浓度 (kg/kg):')
    print(np.round(Ctab, 4))
    print('烘房@1800s: T=%.4f C=%.5f' % (T_air(1800), C_air(1800)))
    print('D(C) 表面值: t=0 -> %.4e, t=1800s -> %.4e' % (DofC(C0), DofC(C[-1])))

    # ---------- 写出 result1.xlsx（与附件3模板一致：每 1 s 一行、每 0.1 cm 一列） ----------
    wb = openpyxl.Workbook()
    wsT = wb.active; wsT.title = '温度'
    wsC = wb.create_sheet('水分浓度')
    xcols = list(range(0, N, 32))     # 0:0.1:2 cm（0.03125 mm 离散上每 32 节点一列）
    for ws in (wsT, wsC):
        ws.cell(1, 1, '时间\\到药材中心的距离')
        for j, i in enumerate(xcols):
            ws.cell(1, 2 + j, round(r[i] * 100, 1))
    T, C = np.full(N, T0), np.full(N, C0)
    Kheat = np.full(N, alpha)
    for n in range(1, 14401):
        T = crank_step(T, Kheat, h / k, T_air((n - 1) * 0.125), T_air(n * 0.125), 0.125)
        Dc = DofC(C)
        C = crank_step(C, Dc, hm / Dc[-1], C_air((n - 1) * 0.125), C_air(n * 0.125), 0.125)
        if n % 8 == 0:                          # 每 1 s 记一行
            m = n // 8
            wsT.cell(1 + m, 1, m); wsC.cell(1 + m, 1, m)   # 时间 1..1800 s
            for j, i in enumerate(xcols):
                wsT.cell(1 + m, 2 + j, round(float(T[i]), 4))
                wsC.cell(1 + m, 2 + j, round(float(C[i]), 4))
    wb.save(OUTXLSX)
    print('已写出 %s' % OUTXLSX)

    # ---------- 验证1：Bessel 解析锚点（常数边界 T_air=41.513，主离散 Δr=0.03125 mm） ----------
    Ta = 41.513
    Bi = h * R / k
    f = lambda b: b * j1(b) - Bi * j0(b)
    bs = [brentq(f, lo, hi) for lo, hi in
          zip(np.arange(0.05, 400, 0.5), np.arange(0.55, 400.5, 0.5)) if f(lo) * f(hi) < 0][:40]
    bs = np.array(bs)
    An = 2 * j1(bs) / (bs * (j0(bs) ** 2 + j1(bs) ** 2))
    def T_exact(rr, t):
        return Ta + (T0 - Ta) * (An * np.exp(-alpha * bs ** 2 * t / R ** 2) * j0(bs * rr / R)).sum()
    Tb, Cb, _, _ = solve(14400, 0.125, lambda t: Ta, lambda t: 0.04986, record=record)
    print('验证1 Bessel: 中心 数值=%.4f 解析=%.4f | 表面 数值=%.4f 解析=%.4f'
          % (Tb[-1][0], T_exact(0.0, 1800), Tb[-1][4], T_exact(0.02, 1800)))

    # ---------- 验证2：离散收敛（0.125 mm 与 0.0625 mm 两套粗离散，Δt=0.125 s 固定，与主离散同位置采样） ----------
    global_dr = dr; global_N = N
    dr = 1.25e-4; N = 161
    ri2 = [0, 40, 80, 120, 160]            # 0.125 mm 离散上的 0 / 0.5 / 1 / 1.5 / 2 cm
    T2, C2, T2f, C2f = solve(14400, 0.125, T_air, C_air, record={tt: ri2 for tt in t_ev})
    dr = 6.25e-5; N = 321
    ri1 = [0, 80, 160, 240, 320]           # 0.0625 mm 离散上的同一物理位置
    T1, C1, T1f, C1f = solve(14400, 0.125, T_air, C_air, record={tt: ri1 for tt in t_ev})
    dr = global_dr; N = global_N
    print('验证2 收敛: 表1 最大差 0.125mm->0.0625mm: %.2e °C, 0.0625mm->0.03125mm: %.2e °C'
          % (np.max(np.abs(T2 - T1)), np.max(np.abs(T1 - Ttab))))
    print('         表2 最大差 0.125mm->0.0625mm: %.2e, 0.0625mm->0.03125mm: %.2e kg/kg'
          % (np.max(np.abs(C2 - C1)), np.max(np.abs(C1 - Ctab))))

    # ---------- 验证3：离散质量守恒（0.125 mm 与 0.03125 mm 两套离散上分别做） ----------
    # 恒等式 Σ_i v_i(C^f_i-C^0_i) = -R∫h_m(C_R-C_air)dt（v_i 为控制体体积，C_R 即表面节点值）
    # 由格式构造逐时间层成立，数值上应闭合至机器精度。
    def cv_volumes(dr, N):
        v = np.zeros(N)
        v[0] = dr ** 2 / 8
        v[1:N - 1] = np.arange(1, N - 1) * dr ** 2
        v[N - 1] = ((N - 1) ** 2 - (N - 1.5) ** 2) * dr ** 2 / 2
        return v

    def check_conserve(Cf, drg, Ng, dt, nt):
        global N, dr
        v = cv_volumes(drg, Ng)
        lhs = (Cf * v).sum() - C0 * v.sum()
        N_save, dr_save = N, dr
        N, dr = Ng, drg                       # crank_step 读取模块级 N, dr
        Cchk = np.full(Ng, C0)
        Fhist = [hm * (C0 - C_air(0.0))]      # F_n = h_m(C^n_R - C_air(t_n))，C_R 即表面节点
        for n in range(1, nt + 1):
            Dc = DofC(Cchk)
            Cchk = crank_step(Cchk, Dc, hm / Dc[-1], C_air((n - 1) * dt), C_air(n * dt), dt)
            Fhist.append(hm * (Cchk[-1] - C_air(n * dt)))
        N, dr = N_save, dr_save
        rhs = -R * np.trapezoid(Fhist, np.arange(0, nt + 1) * dt)
        return lhs, rhs

    lhs2, rhs2 = check_conserve(C2f, 1.25e-4, 161, 0.125, 14400)
    lhs3, rhs3 = check_conserve(C, 3.125e-5, 641, 0.125, 14400)   # 主离散
    print('验证3 离散守恒: 0.125mm 离散 ΣvΔC=%.6e vs -R∫h_m(C_R-C_air)dt=%.6e，差 %.2e'
          % (lhs2, rhs2, lhs2 - rhs2))
    print('              0.03125mm 离散 ΣvΔC=%.6e vs -R∫h_m(C_R-C_air)dt=%.6e，差 %.2e'
          % (lhs3, rhs3, lhs3 - rhs3))

    # ---------- 验证4：渐近行为（常数边界，长时程，1 mm 离散） ----------
    Ta, Ca = 41.513, 0.04986
    ri1mm = [0, 5, 10, 15, 20]            # 1 mm 离散上的 0 / 0.5 / 1 / 1.5 / 2 cm
    rec4 = {20000: ri1mm, 50000: ri1mm, 100000: ri1mm}
    dr = 1e-3; N = 21
    r = np.arange(N) * dr
    T4, C4, _, _ = solve(100000, 1.0, lambda t: Ta, lambda t: Ca, record=rec4)
    dr = global_dr; N = global_N
    r = np.arange(N) * dr
    for k, tt in enumerate(rec4):
        print('验证4 渐近: t=%ds 中心 T=%.4f (目标 %.4f)，表面 C=%.5f (目标 %.5f)'
              % (tt, T4[k][0], Ta, C4[k][4], Ca))

    # ---------- 附件2 收缩数据核对（支撑问题1忽略收缩的假设） ----------
    d2 = np.array([row for row in
                   openpyxl.load_workbook(ANNEX2, data_only=True)['Sheet1'].iter_rows(values_only=True)][1:],
                  float)
    print('附件2: t=0 半径=%.4f cm，t=1800s 半径=%.4f cm，末点 t=%ds 半径=%.4f cm'
          % (d2[0, 1], d2[1, 1], d2[-1, 0], d2[-1, 1]))
