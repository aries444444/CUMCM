# -*- coding: utf-8 -*-
"""2026 CUMCM A 题 表8 参数灵敏度重算（新求解器：Picard<1e-8、t* 根插值）

问题3: 径向离散 0.25 mm（N=81）、dt=2.5 s，附录3 物性、固定 R；
问题4: ζ 离散 1/80（N=81）、dt=5 s，附录4 物性、附件2 收缩 R(t)。
变体: Ta±2 °C、Ca±20%、h±50%、hm±50%（相对题给值），其余参数保持题给值。

依据: 灵敏度取两变体 t* 之差，对空间离散不敏感，用验证离散即可（与 7.7 节说明一致）；
迭代判据改变引起的 t* 系统性偏移（实测 <0.001 h）在差分中相消，故本表只需重算一遍。
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

def _annex(name):
    """定位赛题附件数据:优先脚本同目录(支撑材料布局),其次仓库目录结构。"""
    for p in (os.path.join(ROOT, name),
              os.path.join(ROOT, 'OriginalMaterial', 'A题', '附件', name)):
        if os.path.exists(p):
            return p
    raise FileNotFoundError('未找到 %s:请将赛题附件文件放在脚本同目录下' % name)

ANNEX1 = _annex('附件1.xlsx')
ANNEX2 = _annex('附件2.xlsx')

R0, T0, C0 = 0.02, 28.0, 2.55
TEND, CEND = 50.165, 0.04986
CTARGET = 0.15
PICARD_TOL, PICARD_MAXIT = 1e-8, 20

# ---- 附件数据 ----
d1 = np.array([row for row in
               openpyxl.load_workbook(ANNEX1, data_only=True)['Sheet1'].iter_rows(values_only=True)][1:],
              float)
d2 = np.array([row for row in
               openpyxl.load_workbook(ANNEX2, data_only=True)['Sheet1'].iter_rows(values_only=True)][1:],
              float)
T2END = d2[-1, 0]; REND = d2[-1, 1] / 100.0
_Rp = PchipInterpolator(d2[:, 0], d2[:, 1] / 100.0, extrapolate=False)

# ---- 变体参数（模块级，运行中修改）----
P = dict(dToff=0.0, cscale=1.0, h=25.0, hm=8e-7)
T_air = lambda t: (np.interp(t, d1[:, 0], d1[:, 1]) + P['dToff']) if t <= 14400 else TEND + P['dToff']
C_air = lambda t: (np.interp(t, d1[:, 0], d1[:, 2]) * P['cscale']) if t <= 14400 else CEND * P['cscale']

# 附录3 / 附录4 物性
rho3 = lambda C: 650.0 + 128.0 * C
cp3 = lambda C: 1450.0 + 2736.0 * C / (C + 1.0)
k3 = lambda C: 0.21 + 0.38 * C / (C + 1.0)
D3 = lambda C, T: 2.4e-3 * np.exp(-0.45 / C) * np.exp(-3850.0 / T)
rho4 = lambda C: 760.0 + 90.0 * C
cp4 = lambda C: 1850.0 + 2150.0 * C / (C + 1.0)
k4 = lambda C: 0.12 + 0.20 * C / (C + 1.0)
D4 = lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 0.02)) * np.exp(-3850.0 / T)

def thomas(a, b, c, rhs):
    n = len(rhs)
    ab = np.zeros((3, n))
    ab[0, 1:] = c[:-1]
    ab[1, :] = b
    ab[2, :-1] = a[1:]
    return solve_banded((1, 1), ab, rhs)

def crank_step(u, K, Cp, beta, bc_prev, bc_next, dt, NZ, DZ, R):
    """径向（R 固定）或 ζ 系（R 为当前半径）的 Crank--Nicolson 一步，表面节点在真实表面。"""
    Km = (K[:-1] + K[1:]) / 2
    R2 = R * R
    lo = np.zeros(NZ); di = np.zeros(NZ); up = np.zeros(NZ)
    di[0] = -4 * Km[0] / (DZ ** 2 * R2); up[0] = 4 * Km[0] / (DZ ** 2 * R2)
    for i in range(1, NZ - 1):
        lo[i] = Km[i - 1] * (i - 0.5) / (i * DZ ** 2 * R2)
        di[i] = -(Km[i] * (i + 0.5) + Km[i - 1] * (i - 0.5)) / (i * DZ ** 2 * R2)
        up[i] = Km[i] * (i + 0.5) / (i * DZ ** 2 * R2)
    vhat = (1.0 - ((NZ - 1.5) * DZ) ** 2) / 2.0
    Bc = beta * K[-1]
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

def step(T, C, t_prev, t_cur, dt, NZ, DZ, R, rhoF, cpF, kF, DF):
    """交错 Picard，T、C 相对变化 < 1e-8 或满 20 轮停止。"""
    Tn, Cn = T.copy(), C.copy()
    rel = 0.0
    for it in range(1, PICARD_MAXIT + 1):
        Kk = kF(C); Rc = rhoF(C) * cpF(C)
        T = crank_step(Tn, Kk, Rc, P['h'] * R / Kk[-1], T_air(t_prev), T_air(t_cur), dt, NZ, DZ, R)
        Tk = T + 273.15
        Dc = DF(C, Tk)
        C = crank_step(Cn, Dc, np.ones(NZ), P['hm'] * R / Dc[-1], C_air(t_prev), C_air(t_cur), dt, NZ, DZ, R)
        if it > 1:
            rel = max(np.max(np.abs(T - To)) / max(np.max(np.abs(T)), 1.0),
                      np.max(np.abs(C - Co)) / max(np.max(np.abs(C)), 1.0))
            if rel < PICARD_TOL:
                return T, C, it, rel
        To, Co = T.copy(), C.copy()
    return T, C, PICARD_MAXIT, rel

def run_case(appendix, shrink, NZ, DZ, dt):
    """返回 t*（s）。shrink=False 时 R≡R0；True 时 R(t) 取半步中点。"""
    PF = (rho3, cp3, k3, D3) if appendix == 3 else (rho4, cp4, k4, D4)
    RofT = (lambda t: float(_Rp(t)) if t <= T2END else REND) if shrink else (lambda t: R0)
    T, C = np.full(NZ, T0), np.full(NZ, C0)
    nmax = int(300 * 3600 / dt)
    for n in range(1, nmax + 1):
        tm = (n - 0.5) * dt
        Cprev = C.copy()
        T, C, it, rel = step(T, C, (n - 1) * dt, n * dt, dt, NZ, DZ, RofT(tm), *PF)
        if C.max() <= CTARGET:
            gp = Cprev.max() - CTARGET
            gn = C.max() - CTARGET
            return (n - 1) * dt + dt * gp / (gp - gn)
    raise RuntimeError('300 h 内未达判据')

VARIANTS = [('基准', dict()),
            ('Ta+2', dict(dToff=2.0)), ('Ta-2', dict(dToff=-2.0)),
            ('Ca+20%', dict(cscale=1.2)), ('Ca-20%', dict(cscale=0.8)),
            ('h+50%', dict(h=37.5)), ('h-50%', dict(h=12.5)),
            ('hm+50%', dict(hm=1.2e-6)), ('hm-50%', dict(hm=4e-7))]

if __name__ == '__main__':
    # ===== 表7: 物性×几何 2×2 消融（同一验证离散 Δζ=1/160, Δt=5 s；与 7.5 节一致） =====
    print('==== 表7 消融（Δζ=1/160, Δt=5 s） ====')
    ts = {}
    for name, ap, sh in (('附录3×固定', 3, False), ('附录3×收缩', 3, True),
                         ('附录4×固定', 4, False), ('附录4×收缩', 4, True)):
        ts[name] = run_case(ap, sh, 161, 1.0 / 160.0, 5.0)
        print('%s: t* = %.4f h (%d s)' % (name, ts[name] / 3600, round(ts[name])))
    print('收缩净效应: 附录3 Δt* = %.4f h（%.1f%%），附录4 Δt* = %.4f h（%.1f%%）'
          % ((ts['附录3×收缩'] - ts['附录3×固定']) / 3600,
             (ts['附录3×收缩'] / ts['附录3×固定'] - 1) * 100,
             (ts['附录4×收缩'] - ts['附录4×固定']) / 3600,
             (ts['附录4×收缩'] / ts['附录4×固定'] - 1) * 100))
    print()
    for label, cfg in (('问题3（ζ 离散 1/80 ⇒ 物理 0.25 mm, dt=2.5 s）',
                        dict(appendix=3, shrink=False, NZ=81, DZ=1.0/80.0, dt=2.5)),
                       ('问题4（Δζ=1/80, dt=5 s）',
                        dict(appendix=4, shrink=True, NZ=81, DZ=1.0/80.0, dt=5.0))):
        print('==== %s ====' % label)
        rows = []
        base = None
        for name, v in VARIANTS:
            for k in P:  # 重置为基准再叠加变体
                P[k] = dict(dToff=0.0, cscale=1.0, h=25.0, hm=8e-7)[k]
            P.update(v)
            t = run_case(**cfg)
            rows.append((name, t))
            if name == '基准':
                base = t
        for name, t in rows:
            print('%s: t* = %.4f h, Δt* = %+.4f h（%+.1f%%）'
                  % (name, t / 3600, (t - base) / 3600, (t / base - 1) * 100))
        print()
