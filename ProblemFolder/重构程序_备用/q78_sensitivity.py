# -*- coding: utf-8 -*-
"""2026 CUMCM A 题 表8 参数灵敏度重算（Picard<1e-8、t* 判据事件求根）

问题3: 径向离散 0.25 mm（N=81）、dt=2.5 s，附录3 物性、固定 R；
问题4: ζ 离散 1/80（N=81）、dt=5 s，附录4 物性、附件2 收缩 R(t)。
变体: Ta±2 °C、Ca±20%、h±50%、hm±50%（相对题给值），其余参数保持题给值。
公共数值内核见 common.py。

依据: 灵敏度取两变体 t* 之差，对空间离散不敏感，用验证离散即可（与 7.7 节说明一致）；
迭代判据改变引起的 t* 系统性偏移（实测 <0.001 h）在差分中相消，故本表只需重算一遍。
"""
from functools import partial
import numpy as np
from scipy.interpolate import PchipInterpolator

from common import (crank_zeta, integrate_until_dry, kernel, load_annex,
                    picard_step, setup_console)

setup_console()

R0, T0, C0 = 0.02, 28.0, 2.55
TEND, CEND = 50.165, 0.04986

# ---- 附件数据 ----
d1 = load_annex('附件1.xlsx')
d2 = load_annex('附件2.xlsx')
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

def crank_step(u, K, Cp, beta, bc_prev, bc_next, dt, NZ, DZ, R):
    """径向（R 固定）或 ζ 系（R 为当前半径）的 Crank--Nicolson 一步，见 common.crank_zeta。"""
    return crank_zeta(u, K, Cp, beta, bc_prev, bc_next, dt, NZ, DZ, R)

def step(T, C, t_prev, t_cur, dt, NZ, DZ, R, rhoF, cpF, kF, DF):
    """交错 Picard 耦合步（变体 h、h_m 经模块级 P 传入）。"""
    crk = partial(crank_zeta, dt=dt, NZ=NZ, DZ=DZ, R=R)
    return picard_step(crk, T, C, dt, t_prev, t_cur, rhoF, cpF, kF, DF,
                       T_air, C_air, P['h'] * R, P['hm'] * R)

def run_case(appendix, shrink, NZ, DZ, dt):
    """返回 t*（s）。shrink=False 时 R≡R0；True 时 R(t) 取半步中点。"""
    PF = (rho3, cp3, k3, D3) if appendix == 3 else (rho4, cp4, k4, D4)
    RofT = (lambda t: float(_Rp(t)) if t <= T2END else REND) if shrink else (lambda t: R0)
    T, C = np.full(NZ, T0), np.full(NZ, C0)

    def one_step(T, C, t_prev, t_cur):
        return step(T, C, t_prev, t_cur, dt, NZ, DZ, RofT(0.5 * (t_prev + t_cur)), *PF)[:2]

    tstar, _, _, _ = integrate_until_dry(one_step, T, C, dt)
    return tstar

VARIANTS = [('基准', dict()),
            ('Ta+2', dict(dToff=2.0)), ('Ta-2', dict(dToff=-2.0)),
            ('Ca+20%', dict(cscale=1.2)), ('Ca-20%', dict(cscale=0.8)),
            ('h+50%', dict(h=37.5)), ('h-50%', dict(h=12.5)),
            ('hm+50%', dict(hm=1.2e-6)), ('hm-50%', dict(hm=4e-7))]

if __name__ == '__main__':
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
