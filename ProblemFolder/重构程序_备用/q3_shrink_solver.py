# -*- coding: utf-8 -*-
"""2026 CUMCM A 题 消融对照求解器（审稿意见 6 的 2×2 表）

同一离散设置（ζ 网格 Δζ=1/160、Δt=5 s、Picard 相对残差 < 1e-8）下，
物性（附录3 / 附录4）x 几何（固定半径 R=2 cm / 附件2 收缩 R(t)）四个组合的烘干时长：

             固定半径      收缩半径
  附录3      t3,F          t3,S
  附录4      t4,F          t4,S

依据：问题 3 与问题 4 同时更换了物性与几何，只有同物性比较固定/收缩半径，
才能分离收缩的净效应（审稿意见 6）。固定半径情形 ζ 系方程即径向方程（R≡R0）。
"""
from functools import partial
import numpy as np
from scipy.interpolate import PchipInterpolator

from common import (crank_zeta, integrate_until_dry, kernel, load_annex,
                    picard_step, setup_console)

setup_console()

R0, T0, C0 = 0.02, 28.0, 2.55
h, hm = 25.0, 8e-7
TEND, CEND = 50.165, 0.04986
DZ, NZ = 1.0 / 160.0, 161
zeta = np.arange(NZ) * DZ
DT = 5.0

rho3 = lambda C: 650.0 + 128.0 * C
cp3 = lambda C: 1450.0 + 2736.0 * C / (C + 1.0)
k3 = lambda C: 0.21 + 0.38 * C / (C + 1.0)
D3 = lambda C, T: 2.4e-3 * np.exp(-0.45 / C) * np.exp(-3850.0 / T)
rho4 = lambda C: 760.0 + 90.0 * C
cp4 = lambda C: 1850.0 + 2150.0 * C / (C + 1.0)
k4 = lambda C: 0.12 + 0.20 * C / (C + 1.0)
D4 = lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 0.02)) * np.exp(-3850.0 / T)

d1 = load_annex('附件1.xlsx')
T_air = lambda t: np.interp(t, d1[:, 0], d1[:, 1]) if t <= 14400 else TEND
C_air = lambda t: np.interp(t, d1[:, 0], d1[:, 2]) if t <= 14400 else CEND

d2 = load_annex('附件2.xlsx')
T2END = d2[-1, 0]; REND = d2[-1, 1] / 100.0
_Rp = PchipInterpolator(d2[:, 0], d2[:, 1] / 100.0, extrapolate=False)

def run_case(props, shrink):
    """props: 'A3' 或 'A4'；shrink: True 用附件2 R(t)，False 固定 R=R0。返回 (t*, 统计)。"""
    rhoF, cpF, kF, DF = (rho3, cp3, k3, D3) if props == 'A3' else (rho4, cp4, k4, D4)
    RofT = (lambda t: float(_Rp(t)) if t <= T2END else REND) if shrink else (lambda t: R0)
    T, C = np.full(NZ, T0), np.full(NZ, C0)
    iters, rels = [], []

    def one_step(T, C, t_prev, t_cur):
        R = RofT(0.5 * (t_prev + t_cur))
        crk = partial(crank_zeta, dt=DT, NZ=NZ, DZ=DZ, R=R)
        T, C, it, rel = picard_step(crk, T, C, DT, t_prev, t_cur, rhoF, cpF, kF, DF,
                                    T_air, C_air, h * R, hm * R)
        iters.append(it); rels.append(rel)
        return T, C

    tstar, _, _, _ = integrate_until_dry(one_step, T, C, DT)
    return tstar, np.array(iters), np.array(rels)

if __name__ == "__main__":
    print('2×2 消融对照（Δζ=1/160, Δt=5 s, Picard<1e-8）:')
    t3F, i3F, r3F = run_case('A3', False)
    t3S, i3S, r3S = run_case('A3', True)
    t4F, i4F, r4F = run_case('A4', False)
    t4S, i4S, r4S = run_case('A4', True)
    for name, t, it, rl in (('t3,F', t3F, i3F, r3F), ('t3,S', t3S, i3S, r3S),
                            ('t4,F', t4F, i4F, r4F), ('t4,S', t4S, i4S, r4S)):
        print('%s = %.4f h (%d s), 平均 %.2f 轮 / 最大 %d 轮'
              % (name, t / 3600, round(t), it.mean(), it.max()))
    print('收缩净效应: 附录3 下 Δt* = %.4f h（%.1f%%），附录4 下 Δt* = %.4f h（%.1f%%）'
          % ((t3S - t3F) / 3600, (t3S / t3F - 1) * 100, (t4S - t4F) / 3600, (t4S / t4F - 1) * 100))
