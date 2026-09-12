# -*- coding: utf-8 -*-
"""问题3 验证与增强脚本(评审报告问题2/14 的补充):
A) Picard 两轮耦合残差量化；B) t* 三网格(0.5/0.25/0.125 mm)+ Richardson 外推；
C) 灵敏度:Ta±2°C、Ca±20%、h±50%、hm±50% 对 t*。
复用 q3_solver 的函数(q3_solver 主流程已加 __main__ 守卫)。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else '.')
import numpy as np
import q3_solver as Q

DT3 = 2.5
CTARGET = 0.15


def full_run(drg, ng, dt=DT3, Ta=None, Ca=None, h=None, hm=None, tag=''):
    """全程求解至 max C<=0.15,返回 t*。参数覆盖用模块全局临时替换。"""
    old = (Q.TEND, Q.CEND, Q.h, Q.hm)
    if Ta is not None:
        Q.TEND = Ta
    if Ca is not None:
        Q.CEND = Ca
    if h is not None:
        Q.h = h
    if hm is not None:
        Q.hm = hm
    tstar, _, _, _ = Q.run_drying(drg, ng, dt, sample6h=False)
    Q.TEND, Q.CEND, Q.h, Q.hm = old
    return tstar


if __name__ == '__main__':
    print('===== A) Picard 两轮耦合残差 (0.25mm 网格, dt=2.5) =====')
    Q.N, Q.dr, Q.r = 81, 2.5e-4, np.arange(81) * 2.5e-4
    T, C = np.full(81, Q.T0), np.full(81, Q.C0)
    for n in range(1, int(60000 / DT3) + 1):
        tn, tp = n * DT3, (n - 1) * DT3
        Tn, Cn = T.copy(), C.copy()
        Kk = Q.k23(C); Rc = Q.rho23(C) * Q.cp23(C)
        T1 = Q.crank_step(Tn, Kk, Rc, Q.h / Kk[-1], Q.T_air(tp), Q.T_air(tn), DT3)
        Dc1 = Q.D23(C, T1 + 273.15)
        C1 = Q.crank_step(Cn, Dc1, np.ones(81), Q.hm / Dc1[-1], Q.C_air(tp), Q.C_air(tn), DT3)
        Kk2 = Q.k23(C1); Rc2 = Q.rho23(C1) * Q.cp23(C1)
        T2 = Q.crank_step(Tn, Kk2, Rc2, Q.h / Kk2[-1], Q.T_air(tp), Q.T_air(tn), DT3)
        Dc2 = Q.D23(C1, T2 + 273.15)
        C2 = Q.crank_step(Cn, Dc2, np.ones(81), Q.hm / Dc2[-1], Q.C_air(tp), Q.C_air(tn), DT3)
        dC_pic = np.max(np.abs(C2 - C1)); dC_step = np.max(np.abs(C2 - Cn))
        dT_pic = np.max(np.abs(T2 - T1))
        if tn in (3600, 21600, 129600, 190800):
            print('t=%6ds: Picard |ΔC|=%.3e (步变化 %.3e, 比 %.1e) | |ΔT|=%.3e'
                  % (tn, dC_pic, dC_step, dC_pic / dC_step, dT_pic))
        T, C = T2, C2
        if C.max() <= CTARGET:
            break

    print('===== B) t* 三网格收敛 + 外推 =====')
    t05 = full_run(5e-4, 41, tag='0.5mm')
    t025 = full_run(2.5e-4, 81, tag='0.25mm')
    t0125 = full_run(1.25e-4, 161, tag='0.125mm')
    print('t*(0.5mm)=%.4f h  t*(0.25mm)=%.4f h  t*(0.125mm)=%.4f h'
          % (t05 / 3600, t025 / 3600, t0125 / 3600))
    ext = (4 * t025 - t05) / 3
    ext2 = (4 * t0125 - t025) / 3
    print('Richardson(0.5,0.25): %.4f h ; (0.25,0.125): %.4f h ; |ext2-ext|=%.4f h'
          % (ext / 3600, ext2 / 3600, abs(ext2 - ext) / 3600))

    print('===== C) 灵敏度 (0.25mm, dt=2.5) =====')
    base = t025
    for tag, kw in [('Ta+2', dict(Ta=52.165)), ('Ta-2', dict(Ta=48.165)),
                    ('Ca+20%', dict(Ca=0.04986 * 1.2)), ('Ca-20%', dict(Ca=0.04986 * 0.8)),
                    ('h+50%', dict(h=37.5)), ('h-50%', dict(h=12.5)),
                    ('hm+50%', dict(hm=1.2e-6)), ('hm-50%', dict(hm=4e-7))]:
        ts = full_run(2.5e-4, 81, tag=tag, **kw)
        print('%-8s t* = %.4f h (Δ %+.3f h, %+.1f%%)' % (tag, ts / 3600, (ts - base) / 3600,
              (ts - base) / base * 100))
