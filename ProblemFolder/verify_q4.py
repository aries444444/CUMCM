# -*- coding: utf-8 -*-
"""问题4 验证与增强脚本（评审报告问题4/5/6/14 的复验 + 消融 + 灵敏度）。
复用 q4_solver 的函数（q4_solver 主流程已加 __main__ 守卫,导入不触发求解）。
A) 6h 行网格收敛 N=21/41/81；B) N=81 时间收敛 dt=5 vs 2.5；
C) Picard 两轮耦合残差量化；D) N=81 生产复算(t*、表6、V1 守恒)；
E) t* 网格收敛 21/41/81 + Richardson 外推；F) 消融:附录4 + 固定 R=2cm；
G) 灵敏度:Ta±2、Ca±20%、h±50%、hm±50% 对 t*。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else '.')
import numpy as np
import q4_solver as Q

DT = 5.0
CTARGET = 0.15
DIST = np.array([0.0, 0.005, 0.01, 0.015, 0.02])


def row6h(dz, nz, dt=DT):
    """前 6 h 的表6 行(含表面列),指定网格。"""
    old = (Q.NZ, Q.DZ, Q.zeta)
    Q.NZ, Q.DZ = nz, dz
    Q.zeta = np.arange(nz) * dz
    T, C = np.full(nz, Q.T0), np.full(nz, Q.C0)
    for n in range(1, int(21600 / dt) + 1):
        T, C = Q.step(T, C, (n - 1) * dt, n * dt, dt)
    R = Q.RofT(21600)
    row = [float(np.interp(d / R, Q.zeta, C)) if d <= R else np.nan for d in DIST] + [C[-1]]
    Q.NZ, Q.DZ, Q.zeta = old
    return np.array(row)


def full_run(dz, nz, dt=DT, fixed_R=False, Ta=None, Ca=None, h=None, hm=None,
             tag='', sample=False):
    """全程求解至 max C<=0.15。fixed_R: 消融实验(附录4 物性 + R 恒为 2cm)。
    Ta/Ca/h/hm 传 None 用默认;非 None 时临时覆盖 q4_solver 的模块全局。"""
    old = (Q.NZ, Q.DZ, Q.zeta, Q.TEND, Q.CEND, Q.h, Q.hm)
    Q.NZ, Q.DZ = nz, dz
    Q.zeta = np.arange(nz) * dz
    if Ta is not None:
        Q.TEND = Ta
    if Ca is not None:
        Q.CEND = Ca
    if h is not None:
        Q.h = h
    if hm is not None:
        Q.hm = hm
    R0c = 0.02
    T, C = np.full(nz, Q.T0), np.full(nz, Q.C0)
    tstar = None
    Whist, Cshist, Rhist, Rphist, tt_hist = [], [], [], [], []
    tab6 = []
    nmax = int(300 * 3600 / dt)
    for n in range(1, nmax + 1):
        tn = n * dt
        tp = tn - dt
        tm = 0.5 * (tp + tn)
        R = R0c if fixed_R else Q.RofT(tm)
        Tn, Cn = T.copy(), C.copy()
        for _ in range(2):
            Kk = Q.k4(C); Rc = Q.rho4(C) * Q.cp4(C)
            T = Q.crank_step(Tn, Kk, Rc, R, Q.h * R / Kk[-1], Q.T_air(tp), Q.T_air(tn), dt)
            Dc = Q.D4(C, T + 273.15)
            C = Q.crank_step(Cn, Dc, np.ones(nz), R, Q.hm * R / Dc[-1], Q.C_air(tp), Q.C_air(tn), dt)
        W = Q.water_int(C, R)
        Whist.append(W); Cshist.append(C[-1]); Rhist.append(R); Rphist.append(0.0 if fixed_R else Q.RpofT(tn))
        tt_hist.append(tn)
        if n % int(6 * 3600 / dt) == 0:
            row = [float(np.interp(d / R, Q.zeta, C)) if d <= R else np.nan for d in DIST] + [C[-1]]
            tab6.append(row)
        if C.max() <= CTARGET:
            tstar = tn
            row = [float(np.interp(d / R, Q.zeta, C)) if d <= R else np.nan for d in DIST] + [C[-1]]
            tab6.append(row)
            break
    if tstar is None:
        raise RuntimeError(tag + '300 h 内未达判据')
    # V1 动边界守恒(同 q4_solver 主流程口径)
    W0 = Q.water_int(np.full(nz, Q.C0), R0c)
    lhs = Whist[-1] - W0
    tt = np.array(tt_hist)
    f1 = [2 * rp / r * w for rp, r, w in zip(Rphist, Rhist, Whist)]
    f2 = [-r * Q.hm * (cs - Q.C_air(t)) for r, cs, t in zip(Rhist, Cshist, tt)]
    rhs = np.trapezoid([a + b for a, b in zip(f1, f2)], tt)
    out = dict(tstar=tstar, tab6=np.array(tab6), V1=(lhs, rhs, lhs - rhs))
    if sample:
        out['C'] = C
    Q.NZ, Q.DZ, Q.zeta, Q.TEND, Q.CEND, Q.h, Q.hm = old
    return out


if __name__ == '__main__':
    print('===== A) 6h 行网格收敛 (dt=5s) =====')
    r21 = row6h(1 / 20, 21)
    r41 = row6h(1 / 40, 41)
    r81 = row6h(1 / 80, 81)
    d2141 = np.nanmax(np.abs(r21 - r41))
    d4181 = np.nanmax(np.abs(r41 - r81))
    print('21 vs 41 最大差 %.3e ; 41 vs 81 最大差 %.3e ; 比值 %.2f'
          % (d2141, d4181, d2141 / d4181))
    print('41 网格 6h 行:', np.round(r41, 4))
    print('81 网格 6h 行:', np.round(r81, 4))

    print('===== B) 时间收敛 (N=81, 6h 行, dt=5 vs 2.5) =====')
    t5 = row6h(1 / 80, 81, 5.0)
    t25 = row6h(1 / 80, 81, 2.5)
    print('最大差 %.2e' % np.nanmax(np.abs(t5 - t25)))

    print('===== C) Picard 两轮耦合残差 (N=81, dt=5) =====')
    Q.NZ, Q.DZ = 81, 1 / 80
    Q.zeta = np.arange(81) / 80
    T, C = np.full(81, Q.T0), np.full(81, Q.C0)
    for n in range(1, int(200000 / DT) + 1):
        tn, tp = n * DT, (n - 1) * DT
        tm = 0.5 * (tp + tn)
        R = Q.RofT(tm)
        Tn, Cn = T.copy(), C.copy()
        Kk = Q.k4(C); Rc = Q.rho4(C) * Q.cp4(C)
        T1 = Q.crank_step(Tn, Kk, Rc, R, Q.h * R / Kk[-1], Q.T_air(tp), Q.T_air(tn), DT)
        Dc1 = Q.D4(C, T1 + 273.15)
        C1 = Q.crank_step(Cn, Dc1, np.ones(81), R, Q.hm * R / Dc1[-1], Q.C_air(tp), Q.C_air(tn), DT)
        Kk2 = Q.k4(C1); Rc2 = Q.rho4(C1) * Q.cp4(C1)
        T2 = Q.crank_step(Tn, Kk2, Rc2, R, Q.h * R / Kk2[-1], Q.T_air(tp), Q.T_air(tn), DT)
        Dc2 = Q.D4(C1, T2 + 273.15)
        C2 = Q.crank_step(Cn, Dc2, np.ones(81), R, Q.hm * R / Dc2[-1], Q.C_air(tp), Q.C_air(tn), DT)
        dC_pic = np.max(np.abs(C2 - C1)); dC_step = np.max(np.abs(C2 - Cn))
        dT_pic = np.max(np.abs(T2 - T1)); dT_step = np.max(np.abs(T2 - Tn))
        if tn in (3600, 21600, 72000, 172800):
            print('t=%6ds: Picard |ΔC|=%.3e (步变化 %.3e, 比 %.1e) | |ΔT|=%.3e (步变化 %.3e)'
                  % (tn, dC_pic, dC_step, dC_pic / dC_step, dT_pic, dT_step))
        T, C = T2, C2
        if C.max() <= CTARGET:
            break

    print('===== D) 生产复算 N=81, dt=5 =====')
    res81 = full_run(1 / 80, 81, tag='N=81 生产')
    print('t* = %d s = %.4f h' % (res81['tstar'], res81['tstar'] / 3600))
    print('表6 (N=81):')
    for k, row in enumerate(res81['tab6']):
        rr = ' '.join('  --' if np.isnan(v) else '%6.4f' % v for v in row)
        trow = (k + 1) * 6.0 if k + 1 < len(res81['tab6']) else res81['tstar'] / 3600.0
        print('%7.1fh |%s' % (trow, rr))
    lhs, rhs, diff = res81['V1']
    print('V1 守恒: ΣvΔC=%.6e  ∫(...)dt=%.6e  差 %.2e (相对 %.1e)'
          % (lhs, rhs, diff, abs(diff) / abs(lhs)))

    print('===== E) t* 网格收敛 21/41/81 =====')
    res21 = full_run(1 / 20, 21, tag='N=21')
    res41 = full_run(1 / 40, 41, tag='N=41')
    t21, t41, t81 = res21['tstar'], res41['tstar'], res81['tstar']
    print('t*(21)=%.4f h  t*(41)=%.4f h  t*(81)=%.4f h' % (t21 / 3600, t41 / 3600, t81 / 3600))
    ext = (4 * t81 - t41) / 3
    print('Richardson 外推(41,81): t*≈%.4f h ; (21,41): t*≈%.4f h'
          % (ext / 3600, (4 * t41 - t21) / 3 / 3600))

    print('===== F) 消融实验: 附录4 物性 + 固定 R=2cm (N=81, dt=5) =====')
    resFix = full_run(1 / 80, 81, fixed_R=True, tag='消融固定R')
    print('t*_无收缩 = %.4f h ; t*_有收缩 = %.4f h ; 收缩贡献 Δt = %.4f h (%.1f%%)'
          % (resFix['tstar'] / 3600, res81['tstar'] / 3600,
             (resFix['tstar'] - res81['tstar']) / 3600,
             (resFix['tstar'] - res81['tstar']) / resFix['tstar'] * 100))

    print('===== G) 灵敏度 (N=81, dt=5) =====')
    sens = {}
    for tag, kw in [('基线', {}), ('Ta+2', dict(Ta=52.165)), ('Ta-2', dict(Ta=48.165)),
                    ('Ca+20%', dict(Ca=0.04986 * 1.2)), ('Ca-20%', dict(Ca=0.04986 * 0.8)),
                    ('h+50%', dict(h=37.5)), ('h-50%', dict(h=12.5)),
                    ('hm+50%', dict(hm=1.2e-6)), ('hm-50%', dict(hm=4e-7))]:
        if tag == '基线':
            sens[tag] = res81['tstar']
        else:
            sens[tag] = full_run(1 / 80, 81, tag=tag, **kw)['tstar']
        print('%-8s t* = %.4f h (Δ %.3f h)' % (tag, sens[tag] / 3600,
              (sens[tag] - sens['基线']) / 3600))
