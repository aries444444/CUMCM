# -*- coding: utf-8 -*-
"""问题1 补充验证:表2 各单元格 0.125/0.0625/0.03125 mm 三套离散的差值定位,
判断四位小数可靠性的实际薄弱位置(评审口径:Richardson E(h/2)≈|u_h-u_{h/2}|/3)。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else '.')
import numpy as np
import q1_solver as Q

t_ev = [100, 300, 600, 900, 1200, 1500, 1800]

old = (Q.dr, Q.N)
# 0.125 mm 粗离散（Δt=0.125 s）
Q.dr, Q.N = 1.25e-4, 161
_, C125, _, _ = Q.solve(14400, 0.125, Q.T_air, Q.C_air,
                        record={tt: [0, 40, 80, 120, 160] for tt in t_ev})
# 0.0625 mm 中离散（Δt=0.0625 s，时空同时折半）
Q.dr, Q.N = 6.25e-5, 321
_, C62, _, _ = Q.solve(28800, 0.0625, Q.T_air, Q.C_air,
                       record={tt: [0, 80, 160, 240, 320] for tt in t_ev})
# 0.03125 mm 主离散（Δt=0.0625 s，与生产一致）
Q.dr, Q.N = 3.125e-5, 641
_, C31, _, _ = Q.solve(57600, 0.0625, Q.T_air, Q.C_air,
                       record={tt: [0, 160, 320, 480, 640] for tt in t_ev})
Q.dr, Q.N = old

d1 = np.abs(C62 - C31)   # 0.0625 -> 0.03125
d2 = np.abs(C125 - C62)  # 0.125 -> 0.0625
cols = ['0', '0.5', '1', '1.5', '2']
print('表2 各单元格 |C(0.0625mm)-C(0.03125mm)| (行=时间s, 列=0/0.5/1/1.5/2cm):')
for tt, row in zip(t_ev, d1):
    print('%5d | ' % tt + '  '.join('%.2e' % v for v in row))
i, j = np.unravel_index(d1.argmax(), d1.shape)
print('最大差 %.2e 位于 t=%ds, 列=%scm；Richardson 误差(0.03125mm) ≈ %.2e'
      % (d1.max(), t_ev[i], cols[j], d1.max() / 3))
print()
print('0.125mm->0.0625mm 最大差 %.2e (Richardson 误差(0.0625mm) ≈ %.2e)'
      % (d2.max(), d2.max() / 3))
print()
print('表2 末行(1800s) 0.125/0.0625/0.03125mm:', np.round(C125[-1], 4),
      np.round(C62[-1], 4), np.round(C31[-1], 4))
print('表2 首行(100s)  0.125/0.0625/0.03125mm:', np.round(C125[0], 4),
      np.round(C62[0], 4), np.round(C31[0], 4))
