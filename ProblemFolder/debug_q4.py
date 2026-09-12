# -*- coding: utf-8 -*-
"""q4 短时演化探针:打印前若干步的表面/中心值,检查表面行系数。
仅调试用——所有数值逻辑直接复用 q4_solver.py,不做重复实现
(历史上本文件曾因重复实现表面通量系数多乘 1/DZ 而出错,故改为 import 复用)。
"""
import sys, os
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from q4_solver import (NZ, zeta, T0, C0, R0, step, RofT, C_air, D4)

T, C = np.full(NZ, T0), np.full(NZ, C0)
print('初始: D4(C0,300K)=%.3e  R(0)=%.4f  C_air(0)=%.5f' % (D4(C0, 28 + 273.15), RofT(0), C_air(0)))
for n in range(1, 13):
    T, C = step(T, C, (n - 1) * 5.0, n * 5.0, 5.0)
    if n in (1, 2, 4, 8, 12):
        print('t=%3ds 中心T=%.4f 表面T=%.4f | 中心C=%.5f 表面C=%.5f'
              % (n * 5, T[0], T[-1], C[0], C[-1]))
