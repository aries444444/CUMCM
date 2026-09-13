# -*- coding: utf-8 -*-
"""从 支撑材料/ 的八个旧版程序生成论文附录的完整代码清单（附录源程序_旧版.tex）。

每个程序一个 subsection + 完整 lstlisting；内容与支撑材料中的 .py 文件逐字节一致。
运行：python build_appendix_old.py
"""
import os

BASE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(BASE, '支撑材料')
OUT = os.path.join(BASE, '附录源程序_旧版.tex')

FILES = ['q1_solver.py', 'q2_solver.py', 'q3_solver.py', 'q4_solver.py',
         '表7消融_表8灵敏度.py', 'verify_q1.py', 'verify_q3.py', 'verify_q4.py']

parts = []
for name in FILES:
    path = os.path.join(SUPP, name)
    with open(path, encoding='utf-8') as f:
        code = f.read().rstrip('\n')
    esc = name.replace('_', r'\_')
    parts.append(r'\subsection*{%s}' % esc)
    parts.append(r'\begin{lstlisting}[language=Python]')
    parts.append(code)
    parts.append(r'\end{lstlisting}')
    parts.append('')

with open(OUT, 'w', encoding='utf-8') as f:
    f.write('\n'.join(parts))

print('已生成 %s（%d 行，%d 个程序）' % (OUT, len('\n'.join(parts).splitlines()), len(FILES)))
