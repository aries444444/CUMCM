# -*- coding: utf-8 -*-
"""由 ProblemFolder 下的源程序自动生成论文附录：ProblemFolder/附录源程序.tex

用法: python build_appendix.py
说明: 附录内容与源程序逐字一致，改动代码后重新运行即可同步附录，
      避免手工誊抄漏改，也避免在附录中重复粘贴同一段公共代码。
"""
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, 'ProblemFolder')
OUT = os.path.join(SRC, '附录源程序.tex')

# (附录编号, 标题, 文件名, 功能说明)
FILES = [
    ('附录 1', '公共数值内核', 'common.py',
     '附件定位与读入、三对角求解（纯 Python 追赶法与 LAPACK 带状求解两种等价实现）、'
     '守恒型有限体积 Crank--Nicolson 单步（固定半径 \\texttt{crank\\_radial} 与随体坐标 '
     '\\texttt{crank\\_zeta} 两套）、交错 Picard 耦合步、烘干时长判据事件求根。'
     '四个求解器共用，全文只实现一次。'),
    ('附录 2', '问题 1 求解器：预热平衡阶段', 'q1_solver.py',
     '常物性 + 非线性 $D(C)$，$\\Delta r=0.03125$ mm、$\\Delta t=0.125$ s；'
     '输出表 1、表 2 与 \\texttt{result1.xlsx}，并完成 Bessel 解析锚点、空间收敛、'
     '离散质量守恒、渐近行为四项自检。'),
    ('附录 3', '问题 2 求解器：前 3 h 热质耦合', 'q2_solver.py',
     '附录 3 变物性、热质双向耦合，$\\Delta r=0.125$ mm、$\\Delta t=0.25$ s；'
     '输出表 3、表 4 与 \\texttt{result2.xlsx}（每 1 s、每 0.1 cm），'
     '并做时间收敛、空间收敛、渐近行为检验。'),
    ('附录 4', '问题 3 求解器：全程烘干与时长', 'q3_solver.py',
     '附录 3 变物性、固定半径，$\\Delta r=0.0625$ mm、$\\Delta t=2.5$ s，'
     '以 $\\max C\\le 0.15$ kg/kg 为判据推进至烘干结束；输出 $t^*$、表 5 与 '
     '\\texttt{result3.xlsx}（每 60 s），并做离散质量守恒、时间收敛、'
     '$t^*$ 五点空间收敛链检验。'),
    ('附录 5', '问题 4 求解器：含收缩动边界', 'q4_solver.py',
     '附录 4 变物性，附件 2 半径 $R(t)$ 经 Pchip 保形插值构造动边界，'
     '随体坐标 $\\zeta=r/R(t)$ 下 $\\Delta\\zeta=1/160$、$\\Delta t=5$ s；'
     '输出 $t^*$、表 6 与 \\texttt{result4.xlsx}（每 60 s），'
     '并做动边界质量守恒、时间收敛、$\\zeta$ 五点收敛链检验。'),
    ('附录 6', '物性$\\times$几何消融与参数灵敏度', 'q78_sensitivity.py',
     '表 7 的「附录 3/附录 4 物性 $\\times$ 固定/收缩几何」四组合消融，'
     '与表 8 的 $T_a$、$C_a$、$h$、$h_m$ 各 $\\pm$ 变体下 $t^*$ 重算。'),
    ('附录 7', '验证程序一：表 2 逐格离散差定位', 'verify_q1.py',
     '在 0.125 / 0.0625 / 0.03125 mm 三套离散上逐格比较问题 1 的表 2，'
     '按 Richardson 估计 $E(h/2)\\approx|u_h-u_{h/2}|/3$ 定位四位小数的薄弱格点。'),
    ('附录 8', '验证程序二：问题 3 耦合残差与收敛链', 'verify_q3.py',
     '量化 Picard 两轮耦合残差与单步变化之比；给出 $t^*$ 在 '
     '0.5/0.25/0.125 mm 上的收敛链与 Richardson 外推；重算参数灵敏度。'),
    ('附录 9', '验证程序三：问题 4 复算、消融与灵敏度', 'verify_q4.py',
     '复用求解器完成 6 h 行离散收敛、时间收敛、Picard 残差量化、'
     '$t^*$ 五点收敛链、$\\zeta$ 离散上的守恒闭合，以及消融与灵敏度复核。'),
]

INTRO = r"""
\noindent 本附录给出本文建模与求解用到的\textbf{全部完整、可运行的源程序}。
四个问题共用同一套守恒型有限体积离散与 Crank--Nicolson 时间推进，
其附件读取、三对角求解、离散单步、交错 Picard 耦合步与烘干时长事件求根等公共部分
\textbf{只实现一次}，集中于公共内核 \texttt{common.py}（附录 1）；
各求解器只保留自己的物性经验式、烘房边界数据、网格步长与输出格式。
因此本文任何算法都只有一个实现，附录中\textbf{不存在成段重复的代码}：
同一处改动即可同步全部问题。

\noindent 运行环境：Python 3 + NumPy + SciPy + openpyxl。
运行前将赛题附件 \texttt{附件1.xlsx}、\texttt{附件2.xlsx} 与脚本放在同一目录
（\texttt{common.annex} 亦会自动检索仓库目录结构）。
各程序以 \texttt{python 脚本名.py} 直接运行，结果写入 \texttt{result1--4.xlsx}
并在控制台打印论文中的各张表格与检验量。
"""


def tex(s):
    """转义文件名中的下划线（代码块内不转义）。"""
    return s.replace('_', r'\_')


def main():
    lines = ['% 本文件由 build_appendix.py 自动生成：请勿手工编辑，改动代码后重新生成。',
             '',
             r'\section*{附录：源程序代码}',
             INTRO]
    # 文件清单表
    lines += [r'\begin{table}[htbp]', r'\centering', r'\zihao{-5}',
              r'\setlength{\tabcolsep}{4pt}', r'\renewcommand{\arraystretch}{1.25}',
              r'\begin{tabular}{@{}lll@{}}', r'\toprule',
              r'附录 & 源程序 & 对应内容 \\', r'\midrule']
    for no, title, fn, _ in FILES:
        lines.append('%s & \\texttt{%s} & %s \\\\' % (no, tex(fn), title))
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}', '']

    total = 0
    for no, title, fn, desc in FILES:
        path = os.path.join(SRC, fn)
        with open(path, encoding='utf-8') as fh:
            code = fh.read().rstrip('\n').split('\n')
        total += len(code)
        lines += [r'\subsection*{%s\quad %s（\texttt{%s}）}' % (no, title, tex(fn)),
                  '',
                  r'\noindent %s' % desc,
                  '',
                  r'\begin{lstlisting}']
        lines += code
        lines += [r'\end{lstlisting}', '']
        print('%-22s %4d 行' % (fn, len(code)))

    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('合计 %d 行代码 -> %s' % (total, OUT))


if __name__ == '__main__':
    main()
