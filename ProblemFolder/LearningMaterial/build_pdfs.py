# -*- coding: utf-8 -*-
"""Rebuild both PDFs from the Chinese-named .tex files and place them in ProblemFolder/A题."""
import os, shutil, subprocess

root = r'f:/Desktop/CUMCM/正赛'
dest = os.path.join(root, 'ProblemFolder')  # 文档与 PDF 产物放在 ProblemFolder 根下（A题 子目录只放题目与附件）

# locate source tex files
guide_src = os.path.join(dest, 'A题物理化学补习指南.tex')
lit_src = os.path.join(root, 'A题文献清单.tex')  # lit tex is at root
if not os.path.exists(lit_src):
    lit_src = os.path.join(dest, 'A题文献清单.tex')

jobs = [('guide', guide_src, 'A题物理化学补习指南.pdf'),
        ('lit', lit_src, 'A题文献清单.pdf')]

for tag, src, final_name in jobs:
    tmp = os.path.join(root, tag + '.tex')
    shutil.copy(src, tmp)
    for passn in (1, 2):
        r = subprocess.run(['xelatex', '-interaction=nonstopmode', '-halt-on-error', tmp],
                           cwd=root, capture_output=True, text=True, encoding='utf-8', errors='replace')
        print(f'[{tag}] pass {passn} exit={r.returncode}')
        if r.returncode != 0:
            tail = r.stdout.splitlines()[-15:]
            print('\n'.join(tail))
            break
    pdf = os.path.join(root, tag + '.pdf')
    if os.path.exists(pdf):
        shutil.copy(pdf, os.path.join(dest, final_name))
        print(f'[{tag}] OK -> {final_name} ({os.path.getsize(pdf)} bytes)')

# tidy: remove aux/log/out files at root
for f in os.listdir(root):
    if f.endswith(('.aux', '.log', '.out', '.toc')) or f in ('texput.log',):
        try:
            os.remove(os.path.join(root, f))
        except OSError:
            pass
print('done')
