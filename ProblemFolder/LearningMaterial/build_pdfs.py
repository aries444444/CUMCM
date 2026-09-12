# -*- coding: utf-8 -*-
"""Rebuild all PDFs from the .tex files in this directory (LearningMaterial)."""
import os, shutil, subprocess

here = os.path.dirname(os.path.abspath(__file__))  # LearningMaterial
root = os.path.dirname(here)                        # repo root (CUMCM)
dest = here                                         # PDFs go next to the sources

jobs = [('guide', 'A题物理化学补习指南.tex', 'A题物理化学补习指南.pdf'),
        ('lit', 'A题文献清单.tex', 'A题文献清单.pdf'),
        ('q1lecture', 'A题问题1知识点串讲.tex', 'A题问题1知识点串讲.pdf'),
        ('graddiv', 'A题梯度与散度讲解.tex', 'A题梯度与散度讲解.pdf'),
        ('constraints', 'A题格式与内容硬约束.tex', 'A题格式与内容硬约束.pdf')]

for tag, name, final_name in jobs:
    src = os.path.join(here, name)
    if not os.path.exists(src):
        print(f'[{tag}] skip: {name} not found')
        continue
    tmp = os.path.join(root, tag + '.tex')
    shutil.copy(src, tmp)
    for passn in (1, 2, 3):
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

# tidy: remove temp tex/pdf and aux/log/out/toc at root
for f in os.listdir(root):
    if f.endswith(('.aux', '.log', '.out', '.toc')) or f in ('texput.log',):
        try:
            os.remove(os.path.join(root, f))
        except OSError:
            pass
for f in os.listdir(root):
    if f.endswith('.tex') and f[:-4] in [t for t, _, _ in jobs]:
        os.remove(os.path.join(root, f))
    if f.endswith('.pdf') and f[:-4] in [t for t, _, _ in jobs]:
        os.remove(os.path.join(root, f))
print('done')
