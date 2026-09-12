# 科研可视化任务总结

## 完成内容

1. **安装科研可视化 skill**：`alterlab-scientific-viz`（出版级 matplotlib/seaborn 规范：Okabe-Ito 色盲友好调色板、期刊尺寸、300dpi/矢量导出）。已通过安全审计（无网络调用、无凭据访问、纯绘图工具）。
2. **生成 12 张出版级静态图**（`ProblemFolder/figures_pub/`，每张含 300dpi PNG + 矢量 PDF），数据全部来自 result1-4.xlsx、附件1/2.xlsx 与论文表 7：

| 图 | 内容 |
|---|---|
| fig01 | 附件1 烘房温度/水分浓度边界（双轴） |
| fig02 | 附件2 半径收缩曲线（累计收缩约 40%） |
| fig03 | Q1 表面/中心温度与水分演化（a/b 双面板） |
| fig04 | Q1 温度场时空热图 |
| fig05 | Q2 各位置温度演化 |
| fig06 | Q2 温度场时空热图 |
| fig07 | Q3 含水率演化 + 0.15 判据线 + t\*=57.09 h |
| fig08 | Q3 含水率场时空热图（gamma 增强） |
| fig09 | Q4 含收缩含水率演化 + t\*=50.77 h |
| fig10 | Q4 含水率时空热图（NaN 掩膜呈现代理收缩动边界） |
| fig11 | Q3 vs Q4 中心含水率对比（收缩缩短烘干约 11.1%） |
| fig12 | 表7 参数灵敏度 tornado 图 |

3. **可复现脚本**：`ProblemFolder/figures_pub/make_figures.py`，求解器重跑后可一键再生成。

## 关键核对

- Q1 表面 1800 s = 36.79 °C、Q3 t\*=57.09 h、Q4 t\*=50.77 h 均与论文表 1/5/6 一致。

## 后续建议

- 论文 `A题论文.tex` 目前无插图，可直接 `\includegraphics` 引用 figures_pub 下的 PDF（矢量）。
- smart-charts 已知待办（Q3/Q4 两张 HTML 标题用旧 t\*）仍待更新。
