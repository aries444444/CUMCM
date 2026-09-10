# 2026 CUMCM A题 参考文献清单（教材 + 论文，均经摘要核实）

> 检索时间：2026-09-10。来源：Web 检索（Google/Bing 索引）、Semantic Scholar、ScienceDirect/RSC/MDPI/CNKI 条目页的摘要级内容。
> 标记说明：**[摘要已核]** = 已读到该文摘要原文/摘要级内容并确认与本题相关；**用途** = 建议用在论文哪一部分。
> 注：本机网络策略禁止直连部分出版社网站（pubs.rsc.org、sciencedirect 均被拦截），以下摘要核实基于搜索引擎返回的摘要文本；引用前建议用学校图书馆数据库（知网/Web of Science/Elsevier/RSC）再核对一遍原文。

---

## 一、教材与专著（打基础用，按阅读顺序）

### 1. Crank, J. *The Mathematics of Diffusion*, 2nd ed., Clarendon Press / Oxford University Press, 1975（ISBN 9780198534112）★ 本题圣经
- 覆盖：扩散方程、圆柱中的扩散（第 5 章，含 Bessel 级数解）、**移动边界**（第 7 章）、有限差分法（第 8 章）、变扩散系数、**热与水分同时扩散**（末章）。
- 用途：问题 1~4 的方程、解析验证锚点（冻结 D 时的级数解）、移动边界处理，全部能在本书找到对应章节。有中译本《扩散的数学》。
- 姊妹篇：Carslaw, H.S. & Jaeger, J.C. *Conduction of Heat in Solids*, 2nd ed., Oxford, 1959 —— 热传导侧的同款工具书，两书互相引用，合起来是本题的标准理论底座。

### 2. Incropera, F.P., DeWitt, D.P., Bergman, T.L., Lavine, A.S. *Fundamentals of Heat and Mass Transfer*, 6th ed. (Wiley, 2006) / 8th ed. (Wiley, 2017)
- 第 5 章 Transient Conduction：5.4~5.6 **圆柱径向瞬态导热与对流边界**（5.6.1 精确解、5.6.2 近似解）、5.10 有限差分（显式/隐式）；附录 D Heisler 图。
- 第 14 章 Diffusion Mass Transfer：Fick 定律、14.7 **瞬态扩散**（与第 5 章完全平行的类比）。
- 用途：模块 A/B 的标准教材，第 5 章 5.6 节可直接对着写问题 1 的线性部分。

### 3. Whitaker, S. "Simultaneous Heat, Mass, and Momentum Transfer in Porous Media: A Theory of Drying", *Advances in Heat Transfer*, Vol. 13, pp. 119–203, 1977（DOI 10.1016/S0065-2717(08)70223-5）★ 干燥理论奠基文献
- 体积平均法建立多孔介质干燥的输运方程组；已被引 1183+ 次。1998 年 Whitaker 有续篇 "Coupled Transport in Multiphase Systems: A Theory of Drying"（*Advances in Heat Transfer*, pp. 1–102）。
- 用途：论文"问题重述/模型假设"部分引用，说明为什么可以用连续介质扩散-导热框架描述药材（多孔介质）；也是"为什么题给经验式就够用"的理论台阶。

### 4. Mujumdar, A.S. (ed.) *Handbook of Industrial Drying*, 4th ed., CRC Press, 2015（ISBN 9781466596658）
- 相关章节：Ch.1 Principles/Classification/Selection of Dryers（干燥速率曲线、含水率概念）、"Basic Process Calculations and Simulations in Drying"、"Transport Properties in the Drying of Solids"。
- 用途：干燥三阶段、干基/湿基、临界/平衡含水率的权威出处；引言部分引用。

### 5. 杨世铭, 陶文铨. 《传热学》（第 4/5 版）, 高等教育出版社
- 第 2~3 章稳态/非稳态导热（集总参数法、毕渥数、无限长圆柱）、第 4 章数值解法。
- 用途：中文快速入门传热模块 A（比 Incropera 薄、中文好读）。

### 6. 化工原理教材（任选一本，看"干燥"一章）：陈敏恒《化工原理》下册 / 柴诚敬《化工原理》
- 干燥章：湿空气性质、干燥速率曲线（预热/恒速/降速）、临界与平衡含水率、物料衡算。
- 用途：模块 C 全部概念；中文教材例子最贴合本题背景。

---

## 二、直接相关论文（几何与模型和本题几乎一致）★ 精读

### P1. Simal, S., Rosselló, C., Berna, A., Mulet, A. "Drying of shrinking cylinder-shaped bodies", *Journal of Food Engineering*, 37(4): 423–435, 1998（DOI 10.1016/S0260-8774(98)00095-8）**[摘要已核]**
- 摘要要点：圆柱固体干燥过程中温度、平均水分与水分分布的预测模型；**考虑收缩**；宏观热平衡 + 菲克定律微观质平衡同时求解（Runge-Kutta-Merson + 有限差分）；**有效扩散系数写成温度与局部含水率的函数**；用 90°C 西兰花茎干燥曲线标定，解释方差 99.8%。
- 用途：问题 4（收缩）的**模板论文**，其"热平衡+扩散+移动边界+D(T,C)"结构就是本题问题 2~4 的骨架；模型假设写法（如"导热系数视为无限大→温度均匀"）值得对比讨论。

### P2. Abbasi Souraki, B., Mowla, D. "Axial and radial moisture diffusivity in cylindrical fresh green beans in a fluidized bed dryer with energy carrier: Modeling with and without shrinkage", *Journal of Food Engineering*, 88(1): 9–19, 2008（DOI 10.1016/j.jfoodeng.2007.05.013）**[摘要已核]**
- 摘要要点：圆柱绿豆用**无限/有限圆柱菲克模型**求轴向与径向有效扩散系数；**含收缩与不含收缩**两种建模对比；收缩时扩散系数估计值更小；轴向 D 显著大于径向 D；D 与温度呈 **Arrhenius** 关系。
- 用途：①问题 1 圆柱几何"1D 径向模型"合理性的佐证与反例讨论（长径比大→径向主导；轴向不可忽略时怎么办）；②"忽略收缩会高估 D"可直接引用到问题 4 动机。

### P3. "Heat and moisture transport in okra cylinders with shrinkage effects under solar drying: a multiphysics-based simulation approach", *Sustainable Food Technology* (RSC), 2025（DOI 10.1039/D4FB00343H）**[摘要已核]**
- 要点：秋葵**圆柱**（半径 6 mm）干燥的**热-湿耦合**瞬态模拟；2D 轴对称；COMSOL 耦合"固体传热"+"稀物质传递"；**ALE 动网格处理收缩**；表面对流干燥边界（含蒸发项）。
- 用途：最现代化的同类算例（2025 年），问题 4 的 ALE/移动网格实现细节、边界条件写法都能对照；也可作为"多物理场软件 vs 自编 FDM"的方法对比引用。

### P4. 任迪峰.《中药材干燥过程中质量退化及优化干燥工艺的研究》（博士学位论文，中国农业大学）
- 要点：第四章把常见根类中药材假设为**近似变径轴对称圆柱体**，建立温度传导与水分扩散**耦合的二维传热传质有限元模型**（轴对称圆柱坐标微元体热量/质量平衡）。
- 用途：**中文文献里与本题重合度最高的一份**——同样是"根类药材 + 轴对称圆柱 + 热质耦合 FEM"；其文献综述还梳理了 Chau、Irudayaraj、Liu、Sakai、Hayakawa 等一维/二维热质耦合建模的谱系，可直接当"文献综述"底稿。

### P5. Kim, M.H., Kim, C.S., Park, S.J., Lee, C.H. "Numerical Analysis of Moisture Diffusion for Hot-air Drying of Ginseng"（人参热风干燥水分扩散的数值解析研究），《韩国农业机械学会志》, 1998
- 要点：人参视为**圆柱体**，水分仅径向迁移，控制方程 **∂M/∂t = (1/r)∂/∂r[D_eff(T,M)·r·∂M/∂r]**，数值求解抛物线 PDE 模拟热风干燥。
- 用途：问题 1 水分方程（含 D(C)）的**一字不差的先例**；1999 年同团队在 *Food Engineering Progress* 还有缓苏（tempering）干燥的隐式有限差分续作，可作为时间步长/隐式格式选择参照。

### P6. 庞凌云, 袁志华, 詹丽娟, 等. 《铁棍山药片热风干燥过程中传热传质规律》,《食品与机械》40(8), 2024（DOI 10.13652/j.spjx.1003.5788.2023.60171）
- 要点：山药片按**圆柱体**建模（半径 ~20 mm），测定导热系数/比热容，ANSYS Transient Thermal 模拟温度场；Modified Page 模型拟合传质；中心温度最低、外表面最高，内外温差由 7.43 K 降至 5.82 K。
- 用途：①中文写作范例（结构、术语）；②"圆柱药材热风干燥中心-表面温差"的实测佐证；③问题 1 表 1 温度形态（外高内低、梯度渐小）的实验对应。

### P7. 王学成, 康超超, 伍振峰, 等. 《二至丸热风干燥过程温度均匀性模拟与实验》,《中草药》51(5): 1226–1232, 2020
- 要点：COMSOL 传热传质模型 + Fick 第二定律平板模型求 Deff（0.76×10⁻⁷~3.94×10⁻⁷ m²/s），60/80/100°C 热风，探针测温验证。
- 用途：中药材 + COMSOL + Fick 的中文范式；Deff 量级可与题给公式互证。

### P8. Pacheco-Aguirre, F.M., Ladrón-González, A., Ruiz-Espinosa, H., García-Alvarado, M.A., Ruiz-López, I.I. "A method to estimate anisotropic diffusion coefficients for cylindrical solids: Application to the drying of carrot", *Journal of Food Engineering*, 125: 24–33, 2014（DOI 10.1016/j.jfoodeng.2013.10.015）**[摘要已核]**
- 摘要要点：基于**有限各向异性圆柱**非稳态传质方程**解析解**（含 Bessel 函数）估计轴向/径向/角向扩散系数；胡萝卜 2.2 cm 直径圆柱，80°C；D = 0.53~2.93×10⁻⁹ m²/s；三方向差异显著。
- 用途：①圆柱扩散解析解的现成公式库（验证你的 FDM 代码）；②"各向异性是否必要"的讨论素材（本题 25 cm 长、径 2 cm，长径比 12.5，可用它论证径向 1D 近似的合理性）。

---

## 三、方法与理论支撑论文

### P9. Chandra Mohan, V.P., Talukdar, P. "Three dimensional numerical modeling of simultaneous heat and moisture transfer in a moist object subjected to convective drying", *International Journal of Heat and Mass Transfer*, 53(21–22): 4638–4650, 2010（DOI 10.1016/j.ijheatmasstransfer.2010.06.029）**[摘要已核]**
- 要点：耦合热湿方程 3D 数值解（有限体积法、全隐式）；h 由自编 CFD 给出，h_m 由**热-浓度边界层类比**（Chilton–Colburn）换算；气流 313→353 K 省时 40%。
- 用途：数值方法（隐式、耦合求解）写法；"变表面换热系数 vs 常数系数"的讨论——支持本题"常数 h/h_m"简化假设的合理性论证（原结论：是否必要取决于物料干燥行为）。

### P10. Defraeye, T., Blocken, B., Carmeliet, J. "Analysis of convective heat and mass transfer coefficients for convective drying of a porous flat plate by conjugate modelling", *International Journal of Heat and Mass Transfer*, 55(1–3): 112–124, 2012（DOI 10.1016/j.ijheatmasstransfer.2011.08.047）**[摘要已核]**
- 要点：共轭建模（CFD 空气域+多孔材料域）反推对流系数，发现系数在恒速段/降速段近似常数、**在段间过渡处才显著变化**；热质类比在两段内成立、过渡处偏差大。
- 用途：模型假设处引用"常数 h、h_m 在常规对流干燥中可接受"，堵住评阅人"为什么不考虑系数时空变化"的疑问。

### P11. García del Valle, J., Sierra Pallares, J. "Analytical solution for the coupled heat and mass transfer formulation of one-dimensional drying kinetics", *Journal of Food Engineering*, 230: 99–113, 2018（DOI 10.1016/j.jfoodeng.2018.02.029）**[摘要已核]**
- 要点：耦合热质守恒的**解析解**（内部解耦、边界耦合：内部毛细扩散 + 界面蒸发）；自然涌现**两个干燥阶段**（对流控制段→扩散控制段）；假设恒定 Deff、无收缩。
- 用途：①理解本题"为什么边界条件是耦合的、内部方程是解耦/弱耦合的"；②解析解可作数值代码验证锚点；③论文"模型合理性"部分的引用。

### P12. García-Alvarado, M.A., Pacheco-Aguirre, F.M., Ruiz-López, I.I. "Analytical solution of simultaneous heat and mass transfer equations during food drying", *Journal of Food Engineering*, 142: 39–45, 2014
- 要点：Laplace 变换求食品干燥热质同时传递解析解，考虑温度对界面水分的影响。
- 用途：与 P11 同族，问题 2 两阶段解析讨论备用。

---

## 四、中药材干燥数据与综述（参数标定、量级互证、引言背景）

### P13. 《中草药干燥加工现状及发展趋势》,《南京中医药大学学报》, 2021（DOI 10.14148/j.issn.1672-0482.2021.0786）——中文综述
- 用途：引言"中药材干燥背景与难点"整段可直接取材；经验/半理论/理论模型分类（Weibull、Page、Midilli 等）适合在论文里交代"为什么不直接用经验模型、而要建理论模型"。

### P14. "Drying kinetics of Salviae Miltiorrhizae Radix et Rhizoma (丹参) and dynamics of active components in drying process"（PubMed PMID 39929654, 2025）
- 要点：丹参热风干燥 **Deff = 1.746×10⁻⁸ ~ 9.354×10⁻⁸ m²/s**，活化能 **40.0 kJ/mol**（红外 26.62）。
- 用途：题给公式 D=2.4×10⁻³ e^(−0.45C) e^(−3850/T) 的 Ea≈32 kJ/mol、D~10⁻⁹~10⁻⁸ m²/s 与实测完全同量级 → **论文里用这类数据给题给经验式做"物理合理性背书"**（高分讨论点）。

### P15. 《鲜地黄干燥动力学及其对环烯醚萜类成分影响研究》（《食品研究与开发》等, 2024）
- 要点：热风 40–60°C，Midilli 模型最优；Deff = 6.681×10⁻¹⁰ ~ 8.969×10⁻⁹ m²/s，Ea = 34.94 kJ/mol。
- 用途：同上，量级互证；"变温干燥（60°C 烘 42h 后 40°C）"可作问题 2/3 分段策略的工艺背景。

### P16. "Kinetics and variation of volatile components of Atractylodis Macrocephalae Rhizoma (白术) during hot-air drying",《中国中药杂志》
- 要点：30–70°C 热风，Midilli 模型；Deff = 1.04×10⁻⁹ ~ 6.28×10⁻⁹ m²/s，Ea = 37.47 kJ/mol。
- 用途：量级互证 + 根茎类药材（与本题圆柱形药材同型）干燥行为参考。

### P17. 李坤.《根茎类中药材热风干燥的失水动力学模型及特性机理研究》（硕士论文，云南师范大学, 2020）
- 要点：根茎类药材（占植物类中药材 80% 以上）热风干燥传热传质理论；基于质量/热量守恒建立失水动力学模型；计算 Deff 与活化能；分析温度、风速、相对湿度、厚度的影响。
- 用途：中文"理论模型"写作范例 + 根茎类药材参数影响规律。

---

## 五、检索渠道与检索式（后续自己查补）

- **Google Scholar**（scholar.google.com）：检索式示例
  - `"drying" AND "cylinder" AND ("shrinkage" OR "moving boundary")`
  - `"coupled heat and mass transfer" drying numerical simulation`
  - `"effective moisture diffusivity" cylinder Arrhenius`
  - 中文：`中药材 热风干燥 动力学 有效水分扩散系数`；`圆柱 干燥 数值模拟 传热传质`
- **Semantic Scholar**（semanticscholar.org）：按 P1、P2 的"引用/被引"顺藤摸瓜（P1 被引 83+，P9 被引 130+，Whitaker 1977 被引 1183+）。
- **知网 CNKI**：中文论文与学位论文（P4、P6、P7、P13~P17）；检索式 `中药材 AND 干燥 AND (数值模拟 OR 传热传质 OR 动力学)`。
- **Web of Science / Scopus**：英文权威回溯，按 DOI 定位。
- 获取全文：学校图书馆数据库（Elsevier/RSC/Wiley/知网）均可下载上述文献；部分老文献（P1/P2/P5）也在 ResearchGate 有作者上传版。

---

## 六、文献 ↔ 论文章节映射（写作速查）

| 论文部分 | 建议引用 |
|---|---|
| 问题重述/背景 | P13 综述、Mujumdar 手册、P4 学位论文综述段 |
| 模型假设（连续介质、常数 h/h_m、径向 1D） | Whitaker 1977、P10 Defraeye、P8（长径比论证）、P2（轴向 vs 径向） |
| 控制方程与边界条件 | Crank、Incropera、P5（人参圆柱方程）、P1（热平衡+Fick 结构） |
| 参数来源与合理性 | 附录经验式 + P14~P16 实测 Deff/Ea 量级互证 |
| 数值方法 | Incropera 5.10、Crank 第 8 章、P9（隐式耦合求解） |
| 模型验证 | P11/P12 解析解、Crank 圆柱级数解、P6（中心-表面温差形态） |
| 问题 4 收缩 | P1（模板）、P3（ALE）、P2（收缩对 D 的影响） |
| 灵敏度/讨论 | P9（变系数）、P10（过渡段）、P6/P7（温度均匀性） |
