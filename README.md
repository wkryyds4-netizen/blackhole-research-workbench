# Blackhole Research Workbench · 黑洞研究工作台

通用黑洞场景可视化、公开数值相对论波形对照与群体推断稳健性实验。原创代码采用 [MIT](LICENSE)，第三方代码和科学数据按各自许可分开管理，详见 [NOTICE](NOTICE.md)。

## 研究论文与复现材料

2026-10-27 检查点：研究按用户要求暂停，现有成果随本次 Git 保存。新增[物理形成渠道接口与诊断](research/population-channels/README.md)（[English](research/population-channels/README.en.md)），尚未拟合观测渠道比例；已完成内容与恢复位置见[暂停记录](research/CHECKPOINT_2026-10-27.md)。本轮新增结果尚未并入 v0.3 论文或静态工作台。

2026-10-27追加成果：[质量—红移54模型研究](research/population-mass-redshift/README.md)（[English](research/population-mass-redshift/README.en.md)）。全69来源、216组条件后验经独立复算；幅度族模型差在所测网格中保持10.36–17.14个百分点，质量/红移选择自身可使中位数跨约3.95个百分点。可分离模型、缺失支持及实际选择协议的限制仍保留，尚未并入v0.3或静态工作台。

2026-10-27新增研究：[完整O3混合注入对照](research/population-mixture/README.md)（[English](research/population-mixture/README.en.md)）。固定69事件和12个模型，更换官方注入表示后，g中位数最大移动0.13703个百分点；混合提议仍未覆盖33.305257%的目标概率质量。本轮结果尚未并入下列v0.3论文或静态工作台，详见研究路线。

阶段快照（2026-10-25，研究按用户要求暂停）：[全69来源及后续诊断](research/population-catalog/README.md)已扩展到11180578原始先验后验行，完成O3/事件/联合成对抽稀及未见贡献的条件边界。下列v0.3论文/PDF已纳入这些结果；静态工作台已初步接入同一研究快照，旧81设置视图保留在历史折叠入口。界面尚未完成全部交互与响应式验证。续研所需证据与未完成事项见[研究路线](research/RESEARCH_ROADMAP.md)。

**[中文 LaTeX 论文 PDF v0.3（14页）](output/pdf/blackhole-population-robustness-latex-zh-v0.3.pdf)** · **[English LaTeX PDF v0.3 (14 pages)](output/pdf/blackhole-population-robustness-latex-en-v0.3.pdf)**

[中文 LaTeX 源码](research/paper/v0.3/latex-zh/main.tex) · [English LaTeX source](research/paper/v0.3/latex-en/main.tex) · [双语构建与审阅](research/paper/README.md) · [English guide](research/paper/README.en.md) · [科学计算复现](research/REPRODUCING_POPULATION.md)

The full Chinese and English editions report the same conditional population diagnostics. Both include the methods, results, limitations, references, and reproduction appendix; neither is peer reviewed.

研究稿v0.3《双黑洞群体推断的后验来源、有限银行稳定性与缺失支持》汇总完整来源恢复、两个六模型控制、512次成对人工抽稀及未见贡献的条件边界，含14个公式、5张矢量图、4张表和15条文献。**未同行评审；本稿报告条件诊断，不声称测得形成渠道比例或认证实际遗漏贡献。** 本次快照包含v0.3；原始科学数据和大型模型按来源单独下载。

已公开的历史v0.2仍保留：[中文PDF](output/pdf/blackhole-population-robustness-latex-zh-v0.2.pdf)、[英文PDF](output/pdf/blackhole-population-robustness-latex-en-v0.2.pdf)、[中文Markdown](research/paper/manuscript.zh.md)、[English Markdown](research/paper/manuscript.en.md)。其中81设置、2048次经验重采样属于早期冻结研究，不能代替v0.3的全来源与成对抽稀结果。

## 快速体验：无需安装科学计算环境

在仓库根目录使用 Python 3 启动静态工作台：

```sh
python -m http.server 8765 --bind 127.0.0.1 --directory simulator/public
```

打开 **http://127.0.0.1:8765/workbench.html?view=waveform**。

- **通用场景**：单 Kerr/Schwarzschild 光学、通用轨迹接口与自旋方向；公开版含4份许可明确的SXS场景快照。
- **波形比较**：0389与NRSur7dq4复h22对照、145个配置、PSD/频段/窗函数敏感性及JSON导出。
- **对齐与动力学**：视界轨迹对齐敏感性、代理波形与旋转框架。
- **群体稳健性**：69事件、12个条件模型、10项成对对比；区分可信区间、人工抽稀与假设边界，附来源哈希、中英PDF和快照导入/导出。旧81设置视图保留供历史对照。

静态模式显示仓库中已保存的科研结果；根路径的交互波形计算需另启动 [Ripple后端](simulator/README.md)。WebGL2不可用时采用坐标视图。无需Git LFS。

## 开发与研究复现

Node.js 20或以上用于前端测试：

```sh
npm ci
npm test
```

浏览器发布检查与贡献要求见 [CONTRIBUTING](CONTRIBUTING.md)。科研计算按层独立运行：

| 内容 | 入口 |
|---|---|
| Ripple计算、GW150914处理、渲染重建 | [simulator/README.md](simulator/README.md) |
| 通用场景的单位与能力约定 | [SCENE_CONTRACT.md](simulator/SCENE_CONTRACT.md) |
| NR–代理比较方法、运行与验证 | [第三阶段记录](research/stage3-results/README.md) |
| 原始SXS波形下载与解码 | [波形来源记录](research/wave-source-results/README.md) |
| 代理模型Windows构建路径 | [代理运行记录](research/surrogate-results/README.md) |

发布包不包含Python环境、模型权重、大型原始样本、下载缓存或许可未明确的2327数据。需重算时按相应脚本下载，来源和哈希保留在研究记录中。历史检查记录可能引用未随仓库发布的本地文件，不能把这些记录当作新机器已通过测试的证明。

这些工具用于探索性研究。双中心光学不是完整双黑洞数值时空；复h22失配不是探测器网络失配；训练成员关系未确认，不能声称训练集外验证；现象学倾角混合比例不能直接等同形成渠道比例。

## 群体研究的详细记录

2026-10-25 后续研究公开：

- [阶段三：事件成员与两个PE来源试点](research/population-stage3/README.md)：逐行确认两个事件的IMRPhenomXPHM标签，并发现10/20 Hz的自旋参考频率差异。
- [阶段四：先验恢复与条件替换](research/population-stage4/README.md)：重建504,649条教学先验数值，量化两个事件的样本替换效应，定位教学子集遗漏的高权重点。
- [阶段五：原始先验样本与尾部检查](research/population-stage5/README.md)：恢复GW190805的169,040条原始先验样本；在指定低自旋检验点，使用匹配先验后，权重ESS从完整cosmo表示的30.17提高到998.36，最大单行贡献从15.86%降至1.04%。

上述历史阶段的后续全目录结果已纳入v0.3中英论文/PDF及新的工作台快照；v0.2冻结文件保持不变。更高ESS不证明PE采样已收敛；全局先验支持、共同自旋参考定义和注入执行协议仍有未完成核验。原始HDF、NPZ及第三方文献快照单独下载，仓库提供脚本、来源校验值、派生摘要和SVG。

**English update (2026-10-25):** Stages three through five add event provenance,
reconstruction of stored prior densities, and a single-event comparison of
original-prior and cosmologically reweighted posterior representations. At the
specified low-spin test point, correctly matched prior division raises the
importance-weight ESS from 30.17 to 998.36. This is a numerical representation
diagnostic, not a posterior-convergence proof or a formation-channel measurement.
The bilingual v0.2 papers and historical workbench results remain unchanged.
This snapshot also includes the separate bilingual v0.3 interim papers, the
all-69-source follow-ups and an initial population-view adapter. Research is
paused; physical identification and real-campaign coverage remain unresolved. See the
[reproduction guide](research/REPRODUCING_POPULATION.md) for the ordered workflow.

2026-10-24 科学研究阶段更新：[阶段一：共同先验与数值精度](research/population-stage1/README.md)及[阶段二：注入来源恢复与成对MC诊断](research/population-stage2/README.md)。现已逐行恢复41972条官方注入及保留规则，完成2048次共享样本的经验重采样，区分在该扰动下保持的模型差异与不稳定的低自旋压力案例。原始条件推断仍保留；这些诊断尚不构成形成渠道占比测量。

完成日期：2026-10-22。

2026-10-23 第三阶段：已加入 [NR–代理波形比较记录](research/stage3-results/README.md)，工作台入口 `http://127.0.0.1:8765/workbench.html?view=waveform`。0389 的复 h22 对照、外推阶敏感性及采样检查共 145 个配置；支持 PSD、频段、窗函数筛选和完整导出。2327 波形已解码，进动多模式比较尚待参考框架与记忆效应校准。

2026-10-23 更新：新增[通用场景与研究工作台](simulator/README.md)，入口 `http://127.0.0.1:8765/workbench.html`。第二轮核查10个相关仓库，建立带单位/来源的轨迹契约及光学能力检查；新增后验全域 MC 精度筛查，基准后验质量86.7%通过原探索性阈值。详见[泛化研究记录](research/generalization-review.md)。该诊断不改变原始条件后验，也不是形成渠道的正式测量。

新增模块（2026-10-23）：[双黑洞引力波实验室](simulator/README.md)。支持参数化 IMR 波形、同步轨道与极化示意、双模型比较，以及 GW150914 H1/L1 公开信号对照。Windows 可运行 `simulator/start.cmd` 启动。该模块与下述群体统计实验分别保存，不改变其结果。

这是一个可复现的探索性群体推断：使用真实事件的多维后验样本和真实搜索注入的预处理数据，检验一个**偏对齐倾角分量＋各向同性倾角分量**的混合比例对固定模型假设及选择效应的敏感性。它不是形成渠道的直接分类，也不是LVK官方结果的复现。

## 这次得到的结果

在基准假设下，偏对齐分量比例 f 的后验中位数为 **22.7%**，90%等尾条件可信区间为 **3.7%—49.1%**。忽略选择效应的对照得到中位数 **49.6%**，区间 **22.4%—75.9%**。两者相差约26.9个百分点。

基准模型的两分量相对探测效率为1.1604。忽略选择的曲线只用于展示错误分析造成的影响，不能直接解释成“已探测样本中的实际渠道比例”。

81个固定假设组合的中位数跨越1.3%—81.2%；其中36个组合未通过本实验中位数处的探索性蒙特卡洛精度筛查。其余45个组合仍跨越1.3%—70.4%。这个范围是**不同条件模型之间的敏感性范围**，不是边缘化后的可信区间，不能将81个模型等权混合而得到科学结论。

基准模型中，剔除GW190517_055101使中位数下降约3.5个百分点；剔除GW191109_010717使其上升约3.0个百分点。单事件影响小于一些固定假设的影响。

结论：在这里测试的模型族和数据处理流程中，混合比例尚不能被视为对假设稳健的数值。物理渠道比例的解释需要额外的形成模型映射、完整超参数边缘化、拟合优度检验和数据选择审计。

## 数据来源及样本审计

主要来源是[GWPopulation作者提供的GWTC-3教程](https://colmtalbot.github.io/gwpopulation/examples/GWTC3.html)。该教程提供由gwpopulation-pipe预处理的后验与注入文件；原作者明确提示教程事件选择并不严格，应作定性用途。这一限制保留在本项目中。

- 后验：`data/gwtc-3-samples.pkl`，73组，每组包含源参考系主星质量m₁、质量比q、两自旋幅度a₁/a₂、两自旋倾角余弦c₁/c₂、红移z和原参数估计先验密度prior。
- 注入：`data/gwtc-3-injections.pkl`，41,972条已检出注入；文件元数据记录总生成数73,280,198及分析时长2.0004556年。本实验不推断绝对并合率。
- 使用原文件的prior列做重要性重加权；没有把单事件后验误当成似然，也没有用目录中位数替代后验计算。
- 下载来源、文件大小和SHA-256保存于`data/manifest.json`。
- `prepare_data.py`仅允许指定的NumPy/Pandas数据对象重建，拒绝其他pickle全局对象，转存为可用`allow_pickle=False`读取的NPZ文件。模型之后只读取NPZ与JSON。

教程后验文件不保留事件名称。这里通过GWOSC的GWTC-1、GWTC-2.1、GWTC-3目录，对m₁、m₂、χ_eff、z的中位数做一对一最小代价匹配，并独立核验73个匹配结果的GPS时间严格递增。全部结果及匹配代价在`data/sample_audit.json`，名称属于重建的元数据，不应当作原文件携带的标识。

从73组中排除以下4组，得到69组，共476,997个后验样本：

| 排除对象 | 本实验理由 |
|---|---|
| GW190425 | 伴星质量中位数低于3 M☉ |
| GW190814 | 伴星质量中位数低于3 M☉；其伴星性质本身存在争议 |
| GW191204_110529 | 公开目录FAR为3.3/年 |
| GW200208_222617 | 公开目录FAR为4.8/年 |

其余O3对象使用公开目录FAR<1/年的筛选。O1/O2目录不提供该汇总字段，保留其confident集合中的双黑洞对象；不把缺失FAR当作零。

**尚未完成的选择审计**：教程的已检出注入文件不保留逐管线FAR、搜索阈值、观测运行标签或原始波形标签。因此本实验采用作者提供的found-injection集合，不能独立证明它与这里的事件筛选完全匹配。正式分析必须从原始注入产品重建相同阈值与观测运行的选择函数，并核对后验波形与事件版本。后处理的质量筛选也需要在选择模型中审计。

原始数据入口：

- [GWTC-3数据发布说明](https://gwosc.org/GWTC-3/)
- [GWTC-3后验发布](https://zenodo.org/records/5546663)
- [GWTC-2.1后验发布](https://zenodo.org/records/6513631)
- [教程引用的O1/O2/O3注入发布](https://zenodo.org/records/7890398)
- [GWOSC API说明](https://gwosc.org/api/v1/docs/)

## 统计模型

记θ=(m₁,q,a₁,a₂,c₁,c₂,z)，其中cᵢ=cos θᵢ。群体密度为

\[
p(\theta\mid f,\eta)=f\,p_A(\theta\mid\eta)+(1-f)\,p_I(\theta\mid\eta).
\]

A是偏对齐倾角分量，I是各向同性分量。f是所定义源总体中的内禀混合比例。两个分量共享质量、自旋幅度和红移分布，区别仅在倾角分布。

- A：两个倾角余弦独立，均服从以1为中心、截断到[-1,1]的高斯分布。σₜ取0.3、0.6、1.2。σₜ是余弦空间宽度，不是弧度。
- I：c₁、c₂独立且均匀分布于[-1,1]。
- 自旋幅度独立同分布于Beta(1.5,8.5)、Beta(2,5)或Beta(1,1)。
- 主星质量区间[3,100] M☉；三种模型为α=2幂律加10%高斯峰、α=3幂律加10%峰，以及纯α=2幂律。峰中心35 M☉、宽度5 M☉，截断并归一化。
- 质量比密度p(q|m₁)∝q，范围[3/m₁,1]，归一化后保证伴星质量至少3 M☉。
- p(z)∝(dV_c/dz)(1+z)^(λ−1)，z∈[0,2.3]，λ取0、2、4；宇宙学固定为Astropy Planck15。所有候选模型使用相同支持域。
- f先验可选Beta(1,1)、Beta(2,2)、Beta(1/2,1/2)。

基准：σₜ=0.6，幅度Beta(2,5)，质量α=2加10%峰，λ=2，f均匀先验。

这些是人为指定的现象学模型。孤立双星也可能有倾斜自旋，动力学系统也可能呈现复杂选择后的自旋分布；正χ_eff不能直接证明孤立演化，负χ_eff也不能独立确定动力学形成。

给定已观测事件数N，并对总体率采用对应的对数均匀率先验边缘化形式，群体似然为

\[
\mathcal{L}(f,\eta)\propto
\frac{\prod_i\left[fA_i+(1-f)I_i\right]}
{\left[f\xi_A+(1-f)\xi_I\right]^N},
\qquad
A_i=\frac1{S_i}\sum_s\frac{p_A(\theta_{is}\mid\eta)}{\pi_i(\theta_{is})}.
\]

Iᵢ同理，πᵢ是事件参数估计先验。单事件证据是与f无关的常数，在归一化的f后验中抵消。选择效率由已检出注入计算：

\[
\xi_k=\frac1{N_{\rm generated}}\sum_{j\in\mathrm{found}}
\frac{p_k(\theta_j\mid\eta)}{p_{\rm draw}(\theta_j)}.
\]

实现保留联合后验样本，因此m₁、q、χ_eff等参数相关性被保留。显示用的χ_eff=(a₁c₁+qa₂c₂)/(1+q)只用于散点图，模型本身使用两自旋幅度和倾角。

界面中的“选择效率比压力测试”将ξ_A乘以0.8—1.2的系数，ξ_I不变。它是人为扰动敏感性分析，不是已经测得的校准误差分布。忽略选择效应的对照直接删除似然分母；它并非一个经过正确探测建模的替代推断。

统计形式参考：[GWPopulation层次似然实现](https://github.com/ColmTalbot/gwpopulation/blob/main/gwpopulation/hyperpe.py)、[选择效应实现](https://github.com/ColmTalbot/gwpopulation/blob/main/gwpopulation/vt.py)、[LVK GWTC-3群体论文](https://arxiv.org/abs/2111.03634)。

## 稳健性与数值精度

1. 扫描3×3×3×3=81个固定超参数组合。界面展示当前条件模型与每次只改一个假设的对照，不把不同组合当作联合边缘化样本。
2. 同时计算校正与不校正选择效应的曲线；可改变f先验和相对选择效率。
3. 保持其余设置，逐个剔除事件，诊断单事件影响。该功能是影响分析，不构成重新校准过的事件选择方案。
4. 用权重的一、二阶矩计算事件ESS、注入ESS和近似log似然蒙特卡洛标准差。采用中位数处最小事件ESS≥50、注入ESS≥4N、σ_MC(logL)≤1作为探索性提示标准，**不是发表级收敛证明**。空心点和虚线表示未通过这一提示标准。可信区间本身未边缘化蒙特卡洛积分误差。
5. ESS指重要性权重的有效样本量。基准中位数处最小事件ESS≈150、注入ESS≈6435、σ_MC(logL)≈0.89。
6. 对真实后验样本与已检出注入进行40次独立有放回重采样，估计蒙特卡洛积分波动；基准中位数标准差约0.0254，即2.5个百分点。这个bootstrap没有重采样天体事件，不是总体抽样误差估计。
7. f网格从501点增加到2001点，基准5%、50%、95%分位数的最大变化小于0.00075。
8. 构造已知真值f=0.30的选择偏差模拟：选择校正恢复MAP≈0.3001；错误忽略选择得到≈0.6169。
9. 两分量完全相同时，均匀先验保持为均匀后验，符合不可辨识性要求。所有81组曲线均检查有限性和归一化。

检查结果在`results/validation.json`；逐项结果在`results/model.json`；基准与单事件剔除结果在`results/summary.json`。

## 界面与复现

交互图包括：条件后验对照、逐项假设敏感性、真实事件质量—有效自旋误差图、单事件剔除检验。支持选择效率扰动、先验切换、状态保存和320—1024像素宽度。

Python 3.14环境中的依赖版本见`requirements.txt`。在本目录运行：

```powershell
python -m pip install -r requirements.txt
python download_data.py
python prepare_data.py
python model.py
python validate_model.py
python build_visualization.py
```

`download_data.py`仅下载公开数据，有缓存时不重复下载。`model.py`与验证程序均可离线运行。构建交互图时，`build_visualization.py`默认写入本次聊天的授权可视化目录，可通过`--output`指定其他目录。图表使用固定版本D3 CDN；图中数据已内嵌，不运行网络数据请求。

## 进入正式研究前的下一步

- 从原始LVK文件恢复事件名称、版本、波形和逐观测运行的注入筛选；验证事件样本与选择函数完全一致。
- 对自旋幅度、倾角、质量、红移超参数联合边缘化，计算更完整的参数退化；考虑质量—自旋、质量比—自旋相关性及多次并合等替代模型。
- 在全部相关后验区域检查积分精度，增加关键事件与注入的有效样本量，研究波形及噪声系统误差。
- 做经过选择函数的后验预测检验，确认所用模型能描述数据，再讨论形成渠道映射。
- 将恒星演化和星团动力学模拟给出的可观测量分布引入模型；不能把本实验的倾角混合比例直接改名为孤立双星或动力学形成比例。

本项目没有运行原始应变的单事件参数估计，没有新增观测，没有完成全部正式研究步骤；它完成的是可交互、可复算且明确标注限制的稳健性试验。
