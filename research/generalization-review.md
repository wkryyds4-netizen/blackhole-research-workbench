# 通用场景引擎与科研工作台：第二轮项目调研和实现

日期：2026-10-23。用户选择两条线并重：通用场景引擎＋科研分析工作台。

本轮新增/重新核查 10 个仓库，并与已运行的 algometrix/blackhole-sim 组成 11 项能力清单。仓库默认分支的固定提交、目录、README、许可证与依赖说明保存于 `generalization-snapshots/`。这些候选项目没有全部安装或运行；不要把文档核查写成独立精度验证。

## 结论与优先级

采用可更换的数据适配器和明确的能力契约，而不是将某个浏览器模拟器扩大成所有研究问题的物理引擎。当前 WebGL 保留；进动、NR 数据、后验与层次推断分别接入专门工具。

| 项目 | 本轮确认的价值 | 对本课题的接入方式 | 许可与限制 |
|---|---|---|---|
| [GYOTO](https://gyoto.obspm.fr/) | Metric / Astrobj / Screen / Spectrum 分层，解析和数值度规，Python/C++ 扩展 | 优先借鉴架构；作为离线光学交叉验证候选 | GPL-3.0-or-later；官方主要支持 Linux/macOS，本机尚未运行；有完整数值度规才可追踪相应时空 |
| [SXS](https://github.com/sxs-collaboration/sxs) | 波形模式和 horizons A/B/C 的位置、质量、自旋时序；可固定目录版本 | 第一批真实轨迹适配器，保留质量归一化、坐标规范和时间对齐 | 代码 MIT；数据须引用其版本及适用条款。坐标轨迹不是可观测轨迹，视界文件不是完整时空度规 |
| [GWSurrogate](https://github.com/sxs-collaboration/gwsurrogate) | 部分进动模型返回动力学，模式波形与极化接口 | 与 SXS 并行的模型适配器候选；先选 NRSur7dq4 再逐模型扩展 | MIT/X Consortium；模型文件、C 扩展、质量比/自旋/时长/参考频率范围需分别核对。当前未加载模型 |
| [PESummary](https://pesummary.readthedocs.io/en/stable/gw/tutorials/plotting_from_metafile.html) | 多分析结果容器，后验转换、corner、violin、自旋图、PSD | 借鉴统一“事件—分析—样本—先验”层及模型对照 | LICENSE.md 实为 MIT，不能因 GitHub 标为 NOASSERTION 就误判无许可证；不是推断引擎 |
| [GWPopulation](https://colmtalbot.github.io/gwpopulation/) | 标准群体模型、层次似然、选择效应、Bilby 与加速后端 | 独立核对现有重要性重加权；正式超参数分析时使用标准接口 | MIT；不能补回精简教程缺失的搜索阈值、运行标签与波形版本 |
| [Bilby](https://bilby-dev.github.io/bilby/index.html) | 先验、似然、采样器与单事件/超参数研究接口 | 下一阶段完整超参数边缘化与采样诊断 | MIT；耗时由采样问题决定。此轮没有启动新的参数估计 |
| [qnm](https://github.com/duetosymmetry/qnm) | Kerr 复准正常模频率与球/扁球谐混合 | 后续独立遗迹振铃面板 | MIT；本征频率不决定激发振幅，也不提供并合轨道。未运行 |
| [EinsteinPy](https://github.com/einsteinpy/einsteinpy) | 度规与测地线的模块化教学/计算接口 | 低维轨迹交叉检查候选 | MIT；不提供通用双黑洞并合解，未运行 |
| [kyleyhw/black_hole](https://github.com/kyleyhw/black_hole) | 相机标架、float64 验证目录、CPU/GPU 检查组织 | 参考守恒量、弱偏折、临界冲量参数和收敛测试思路 | 当前固定提交未发现明确许可证；未复制实现，未执行其测试；不采信未经验证的 exact/accurate 宣称 |
| [samrht/blackhole-sim](https://github.com/samrht/blackhole-sim) | WebGPU、shader/CPU 一致性测试、步数与高阶光子环误差讨论 | 参考精度诊断；暂不替换兼容性更好的现有 WebGL2 路径 | 未发现明确许可证；未复制实现，未运行；浏览器和 GPU 兼容性仍需检查 |

本轮确认了“相似名字不代表同一项目”：algometrix、samrht 的 blackhole-sim 分属不同仓库，许可和物理实现不能互相套用。代码许可证也不能自动覆盖模拟模型文件和观测数据。

## 固定提交清单

| 仓库 | 本轮读取提交 |
|---|---|
| bilby-dev/bilby | `2a3306a08d81aa4caf264a7229ee08275a119a75` |
| ColmTalbot/gwpopulation | `b3a34f9dbae5815f13d9c2198f203ac390c05cb7` |
| duetosymmetry/qnm | `c148c7b3e7705415b6b4afbfd512c8751b52a778` |
| einsteinpy/einsteinpy | `b9e3fc13cb5bdcb75accbb0b796489ced1476a2d` |
| gyoto/Gyoto | `b77286cab809997d179506c5b871b3fd6f49ef45` |
| kyleyhw/black_hole | `96a6835f2624bc6425682f723a050f397893cb6d` |
| pesummary/pesummary | `ff4e13428bb2d00fc893e40e447318f92995ecef` |
| samrht/blackhole-sim | `7b0a4d32a01a2367aabffb57ddf48117fd91230b` |
| sxs-collaboration/gwsurrogate | `13cd6a5a4667e0dcae626ddd11677a0d140b54eb` |
| sxs-collaboration/sxs | `b7d25d5aadc2c15ffa4e989fe1c2f661040660b6` |

上述元数据中没有 archived=true。没有将最近提交日期当作正确性或长期维护承诺，也没有执行候选仓库的安装脚本。

## 已实现的泛化基础

入口：<http://127.0.0.1:8765/workbench.html>。原有波形实验保留在根路径，并增加工作台入口。

1. `blackhole-scene/1` 契约：显式时间单位、长度单位、参考质量、质量参考系、来源类别/版本、坐标说明与逐帧天体 ID、质量、位置、自旋向量。最多 30000 帧、每帧 12 天体；浏览器导入上限 5 MB。
2. `sampleScene` 使用非均匀时间轴查找与显示插值，不假设所有模型都是 4096 Hz。天体 ID/数量改变时不插值拓扑；同 ID 位置和自旋的线性插值仅为回放方法。
3. `opticalCapability` 明确拒绝两个自旋黑洞、任意轴单 Kerr、多于两天体、超出当前积分尺度的场景，保留原始数据并切到坐标/自旋视图。不会把自旋清零以假装完整支持。
4. 单 Kerr 与 Schwarzschild 复用已有 MIT 渲染源码；Kerr 内盘边界随顺行 ISCO 改变。通用 `setFrame` 不要求有波形数组；原有 `update` 接口仍服务 Ripple 同步动画。
5. 内置牛顿偏心轨道是接口示例，含明确说明：e=0.25，无辐射反作用，不是偏心 GR 波形。任意自旋示例只展示向量，不声称提供进动动力学。
6. 导入不执行脚本、不读取文件中的远程 URL；来源字符串使用 textContent。格式或能力检查通过均不代表数据经过科学验证。
7. 导出场景与群体研究状态，后者记录输入快照 SHA-256、假设、先验、选择扰动和结果。

协议说明与示例见 `../simulator/SCENE_CONTRACT.md`。现阶段尚未完成 SXS/GWSurrogate/PESummary 的生产适配器，也没有将任意波形自动转换为三维强场轨迹。

## 本轮新增研究：后验全域的蒙特卡洛精度

输入仍是原有 `results/model.json` 的 81 个固定模型权重矩，没有更换事件、注入、选择阈值或形成渠道假设。新程序为 `generalization_study.py`，结果为 `results/generalization.json` 和 `generalization-summary.json`；浏览器使用同一输出快照。

对每个 f 网格点重新计算最小事件 ESS、注入 ESS、近似 σMC(log L)。沿用事件 ESS≥50、注入 ESS≥4N、σMC≤1 的探索性筛查，并计算通过点占归一化条件后验的质量。

- 基准通过质量：**0.8669536553（86.7%）**。
- **39 / 81** 个条件模型的通过质量达到90%。90%是新增汇总阈值，并非发表标准。
- 该比例是积分质量覆盖，不是模型正确概率、证据或新的可信区间。不能在失败区间删点后重新归一化，程序没有这样做。
- 选择效率 k 改变时，界面同时更新选择矩的一、二阶权重与相应 MC 估计；不是沿用 k=1 的诊断假装重新计算。

记 ξ(f,k)=f k ξA+(1−f)ξI，条件事件数为 N。对选择扰动的对数后验导数（忽略与 f 无关的归一化项）为

`s(f,k) = −N f k ξA / ξ(f,k)`。

于是 `d E[f] / d ln k = Cov(f, s)`。基准值为 **−1.4206707316**，含义是在局部一阶近似下，效率比增加约1%，后验均值约下降1.42个百分点。它不是中位数导数，也不是已经测得的选择系统误差。

81组解析协方差导数与 `ln k` 上步长1e−5的中心差分对照，最大绝对差 **2.42e−8以内**。最初1e−4的差分步长未满足既定1e−7阈值，减小步长后达标；没有放宽阈值。

新结果强化的是精度与选择函数审计的必要性；没有解决原输入的搜索选择不匹配风险，也没有获得可靠的物理形成渠道比例。

## 检查与保留限制

- 8 项数据契约与统计检查：单位/来源、非均匀时间、质心、拓扑、能力回退、81模型原始区间复现、导数交叉验证、忽略选择时的校准不变性。
- 27 项工作台浏览器检查：Kerr自旋、动态轨迹、回退、导入/错误保留、标签转义、研究交互/导出、11项目卡片、320–1024布局。
- 原有18项透镜检查和23项波形/真实数据检查通过。
- 已检查桌面、研究面板和移动端截图。上游 shader 的 potentially uninitialized `f_photonDirection` 警告仍存在；没有观察到页面脚本或 WebGL error。此结果不是光学定量误差验证。

## 下一步建议及完成标准

1. **SXS adapter**：先选一个非进动和一个进动数据集，固定目录 tag、模拟 ID、分辨率和数据版本；输出真实坐标轨迹及自旋，并记录波形与近区轨迹的时间对齐方式。
2. **GWSurrogate adapter**：固定模型文件校验值与支持域；检查 q 的方向、总质量时间尺度、模式/极化约定和参考频率，验证返回动力学的坐标框架；不能把训练域之外的外推默认为安全。
3. **光学参考**：优先比较 Schwarzschild 临界冲量参数、弱偏折、Kerr 暗影边界，再检查不同步长/步数与视距的收敛；通过后再考虑相机标架和任意自旋轴。
4. **研究输入升级**：从原始发布恢复事件分析标签、先验、观测运行、注入阈值和选择掩码，再用 GWPopulation/Bilby 联合边缘化超参数并做后验预测检验。先解决样本/选择一致性，再扩大模型库。
