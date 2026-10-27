# 物理形成渠道研究

当前阶段：[24组完整条件对照](RESULTS.md)（[English](RESULTS.en.md)）已完成，包含全69事件、完整选择积分、两算法后验核对及中英LaTeX/PDF。工作台已接入对应的可切换面板。模型条件、贡献集中与未见支持的限制仍保留，不能将这些权重称为宇宙形成比例。

## 已公开暂停检查点的历史背景

以下保留暂停时的记录；当时的未完成状态不代表上方新阶段进度。

2026-10-27，按用户要求暂停并保存当前版本。[English](README.en.md)。
这是模型接入与数值诊断记录，尚未用观测事件拟合形成渠道比例。

## 已保存内容

- 从 [Colloms et al. 的公开模型](https://zenodo.org/records/14967688)
  接入 CE、CHE、GC、NSC、SMT 五渠道的四维流密度，变量为源系啁啾质量、
  质量比、有效自旋和红移。原始产品与权重留在忽略目录；公开文件保存来源、
  校验和、原创适配器及汇总诊断。
- `channel_flow.py` 使用 Glasflow 和严格匹配的数值 state_dict；
  不执行下载的 AMAZE 源码。`runtime-results.json` 保留运行检查与限制。
- `channel-preview.json` 记录三种出生自旋、每渠道每例 32,768 行的模型预览，
  共 491,520 行。这些是模拟器抽样，不是观测后验或渠道率测量。
- `projected_measure.py` 与 `measure-validation.json` 实现并检查四维渠道密度
  和七维自旋先验之间的测度桥。它需要明确选择未建模自旋的条件分布；
  尚未应用到全 69 事件的渠道似然。
- `sign-threshold-results.json` 比较固定预览样本在不同零附近阈值下的符号比例。
  这是事后诊断，不是预注册检验。

## 必须保留的限制

1. 原始运行门槛并非全部通过：NSC 原 sampled-latent 密度误差
   3.58e-8 超过预设 1e-8。以实际存储物理坐标为参考的补充检查通过；
   原失败保留在 `runtime-failures.json`，没有放宽门槛。
2. CE 使用 `[chi_birth, ln(alpha_CE)]` 的解释得到公开参考样本诊断支持，
   但尚未恢复历史训练执行配置。不能把当前示例配置直接当历史事实。
3. GC 的公开坐标上界 q=1.001，预览有约 0.24%–0.94% 的 q>1 行。
   这些行保留并计数，尚未实施物理截断和重新归一化。
4. 出生自旋为零时，GC 严格负有效自旋比例约 49.51%，
   换成 chi_eff<-0.01 后约 6.92%，零附近带内约 86.61%。
   连续流的展宽会影响符号统计；这些数值不能当真实零自旋原子质量，
   也不能从既有 g 结果直接推成渠道比例。
5. 测度检查最初未分段积分的门槛失败保留在 `measure-history.json`；
   在零点分段后的计算通过原门槛。数学单元检查不等于实际观测分析完成。

## 复现入口与恢复位置

现有本机环境为 Python 3.14、CPU Torch 2.13.0、Glasflow 0.4.1，
并使用 NumPy、SciPy、h5py 和 Requests。来源获取工具依次为
`fetch_sources.py`、`fetch_products.py` 和 `inspect_products.py`；
固定来源 URL、版本与校验和见 `sources/manifest.json` 和 `product-inspection.json`。
这不是洁净环境一键复现承诺；原始产品和抽样缓存不随 Git 发布。

取得并核验上述来源后，可运行 `check_runtime.py`、`sample_contexts.py`、
`check_measure.py`、`diagnose_sign_threshold.py` 重建诊断。
重跑会改写本目录同名输出，应先保留本检查点。
`check_runtime.py` 的 CE raw-context 对照是映射诊断，物理预览采用
`preview-design.json` 中的 log-alpha 约定。

恢复研究时，优先解决 GC 物理域归一化、快速且经验证的自旋投影密度、
事件与注入的一致条件分布和共同支持，再建立完整渠道似然。
密度正则化须同时进入事件项和选择项。不要把本检查点标为渠道推断已完成。

来源对照包括 [Zevin et al.](https://arxiv.org/abs/2011.10057v2)、
[Colloms et al.](https://arxiv.org/abs/2503.03819v1) 及固定的
[AMAZE 源码版本](https://github.com/michaelzevin/AMAZE/tree/4b3a0c78eec8ca4780f6e78fd50f0146ced33878)。
该源码快照未见明确许可证，留在本机，不随本项目再发布。
