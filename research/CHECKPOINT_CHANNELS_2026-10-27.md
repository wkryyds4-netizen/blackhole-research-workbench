# 五渠道阶段暂停与 Git 检查点

按用户要求暂停研究并保存现有版本。恢复研究须等待用户明确指示；整体科学目标尚未完成。先前检查点 `68ab257` 及冻结 v0.3 文稿保留。

本次保存：

- [24 组五渠道条件后验](population-channels/RESULTS.md)，全 69 事件、11,180,578 行，五渠道 CE/CHE/GC/NSC/SMT；两隐藏自旋条件、两 O3 银行、两搜索口径、三正则化概率。
- 最大分支中位数变化：正则化 8.65、隐藏自旋条件 3.41、银行替换 0.61、搜索口径 0.43 个百分点。独立 Metropolis 与 QMC 的最大分位数差约 0.17 个百分点。结果仅适用于记录的模型条件、参数域与保存银行。
- [中英 LaTeX 研究稿](paper/channels-2026-10-27/README.md)，各 7 页，及 `output/pdf/` 和工作台中的两份 PDF。
- 通用混合模型数据合同、可切换面板、含科学量定义和数据 SHA256 的会话导出。`npm test` 全部通过，包括 24 组原始结果对照。
- 原始数值失败、修复、独立验证和历史积分协议全部保留；旧直接权重缓存不可混入 `channel-log-*` 协议。

暂停前没有启动新研究。原始 PE 行、注入银行、模型权重、下载源码、环境和本机会话记录不进入公开 Git。公开范围沿用根 `NOTICE.md` 中的归属与许可证边界。

恢复时先读 `population-channels/RESULTS.md` 与 `measure-method.md`。下一项实质研究是对本批渠道积分做成对有限银行传播；当前事件贡献集中最低约 5.24，选择贡献集中最低约 218.21，尚未获得实际误差校准。现有汇总矩没有保存跨隐藏自旋族的全部配对协方差，不能假装已有逐行配对缓存，也不能套用旧 g/r 抽稀范围。历史训练/选择参考、缺失 PE 支持、真实渠道隐藏自旋和波形误差仍待解决。

This is a user-requested pause and public Git checkpoint, not completion of the overall research goal. Resume only on a new user instruction.
