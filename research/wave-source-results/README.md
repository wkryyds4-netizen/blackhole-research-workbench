# 第三阶段：公开数值相对论辐射波形读取

2026-10-23。此目录是源数据适配层；不把解码成功等同于代理模型精度验证。实际对照由上层分析保存自己的窗口、参考参数、旋转与时间平移。

## 已获取的固定版本

| 数据 | 原始产品 | 体积 | 本阶段用途 | 许可 |
|---|---|---:|---|---|
| SXS:BBH:0389，Lev3，[Zenodo 3274661](https://doi.org/10.5281/zenodo.3274661) | rhOverM_Asymptotic_GeometricUnits_CoM.h5 | 88,045,582 B | N2/N3/N4 模式对照与外推阶敏感性 | 记录声明 CC-BY-4.0；使用须归属 SXS 原记录作者并保留 DOI |
| SXS:BBH:2327 v3.0，Lev2，[CaltechDATA f70bw-gcd50](https://doi.org/10.26138/SXS:BBH:2327v3.0) | Strain_N2.h5 + JSON | 4,425,307 + 2,231 B | 解码及进动框架研究准备 | 记录未声明许可；仅本地研究快照，不推定可公开再发布 |
| 同上，Lev3 | Strain_N2.h5 + JSON | 4,253,625 + 2,233 B | 解码及进动框架研究准备 | 同上 |

原始记录保存在 `../adapter-sources/`。`0389-download-manifest.json` 与 `download-manifest.json` 保存 URL、原始字节数、官方 MD5 和本地 SHA256；下载器已逐个校验大小和 MD5。NPZ 同名 JSON 保存原始输入与导出 NPZ 的 SHA256。2327 的旧编号关联/跑次日期仅说明其来源沿革，不能据此证明是代理训练集外的数据。

## 接口与无损读取

`../load_nr_waveform.py` 提供 `load_nr_waveform(path, group, metadata_path)`，返回 `(t, modes, info)`；`modes[(ell,m)]` 是复数数组。NPZ keys 为 `t`、`h_2_-2` 等。所有产品均保存完整 l=2..8 的 77 个模式，无天区投影、不裁去前期瞬态、不平移时钟、不作额外空间旋转。时序严格递增且全部有限。主对照应自行限定有效时段、频段与 l≤4 模型公共模式。

0389 是传统 NRAR 的 `[time, real, imag]` 三列。逐点验证全部 3×77 个导出模式及时间与原 HDF5 完全相等。2327 是 `rotating_paired_diff_multishuffle_bzip2`；用与数据写入版本一致的官方 `sxs==2024.0.25` 读取器，执行其解压、解配对及共旋到惯性框架转换。没有自行实现或近似 Wigner 旋转。

独立检查共旋与惯性各 l 多重态的功率不变量，最大相对误差（按该 l 的峰值功率归一化）为 Lev2 `3.14e-15`、Lev3 `3.57e-15`，四元数范数与原时序也通过。验证详情见 `validation.json`。这检验读取与旋转一致性，不等于物理误差评估。

## 0389 的物理约定与证据

- HDF5 N2/N3/N4 属性 `FrameType=1` 为惯性系，`DataType=1` 按该文件 GWFrames 枚举为 strain，`MIsScaledOut=RIsScaledOut=1`。输出是 spin weight −2 的 `r h_lm/M`，复数 convention 为 `h_plus - i h_cross`。不是某个探测器的 strain。
- `History.txt` 给出 `ChMass=0.999988106301`，与该 Lev3 的参考质量之和完全一致。固定生成源码 [GWFrames commit 188298f 的 Extrapolation.py](https://github.com/moble/GWFrames/blob/188298fa0caf9a889a533cacfbe29cd52b026cac/Code/GWFrames/Extrapolation.py) 的 `ReadFiniteRadiusWaveform` 明确执行 `SetTime(T/ChMass)` 与 strain 的 `UnitScaleFactor=1/ChMass`，所以**不得再次用参考总质量除以已归一化波形**。
- 同一源码先以提取球的平均 lapse 修正时间积分，再减去以 ADM 能量构造的 tortoise radius，之后外推；它是辐射迟滞时间，不是近区视界坐标时钟。源码快照在 `GWFrames-Extrapolation.py`，生成历史在 `0389-history.txt`。
- 本数据是有限半径提取后 N 阶外推到无穷远，并非 CCE。History 保存原提取半径 100..695。N2/N3/N4 比较只测量外推处理敏感性，不是三个网格分辨率，不构成收敛阶或严格误差上界。
- 质心空间平移 `[-2.12387684e-5,-1.38374957e-5,-2.45055831e-6]` 与 boost `[1.61350735e-7,9.24637820e-8,2.18221275e-9]` 已经在存储产品中应用。详细精度以同名 JSON 属性为准，不再重复应用。
- 原始时间范围为 `[-107.76155555281474, 5223.281358314723] M`；离散 l≤4 总振幅峰值位于 `5065.402632488271 M`。该峰值仅供对齐起点参考，后续可在明确的方法下插值求峰；不是共同视界时刻。
- 该旧产品 history 没有新目录的 `add_memory` 后处理记录；这不意味着所有 m=0 数据严格为零，也不应将其称为物理上不存在记忆效应。

## 2327：已解码，尚未完成模型对照

2327 JSON 明确：N2 外推、质心平移/boost、从各自参考时刻加入 memory（Lev3 为 1027 M）、精度截断容差 1e-10。l=2..8，`m_is_scaled_out` 与 `r_is_scaled_out` 为 true。Lev2 有 24233 帧，Lev3 有 24280 帧；源时间与解压时间一致。

相较旧的0389产品，本产品附加了 memory；不能直接将全模式差异归因于 NRSur7dq4 的拟合误差。读取器保留全部模式和后处理信息，并没有静默去除 m=0 或 memory。当前简短 sidecar 不给出旧外推历史的 `ChMass` 标量，比较前还应独立核实 v3 具体质量归一化约定；不得未经证据把视界参考总质量再施加给已经 scaled-out 的模式。进动比较还需匹配辐射共轨框架、参考频率/时刻、两颗黑洞的标签与自旋矢量。目前只完成数据和框架准备，没有公布2327的代理误差。

## “独立”的范围

这里的独立源是从公开 NR 档案独立下载的数值数据，不表示未参与 NRSur7dq4 训练。模型论文 [Varma et al. 2019](https://arxiv.org/abs/1905.09300) 使用早期 890 个进动算例、64 个对齐自旋算例及574个新算例；当前模型HDF5没有可直接查阅的训练 ID 清单。0389、2327 的成员身份尚未证实，所以所有导出 metadata 均写 `Not established`，不能称 held-out validation 或泛化误差。

该论文 III B 脚注还说明训练处理中把辐射 t 与近场视界时间识别是 gauge-dependent 约定；III C 先平移到总模式振幅峰值，再在 −4300 M 将惯性框架对齐辐射共轨框架并旋转自旋。它不是任意使用 metadata.reference_time 的自旋就能复制的处理。低自旋0389可作为受控起点，但仍须报告参考时刻差异和剩余参数假设。

## 复现

```powershell
simulator/.adapters-venv/Scripts/python.exe research/download_nr_waveform.py --0389
simulator/.adapters-venv/Scripts/python.exe research/download_nr_waveform.py
simulator/.adapters-venv/Scripts/python.exe research/load_nr_waveform.py
simulator/.waveform-venv/Scripts/python.exe research/load_nr_waveform.py --dataset 2327
simulator/.waveform-venv/Scripts/python.exe research/validate_nr_waveform.py
```

新读取环境为 `simulator/.waveform-venv`，未改动主服务器 `.venv` 或代理 `.adapters-venv`。SXS wheel 按 `--no-deps` 安装，再显式安装 numpy/scipy/h5py/numba/quaternionic/spherical/inflection/pandas/pytz/tqdm/requests；**没有安装或导入 Julia**。精确已安装版本与官方 wheel SHA256 在 validation.json。SXS读取器为MIT代码；归档研究用外推源码来自其原项目，不作为本工作原创实现发布。
