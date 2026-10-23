# blackhole-scene/1 接入协议

本协议用于数据回放和视图选择，不是通用引力求解器。外部模型先将数据转换成此格式，再使用工作台的“导入场景”。旧波形实验导出的 JSON 不是此协议，不能直接作为数值轨迹导入。

```json
{
  "schema": "blackhole-scene/1",
  "title": "单黑洞接口示例",
  "units": {
    "time": "GM/c^3",
    "length": "GM/c^2",
    "reference_mass_msun": 65,
    "mass_frame": "geometric"
  },
  "provenance": {
    "kind": "illustration",
    "source": "接口示例，不是观测重建",
    "version": "example/1",
    "frame": "Cartesian coordinates, spin along -y"
  },
  "frames": [
    {
      "t": 0,
      "bodies": [
        {"id": "BH", "mass": 1, "position": [0, 0, 0], "spin": [0, -0.7, 0]}
      ]
    }
  ]
}
```

- `mass` 是 `reference_mass_msun` 的倍数，不是太阳质量数值；位置单位 GMref/c²。
- 时间单位为 `s` 或 `GM/c^3`。几何单位转秒乘 `reference_mass_msun * 4.9254909476412675e-6`。
- `mass_frame` 必须为 `source` / `detector` / `geometric`。本页面不自动乘红移，也不推断未提供的红移。
- `provenance.kind` 为 `analytic` / `numerical` / `illustration` / `waveform-derived`。来源、版本、坐标说明不能为空。
- `spin` 是笛卡尔无量纲自旋向量，模长≤1；通用坐标视图保留其方向。光学单 Kerr 只接受 −y 轴、χ≤0.998。
- 时间严格递增，不要求均匀。1–30000帧，每帧1–12天体；同一帧 ID 唯一。拓扑改变时使用离散切换，不创造中间天体。
- 光学积分当前只接受质量比参考尺度内 `0.15≤mass≤1.01`，位置半径≤22；双体分离度须大于等于显示半径和且≤25。超范围仍可在坐标视图查看。这些是当前渲染器的显示能力约束，不是天体物理边界。
- 同一拓扑的质量、位置、自旋采用线性显示插值；若模型需要四元数/旋转框架插值，应在适配器中采样到足够密度，并说明坐标约定。
- 视界轨迹不能替代完整时空度规；导入数值相对论坐标不会使两中心光学变成精确 NR 光学。
- JSON上限5MB，本地解析。文件格式合法不代表科学数据有效。不要把未知外禀参数、轨道相位或参考频率补成“测量值”。

代码接口：`public/scene-contract.mjs` 提供 `validateScene`、`sampleScene`、`opticalCapability`、`opticalFrame`。`renderer-src/adapter.js` 的 `setFrame` 消费已检查的光学帧；旧 `update` 继续接受现有 Ripple 波形状态。

## SXS 数据适配（2026-10-23）

`../research/sxs_adapter.py` 直接读取已校验的 SpEC `Horizons.h5`，不依赖完整 SXS Python 包。七份本地快照通过 `public/datasets/catalog.json` 注册；`?dataset=sxs-2327-lev3` 可打开未被目录标记为弃用的进动基准，0003 保留为历史演示。

参考总质量取原始元数据的 `reference_mass1+reference_mass2` 或旧版 `relaxed_mass1+relaxed_mass2`。导出时间为 `(t_source−t_common)/Mref`，位置和视界 Christodoulou 质量均除以 Mref，自旋保持无量纲。位置和自旋同时经过正旋转 `(x,y,z)→(x,−z,y)`，无坐标质心平移。65 M☉ 只是一种尺度选择。

可选 provenance 字段包括 `source_common_horizon_time`、`rotation_matrix`、`source_files`（原文件 URL 和 SHA-256）、`license`、`deprecated`、`superseded_by` 与 `gaps`。`gaps` 是 `{start,end,behavior}` 数组，时刻使用场景单位。当前 SXS 缺口仅出现在 A/B 结束至 C 开始之间；保持最后 A/B 帧并显示缺口说明，不跨拓扑插值。其他模型若有同拓扑缺口，适配器必须另外处理，不能假定该字段自动改变插值器。

每个量在自身时间网格上插值，禁止外推。自适应压缩相对原始分段线性数据的每分量容限为位置 0.003 Mref、质量 1e-5 Mref、自旋 1e-4；这些不是数值相对论模型本身的误差界。共同视界首个样本不等于波形峰值，未与波形迟滞时间同步。0003 已被官方弃用并由 2327 替代，作为历史接口演示保留。

原始文件、来源清单、检查报告与复现说明位于 `../research/sxs-data/` 和 `../research/adapter-results/`。

2327 v3.0 的 Lev1/2/3 来自 CaltechDATA，记录未声明许可。适配器 v2 保留每级 `pbj_info`、原始参考时刻和参考轨道频率，按清单填写 repository/license；不把各分辨率的不同参考时刻等同，不把0003与2327的差异作为分辨率误差。重新下载时同时运行 `research/fetch_sxs.py` 与 `research/fetch_reference.py`，再运行适配器。

## 独立研究数据类型

`blackhole-alignment-study/1` 用于保存比较协议。`time_M` 为共同的 Lev3 原始坐标时钟/Mref，`methods` 各含固定 Δt、φ、训练/验证指标，以及三类验证残差序列。相对向量差记录非负范数，分离度差为绝对值，平面相位残差有符号。`training_window_M`、`validation_window_M` 与原文件校验必须随导出保留；它不是场景轨迹输入。

`blackhole-dynamics/1` 用于没有位置输出的代理模型。`t` 是几何时间数组，`q_copr` 是 N×4 scalar-first 单位四元数，以主动旋转 q v q⁻¹ 从共进动到惯性系；`chiA/chiB` 是惯性系 N×3 无量纲自旋；`orbphase` 是共进动框架轨道相位；`h22_real/imag` 是惯性球谐模式，单位 r h/M。共轨道 +x 方向由 q 旋转 [cos(orbphase),sin(orbphase),0] 获得。它仅表示方向，没有距离量纲。

代理面板选择最近的原始采样帧，不插值四元数，也不推算位置、分离度、物理质量或距离。当前 NRSur7dq4 的 `dynamics_valid_until_M=-100`，之后隐藏自旋和框架方向，保留波形。截止比较允许 1e-7 M 的浮点容差，以容纳模型网格的约3e-9 M数值偏移。t=0采用模型多模式总振幅峰值约定，不能与SXS共同视界时钟直接联动。固定参数、参考频率的特殊标记、模型/源码哈希与适用范围在 metadata 中保存。
