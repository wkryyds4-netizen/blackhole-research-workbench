# 双黑洞引力波实验室

这是现有黑洞研究项目的交互模块。物理波形由 Ripple 计算，透镜渲染复用 MIT 的 algometrix/blackhole-sim；另保留独立实现的 Canvas 轨道示意。它支持非进动双黑洞的参数实验、保存对照、波形声化和 GW150914 真实数据的教学形状拟合。

## 启动

已接入[波形比较](http://127.0.0.1:8765/workbench.html?view=waveform)：SXS0389 / NRSur7dq4 复h22、PSD/频带/渐消窗敏感性、NR外推阶参照与完整导出。第三阶段结果和复现步骤见 [研究记录](../research/stage3-results/README.md)。

新增[通用研究工作台](http://127.0.0.1:8765/workbench.html)：单 Kerr / Schwarzschild、通用 JSON 轨迹回放、自旋能力回退、81 条件模型的先验与选择效应诊断，以及相关项目接入清单。协议见 `SCENE_CONTRACT.md`，研究与项目结论见 `../research/generalization-review.md`。原波形实验仍在根路径。

本工作区已准备隔离环境 `.venv` 和公开数据。Windows 双击 `start.cmd`，待出现 `Local:` 后打开 http://127.0.0.1:8765 。关闭终端或按 Ctrl+C 停止服务。只监听本机回环地址。

重新安装时使用 Python 3.12 环境：

```powershell
cd simulator
python -m venv .venv
.venv/Scripts/python.exe -m pip install -r requirements.txt
.venv/Scripts/python.exe prepare_event.py
.venv/Scripts/python.exe server.py
```

如果端口占用，执行 `.venv/Scripts/python.exe server.py --port 8766` 并打开相应端口。数据下载只发生于手动运行 `prepare_event.py` 时，服务器不向外部上传参数。

## 功能和借鉴

新增透镜视图：程序化星空、弯曲光线、盘的多重像、HDR 与 bloom；提供三档画质、环境盘开关、镜头旋转/缩放和轨道示意对照。画面与应变共用时间轴，按质量比绕质心运动。暂停冻结纹理时钟并在 128 帧累积完成后停止渲染，隐藏页面/切换视图停止，WebGL 不可用或丢失时回退 Canvas。

**适用范围：**双黑洞光学是无自旋两中心近似，自旋输入仍作用于 Ripple 波形。仅主黑洞有示意盘，不代表 GW150914 的观测证据；并合最后 12 ms 和遗迹半径采用视觉过渡，不能测量真实遗迹参数。未接入上游独立 Peters 轨道、艺术网格波或音频 chirp。详见 `../research/blackhole-sim-analysis/README.md`。

构建产物已放在 `public/lensing/`，正常启动不需要 Node。修改渲染源后重建：

```powershell
cd renderer-src
npm ci --ignore-scripts --no-audit --no-fund
npm test
npm run build
```

源码：`renderer-src/adapter.js`、`scene-state.mjs`；上游许可证和固定提交见 `THIRD_PARTY.md`。光学浏览器验证：`node simulator/qa-lensing.cjs`（从工作区根目录，服务器须已启动；脚本中的 Playwright 路径适用于当前机器）。

| 原项目的经验 | 本项目改造 |
|---|---|
| binaryBHexp 的轨道、自旋与波形联动 | 一条物理时间轴驱动相位示意、极化粒子环和波形；播放倍率独立 |
| Cardiff 的模板与观测叠加 | 使用两台探测器的数据；加入白化、残差、时频图和可追溯处理 |
| 黑洞浏览器模拟器的参数与场景组织 | 可拖动镜头，独立控制物理倾角；参数自动重算，保存对照 |
| 科学波形库与显示分层 | Ripple 双精度 IMRPhenomD / IMRPhenomXAS，后端缓存；导出包含模型与版本 |

`public/index.html`、`style.css`、`app.js` 是界面；`physics.py` 是计算适配器；`server.py` 提供本地接口；`prepare_event.py` 保存数据和处理记录。切换 PyCBC 时只需实现同一波形数据接口，不应改动前端物理时间轴。

## 物理约定

- 质量为探测器参考系，单位太阳质量；源参考系质量需乘以 `(1+z)`。光度距离单位 Mpc。
- 自旋参数为沿轨道角动量方向的无量纲分量，限于 −0.8 至 +0.8；本工具范围 `m1>=m2`、`m1/m2<=4`。
- 两个模型都只在本界面使用准圆、非进动、主导四极波形。不支持偏心、进动或更高阶模式。
- 频域计算使用32秒长度对应的频率网格，4096 Hz采样。20–25 Hz以及900–1000 Hz余弦型平滑收边避免硬截断振铃。收边是信号处理，不是物理源的突然开始或结束。
- 逆傅里叶变换按连续变换单位乘采样率，时间原点位于正面观测双极化包络的最大值。导出的 `hp` 和 `hc` 数值单位为 `10^-21`。
- 粒子环按该波形峰值归一化放大；轨道相位取自波形正交分量，间距由 Kepler 关系估计。并合前12毫秒进入显示过渡，球体尺寸与遗迹外观均不作定量解释。
- 画面的旋近频率来自相位导数，只在幅度和时间筛选有效时显示。没有把并合后的振铃频率称为轨道频率。
- 声音用应变归一化后按原采样率播放；振幅不是探测器的真实声压。动画可减速，导出的物理信号保持不变。
- 默认“GW150914类参数”是教学配置，不是从完整联合后验选取的事件重建。

## 真实数据对照

来源是 [GWOSC 的 GW150914 v3](https://gwosc.org/eventapi/html/GWTC-1-confident/GW150914/v3/)，每台探测器32秒、4096 Hz。下载URL、字节数和SHA-256保存在 `data/manifest.json`；原始压缩文本和元数据均保留。

处理使用首尾各10秒数据独立估计 Welch PSD（2秒片段，中位数平均），将两估计平均并插值到完整频率网格。整段去趋势和 Tukey 收边后白化，再以四阶 Butterworth 35–350 Hz 双向滤波。图示窗口避开整段数据边缘。32秒的噪声谱只适合此教学对照，不替代正式事件分析的长时噪声与数据质量审计。

拟合使用两个正交模板分量，在相对GPS 1126259462.4的 −0.20 至 +0.08秒窗口中解振幅系数，搜索约 ±40毫秒平移。模板与数据经过同一白化和带通流程。相关系数是该窗口内的形状指标，受到所选窗口、PSD、自由对齐和噪声影响。

每台探测器独立拟合时间、相位、振幅。没有强制两个探测器共同天空位置和传播时延；这不是物理探测器投影，不是相干分析，也不是事件参数的后验估计。由于振幅被自由拟合，改变距离和倾角不能在此模式中用于推断距离或倾角。

## 验证

运行 `.venv/Scripts/python.exe validate_physics.py`，结果写入 `validation.json`。检查距离反比缩放、倾角响应、真空质量缩放、啁啾时间方向、自旋的实际波形影响、双模型输出、参数边界与真实数据残差恒等式。

浏览器验证脚本 `qa.cjs` 检查参数交互、模型切换、时间游标、播放、对照、导出、真实数据拟合、手机布局和错误日志。截图与报告保存在 `qa-output`。这些检查不等同于独立 LALSuite 或数值相对论精度交叉验证；正式研究使用前需要补充后者。

## 后续研究接口

现有 `../data/posteriors.npz` 的精简事件后验缺少进动方位角、参考频率和完整外禀参数，不足以直接重建事件探测器信号。接入联合后验后，可把多条波形叠成不确定性带；接入与实际搜索阈值匹配的注入集后，再展示探测前后总体分布。偏对齐与各向同性代理分量不能直接标记为已确认的形成渠道。

参考与许可见 `THIRD_PARTY.md`。相关项目调研见 `../research/gravitational-wave-project-review.md`。

## 通用工作台与 SXS 数据

打开 `/workbench.html` 可在真实 SXS 轨迹、解析/示意场景和群体稳健性分析间切换。`/workbench.html?dataset=sxs-0003-lev5` 直接载入历史进动算例。下拉菜单提供 0001、0003 与 0389 的两个分辨率；质量、自旋、分离度曲线与场景时间联动，支持来源完整导出。

0003 已被官方标记为旧版，只用于接口演示。当前新进动基准为 2327 v3.0，提供三个分辨率，总计七份 SXS 场景。数值轨迹是坐标数据，不是直接观测图像；单 Kerr 渲染只是定常参考。

`?view=alignment` 展示 0389 三种对齐规则与独立晚期验证结果。`?view=dynamics` 展示隔离 Zig 环境中实际生成的 NRSur7dq4 方向、自旋及 h22 模式；−100 M 后不展示非物理延伸自旋。代理输出没有位置或分离度。最新研究、误差与复现汇总见 `../research/stage2-results/README.md`；数据协议见 `SCENE_CONTRACT.md`。
