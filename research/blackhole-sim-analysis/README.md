# blackhole-sim 实测与改造记录

记录日期：2026-10-23。上游：https://github.com/algometrix/blackhole-sim
固定版本：`0ee5541ec4ffcd9ca708f37c080bf00b7f3ef6d1`，MIT，Copyright 2026 algometrix。

## 已完成

- 原始源码保留在 `../blackhole-sim-upstream/`，未改动。
- `npm ci --ignore-scripts --no-audit --no-fund` 安装成功。
- `npm test`：19 个测试文件，152 项通过。README 的 100 项徽标落后于实际测试数。
- `npm run build`：TypeScript 检查与 Vite 构建成功；主包约 710 kB，存在大包提示。
- 原项目运行地址：<http://127.0.0.1:5174/>。进程退出后需重新启动。
- Edge / WebGL2 实测正常，GPU 为 ANGLE / NVIDIA RTX 5070 Laptop / Direct3D11。HUD 的约 165 fps 仅是一次读数，不是性能基准。
- 检查默认、Kerr、双黑洞并合、曲率网格预设，截图保存在本目录。预设检查没有页面异常或 HTTP 失败记录。
- 首次运行见 shader 的 potentially uninitialized `f_photonDirection` 警告及一次未记录地址的 404；后续预设检查未复现 HTTP 失败。尚未证明警告根因，不能当作定量准确性验证。

脚本：`../run-blackhole-inspection.cjs`、`../inspect-presets.cjs`。浏览器结果：`initial-browser.json`、`presets.json`。

## 架构与可借鉴之处

| 模块 | 工作 | 改造价值 |
|---|---|---|
| `render/blackHolePass.ts` + `shaders/geodesic.frag` | 每像素积分光线，计算暗影、盘的多重像、背景透镜 | 直接复用渲染核心 |
| `shaders/kerr.glsl` | 单 Kerr 黑洞的 Hamiltonian 光线积分 | 后续独立单黑洞实验；不能直接用作双 Kerr 解 |
| `render/starfield.ts` | 启动时烘焙程序化 HDR 星图立方体 | 本地资源、无需外部纹理 |
| `render/pipeline.ts` | 光线追踪、覆盖场景、HDR 合成、bloom、色调映射 | 复用画质控制与暂停累积 |
| `render/bloom.ts` | 六级双滤波 bloom | 控制光晕，避免过曝遮住双暗影 |
| `main.ts` | 固定步进，最大发散补帧限制，性能降级 | 借鉴资源与性能管理，不接入第二套物理时钟 |
| `ui/presets.ts` | 数据驱动场景参数与相机 | 借鉴可重复预设 |

低/中/高：内部比例 0.5 / 0.75 / 1；光线步数 128 / 256 / 420。Schwarzschild 低档使用 RK2，中高档 RK4。静止画面可累积至 128 帧。积分步长由距离启发式控制，并非误差容限控制。

## 必须保留的物理边界

1. 双黑洞光学为两个 Schwarzschild 中心的加速度叠加近似，不是双黑洞时变度规。自旋 Kerr 与第二黑洞不能同时启用。
2. 上游轨道采用圆轨道 Peters 衰减和 Kepler 角速度，显式步进；主黑洞固定在原点，伴星绕其运动。注释里的 exact 不代表完整相对论双黑洞轨迹。
3. 网格基于单黑洞 Flamm 嵌入曲面，再加艺术化传播速度与正弦波纹；它不是探测器的 h+ / h×，也不是三维空间实际下陷。
4. 并合后半径的阻尼呼吸是视觉过渡，不是求解出的准正常模。
5. 吸积盘是薄盘形状、噪声纹理、近似色温与 Doppler/红移着色，不是辐射流体动力学；只环绕主黑洞。不能暗示 GW150914 有已证实的明亮吸积盘。
6. CPU 测试通过并不构成 GLSL 渲染的定量误差验证；Kerr 极限步数耗尽时的捕获判定也可能影响边界。

## 本地已有模拟器

`../../simulator/`：Ripple 0.2.1、IMRPhenomD / IMRPhenomXAS，统一波形、轨道示意、极化、音频、真实 H1/L1 GW150914 对照。
地址：<http://127.0.0.1:8765/>；重启：`simulator/start.cmd`。
此前浏览器验证 23 项通过。波形算法和 GWOSC 数据处理本次保持原样。

## 接入方案

- 复制 MIT 渲染模块的依赖闭包到 `simulator/renderer-src/vendor/`，保留许可证、版本与来源；原始仓库不动。
- 浏览器用本地构建包；提供透镜渲染和现有 Canvas 轨道示意两种模式。
- 不调用上游 `stepWorld`、Peters 轨道或合成 chirp。渲染只读取我们的波形采样时间、相位与 Kepler 间距。
- 取长度单位 GM_total/c²，两个 Schwarzschild 半径分别 2m₁/M、2m₂/M。通过相机与第二黑洞相对主黑洞的坐标变换，实现绕质心的画面。
- 自旋参数继续影响波形，但双黑洞光学保持无自旋近似。并合视觉过渡不能输出未经计算的遗迹自旋。
- 吸积盘明确标注为可关闭的环境示意。画质、相机与吸积盘控制不改变科学数据。
- 时间跳转、相机、参数变化重置累积；暂停冻结纹理时钟；隐藏页面停止渲染；WebGL 不可用回退 Canvas。
- 验证时间同步、绕质心、暂停/倒放、相机倾角分离、画质/吸积盘不改波形、移动端布局与 WebGL 错误。

## 恢复工作

用户曾要求保存，随后明确“继续研究”。以下首轮改造现已完成，不需重新下载、重复基础测试或改动已有波形后端。

## 已落地的首轮改造

- 完整上游版本另存 `upstream-0ee5541.zip`，原始 checkout 的 `git status --short` 为空。
- MIT 依赖闭包 25 个文件已复制到 `simulator/renderer-src/vendor/blackhole-sim/`，保留逐文件 SHA-256 与许可证；上游源码未改。
- `adapter.js` 接管质心相机变换、画质、可选环境盘、暂停/隐藏节能与异常回退。`scene-state.mjs` 为无状态时间采样函数，可确定性反向拖动。
- 默认开启透镜视图；可切回原 Canvas 轨道示意。窄屏扩大垂直视场，减轻双黑洞被裁切。
- 暂停后累积 128 帧再停止渲染；场景/相机变化重新累积。切换页面或关闭透镜停止渲染。模拟 WebGL 上下文丢失可回退轨道示意。
- 波形后端、公开数据与拟合算法未改。镜头不会改变物理倾角，显示选项不改变波形数组。
- 构建包 `simulator/public/lensing/renderer.js` 约 686 kB（gzip 158 kB），较上游完整应用的 gzip 188 kB 小；未引入上游交互物体、音频、GUI 与独立动力学主循环。包大小不是性能测试结果。
- 未移植网格波纹、恒星撕裂或单 Kerr 面板。这些功能需要另设明确的物理问题与时间尺度，避免混入当前双黑洞波形实验。

### 最终验证

| 检查 | 结果 | 保存记录 |
|---|---|---|
| 上游测试 | 19 文件 / 152 项通过 | 本文记录，本地可 `npm test` 重现 |
| 新增质心与时间采样测试 | 3 项通过 | `simulator/renderer-src/scene-state.test.mjs` |
| 透镜浏览器验证 | 18 项通过 | `integrated-lensing-validation.json` |
| 原有功能回归 | 23 项通过 | `integrated-regression-validation.json` |
| 最终页面脚本、HTTP、WebGL error | 未记录到 | 透镜验证 JSON |
| 原项目继承的 shader warning | 仍存在，三档画质各一次 | `f_photonDirection` 的 potentially uninitialized 警告 |

最新截图：`integrated-desktop.png`、`integrated-mobile.png`、`integrated-remnant.png`。运行入口 <http://127.0.0.1:8765/>；原项目对照入口 <http://127.0.0.1:5174/>。它们是本机进程，进程停止后需重启。

### 下一步研究建议

1. 优先测量光学误差：固定单黑洞相机，比较临界冲量参数、暗影边界与高精度离线光线积分；量化步数、步长和有限逃逸半径的影响。先建立误差预算，再调整 shader。
2. 为双黑洞显示限制可解释范围：按分离度比较两中心近似与更高阶近似/数值时空的光线追踪，不将视觉一致性当成物理验证。
3. 若研究轨道与自旋，单独接入经验证的进动波形或 NR surrogate 轨迹；当前两种波形后端都未提供完整进动轨迹。
4. 环境盘应独立演进为“可视环境”模块；若要解释电磁辐射，需要有明确参数与校准的双黑洞盘/辐射模型。
5. 形成渠道稳健性研究继续在群体后验与探测选择函数层进行；可视化只负责展示假设、样本和结果之间的变化。
