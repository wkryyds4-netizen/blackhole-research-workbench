# 来源与许可

- [Ripple](https://github.com/GW-JAX-Team/ripple)：调用安装的 rippleGW 0.2.1，MIT，版权归其作者。依赖包的许可证随安装包保留。论文：[Edwards et al. 2024](https://arxiv.org/abs/2302.05329)。
- [IMRPhenomD](https://arxiv.org/abs/1508.07253) 和 [IMRPhenomXAS](https://arxiv.org/abs/2001.11412)：波形模型文献；具体调用 Ripple 实现。
- [GWOSC GW150914 v3](https://gwosc.org/eventapi/html/GWTC-1-confident/GW150914/v3/)：公开 H1/L1 应变，原始地址和校验值见 `data/manifest.json`。使用研究成果时应引用相应LVK论文与数据发布。
- [binaryBHexp](https://github.com/vijayvarma392/binaryBHexp)：MIT，借鉴轨道、自旋、波形的联动表达；未复制代码。
- [Cardiff Waveform Fitter](https://github.com/cardiffgravity/waveform-fitter)：借鉴滑块调参和观测叠加交互；未复制代码或缓存模板。
- [algometrix/blackhole-sim](https://github.com/algometrix/blackhole-sim)：MIT，Copyright 2026 algometrix。固定提交 `0ee5541ec4ffcd9ca708f37c080bf00b7f3ef6d1`。复制光线追踪、星图、HDR 渲染管线及其依赖闭包共 25 个文件（含 LICENSE），保存在 `renderer-src/vendor/blackhole-sim/`；`PROVENANCE.json` 记录逐文件原始 SHA-256。上游文件保持原样，新增适配器负责波形时间、质心坐标、暂停累积、相机、画质与回退。
- [Three.js](https://github.com/mrdoob/three.js)：0.185.1，MIT；OrbitControls 与渲染基础库。发布包附带许可证 `public/licenses/THIRD-PARTY.txt`。

本项目的 Canvas 轨道示意、界面与信号处理连接代码为独立实现；新增 WebGL 透镜部分明确复用上述 MIT 源码。仅使用本地资源，运行时不从第三方 CDN 加载脚本。Vite 和 GLSL 插件只用于构建，版本由 `renderer-src/package-lock.json` 锁定。
