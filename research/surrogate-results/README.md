# NRSur7dq4 Windows 运行与核验

2026-10-23，实际运行成功。`dynamics.json` 为可视化数据，`nrsur7dq4-full.npz` 保留 21 个球谐模式和全部模型动力学。它们均来自 GWSurrogate 1.2.0 与已校验的 NRSur7dq4 模型文件，不是自行构造的轨道。

## 运行环境与精确改动

主 `simulator/.venv` 未改变。运行使用独立 `simulator/.adapters-venv`，从官方 1.2.0 源码包重新解压至忽略目录 `research/adapter-build/zig-runtime/`，通过局部 `sys.path` 导入；没有宣称完成官方 Windows wheel 安装。

便携 Zig 0.14.1 从 [官方发布目录](https://ziglang.org/download/0.14.1/zig-x86_64-windows-0.14.1.zip) 下载，SHA-256 为 `554f5378228923ffd558eac35e21af020c73789d87afeabf4bfd16f2e6feed2c`。使用其 x86_64-windows-gnu 目标编译 C11/C++14，支持复数和变长数组。源码公式保持不变，三项兼容调整是：

- Windows LLP64 的 `long` 为 32 位，而模型 `bfOrders` 数据为 int64。五处阶数数组指针明确改为 `npy_int64`，补丁见 `windows-int64.patch`。
- 缺失的 `CMPLX` 宏由 Clang 的 `__builtin_complex` 构造器补齐。
- 三个供 ctypes 调用的样条函数通过 `.def` 文件显式导出。

`runtime-build.json` 记录命令、源码和二进制哈希。该适配运行环境只验证了本次 NRSur7dq4 路径；不代表所有 GWSurrogate 模型均已兼容 Windows。

## 数据与物理约定

输入 q=2，chiA=(0.4,0.1,0.2)，chiB=(-0.1,0.2,-0.1)，位于 q≤4、|chi|≤0.8 的训练范围内。`f_low=0` 保留完整波形，`f_ref=0` 是使用模型初始时刻的特殊标记，不能解释为物理零频率。初始时刻约 -4300 M，dt=2 M，共 2201 帧至约 +100 M。M 是总质量的几何时间尺度；没有选择物理总质量或距离。

参考惯性系：初始轨道角动量沿 +z，从轻黑洞到重黑洞的方向沿 +x，y 完成右手系。`q_copr` 是 scalar-first [w,x,y,z]，以 `q v q^-1` 将共进动系向量旋转到惯性系。chiA/chiB 为惯性系无量纲自旋；orbphase 是共进动系轨道相位。h22 是惯性系球谐模式的复数值，单位 r h/M，不是任意观测者或探测器的应变。

模型 t=0 采用多模式总振幅峰值时间约定，不必等于 h22 模式峰值。它与 SXS 共同视界形成时刻、坐标时间或无穷远波形迟滞时间均未对齐。返回量没有真实位置或分离度，JSON 明确标记缺失，不生成虚构轨道。

[模型论文](https://arxiv.org/abs/1905.09300) 第 III C 节说明，自 t/M=-100 起，训练自旋通过 PN 方程延伸，延伸自旋不具物理意义，仅作为晚期拟合参数。因此 `dynamics_valid_until_M=-100`；可视化方向/自旋面板保守地在此截止。波形可继续显示并合及振铃，不能将晚期双自旋箭头解释成两个仍独立存在的黑洞。该截止点也不是共同视界形成的测量值。

## 检查结果与范围

`validation.json` 保存实测量。全部 2345 组模型基函数拟合独立用 NumPy 多项式求和对照，最大绝对差 2.78×10^-16；四元数范数、自旋范数及与 SciPy Rotation 比较的最大误差均低于 3.34×10^-16。以显式 times 重复采样得到完全相同的对应模式值。初始自旋与给定值一致，所有输出有限。

官方源码中的六组测试文件共 **79 项通过**，其中自然样条独立对照 SciPy，Wigner 矩阵有 Python 参考和幺正性检验，固定波形值回归采用上游硬编码预期值。测试命令及结果见 `upstream-tests.json`。这些检查支持数值实现和适配一致性，尚不是与独立 NR 数据的波形失配评估，也未得到新的形成渠道后验结论。导入时的 `cannot import LAL` 是可选 LAL 模块未安装的提示，本次运行和上述测试不依赖 LAL。

## 复现

项目根目录 PowerShell：

```powershell
simulator/.adapters-venv/Scripts/python.exe -m pip install -r research/surrogate-results/requirements-lock.txt
simulator/.adapters-venv/Scripts/python.exe research/fetch_surrogate_compiler.py
simulator/.adapters-venv/Scripts/python.exe research/build_surrogate_runtime.py
simulator/.adapters-venv/Scripts/python.exe -X utf8 research/run_surrogate.py
```

模型仍位于 `research/surrogate-data/NRSur7dq4.h5`，来源、许可、远程校验见其 manifest。重新运行上游固定回归时，将此模型复制到运行源码的 `gwsurrogate/surrogate_downloads/NRSur7dq4.h5`（测试固定查找此处），不需要重新下载其他模型。

