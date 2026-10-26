# 液膜代码版本交接（2026-10-26）

**图中的 `Basilisk .3` 是移植后的保守随机中点版**，代码在 [film_rupture/midpoint/](film_rupture/midpoint/README.md)。它确实使用 Basilisk/qcc 的网格、面通量、周期边界和 MPI；Python 负责布朗路径、步长控制和输出。它与最早 dt=.1 的分步 Basilisk 程序不同，也不等于“破碎时间已经收敛”。

## 一张版本表

| 版本及固定称呼 | 入口 | 实现与主要区别 | 适用范围/证据边界 |
|---|---|---|---|
| TCW Basilisk 基线 | [tcw/solvers/](tcw/solvers/) | 原 Basilisk 薄膜程序；修正 y 邻居索引 | 既有 TCW 对照；不是随机中点格式 |
| 守恒/yfix 无排斥 v3 | [film_rupture/](film_rupture/README.md) 的 `source_yfix_h*.c` | 周期噪声坐标回绕、体积归约和 y 索引修复；hstar=0；仍为旧时间格式 | 与 `source_yold_h1.05.c` 的 y 索引实验对照分开；短测守恒不证明统计正确 |
| **旧 v3 分档版（最早 dt=.1 那组）** | [film_rupture/v3_tiered/](film_rupture/v3_tiered/README.md) | hmin≥.3 用 cap .1，<.3 用 .01，<.1 用 .003；旧分步时间格式、旧非线性吸引项离散；有 .005 厚度截断 | 8 高度×5 seed；旧 h1.2 “t300”实际 t300.1 的高波数 TCW 偏移 +178.113%；不再称 dt=.1 安全/成核无偏 |
| L6 噪声双线性插值版 | [film_rupture/l6_noise_bilinear/](film_rupture/l6_noise_bilinear/README.md) | 仍是旧时间格式；把两面分量 L6 交错噪声周期插值到 L7/L8，固定 Delta6 归一化 | 固定空间噪声截止测试；90 条均达阈值，但 4 条触及 .005 截断；不能称网格收敛 |
| **Python 中点参考版** | [midpoint_candidate.py](film_rupture/midpoint/midpoint_candidate.py) + [adaptive_runner.py](film_rupture/midpoint/adaptive_runner.py) | NumPy/SciPy 有限体积、离散能量梯度、隐式随机中点；Newton/GMRES/FFT 均在 Python 栈 | 用于公式核验与实现参照；不是 Basilisk；拒绝步沿同一布朗路径细分 |
| **Python + native C 算子版** | [accelerated_runner.py](film_rupture/midpoint/accelerated_runner.py) + [native_midpoint_candidate.py](film_rupture/midpoint/native_midpoint_candidate.py) + [native_midpoint_kernel.c](film_rupture/midpoint/native_midpoint_kernel.c) | 只把 FV 漂移/噪声/Jacobian 算子放进 C；Newton/GMRES/FFT 仍由 Python/SciPy/NumPy 完成 | 旧 135 身份批次使用的独立求解器；不是 Basilisk，也不是“全 C 求解器” |
| **Basilisk 随机中点版（新版正式 40 条）** | [basilisk_midpoint.c](film_rupture/midpoint/basilisk_midpoint.c) + [basilisk_midpoint_candidate.py](film_rupture/midpoint/basilisk_midpoint_candidate.py) | qcc/Basilisk 均匀周期网格；Newton/GMRES 和 FFT 预条件器也在 C；Python 提供相同布朗增量 | N64、cap .3、8 高度×5 seed、16 MPI×1 线程；实现/TCW 门槛通过，FPT/MFPT 收敛尚未证明 |
| 液滴、钉扎与铺展 | [drop/](drop/README.md) | 另一类几何、初态与边界问题 | 保持原有说明和 notebook，不与周期平膜破碎直接混算 |
| 历史 AMS/DMC | [AMS/](AMS/README.md)、[DMC/](DMC/README.md)、[rare_events/](rare_events/README.md) | 仓库保留的其他采样框架 | **本次中点批次全部直接 MC，不调用分裂、重采样、权重方法** |

## 两层修正，不能混称

1. **空间索引/守恒修正**：y 邻居、周期面噪声一致性、避免体积重复全局归约。旧 v3 已包含其中相应修正；这不会自动消除时间离散造成的 TCW 畸变。
2. **随机中点修正**：漂移与噪声均在新旧场中点求解；吸引力改为离散能量的化学势差，面迁移率与噪声协方差匹配。拒绝步和阈值事件沿原布朗路径细分，不夹值、不补体积、不重抽噪声。它同时改变了时间格式和非线性空间离散，旧/新差异不能全部归因于 dt。

中点核心的七份源码保持部署时的字节；[provenance.json](film_rupture/midpoint/provenance.json) 记录 SHA256。新增的本地示例、构建脚本和测试不是服务器生产启动器。没有把库二进制、服务器配置、原场、检查点或未运行候选打进本包。

## 共同物理与比较口径

平膜初值 H=1，N=64，周期正方形；外部高度参数记为 h0：

```text
A = 0.0013159 / h0^2       phi = 0.01 / h0^2
L = 100 / h0              ell_code = 0.2 / h0
M(H) = H^3 + ell_code H^2
hstar = 0                破碎阈值 hmin <= 0.01
```

ell_code 是 h² 的系数；采用 M=H³+3bH² 时 ell_code=3b。理论基准 ell_m=.2 对应 Navier b/H*=1/15。不能把 .2 直接称为该约定下的 Navier b。

新中点生产布朗根块 .064、depth18；hmin<.3/.1 时将上限降到 .01/.003，失败与事件另作同路径细分。同 seed、同根块和同版本 RNG 的不同 dt/实现可共享路径；它们不是独立重复。跨高度复用 seed 标签，也不应把五个高度当成 25 个独立 seed。

时间图统一用 **tau = h0 × t_code / 1e8**；原始 t_code 保留。EK 参考是全迁移率谱 N49 的鞍点首达预测；这里的 hmin 阈值破碎不是同一事件，N49 谱与 N64 FV 也不同。不外推 h0≤1，不用理论填补未破碎轨迹。

## 已知结果与没有证明的事

- 新 Basilisk t300：8 高度、每组 5 seed，cap .3 高波数 TCW 组均偏移约 −0.021%～+1.695%；相对同路径 .03 的个体功率变化最大 4.330%。没有人为调小噪声来达到门槛。详见 [TCW 与速度表](film_rupture/midpoint/validation/README.md)。
- 实现门槛：原部署的算子/步进对照、续算、MPI 1/4/16/64 对照通过；这是实现验收。复制到本仓库后另做的本地测试有独立日期、运行环境和结果，不冒称重新完成集群验收。
- **TCW 接近不保证首次破碎时间接近**：native 中点 h.9 的 5-seed 样本算术均值，cap .03 为 109315.767，.3 为 64126.063（低 41.34%），.5 为 61214.054（低 44.00%）。不能因为 .3/.5 相近便宣称 dt 极限收敛。
- 本地冻结快照 2026-10-26 14:04:54：Basilisk 40 身份为 33 阈值事件+7 运行；h.6/.7/.8/.9/1/1.1 各5/5，h1.2 为3/5，h1.3 为0/5。**这不是在线状态服务**，后续以新的真实审计为准。
- 旧 native 135 身份已终止：81 事件、29 固定窗右删失、25 用户停止；另 4 条旧参考为 1 事件+3 墙时软暂停，单列。Basilisk/native 复用路径，不能把两实现身份相加当独立 MC 样本。
- 只有组内五条都真破碎才报五条样本算术均值；运行、固定窗删失、用户停止、软暂停及数值失败分别保留。调度 COMPLETED 不等于破碎；五条均值也不等于收敛的总体 MFPT。

## 从哪里开始复用

先看 [midpoint/README.md](film_rupture/midpoint/README.md) 的文件表、构建和有限时间示例；纯 Python 参考不需要编译器。旧 dt=.1 源码保留在 v3_tiered，L6 插值源单列。现有 droplet notebook 未改动；本次不发布生产原场和完整轨迹数据。
