# 保守随机中点：Python、native C 和 Basilisk

这是同一离散模型的三种实现；[总版本表](../../CODE_VERSIONS.md) 解释它们与旧 dt=.1、yfix 和 L6 噪声插值版的关系。核心源码原样复制，SHA256 见 [provenance.json](provenance.json)。本目录不包含生产原场、检查点、编译库或服务器连接配置。

## 文件与职责

| 文件 | 职责 |
|---|---|
| `midpoint_candidate.py` | 纯 NumPy/SciPy FV 漂移、噪声、解析 Jacobian、Newton/GMRES、FFT 预条件；`--self-test` 验证离散结构 |
| `adaptive_runner.py` | BrownianTape、同路径拒绝/事件细分、hmin 分档、SHA/运行版本检查点；历史 schema1 默认根块 .001/depth12 |
| `accelerated_runner.py` | reference/native 统一轨迹 schema2；可显式选生产根块 .064/depth18 |
| `native_midpoint_kernel.c` | FV 力、面噪声与 Jacobian 算子 C 加速；不负责 Newton/GMRES |
| `native_midpoint_candidate.py` | ctypes 绑定、GCC 构建和编译清单；Newton/GMRES 仍用 Python 参考实现 |
| `basilisk_midpoint.c` | Basilisk 多重网格、周期面算子、Newton/GMRES 与自带 radix-2 FFT；可编译 MPI |
| `basilisk_midpoint_candidate.py` | qcc 库绑定与 schema3；Python 驱动相同布朗路径 |
| `example_run.py` | 新增的有限时间、本地直接 MC 示例；保留阈值事件或右删失；不是正式批量启动器 |
| `build_basilisk.sh` | 新增的本地 qcc 构建命令，默认串行，`MPI=1` 编译 MPI 库 |
| `test_replay.py` | 新增的参考/native、布朗可加性、逐字节续算和错误 seed 拒收检查 |
| `validate_and_benchmark.py` | 原集群实现门槛和同节点固定窗测速工具；使用 `compiled/`，测速需 Linux CPU affinity |
| `build_and_smoke.py` | 原 qcc/MPI 编译与算子/单步门槛；依赖实际分配节点和兼容环境 |
| `long_window_pilot.py` | 原 reference/native 固定窗驱动，定时原场、检查点、显式恢复与软停 |
| `run_formal.py` | 原 Basilisk 事件终止驱动；依赖本机重新验收后的 case card，默认不提供生产选择配置 |
| `tcw_probe.py` | 原 t300 配对 .03/.3 TCW 工具；读取 task index，详见下文 |
| `validation/` | 脱离主机配置的组统计与验证范围；不是完整原场复现包 |

## 离散模型与改变

```text
mu(h) = -lap_FD(h) + A / h^3
M_face = h_face^3 + ell_code * h_face^2
h_face = (h_left + h_right) / 2
F(h) = -K(h) mu(h)
h_new - h_old = dt F((h_new+h_old)/2) + G((h_new+h_old)/2) dW
```

两个方向使用同一面迁移率的平方根噪声，协方差满足有限维 FDT。面通量周期相消；没有厚度夹值/体积修补。Newton/GMRES 残差容差 1e-10，正性线搜索；不正/不收敛则拒绝并细分原布朗增量。事件阈值 hmin≤.01，事件步细分至约 1e-5 代码时间。源码中 `legacy_chain` 只是旧吸引项链式离散诊断，不是生产默认值。

Basilisk 此处解的是标量薄膜方程，**没有使用 VOF 两相 Navier–Stokes**；固定 N64，无 AMR。MPI 面算子在本地单元计算并交换边界，4096 维向量与 FFT 预条件器在各 rank 复制，不能推断无限加核会加速。

## 安装与最小运行

下列命令均在本目录执行；Python 3.6.8/NumPy1.19.5/SciPy1.5.4 是原集群冻结组合。现代 Python 的独立复用可安装 `numpy`、`scipy`，但随机数/检查点不保证跨版本相同；不要把新版环境直接当成旧轨迹续算环境。

```sh
python -m pip install -r requirements.txt
python midpoint_candidate.py --self-test
python example_run.py --kernel reference --height 1.1 --time 3 --out runs/reference
```

**native C 算子版**需要 GCC；构建目录必须尚不存在。动态库与 build.json 一同保留：

```sh
python -c "from native_midpoint_candidate import build; print(build('compiled'))"
# Linux: compiled/kernel.so；Windows: compiled/kernel.dll
python example_run.py --kernel native --library compiled/kernel.so --height 1.1 --time 30 --out runs/native
python test_replay.py --out validation/local_replay.json
```

**Basilisk**需要已有兼容的 Basilisk/qcc、C 编译器；MPI 版另需匹配的 MPI 编译器/运行库，不需要 mpi4py。构建脚本不会覆盖已有库。qcc 使用相对源文件名，并关闭 fast-math/浮点乘加收缩：

```sh
export BASILISK=/path/to/basilisk/src
bash build_basilisk.sh
python example_run.py --kernel basilisk --height 1.1 --time 30 --out runs/basilisk

# 在另一个干净的代码副本中构建 MPI 版：
export CC99=mpicc
MPI=1 bash build_basilisk.sh
mpiexec -n 16 python example_run.py --kernel basilisk --height 1.1 --time 30 --out runs/basilisk_mpi
```

同一副本一次只使用一个 Basilisk 库/活动网格配置。MPI 的每个 rank 都运行 Python 控制器和同一路径；仅 rank0 写文件。`example_run.py` 使用 N64、cap .3、根块 .064/depth18，默认到 t=30 即停止；输出 `RIGHT_CENSORED_AT_FIXED_HORIZON` 表示随访到该时间仍未破碎，不是破碎完成。输出目录必须为空且不存在，防止覆盖。

示例只写最后检查点、hmin 轨迹和 summary，不替代原生产的每300时间原场/每1000步检查点策略。恢复 API：`accelerated_runner.Trajectory.load(path, library=...)` 或 `basilisk_midpoint_candidate.Trajectory.load(path)`；后者记录原绝对库路径和 MPI 核数。同源、同实际库、同运行环境才能续算。schema1/2/3 不混读；没有跨平台移植检查点的默认保证。

## 原长窗与正式驱动

`long_window_pilot.py` 保留原 native 批次的固定窗/软停协议，可在新目录运行：

```sh
python long_window_pilot.py --out runs/window --height 1.1 --seed 2026102201 --caps .03 .3 --horizon 300 --interval 300 --native-library compiled/kernel.so
```

有 `--native-library` 时显式采用根块 .064/depth18；**没有该选项的历史 reference CLI 仍默认 schema1/.001/depth12**。要做共享路径的 reference/native 对照，使用上文 example_run 或显式 `accelerated_runner.Trajectory(kernel="reference", base_dt=.064, depth=18)`，不要混用不同布朗协议。

`run_formal.py` 是已部署的原 Basilisk 驱动，使用 `SLURM_ARRAY_TASK_ID`、`FORMAL_SELECTION.json` 的 PASS/cases，以及本机 TCW_PASS.json；它每300时间输出原场，每1000步保存状态，事件终止且无固定模拟终点。**已有 formal/task_XX/state.npz 时该驱动自动尝试同指纹恢复**，并非总是新起轨迹。仓库不提供实际生产选择文件、调度器提交器或继续任务脚本，不能把 validation/ 的历史 PASS 复制为新平台验收；先在独立目录重新构建、验收、绑定源/库/案例再接入自己的调度器。源码保存是为了复核原生产行为，最小复用从有限时间示例开始。

## 复用原验收工具

构建 `compiled/kernel.so` 和 `compiled/basilisk_midpoint.so` 后：

```sh
python validate_and_benchmark.py gate
mpiexec -n 16 python validate_and_benchmark.py mpi
python validate_and_benchmark.py benchmark --solver native --height 1.1 --repeat 0
mpiexec -n 16 python validate_and_benchmark.py benchmark --solver basilisk --height 1.1 --repeat 0
```

同节点、同核数、同绑定才用于实现速度比较。`benchmark` 固定 seed、暖机至 t3、测速 t3→33、不含异步输出；比较 1CPU/native 与 16CPU/Basilisk 时必须同时报告核数，不能称同资源加速。

`tcw_probe.py` 的位置参数为 0–39（8 高度依次 .6/.7/.8/.9/1/1.1/1.2/1.3，每高度5seed），使用串行 qcc 库；每身份沿共同布朗路径跑 cap .03 和 .3 到 t300。它会写原场，需要另行空间；结果不自动提交 Git。

独立原型 `midpoint_candidate.py --out ...` 的 probe 使用 .001 根时间网格和 .15 宏块，输出标记 `restartable=False`。它与生产 .064/depth18 轨迹不同，不能混作额外独立重复。

## 验证边界

见 [validation/README.md](validation/README.md)。TCW 是同 N64 的有限时间线性 FD 参考；高段 24<sqrt(mx²+my²)≤46。固定 t300 与不同末时刻动态谱不能混称平衡测度误差。门槛通过不证明首次破碎、总体 MFPT、网格或 dt 收敛。

本目录仅直接 MC。状态类里保留的 `fork()` 工具不是这个示例的采样步骤；示例/正式中点批次不调用它、不分裂/重采样/加权。
