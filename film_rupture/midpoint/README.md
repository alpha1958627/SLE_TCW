# 随机中点源码：先看C / CUDA

**四条中点实现，同一时间格式；MPI核数和dt是运行配置。** 全部直接Monte Carlo。总版本表见[CODE_VERSIONS.md](../../CODE_VERSIONS.md)，同16条速度见[耗时表](BENCHMARK.md)。旧分步格式保留在[历史源码](../v3_tiered/)。

| 稳定名称 | 优先阅读源码 | 谁完成求解 | 并行方式 |
|---|---|---|---|
| SLE-SMID-FV-NUMPY | [运行必需的参考核心](runtime/midpoint_candidate.py) | NumPy/SciPy算子、Newton/GMRES/FFT | 单轨迹；也可由MPI分发多轨迹 |
| **SLE-SMID-FV-C99** | [sle_smid_fv_c99.c](sle_smid_fv_c99.c) | C算通量/Jacobian；Python/SciPy解Newton/GMRES/FFT | CPU轨迹MPI，1/4/16核共用这一个C内核 |
| **SLE-SMID-FV-BASILISK** | [sle_smid_fv_basilisk.c](sle_smid_fv_basilisk.c) | Basilisk面算子，C Newton/GMRES及FFT | 单场空间MPI；小网格向量/FFT各rank复制 |
| **SLE-SMID-FV-CUDA** | [sle_smid_fv_cuda.cu](sle_smid_fv_cuda.cu) | 显存常驻批量FP64算子、Krylov正交化/旋转/回代 | 单GPU多轨迹批量；CPU提供兼容布朗增量 |

`.cu`是CUDA源码，保留真实扩展名。三个C/CUDA内核均与其已验收来源**逐字节相同**，仅换描述性文件名。这里没有把Python驱动伪装成C，也没有重写或修改正在运行的服务器源码。

## 为什么仍有runtime目录

[runtime/](runtime/README.md)只放运行必需的控制器、绑定与布朗树：它们负责随机路径、拒绝步/事件细分、检查点、MPI分工和CUDA调度。去掉后C文件不能独立复现这些轨迹。绘图、看板、SSH、竞跑、部署和历史门槛脚本已从这一展示目录移除，旧文件仍在[整理前提交](https://github.com/alpha1958627/SLE_TCW/tree/7e1d0959a22cc1a69ed65da2504914275ace0d2b/film_rupture/midpoint)。

## 物理和格式

周期二维均匀FV，`mu=-lap(h)+A/h^3`，面迁移率`M=h_face^3+ell*h_face^2`，噪声用同一`M`的平方根。漂移和噪声在`(h_new+h_old)/2`共同求解。Newton物理残差默认<=1e-10；不夹值、不修补体积、不重抽噪声。

固定phi接口：`A=.0013159/h0², phi=.01/h0², L=100/h0, ell=.2/h0, H=1, hstar=0`。布朗根块.064/depth18，hmin<.3/.1将dt上限降至.01/.003，事件hmin<=.01沿同一路径细分。`ell=.2`是迁移率h²系数，采用`h³+3b*h²`约定时对应Navier `b/H*=1/15`。

`phi_ref=.013345`扩展仍用同一C99内核，驱动独立放在[runtime/phi_parameterized/](runtime/phi_parameterized/README.md)，不能把该物理参数下的轨迹当成上面phi=.01测速或旧135条的重复。

## 最小运行

从本目录执行；Linux CPU MPI需GCC、MPI、mpi4py、NumPy/SciPy。先限制每rank线程数：

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# CPU单条有限窗（纯NumPy/SciPy参考）
python runtime/run_single.py --kernel reference --height 1.2 --time 3 --out runs/reference
# CPU轨迹MPI；首次由rank0编译上面的C99源码
mpiexec -n 4 python runtime/run_cpu_mpi.py --heights 1.2 \
  --seeds 2026102201 2026102202 2026102203 2026102204 --until 300 --out runs/cpu4
# GPU批量；安装与本机CUDA匹配的CuPy，显式选择可用显卡
CUDA_VISIBLE_DEVICES=0 python runtime/run_gpu_batch.py --heights 1.2 \
  --seeds 2026102201 2026102202 2026102203 2026102204 --until 300 --out runs/gpu4
```

Basilisk需兼容qcc；串行或MPI库分别在干净副本构建：

```sh
export BASILISK=/path/to/basilisk/src
bash build_basilisk.sh
python runtime/run_single.py --kernel basilisk --height 1.2 --time 3 --out runs/basilisk
# 空编译目录内：export CC99=mpicc; MPI=1 bash build_basilisk.sh
# mpiexec -n 16 python runtime/run_single.py --kernel basilisk ...
```

这里的命令是功能示例，**不是上述16条预热测速程序**，也不自动向集群提交任务。GPU与CPU MPI完整入口使用`--until 0`等事件，墙时到限保存，`--resume`显式恢复；有限窗`CENSORED`不等于破碎。单条示例只提供有限窗。

源码搬移后路径和Python指纹改变：**新目录不能接管旧生产检查点**。旧任务继续使用原冻结目录/库/运行环境。原核SHA、来源与路径整理记录见[source_manifest.json](source_manifest.json)；本次本地包装检查见[packaging_checks.json](packaging_checks.json)。

## 使用结论

这个N64、16条基准优先选C99轨迹MPI，16核12.23秒；GPU批量34.00秒是可用实现，仍受CPU兼容噪声约束。Basilisk适合核对空间MPI路径；不能把加核承诺成无限加速。短窗一致性、TCW以及5seed事件均值均不等于首次破碎/MFPT或dt/网格收敛。
