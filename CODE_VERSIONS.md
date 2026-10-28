# 液膜求解器版本总表

**按周期平膜破碎这条研究线，共2种时间格式、5条主要实现路线。** 1核、4核MPI、16核MPI是同一CPU实现的运行配置；h、dt、网格、噪声插值、phi和seed不是独立算法版本。仓库的液滴/TCW历史/AMS/DMC属于其他问题或框架，不加入这5条计数。

## 固定命名法

`SLE-时间格式-空间离散-计算后端`，精确源码版本由**Git提交SHA + 源SHA256**指定，不用“优化版、修复版、最新版”。

| 字段 | 取值与含义 |
|---|---|
| 时间格式 | `SPLIT`：历史分步；`SMID`：保守随机隐式中点 |
| 空间离散 | `FV`：周期均匀有限体积 |
| 后端 | `BASILISK`、`NUMPY`、`C99`、`CUDA` |
| 运行配置 | 单独记`N64 / cap0.3 / phi0.01 / MPI16 / B16`；不写进算法版本号 |

例如 `SLE-SMID-FV-C99 @ <commit> ; N64 ; cap0.3 ; phi0.01 ; MPI16 ; seeds2026102201-216` 唯一说明“算法/源版本/参数/资源”。`C99`指C算子与Python求解控制的混合实现，并非全C；CUDA保留`.cu`扩展名。

## 五条主线，从这里认代码

| 固定名称 | 核心源码 | 数值与并行职责 | 本项目中的用途 |
|---|---|---|---|
| **SLE-SPLIT-FV-BASILISK** | [旧分步C源码，含h1.2](film_rupture/v3_tiered/source_v3_h1.2.c) | 旧Basilisk分步更新、旧吸引项离散，.005厚度截断 | 最早cap=.1分档组；只作历史对照 |
| **SLE-SMID-FV-NUMPY** | [NumPy参考核心](film_rupture/midpoint/runtime/midpoint_candidate.py) | NumPy/SciPy算子、Newton/GMRES/FFT | 中点公式与实现参照；不是Basilisk |
| **SLE-SMID-FV-C99** | [C99面算子](film_rupture/midpoint/sle_smid_fv_c99.c) | C通量/噪声/Jacobian，Python/SciPy隐式解；CPU轨迹MPI | 原native批次及1/4/16核同工作量测速；同一内核 |
| **SLE-SMID-FV-BASILISK** | [Basilisk中点C源码](film_rupture/midpoint/sle_smid_fv_basilisk.c) | qcc周期FV，C Newton/GMRES/FFT，单场空间MPI | 正式40条所用中点核心；不能与旧分步程序混称 |
| **SLE-SMID-FV-CUDA** | [CUDA批量内核](film_rupture/midpoint/sle_smid_fv_cuda.cu) | FP64显存常驻多轨迹、融合Krylov；CPU兼容布朗路径 | A100批量16条34.00秒；仍是直接MC |

## 同16条轨迹总耗时

| 配置 | 总墙时 | 相对CPU1核 |
|---|---:|---:|
| C99，CPU1核 | 192.36秒 | 1.00× |
| C99，CPU4核MPI | 47.99秒 | 4.01× |
| CUDA，A100批量16条+1主机线程 | **34.00秒** | 5.66× |
| C99，CPU16核MPI | **12.23秒** | 15.72× |

同一g1节点、N64、h0=1.2、cap=.3、phi_ref=.01、同16seed；预热0→3，计时3→300，包括兼容噪声、上传和控制，末场写盘在计时外。一次测量/配置。**本基准CPU16核快于GPU约2.78倍**；不是GPU全面领先，也不是CPU单条轨迹需要192秒。完整条件、精确表与验收见[中点耗时说明](film_rupture/midpoint/BENCHMARK.md)。

## 这些区别不另造版本名

- `noise=L6-bilinear`：[L6噪声插值源码](film_rupture/l6_noise_bilinear/README.md)，仍属SPLIT-BASILISK，固定空间噪声截止；不是新的中点或白噪声网格极限。
- `y-index`历史实验：[y索引对照](film_rupture/source_yold_h1.05.c)与[周期守恒源](film_rupture/source_yfix_h1.05.c)，属旧格式的索引控制；不作为当前推荐的独立主线。
- `phi_ref=.013345`：[CPU统一参数入口](film_rupture/midpoint/runtime/README.md)，C99内核SHA完全相同；是物理参数扩展。统一入口默认phi=.01；旧冻结接口可按下方快照恢复，四项耗时仍只属于原版本。
- GPU早期逐步CuPy、未融合Krylov及4噪声线程探索属于性能试验历史；当前展示已验收的显存常驻融合Krylov源码，历史数据留本地档案，不重复命名为新算法。

## 算法差别与证据边界

SMID将漂移与噪声在`(h_new+h_old)/2`联合求解，以离散能量差构造化学势，面迁移率与噪声匹配。Newton物理残差<=1e-10，拒绝和事件沿同一布朗路径细分，没有厚度夹值、体积修补或噪声重抽。它和旧SPLIT同时改变时间及非线性空间离散，不能只称“dt调大”。

旧h1.2 cap=.1在实际t300.1的高段TCW相对同N64离散参考偏移约+178.113%。中点对同N64离散参考接近，但N64高波数对连续正则模仍有明显空间误差；参考不同必须明确。native h.9五seed均值cap.03=109315.767、.3=64126.063，差−41.34%，所以TCW接近也不能证明首次破碎/MFPT或dt收敛。

CPU/GPU配对16seed最大场差9.20e-14；GPU独立末步物理残差最大7.08e-12。基准到T300未破碎，不能把短窗性能验收当破碎统计。旧native135条与Basilisk40条复用路径，不能相加成175个独立seed。只有全组真实事件才给全组样本算术均值；其余保留删失/暂停/停止/失败。

图中时间保持`tau=h0*t_code/1e8`。EK全迁移率N49谱的鞍点首达和hmin阈值破碎不同；ell_m=.2对应Navier b/H*=1/15，不以理论填未破碎数据，不外推h<=1。

## GitHub展示与复用

中点目录优先展示三个[源码文件](film_rupture/midpoint/README.md)，运行必需的Python集中于runtime；绘图、部署、看板、SSH和一次性门槛脚本从这次液膜展示范围移出，Git历史仍可追溯。没有删除本地证据、改服务器运行源或改其他项目的Python脚本。[整理前完整目录](https://github.com/alpha1958627/SLE_TCW/tree/7e1d0959a22cc1a69ed65da2504914275ace0d2b/film_rupture)保留原路径与原启动历史。

文件重命名与Python路径绑定使新包指纹不同；**不要用新包恢复旧检查点**。C/CUDA内核逐字节保留原核，SHA见[清单](film_rupture/midpoint/source_manifest.json)。状态监测仍由各批独立审计负责，本页是代码目录，不宣称实时生产状态。

## 仓库整理记录

2026-10-28：日常分支统一为 `main`；独立实验只作为各项目的 `experiments/` 内容保留，不替换默认模型。四个研究入口为 `film_rupture/`、`drop/`、`tcw/`、`rare_events/`。

| 原文件/目录 | 当前入口 |
|---|---|
| `AMS/run.py`、`DMC/run.py` | `rare_events/ams.py`、`rare_events/dmc.py`；共用引擎与CLI，文档合入同目录README |
| `tools/` | `rare_events/check_fields.py`、`rare_events/plot_progress.py` |
| `runtime/phi_parameterized/` 四份重复Python | 共用 `runtime/`，CPU MPI入口增加 `--phi-ref`，默认.01；.013345需显式指定 |
| 两份TCW绘图脚本 | `tcw/plot_tcw_H0_reference.py`，纯驱动选择 `--pure-driven`；公式、时钟与输出名分别保持 |
| 日期实验中的4个submit/monitor脚本 | 原分支快照；当前保留参数、运行协议、独立数值验证 |
| 液膜旧部署/看板/源码生成脚本 | 已沿Notebook分支整理移出；见原main快照 |

整理前的完整分支固定为以下标签，保留精确源码、路径和历史：

| 原分支 | 快照标签 | 原HEAD |
|---|---|---|
| `main` | [archive/20261028/main](https://github.com/alpha1958627/SLE_TCW/tree/archive/20261028/main) | `cff0e3b` |
| `droplet-notebook-20261025` | [archive/20261028/droplet-notebook](https://github.com/alpha1958627/SLE_TCW/tree/archive/20261028/droplet-notebook) | `fa00372` |
| `codex/capillary-midpoint-20261027` | [archive/20261028/capillary-midpoint](https://github.com/alpha1958627/SLE_TCW/tree/archive/20261028/capillary-midpoint) | `cddbfcc` |
| `codex/q17-observer-20261021` | [archive/20261028/q17-observer](https://github.com/alpha1958627/SLE_TCW/tree/archive/20261028/q17-observer) | `78496ed` |
| `droplet-spreading-noise-modes-20261021` | [archive/20261028/droplet-spreading](https://github.com/alpha1958627/SLE_TCW/tree/archive/20261028/droplet-spreading) | `92f3eb6` |

例如 `git switch --detach archive/20261028/droplet-notebook` 可读回完整原版本。生产恢复仍使用原冻结目录、库及环境；新源码指纹拒绝接管旧检查点。C/CUDA物理核、原始结果和验收数值保持原字节；历史JSON只说明其对应源码版本。
