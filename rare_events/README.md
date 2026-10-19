# rupture基准的共用罕见事件实现

此实现于2026-10-19从仓库 `fed7bb3bf963afe44466017e61199c61613b417d` 的 `film_rupture/source_yfix_h1.05.c` 重建。本线程旧AMS/DMC代码、场、种子和权重不导入。原 `tcw/`、`film_rupture/` 保持不变。

发布时合入上游 `77ce5cb` 的无排斥v3更新：它仅将参考源码的hstar置零。生成器兼容这两份精确SHA，显式消去零项后，`engine.c`、`physics.h` 和完整参考程序的字节与已验收/正在运行版本保持一致，物理源码签名不变。兼容性证据见 `validation/upstream_compatibility_20261019.json`；不为这一生成器兼容更新重复计算已有路径。今后克隆构建仍须自行跑数值门槛，不能把本仓库的历史验收JSON当作本机门槛。

## 物理与差异

用户已确认保留ratio=1.01、ls_ref=.2、无排斥。初始归一化平膜h=1；`mu=-lap(h)+A/h^3`、`M=h^3+ell_code*h^2`。

|量|值/约定|
|---|---|
|A、phi|.0013159/ratio²、.01/ratio²|
|周期域|L=100/ratio，64×64均匀网格|
|滑移|ell_code=.2/ratio；若写M=h³+3b h²，则ell_code=3b|
|时间步|固定请求dt=.003；半dt=.0015不会爬回.003|
|事件与数值下限|hmin<=.01破碎；.005数值截断另记录原始最小值和增质量|
|时间|已完成步数×dt；乘ratio为共同参考，无ns映射|

`prepare_sources.py`核对原文件SHA，提取隐式求解器与周期守恒面噪声，保留修复的y索引。按用户要求只删除beta/gamma中的h^-5排斥贡献，完整补丁在 `no_repulsion.diff`。核的生成可用 `--check` 字节复验。均匀网格使用multigrid后端、Jacobi迭代，MPI只用1/4/16等4的幂；这些数值选择也在provenance中标明。

每次更新先复制旧场，再求隐式确定性演化，施加以旧场迁移率计算的周期随机面通量，记录原始最小值并应用原下限，最后判断事件、评分与守恒。首次tau是离散检测的完成步时刻，所在步区间可由dt给出；不伪造连续时间精确穿越。超过质量漂移/高度门槛时保存失败全场并非零退出，不给物理结果。

## 构建和验收

兼容Python3.6及Basilisk/qcc。需要对应MPI与Basilisk库；不依赖旧目录或旧二进制。从仓库根目录：

```sh
python3 rare_events/prepare_sources.py --check
python3 -m rare_events.test_estimators
export BASILISK=/path/to/compatible/basilisk/src
export CC99=mpicc
bash rare_events/build.sh
sbatch rare_events/validate.slurm
```

验收包括：

- 两状态两步过程完整枚举：真值P=.0494、RMST=1.98、随机质量期望1；演示不补偿的克隆比例有偏。另以解析Uniform条件采样检验AMS，p=.05、2000独立重复。
- 从完整 `film_rupture` 程序生成独立事件驱动参考，作相同无排斥/固定dt/multigrid/Jacobi修改，与分段引擎比较2000次真实更新。参考原事件i=0已更新一次，故终止i=1999对应2000次更新；新引擎记录完成步时钟，避免原标签偏一拍。
- 两种子真串行/MPI1/4/16、重复、观察非干扰、1000+1000重启、独立未来噪声和质量漂移；执行真实半dt并拒绝外国/旧格式检查点。
- 当前控制器在正式试点启动前再运行独立小烟测，含AMS、DMC、QDMC及直接MC。烟测文件放在 `_runtime_gate`，不混入试点统计。

检查点 `RESTFE26` 明确保存版本、64字符物理源码签名、网格、已完成步、源种子、dt/ratio/slip/phi/A/L/事件/下限、初始质量、历史评分和累计截断增质量，以及64²全场；Python与C都拒绝参数/签名不符的文件。旧AMS2/3不能恢复。

数值门槛只覆盖所列短时操作。通过不代表长期、临近破碎、半dt收敛、空间/噪声截断或物理寿命已经验收。原始V由场积分得到，不减均值、不归一化掩盖漂移。

## 有限预算试点与监测

验收通过后（SIPF示例环境脚本在slurm文件中；其他站须适配环境和分区）：

```sh
export VALIDATED_GATE=validation/gate_<verified_job>
sbatch --export=ALL,VALIDATED_GATE rare_events/pilot.slurm
python3 -m rare_events.monitor --workdir runs/job_<actual_job> --job <actual_job>
```

首轮各N12、T3000、每组3条MPI4并行：AMS D1最多16轮；DMC/QDMC每150同步、末端强度4/.5线性爬升；独立直接MC也12条。四组共48核，初始48种子各不相同，namespace含真实作业号和方法。采用共同q3，只是有预算的初始可行性试点，不是完成时间窗选择或优化后的效率对照。

完成后必须检查事件、场结构、评分的持续性、家系和成本，再决定独立新种子/新窗口/评分试验和64条直接MC及半dt验证；不能无限分裂至极小W，也不能止于把试点当完整寿命。跨重复比较需固定物理、评分和窗口；三方法选择强度不是等效或最优的已知值。所有MFPT字段保持null，直到尾部和物理验收支持结论。

日志/检查点/执行清单、每轮父索引和权重保存于显式workdir。完成求解核时=wall_seconds×实际MPI秩数；活动调用经过时间估计另列，Slurm分配核时还包括环境门槛、控制器和空闲。运行/下载观察出错时不自动重投；未确认的物理调用保留待审查。

GitHub不保存运行场、二进制、实际服务器凭据或可执行文件。验收结果应另保存带日期的JSON与源码/二进制SHA，并分别报告调度、数值、统计和物理状态。
