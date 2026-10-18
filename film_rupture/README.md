# 液膜破碎：ls 等效 0.2 / 21 条轨迹

2026-10-18 批次的源码、参数与监测工具。20 条 y 索引修正版轨迹（4 种厚度 × 5 seeds），另加 1 条旧版 h1.05/run1 同种子对照，共 **21 run**。这里的“罕见事件”是待观察的失稳事件；本批是直接 SLE 模拟，没有实现 AMS/DMC 等采样算法。

## 参数

令 r 为相对参考膜厚倍率，计算场初值统一 h=1：A=0.0013159/r²，varphi=0.01/r²，L0=100/r，ell_code=0.2/r。滑移迁移率为 M(h)=h³+ell_code h²；若用 Navier 参数 b 写作 h³+3b h²，则 ell_code=3b。此处参数定义与 `../tcw/solvers/` 的模型不能直接互换。

| r | A | varphi | L0 | ell_code | yfix runs |
|---|---|---|---|---|---|
| 1.05 | 0.0011935600907029478 | 0.0090702947845804989 | 95.238095238095238 | 0.19047619047619048 | 5 |
| 1.10 | 0.0010875206611570248 | 0.0082644628099173554 | 90.909090909090909 | 0.18181818181818182 | 5 |
| 1.15 | 0.00099500945179584121 | 0.0075614366729678639 | 86.956521739130435 | 0.17391304347826087 | 5 |
| 1.20 | 0.00091381944444444444 | 0.0069444444444444444 | 83.333333333333333 | 0.16666666666666667 | 5 |

各厚度 run1–5 分别用种子 **2026101801–2026101805**，跨厚度配对；task20 是旧版 r=1.05/run1/seed2026101801，匹配 task0。完整 21 行参数见 [runs.csv](runs.csv)，实际启动映射和源码 SHA-256 见 [tasks.tsv](tasks.tsv)。均为 level6（64×64）、dt=0.003、步数上限 200000000、8 MPI ranks、8GB、14 天 walltime。hstar=0.01，数值厚度下限0.005，事件停止阈值 hmin≤0.01。

## 新旧区别与证据边界

[reference/sipf_old_h105.c](reference/sipf_old_h105.c) 保留原始 SIPF 源码；[y_index_fix.diff](y_index_fix.diff) 只有五行改动：两处 `alpha.y[1]`、两处 `gamma.y[1]` 改成 `[0,1]`，显式 y 面中的 `f_old[-1]` 改成 `[0,-1]`。h1.05 的其他源码及参数字面量不变，其他厚度仅进一步替换 A、varphi、L0，滑移由命令行传入。参考索引修复来自 [SLE_TCW b5898c3](https://github.com/alpha1958627/SLE_TCW/commit/b5898c3a3a4eb8a350264c21dd004dd67a6d346e)，本批并未整体替换成 TCW 求解器。

**V 守恒尚未通过。** 本批原始 `volume_stats` 在 Basilisk reduction 后重复 MPI_Allreduce，8 ranks 的日志绝对 V 约为真实体积的8倍；相对首条日志值的漂移不受常数因子影响，但绝对 V 需独立积分原始场。2026-10-18 11:10（北京时间）的运行快照记录最大日志相对偏差约1.87%，最近保存场相对初始体积的最大偏差约1.68%。这些是该时刻观察，非当前实时状态。为保留新旧对照，此次发布不修复归约问题或改写在运行源码；不可据此声称体积收敛、时间步/网格收敛、破碎寿命或物理验收通过。

## 本地静态检查（不运行模拟）

```sh
python build_plan.py
python monitor_volume.py --self-check
```

第一条从冻结的旧源码重建四份 yfix 源码、21 行 runs.csv、tasks.tsv、补丁与 sources.sha256；`verification.json` 仅描述此静态检查。第二条用合成场验证二进制读取、全场积分和日志解析。

## 编译与提交（SIPF 模板）

要求兼容本源码的 Basilisk/qcc、MPI C 编译器及图形链接库 glutils/fb_tiny。具体 Basilisk 版本和 MPI 兼容性需在目标站点验证；源码发布不等于任意 qcc 版本均可编译。

```sh
cd film_rupture
python build_plan.py
export BASILISK=/path/to/compatible/basilisk/src
export CC99=/path/to/mpicc
export PATH="$(dirname "$CC99"):$PATH"
bash build_sipf.sh
# 核对 C64 分区、14天 walltime 和 pmix_v3 是否适合当前账号后提交。
sbatch run_sipf.slurm
```

数组task0–20；启动参数顺序是 `seed dt step_limit ell_code initial_h`。每个算例目录用 `.started` 防止意外重复启动；不要删除标记来覆盖已开始轨迹。`BUILD_OK` 和源码哈希检查只是启动条件，不是数值验收。若换集群，应先确认同一真实作业的排队/运行状态、编译器/MPI 与调度器接口，再适配分区和启动器；此模板不自动提交跨站副本或清理作业。

## 监测与看板

在已部署目录运行 `python monitor_volume.py > audit.json`，读取本地算例；可选 `SUBMITTED_JOB_ID` 只用于查询对应 Slurm 数组。独立体积为 `sum(h)*(L0/64)^2`，参考初始体积 L0²；不缩放场或补偿质量。0.1% 是监测警报线，不能替代 dt/网格收敛标准。有限运行或被停止的未破碎轨迹按右删失处理。

看板组件、只读采集和数据构建说明见 [dashboard/README.md](dashboard/README.md)。公开包不包含二进制、场文件、运行日志、实际 SSH 配置和账号目录。当前服务器计算及已部署看板继续使用原路径，仓库整理不触发重提计算。
