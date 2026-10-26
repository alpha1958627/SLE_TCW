# 液膜破碎：无排斥版 v3 / ls等效0.2 / 21条轨迹

> **2026-10-26 版本交接**：本页下文是原21条的历史记录，不代表新版正式批次。[总版本表](../CODE_VERSIONS.md)区分空间索引/守恒修复与随机中点修正。新版 Python/native/Basilisk 源码及构建见 [midpoint/](midpoint/README.md)；固定 L6 噪声网格测试见 [l6_noise_bilinear/](l6_noise_bilinear/README.md)。新版全为直接MC，TCW验收不等于首次破碎时间收敛。

2026-10-19 按用户要求关闭 h^-5 排斥项，从相同初态和五个种子重新运行。保留 h^-3 吸引项；Pi(h)=-A/h^3，排斥系数 hstar=0。除这一物理项外，其余参数与守恒v2一致。此批是直接SLE模拟，不包含AMS/DMC等罕见事件加速采样。

> **2026-10-20 扩展**：[v3_tiered/](v3_tiered/README.md) 在同一物理（hstar=0）上加入
> **分档 dt cap**（hmin≥0.3→1e-1 / <0.3→1e-2 / <0.1→3e-3，依据 dt 上限测试）与
> **剖面/液面输出**（hmin 每降 0.01 一档 + 每 1000 t + 每 100 t 全场），
> 并把扫描扩展到 r=0.5–1.2（步长 0.1）× 5 seed = 40 条轨迹（SIPF 作业 29333）。

## 修复与验收

所有当前源码都包含：

1. 随机面通量的两个坐标哈希索引按 N=1<<LEVEL_MIN 做非负周期回绕，消除MPI周期两端随机通量不一致造成的真实质量漂移。
2. 体积统计只使用一次全局foreach(reduction)；移除随后重复MPI_Allreduce。当前日志绝对V不再多计8倍。
3. 达到步数上限时输出stop.hfield，便于独立体积核验。

不守恒的旧参考源码和旧生成流程已从当前项目移除。`source_yold_h1.05.c` 是用户要求的y索引实验对照：它与h1.05/yfix恰有五行y索引差异，**也包含上述两项守恒修复**。它不是旧MPI漏量程序；此实验对照不视作已验证生产模型。

先前15个分步测试中，MPI8纯噪声五种子的旧周期索引产生1.9e-5至1.85e-4量级的相对体积变化，周期回绕后≤7.1e-15。最终五份生产程序覆盖全部21个参数/种子组合，以MPI8、dt=0.003分别运行1024步：日志最大相对漂移5.788703e-13，最终保存场独立积分最大偏差6.008527e-13。[验收结果](validation/smoke_verification.json)、[分步测试](validation/periodic_noise_tests.csv)。短测通过不代表长期、临近破碎、dt/网格或物理验收通过。

## 参数与重复

r=1.05、1.10、1.15、1.20，各5 run；另加yold h1.05/run1，共21条。run1–5种子为2026101801–2026101805，跨厚度配对；对照与yfix run1同种子。

计算场初值h=1；A=0.0013159/r²、varphi=0.01/r²、L0=100/r、ell_code=0.2/r。迁移率M(h)=h³+ell_code h²；若采用h³+3b h²约定，则ell_code=3b，与TCW目录参数不能直接互换。

|r|A|varphi|L0|ell_code|yfix run数|
|---|---:|---:|---:|---:|---:|
|1.05|0.0011935600907029478|0.0090702947845804989|95.238095238095238|0.19047619047619048|5|
|1.10|0.0010875206611570248|0.0082644628099173554|90.909090909090909|0.18181818181818182|5|
|1.15|0.00099500945179584121|0.0075614366729678639|86.956521739130435|0.17391304347826087|5|
|1.20|0.00091381944444444444|0.0069444444444444444|83.333333333333333|0.16666666666666667|5|

level6（64²）、dt=0.003、上限2e8步、8 MPI ranks、8GB、14天墙钟；hstar=0（h^-5排斥关闭）、厚度截断下限0.005、hmin≤0.01停止。完整表见[runs.csv](runs.csv)，启动映射和源码SHA256见[tasks.tsv](tasks.tsv)。

## 构建、短测、正式运行

需要兼容的Basilisk/qcc、MPI、glutils和fb_tiny，先在SIPF环境核对分区、账号及MPI插件。以下命令从film_rupture目录执行：

```sh
python build_plan.py
python monitor_volume.py --self-check
export BASILISK=/path/to/compatible/basilisk/src
export CC99=/path/to/mpicc
export PATH="$(dirname "$CC99"):$PATH"
bash build_sipf.sh
sbatch --wait smoke.slurm
sbatch run_sipf.slurm
```

短测逐一使用全部21条源/种子组合，`check_smoke.py` 验证日志绝对V约为L0²及stop.hfield独立积分误差<1e-9，成功后写SMOKE_OK；正式启动器同时要求BUILD_OK和SMOKE_OK。`srun --input=none` 避免吞读批量任务表。每条目录的`.started`阻止覆盖已有轨迹；从新的空目录启动，不复制旧场或checkpoint续跑。

## 监测与证据边界

`monitor_volume.py` 仅读取已有log_restart、stderr和原始64²场。独立V=sum(h)*(L0/64)²，不归一化场或减均值；0.1%是警报线。厚度截断可能增质量，长期仍需观察。未结束轨迹保持右删失；调度COMPLETED或阈值事件不是物理验收。

看板说明见[dashboard/README.md](dashboard/README.md)。复制deployment.example.json并填写真实SSH及精确作业映射后，用refresh_dashboard.py只读采集；支持每站shell_prefix（如SIPF的环境初始化）。hmin曲线包含stderr中的精确事件终点，终点不补填缺失V。追加中的日志只接受精确白名单tar警告，其他采集失败保留上一成功看板。
