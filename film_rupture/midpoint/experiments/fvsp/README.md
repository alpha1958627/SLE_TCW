# SLE-SMID-FVSP-BASILISK：FV面通量＋谱曲率

**这就是h1.2 TCW对照图的红线，是实际求解器的空间离散变体，不是后处理。** 黑线仍为[SLE-SMID-FV-BASILISK原FV核](../../sle_smid_fv_basilisk.c)。本实验源码：[sle_smid_fvsp_basilisk.c](sle_smid_fvsp_basilisk.c)。

`FVSP`固定表示“FV面通量＋谱能量曲率”。它沿用随机中点时间格式，同时在每个实际积分步的化学势、Jacobian和FFT预条件器中用谱q²替换原FV曲率符号lambda。迁移率、噪声面通量和散度仍是FV；不是完整谱求解器。两者后处理都对原场独立FFT，并与同一连续正则模比较；没有对黑线功率进行画图缩放。

| 项目 | 黑线FV | 红线FVSP |
|---|---|---|
| 化学势 | -lap_FV(h)+A/h³ | IFFT(q² FFT(h))+A/h³ |
| Jacobian曲率项 | -lap_FV(v) | IFFT(q² FFT(v)) |
| 常系数预条件符号 | 1+.5*dt*M*lambda*(lambda-3A) | 1+.5*dt*M*lambda*(q²-3A) |
| 迁移率/面噪声/散度 | FV | 同一FV |
| 中点/布朗树/拒绝/事件 | 同路径，无夹值或修补 | 同一协议 |

## 物理和数值边界

本试验h0=1.2、H=1、N64、cap=.3、phi_ref=.01、ls_ref=.2、hstar0，seed2026102201–205。A=.0013159/h0²、phi=.01/h0²、L=100/h0、ell=.2/h0。布朗根块.064/depth18，hmin<.3/.1降dt至.01/.003，事件hmin<=.01同路径细分。没有降低噪声、夹值、质量投影或重抽。新增bm_curvature()在推进场上运算；并非只改保存谱。

此变体改变了有限网格能量：梯度能改为.5*h*Q*h，其中Q的谱符号q²；其余势能和耗散保持。稳定线性模平衡谱因此更接近连续理论，但有限时间松弛仍含FV lambda，不能承诺长波/非线性平衡/首达或MFPT全面准确。旧FV的离散鞍点/前因子也不能直接当本模型精确EK。默认主程序仍使用原FV；本文件在experiments下单列。

## 固定t300证据

同五个seed、同参数，相对连续正则模，功率偏移为100*(P/P_NM−1)：

|实现|低段0<r<=8|中段8<r<=24|高段24<r<=46|
|---|---:|---:|---:|
|原FV|+4.7693%|+13.4953%|+71.5523%|
|FVSP|+3.3363%|−1.9547%|−.3876%|

原FV相对同N64离散参考的高段+.19345%是另一口径，不能称相对连续理论仅偏+.193%。[逐频段汇总](tcw_t300.csv)保留均值和五seed最小–最大，不是CI；配对实现不是10个独立seed。固定t300前缀也不是首次破碎结果。原FV数据来源2026-10-28 10:12审计，FVSP来源07:47审计；不当在线状态。原始场、库和全量日志留本地，不包含在Git。

实际选中平台为WUZH串行1CPU/条，h1.2五条；这是原FV正式40条之外的独立计算身份，不能混入原40统计。门槛见[gates.json](gates.json)：虽然历史文件名叫MPI_PASS，实际选中平台ranks=1，不能称本包重测过MPI16。四项192.36/47.99/34.00/12.23秒的CPU/GPU速度表属于原FV实现，不代表本FVSP的速度。

## 最小复用

在中点主目录执行；需兼容Basilisk/qcc、Python/NumPy/SciPy：

```sh
export BASILISK=/path/to/basilisk/src
CURVATURE=fvsp bash build_basilisk.sh
python runtime/run_single.py --kernel basilisk_fvsp --height 1.2 --cap .3 \
  --time 300 --seed 2026102201 --out runs/fvsp_t300
```

同一构建脚本不设CURVATURE仍编译原FV到basilisk_midpoint.so；本试验单独生成basilisk_fvsp.so，不覆盖原库。新增必要绑定在runtime/basilisk_fvsp_candidate.py，复用既有布朗树和有限窗入口，不新增绘图/SSH/调度器脚本。

检查点schema4含correction开关、源/库/运行版本指纹；与原FV schema3分开。此仓库的目录包装及共用Python指纹已变化，**不接管旧服务器检查点，不自动恢复任务**。C核与实际试验源逐字节相同，SHA见[provenance.json](provenance.json)。本次发布未重新进行qcc集群计算，门槛是有日期的既有证据。
