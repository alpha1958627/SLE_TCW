# DMC / QDMC（2026-10-19重建）

这是按时间同步的带权重复制/淘汰算法，物理核、破碎事件和q3评分与 [AMS](../AMS/README.md) 相同，不修改薄膜动力学。构建和验收见 [共用说明](../rare_events/README.md)。

从仓库根目录运行：

```sh
python3 -m DMC.run --mode dmc --workdir runs/my_dmc --namespace my_new_campaign:dmc:R0 \
  --N 12 --horizon 3000 --block 150 --strength 4 --ranks 4 --workers 3
python3 -m DMC.run --mode qdmc --workdir runs/my_qdmc --namespace my_new_campaign:qdmc:R0 \
  --N 12 --horizon 3000 --block 150 --strength .5 --ranks 4 --workers 3
```

每一步：

1. 从平膜启动N个独立种子，所有粒子前进到同一个物理时刻；已经破碎的粒子吸收，首次tau不覆盖。
2. DMC以同一个q3构造势 `V_k=c*(t_k/T)*(1-q3)`。QDMC先用逆上一轮倾斜的权重 `exp(-V_old)` 建立加权中点CDF，再作逆正态映射；不是无权排名，也不是Z-score。
3. 增量权重 `G_i=exp(V_new_i-V_old_i)`，累计 `logZ += log(mean(G))`，按归一化G系统重采样。一个父粒子分出多个孩子时，未来噪声种子全部独立。
4. 保存父索引、祖先家系、势、logZ、选择ESS、完整检查点和首次tau。继承的是被选父粒子的势与tau，不能丢失之前的选择补偿。
5. 后选择粒子贡献 `w_i=exp(logZ-V_i)/N`。计算 `P(T)=sum(w_i*I(tau_i<=T))`、`RMST(T)=T-sum(w_i*(T-tau_i)+)`；17点CDF使用同样权重。

`sum(w)`是期望为1的随机质量诊断，不是物理体积。**不再除sum(w)、不截断越界P、不剔除不利重复。** 家系坍缩、选择ESS和补偿权重波动都要监测；粒子ESS不等于独立重复数。

代码中的强度是探索参数，普通DMC和QDMC的数值不能解释成等效选择强度；试点强度4/.5不是已优化参数。每组完整重复才可作为误差样本。局部时间乘ratio转换到共同参考，无ns映射。完整MFPT必须另验尾部、步长和空间噪声截断。
