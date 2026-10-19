# 固定时间窗 AMS（2026-10-19重建）

物理核来自本仓库当前修复的 `film_rupture`，与DMC/QDMC、直接模拟共用。本线程已确认的条件是ratio=1.01、ls_ref=.2、无h^-5排斥；原 `film_rupture/` 未修改。源码来源、物理差异和数值门槛见 [共用说明](../rare_events/README.md)。

从仓库根目录运行（先构建并完成数值门槛）：

```sh
python3 -m AMS.run --workdir runs/my_ams --namespace my_new_independent_campaign:ams:R0 \
  --N 12 --D 1 --max-rounds 16 --horizon 3000 --ratio 1.01 --slip .2 --dt .003 --ranks 4 --workers 3
```

每一步：

1. 从平膜启动N个独立噪声种子，保存全场、完成步数、物理参数、源码签名以及沿路更高分状态。
2. 统一评分为 `progress=1-q3`，`q3=hmin+clip((hmin-.01)/.19,0,1)*(global_min_3x3_mean-hmin)`。它只决定采样，不改变力；q3<=.01与原破碎事件一致。
3. 淘汰历史最高分最差的D条；并列分数一起淘汰，实际淘汰数进入 `logW += log(1-k/N)`。
4. 从幸存路径**第一次严格超过淘汰等级**的完整检查点分支，保留之前路径和绝对时钟，之后改用独立噪声。不会从幸存者未来终点倒退拼接。
5. 目标等级达到后，`P(T)=W*事件比例`，CDF及RMST由相同路径权重和首次tau计算。预算用完则 `P/RMST=null`，只保留logW；不把它当概率或概率下界。

`state.json`、逐轮 `selection_*.json`、完整检查点、轨迹和原始质量漂移都保留。源/参数发生变化就拒绝续算；只有确认原控制器停止后才显式 `--resume`。未完成求解调用不会自动重放。

AMS路径/克隆不是独立误差样本。效率比较需要独立全程重复、相同物理/评分/窗口及方差×实际求解核时；全零估计不能拿零样本方差证明效率。有限窗CDF/RMST不是完整MFPT。
