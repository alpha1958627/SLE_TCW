# SLE_TCW

随机润滑方程的物理研究与罕见事件采样：

**代码入口：[五条实现路线与标准命名](CODE_VERSIONS.md)**。先看[中点C/CUDA源码](film_rupture/midpoint/README.md)和[同16条耗时表](film_rupture/midpoint/BENCHMARK.md)：CPU1核192.36秒、CPU4核MPI47.99秒、A100批量34.00秒、CPU16核MPI12.23秒。核数是配置，不是算法版本。

| 目录 | 内容 |
|---|---|
| [tcw/](tcw/README.md) | 热毛细波谱、线性理论对照、绘图脚本、已有图表与 Basilisk TCW 求解器 |
| [film_rupture/](film_rupture/README.md) | 旧v3与分档版、L6噪声插值网格对照，以及Python/native/Basilisk保守随机中点直接MC；各批身份与验证范围分开 |
| [AMS/](AMS/README.md) | 2026-10-19重建的固定时间窗AMS；严格首次越级、完整状态分支、预算上限 |
| [DMC/](DMC/README.md) | 2026-10-19重建的DMC/QDMC；同物理核和评分，谱系及概率补偿 |
| [rare_events/](rare_events/README.md) | 两种采样框架共用的rupture参考物理核、构建、验收、独立直接模拟和有限试点 |
| [drop/](drop/README.md) | 液滴算例集合：钉扎/接触角与自由铺展（含可选涨落模式） |

TCW 原有内容整体移入 `tcw/`；在该目录内执行原 README 的命令，已有相对路径保持成立。绘图脚本中的外部原始数据目录仍需按本机数据位置设置。液膜破碎无排斥版v3已通过21组合1024步短测；长期、dt/网格及物理收敛仍未验收；两目录的滑移参数约定见各自源码与说明。
