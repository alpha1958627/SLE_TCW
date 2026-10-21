# SLE_TCW

随机润滑方程的物理研究与罕见事件采样：

| 目录 | 内容 |
|---|---|
| [tcw/](tcw/README.md) | 热毛细波谱、线性理论对照、绘图脚本、已有图表与 Basilisk TCW 求解器 |
| [film_rupture/](film_rupture/README.md) | 液膜破碎无排斥版v3，ls等效0.2：4种厚度×5种子及1条旧版对照，共21条轨迹；参数、源码、批处理与 hmin–t / V 监测 |
| [AMS/](AMS/README.md) | 2026-10-19重建的固定时间窗AMS；严格首次越级、完整状态分支、预算上限 |
| [DMC/](DMC/README.md) | 2026-10-19重建的DMC/QDMC；同物理核和评分，谱系及概率补偿 |
| [rare_events/](rare_events/README.md) | 两种采样框架共用的rupture参考物理核、构建、验收、独立直接模拟和有限试点 |
| [droplet_contact_angle/](droplet_contact_angle/README.md) | 前驱膜0.01的无涨落液滴接触角；保守AMR、方向索引修复、[控制与测量公式](droplet_contact_angle/METHODS.md)、[2026-10-21曲线/截面/数值结果](droplet_contact_angle/results/2026-10-21/README.md) |

TCW 原有内容整体移入 `tcw/`；在该目录内执行原 README 的命令，已有相对路径保持成立。绘图脚本中的外部原始数据目录仍需按本机数据位置设置。液膜破碎无排斥版v3已通过21组合1024步短测；长期、dt/网格及物理收敛仍未验收；两目录的滑移参数约定见各自源码与说明。
