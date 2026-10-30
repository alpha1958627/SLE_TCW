# SLE_TCW

随机润滑方程：热毛细波、液膜破碎、液滴与罕见事件采样。

| 研究内容 | 入口 |
|---|---|
| 液膜破碎 | [film_rupture/](film_rupture/README.md) · [中点求解器](film_rupture/midpoint/README.md) |
| 液滴 | [drop/](drop/README.md) · [铺展 Notebook](drop/droplet-spreading/analysis/phi0001_20261025/drop.ipynb) · [液滴聚合](drop/droplet-coalescence/README.md) · [接触角](drop/pinned/README.md) |
| 热毛细波 | [tcw/](tcw/README.md) · [谱分析 Notebook](tcw/tcw.ipynb) |
| AMS / DMC / QDMC | [rare_events/](rare_events/README.md) |

平膜求解器的格式、后端与历史测速见 [CODE_VERSIONS.md](CODE_VERSIONS.md)。

日常使用 `main`。独立实验放在各研究目录的 `experiments/` 中：
[液滴毛细中点](drop/droplet-spreading/experiments/capillary_midpoint_20261027/README.md)、
[q17 评分](rare_events/experiments/width_q17_20261021/README.md)。实验代码合入目录不改变默认模型，也不表示物理验收。

2026-10-28 整理：合并开发分支，移除重复 Python 与一次性提交/监控脚本；
[文件去向与原分支快照](CODE_VERSIONS.md#仓库整理记录)。已有生产任务继续使用其原冻结源码和环境；整理后的源码指纹不接管旧检查点。
