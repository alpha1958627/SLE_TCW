# 历史hmin–t看板组件

这里的React/CSS对应2026-10-19旧21身份看板，不是当前中点批次状态页。当前代码入口见[版本总表](../../CODE_VERSIONS.md)和[中点源码](../midpoint/README.md)。

组件依赖原Data app的公共API；不是独立Vite项目。一次性JSON生成、SSH采集和刷新脚本已从当前源码展示移出，可在[整理前提交的film_rupture目录](https://github.com/alpha1958627/SLE_TCW/tree/7e1d0959a22cc1a69ed65da2504914275ace0d2b/film_rupture)查看。原本地证据和服务器监测未修改；不应把旧看板部署说明当成今天仍在线的运行状态。

历史横轴为各工况原始代码时间，未做跨高度物理时间换算。当前时间比较使用tau=h0*t_code/1e8，并保留删失、停止与失败。
