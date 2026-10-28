# 周期液膜破碎：先认时间格式与后端

**版本入口：[五条主线与命名规则](../CODE_VERSIONS.md)。** 2种时间格式，5条实现路线；1/4/16核MPI是同一C99中点代码的资源配置。所有本次中点计算均直接MC。

| 想找什么 | 源码入口 |
|---|---|
| 同16条：CPU1/4/16核与A100耗时 | [192.36 / 47.99 / 12.23 / 34.00秒](midpoint/BENCHMARK.md) |
| 当前中点C / Basilisk C / CUDA | [midpoint/](midpoint/README.md)：三个源码优先展示 |
| 最早h1.2、cap=.1那组 | [SLE-SPLIT-FV-BASILISK](v3_tiered/source_v3_h1.2.c) |
| L6噪声双线性插值到L7/L8 | [噪声协议变体](l6_noise_bilinear/README.md) |
| 历史y索引实验 | [source_yold_h1.05.c](source_yold_h1.05.c) / [source_yfix_h1.05.c](source_yfix_h1.05.c) |

旧分步格式带.005截断；其守恒和索引检查不消除TCW/首达误差。中点内核没有夹值/体积修补/重抽；同路径拒绝和事件细分由runtime负责。短窗TCW、数值残差和MC首次破碎统计分别验收。

本目录根部`source_yfix_h*.c`、runs.csv、tasks.tsv及validation为2026-10-19原21身份历史证据，不代表当前生产状态。旧40分档和90噪声插值扫描各自保留说明。原部署、生成、监测和看板脚本集中留在[整理前提交](https://github.com/alpha1958627/SLE_TCW/tree/7e1d0959a22cc1a69ed65da2504914275ace0d2b/film_rupture)，当前展示不再让它们占据代码入口；本地完整档案与服务器文件未动。

共同外部高度记h0，计算初值H=1，`A=.0013159/h0², L=100/h0, ell=.2/h0`；固定phi版本用`.01/h0²`，显式phi接口另列。采用迁移率`H³+3bH²`约定时ell=3b。时间`tau=h0*t_code/1e8`。具体命令和运行依赖见[中点目录](midpoint/README.md)。
