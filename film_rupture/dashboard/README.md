# hmin–t 看板源码

`src/content/dashboard/` 是已部署看板的自定义 React 组件和 CSS。组件依赖 Codex Data Analytics 1.0.11 的 Data app 公共 API `src/data-app-public.jsx`，不能当作独立 Vite 项目启动。此处不复制该插件的受保护运行时或打包依赖。

在已安装该插件的环境生成规范 Data app，将这两个文件放进相同相对路径，设置 `DASHBOARD_APP` 为该 app 的绝对目录。`build_dashboard_data.py` 读取 `deployment.json` 与 `audits/<site>_final.json`，写出 `dashboard_snapshot.json`；若 app 的 `src/data.json` 已存在，也原位更新它。三个查询是 `film_history`、`run_status`、`field_volume`。

复制 `deployment.example.json` 为 `deployment.json`，填实际 SSH 别名、远端根目录、Python 路径、数组作业号和每条轨迹的保留位置。远端必须已有本项目 `monitor_volume.py`、`runs.csv` 和算例输出。然后从 `film_rupture/` 运行：

```sh
python refresh_dashboard.py
```

它只读远端已有文件与 Slurm 状态，不提交、取消、重启或上传文件。全部站点采集成功才替换看板数据；采集错误保存到 `audits/refresh_errors.json`。日志较长时保留每桶首尾及 hmin 极值，原始下载日志保留在本地 `audits/raw_logs/`。

设置 `DATA_APP_SCRIPT` 为插件 `scripts/data-app.mjs` 的路径，可在刷新后构建；`NODE` 可指定 Node 可执行文件。没有插件时仍可生成并检查 JSON 快照。每 12 小时调用一次刷新由外部调度器负责，仓库源码本身不会创建定时任务；当前已部署看板仍使用原本的 Codex 12h 自动化。

横轴是各工况的原始数值时间，没有跨厚度物理时间换算。紫色虚线是旧版 h1.05/run1 对照；全场体积图使用原始场直接积分。下载、构建与展示成功均不代表数值收敛。
