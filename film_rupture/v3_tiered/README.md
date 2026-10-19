# v3 tiered-dt 扫描（ls=0.2 × r=0.5–1.2 × 5 seed）

2026-10-20 由 lead 会话生成并部署（SIPF 作业 29333）。基于本仓库"无排斥版 v3"物理
（hstar=0，SLE_TCW 77ce5cb）+ 全部守恒修复（ee3646d）+ **分档 dt cap** + 剖面/液面输出。

## 本目录新增内容（相对 film_rupture 当前版）

1. **分档 dt cap**（`slip_prod_L6_v3.c`，CONF 标记 `v3=1`）：

   ```c
   double cap = hmin < 0.1 ? 3e-3 : (hmin < 0.3 ? 1e-2 : 1e-1);
   ```

   依据 2026-10-19 dt 上限测试（XH，MPI8，本代码族）：
   - 悬挂段（hmin≥0.3）dt=1e-1 完全安全（守恒 ~1e-12，成核统计无偏，t=1e5 仅 ~1–2 小时）；
   - 俯冲段 dt=1e-1 有 1/3 SIGFPE（hmin≈0.04–0.22 处浮点异常），dt≤3e-2 全生命周期 100% 干净；
   - 破碎时间在 3e-3→1e-1 无系统偏移（GM 576/651/673/606，3-seed 散布内）。
   分档恰好覆盖观测到的崩溃区（<0.3 与 <0.1 两档）。

2. **剖面与液面输出**（供复现）：
   - `sect_hmin_%06.2f_r%d.dat`：hmin 每跌破一档 0.01（0.99→0.01，~99 帧，while 防快潜跳档）；
   - `sect_t_%09.0f_r%d.dat`：定时（每 1000 t，~100 帧）；
   - `hfield_t%09.1f.bin`：全场每 100 t；
   - 破碎瞬间 per-rank `rupture.hfield/rupture_field_r*/rupture_section_r*`（22712 修复版）。
   均不消耗随机数、不改变轨迹；MPI 按进程分文件，本地合并即可。

3. **扫描范围**：r=0.5,0.6,…,1.2（8 组）× 5 seed = 40 条轨迹，换算口径
   （A=0.0013159/r²，varphi=0.01/r²，L0=100/r，ell_code=0.2/r，代码初值 h=1）。
   seed=2026112100+hidx×10+run。level 6、MPI8、步限 2e8、hmin≤0.01 停车。

## 文件

| 文件 | 说明 |
|---|---|
| `slip_prod_L6_v3.c` | 模板（分档 dt + 全部修复；物理参数为 r=1.0 默认） |
| `source_v3_h{r}.c` | 8 份预补丁源（常数按 r 烘入，sha256 见 `sources.sha256`） |
| `make_v3_template_20261020.py` / `gen_v3_sources_20261020.py` | 生成链（Decimal 40 位精度，与 build_plan.py 同口径） |
| `run_v3_sipf.slurm` | SIPF C64 提交脚本（array 0–39%40） |

## 口径备注

- 滑移迁移率沿用 `m(h)=h^3(1+ell/h)`（= h^3+ell_code·h^2，与本目录 ell_code 约定一致）；
  噪声通量乘 sqrt(1+ell/h)（FDT 一致）。
- 短测（dt 分档切换、输出文件齐全）在提交前经 CONF/日志抽查；**长期与 dt/网格收敛未验收**，
  与主批同等的证据边界适用。
