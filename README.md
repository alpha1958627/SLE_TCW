# SLE_TCW

随机润滑方程（stochastic lubrication equation, SLE）在热毛细波（thermal capillary
wave, TCW）标度下的**集合谱作图脚本**。脚本读取求解器输出的高度场 `*.hfield`，做二维
FFT、按整数壳层做旋转平均，然后与线性解析解对照。

## 物理与公式

线性 TCW 方程每个波数模满足 `dH/dt = -k^4 H + 噪声`，本文这套求解器的噪声是按**面通量**
注入的（`f += dt * (flux.x[1] - flux.x[]) / Delta`，即 `d_t h = div(flux)`），傅里叶空间
多一个 `ik`，所以驱动功率 `∝ k^2`。以实测初谱 `H0`（第一帧 `t = 0.001`）为初值的解析演化是

```
|H(k,t)|^2 = |H0(k)|^2 exp(-2 k^4 tau) + varphi L^2 k^2 / k^4 * (1 - exp(-2 k^4 tau)),   tau = t - 0.001
```

- 实线：上式（`varphi = 0.001`）——`t = 0.001` 那条按构造就等于 `H0` 本身；
- 虚线：同样的 `H0` 衰减但 `varphi = 0`；
- 圆圈 + 点划线：SLE 集合数值谱。

`plot_tcw_H0_reference_pure_driven.py` 画的是另一个量——`H0 = 0` 的纯驱动解析解

```
|H(k,t)| = sqrt( varphi L^2 k^2 / k^4 * (1 - exp(-2 k^4 t)) )
```

用**绝对物理时刻** `t` 求值；因为三个时刻在 `k >~ 6` 都已饱和，所以它三条实线在高 k
**重合**（教科书纯理论图的样子）。代价是低 k / 早期它低于数据（那些模还保留着初值）。

## 算例参数

| 量 | continuous | deterministic / initial-noise |
|---|---|---|
| 网格 | L10, 1024 x 1024 | 同 |
| 域长 | `L0 = L = 100` | 同 |
| 噪声强度 | `varphi = 0.001` | `0`（初值 `h = 1 + 0.01 * noise()`） |
| 步长 | `dt = 2e-5`，5000 步 | 同 |
| 输出帧 | `t = 0.001, 0.01, 0.1` | 同 |
| realization | 10 | 30 |

数据路径在脚本顶部的 `RESULTS` / `base.L10` 常量里，指向本机结果目录；换机器时改这两处即可。

## 文件

| 文件 | 说明 |
|---|---|
| `tcw_common.py` | 公共模块：读 `hfield`、FFT 壳层平均、解析解函数 |
| `plot_tcw_H0_reference.py` | 主图：`H0` 锚定的解析演化（16:9 大图 + 5:3 小图） |
| `plot_tcw_H0_reference_pure_driven.py` | 纯驱动版（`H0 = 0`，绝对时刻），高 k 三条实线重合 |
| `diagnose_highk_dt_residual.py` | 诊断脚本：对账 `k > 8` 残差的两个候选机制 |
| `figures/` | 输出：无噪声版大图 / 小图（PNG + PDF）与 CSV 数据表 |

大图 16:9（3520 x 1980）、小图 5:3（1650 x 990），保存时不裁边，比例精确，便于拼版。

## 运行

```
pip install numpy scipy matplotlib
python plot_tcw_H0_reference.py
python plot_tcw_H0_reference_pure_driven.py
```

文字用 LaTeX 渲染（`text.usetex = True`），需要本机有 LaTeX。

## 已知残差：`k > 8` 数值高于理论 6~15%

三个时刻一致、静态、随 k 增大（k=5 约 +1%，k=8 约 +8%，k=10 约 +15%），不是公式错误，
而是求解器**有限时间步**造成的：模的弛豫率 `omega = k^4`，`dt = 2e-5` 时

```
k = 5   ->  omega*dt = 0.0125
k = 8   ->  omega*dt = 0.082
k = 10  ->  omega*dt = 0.20
```

离散格式的稳态方差是 `q / (1 - g^2)` 而不是连续的 `q / (2 omega dt)`。本求解器的确定性
部分是隐式的（`g ≈ 1/(1 + omega*dt)`），于是振幅偏高一个因子

```
(1 + omega*dt) / sqrt(1 + omega*dt/2)      ->  k=8 约 +6%，k=10 约 +14%
```

把观测比值除以这个因子后，`k = 1..10` 全部回到 `1.00 ± 0.03`（三个时刻一致）。
`diagnose_highk_dt_residual.py` 复现这套对账；同时排除了乘性噪声 `pow(ff, 1.5)` 的
候选解释（实测只贡献 +0.15%，因为 `var(h) = 0.0022`）。

推论：`k > 8` 那一小段的偏差会随 `dt -> 0` 消失；若要理论线与数据在高 k 也严格重合，
给驱动项乘上上面那个离散因子即可。
