# 控制方程、理论接触角与测量公式

以下均为程序使用的无量纲量；本组结果取 $h_*=0.01$、$\gamma=1$、噪声幅度为零。代码中的 $k$ 是分离压的强度系数。它控制短程排斥与长程吸引的共同幅度，不是接触角本身。

## 1. 确定性薄膜方程

二维周期底面上的高度 $h(x,y,t)$ 满足

```math
\partial_t h=\nabla\cdot\left[h^3\nabla p\right],
\qquad p=-\gamma\nabla^2h-\Pi(h).
```

本源的毛细项使用小斜率曲率 $\nabla^2h$；没有使用完整曲率
$\nabla\cdot[\nabla h/\sqrt{1+|\nabla h|^2}]$。
无量纲迁移率是 $h^3$，时间尺度已吸收黏性前因子。

分离压为

```math
\Pi(h)=k\left[\left(\frac{h_*}{h}\right)^n-
                   \left(\frac{h_*}{h}\right)^m\right],
\qquad n=5,\quad m=3,\quad n>m>1.
```

$h<h_*$ 时 $\Pi>0$，短程排斥占优；$h>h_*$ 时 $\Pi<0$，长程吸引占优；$\Pi(h_*)=0$。
前驱膜参数 $h_*$ 是势能最小点，不是每个时间步强制实际膜厚。
与源码系数的对应关系为

```math
A=k h_*^3,\qquad
\Pi(h)=\frac{A h_*^2}{h^5}-\frac{A}{h^3},\qquad
p=-\nabla^2 h+\frac{A}{h^3}-\frac{A h_*^2}{h^5}
\quad(\gamma=1).
```

其中 $A$ 对应 `COEFF_A`，是本代码的无量纲吸引系数。
可用 $\Phi'(h)=-\Pi(h)$ 定义界面势：

```math
\Phi(h)=k\left[\frac{h_*^n}{(n-1)h^{n-1}}
                    -\frac{h_*^m}{(m-1)h^{m-1}}\right].
```

## 2. 完整曲率标定与小斜率参考必须分开

定义正数

```math
C=\frac{n-m}{(n-1)(m-1)},\qquad
\int_{h_*}^{\infty}\Pi(h)\,dh=-Ckh_*.
```

在局部宏观界面极限，完整曲率的Young关系给出

```math
\gamma(1-\cos\theta_{\rm cal})=-\int_{h_*}^{\infty}\Pi(h)\,dh,
\qquad
\theta_{\rm cal}=\arccos\left(1-\frac{Ckh_*}{\gamma}\right).
```

本说明取 $\Pi>0$ 为短程排斥，因此上面积分为负。对 $n=5,m=3$，$C=1/4$，故

```math
\theta_{\rm cal}=\arccos\left(1-\frac{kh_*}{4\gamma}\right),
\qquad
k=\frac{4\gamma(1-\cos\theta_{\rm cal})}{h_*}.
```

这就是15°、25°、35°**标定标签**的来源。当前程序采用小斜率曲率，相应的宏观坡度与几何参考角为

```math
s_{\infty}^2=-\frac{2}{\gamma}\int_{h_*}^{\infty}\Pi(h)\,dh
            =\frac{2Ckh_*}{\gamma},
\qquad
\theta_{\rm model}=\arctan(s_{\infty})
                  =\arctan\sqrt{\frac{kh_*}{2\gamma}}.
```

表中的角度均从弧度乘 $180/\pi$ 转为度。`theta_lub_degrees` 是坡度 $s_\infty$ 的角度单位表达，不能代替 $\arctan(s_\infty)$。

| $k$ | $A$ | 完整曲率标定 | 小斜率几何参考 |
|---:|---:|---:|---:|
| 13.62966948 | 1.362966948e-5 | 15° | 14.63068090° |
| 37.47688519 | 3.747688519e-5 | 25° | 23.40678414° |
| 72.33918228 | 7.233918228e-5 | 35° | 31.02318910° |

这些是宏观界面极限参考，不能强迫有限液滴在尚有压力差、实际膜厚偏移或时间/网格误差时恰好等于该角度。

## 3. 初场与相同体积

底面周期边长 $L=32$，$r=\sqrt{x^2+y^2}$：

```math
h(r,0)=h_*+H_0\cos^4\left(\frac{\pi r}{2R_0}\right)
\quad(r<R_0),\qquad h(r,0)=h_*\quad(r\ge R_0),
\qquad H_0=1,\ R_0=5.
```

三 $k$、两网格共享此初场。修复版实际原生初始总体积约23.776936818258157；以每例初始化保存的原 $V_0$ 验收，不能在续算时重置。

## 4. 主体外推接触角

对原生网格采样的 $x^\pm/y^\pm$ 四条截面，取主体超额高度窗口

```math
0.4<\frac{h-h_*}{h_{\max}-h_*}<0.8,
\qquad h(r)\simeq a+b r^2,\quad b<0.
```

外推至 $h=h_*$，得到

```math
R_{\rm fit}=\sqrt{\frac{h_*-a}{b}},\qquad
s_{\rm fit}=2(-b)R_{\rm fit}=2\sqrt{(-b)(a-h_*)},\qquad
\theta_{\rm fit}=\arctan(s_{\rm fit}).
```

已发布散点是各组**最新完整帧的四方向均值**，与截面角度标签一致；误差棒是四方向min/max，不是置信区间或时间不确定性。另保存20–60%、60–90%窗口与圆弧拟合，便于审查测法偏差。

局部图以 $(r-R_{\rm fit},h-h_*)$ 等比例画出外推切线和角弧。$R_{\rm fit}$ 不等于监控中四方向首次 $h=0.011$ 交点半径，也不等于包含外环的全阈值面积半径。

圆弧对照拟合 $r^2+(h-c)^2=R_c^2$，角度为 $\arccos[(h_*-c)/R_c]$。本模型的主法仍是抛物线；应同时检查拟合残差与窗口依赖，不能挑选最接近标定标签的拟合法。

## 5. 体积与AMR验收

```math
V(t)=\int h\,dx\,dy\simeq\sum_{i\in\mathrm{native\ leaves}}h_i\Delta_i^2,
\quad V_{\rm signed}=V-h_*L^2,
\quad V_{+}=\int\max(h-h_*,0)\,dx\,dy.
```

```math
\epsilon_V^{\max}=\max_{\tau\le t}\frac{|V(\tau)-V_0|}{V_0}<10^{-8},
\qquad
\epsilon_{\rm AMR}=\frac{|V_{\rm after}-V_{\rm before}|}{V_{\rm before}}<10^{-10}.
```

包含初始化、每次AMR前后和所有续算历史；原生积分与保守导出独立补偿积分相对差要求小于 $10^{-12}$。不以高度裁剪或体积强制归一化掩盖漂移。

稳定验收还需三个连续50模型时间窗面积及四向半径范围/均值≤0.1%，末50完整场及双向截面变化≤峰超额高度0.1%，再检查末三帧四向角度和同时间L13/L14：半径差≤0.2%、角度差≤0.2°。压力、物理RHS和必要的时间步一致性单独检查。MG/Picard≤ $10^{-8}$、作业COMPLETED或达到 $t=400$ 均不能代替这些条件。

当前数据、图及逐组状态见 [2026-10-21结果](results/2026-10-21/README.md)。
