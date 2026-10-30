# 双液滴聚合：无涨落与 OU–SE

固定预加密的二维 Basilisk 双液滴试验；`coalescence.ipynb` 读取本地计算输出，绘制液桥高度、中心连线截面和俯视图。一个单元的膜厚首次严格超过 `0.03` 后永久允许随机源；为保持离散守恒，随机面通量只有在其全部非零梯度模板单元均已激活时才开启。OU 潜在场连续演化，单元激活历史不因高度降低而消失。

当前为实验性串行、固定网格求解器：动态 AMR 历史转移和包含 OU 隐状态的断点续算尚不支持。噪声门控没有同步修改确定性耗散，不能自动视为满足热平衡涨落–耗散定理。前驱膜始终存在，阈值连通是操作定义而非分子首次接触；一维中心线最小值是液桥高度代理，不等同于二维鞍点。

## 编译及运行

设置 `BASILISK` 指向已安装的 Basilisk 源目录，提供 `gcc` 和 `qcc`：

```sh
export BASILISK=/path/to/basilisk/src
export CC99=gcc
sh build.sh
mkdir -p data/det data/ou
(cd data/det && ../../droplet -1 103203 0 0.003125 32000 8 0.002 0.001 0.03)
(cd data/ou  && ../../droplet  2 103203 0 0.003125 32000 8 0.002 0.001 0.03)
```

最后三个参数分别为 `phi`、保守膜厚下界和永久激活阈值；无涨落模式会把 phi 设为零。示例是 `t=100`，每个完整轨迹可能耗费数小时；先做短测再运行。每个输出目录必须独立，保留 `runtime_config.json`、`activation_config.json`、`metrics.csv`、`activation_events.csv`、`activation_metrics.csv`、`activation_final.bin` 和 `*.native`。每 0.1 保存一个原生高度场。运行期间不可把不完整的场当成结果。

参数：方域 `L=64`，两滴中心 `(27.25,32)` 和 `(36.75,32)`，半径均为 5，初始超额高度为 `cos^4`、峰高 1，前驱膜 `hp=0.01`；最细 `dx=0.25`，初始静态细网格保护半径 14，过渡半径 18。OU `sigma=0.05`、`tau=0.01`；其相关长度小于流体网格，需另做分辨率研究。二元格式为 20 字节头 `<iQd` 加各叶单元 `[x,y,Delta,h,p,level]` 的 6 个 float64。请确认输入实际文件哈希、模拟身份和可比物理时刻。

## 画图

安装 NumPy、Matplotlib 和 Jupyter，再从此目录启动 `coalescence.ipynb`。默认读取 `data/det` 和 `data/ou`；也可通过环境变量 `COALESCENCE_DATA_DIR` 指向外部包含 `det/` 和 `ou/` 的目录。notebook 比较实际共同时间步，不外推或填补缺失数据。需要 LaTeX、dvipng、dvips 生成论文风格 PDF 和 600 dpi PNG。

输出数据和计算日志默认不进入 Git；数据由使用者在自己的工作环境中保留。不要将同 seed 的回放当作独立样本。
