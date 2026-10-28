# 运行必需的Python，集中放在这里

主目录优先展示三个C/CUDA内核。本目录不含绘图、服务器账号、SSH、调度提交、竞跑或生产原场。

| 文件 | 必须保留的职责 |
|---|---|
| `midpoint_candidate.py` | NumPy参考、Newton/GMRES/FFT；C99类继承它 |
| `native_midpoint_candidate.py` | C99 ctypes绑定、构建、库SHA |
| `adaptive_runner.py` / `accelerated_runner.py` | 同路径BrownianTape、拒绝/事件细分、检查点/指纹 |
| `run_cpu_mpi.py` | 一个rank处理不同轨迹；完整CLI/输出/显式恢复 |
| `basilisk_fvsp_candidate.py` | FVSP谱曲率实验绑定和schema4，使用独立C库；不替换原FV |
| `basilisk_midpoint_candidate.py` | qcc库绑定及MPI路径；不是独立时间格式 |
| `run_single.py` | 最小reference/native/Basilisk有限窗入口 |
| `batched_film.py` | CuPy设备数组、FFT及CUDA Newton/GMRES控制 |
| `batch_runner.py` | GPU各lane独立时间/拒绝/事件与检查点 |
| `run_gpu_batch.py` | GPU批量CLI、输出和显式恢复 |
| `run_cpu_mpi.py --phi-ref` | CPU噪声强度参数，默认.01；替代原phi_parameterized四份重复代码 |

CPU扩展示例（从中点主目录执行）：

```sh
mpiexec -n 5 python runtime/run_cpu_mpi.py --heights 1.2 \
  --seeds 2026102701 2026102702 2026102703 2026102704 2026102705 \
  --phi-ref .013345 --dt-max .03 --until 300 --out runs/phi013345
```

CPU仅用NumPy/SciPy；CPU MPI另需mpi4py/MPI，GPU另需CuPy/CUDA。GPU没有调用CPU C99库做场求解；共享CPU模块为布朗协议和参考控制器提供定义。所有输出目录必须独立；不要用此包装目录恢复旧部署检查点。
