# 显式phi参数接口：仍是SLE-SMID-FV-C99

这里不是第六种数值格式。C内核与主目录`sle_smid_fv_c99.c`逐字节相同；四个必要Python文件给Film/布朗检查点/MPI入口增加显式`phi_ref`并核恢复指纹。共用上一层的C绑定，不复制C内核。

从中点主目录执行Linux命令：

```sh
export PYTHONPATH="$PWD/runtime${PYTHONPATH:+:$PYTHONPATH}"
mpiexec -n 5 python runtime/phi_parameterized/run_cpu_mpi.py --heights 1.2 \
 --seeds 2026102701 2026102702 2026102703 2026102704 2026102705 \
 --phi-ref .013345 --dt-max .03 --until 300 --out runs/phi013345
```

不要在运行中的轨迹或检查点中修改phi；.013345与.01是不同物理强度。四项速度表测的是固定phi_ref=.01接口，不宣称本扩展重新完成相同测速。
