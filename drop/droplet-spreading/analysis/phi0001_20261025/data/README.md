# 本地输入数据（尚未授权公开）
+
+公开仓库仅包含notebook、图像和说明，不包含模拟原场或完整半径时序。GitHub可直接查看已执行输出；完整复用须将本地数据包的以下文件放到此目录，或修改notebook的DATA_DIR：
+
+- native_fields_t1000.npz
+- radii_noise.csv
+- radii_deterministic.csv
+- metadata.json
+- manifest.sha256.json
+
+notebook会校验输入SHA、原场SHA及半径数值。数据缺失时不会生成替代或虚构数据。
+