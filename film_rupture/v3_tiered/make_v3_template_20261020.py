# -*- coding: utf-8 -*-
"""2026-10-20 v3 模板: slip_prod_L6_dt3e-3_yfix2.c + 分档 dt cap
档位依据 2026-10-19 dt 上限测试 (viz/dt_limit_20261019.png):
  hmin >= 0.3 -> 1e-1 (悬挂段, 守恒机器精度, 提速 33x)
  hmin <  0.3 -> 1e-2 (俯冲段保护, 测试中 3e-2 也 100% 干净, 取 1e-2 更稳)
  hmin <  0.1 -> 3e-3 (近破碎精度段, 步数很少不耗时)
CONF 标记 yfix2=1 -> v3=1。
"""
from pathlib import Path
import re

SRC = Path('tools/slip_linear/slip_prod_L6_dt3e-3_yfix2.c')
OUT = Path('tools/slip_linear/slip_prod_L6_v3.c')
t = SRC.read_text(encoding='utf-8')

# --- 无排斥物理 (对齐 SLE_TCW 77ce5cb: hstar=0, 保留 h^-3 吸引) ---
oldhs = 'double hstar = pow(0.1, p.dd - 3.);'
assert t.count(oldhs) == 2, t.count(oldhs)
t = t.replace(oldhs, 'double hstar = 0.0;   /* v3: h^-5 repulsion disabled; h^-3 attraction retained (SLE_TCW 77ce5cb) */')
old = 'double cap = 3e-3;   /* dt 封顶 3e-3: 抑制长亚稳期质量漂移 */'
assert t.count(old) == 1, t.count(old)
new = ('double cap = hmin < 0.1 ? 3e-3 : (hmin < 0.3 ? 1e-2 : 1e-1);'
       '   /* v3: 分档 dt (dt-limit test 2026-10-19) */')
t = t.replace(old, new)
t = t.replace('fprintf(stdout, "CONF yfix2=1 slip=%g', 'fprintf(stdout, "CONF v3=1 slip=%g')

OUT.write_text(t, encoding='utf-8', newline='\n')
assert 'CONF v3=1' in t and 'hmin < 0.3 ? 1e-2 : 1e-1' in t and t.count('hstar = 0.0;') == 2
assert '(ix%nidx+nidx)%nidx' in t and 'alpha.y[0,1]' in t  # 守恒+y索引修复保留
for m in re.finditer(r'\.y\[\s*-?\d\s*\]', t):
    raise AssertionError(m.group(0))
print('v3 template OK:', OUT, len(t), 'chars')
