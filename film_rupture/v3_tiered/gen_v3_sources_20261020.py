# -*- coding: utf-8 -*-
"""2026-10-20 v3 批源生成: ℓ=0.2 × h0=0.5-1.2 (步长0.1) × 5 run, 换算口径, 分档 dt
输出 stage_v3_20261020/source_v3_h{0.5,...,1.2}.c (含 hfield 每100t + 截面 hmin每0.01/每1000t)
"""
import re
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 40
src = Path('tools/slip_linear/slip_prod_L6_v3.c').read_text(encoding='utf-8')

periodic = (
    'event hfield_periodic (t = 0.; t += 100.) {\n'
    '    char name[128];\n'
    '    snprintf(name, sizeof(name), "hfield_t%09.1f.bin", t);\n'
    '    write_hfield_binary(name);\n'
    '    if (pid() == 0) fprintf(stderr, "HFIELD t=%.12g file=%s\\n", t, name);\n'
    '}\n\n'
    'double sect_hmin_next = 0.99;\n\n'
    'event section_hmin(i++) {\n'
    '    while (hmin < sect_hmin_next && sect_hmin_next > 0.004) {\n'
    '        char nms[64];\n'
    '        snprintf(nms, sizeof nms, "sect_hmin_%06.2f_r%d.dat", sect_hmin_next, pid());\n'
    '        FILE *fs = fopen(nms, "w");\n'
    '        double ymid = L0 / 2.;\n'
    '        foreach()\n'
    '            if (fabs(y - ymid) < Delta)\n'
    '                fprintf(fs, "%.6g %.6g\\n", x, h[]);\n'
    '        fclose(fs);\n'
    '        sect_hmin_next -= 0.01;\n'
    '    }\n'
    '    return 0;\n'
    '}\n\n'
    'event section_time (t = 1000.; t += 1000.) {\n'
    '    char nms[64];\n'
    '    snprintf(nms, sizeof nms, "sect_t_%09.0f_r%d.dat", t, pid());\n'
    '    FILE *fs = fopen(nms, "w");\n'
    '    double ymid = L0 / 2.;\n'
    '    foreach()\n'
    '        if (fabs(y - ymid) < Delta)\n'
    '            fprintf(fs, "%.6g %.6g\\n", x, h[]);\n'
    '    fclose(fs);\n'
    '    if (pid() == 0) fprintf(stderr, "SECT_T t=%.12g\\n", t);\n'
    '}\n\n'
)
assert 'event hfield_periodic' not in src
assert src.count('event rupture_stop') == 1
OUT = Path('stage_v3_20261020')
OUT.mkdir(exist_ok=True)

H0S = ['0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1.1', '1.2']
slips = {}
for h0s in H0S:
    h0 = Decimal(h0s)
    text = src

    def replace_once(pattern, replacement, value, label):
        value, count = re.subn(pattern, replacement, value, count=1)
        assert count == 1, (label, count)
        return value

    a = format(Decimal('0.0013159') / (h0 * h0), '.17g')
    phi = format(Decimal('0.01') / (h0 * h0), '.17g')
    l0 = format(Decimal(100) / h0, '.17g')
    text = replace_once(r'double tcw_A = [-+0-9.eE]+;', f'double tcw_A = {a};', text, 'tcw_A')
    text = replace_once(r'double tcw_varphi = [-+0-9.eE]+;', f'double tcw_varphi = {phi};', text, 'tcw_varphi')
    text = replace_once(r'L0\s*=\s*100\s*;', f'L0 = {l0};', text, 'L0')
    idx = text.find('event rupture_stop')
    text = text[:idx] + periodic + text[idx:]
    for marker in ('CONF v3=1', 'sqrt(1.0 + tcw_slip/ff)', 'rupture_field_r%d.dat',
                   'hmin < 0.3 ? 1e-2 : 1e-1', f'double tcw_A = {a};',
                   f'double tcw_varphi = {phi};', f'L0 = {l0};',
                   'event hfield_periodic', 'event section_hmin(i++)',
                   'sect_hmin_%06.2f_r%d.dat', '(ix%nidx+nidx)%nidx',
                   'alpha.y[0,1]', 'f_old[0,-1]'):
        assert marker in text, marker
    p = OUT / f'source_v3_h{h0s}.c'
    p.write_text(text, encoding='utf-8', newline='\n')
    slips[h0s] = format(Decimal('0.2') / h0, '.17g')
    print(h0s, 'A=' + a, 'phi=' + phi, 'L0=' + l0, 'slip=' + slips[h0s], len(text))
print('staged in', OUT)
