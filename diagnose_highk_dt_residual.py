"""诊断 k>8 残差：乘性噪声温度重整化 vs 有限时间步稳态方差膨胀。"""
from pathlib import Path

import numpy as np

import tcw_common as base

DT = 2.0e-5
PHI, L = base.PHI, base.LENGTH
cont = base.ring_spectra(base.L10, "repeat_run_*", base.L10_FILES)
runs = sorted(base.L10.glob("repeat_run_*"))

print(f"dt = {DT:g}   phi = {PHI}   L = {L}")
print("\n--- (a) 乘性噪声因子 sqrt(<ff^3>)，ff = 面平均 h ---")
for t in base.TIMES:
    vals = []
    for run in runs:
        h = base.read_hfield(run / base.L10_FILES[t])
        ff = 0.5 * (h + np.roll(h, 1, axis=0))
        vals.append(np.mean(ff**3))
        var = float(np.var(h))
    factor = float(np.sqrt(np.mean(vals)))
    print(f"  t={t:<6g} mean(h)={np.mean(h):.6f}  var(h)={var:.6g}  sqrt(<ff^3>)={factor:.5f}")

print("\n--- (b) 有限时间步稳态方差膨胀 ---")
print("  后向 Euler: V/V_cont = 2*w*dt/(1-g^2), g=1/(1+w*dt)  =>  振幅因子 (1+w*dt)/sqrt(1+w*dt/2)")
print("  精确分裂  : g=exp(-w*dt)                            =>  振幅因子 sqrt(2*w*dt/(1-exp(-2*w*dt)))")
print(f"  {'k':>7} {'w*dt':>9} {'后向Euler':>10} {'分裂':>9}")

print("\n--- 观测比值 numeric/theory ---")
print(f"  {'k':>7} " + " ".join(f"{f't={t:g}':>9}" for t in base.TIMES))
ks = (1.0, 2.0, 3.0, 5.0, 6.0, 8.0, 10.0)
ratios = {}
for kk in ks:
    wdt = kk**4 * DT
    be = (1 + wdt) / np.sqrt(1 + wdt / 2)
    sp = np.sqrt(2 * wdt / (1 - np.exp(-2 * wdt)))
    print(f"  {kk:>7.1f} {wdt:>9.4f} {be:>10.4f} {sp:>9.4f}")
for t in base.TIMES:
    k, y = cont[t]
    theory = np.sqrt(PHI * L**2 * k**2 / k**4 * (1.0 - np.exp(-2.0 * k**4 * t)))
    row = []
    for kk in ks:
        i = int(np.argmin(np.abs(k - kk)))
        row.append(y[i] / theory[i])
    ratios[t] = row
    print(f"  {t:>7g} " + " ".join(f"{r:>9.3f}" for r in row))

print("\n--- 观测比值 / (a)因子 / (b)后向Euler  (应接近 1) ---")
for t in base.TIMES:
    h = base.read_hfield(runs[0] / base.L10_FILES[t])
    ff = 0.5 * (h + np.roll(h, 1, axis=0))
    fac = float(np.sqrt(np.mean(ff**3)))
    out = []
    for kk, r in zip(ks, ratios[t]):
        wdt = kk**4 * DT
        be = (1 + wdt) / np.sqrt(1 + wdt / 2)
        out.append(r / fac / be)
    print(f"  t={t:<6g} " + " ".join(f"{o:>7.3f}" for o in out))
