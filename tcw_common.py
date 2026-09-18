"""Plot TCW spectra using each numerical first saved frame as H0."""
from pathlib import Path
import struct

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.ndimage import uniform_filter1d

L10 = Path(r"F:\alpha_cx\罕见事件理论\results\tcw_L10_sqrt2_phi0p001_10runs_20260728\raw")
ONE_KICK = Path(r"F:\alpha_cx\罕见事件理论\results\firstkick_30run_20260804")
OUT = Path(r"F:\MD\jingbang\tcw_L10_H0_reference_20260804")
TIMES = (0.001, 0.01, 0.1)
T_REF = TIMES[0]
L10_FILES = {0.001: "tcw_t0000.hfield", 0.01: "tcw_t0001.hfield", 0.1: "tcw_t0002.hfield"}
KICK_FILES = {0.001: "tcw_t0000.hfield", 0.01: "tcw_t0003.hfield", 0.1: "tcw_t0012.hfield"}
COLORS = {0.001: "#000000", 0.01: "#d62828", 0.1: "#2455d6"}
LENGTH = 100.0
PHI = 1.0e-3

mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 20,
    "axes.labelsize": 29,
    "legend.fontsize": 19
})


def read_hfield(path: Path) -> np.ndarray:
    raw = path.read_bytes()
    n = struct.unpack_from("i", raw, 0)[0]
    return np.frombuffer(raw, dtype="<f8", count=n*n, offset=4).reshape(n, n)


def theory_from_h0(h0: np.ndarray, k: np.ndarray, tau: float, driven: bool) -> np.ndarray:
    omega = k**4
    power = h0**2 * np.exp(-2.0 * omega * tau)
    if driven:
        eps = 1e-12
        power += PHI * LENGTH**2 * k**2 / (omega + eps) * (1.0 - np.exp(-2.0 * omega * tau))
    return np.sqrt(power)

def pure_driven_spectrum(k: np.ndarray, t_phys: float):
    """纯驱动解析解，h0=0，t=0开启噪声，t_phys为真实物理时刻，不是相对tau"""
    omega = k**4
    eps = 1e-12
    power = PHI * LENGTH**2 * k**2 / (omega + eps) * (1.0 - np.exp(-2.0 * omega * t_phys))
    return np.sqrt(power)


def ring_spectra(root: Path, pattern: str, files: dict[float, str]) -> dict[float, tuple[np.ndarray, np.ndarray]]:
    runs = sorted(root.glob(pattern))
    if not runs:
        raise RuntimeError(f"no runs found in {root}")
    output = {}
    for time in TIMES:
        fields = [read_hfield(run / files[time]) for run in runs]
        n = fields[0].shape[0]
        if any(field.shape != (n, n) for field in fields):
            raise ValueError(f"inconsistent grid at t={time:g}")
        dx = LENGTH / n
        power = np.mean([
            np.abs(np.fft.fftshift(np.fft.fft2(field - field.mean())))**2
            for field in fields
        ], axis=0)
        rms = np.sqrt(power) * dx**2
        k2d = np.fft.fftshift(2.0 * np.pi * np.fft.fftfreq(n, d=dx))
        kx, ky = np.meshgrid(k2d, k2d)
        shell = np.rint(np.hypot(kx, ky) / (2.0 * np.pi / LENGTH)).astype(np.int64)
        count = np.bincount(shell.ravel())
        shell_power = np.bincount(shell.ravel(), weights=(rms * rms).ravel())
        valid = np.arange(1, len(count))
        output[time] = valid * 2.0 * np.pi / LENGTH, np.sqrt(shell_power[valid] / count[valid])
    return output


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    continuous = ring_spectra(L10, "repeat_run_*", L10_FILES)
    one_kick = ring_spectra(ONE_KICK, "run_*", KICK_FILES)

    k_cont0, h0_cont = continuous[T_REF]
    k_kick0, h0_kick = one_kick[T_REF]

    # 平滑用于理论计算的h0，消除初场统计噪声
    h0_cont = uniform_filter1d(h0_cont, size=4, mode="nearest")
    h0_kick = uniform_filter1d(h0_kick, size=4, mode="nearest")

    fig, ax = plt.subplots(figsize=(11.0, 6.0), constrained_layout=True)

    rows = []

    for time in TIMES:
        tau = time - T_REF
        k_cont, y_cont = continuous[time]
        k_kick, y_kick = one_kick[time]

        cont_theory = theory_from_h0(h0_cont, k_cont0, tau, driven=True)
        kick_theory = theory_from_h0(h0_kick, k_kick0, tau, driven=False)
        pure_theory = pure_driven_spectrum(k_cont0, time)

        color = COLORS[time]

        # 数值圆圈，降采样
        ax.plot(
            k_cont[::3], y_cont[::3],
            linestyle="none",
            marker="o",
            markersize=9,
            markerfacecolor="none",
            markeredgecolor=color,
            markeredgewidth=1.3,
            label=rf"numerical continuous (L10 $\sqrt{{2}}$), $t={time:g}$"
        )

        # 数值菱形，降采样
        ax.plot(
            k_kick[::3], y_kick[::3],
            linestyle="none",
            marker="D",
            markersize=8,
            markerfacecolor="none",
            markeredgecolor=color,
            markeredgewidth=1.3,
            label=rf"numerical one kick (L9 Gaussian), $t={time:g}$"
        )

        # 带模拟H0的理论实线（主曲线）
        ax.plot(
            k_cont0, cont_theory,
            color=color,
            linestyle="-",
            linewidth=2.2,
            label=rf"theory continuous from $H_0$, $t={time:g}$"
        )

        # 【新增】纯驱动解析解 h0=0，浅灰细实线做对比
        ax.plot(
            k_cont0, pure_theory,
            color="#aaaaaa",
            linestyle="-",
            linewidth=1.4,
            alpha=0.75
        )

        ax.plot(
            k_kick0, kick_theory,
            color=color,
            linestyle="--",
            linewidth=2.2,
            label=rf"theory decay from $H_0$, $t={time:g}$"
        )

        rows.extend((time, tau, k, y, "L10_continuous_numeric") for k, y in zip(k_cont, y_cont))
        rows.extend((time, tau, k, y, "L9_one_kick_numeric") for k, y in zip(k_kick, y_kick))
        rows.extend((time, tau, k, y, "continuous_theory_from_H0") for k, y in zip(k_cont0, cont_theory))
        rows.extend((time, tau, k, y, "one_kick_theory_from_H0") for k, y in zip(k_kick0, kick_theory))
        rows.extend((time, time, k, y, "pure_driven_h0_eq_0") for k, y in zip(k_cont0, pure_theory))

    ax.set(xlim=(0.0, 8.1), xlabel=r"$|k|$", ylabel=r"$|H|_{\rm rms}$")
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=22, length=7)

    legend_handles = [
        Line2D([], [], color="black", linewidth=2.4, linestyle="-", label=r"theory, $\varphi=0.001$"),
        Line2D([], [], color="black", linewidth=2.4, linestyle="--", label=r"theory, $\varphi=0$"),
        Line2D([], [], color="black", marker="o", markersize=9, markerfacecolor="none",
               markeredgewidth=1.3, linestyle="none", label=r"SLE, $\varphi=0.001$"),
        Line2D([], [], color="black", marker="D", markersize=8, markerfacecolor="none",
               markeredgewidth=1.3, linestyle="none", label=r"SLE, $\varphi=0$"),
        Line2D([], [], color="#aaaaaa", linewidth=1.4, linestyle="-", label=r"pure driven $H_0=0$")
    ]

    ax.legend(
        handles=legend_handles,
        ncol=1,
        loc="upper right",
        frameon=False,
        columnspacing=1.8,
        handlelength=2.7,
        fontsize=22
    )

    fig.savefig(OUT / "tcw_H0_reference_continuous_and_onekick.png", dpi=220)
    fig.savefig(OUT / "tcw_H0_reference_continuous_and_onekick.pdf", bbox_inches="tight")

    np.savetxt(
        OUT / "tcw_H0_reference_continuous_and_onekick.csv",
        rows,
        delimiter=",",
        fmt="%s",
        header="physical_time,tau_from_H0,k,rms,series",
        comments=""
    )


if __name__ == "__main__":
    main()