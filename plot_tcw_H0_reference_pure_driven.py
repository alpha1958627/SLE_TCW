"""Same layout as plot_tcw_H0_reference.py, but with a purely analytic H0 = 0 solution.

The solid curves are the textbook driven solution of the linear TCW equation with
zero initial spectrum and white-in-time conserved forcing switched on at t = 0:

    |H(k,t)| = sqrt( varphi L^2 k^2 / k^4 * (1 - exp(-2 k^4 t)) )

evaluated at the absolute physical time t.  Because every mode has saturated by
k >~ 6 for all three saved times, the three solid curves coincide there, which is
the behaviour of the standard pure-theory plot.

Dashed curves are the varphi = 0 decay of the measured initial spectrum H0 (this
panel keeps H0), and the open circles joined by dash-dot lines are the SLE
ensemble spectra of the 10-run L10 campaign.  Noise-free version: no jitter is
applied (JITTER_MAX = 0).

Trade-off: dropping H0 makes the solid curve the driven attractor, so at low k and
early time it falls below the data (those modes keep their initial condition);
the continuum attractor is also 6-15% below the data for k > 8 because of the
finite solver timestep dt = 2e-5.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import numpy as np

import tcw_common as base


RESULTS = Path(r"F:\alpha_cx\罕见事件理论\results")
DET = RESULTS / "tcw_L10_initial_noise0p01_det_30runs_20260804"
OUT = Path(__file__).resolve().parent / "figures"
DET_FILES = {0.001: "tcw_t0000.hfield", 0.01: "tcw_t0001.hfield", 0.1: "tcw_t0002.hfield"}

MAIN_K_RANGE = (0.0, 10.0)
ZOOM_K_RANGE = (0.0, 8.0)
ZOOM_K_STEP = 2.0
MAIN_FIGSIZE = (16.0, 9.0)
ZOOM_FIGSIZE = (7.5, 4.5)
MAIN_TICK_SIZE, MAIN_LABEL_SIZE, MAIN_LEGEND_SIZE = 32, 38, 32
ZOOM_TICK_SIZE, ZOOM_LABEL_SIZE = 36, 42
LINE_WIDTH = 3.4
DATA_LINE_WIDTH = 2.4
MARKER_SIZE_O, MARKER_EDGE = 11, 1.8
ZOOM_DECADE = 0.01
ZOOM_Y_STEP = 2.0 * ZOOM_DECADE
JITTER_MAX = 0.0
JITTER_SEED = 20260917

mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 24,
    "axes.labelsize": MAIN_LABEL_SIZE,
    "legend.fontsize": MAIN_LEGEND_SIZE,
})


def jitter_halfwidth(y: np.ndarray, y_max: float) -> np.ndarray:
    """Jitter grows with the value of the point (JITTER_MAX at the global peak)."""
    sigma = JITTER_MAX * y / y_max
    sigma[0] = 0.0
    return sigma


def style(axis, xlim, tick_size: float) -> None:
    axis.set(xlim=xlim)
    axis.tick_params(axis="both", which="both", direction="in", top=True, right=True,
                     labelsize=tick_size, length=9, width=1.6)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(JITTER_SEED)
    continuous = base.ring_spectra(base.L10, "repeat_run_*", base.L10_FILES)
    deterministic = base.ring_spectra(DET, "repeat_run_*", DET_FILES)
    k_cont0, _ = continuous[base.T_REF]
    k_det0, h0_det = deterministic[base.T_REF]
    h0_det = base.uniform_filter1d(h0_det, size=4, mode="nearest")
    y_max = max(float(continuous[t][1][::3].max()) for t in base.TIMES)

    fig_main, ax_main = plt.subplots(figsize=MAIN_FIGSIZE, constrained_layout=True)
    fig_zoom, ax_zoom = plt.subplots(figsize=ZOOM_FIGSIZE, constrained_layout=True)

    rows, zoom_values, main_values = [], [], []
    for time in base.TIMES:
        tau = time - base.T_REF
        k_cont, y_cont = continuous[time]
        k_det, y_det = deterministic[time]
        cont_theory = base.pure_driven_spectrum(k_cont0, time)
        det_theory = base.theory_from_h0(h0_det, k_det0, tau, driven=False)
        color = base.COLORS[time]

        k_show = k_cont[::3]
        y_raw = y_cont[::3]
        y_show = y_raw + rng.uniform(-1.0, 1.0, size=k_show.size) * jitter_halfwidth(y_raw, y_max)
        y_show = np.maximum(y_show, 0.0)
        ax_main.plot(k_show, y_show, color=color, linestyle="-.", linewidth=DATA_LINE_WIDTH,
                     marker="o", markersize=MARKER_SIZE_O, markerfacecolor="none",
                     markeredgecolor=color, markeredgewidth=MARKER_EDGE)
        ax_main.plot(k_cont0, cont_theory, color=color, linewidth=LINE_WIDTH)
        ax_main.plot(k_det0, det_theory, color=color, linestyle="--", linewidth=LINE_WIDTH)

        ax_zoom.plot(k_det0, det_theory, color=color, linestyle="--", linewidth=LINE_WIDTH)

        main_values.extend(y_show)
        main_values.extend(cont_theory[k_cont0 <= MAIN_K_RANGE[1] + 1e-9])
        zoom_values.extend(det_theory[k_det0 <= ZOOM_K_RANGE[1] + 1e-9])

        rows.extend((time, tau, k, y, "L10_continuous_numeric") for k, y in zip(k_cont, y_cont))
        rows.extend((time, tau, k, y, "initial_noise_numeric") for k, y in zip(k_det, y_det))
        rows.extend((time, time, k, y, "pure_driven_h0_eq_0") for k, y in zip(k_cont0, cont_theory))
        rows.extend((time, tau, k, y, "initial_noise_theory_from_H0") for k, y in zip(k_det0, det_theory))

    style(ax_main, MAIN_K_RANGE, MAIN_TICK_SIZE)
    ax_main.set(ylim=(0.0, float(np.max(main_values)) * 1.10),
                xlabel=r"$|k|$", ylabel=r"$|H|_{\rm rms}$")
    ax_main.xaxis.label.set_size(MAIN_LABEL_SIZE)
    ax_main.yaxis.label.set_size(MAIN_LABEL_SIZE)

    style(ax_zoom, ZOOM_K_RANGE, ZOOM_TICK_SIZE)
    zoom_top = float(np.ceil(np.max(zoom_values) / ZOOM_Y_STEP)) * ZOOM_Y_STEP
    ax_zoom.set(ylim=(0.0, zoom_top + 0.5 * ZOOM_Y_STEP), xlabel=r"$|k|$")
    ax_zoom.xaxis.label.set_size(ZOOM_LABEL_SIZE)
    ax_zoom.xaxis.set_major_locator(MultipleLocator(ZOOM_K_STEP))
    ax_zoom.yaxis.set_major_locator(MultipleLocator(ZOOM_Y_STEP))
    ax_zoom.ticklabel_format(axis="y", style="sci", scilimits=(-2, -2), useMathText=True)
    ax_zoom.yaxis.get_offset_text().set_size(ZOOM_TICK_SIZE)

    handles = [
        Line2D([], [], color="black", linewidth=LINE_WIDTH,
               label=r"theory, $\varphi=0.001$, $H_0=0$"),
        Line2D([], [], color="black", linewidth=LINE_WIDTH, linestyle="--",
               label=r"theory, $\varphi=0$, from $H_0$"),
        Line2D([], [], color="black", marker="o", markersize=MARKER_SIZE_O, markerfacecolor="none",
               markeredgewidth=MARKER_EDGE, linestyle="-.", linewidth=DATA_LINE_WIDTH,
               label=r"SLE, $\varphi=0.001$"),
    ]
    ax_main.legend(handles=handles, loc="upper right", frameon=False, fontsize=MAIN_LEGEND_SIZE,
                   handlelength=3.0, labelspacing=0.35, borderaxespad=0.5)

    for fig, tag in ((fig_main, "pure_main16x9"), (fig_zoom, "pure_zoom5x3")):
        fig.savefig(OUT / f"tcw_H0_reference_{tag}.png", dpi=220)
        fig.savefig(OUT / f"tcw_H0_reference_{tag}.pdf")

    np.savetxt(OUT / "tcw_H0_reference_pure_data.csv", rows, delimiter=",", fmt="%s",
               header="physical_time,tau_from_H0,k,rms,series", comments="")


if __name__ == "__main__":
    main()
