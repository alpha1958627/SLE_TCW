#!/usr/bin/env python3
"""Deterministic Matplotlib house style for academic figures."""

from __future__ import annotations

import argparse
import shutil
import tempfile
from pathlib import Path

import matplotlib as mpl


FIGSIZE = (8.0, 6.0)
FONT_SIZE = 20.0
COLORS = ("#000000", "#D62728", "#1F77B4")


def require_latex() -> None:
    """Require the executables used by Matplotlib's LaTeX rendering path."""
    missing = [name for name in ("latex", "dvipng", "dvips") if shutil.which(name) is None]
    if missing:
        raise RuntimeError("Missing LaTeX rendering tools: " + ", ".join(missing))


def apply_academic_style() -> None:
    """Apply the user's non-negotiable academic Matplotlib defaults."""
    require_latex()
    mpl.rcParams.update(
        {
            "figure.figsize": FIGSIZE,
            "figure.constrained_layout.use": False,
            "figure.facecolor": "white",
            "font.family": "serif",
            "font.size": FONT_SIZE,
            "text.usetex": True,
            "text.latex.preamble": r"\usepackage{amsmath}",
            "axes.titlesize": FONT_SIZE,
            "axes.labelsize": FONT_SIZE,
            "axes.linewidth": 1.2,
            "axes.grid": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.prop_cycle": mpl.cycler(color=COLORS),
            "xtick.labelsize": FONT_SIZE,
            "ytick.labelsize": FONT_SIZE,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.size": 6.0,
            "ytick.major.size": 6.0,
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
            "legend.fontsize": FONT_SIZE,
            "legend.frameon": False,
            "figure.titlesize": FONT_SIZE,
            "lines.linewidth": 2.0,
            "lines.markersize": 8.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
            "savefig.facecolor": "white",
            "savefig.bbox": None,
        }
    )


def validate_active_style() -> None:
    """Raise if the active rcParams no longer match the house style."""
    checks = {
        "figure.figsize": tuple(float(value) for value in mpl.rcParams["figure.figsize"])
        == FIGSIZE,
        "text.usetex": bool(mpl.rcParams["text.usetex"]),
        "font.size": float(mpl.rcParams["font.size"]) == FONT_SIZE,
        "axes.labelsize": float(mpl.rcParams["axes.labelsize"]) == FONT_SIZE,
        "axes.titlesize": float(mpl.rcParams["axes.titlesize"]) == FONT_SIZE,
        "xtick.labelsize": float(mpl.rcParams["xtick.labelsize"]) == FONT_SIZE,
        "ytick.labelsize": float(mpl.rcParams["ytick.labelsize"]) == FONT_SIZE,
        "legend.fontsize": float(mpl.rcParams["legend.fontsize"]) == FONT_SIZE,
        "axes.grid": not bool(mpl.rcParams["axes.grid"]),
        "axes.prop_cycle": tuple(mpl.rcParams["axes.prop_cycle"].by_key()["color"][:3])
        == COLORS,
    }
    failed = [name for name, ok in checks.items() if not ok]
    if failed:
        raise AssertionError("Academic style validation failed: " + ", ".join(failed))


def smoke_test() -> None:
    """Render temporary PNG/PDF files to exercise the full LaTeX pipeline."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    apply_academic_style()
    validate_active_style()
    x = np.linspace(0.0, 2.0 * np.pi, 80)
    fig, ax = plt.subplots()
    for index, (phase, style) in enumerate(zip((0.0, 0.4, 0.8), ("-", "--", ":"))):
        ax.plot(x, np.sin(x + phase), linestyle=style, label=rf"$s_{index + 1}$")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$f(x)$")
    ax.legend()
    ax.grid(False, which="both")
    with tempfile.TemporaryDirectory() as directory:
        base = Path(directory) / "academic-style-smoke"
        fig.savefig(base.with_suffix(".png"), dpi=600)
        fig.savefig(base.with_suffix(".pdf"))
        if not all(base.with_suffix(suffix).stat().st_size > 0 for suffix in (".png", ".pdf")):
            raise RuntimeError("Academic style smoke-test export failed")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="run the LaTeX render smoke test")
    args = parser.parse_args()
    if not args.check:
        parser.error("use --check")
    smoke_test()
    print("Academic Matplotlib style: OK")


if __name__ == "__main__":
    main()
