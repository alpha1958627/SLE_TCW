# Capillary midpoint experiment

Independent experimental comparison, based on SLE_TCW main cff0e3b and the locally validated Gaussian/time-bin OU–SE implementation. It does not replace the default solver.

`baseline.c` uses the existing lagged-mobility implicit Euler step. `midpoint.c` evaluates the capillary pressure and h^3 mobility at the implicit time midpoint by Picard iteration. Both use identical original old-state stochastic face-flux construction; the same 10 noise/amplitude/limiter functions are unchanged. This is **not** a stochastic midpoint/Stratonovich conversion and is not a claim of second-order stochastic accuracy.

Both sources contain a corrected selftest threshold: the configured `white_lower`, rather than precursor height `hp`, is compared against synthetic strong-forcing minima. No physical lower bound or original absolute rounding tolerance is changed.

## Build

Use Basilisk qcc and the existing MPI/graphics link libraries:

```sh
export CC99=mpicc
qcc -disable-dimensions -D_GNU_SOURCE -D_MPI=1 -O2 baseline.c -o baseline_mpi -L"$BASILISK/gl" -lglutils -lfb_tiny -lm
qcc -disable-dimensions -D_GNU_SOURCE -D_MPI=1 -O2 midpoint.c -o midpoint_mpi -L"$BASILISK/gl" -lglutils -lfb_tiny -lm
```

Run each binary with `--operator-selftest`, `--nonfinite-selftest`, `--face-index-selftest`, `--implicit-selftest`, and `--coupled-selftest` before simulation. The Fourier oracle in each source matches its own time integrator.

Example OU short test (256², L64, dt=.003125, t1, phi=.003, sigma=.05, tau=.01):

```sh
mpiexec -n 8 ./midpoint_mpi 102901 .003125 320 .01 2 320 0 0 .05 8 0 1 8 1e-8 1e-8 60 0 1 200 1 64 .003 0 1
```

For the deterministic control set mode to -1. White mode is 0, uses seed103001 and additionally accepts `1e-6 0 -` for lower-bound/checkpoint arguments. Keep the actual physical core count, output cadence, tolerances, and physical parameters matched when timing. Timing replays are not independent physical samples; different native dt at the same seed is not automatically the same random path.

AMR is disabled. The fixed FFT sampling lattice and fourth-order custom smoother are not yet validated for conservative/statistically consistent coarse-fine coupling. No raw fields or private trajectory tables belong in this code package.
