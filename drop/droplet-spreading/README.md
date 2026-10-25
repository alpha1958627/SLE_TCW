# Droplet spreading with optional stochastic fluxes

This folder contains an English-documented, single-source droplet-spreading solver. The default mode is deterministic (`noise_mode = -1`). The same source exposes two optional modes: Gaussian white-in-time/white-in-space face noise (`0`) and periodic Ornstein–Uhlenbeck noise with a squared-exponential spatial spectrum (OU-SE, `2`). No hard height floor is applied by default.

## Model and initial condition

The code reuses the thin-film equation and the cos^4 droplet used in the OU-SE pilot:

```text
partial_t h = div(h^3 grad(p)),   p = -laplacian(h) - Pi(h)
h(x,y,0) = h_p + H_d cos^4(pi r/(2 R_0))  for r < R_0; h = h_p otherwise
```

Defaults for the droplet case are `L=32`, periodic boundaries, `h_p=0.01`, `R_0=5`, `H_d=1`, `A=0`, and `varphi=0.03` (used only when a noise mode is explicitly selected). The potential coefficients and discretization are those of the included OU-SE source; this is an exploratory solver, not a claim of long-time validation.

## Noise choices

The comments beside `TCW_NOISE_MODE_DEFAULT` and the command line show the three modes. The checked-in default is `-1` (no fluctuations).

| Mode | Meaning | Precursor film | Time scaling |
|---:|---|---|---|
| `-1` | Deterministic reference | No stochastic forcing | N/A |
| `0` | Independent Gaussian face noise each step | Exactly zero at `h <= h_p`; smooth ramp over `0.005` above `h_p` | White-noise `1/sqrt(dt)` scaling |
| `2` | Periodic OU-SE face fields, generated spectrally | Full-field, matching the prior OU-SE pilot | Eq. (52) by default (`legacy_scale=0`) |

White-noise precursor masking multiplies only the stochastic flux amplitude; it does not clip or raise the height. The OU-SE branch intentionally retains the pilot's full-field forcing and does not apply that white-noise mask. The numerical floor remains disabled (`height_floor=0`) so that any positivity or rupture behavior is visible rather than hidden by clipping.

## Build

Use the same Basilisk/qcc environment as the source snapshot. For example:

```bash
qcc -O2 -fopenmp -D_GNU_SOURCE droplet_spreading.c -o droplet_spreading \
  -L"$BASILISK/gl" -lglutils -lfb_tiny -lm
```

The source includes Basilisk `view.h`; if graphics libraries are unavailable on a compute node, use the project's established headless build configuration.

## Run examples

The command-line argument order is:

```text
seed dt steps tau mode sample_every noise_snapshot_every noise_probe_every
sigma OU_grid_level legacy_scale OU_amplitude_multiplier grid_level
MG_tolerance Picard_tolerance Picard_max height_floor varphi_scale
MG_nitermax droplet_case
```

All arguments should be supplied when selecting the droplet case (`droplet_case=1`, argument 20). Example runs to `t=1` on a uniform level-8 grid:

```bash
# Deterministic baseline (noise mode -1)
./droplet_spreading 314159 0.01 100 0.01 -1 10 0 10 0.5 8 0 1.0 8 1e-8 1e-8 30 0 1.0 100 1

# Gaussian white noise, switched off at the precursor film (mode 0)
./droplet_spreading 314159 0.01 100 0.01 0 10 0 10 0.5 8 0 1.0 8 1e-8 1e-8 30 0 1.0 100 1

# OU-SE noise with Eq. (52) scaling (mode 2)
./droplet_spreading 314159 0.01 100 0.01 2 10 0 10 0.5 8 0 1.0 8 1e-8 1e-8 30 0 1.0 100 1
```

These are short smoke examples only. For a longer physical time, increase `steps` consistently with the chosen `dt`; keep the same mesh, physical parameters, and seed when comparing noise modes. The source's `varphi` is set by `COEFF_VARPHI` at compile time (default `0.03`).

## Provenance and limitations

`droplet_spreading.c` is adapted from the OU-SE solver snapshot used for the local `t=1` pilot (source SHA-256 before this packaging: `E27533665A993F017FCEAA4602CADD4ADBA8722257B277CB55011977E024E514`). The OU-SE pilot reached `t=1` with finite positive height and small measured mass drift on a uniform level-8 grid; that does not demonstrate long-time stability or convergence. The white-noise precursor mask is a new amplitude taper and requires its own short-step, mass, timestep, and resolution checks before production use.

## Reusable plotting notebook

[`analysis/phi0001_20261025/drop.ipynb`](analysis/phi0001_20261025/drop.ipynb) packages the existing L64 low-noise droplet data, radius plots, raw/annular/ensemble profile comparisons, and native-cell zooms. All plotting code and required data are included; see its [README](analysis/phi0001_20261025/README.md) for definitions and the older deterministic baseline's different timestep. This post-processing package does not change or validate the solver defaults above.
