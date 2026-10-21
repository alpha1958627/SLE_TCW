# Independent q17 observer validation (2026-10-21)

This experiment prepares a separate source directory from the frozen q3 rupture
engine. It changes sampling observations, labels, and checkpoint identity only.
`physics.h` must retain SHA256
`1f108bce15a106991904c0d4744ef142b4fd6fc6bbf50e01796635aaa1d46290`.
Never run the preparation or a rebuild inside a running q3 campaign directory.

The score is

```
q17 = hmin + clip((hmin-.01)/.19, 0, 1) * (global_min_17x17_mean - hmin)
```

The global mean uses periodic indices, not an assumed 9-cell MPI halo. The
original event `hmin <= .01` is unchanged. RESTFW17 checkpoints reject q3 score
histories. New production trajectories must start from flat fields with new seeds.

From a compatible frozen checkout (baseline published at `1c5f20a`):

```
python prepare_variant.py --base /path/to/frozen/SLE_TCW --output /new/isolated/source
```

Copy `validate.slurm` to that source root. On a compatible Basilisk installation,
set `BASILISK` and `BASELINE_ROOT` (the validated q3 source root), then submit the
validation script from the isolated root. It builds actual serial/MPI binaries,
runs the full existing two-seed 2000-step reference/MPI1/4/16/repeat/restart and
controller gates, compares fields to the saved q3 baseline, checks static score
endpoints and periodic pits, and rejects q3 checkpoints. Synthetic static score
tests do not validate near-rupture integration or time-step convergence.

`check_d3_oracle.py` independently tests ideal AMS N12/D1 and N12/D3 with exactly
5000 toy repetitions each, retaining the terminal event fraction. This does not
replace an actual D3 startup smoke test for a subsequent production job.

The next parameter plan is exploratory: same q17/physics/event/T12000 across
AMS/DMC/QDMC; N12, AMS D3/max16, DMC terminal strength14, QDMC .5, block150.
DMC14 targets potential SD about .4 using the 12 completed independent direct
baseline terminal fields (q17 SD .02814). These are neither optimal nor exactly
equivalent strengths. D3 changes the budget and concurrent work as well as the
score, so improvements cannot all be attributed to the score or framework.
Numerical gate30902 passed (paired baseline field difference0, score error
2.23e-16, MPI error8.66e-15). Independent candidate tests were submitted as
AMS30927/DMC30929/QDMC30931. Each runs the actual new-controller startup checks
before its long-window production. This establishes no physical lifetime.

The immutable parameter contract retains its preparation-time status field.
Actual deployment status comes only from `submission_production.json` and Slurm;
do not infer it from the plan's historical `status` value. `submit_candidate.py`
uses a single-use claim and never replays unknown attempts. Its scripts belong
in `admin/` beside `source/` and `plan.json`; `run_candidate.py` accepts that
experiment root. `monitor_candidate.py` reads one source-checked snapshot.

All finite-window P estimates must retain inverse-potential compensation and
inherited first passage. Unresolved AMS W is not P. No MFPT, acceleration factor,
or efficiency winner is established by this numerical validation.
