"""Independent native/export/AMR/real-step audit. No physical steady-state claim."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import struct
from native_dump import parse_dump, self_check


def meta(path):
    return dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def field(path):
    raw = path.read_bytes()
    n, = struct.unpack_from("<i", raw)
    assert 2 <= n <= 4096 and len(raw) == 4+8*n*n
    vals = struct.unpack_from("<%dd" % (n*n), raw, 4)
    assert all(math.isfinite(h) and h > 0 for h in vals)
    return n, vals


def rows(path):
    with path.open() as stream:
        first = stream.readline().strip()
        stream.seek(0)
        if path.name == 'volume.csv' and first and not first.startswith('step,'):
            # The legacy restore branch appends numeric rows without a CSV header.
            # Parse its documented schema without changing the raw file.
            names = ('step,time,V_total,V_excess,rel_drift,max_rel_drift,picard_iter,'
                     'picard_res,mg_i,mg_resa_after,mg_resb_before').split(',')
            result = list(csv.DictReader(stream, fieldnames=names))
            assert all(None not in row and all(v is not None for v in row.values()) for row in result)
            return result
        return list(csv.DictReader(stream))


def validate(case, steps, end):
    self_check()
    assert (case/"solver_exit.txt").read_text().strip() == "exit_code=0"
    err = (case/"stderr.log").read_text()
    assert not any(s in err for s in ("FATAL:", "SOLVER_FAILURE", "convergence for"))
    assert not (case/"native_invalid").exists()
    resume = (case/"checkpoint_input").exists()
    start_path = case/("checkpoint_input" if resume else "native_initial")
    initial, final = (parse_dump(p.read_bytes()) for p in (start_path, case/"native_final"))
    start_step = initial["step"]
    assert final["step"] == start_step+steps, (initial["step"], final["step"], steps)
    assert abs(final["time"]-end) < 1e-8
    # Diagnostics use the REAL coupled invocation index, including restored first step.
    indices = list(map(int, re.findall(r"HP001_COUPLED step=(\d+)", err)))
    assert indices == list(range(start_step, start_step+steps)), indices
    start_tag = "resume" if resume else "t0000"
    export_errors = []
    for tag, native in ((start_tag, initial), ("final", final)):
        smeta = meta(case/("snapshot_"+tag+".meta"))
        assert abs(float(smeta["time"])-native["time"]) < 1e-8
        assert int(smeta["step"]) == native["step"]
        assert float(smeta["noise"]) == 0 and float(smeta["hstar"]) == .01
        assert float(smeta["n_repulsive"]) == 5 and float(smeta["m_attractive"]) == 3
        assert float(smeta["gamma"]) == 1
        n, vals = field(case/("hfield_"+tag+".bin"))
        export_errors.append(abs(math.fsum(vals)*(32/n)**2-native["V_native"])/native["V_native"])
    assert max(export_errors) < 1e-12
    state = meta(case/"volume_state.txt")
    v0 = float(state["V0"])
    volumes = rows(case/"volume.csv")
    amr = rows(case/"amr_native.csv") if (case/"amr_native.csv").exists() else []
    drift = max([float(state["max_rel_volume_drift"]),
                 abs(final["V_native"]-v0)/v0] +
                [abs(float(r["V_total"])-v0)/v0 for r in volumes] +
                [abs(float(r["V_after"])-v0)/v0 for r in amr])
    transfer = max([0.] + [float(r["rel_transfer"]) for r in amr])
    assert drift < 1e-8 and transfer < 1e-10, (drift, transfer)
    for r in volumes:
        for key in ("picard_res", "mg_resa_after"):
            assert math.isfinite(float(r[key])) and abs(float(r[key])) <= 1e-8
    physical = [list(map(float, m)) for m in re.findall(
        r"HP001_COUPLED step=\d+ picard=\d+ dh=([^ ]+) mass=([^ ]+) pressure=([^ ]+)", err)]
    assert len(physical) == steps
    max_residual = max(max(r) for r in physical)
    assert max_residual <= 1e-8
    pressure = rows(case/"pressure_native.csv")[-1]
    return dict(status="PASS", stable_angle=False, real_coupled_steps=len(indices),
                start_step=start_step, final_step=final["step"], start_time=initial["time"],
                final_time=final["time"], V0=v0, V_final_native=final["V_native"],
                max_cumulative_volume_drift=drift, max_amr_transfer=transfer,
                max_native_export_error=max(export_errors), max_solver_residual=max_residual,
                final_native_finest=final["max_level"], final_leaf_count=final["leaves"],
                final_hmin=final["hmin"], final_hmax=final["hmax"],
                final_pressure=pressure,
                source_sha256=hashlib.sha256((case/"contact_angle.c").read_bytes()).hexdigest(),
                executable_sha256=hashlib.sha256((case/"solver").read_bytes()).hexdigest(),
                checkpoint_sha256=hashlib.sha256(start_path.read_bytes()).hexdigest())


def compare(first, second, transpose=False):
    n, a = field(first/"hfield_final.bin")
    nb, b = field(second/"hfield_final.bin")
    assert n == nb
    delta = max(abs(a[y*n+x]-b[x*n+y if transpose else y*n+x])
                for y in range(n) for x in range(n))
    peak = max(a)-.01
    return dict(max_height_difference=delta, relative_to_peak_excess=delta/peak if peak else 0,
                transposed=transpose)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--end", type=float, required=True)
    args = parser.parse_args()
    try:
        result = validate(args.case, args.steps, args.end)
    except Exception as exc:
        result = dict(status="FAIL", reason=repr(exc), stable_angle=False)
    (args.case/"verification.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    raise SystemExit(0 if result["status"] == "PASS" else 1)
