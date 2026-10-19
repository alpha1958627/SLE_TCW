"""Independent geometric contact-angle audit for the thin-film pilot fields.

This script deliberately separates four measurements:

* the existing h=a+b*r^2 extrapolation;
* the same parabola fit using only the centerline profile;
* a one-dimensional circular-arc (sphere-cap section) fit;
* the largest finite-difference slope, reported only as a local diagnostic.

The field/profile files are output data, so no solver or noise setting is
modified.  Angles are transient geometric estimates; this script does not
declare equilibrium.
"""
from __future__ import annotations

import argparse
import array
import csv
import json
import math
import re
import struct
from pathlib import Path
from typing import Iterable, Sequence


WINDOWS = ((0.2, 0.6), (0.4, 0.8), (0.6, 0.9))


def read_field(path: Path) -> tuple[int, list[float]]:
    raw = path.read_bytes()
    if len(raw) < 4:
        raise ValueError(f"incomplete hfield: {path}")
    n = struct.unpack_from("<i", raw, 0)[0]
    if not 2 <= n <= 4096 or len(raw) != 4 + 8 * n * n:
        raise ValueError(f"invalid hfield dimensions: {path}")
    vals = array.array("d")
    vals.frombytes(raw[4:])
    if vals.itemsize != 8:
        raise ValueError("unexpected host double size")
    if not all(math.isfinite(v) and v > 0.0 for v in vals):
        raise ValueError(f"nonpositive/nonfinite hfield: {path}")
    return n, list(vals)


def read_meta(path: Path) -> dict[str, float]:
    result: dict[str, float] = {}
    if not path.exists():
        return result
    for key, value in re.findall(
        r"([A-Za-z_][A-Za-z_0-9]*)\s*(?:=|:)\s*([-+0-9.eE]+)",
        path.read_text(encoding="utf-8"),
    ):
        result[key] = float(value)
    return result


def read_profile(path: Path) -> tuple[list[float], list[float], list[float]]:
    s: list[float] = []
    hx: list[float] = []
    hy: list[float] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"malformed profile line in {path}: {line!r}")
        a, b, c = map(float, fields)
        if not all(math.isfinite(v) for v in (a, b, c)):
            raise ValueError(f"nonfinite profile in {path}")
        s.append(a)
        hx.append(b)
        hy.append(c)
    if len(s) < 8 or any(b <= 0.0 or c <= 0.0 for b, c in zip(hx, hy)):
        raise ValueError(f"profile is empty or invalid: {path}")
    if any(b >= a for a, b in zip(s[1:], s[:-1])):
        raise ValueError(f"profile coordinates are not increasing: {path}")
    return s, hx, hy


def linear_fit(x: Sequence[float], y: Sequence[float]) -> tuple[float, float, float]:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("at least two paired points are required")
    xm = math.fsum(x) / len(x)
    ym = math.fsum(y) / len(y)
    den = math.fsum((v - xm) ** 2 for v in x)
    if den <= 0.0:
        raise ValueError("singular linear fit")
    slope = math.fsum((a - xm) * (b - ym) for a, b in zip(x, y)) / den
    intercept = ym - slope * xm
    rms = math.sqrt(math.fsum((b - (intercept + slope * a)) ** 2 for a, b in zip(x, y)) / len(x))
    return intercept, slope, rms


def parabola_angle(points: Sequence[tuple[float, float]], hstar: float) -> dict[str, float | int | None]:
    """Fit h=a+b*r^2 and extrapolate the slope at h=hstar."""
    r2 = [r * r for r, _ in points]
    h = [v for _, v in points]
    a, b, rms = linear_fit(r2, h)
    if b >= 0.0 or a <= hstar:
        return {"status": "invalid_extrapolation", "points": len(points), "rms_height_residual": rms}
    rc = math.sqrt((hstar - a) / b)
    slope = -2.0 * b * rc
    return {
        "status": "ok",
        "points": len(points),
        "intercept": a,
        "quadratic_coefficient": b,
        "contact_radius": rc,
        "slope": slope,
        "geometric_angle_deg": math.degrees(math.atan(abs(slope))),
        "small_slope_parameter_deg": math.degrees(abs(slope)),
        "rms_height_residual": rms,
    }


def solve_2x2(a11: float, a12: float, a21: float, a22: float, b1: float, b2: float) -> tuple[float, float]:
    det = a11 * a22 - a12 * a21
    if abs(det) <= 1e-18:
        raise ValueError("singular circle fit")
    return ((b1 * a22 - a12 * b2) / det, (a11 * b2 - b1 * a21) / det)


def circular_arc_angle(points: Sequence[tuple[float, float]], hstar: float) -> dict[str, float | int | None]:
    """Fit r^2+h^2-2*c*h+F=0, a circle centred on the symmetry axis."""
    if len(points) < 3:
        return {"status": "too_few_points", "points": len(points)}
    # 2*c*h - F = r^2+h^2.  This is a linear least-squares fit for c,F.
    h = [v for _, v in points]
    rhs = [r * r + v * v for r, v in points]
    hmean = math.fsum(h) / len(h)
    rmean = math.fsum(rhs) / len(rhs)
    den = math.fsum((v - hmean) ** 2 for v in h)
    if den <= 1e-24:
        return {"status": "singular_circle_fit", "points": len(points)}
    # Regress rhs = p*h + q; p=2c and q=-F.
    p = math.fsum((v - hmean) * (z - rmean) for v, z in zip(h, rhs)) / den
    q = rmean - p * hmean
    c = 0.5 * p
    F = -q
    radius2 = c * c - F
    if radius2 <= 0.0 or not math.isfinite(radius2):
        return {"status": "invalid_circle", "points": len(points)}
    radius = math.sqrt(radius2)
    residuals: list[float] = []
    for r, value in points:
        radicand = radius2 - r * r
        if radicand < -1e-10:
            return {"status": "circle_domain_failure", "points": len(points)}
        predicted = c + math.sqrt(max(0.0, radicand))
        residuals.append(value - predicted)
    rms = math.sqrt(math.fsum(v * v for v in residuals) / len(residuals))
    base_cos = (hstar - c) / radius
    if not -1.0 <= base_cos <= 1.0:
        return {
            "status": "base_outside_circle",
            "points": len(points),
            "circle_center_height": c,
            "circle_radius": radius,
            "rms_height_residual": rms,
            "base_cosine_unclipped": base_cos,
        }
    rc = radius * math.sqrt(max(0.0, 1.0 - base_cos * base_cos))
    slope = math.inf if hstar == c else rc / abs(hstar - c)
    return {
        "status": "ok",
        "points": len(points),
        "circle_center_height": c,
        "circle_radius": radius,
        "contact_radius": rc,
        "slope": slope,
        "geometric_angle_deg": math.degrees(math.acos(max(-1.0, min(1.0, base_cos)))),
        "small_slope_parameter_deg": math.degrees(slope),
        "rms_height_residual": rms,
        "base_cosine": base_cos,
    }


def moving_average(values: Sequence[float], radius: int = 2) -> list[float]:
    out: list[float] = []
    for i in range(len(values)):
        lo = max(0, i - radius)
        hi = min(len(values), i + radius + 1)
        out.append(math.fsum(values[lo:hi]) / (hi - lo))
    return out


def local_max_slope(r: Sequence[float], h: Sequence[float], hstar: float) -> dict[str, float | int | None]:
    """Return the strongest resolved local slope; not a macroscopic angle."""
    if len(r) < 5:
        return {"status": "too_few_points"}
    smooth = moving_average(h, 2)
    candidates: list[tuple[float, int, float]] = []
    peak = max(h)
    for i in range(1, len(r) - 1):
        if r[i] <= 0.0 or h[i] <= hstar + 1e-5 or h[i] >= peak - 1e-10:
            continue
        dr = r[i + 1] - r[i - 1]
        if dr <= 0.0:
            continue
        slope = (smooth[i + 1] - smooth[i - 1]) / dr
        candidates.append((abs(slope), i, slope))
    if not candidates:
        return {"status": "no_candidate"}
    magnitude, i, signed = max(candidates)
    return {
        "status": "diagnostic_only",
        "index": i,
        "radius": r[i],
        "height": h[i],
        "signed_slope": signed,
        "slope": magnitude,
        "local_angle_deg": math.degrees(math.atan(magnitude)),
        "small_slope_parameter_deg": math.degrees(magnitude),
        "smoothing_halfwidth_points": 2,
    }


def profile_symmetry(s: Sequence[float], hx: Sequence[float], hy: Sequence[float], hstar: float) -> dict[str, float | int]:
    n = len(s)
    pairs = []
    for i in range(n // 2):
        j = n - 1 - i
        pairs.append((hx[i], hx[j], hy[i], hy[j]))
    all_values = list(hx) + list(hy)
    scale = max(max(all_values) - hstar, 1e-15)
    lr = [abs(a - b) for a, b, _, _ in pairs]
    xy = [abs((a + b) * 0.5 - (c + d) * 0.5) for a, b, c, d in pairs]
    return {
        "pairs": len(pairs),
        "max_left_right_abs": max(lr),
        "rms_left_right_abs": math.sqrt(math.fsum(v * v for v in lr) / len(lr)),
        "max_horizontal_vertical_abs": max(xy),
        "rms_horizontal_vertical_abs": math.sqrt(math.fsum(v * v for v in xy) / len(xy)),
        "max_left_right_relative_to_excess": max(lr) / scale,
        "max_horizontal_vertical_relative_to_excess": max(xy) / scale,
    }


def field_centerline_symmetry(n: int, values: Sequence[float], hstar: float) -> dict[str, float | int]:
    """Check the exported field on the two central rows and columns."""
    mid0, mid1 = n // 2 - 1, n // 2
    horizontal = []
    vertical = []
    for i in range(n):
        horizontal.append(0.5 * (values[mid0 * n + i] + values[mid1 * n + i]))
        vertical.append(0.5 * (values[i * n + mid0] + values[i * n + mid1]))
    lr = [abs(horizontal[i] - horizontal[n - 1 - i]) for i in range(n // 2)]
    xy = [abs(horizontal[i] - vertical[i]) for i in range(n)]
    scale = max(max(values) - hstar, 1e-15)
    return {
        "N": n,
        "max_left_right_abs": max(lr),
        "rms_left_right_abs": math.sqrt(math.fsum(v * v for v in lr) / len(lr)),
        "max_horizontal_vertical_abs": max(xy),
        "rms_horizontal_vertical_abs": math.sqrt(math.fsum(v * v for v in xy) / len(xy)),
        "max_left_right_relative_to_excess": max(lr) / scale,
        "max_horizontal_vertical_relative_to_excess": max(xy) / scale,
    }


def analyze_one(profile_path: Path, field_path: Path, meta_path: Path) -> dict:
    meta = read_meta(meta_path)
    s, hx, hy = read_profile(profile_path)
    hstar = meta.get("hstar", 0.1)
    # Positive radius contains one radial half of each centerline.  Averaging
    # x/y suppresses the tiny orientation difference without changing geometry.
    positive = [(abs(q), 0.5 * (a + b)) for q, a, b in zip(s, hx, hy) if q >= 0.0]
    positive.sort()
    r = [q for q, _ in positive]
    h = [v for _, v in positive]
    far = [v for q, v in positive if q > (meta.get("L0", 32.0) / 3.0)]
    far_height = math.fsum(far) / len(far)
    peak = max(h)
    windows: list[dict] = []
    for lo, hi in WINDOWS:
        selected = [(q, v) for q, v in positive if far_height + lo * (peak - far_height) < v < far_height + hi * (peak - far_height)]
        row = {
            "height_fraction_window": [lo, hi],
            "points": len(selected),
            "parabola_centerline": parabola_angle(selected, hstar) if len(selected) >= 4 else {"status": "too_few_points", "points": len(selected)},
            "circular_arc_centerline": circular_arc_angle(selected, hstar) if len(selected) >= 4 else {"status": "too_few_points", "points": len(selected)},
        }
        windows.append(row)
    n, field = read_field(field_path)
    return {
        "case": profile_path.parent.name,
        "time": meta.get("time"),
        "step": meta.get("step"),
        "N": n,
        "level": meta.get("level"),
        "k": meta.get("k"),
        "target_theta_deg": meta.get("theta_exact_degrees"),
        "profile_file": str(profile_path.resolve()),
        "field_file": str(field_path.resolve()),
        "symmetry_profile": profile_symmetry(s, hx, hy, hstar),
        "symmetry_field": field_centerline_symmetry(n, field, hstar),
        "farfield_height": far_height,
        "peak_height": peak,
        "windows": windows,
        "local_max_slope": local_max_slope(r, h, hstar),
        "profile_centerline_points": len(r),
    }


def discover_case(case: Path, wanted_time: float | None = None) -> dict:
    profiles = sorted(case.glob("profile_t*.dat"))
    records = []
    for profile in profiles:
        tag = profile.stem.removeprefix("profile_")
        field = case / f"hfield_{tag}.bin"
        meta = case / f"snapshot_{tag}.meta"
        if not field.exists() or not meta.exists():
            continue
        values = read_meta(meta)
        time = values.get("time")
        if wanted_time is not None and (time is None or abs(time - wanted_time) > 1e-8):
            continue
        records.append(analyze_one(profile, field, meta))
    if not records:
        raise ValueError(f"no complete profile/field/meta triplets in {case}")
    records.sort(key=lambda row: (row["time"] is None, row["time"] or 0.0))
    return {"case": case.name, "records": records}


def self_check() -> None:
    hstar = 0.1
    theta = math.radians(25.0)
    radius = 4.0
    center = hstar - radius * math.cos(theta)
    contact_r = radius * math.sin(theta)
    pts = []
    for i in range(1, 101):
        r = contact_r * i / 101.0
        pts.append((r, center + math.sqrt(radius * radius - r * r)))
    circle = circular_arc_angle(pts, hstar)
    assert circle["status"] == "ok"
    assert abs(float(circle["geometric_angle_deg"]) - 25.0) < 1e-8
    parabola = parabola_angle([(r, 1.0 - 0.1 * r * r) for r in (.5, 1.0, 1.5, 2.0)], hstar)
    assert parabola["status"] == "ok"
    assert abs(float(parabola["geometric_angle_deg"]) - math.degrees(math.atan(.6))) < 1e-10
    s = [-1.5, -0.5, 0.5, 1.5]
    symmetry = profile_symmetry(s, [1.0, .8, .8, 1.0], [1.0, .8, .8, 1.0], hstar)
    assert symmetry["max_left_right_abs"] == 0.0


def synthetic_calibration() -> dict:
    """Quantify method bias on exact spherical caps and exact parabolas."""
    caps = []
    for theta_deg in (15.0, 25.0, 35.0):
        hstar = 0.1
        radius = 4.0
        theta = math.radians(theta_deg)
        center = hstar - radius * math.cos(theta)
        contact_r = radius * math.sin(theta)
        radii = [contact_r * (i + 0.5) / 512.0 for i in range(512)]
        heights = [center + math.sqrt(radius * radius - r * r) for r in radii]
        windows = []
        peak = max(heights)
        for lo, hi in WINDOWS:
            selected = [(r, h) for r, h in zip(radii, heights) if hstar + lo * (peak - hstar) < h < hstar + hi * (peak - hstar)]
            windows.append({
                "height_fraction_window": [lo, hi],
                "parabola_angle_deg": parabola_angle(selected, hstar)["geometric_angle_deg"],
                "circle_angle_deg": circular_arc_angle(selected, hstar)["geometric_angle_deg"],
            })
        caps.append({"synthetic": "exact_spherical_cap", "target_angle_deg": theta_deg, "windows": windows})
    # This parabola is the analytic form assumed by the existing estimator.
    selected = [(r, 1.0 - 0.1 * r * r) for r in (.5, 1.0, 1.5, 2.0)]
    return {"exact_caps": caps, "exact_parabola": parabola_angle(selected, 0.1)}


def write_csv(result: dict, path: Path) -> None:
    rows = []
    for case in result["cases"]:
        for record in case["records"]:
            base = {key: record.get(key) for key in ("case", "time", "step", "N", "level", "k", "target_theta_deg", "farfield_height", "peak_height")}
            local = record["local_max_slope"]
            base.update({"local_max_status": local.get("status"), "local_max_radius": local.get("radius"), "local_max_height": local.get("height"), "local_max_angle_deg": local.get("local_angle_deg")})
            base.update({f"profile_{key}": value for key, value in record["symmetry_profile"].items()})
            base.update({f"field_{key}": value for key, value in record["symmetry_field"].items()})
            for window in record["windows"]:
                row = dict(base)
                lo, hi = window["height_fraction_window"]
                row["window_lo"], row["window_hi"] = lo, hi
                row["points"] = window["points"]
                for prefix, fit in (("parabola", window["parabola_centerline"]), ("circle", window["circular_arc_centerline"])):
                    row[f"{prefix}_status"] = fit.get("status")
                    row[f"{prefix}_angle_deg"] = fit.get("geometric_angle_deg")
                    row[f"{prefix}_radius"] = fit.get("contact_radius")
                    row[f"{prefix}_rms_height"] = fit.get("rms_height_residual")
                rows.append(row)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dirs", nargs="+", type=Path)
    parser.add_argument("--time", type=float, default=50.0, help="analyze one output time; use --time -1 for all")
    parser.add_argument("--output", type=Path, default=Path("angle_measure_methods_20261011.json"))
    args = parser.parse_args()
    self_check()
    wanted = None if args.time < 0.0 else args.time
    cases = [discover_case(path, wanted) for path in args.case_dirs]
    result = {
        "self_check": "PASS",
        "purpose": "geometry-only contact-angle method audit; no equilibrium claim",
        "angle_definitions": {
            "parabola": "atan(|dh/dr|) from h=a+b*r^2 extrapolated to hstar",
            "circular_arc": "acos((hstar-c)/R) from an axis-centred circle fit",
            "local_max": "atan(max resolved |dh/dr|), diagnostic only, not the macroscopic angle",
        },
        "synthetic_calibration": synthetic_calibration(),
        "cases": cases,
    }
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    write_csv(result, args.output.with_suffix(".csv"))
    print(f"Self-check PASS; {sum(len(c['records']) for c in cases)} records; {args.output.resolve()}")


if __name__ == "__main__":
    main()
