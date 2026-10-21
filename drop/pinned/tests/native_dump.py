"""Integrate native Cartesian 2-D, single-process Basilisk 170901 dumps.

Uses actual AMR leaf areas, not interpolated uniform-grid samples.
Run: python read_native_volume_20261011.py DUMP [DUMP ...]
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct

HEADER = struct.Struct("<dqiiii3d")  # Linux x86-64 DumpHeader, sizeof=56.


def parse_dump(raw):
    if len(raw) < HEADER.size:
        raise ValueError("Truncated dump header")
    time, count, step, depth, npe, version, nx, ny, nz = HEADER.unpack_from(raw)
    if version != 170901 or npe != 1 or (nx, ny, nz) != (1., 1., 0.):
        raise ValueError("Unsupported dump version, process count, or root geometry")
    if not (1 <= count <= 256 and 0 <= depth <= 16 and math.isfinite(time)):
        raise ValueError("Invalid dump header")
    offset, names = HEADER.size, []
    for _ in range(count):
        if offset + 4 > len(raw):
            raise ValueError("Truncated field name")
        size, = struct.unpack_from("<I", raw, offset)
        offset += 4
        if not 1 <= size <= 256 or offset + size > len(raw):
            raise ValueError("Invalid field name length")
        names.append(raw[offset:offset + size].decode("utf-8"))
        offset += size
    if names[0] != "size" or "h" not in names or "cm" in names:
        raise ValueError("Expected size/h fields and Cartesian unit metric")
    if offset + 32 > len(raw):
        raise ValueError("Truncated origin")
    x0, y0, z0, length = struct.unpack_from("<4d", raw, offset)
    offset += 32
    if not all(map(math.isfinite, (x0, y0, z0, length))) or length <= 0:
        raise ValueError("Invalid domain")
    record = struct.Struct("<I" + "d" * count)
    if (len(raw) - offset) % record.size:
        raise ValueError("Incomplete node record")
    levels, areas, masses, heights = [0], [], [], []
    h_index, nodes = names.index("h") + 1, 0
    for row in record.iter_unpack(memoryview(raw)[offset:]):
        if not levels:
            raise ValueError("Trailing tree nodes")
        level, flags = levels.pop(), row[0]
        if level > depth or flags not in (0, 2):
            raise ValueError("Invalid quadtree depth or leaf flag")
        nodes += 1
        if flags == 0:
            levels.extend([level + 1] * 4)
        else:
            height = row[h_index]
            if not math.isfinite(height) or height <= 0:
                raise ValueError("Nonfinite/nonpositive leaf height")
            area = (length / (1 << level)) ** 2
            areas.append(area)
            masses.append(height * area)
            heights.append(height)
    if levels or abs(math.fsum(areas) - length ** 2) > 1e-12 * length ** 2:
        raise ValueError("Incomplete tree or wrong covered area")
    return dict(time=time, step=step, max_level=depth, L0=length, fields=names,
                nodes=nodes, leaves=len(areas), area=math.fsum(areas),
                V_native=math.fsum(masses), hmin=min(heights), hmax=max(heights))


def self_check():
    raw = HEADER.pack(2., 2, 10, 1, 1, 170901, 1., 1., 0.)
    for name in (b"size", b"h"):
        raw += struct.pack("<I", len(name)) + name
    raw += struct.pack("<4d", 0, 0, 0, 2)
    raw += struct.pack("<I2d", 0, 5, 2.5)
    for height in (1., 2., 3., 4.):
        raw += struct.pack("<I2d", 2, 1, height)
    result = parse_dump(raw)
    assert result["V_native"] == 10 and result["area"] == 4 and result["leaves"] == 4
    try:
        parse_dump(raw[:-1])
    except ValueError:
        pass
    else:
        raise AssertionError("Truncated dump accepted")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dumps", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=Path("native_volume_audit_20261011.json"))
    args = parser.parse_args()
    self_check()
    results = []
    for path in args.dumps:
        raw = path.read_bytes()
        result = parse_dump(raw)
        result.update(file=str(path.resolve()), sha256=hashlib.sha256(raw).hexdigest())
        results.append(result)
    args.output.write_text(json.dumps(dict(self_check="PASS", dumps=results), indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))
