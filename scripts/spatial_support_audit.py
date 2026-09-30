"""Summarize sample support and nearest-neighbor coverage for the Cu dataset."""
import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "drillholes.csv"
OUT = ROOT / "reports" / "spatial_support_audit.json"
TILE = 300.0


def quantiles(values):
    values = sorted(values)
    return {"min": values[0], "q1": statistics.quantiles(values, n=4)[0],
            "median": statistics.median(values), "q3": statistics.quantiles(values, n=4)[2],
            "max": values[-1]}


def main():
    with DATA.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    pts = [(r["HOLEID"], tuple(float(r[c]) for c in ("X", "Y", "Z")), float(r["Cu ppm"])) for r in rows]
    x0, y0 = min(p[1][0] for p in pts), min(p[1][1] for p in pts)
    tiles = Counter((math.floor((p[1][0]-x0)/TILE), math.floor((p[1][1]-y0)/TILE)) for p in pts)
    nearest_3d, nearest_xy, nearest_other_hole_3d = [], [], []
    for i, (hole, xyz, _) in enumerate(pts):
        others = [(math.dist(xyz, p[1]), math.hypot(xyz[0]-p[1][0], xyz[1]-p[1][1]), p[0]) for j,p in enumerate(pts) if i != j]
        nearest_3d.append(min(v[0] for v in others))
        nearest_xy.append(min(v[1] for v in others))
        other_hole = [v[0] for v in others if v[2] != hole]
        if other_hole:
            nearest_other_hole_3d.append(min(other_hole))
    per_hole = defaultdict(list)
    for hole, xyz, _ in pts:
        per_hole[hole].append(xyz)
    hole_depth = [max(p[2] for p in ps)-min(p[2] for p in ps) for ps in per_hole.values()]
    out = {
        "dataset": "drillholes.csv", "target": "Cu ppm", "rows": len(pts), "distinct_holes": len(per_hole),
        "coordinate_units": "Not documented by source; values retained in supplied local coordinates.",
        "extent": {axis: {"min": min(p[1][j] for p in pts), "max": max(p[1][j] for p in pts)} for j,axis in enumerate(("X","Y","Z"))},
        "nearest_sample_spacing_units": {"3d": quantiles(nearest_3d), "xy": quantiles(nearest_xy), "nearest_sample_from_different_hole_3d": quantiles(nearest_other_hole_3d)},
        "samples_per_hole": quantiles([len(ps) for ps in per_hole.values()]),
        "within_hole_z_span_units": quantiles(hole_depth),
        "xy_tile_occupancy": {"tile_size": TILE, "occupied_tiles": len(tiles), "empty_tiles_within_bounding_grid": (math.ceil((max(p[1][0] for p in pts)-x0)/TILE)*math.ceil((max(p[1][1] for p in pts)-y0)/TILE))-len(tiles), "samples_per_occupied_tile": quantiles(list(tiles.values())), "occupied_tile_counts": [{"ix": k[0], "iy": k[1], "n": v} for k,v in sorted(tiles.items())]},
        "interpretation": ["XY tile occupancy describes sampling footprint only and does not establish geological continuity.", "Nearest-sample distances measure local geometric support, not kriging uncertainty or resource confidence.", "No extrapolation beyond sampled support, tonnage, classification, or economic interpretation is supported without CRS, survey, interval, geology, density, and QA metadata."],
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k:v for k,v in out.items() if k != "xy_tile_occupancy"}, indent=2))
    print("XY tiles:", out["xy_tile_occupancy"]["occupied_tiles"], "occupied of", out["xy_tile_occupancy"]["occupied_tiles"]+out["xy_tile_occupancy"]["empty_tiles_within_bounding_grid"])


if __name__ == "__main__":
    main()
