"""Quantify ambiguity in spatial candidate links; never assert sample identity."""
import csv
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def xyz(row):
    return tuple(float(row[k]) for k in ("X", "Y", "Z"))


def main():
    drill = read(RAW / "drillholes.csv")
    by_hole = {}
    for i, row in enumerate(drill, 1):
        by_hole.setdefault(row["HOLEID"], []).append((i, xyz(row)))

    out = {}
    for filename in ("comminution.csv", "flotation.csv"):
        rows = read(RAW / filename)
        records = []
        for i, row in enumerate(rows, 1):
            point = xyz(row)
            same = sorted((math.dist(point, p), j) for j, p in by_hole.get(row["HOLEID"], []))
            any_hole = sorted((math.dist(point, xyz(d)), j, d["HOLEID"]) for j, d in enumerate(drill, 1))
            d1, j1 = same[0] if same else (None, None)
            d2, j2 = same[1] if len(same) > 1 else (None, None)
            records.append({
                "companion_row": i, "HOLEID": row["HOLEID"],
                "same_hole_candidate_count": len(same),
                "nearest_same_hole_row": j1, "nearest_same_hole_distance": d1,
                "second_same_hole_row": j2, "second_same_hole_distance": d2,
                "nearest_any_row": any_hole[0][1], "nearest_any_HOLEID": any_hole[0][2],
                "nearest_any_distance": any_hole[0][0],
                "same_hole_distance_margin": d2 - d1 if d2 is not None else None,
                "same_hole_distance_ratio": d2 / d1 if d2 is not None and d1 > 0 else None,
            })
        distances = [r["nearest_same_hole_distance"] for r in records if r["nearest_same_hole_distance"] is not None]
        margins = [r["same_hole_distance_margin"] for r in records if r["same_hole_distance_margin"] is not None]
        out[filename] = {
            "n": len(rows), "rows_with_matching_HOLEID": len(distances),
            "rows_without_matching_HOLEID": len(rows) - len(distances),
            "same_hole_nearest_distance_units": {
                "median": statistics.median(distances), "max": max(distances),
                "within_5": sum(d <= 5 for d in distances),
                "within_10": sum(d <= 10 for d in distances),
                "within_20": sum(d <= 20 for d in distances),
            },
            "second_candidate_separation_units": {
                "n": len(margins), "median_margin": statistics.median(margins),
                "min_margin": min(margins),
                "nearest_at_least_2x_closer": sum(r["same_hole_distance_ratio"] >= 2 for r in records if r["same_hole_distance_ratio"] is not None),
                "nearest_within_1_unit_of_second": sum(r["same_hole_distance_margin"] <= 1 for r in records if r["same_hole_distance_margin"] is not None),
            },
            "interpretation": "Coordinate ranking measures geometric uniqueness only; it cannot establish sample identity or interval support.",
            "records": records,
        }
    path = ROOT / "reports" / "companion_link_sensitivity.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: {a: b for a, b in v.items() if a != "records"} for k, v in out.items()}, indent=2))


if __name__ == "__main__":
    main()
