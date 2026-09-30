"""Whole-hole cross-validation baselines for Cu ppm (standard library only).

Run after initial_qaqc.py: python scripts/baseline_cv.py
Uses 3D Euclidean distance in the source's local Cartesian coordinates. Results
are exploratory interpolation benchmarks, not a resource estimate.
"""
from __future__ import annotations

import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "drillholes.csv"
OUT = ROOT / "reports" / "baseline_cv.json"
K = 12
POWER = 2.0


def metrics(actual: list[float], predicted: list[float]) -> dict[str, float]:
    errors = [p - y for y, p in zip(actual, predicted)]
    return {
        "n": len(actual),
        "mae_ppm": statistics.fmean(abs(e) for e in errors),
        "rmse_ppm": math.sqrt(statistics.fmean(e * e for e in errors)),
        "mean_error_ppm": statistics.fmean(errors),
    }


def main() -> None:
    with DATA.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    required = {"HOLEID", "X", "Y", "Z", "Cu ppm"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit(f"Missing required fields: {sorted(required - set(rows[0]) if rows else required)}")
    samples = []
    for row in rows:
        try:
            samples.append((row["HOLEID"], tuple(float(row[c]) for c in ("X", "Y", "Z")), float(row["Cu ppm"])))
        except (ValueError, TypeError):
            continue
    groups: dict[str, list[int]] = defaultdict(list)
    for i, (hole, _, _) in enumerate(samples):
        groups[hole].append(i)

    truth: list[float] = []
    predictions = {"training_mean": [], "nearest_neighbor_3d": [], f"idw_p{POWER:g}_k{K}": []}
    for heldout in groups.values():
        held = set(heldout)
        train = [s for i, s in enumerate(samples) if i not in held]
        mean_grade = statistics.fmean(s[2] for s in train)
        for i in heldout:
            _, point, grade = samples[i]
            nearest = sorted((math.dist(point, xyz), value) for _, xyz, value in train)
            truth.append(grade)
            predictions["training_mean"].append(mean_grade)
            predictions["nearest_neighbor_3d"].append(nearest[0][1])
            neighbors = nearest[:K]
            exact = [v for d, v in neighbors if d == 0]
            if exact:
                estimate = statistics.fmean(exact)
            else:
                weights = [1.0 / (d ** POWER) for d, _ in neighbors]
                estimate = sum(w * v for w, (_, v) in zip(weights, neighbors)) / sum(weights)
            predictions[f"idw_p{POWER:g}_k{K}"].append(estimate)

    report = {
        "target": "Cu ppm",
        "validation": "Leave-one-HOLEID-out; all samples from each held-out hole are excluded from training.",
        "distance": "3D Euclidean on X/Y/Z in source local Cartesian coordinates; no CRS transformation.",
        "idw": {"power": POWER, "neighbors": K},
        "n_samples": len(samples),
        "n_holes": len(groups),
        "methods": {name: metrics(truth, pred) for name, pred in predictions.items()},
        "caveats": [
            "Holes may be spatially interspersed; this grouped split does not guarantee spatially separated test areas.",
            "Distance isotropy is only a baseline assumption; directional continuity and variography remain unmodeled.",
            "Metrics do not support resource classification, tonnage, or economic claims.",
        ],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"Report saved: {OUT}")


if __name__ == "__main__":
    main()
