"""Check fixed spatial-CV baselines against XY grid-origin choices.

The repeated partitions are a design-sensitivity analysis, not independent
replicates or block/resource validation.
"""
from __future__ import annotations

import csv
import heapq
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/raw/drillholes.csv"
REPORT = ROOT / "reports/spatial_origin_sensitivity.json"
FIGURE = ROOT / "outputs/spatial_origin_sensitivity.svg"
TILE = 300.0
BUFFER = 50.0
OFFSETS = (0.0, TILE / 3.0, 2.0 * TILE / 3.0)
K = 12
POWER = 2.0


def metrics(actual: list[float], predicted: list[float]) -> dict[str, float | int]:
    error = [p - y for y, p in zip(actual, predicted)]
    return {
        "n": len(error),
        "mae_ppm": statistics.fmean(abs(e) for e in error),
        "rmse_ppm": math.sqrt(statistics.fmean(e * e for e in error)),
        "mean_error_ppm": statistics.fmean(error),
    }


def evaluate(samples, offset_x: float, offset_y: float):
    x_min = min(s[1][0] for s in samples)
    y_min = min(s[1][1] for s in samples)
    origin_x, origin_y = x_min + offset_x, y_min + offset_y
    tiles = defaultdict(list)
    for i, (_, point, _) in enumerate(samples):
        tile = (math.floor((point[0] - origin_x) / TILE),
                math.floor((point[1] - origin_y) / TILE))
        tiles[tile].append(i)

    actual, predicted_idw, predicted_mean, folds = [], [], [], []
    skipped = 0
    for (bx, by), test in sorted(tiles.items()):
        left, bottom = origin_x + bx * TILE, origin_y + by * TILE
        right, top = left + TILE, bottom + TILE
        test_set = set(test)
        train = []
        for i, (_, point, _) in enumerate(samples):
            if i in test_set:
                continue
            dx = max(left - point[0], 0.0, point[0] - right)
            dy = max(bottom - point[1], 0.0, point[1] - top)
            if math.hypot(dx, dy) > BUFFER:
                train.append(i)
        if len(train) < K:
            skipped += 1
            continue

        train_mean = statistics.fmean(samples[i][2] for i in train)
        fold_actual, fold_idw, fold_mean = [], [], []
        for i in test:
            point, grade = samples[i][1], samples[i][2]
            nearest = heapq.nsmallest(
                K,
                ((math.dist(point, samples[j][1]), j) for j in train),
            )
            if nearest[0][0] == 0:
                estimate = statistics.fmean(
                    samples[j][2] for distance, j in nearest if distance == 0
                )
            else:
                weights = [1.0 / (distance ** POWER) for distance, _ in nearest]
                estimate = math.fsum(
                    weight * samples[j][2]
                    for weight, (_, j) in zip(weights, nearest)
                ) / math.fsum(weights)
            actual.append(grade)
            predicted_idw.append(estimate)
            predicted_mean.append(train_mean)
            fold_actual.append(grade)
            fold_idw.append(estimate)
            fold_mean.append(train_mean)
        folds.append({
            "tile": [bx, by],
            "n_test": len(test),
            "n_train": len(train),
            "idw": metrics(fold_actual, fold_idw),
            "training_mean": metrics(fold_actual, fold_mean),
        })
    return {
        "origin_offset_xy_coordinate_units": [offset_x, offset_y],
        "tile_width_coordinate_units": TILE,
        "buffer_coordinate_units": BUFFER,
        "occupied_tiles": len(tiles),
        "evaluated_folds": len(folds),
        "skipped_folds": skipped,
        "evaluated_samples": len(actual),
        "coverage_percent": 100.0 * len(actual) / len(samples),
        "idw": metrics(actual, predicted_idw) if actual else None,
        "training_mean": metrics(actual, predicted_mean) if actual else None,
        "fold_mae_summary": {
            method: {
                "min": min(f[method]["mae_ppm"] for f in folds),
                "median": statistics.median(f[method]["mae_ppm"] for f in folds),
                "max": max(f[method]["mae_ppm"] for f in folds),
            }
            for method in ("idw", "training_mean")
        } if folds else {},
    }


def make_figure(results) -> None:
    width, height = 820, 390
    left, right, top, bottom = 105, 760, 70, 320
    vals = [r[m]["mae_ppm"] for r in results for m in ("idw", "training_mean")]
    lo, hi = min(vals) - 50, max(vals) + 50
    def y(value):
        return bottom - (value - lo) / (hi - lo) * (bottom - top)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font:13px Arial,sans-serif;fill:#202124}.title{font-size:18px;font-weight:bold}.sub{font-size:12px;fill:#555}.axis{stroke:#777;stroke-width:1}.grid{stroke:#ddd;stroke-width:1}.idw{fill:#276FBF}.mean{fill:#E07A2D}</style>',
        '<text x="40" y="30" class="title">Spatial CV baseline sensitivity to XY grid origin</text>',
        '<text x="40" y="50" class="sub">300-unit tiles; 50-unit buffer; source coordinate units and CRS are undocumented</text>',
    ]
    for i in range(5):
        value = lo + i * (hi - lo) / 4
        yy = y(value)
        parts.append(f'<line x1="{left}" y1="{yy:.1f}" x2="{right}" y2="{yy:.1f}" class="grid"/>')
        parts.append(f'<text x="{left-10}" y="{yy+4:.1f}" text-anchor="end">{value:.0f}</text>')
    parts.extend([
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" class="axis"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" class="axis"/>',
        f'<text transform="translate(24 {(top+bottom)//2}) rotate(-90)" text-anchor="middle">Pooled MAE (ppm)</text>',
    ])
    for i, result in enumerate(results):
        x = left + (i + 0.5) * (right - left) / len(results)
        for method, cls, dx in (("training_mean", "mean", -5), ("idw", "idw", 5)):
            yy = y(result[method]["mae_ppm"])
            parts.append(f'<circle cx="{x+dx:.1f}" cy="{yy:.1f}" r="5" class="{cls}"/>')
        ox, oy = result["origin_offset_xy_coordinate_units"]
        parts.append(f'<text x="{x:.1f}" y="{bottom+20}" text-anchor="middle">{ox:.0f},{oy:.0f}</text>')
    parts.extend([
        f'<text x="{(left+right)//2}" y="{height-18}" text-anchor="middle">Grid-origin shift (X,Y; source coordinate units)</text>',
        f'<circle cx="{right-195}" cy="35" r="5" class="mean"/><text x="{right-184}" y="39">Training mean</text>',
        f'<circle cx="{right-85}" cy="35" r="5" class="idw"/><text x="{right-74}" y="39">3D IDW</text>',
        '</svg>',
    ])
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    FIGURE.write_text("".join(parts), encoding="utf-8")


def main() -> None:
    with DATA.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"HOLEID", "X", "Y", "Z", "Cu ppm"}
    missing = required - set(rows[0]) if rows else required
    if missing:
        raise SystemExit(f"Missing required fields: {sorted(missing)}")
    samples = [
        (r["HOLEID"], tuple(float(r[c]) for c in ("X", "Y", "Z")), float(r["Cu ppm"]))
        for r in rows
    ]
    results = [evaluate(samples, ox, oy) for ox in OFFSETS for oy in OFFSETS]
    report = {
        "target": "Cu ppm",
        "purpose": "Assess sensitivity of fixed spatial-CV baselines to the arbitrary XY grid origin.",
        "design": "Nine buffered XY tile partitions with origin offsets at 0, 1/3, and 2/3 of the 300-coordinate-unit tile width in X and Y.",
        "methods": {"idw": {"power": POWER, "neighbors": K, "distance": "3D Euclidean"},
                    "baseline": "Training-sample mean, recomputed within each fold"},
        "coordinates": "Supplied local X/Y/Z values; CRS, datum, units, and vertical datum are undocumented.",
        "results": results,
        "figure": "outputs/spatial_origin_sensitivity.svg",
        "caveats": [
            "Origin variants are sensitivity designs, not independent validation replicates; folds overlap in their training data.",
            "The coordinate-unit offsets are relative design settings, not real-world mine block dimensions.",
            "These sample-level errors do not validate any unsampled block or resource estimate.",
        ],
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    make_figure(results)
    for result in results:
        ox, oy = result["origin_offset_xy_coordinate_units"]
        print(
            f"origin=({ox:.0f},{oy:.0f}) n={result['evaluated_samples']} "
            f"folds={result['evaluated_folds']} "
            f"IDW MAE={result['idw']['mae_ppm']:.1f} "
            f"mean MAE={result['training_mean']['mae_ppm']:.1f}"
        )
    print(f"Saved {REPORT} and {FIGURE}")


if __name__ == "__main__":
    main()
