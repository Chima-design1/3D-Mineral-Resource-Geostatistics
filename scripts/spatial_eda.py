"""Exploratory spatial maps, directional Cu variograms, and XY-blocked CV."""
from __future__ import annotations

import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "drillholes.csv"
REPORT = ROOT / "reports" / "spatial_eda.json"
FIGURE = ROOT / "reports" / "spatial_eda.svg"
LAG_WIDTH = 20.0
MAX_LAG = 500.0
AZ_TOL = math.radians(22.5)
VERT_TOL = math.sin(AZ_TOL)
BLOCK_SIZE = 300.0
SPATIAL_BUFFER = 50.0
IDW_K = 12
IDW_POWER = 2.0


def read_samples():
    with DATA.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    return [(r["HOLEID"], tuple(float(r[c]) for c in ("X", "Y", "Z")), float(r["Cu ppm"])) for r in rows]


def directional_variograms(samples, value_transform="raw", estimator="matheron"):
    directions = [0, 45, 90, 135]
    bins = int(MAX_LAG / LAG_WIDTH)
    sums = {f"horizontal_{a}_deg": [0.0] * bins for a in directions}
    counts = {name: [0] * bins for name in sums}
    sums["vertical"] = [0.0] * bins
    counts["vertical"] = [0] * bins
    cos_tol = math.cos(AZ_TOL)
    coords = [s[1] for s in samples]
    if value_transform == "log1p":
        values = [math.log1p(s[2]) for s in samples]
    else:
        values = [s[2] for s in samples]
    for i in range(len(samples) - 1):
        x1, y1, z1 = coords[i]
        for j in range(i + 1, len(samples)):
            x2, y2, z2 = coords[j]
            dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
            h = math.hypot(dx, dy)
            d = math.hypot(h, dz)
            if d <= 0 or d > MAX_LAG:
                continue
            b = min(int(d / LAG_WIDTH), bins - 1)
            delta = abs(values[j] - values[i])
            gamma = 0.5 * delta**2 if estimator == "matheron" else math.sqrt(delta)
            if h / d <= VERT_TOL:
                sums["vertical"][b] += gamma
                counts["vertical"][b] += 1
            if h == 0 or abs(dz) / d > VERT_TOL:
                continue
            az = math.atan2(dy, dx) % math.pi
            for a in directions:
                center = math.radians(a) % math.pi
                diff = abs((az - center + math.pi / 2) % math.pi - math.pi / 2)
                if diff <= AZ_TOL:
                    name = f"horizontal_{a}_deg"
                    sums[name][b] += gamma
                    counts[name][b] += 1
    result = {}
    for name in sums:
        result[name] = [
            {"lag_center": (b + 0.5) * LAG_WIDTH, "pairs": n,
             "semivariance": (sums[name][b] / n if estimator == "matheron" else (sums[name][b] / n) ** 4 / (0.457 + 0.494 / n + 0.045 / (n * n))) if n else None}
            for b, n in enumerate(counts[name])
        ]
    return result


def idw_predict(point, train):
    nearest = sorted((math.dist(point, xyz), val) for xyz, val in train)[:IDW_K]
    exact = [v for d, v in nearest if d == 0]
    if exact:
        return statistics.fmean(exact)
    weights = [1.0 / (d ** IDW_POWER) for d, _ in nearest]
    return sum(w * pair[1] for w, pair in zip(weights, nearest)) / sum(weights)


def metrics(actual, predicted):
    errors = [p - y for y, p in zip(actual, predicted)]
    return {"n": len(actual), "mae_ppm": statistics.fmean(map(abs, errors)), "rmse_ppm": math.sqrt(statistics.fmean(e * e for e in errors)), "mean_error_ppm": statistics.fmean(errors)}


def blocked_cv(samples, buffer=0.0):
    xs, ys = [s[1][0] for s in samples], [s[1][1] for s in samples]
    x0, y0 = min(xs), min(ys)
    blocks = defaultdict(list)
    for i, (_, (x, y, _), _) in enumerate(samples):
        blocks[(math.floor((x - x0) / BLOCK_SIZE), math.floor((y - y0) / BLOCK_SIZE))].append(i)
    actual, pred_idw, pred_mean = [], [], []
    for (bx, by), test_indices in blocks.items():
        test_set = set(test_indices)
        left, bottom = x0 + bx * BLOCK_SIZE, y0 + by * BLOCK_SIZE
        right, top = left + BLOCK_SIZE, bottom + BLOCK_SIZE
        train = []
        for i, s in enumerate(samples):
            if i in test_set:
                continue
            x, y, _ = s[1]
            dx = max(left - x, 0.0, x - right)
            dy = max(bottom - y, 0.0, y - top)
            if math.hypot(dx, dy) > buffer:
                train.append((s[1], s[2]))
        if not train:
            continue
        mean = statistics.fmean(v for _, v in train)
        for i in test_indices:
            point, value = samples[i][1], samples[i][2]
            actual.append(value)
            pred_mean.append(mean)
            pred_idw.append(idw_predict(point, train))
    return {"block_size_xy": BLOCK_SIZE, "folds": len(blocks), "buffer_xy": buffer, "split": "Leave one occupied 300-unit XY grid tile out; training points within buffer distance of held-out tile boundaries excluded.", "idw_parameters": {"power": IDW_POWER, "neighbors": IDW_K}, "training_mean": metrics(actual, pred_mean), "idw": metrics(actual, pred_idw)}


def make_svg(samples):
    w, h, pad = 1100, 500, 70
    panels = [("XY plan view", 0, 1), ("X-Z section", 0, 2)]
    coords = [s[1] for s in samples]
    allvals = [s[2] for s in samples]
    logvals = [math.log1p(v) for v in allvals]
    vmin, vmax = min(logvals), max(logvals)
    body = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">', '<rect width="100%" height="100%" fill="#fbfcfe"/>', '<style>text{font-family:Arial,sans-serif;fill:#263445}.title{font-size:22px;font-weight:600}.axis{stroke:#64748b;stroke-width:1.3}.dot{stroke:#fff;stroke-width:.5}</style>']
    for panel, (title, ix, iy) in enumerate(panels):
        left = pad + panel * 510
        top = 95
        pw, ph = 430, 320
        lo_x, hi_x = min(c[ix] for c in coords), max(c[ix] for c in coords)
        lo_y, hi_y = min(c[iy] for c in coords), max(c[iy] for c in coords)
        body.append(f'<text class="title" x="{left}" y="48">{title} • Cu ppm (log color)</text>')
        body.append(f'<line class="axis" x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top+ph}"/><line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+ph}"/>')
        body.append(f'<text x="{left+pw/2}" y="{top+ph+36}" text-anchor="middle">{"X" if ix == 0 else "Y"} (local coordinate)</text>')
        body.append(f'<text x="{left-40}" y="{top+ph/2}" text-anchor="middle" transform="rotate(-90 {left-40} {top+ph/2})">{"Y" if iy == 1 else "Z"} (local coordinate)</text>')
        for (_, c, value), lv in zip(samples, logvals):
            px = left + (c[ix] - lo_x) / (hi_x - lo_x) * pw
            py = top + ph - (c[iy] - lo_y) / (hi_y - lo_y) * ph
            t = (lv - vmin) / (vmax - vmin) if vmax > vmin else 0
            red, green, blue = int(40 + 215*t), int(105 - 65*t), int(180 - 145*t)
            body.append(f'<circle class="dot" cx="{px:.2f}" cy="{py:.2f}" r="3.2" fill="rgb({red},{green},{blue})"><title>Cu: {value:g} ppm</title></circle>')
        body.append(f'<text x="{left}" y="{top+ph+58}" font-size="12">Range: {lo_x:.1f} to {hi_x:.1f}; {lo_y:.1f} to {hi_y:.1f}</text>')
    body.append('<text x="550" y="486" text-anchor="middle" font-size="12">2,000 sample centroids; color is log1p(Cu ppm). Coordinates are local Cartesian; CRS not provided.</text></svg>')
    FIGURE.write_text("\n".join(body), encoding="utf-8")


def main():
    samples = read_samples()
    report = {
        "n_samples": len(samples),
        "target": "Cu ppm",
        "coordinate_context": "Source documentation describes local Cartesian sample-centroid coordinates; CRS/datum not provided.",
        "variogram_parameters": {"max_lag": MAX_LAG, "lag_width": LAG_WIDTH, "horizontal_azimuth_tolerance_degrees": 22.5, "vertical_tolerance_degrees_from_vertical": 22.5, "estimators": ["Matheron raw Cu", "Matheron log1p(Cu)", "Cressie-Hawkins robust raw Cu"]},
        "directional_variograms_raw": directional_variograms(samples),
        "directional_variograms_log1p": directional_variograms(samples, value_transform="log1p"),
        "directional_variograms_robust_raw": directional_variograms(samples, estimator="cressie_hawkins"),
        "xy_blocked_validation": blocked_cv(samples),
        "xy_blocked_validation_buffered": blocked_cv(samples, buffer=SPATIAL_BUFFER),
        "figure": FIGURE.name,
        "limitations": ["Full-data experimental variograms are exploratory only; directional variograms must be refit within each training fold for unbiased kriging validation.", "Raw Cu is strongly right-skewed; sensitivity variants are exploratory and do not themselves establish a valid transformation or variogram model.", "The buffered XY folds may discard substantial training data in some folds.", "No geological domains, collars/surveys, density, or resource-estimation support are available."],
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    make_svg(samples)
    print(json.dumps({"n_samples": len(samples), "xy_blocked_validation": report["xy_blocked_validation"], "xy_blocked_validation_buffered": report["xy_blocked_validation_buffered"], "variogram_pair_counts_at_100_units": {name: next((r["pairs"] for r in series if r["lag_center"] == 110.0), 0) for name, series in report["directional_variograms_raw"].items()}, "reports": [REPORT.name, FIGURE.name]}, indent=2))


if __name__ == "__main__":
    main()
