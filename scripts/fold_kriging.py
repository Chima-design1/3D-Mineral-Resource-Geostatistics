"""Experimental isotropic ordinary-kriging CV with fold-only variogram fitting.

Models: ordinary kriging on raw Cu; ordinary kriging on log1p(Cu) with
lognormal variance correction. Variogram family/range/sill/nugget are selected
from each training fold's sampled empirical variogram only.
"""
from __future__ import annotations

import csv
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "drillholes.csv"
OUT = ROOT / "reports" / "fold_kriging.json"
BLOCK = 300.0
BUFFER = 50.0
LAG_WIDTH = 20.0
MAX_LAG = 500.0
PAIR_DRAWS = 40000
K = 12
RANGES = (40, 60, 80, 100, 125, 150, 200, 250, 300, 400, 500)
FAMILIES = ("spherical", "exponential", "gaussian")


def load():
    with DATA.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    samples = []
    for r in rows:
        samples.append({"hole": r["HOLEID"], "xyz": tuple(float(r[c]) for c in ("X", "Y", "Z")), "cu": float(r["Cu ppm"])})
    return samples


def grid_folds(samples):
    x0 = min(s["xyz"][0] for s in samples)
    y0 = min(s["xyz"][1] for s in samples)
    folds = defaultdict(list)
    for i, s in enumerate(samples):
        x, y, _ = s["xyz"]
        folds[(math.floor((x - x0) / BLOCK), math.floor((y - y0) / BLOCK))].append(i)
    output = []
    for (bx, by), tests in sorted(folds.items()):
        left, bottom = x0 + bx * BLOCK, y0 + by * BLOCK
        right, top = left + BLOCK, bottom + BLOCK
        train = []
        for i, s in enumerate(samples):
            if i in tests:
                continue
            x, y, _ = s["xyz"]
            dx, dy = max(left - x, 0.0, x - right), max(bottom - y, 0.0, y - top)
            if math.hypot(dx, dy) > BUFFER:
                train.append(i)
        output.append((tests, train))
    return output


def semivar(d, family, nugget, psill, a):
    if d <= 1e-12:
        return 0.0
    if family == "pure_nugget":
        return nugget
    r = d / a
    if family == "spherical":
        shape = 1.5 * r - 0.5 * r**3 if r < 1 else 1.0
    elif family == "exponential":
        shape = 1.0 - math.exp(-3.0 * r)
    else:
        shape = 1.0 - math.exp(-3.0 * r * r)
    return nugget + psill * shape


def empirical(train, samples, transformed, seed):
    rng = random.Random(seed)
    n = len(train)
    sums, counts = [0.0] * int(MAX_LAG / LAG_WIDTH), [0] * int(MAX_LAG / LAG_WIDTH)
    vals = [math.log1p(samples[i]["cu"]) if transformed else samples[i]["cu"] for i in train]
    pts = [samples[i]["xyz"] for i in train]
    for _ in range(PAIR_DRAWS):
        i, j = rng.randrange(n), rng.randrange(n - 1)
        if j >= i:
            j += 1
        d = math.dist(pts[i], pts[j])
        if 0 < d <= MAX_LAG:
            b = min(int(d / LAG_WIDTH), len(counts) - 1)
            sums[b] += 0.5 * (vals[i] - vals[j]) ** 2
            counts[b] += 1
    bins = []
    for b, count in enumerate(counts):
        if count:
            bins.append(((b + 0.5) * LAG_WIDTH, sums[b] / count, count))
    return bins


def fit_variogram(bins, allowed_families=None, allowed_ranges=None):
    if len(bins) < 4:
        raise ValueError("Too few populated empirical variogram bins")
    best = None
    for family in (allowed_families or FAMILIES):
        for a in (allowed_ranges or RANGES):
            shapes = []
            for h, _, _ in bins:
                shapes.append(semivar(h, family, 0, 1, a))
            # Weighted regression: gamma = nugget + partial_sill * shape.
            ws = [math.sqrt(n) for _, _, n in bins]
            sw = sum(w*w for w in ws)
            sx = sum(w*w*x for w, x in zip(ws, shapes))
            sy = sum(w*w*y for w, (_, y, _) in zip(ws, bins))
            sxx = sum(w*w*x*x for w, x in zip(ws, shapes))
            sxy = sum(w*w*x*y for w, x, (_, y, _) in zip(ws, shapes, bins))
            den = sw * sxx - sx * sx
            if den <= 1e-12:
                continue
            psill = (sw * sxy - sx * sy) / den
            nugget = (sy - psill * sx) / sw
            if psill <= 0:
                continue
            if nugget < 0:
                nugget = 0.0
                psill = sxy / max(sxx, 1e-12)
            error = sum(n * (y - (nugget + psill * shape)) ** 2 for (h, y, n), shape in zip(bins, shapes))
            if best is None or error < best["weighted_sse"]:
                best = {"family": family, "range": a, "nugget": nugget, "partial_sill": psill, "sill": nugget + psill, "weighted_sse": error, "fit_bins": len(bins)}
    if best is None:
        # A pure-nugget fallback is valid when the constrained structured fit
        # has no positive partial sill; it intentionally carries no continuity.
        sw = sum(n for _, _, n in bins)
        nugget = sum(y * n for _, y, n in bins) / sw
        best = {"family": "pure_nugget", "range": 0.0, "nugget": max(nugget, 1e-9), "partial_sill": 0.0, "sill": max(nugget, 1e-9), "weighted_sse": sum(n * (y-nugget)**2 for _, y, n in bins), "fit_bins": len(bins), "fallback": True}
    return best


def solve(matrix, rhs):
    n = len(rhs)
    a = [matrix[i][:] + [rhs[i]] for i in range(n)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-12:
            raise ArithmeticError("Singular ordinary-kriging system")
        a[col], a[pivot] = a[pivot], a[col]
        div = a[col][col]
        for j in range(col, n + 1):
            a[col][j] /= div
        for r in range(n):
            if r == col:
                continue
            factor = a[r][col]
            if factor:
                for j in range(col, n + 1):
                    a[r][j] -= factor * a[col][j]
    return [a[i][n] for i in range(n)]


def ordinary_kriging(point, train, samples, fit, transformed, k=None):
    nearest = sorted((math.dist(point, samples[i]["xyz"]), i) for i in train)[:(k or K)]
    coords = [samples[i]["xyz"] for _, i in nearest]
    vals = [math.log1p(samples[i]["cu"]) if transformed else samples[i]["cu"] for _, i in nearest]
    n = len(nearest)
    matrix = [[0.0] * (n + 1) for _ in range(n + 1)]
    gamma0 = []
    for i in range(n):
        for j in range(n):
            matrix[i][j] = semivar(math.dist(coords[i], coords[j]), fit["family"], fit["nugget"], fit["partial_sill"], fit["range"])
        matrix[i][n] = matrix[n][i] = 1.0
        gamma0.append(semivar(nearest[i][0], fit["family"], fit["nugget"], fit["partial_sill"], fit["range"]))
    solution = solve(matrix, gamma0 + [1.0])
    weights, lagrange = solution[:n], solution[n]
    estimate = sum(w * v for w, v in zip(weights, vals))
    variance = max(0.0, sum(w * g for w, g in zip(weights, gamma0)) + lagrange)
    if transformed:
        # Conditional lognormal mean correction from ordinary-kriging variance.
        estimate = math.expm1(min(estimate + 0.5 * variance, 700.0))
    return max(0.0, estimate)


def metrics(actual, predicted):
    error = [p-y for y, p in zip(actual, predicted)]
    return {"n": len(actual), "mae_ppm": statistics.fmean(abs(e) for e in error), "rmse_ppm": math.sqrt(statistics.fmean(e*e for e in error)), "mean_error_ppm": statistics.fmean(error)}


def main():
    samples = load()
    folds = grid_folds(samples)
    predictions = {"ok_raw_cu": [], "ok_log1p_lognormal_mean": []}
    actual = []
    fold_details = []
    for fold_id, (test_ids, train_ids) in enumerate(folds):
        fits = {}
        for transformed, name in ((False, "raw"), (True, "log1p")):
            bins = empirical(train_ids, samples, transformed, seed=1717 + fold_id)
            fits[name] = fit_variogram(bins)
        for i in test_ids:
            point = samples[i]["xyz"]
            actual.append(samples[i]["cu"])
            for transformed, name, out_name in ((False, "raw", "ok_raw_cu"), (True, "log1p", "ok_log1p_lognormal_mean")):
                predictions[out_name].append(ordinary_kriging(point, train_ids, samples, fits[name], transformed))
        fold_details.append({"fold": fold_id + 1, "test_n": len(test_ids), "train_n": len(train_ids), "raw_fit": fits["raw"], "log1p_fit": fits["log1p"]})
        print(f"completed spatial fold {fold_id + 1}/{len(folds)}", flush=True)
    report = {
        "target": "Cu ppm",
        "validation": {"type": "XY tile CV with training-only buffer", "tile_xy": BLOCK, "buffer_xy": BUFFER, "folds": len(folds)},
        "variogram_fitting": {"pairs_randomly_drawn_per_training_fold": PAIR_DRAWS, "max_lag": MAX_LAG, "lag_width": LAG_WIDTH, "families_compared": list(FAMILIES), "selection": "minimum pair-count-weighted SSE across training-fold experimental bins", "neighbor_count": K},
        "methods": {name: metrics(actual, preds) for name, preds in predictions.items()},
        "folds": fold_details,
        "caveats": ["All model and variogram fitting uses training samples only; holdout grades are used only for final scores.", "Local isotropic variograms are selected from random-pair empirical bins and may be unstable in sparse folds.", "The log1p model uses kriging-variance lognormal back-transformation; this is an exploratory assumption.", "No nested spatial validation was used to tune the design choices, and no geological domain model is available."],
    }
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"methods": report["methods"], "variogram_selection_counts": {name: {fam: sum(f["raw_fit" if name == "raw" else "log1p_fit"]["family"] == fam for f in fold_details) for fam in FAMILIES} for name in ("raw", "log1p")}, "report": OUT.name}, indent=2))


if __name__ == "__main__":
    main()
