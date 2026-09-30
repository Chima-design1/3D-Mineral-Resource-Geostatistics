"""Nested XY-block CV for choosing the local ordinary-kriging neighborhood."""
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import fold_kriging as fk

KS = (8, 16, 24, 32)
FAMILIES = ("spherical", "gaussian")
RANGES = (80, 200, 400)


def inner_splits(pool, samples, x0, y0):
    tiles = defaultdict(list)
    for i in pool:
        x, y, _ = samples[i]["xyz"]
        tiles[(math.floor((x-x0)/fk.BLOCK), math.floor((y-y0)/fk.BLOCK))].append(i)
    out = []
    pool_set = set(pool)
    for (bx, by), test in sorted(tiles.items()):
        left, bottom = x0+bx*fk.BLOCK, y0+by*fk.BLOCK
        right, top = left+fk.BLOCK, bottom+fk.BLOCK
        train = []
        for i in pool_set.difference(test):
            x, y, _ = samples[i]["xyz"]
            dx, dy = max(left-x, 0.0, x-right), max(bottom-y, 0.0, y-top)
            if math.hypot(dx, dy) > fk.BUFFER:
                train.append(i)
        if len(test) and len(train) >= max(KS):
            out.append((test, train))
    # Use two geographically distributed inner tiles per outer fold to keep
    # the nested experiment tractable while preserving spatial separation.
    if len(out) > 2:
        positions = [0, len(out)-1]
        out = [out[i] for i in positions]
    return out


def mae(actual, pred):
    return statistics.fmean(abs(p-y) for y, p in zip(actual, pred))


def main():
    # Nested CV repeats variogram fitting many times; use a deterministic, lower
    # pair sample here while retaining the same max lag and binning design.
    fk.PAIR_DRAWS = 10000
    samples = fk.load()
    outer = fk.grid_folds(samples)
    x0 = min(s["xyz"][0] for s in samples)
    y0 = min(s["xyz"][1] for s in samples)
    actual_all, predicted_all = [], []
    result_folds = []
    outer_predictions = []
    for ofold, (test_ids, outer_train) in enumerate(outer):
        splits = inner_splits(outer_train, samples, x0, y0)
        errors = {(family, ran, k): [] for family in FAMILIES for ran in RANGES for k in KS}
        for ifold, (ival, itrain) in enumerate(splits):
            bins = fk.empirical(itrain, samples, False, seed=50000 + ofold*100 + ifold)
            fits = {(family, ran): fk.fit_variogram(bins, [family], [ran]) for family in FAMILIES for ran in RANGES}
            for i in ival:
                point = samples[i]["xyz"]
                y = samples[i]["cu"]
                for (family, ran), fit in fits.items():
                    for k in KS:
                        pred = fk.ordinary_kriging(point, itrain, samples, fit, False, k=k)
                        errors[(family, ran, k)].append(abs(pred-y))
        candidate_mae = {"|".join(map(str, key)): statistics.fmean(v) for key, v in errors.items() if v}
        selected_key = min(errors, key=lambda key: statistics.fmean(errors[key]))
        selected = selected_key[2]
        outer_fit = fk.fit_variogram(fk.empirical(outer_train, samples, False, seed=90000 + ofold), [selected_key[0]], [selected_key[1]])
        ys, ps = [], []
        for i in test_ids:
            ys.append(samples[i]["cu"])
            ps.append(fk.ordinary_kriging(samples[i]["xyz"], outer_train, samples, outer_fit, False, k=selected))
            outer_predictions.append({"outer_fold": ofold+1, "row": i+1, "x": samples[i]["xyz"][0], "y": samples[i]["xyz"][1], "z": samples[i]["xyz"][2], "observed_cu_ppm": samples[i]["cu"], "predicted_cu_ppm": ps[-1], "residual_ppm": ps[-1]-ys[-1]})
        actual_all.extend(ys)
        predicted_all.extend(ps)
        result_folds.append({"outer_fold": ofold+1, "n_inner_folds": len(splits), "selected_family": selected_key[0], "selected_range": selected_key[1], "selected_k": selected, "inner_best_candidates_mae_ppm": sorted(({"family": key.split("|")[0], "range": int(key.split("|")[1]), "k": int(key.split("|")[2]), "mae_ppm": val} for key,val in candidate_mae.items()), key=lambda r:r["mae_ppm"])[:8], "outer_test_n": len(ys), "outer_test_mae_ppm": mae(ys, ps), "outer_test_rmse_ppm": math.sqrt(statistics.fmean((p-y)**2 for y,p in zip(ys,ps))), "outer_variogram": outer_fit})
        print(f"nested outer fold {ofold+1}/{len(outer)}; selected {selected_key}", flush=True)
    report = {"design": "Outer 300-unit XY tile CV with 50-unit horizontal buffer; within each outer training set, two geographically distributed occupied tiles are used as buffered inner folds. Inner tuning jointly selects neighborhood count and constrained variogram family/range.", "target": "Cu ppm", "candidate_k": list(KS), "candidate_families": list(FAMILIES), "candidate_ranges": list(RANGES), "pair_draws_per_fit": fk.PAIR_DRAWS, "pooled_outer_mae_ppm": mae(actual_all, predicted_all), "pooled_outer_rmse_ppm": math.sqrt(statistics.fmean((p-y)**2 for y,p in zip(actual_all,predicted_all))), "folds": result_folds, "outer_predictions": outer_predictions, "caveats": ["Inner tuning uses two distributed inner tiles rather than every available tile; results can be sensitive to this small inner sample.", "The grid varies only family, range, and neighborhood size; nugget/partial sill are fitted on each training variogram.", "Occupied XY tiles share training data and are descriptive spatial partitions, not independent replications.", "The grid uses supplied local coordinates and no geologic domains; results are exploratory, not resource estimates."]}
    out = ROOT / "reports" / "nested_kriging_sensitivity.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"pooled_outer_mae_ppm": report["pooled_outer_mae_ppm"], "pooled_outer_rmse_ppm": report["pooled_outer_rmse_ppm"], "selected_k_counts": {str(k): sum(f["selected_k"] == k for f in result_folds) for k in KS}, "report": out.name}, indent=2))


if __name__ == "__main__":
    main()
