"""Summarize outer spatial-fold errors and training-coverage diagnostics."""
import json
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import fold_kriging as fk


def scores(actual, pred):
    err = [p - y for y, p in zip(actual, pred)]
    return {"n": len(err), "mae_ppm": statistics.fmean(abs(e) for e in err), "rmse_ppm": math.sqrt(statistics.fmean(e*e for e in err)), "mean_error_ppm": statistics.fmean(err)}


def main():
    samples = fk.load()
    folds = fk.grid_folds(samples)
    details, all_actual = [], []
    predictions = {"training_mean": [], "idw_p2_k12": [], "ok_raw_k12": [], "ok_log1p_k12": []}
    for fno, (test_ids, train_ids) in enumerate(folds):
        fit_raw = fk.fit_variogram(fk.empirical(train_ids, samples, False, 1717 + fno))
        fit_log = fk.fit_variogram(fk.empirical(train_ids, samples, True, 1717 + fno))
        train_mean = statistics.fmean(samples[i]["cu"] for i in train_ids)
        ys, preds = [], {k: [] for k in predictions}
        nearest = []
        for i in test_ids:
            s = samples[i]
            point = s["xyz"]
            ys.append(s["cu"])
            ds = sorted((math.dist(point, samples[j]["xyz"]), j) for j in train_ids)
            nearest.append(ds[0][0])
            preds["training_mean"].append(train_mean)
            chosen = ds[:fk.K]
            w = [1 / max(d, 1e-9) ** 2 for d, _ in chosen]
            preds["idw_p2_k12"].append(sum(ww * samples[j]["cu"] for ww, (_, j) in zip(w, chosen)) / sum(w))
            preds["ok_raw_k12"].append(fk.ordinary_kriging(point, train_ids, samples, fit_raw, False))
            preds["ok_log1p_k12"].append(fk.ordinary_kriging(point, train_ids, samples, fit_log, True))
        all_actual.extend(ys)
        for key in predictions:
            predictions[key].extend(preds[key])
        detail = {"fold": fno + 1, "n_test": len(ys), "n_train": len(train_ids), "nearest_train_distance_median": statistics.median(nearest), "nearest_train_distance_max": max(nearest), "methods": {k: scores(ys, v) for k, v in preds.items()}, "raw_variogram": fit_raw, "log1p_variogram": fit_log}
        details.append(detail)
        print(f"diagnosed fold {fno + 1}/{len(folds)}", flush=True)
    report = {"validation": "same 22 occupied 300-unit XY tiles, 50-unit horizontal buffer", "overall": {k: scores(all_actual, v) for k, v in predictions.items()}, "folds": details, "interpretation": "Nearest training distance is in the supplied local XYZ coordinate units. Fold scores are descriptive; folds share training observations and are not independent replicates."}
    out = ROOT / "reports" / "fold_diagnostics.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    worst = sorted(details, key=lambda d: d["methods"]["ok_raw_k12"]["mae_ppm"], reverse=True)[:3]
    print(json.dumps({"overall": report["overall"], "worst_raw_kriging_folds": [{"fold": f["fold"], "mae": f["methods"]["ok_raw_k12"]["mae_ppm"], "nearest_distance_median": f["nearest_train_distance_median"]} for f in worst], "report": out.name}, indent=2))


if __name__ == "__main__":
    main()
