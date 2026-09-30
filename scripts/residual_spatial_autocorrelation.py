"""Descriptive Moran's I diagnostics on buffered-CV out-of-fold residuals."""
import csv
import json
import math
import random
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
N_PERMUTATIONS = 499
K_VALUES = (4, 8, 12)
SEED = 20260930


def morans_i(values, edges):
    n = len(values)
    mean = statistics.fmean(values)
    centered = [v - mean for v in values]
    denom = math.fsum(v * v for v in centered)
    if denom == 0:
        return 0.0
    numerator = math.fsum(centered[i] * centered[j] for i, j in edges)
    return (n / len(edges)) * (numerator / denom)


def main():
    with (ROOT / "data/raw/drillholes.csv").open(encoding="utf-8-sig", newline="") as f:
        source = list(csv.DictReader(f))
    xyz = [(float(r["X"]), float(r["Y"]), float(r["Z"])) for r in source]
    n = len(xyz)
    if n < max(K_VALUES) + 1:
        raise ValueError("Too few sample points for the configured nearest-neighbor graph")

    neighbors = []
    for i, point in enumerate(xyz):
        distances = []
        for j, other in enumerate(xyz):
            if i == j:
                continue
            d2 = sum((a - b) ** 2 for a, b in zip(point, other))
            distances.append((d2, j))
        neighbors.append([j for _, j in sorted(distances)[:max(K_VALUES)]])
    edge_sets = {k: [(i, j) for i, js in enumerate(neighbors) for j in js[:k]]
                 for k in K_VALUES}

    kriging = json.loads((ROOT / "reports/nested_kriging_sensitivity.json").read_text(encoding="utf-8"))
    ridge_low = json.loads((ROOT / "reports/geochem_covariate_cv.json").read_text(encoding="utf-8"))
    feature = json.loads((ROOT / "reports/geochem_feature_sensitivity.json").read_text(encoding="utf-8"))
    model_rows = {
        "nested_ordinary_kriging": kriging["outer_predictions"],
        "conditional_ridge_all17_alpha_0.1_to_100": ridge_low["outer_predictions"],
        "conditional_ridge_all17_alpha_100_to_10000": feature["models"]["all_17"]["outer_predictions"],
        "conditional_ridge_S_Ag_Au": feature["models"]["S_Ag_Au"]["outer_predictions"],
    }

    rng = random.Random(SEED)
    results = {}
    for name, rows in model_rows.items():
        by_row = {int(r["row"]): r for r in rows}
        if len(by_row) != n:
            raise ValueError(f"{name}: expected {n} unique out-of-fold rows, found {len(by_row)}")
        residual = [float(by_row[i + 1]["residual_ppm"]) for i in range(n)]
        fold = [int(by_row[i + 1]["outer_fold"]) for i in range(n)]
        fold_groups = {}
        for i, fold_id in enumerate(fold):
            fold_groups.setdefault(fold_id, []).append(i)
        results[name] = {}
        for k, edges in edge_sets.items():
            model_stats = {}
            for series_name, values in (("signed_residual_ppm", residual),
                                        ("absolute_residual_ppm", [abs(x) for x in residual])):
                observed = morans_i(values, edges)
                permuted = []
                for _ in range(N_PERMUTATIONS):
                    shuffled = list(values)
                    for group in fold_groups.values():
                        group_values = [values[i] for i in group]
                        rng.shuffle(group_values)
                        for i, value in zip(group, group_values):
                            shuffled[i] = value
                    permuted.append(morans_i(shuffled, edges))
                exceed = sum(abs(v) >= abs(observed) for v in permuted)
                ordered = sorted(permuted)
                model_stats[series_name] = {
                    "moran_i": observed,
                    "permutation_mean_i": statistics.fmean(permuted),
                    "permutation_95pct_range_i": [
                        ordered[math.floor(0.025 * N_PERMUTATIONS)],
                        ordered[math.floor(0.975 * N_PERMUTATIONS)],
                    ],
                    "two_sided_permutation_p": (exceed + 1) / (N_PERMUTATIONS + 1),
                    "permutations": N_PERMUTATIONS,
                }
            results[name][str(k)] = model_stats

    report = {
        "design": "Global Moran's I over 3D Euclidean k-nearest-neighbor graph built from supplied sample-centroid XYZ coordinates.",
        "n_samples": n,
        "k_neighbors_sensitivity": list(K_VALUES),
        "directed_edges_by_k": {str(k): len(edges) for k, edges in edge_sets.items()},
        "residual_definition": "predicted Cu minus observed Cu; all predictions are out-of-fold under the 300-unit XY tile / 50-unit buffer design.",
        "permutation": {"type": "random-label permutation restricted within outer validation fold; two-sided statistic", "n": N_PERMUTATIONS, "seed": SEED},
        "results": results,
        "caveats": [
            "Coordinates are sample centroids; CRS and coordinate units are undocumented.",
            "Moran's I is a descriptive diagnostic on cross-validated residuals, not a formal inferential test for a target population.",
            "Within-fold permutation preserves each fold's residual distribution but does not make the folds independent or account for their spatial design.",
            "A 3D Euclidean graph assumes comparable scales across X, Y, and Z and cannot substitute for geological domains or anisotropy analysis.",
        ],
    }
    out = ROOT / "reports/residual_spatial_autocorrelation.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"report": out.name, "results": results}, indent=2))


if __name__ == "__main__":
    main()
