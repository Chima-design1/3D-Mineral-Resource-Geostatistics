"""Map paired out-of-fold absolute-error differences for Cu models."""
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def summarize(rows):
    y = [r["observed_cu_ppm"] for r in rows]
    krig_abs = [abs(r["kriging_residual_ppm"]) for r in rows]
    ridge_abs = [abs(r["ridge_residual_ppm"]) for r in rows]
    delta = [a - b for a, b in zip(krig_abs, ridge_abs)]
    return {
        "n": len(rows),
        "kriging_mae_ppm": statistics.fmean(krig_abs),
        "ridge_mae_ppm": statistics.fmean(ridge_abs),
        "paired_mae_gain_ppm_positive_means_ridge_better": statistics.fmean(delta),
        "ridge_lower_absolute_error_fraction": sum(v > 0 for v in delta) / len(delta),
        "ridge_tie_fraction": sum(v == 0 for v in delta) / len(delta),
        "ridge_higher_absolute_error_fraction": sum(v < 0 for v in delta) / len(delta),
        "paired_delta_median_ppm": statistics.median(delta),
    }


def percentile(values, p):
    ordered = sorted(values)
    index = max(0, min(len(ordered) - 1, math.ceil(p * len(ordered)) - 1))
    return ordered[index]


def main():
    krig = json.loads((ROOT / "reports/nested_kriging_sensitivity.json").read_text(encoding="utf-8"))["outer_predictions"]
    chemistry = json.loads((ROOT / "reports/geochem_covariate_cv.json").read_text(encoding="utf-8"))["outer_predictions"]
    high_alpha = json.loads((ROOT / "reports/geochem_feature_sensitivity.json").read_text(encoding="utf-8"))["models"]["all_17"]["outer_predictions"]
    three_assays = json.loads((ROOT / "reports/geochem_feature_sensitivity.json").read_text(encoding="utf-8"))["models"]["S_Ag_Au"]["outer_predictions"]
    k_by_row = {int(r["row"]): r for r in krig}
    low_by_row = {int(r["row"]): r for r in chemistry}
    high_by_row = {int(r["row"]): r for r in high_alpha}
    three_by_row = {int(r["row"]): r for r in three_assays}
    if not (len(k_by_row) == len(low_by_row) == len(high_by_row) == len(three_by_row) == 2000):
        raise ValueError("Expected 2,000 unique out-of-fold rows for each compared model")

    rows = []
    for row_id, k in k_by_row.items():
        low, high, three = low_by_row[row_id], high_by_row[row_id], three_by_row[row_id]
        if len({k["outer_fold"], low["outer_fold"], high["outer_fold"], three["outer_fold"]}) != 1:
            raise ValueError(f"Model fold mismatch at row {row_id}")
        rows.append({
            "row": row_id, "outer_fold": int(k["outer_fold"]),
            "x": float(k["x"]), "y": float(k["y"]),
            "observed_cu_ppm": float(k["observed_cu_ppm"]),
            "kriging_residual_ppm": float(k["residual_ppm"]),
            "ridge_residual_ppm": float(low["residual_ppm"]),
            "ridge_high_alpha_residual_ppm": float(high["residual_ppm"]),
            "ridge_three_assay_residual_ppm": float(three["residual_ppm"]),
        })

    grades = sorted(r["observed_cu_ppm"] for r in rows)
    q75 = percentile(grades, .75)
    q90 = percentile(grades, .90)
    groups = {
        "all": rows,
        "observed_cu_below_q75": [r for r in rows if r["observed_cu_ppm"] < q75],
        "observed_cu_q75_to_q90": [r for r in rows if q75 <= r["observed_cu_ppm"] < q90],
        "observed_cu_at_or_above_q90": [r for r in rows if r["observed_cu_ppm"] >= q90],
    }
    summaries = {}
    for group_name, group in groups.items():
        values = summarize(group)
        for model_name, key in (("all17_alpha_100_to_10000", "ridge_high_alpha_residual_ppm"),
                                ("S_Ag_Au", "ridge_three_assay_residual_ppm")):
            errors = [abs(r[key]) for r in group]
            values[model_name + "_mae_ppm"] = statistics.fmean(errors)
        summaries[group_name] = values
    fold_summary = {str(fold_id): summarize([r for r in rows if r["outer_fold"] == fold_id])
                    for fold_id in sorted({r["outer_fold"] for r in rows})}
    report = {
        "design": "Paired row-level comparison of out-of-fold Cu absolute errors on the same buffered XY holdout rows.",
        "conditional_model": "Ridge regression using the 17 other co-located assays; alpha selected within nested buffered spatial CV.",
        "comparator": "Nested ordinary kriging using only training-fold Cu and sample-centroid geometry.",
        "grade_thresholds_ppm": {"q75": q75, "q90": q90},
        "groups": summaries,
        "by_outer_fold": fold_summary,
        "map_delta_definition": "absolute kriging residual minus absolute alpha 0.1–100 ridge residual; positive green values mean the conditional ridge has lower absolute error.",
        "map_color_cap_ppm": percentile([abs(abs(r["kriging_residual_ppm"]) - abs(r["ridge_residual_ppm"])) for r in rows], .95),
        "caveats": [
            "The conditional ridge requires all 17 non-Cu assays at each target location.",
            "This is sample-level paired prediction comparison, not block-support resource estimation.",
            "XY coordinates are supplied sample-centroid coordinates with undocumented CRS and units; the map is not a geographic map.",
            "Local patterns do not identify geological cause and must not be interpreted as domains or an economic priority map.",
        ],
    }
    (ROOT / "reports/paired_error_comparison.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    w, h = 1060, 760
    x0, x1, y0, y1 = min(r["x"] for r in rows), max(r["x"] for r in rows), min(r["y"] for r in rows), max(r["y"] for r in rows)
    px0, px1, py0, py1 = 90, 980, 100, 645
    cap = report["map_color_cap_ppm"]
    elements = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
                '<rect width="100%" height="100%" fill="white"/>',
                '<style>text{font:15px Arial,sans-serif}.title{font-size:22px;font-weight:bold}.small{font-size:12px;fill:#555}</style>',
                '<text x="40" y="38" class="title">Paired Cu prediction error difference at held-out samples</text>',
                '<text x="40" y="64" class="small">3D XY projection · positive difference means assay-assisted ridge has lower absolute error · color clipped at 95th percentile</text>',
                f'<rect x="{px0}" y="{py0}" width="{px1-px0}" height="{py1-py0}" fill="#fafafa" stroke="#777"/>']
    for r in rows:
        x = px0 + (r["x"] - x0) / (x1 - x0) * (px1 - px0)
        y = py1 - (r["y"] - y0) / (y1 - y0) * (py1 - py0)
        delta = abs(r["kriging_residual_ppm"]) - abs(r["ridge_residual_ppm"])
        t = min(1.0, abs(delta) / cap)
        if delta > 0:
            color = f'rgb({round(245-170*t)},{round(245-105*t)},{round(245-190*t)})'
        elif delta < 0:
            color = f'rgb({round(245-10*t)},{round(245-135*t)},{round(245-205*t)})'
        else:
            color = "#eeeeee"
        elements.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.5" fill="{color}" fill-opacity=".82"/>')
    elements += [f'<text x="{px0}" y="{py1+25}" class="small">X: {x0:.1f} to {x1:.1f} (source coordinate units)</text>',
                 f'<text x="{px0}" y="{py1+48}" class="small">Y: {y0:.1f} to {y1:.1f} (source coordinate units)</text>',
                 '<rect x="90" y="710" width="18" height="14" fill="#4b8c37"/><text x="115" y="722" class="small">Ridge improves</text>',
                 '<rect x="270" y="710" width="18" height="14" fill="#e76e28"/><text x="295" y="722" class="small">Kriging improves</text>',
                 f'<text x="515" y="722" class="small">Color scale clipped at |error difference| = {cap:,.0f} ppm. Conditional model assumes 17 assays are available.</text>',
                 '</svg>']
    (ROOT / "outputs/paired_error_comparison.svg").write_text("".join(elements), encoding="utf-8")
    print(json.dumps({"report": "paired_error_comparison.json", "figure": "paired_error_comparison.svg", "summary": summaries}, indent=2))


if __name__ == "__main__":
    main()
