"""Summarize nested outer residuals and draw XY/XZ residual maps."""
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports" / "nested_kriging_sensitivity.json"


def metric(rows):
    e = [r["residual_ppm"] for r in rows]
    return {"n": len(e), "mae_ppm": statistics.fmean(abs(x) for x in e), "rmse_ppm": math.sqrt(statistics.fmean(x*x for x in e)), "mean_error_ppm": statistics.fmean(e)}


def main():
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    rows = report["outer_predictions"]
    absr = sorted(abs(r["residual_ppm"]) for r in rows)
    cap = absr[math.ceil(.95*len(absr))-1]
    actual_sorted = sorted(r["observed_cu_ppm"] for r in rows)
    edges = [actual_sorted[math.ceil(q*len(rows))-1] for q in (.25,.5,.75)]
    bands = []
    for i, (lo, hi) in enumerate(zip([float("-inf")]+edges, edges+[float("inf")])):
        group = [r for r in rows if lo < r["observed_cu_ppm"] <= hi]
        bands.append({"band": i+1, "observed_cu_ppm_min": min(r["observed_cu_ppm"] for r in group), "observed_cu_ppm_max": max(r["observed_cu_ppm"] for r in group), **metric(group)})
    summary = {"all": metric(rows), "residual_cap_95_abs_ppm": cap, "by_observed_cu_quartile": bands, "largest_absolute_residuals": sorted(rows, key=lambda r: abs(r["residual_ppm"]), reverse=True)[:20], "selected_hyperparameters": {"k": {}, "family": {}, "range": {}}}
    for f in report["folds"]:
        for key, val in (("k", f["selected_k"]), ("family", f["selected_family"]), ("range", f["selected_range"])):
            summary["selected_hyperparameters"][key][str(val)] = summary["selected_hyperparameters"][key].get(str(val), 0)+1
    (ROOT / "reports" / "residual_diagnostics.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    w,h=960,690
    panels=[("Plan view (X, Y)","x","y",40), ("Vertical section (X, Z)","x","z",370)]
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="white"/><style>text{{font:14px Arial,sans-serif}} .title{{font-size:19px;font-weight:bold}} .sub{{fill:#555;font-size:12px}}</style><text x="40" y="26" class="title">Nested spatial CV: Cu prediction residuals (prediction − observed)</text>']
    for title,cx,cy,oy in panels:
        xs=[r[cx] for r in rows]; ys=[r[cy] for r in rows]
        xlo,xhi=min(xs),max(xs); ylo,yhi=min(ys),max(ys)
        px0,px1=70,900; py0,py1=oy+35,oy+250
        parts.append(f'<text x="40" y="{oy+18}" class="title">{title}</text><rect x="{px0}" y="{py0}" width="{px1-px0}" height="{py1-py0}" fill="#fafafa" stroke="#777"/>')
        for r in rows:
            x=px0+(r[cx]-xlo)/(xhi-xlo)*(px1-px0)
            y=py1-(r[cy]-ylo)/(yhi-ylo)*(py1-py0)
            ratio=min(1,abs(r["residual_ppm"])/cap)
            if r["residual_ppm"] < 0:
                color=f'rgb({round(245-150*ratio)},{round(245-190*ratio)},{round(245-20*ratio)})'
            else:
                color=f'rgb({round(245-15*ratio)},{round(245-150*ratio)},{round(245-190*ratio)})'
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.3" fill="{color}" fill-opacity=".8"/>')
        parts.append(f'<text x="{px0}" y="{py1+20}" class="sub">{cx.upper()} {xlo:.0f} to {xhi:.0f}; {cy.upper()} {ylo:.0f} to {yhi:.0f}</text>')
    parts.append(f'<text x="40" y="675" class="sub">Blue = underprediction; red = overprediction. Color scale clipped at 95th percentile |residual| = {cap:,.0f} ppm. Not a resource estimate.</text></svg>')
    (ROOT / "reports" / "nested_residuals.svg").write_text("".join(parts), encoding="utf-8")
    print(json.dumps({"all": summary["all"], "quartiles": bands, "selected": summary["selected_hyperparameters"], "residual_cap": cap, "outputs": ["residual_diagnostics.json", "nested_residuals.svg"]}, indent=2))


if __name__ == "__main__":
    main()
