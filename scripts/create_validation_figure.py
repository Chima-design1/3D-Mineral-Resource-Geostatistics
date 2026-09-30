"""Render observed-versus-predicted Cu diagnostics for buffered spatial holdouts."""
import json
import math
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REPORTS=ROOT/"reports"
FIG=ROOT/"outputs"/"model_validation_figure.svg"
META=ROOT/"reports"/"model_validation_figure.json"
Q75,Q90=8900.0,18100.0


def mapping(filename,key):
    d=json.loads((REPORTS/filename).read_text(encoding="utf-8"))
    return {int(r["row"]):float(r[key]) for r in d["outer_predictions"]}


def main():
    grid=json.loads((REPORTS/"grid_support_validation.json").read_text(encoding="utf-8"))["records"]
    actual={int(r["row"]):float(r["actual"]) for r in grid}
    models={"Raw Cu IDW, p=2, k=12":{int(r["row"]):float(r["idw"]) for r in grid},"Nested ordinary kriging":mapping("nested_kriging_sensitivity.json","predicted_cu_ppm"),"Conditional all-assay ridge":mapping("geochem_covariate_cv.json","predicted_cu_ppm")}
    feature=json.loads((REPORTS/"geochem_feature_sensitivity.json").read_text(encoding="utf-8"))["models"]["all_17"]["outer_predictions"]
    models["Conditional ridge, alpha 100–10,000"]={int(r["row"]):float(r["predicted_cu_ppm"]) for r in feature}
    models["Conditional ridge, alpha 0.1–100"]=models.pop("Conditional all-assay ridge")
    sub=json.loads((REPORTS/"subcomposition_spatial_sensitivity.json").read_text(encoding="utf-8"))["designs"]["buffered_XY_tile"]["predictions_by_row"]
    models["Cu-Au-Ag-S subcomposition IDW"]={int(r["row"]):float(r["clr_idw_total_idw_zero_0.1x_min"]) for r in sub}
    rows=sorted(actual); vmax=max(actual.values());fmax=math.log1p(vmax);ticks=(0,5000,10000,20000,40000,60000)
    nrows=math.ceil(len(models)/2);W,H=1200,360+350*nrows;pw,ph=470,300;body=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">','<rect width="100%" height="100%" fill="#fff"/>','<style>text{font-family:Arial,sans-serif;fill:#243447}.main{font-size:25px;font-weight:700}.sub{font-size:14px;fill:#475569}.title{font-size:16px;font-weight:600}.axis{stroke:#526174;stroke-width:1}.grid{stroke:#e2e8f0;stroke-width:1}.diag{stroke:#334155;stroke-width:1.5;stroke-dasharray:7 5}.thresh{stroke:#c2410c;stroke-width:1;stroke-dasharray:4 4}.dot{stroke:#fff;stroke-width:.2}</style>']
    body.extend(['<text class="main" x="600" y="38" text-anchor="middle">Buffered spatial holdouts: observed vs predicted Cu</text>','<text class="sub" x="600" y="62" text-anchor="middle">22 XY tiles • 50-unit buffer • log1p axes • observed top-decile samples highlighted</text>'])
    stats={}
    for k,(name,pred) in enumerate(models.items()):
        col=k%2;row=k//2;left=95+col*555;top=110+row*350
        es=[pred[i]-actual[i] for i in rows];mae=statistics.fmean(abs(e) for e in es);rmse=math.sqrt(statistics.fmean(e*e for e in es));high=[i for i in rows if actual[i]>=Q90];rec=sum(pred[i]>=Q90 for i in high)/len(high)
        stats[name]={"n":len(rows),"mae_ppm":mae,"rmse_ppm":rmse,"top_decile_recall_at_prediction_cutoff":rec}
        body.append(f'<text class="title" x="{left}" y="{top-14}">{name}</text>')
        body.append(f'<rect x="{left}" y="{top}" width="{pw}" height="{ph}" fill="#fbfdff" stroke="#cbd5e1"/>')
        for val in ticks:
            pos=left+math.log1p(val)/fmax*pw;ypos=top+ph-math.log1p(val)/fmax*ph
            body.append(f'<line class="grid" x1="{pos:.2f}" y1="{top}" x2="{pos:.2f}" y2="{top+ph}"/><line class="grid" x1="{left}" y1="{ypos:.2f}" x2="{left+pw}" y2="{ypos:.2f}"/>')
            body.append(f'<text x="{pos:.2f}" y="{top+ph+17}" text-anchor="middle" font-size="10">{val//1000}k</text>')
            body.append(f'<text x="{left-8}" y="{ypos+4:.2f}" text-anchor="end" font-size="10">{val//1000}k</text>')
        body.append(f'<line class="diag" x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top}"/>')
        for val in (Q75,Q90):
            pos=left+math.log1p(val)/fmax*pw;ypos=top+ph-math.log1p(val)/fmax*ph
            body.append(f'<line class="thresh" x1="{pos:.2f}" y1="{top}" x2="{pos:.2f}" y2="{top+ph}"/><line class="thresh" x1="{left}" y1="{ypos:.2f}" x2="{left+pw}" y2="{ypos:.2f}"/>')
        for i in rows:
            x=left+math.log1p(actual[i])/fmax*pw;y=top+ph-math.log1p(max(0,pred[i]))/fmax*ph
            colr="#c2410c" if actual[i]>=Q90 else "#2563eb"
            body.append(f'<circle class="dot" cx="{x:.2f}" cy="{y:.2f}" r="2.2" fill="{colr}" fill-opacity=".42"><title>row {i}: observed {actual[i]:.0f} ppm; predicted {pred[i]:.0f} ppm</title></circle>')
        body.append(f'<text x="{left+pw/2}" y="{top+ph+37}" text-anchor="middle" font-size="12">Observed Cu (ppm; log1p)</text>')
        body.append(f'<text x="{left-52}" y="{top+ph/2}" text-anchor="middle" transform="rotate(-90 {left-52} {top+ph/2})" font-size="12">Predicted Cu (ppm; log1p)</text>')
        body.append(f'<text x="{left+pw-6}" y="{top+19}" text-anchor="end" font-size="11">MAE {mae:,.0f} • RMSE {rmse:,.0f} ppm • Q90 recall {rec:.1%}</text>')
    foot=110+350*nrows+10
    body.extend([f'<circle cx="370" cy="{foot}" r="4" fill="#2563eb"/><text x="380" y="{foot+4}" font-size="12">Other observed samples</text>',f'<circle cx="560" cy="{foot}" r="4" fill="#c2410c"/><text x="570" y="{foot+4}" font-size="12">Observed Cu ≥ 18,100 ppm</text>',f'<line x1="790" y1="{foot}" x2="820" y2="{foot}" class="diag"/><text x="830" y="{foot+4}" font-size="12">1:1 line</text>',f'<line x1="925" y1="{foot}" x2="955" y2="{foot}" class="thresh"/><text x="965" y="{foot+4}" font-size="12">Q1/Q9 threshold lines</text>',f'<text x="600" y="{foot+33}" text-anchor="middle" font-size="11" fill="#475569">Thresholds are descriptive; both chemistry-ridge variants assume non-Cu assays are available at target locations.</text>','</svg>'])
    FIG.write_text("\n".join(body),encoding="utf-8")
    meta={"figure":str(FIG.relative_to(ROOT)),"validation":"Same 22 buffered XY tiles, 300-unit tiles and 50-unit horizontal buffer.","axes":"log1p(Cu ppm)","thresholds_ppm":{"top_quartile":Q75,"top_decile":Q90},"models":stats,"caveat":"Pooled predicted-versus-observed values summarize out-of-fold calibration and show model shrinkage; fold-specific ranking metrics are in top_rank_screening.json. The all-assay ridge is conditional on co-located non-Cu assays."}
    META.write_text(json.dumps(meta,indent=2),encoding="utf-8")
    print(json.dumps(meta,indent=2))


if __name__=="__main__":main()
