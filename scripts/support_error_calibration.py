"""Describe buffered IDW errors across the full holdout distance distribution."""
import json
import math
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"reports"/"grid_support_validation.json"
OUT=ROOT/"reports"/"support_error_calibration.json"
FIG=ROOT/"reports"/"support_error_calibration.svg"


def decile_rows(records, key, groups=10):
    ordered=sorted(records,key=lambda r:r[key])
    out=[]
    for b in range(groups):
        chunk=ordered[b*len(ordered)//groups:(b+1)*len(ordered)//groups]
        if not chunk: continue
        def score(pred):
            e=[r[pred]-r["actual"] for r in chunk]
            return {"mae_ppm":statistics.fmean(map(abs,e)),"rmse_ppm":math.sqrt(statistics.fmean(v*v for v in e))}
        out.append({"bin":b+1,"n":len(chunk),"distance_min":chunk[0][key],"distance_median":statistics.median(r[key] for r in chunk),"distance_max":chunk[-1][key],"idw":score("idw"),"training_mean":score("mean")})
    return out


def make_svg(results):
    w,h=1100,500; panels=[("Nearest training sample (3D)",results["nearest_distance_deciles"]),("12th-nearest training sample (3D)",results["12th_neighbor_distance_deciles"])]
    body=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">','<rect width="100%" height="100%" fill="#fbfcfe"/>','<style>text{font-family:Arial,sans-serif;fill:#263445}.title{font-size:19px;font-weight:600}.axis{stroke:#64748b;stroke-width:1}.idw{fill:none;stroke:#1769aa;stroke-width:3}.mean{fill:none;stroke:#d97706;stroke-width:3}.grid{stroke:#dbe3ec;stroke-width:1}</style>']
    for pi,(title,bins) in enumerate(panels):
        x0=72+pi*530; y0=110; pw=445; ph=290
        ymax=max(max(r["idw"]["mae_ppm"],r["training_mean"]["mae_ppm"]) for r in bins)*1.08
        body.append(f'<text class="title" x="{x0}" y="55">{title}</text>')
        for tick in range(5):
            y=y0+ph-tick*ph/4; val=ymax*tick/4
            body.append(f'<line class="grid" x1="{x0}" y1="{y:.1f}" x2="{x0+pw}" y2="{y:.1f}"/><text x="{x0-8}" y="{y+4:.1f}" text-anchor="end" font-size="11">{val:.0f}</text>')
        def coords(method):
            return [(x0+(i+0.5)*pw/len(bins),y0+ph-r[method]["mae_ppm"]/ymax*ph) for i,r in enumerate(bins)]
        for name,cls,col in (("idw","idw","#1769aa"),("training_mean","mean","#d97706")):
            xy=coords(name); body.append(f'<polyline class="{cls}" points="'+' '.join(f'{x:.1f},{y:.1f}' for x,y in xy)+'"/>')
            for x,y in xy: body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{col}"/>')
        for i,r in enumerate(bins):
            x=x0+(i+0.5)*pw/len(bins); body.append(f'<text x="{x:.1f}" y="{y0+ph+19}" text-anchor="middle" font-size="10">{r["n"]}</text>')
        body.append(f'<text x="{x0+pw/2}" y="{y0+ph+43}" text-anchor="middle" font-size="12">Distance decile (number below)</text>')
        body.append(f'<text x="{x0-48}" y="{y0+ph/2}" text-anchor="middle" transform="rotate(-90 {x0-48} {y0+ph/2})" font-size="12">MAE (Cu ppm)</text>')
    body.extend(['<line x1="430" y1="462" x2="460" y2="462" stroke="#1769aa" stroke-width="3"/><text x="467" y="466" font-size="12">IDW p=2, k=12</text>','<line x1="610" y1="462" x2="640" y2="462" stroke="#d97706" stroke-width="3"/><text x="647" y="466" font-size="12">Training mean baseline</text>','<text x="550" y="492" text-anchor="middle" font-size="11">Buffered 300-unit XY tile holdouts; deciles have near-equal sample counts. Descriptive calibration, not an independent test.</text>','</svg>'])
    FIG.write_text("\n".join(body),encoding="utf-8")


def main():
    data=json.loads(SOURCE.read_text(encoding="utf-8")); records=data["records"]
    out={"source_report":SOURCE.name,"validation_design":data["validation"],"n":len(records),"binning":"Ten equal-count bins, ordered separately by each distance measure.","nearest_distance_deciles":decile_rows(records,"nearest"),"12th_neighbor_distance_deciles":decile_rows(records,"kth"),"interpretation":"Binned results are descriptive because the same buffered holdouts are used to display all distance bands; do not choose and report a cutoff as independently validated from these results.","figure":FIG.name}
    OUT.write_text(json.dumps(out,indent=2),encoding="utf-8"); make_svg(out)
    print(json.dumps({k:v for k,v in out.items() if "deciles" not in k},indent=2))
    for key in ("nearest_distance_deciles","12th_neighbor_distance_deciles"):
        print(key)
        for row in out[key]: print(row["bin"],row["n"],round(row["distance_min"],1),round(row["distance_max"],1),round(row["idw"]["mae_ppm"]),round(row["training_mean"]["mae_ppm"]))


if __name__=="__main__": main()
