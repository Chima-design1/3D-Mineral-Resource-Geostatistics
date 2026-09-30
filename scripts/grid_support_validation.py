"""Measure buffered spatial-CV IDW error by geometric support and cutoff."""
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "drillholes.csv"
SUPPORT = ROOT / "reports" / "spatial_support_audit.json"
OUT = ROOT / "reports" / "grid_support_validation.json"
TILE, BUFFER, K, POWER = 300.0, 50.0, 12, 2.0


def metrics(items):
    if not items: return {"n":0,"mae_ppm":None,"rmse_ppm":None}
    es=[p-y for y,p in items]
    return {"n":len(es),"mae_ppm":statistics.fmean(map(abs,es)),"rmse_ppm":math.sqrt(statistics.fmean(e*e for e in es))}


def main():
    with DATA.open(encoding="utf-8-sig",newline="") as f: rows=list(csv.DictReader(f))
    pts=[(r["HOLEID"],tuple(float(r[c]) for c in ("X","Y","Z")),float(r["Cu ppm"])) for r in rows]
    sup=json.loads(SUPPORT.read_text(encoding="utf-8"))
    near0=sup["nearest_sample_spacing_units"]["3d"]["q3"]
    kth0=2*sup["nearest_sample_spacing_units"]["nearest_sample_from_different_hole_3d"]["median"]
    xs=[p[1][0] for p in pts]; ys=[p[1][1] for p in pts]; x0,y0=min(xs),min(ys)
    tiles=defaultdict(list)
    for i,p in enumerate(pts): tiles[(math.floor((p[1][0]-x0)/TILE),math.floor((p[1][1]-y0)/TILE))].append(i)
    records=[]
    for (bx,by), test in sorted(tiles.items()):
        left,bottom=x0+bx*TILE,y0+by*TILE; right,top=left+TILE,bottom+TILE; testset=set(test); train=[]
        for i,p in enumerate(pts):
            if i in testset: continue
            x,y,_=p[1]; dx=max(left-x,0,x-right); dy=max(bottom-y,0,y-top)
            if math.hypot(dx,dy)>BUFFER: train.append(i)
        if len(train)<K: continue
        mean=statistics.fmean(pts[i][2] for i in train)
        for i in test:
            hole,xyz,y=pts[i]
            near=sorted((math.dist(xyz,pts[j][1]),j) for j in train)[:K]
            ds=[v[0] for v in near]
            if ds[0]==0: pred=pts[near[0][1]][2]
            else:
                ws=[1/(d**POWER) for d in ds]
                pred=sum(w*pts[j][2] for w,(_,j) in zip(ws,near))/sum(ws)
            records.append({"row":i+1,"fold":[bx,by],"actual":y,"idw":pred,"mean":mean,"nearest":ds[0],"kth":ds[-1]})
    combos=[]
    near_limits=[near0,2*near0,max(sup["nearest_sample_spacing_units"]["3d"]["max"],2*near0)]
    kth_limits=[kth0/2,kth0,1.5*kth0]
    for a in near_limits:
        for b in kth_limits:
            subset=[r for r in records if r["nearest"]<=a and r["kth"]<=b]
            combos.append({"max_nearest_distance":a,"max_12th_neighbor_distance":b,"idw":metrics([(r["actual"],r["idw"]) for r in subset]),"training_mean":metrics([(r["actual"],r["mean"]) for r in subset]),"heldout_fraction_percent":100*len(subset)/len(records)})
    out={"validation":"Each occupied 300-unit XY tile held out in turn; training points within 50 horizontal units of that tile are excluded. IDW p=2, k=12.","n_heldout_predictions":len(records),"folds":len(tiles),"distance_cutoffs_selected_without_grade_error_tuning":True,"baseline_cutoffs_used_for_current_grid":{"nearest":near0,"12th_neighbor":kth0},"all_holdout_metrics":{"idw":metrics([(r["actual"],r["idw"]) for r in records]),"training_mean":metrics([(r["actual"],r["mean"]) for r in records])},"support_band_metrics":combos,"interpretation":"Support-conditioned scores describe interpolation accuracy at held-out sample locations under these splits. They do not certify unsampled nodes, correct geology, or resource confidence.","records":records}
    OUT.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in out.items() if k!="records"},indent=2))


if __name__=="__main__": main()
