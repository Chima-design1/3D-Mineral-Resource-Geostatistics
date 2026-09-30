"""Check baseline validation sensitivity to XY tile width and buffer."""
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"/"raw"/"drillholes.csv"
OUT=ROOT/"reports"/"spatial_cv_design_sensitivity.json"
TILE_WIDTHS=(200.0,300.0,500.0)
BUFFERS=(0.0,50.0,100.0)
K,POWER,Q75,Q90=12,2.0,8900.0,18100.0


def metrics(actual,pred):
    e=[p-y for y,p in zip(actual,pred)]
    hi=[i for i,y in enumerate(actual) if y>=Q75]
    return {"n":len(actual),"mae_ppm":statistics.fmean(abs(v) for v in e),"rmse_ppm":math.sqrt(statistics.fmean(v*v for v in e)),"mean_error_ppm":statistics.fmean(e),"high_grade_n":len(hi),"high_grade_mae_ppm":statistics.fmean(abs(pred[i]-actual[i]) for i in hi) if hi else None}


def evaluate(samples,width,buffer):
    xs=[r[1][0] for r in samples];ys=[r[1][1] for r in samples];x0,y0=min(xs),min(ys);tiles=defaultdict(list)
    for i,(_,p,_) in enumerate(samples):tiles[(math.floor((p[0]-x0)/width),math.floor((p[1]-y0)/width))].append(i)
    actual=[];idw=[];mean=[];folds=[];skipped=0
    for (bx,by),test in sorted(tiles.items()):
        left,bottom=x0+bx*width,y0+by*width;right,top=left+width,bottom+width;ts=set(test);train=[]
        for i,(_,p,_) in enumerate(samples):
            if i in ts:continue
            dx=max(left-p[0],0,p[0]-right);dy=max(bottom-p[1],0,p[1]-top)
            if math.hypot(dx,dy)>buffer:train.append(i)
        if len(train)<K:
            skipped+=1;continue
        ys=[samples[i][2] for i in train];m=statistics.fmean(ys);fy=[];fp=[];fm=[]
        for i in test:
            y=samples[i][2];nearest=sorted((math.dist(samples[i][1],samples[j][1]),j) for j in train)[:K]
            if nearest[0][0]==0:p=samples[nearest[0][1]][2]
            else:
                w=[1/(d**POWER) for d,_ in nearest];p=sum(a*samples[j][2] for a,(_,j) in zip(w,nearest))/sum(w)
            actual.append(y);idw.append(p);mean.append(m);fy.append(y);fp.append(p);fm.append(m)
        folds.append({"tile":[bx,by],"test_n":len(test),"train_n":len(train),"idw":metrics(fy,fp),"mean":metrics(fy,fm)})
    return {"tile_width_coordinate_units":width,"buffer_coordinate_units":buffer,"occupied_tiles":len(tiles),"evaluated_folds":len(folds),"skipped_folds_fewer_than_12_train_samples":skipped,"evaluated_n":len(actual),"coverage_percent":100*len(actual)/len(samples),"idw":metrics(actual,idw),"training_mean":metrics(actual,mean),"folds":folds}


def main():
    with DATA.open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
    samples=[(r["HOLEID"],tuple(float(r[c]) for c in ("X","Y","Z")),float(r["Cu ppm"])) for r in rows]
    results=[evaluate(samples,t,b) for t in TILE_WIDTHS for b in BUFFERS]
    out={"target":"Cu ppm","method":"3D IDW (p=2, k=12) versus sample-weighted training mean","coordinates":"Source local Cartesian values; units and CRS not documented in CSV.","designs":results,"caveats":["This varies the validation design for fixed, untuned baselines; it does not retune kriging or conditional chemistry models.","Tile width and buffer values are in undocumented source coordinate units.","Occupied XY tiles may not be independent, and this sensitivity is not an uncertainty or resource-classification analysis."]}
    OUT.write_text(json.dumps(out,indent=2),encoding="utf-8")
    for r in results:print(r["tile_width_coordinate_units"],r["buffer_coordinate_units"],r["evaluated_n"],"folds",r["evaluated_folds"],"IDW MAE",round(r["idw"]["mae_ppm"]),"mean MAE",round(r["training_mean"]["mae_ppm"]),"IDW RMSE",round(r["idw"]["rmse_ppm"]),"mean RMSE",round(r["training_mean"]["rmse_ppm"]))


if __name__=="__main__":main()
