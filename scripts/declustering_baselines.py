"""Compare sample-, hole-, and occupied-XY-tile-weighted mean baselines."""
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"/"raw"/"drillholes.csv"
OUT=ROOT/"reports"/"declustering_baselines.json"
BLOCK,BUFFER,TOP_CU=300.0,50.0,8900.0


def load():
    with DATA.open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
    return [{"hole":r["HOLEID"],"xyz":tuple(float(r[c]) for c in ("X","Y","Z")),"cu":float(r["Cu ppm"])} for r in rows]


def splits_loho(samples):
    groups=defaultdict(list)
    for i,s in enumerate(samples):groups[s["hole"]].append(i)
    return [(idx,[i for i,s in enumerate(samples) if s["hole"]!=hole]) for hole,idx in groups.items()]


def splits_xy(samples):
    x0=min(s["xyz"][0] for s in samples);y0=min(s["xyz"][1] for s in samples); groups=defaultdict(list)
    for i,s in enumerate(samples):groups[(math.floor((s["xyz"][0]-x0)/BLOCK),math.floor((s["xyz"][1]-y0)/BLOCK))].append(i)
    out=[]
    for (bx,by),test in sorted(groups.items()):
        left,bottom=x0+bx*BLOCK,y0+by*BLOCK;right,top=left+BLOCK,bottom+BLOCK;ts=set(test);train=[]
        for i,s in enumerate(samples):
            if i in ts:continue
            x,y,_=s["xyz"];dx=max(left-x,0,x-right);dy=max(bottom-y,0,y-top)
            if math.hypot(dx,dy)>BUFFER:train.append(i)
        if train:out.append((test,train))
    return out


def weighted_means(samples,train,tile_origin):
    ys=[samples[i]["cu"] for i in train]
    holes=defaultdict(list);tiles=defaultdict(list)
    for i in train:
        s=samples[i];holes[s["hole"]].append(s["cu"])
        bx=math.floor((s["xyz"][0]-tile_origin[0])/BLOCK);by=math.floor((s["xyz"][1]-tile_origin[1])/BLOCK)
        tiles[(bx,by)].append(s["cu"])
    return {"sample_weighted":statistics.fmean(ys),"equal_hole":statistics.fmean(statistics.fmean(v) for v in holes.values()),"equal_xy_tile":statistics.fmean(statistics.fmean(v) for v in tiles.values())}


def metrics(actual,pred):
    e=[p-y for y,p in zip(actual,pred)]
    high=[i for i,y in enumerate(actual) if y>=TOP_CU]
    return {"n":len(actual),"mae_ppm":statistics.fmean(abs(x) for x in e),"rmse_ppm":math.sqrt(statistics.fmean(x*x for x in e)),"mean_error_ppm":statistics.fmean(e),"high_grade_cutoff_ppm":TOP_CU,"high_grade_n":len(high),"high_grade_mae_ppm":statistics.fmean(abs(pred[i]-actual[i]) for i in high) if high else None}


def evaluate(samples,splits,name,origin):
    actual=[];preds={k:[] for k in ("sample_weighted","equal_hole","equal_xy_tile")};folds=[]
    for n,(test,train) in enumerate(splits,1):
        estimates=weighted_means(samples,train,origin);ys=[]
        for i in test:
            y=samples[i]["cu"];ys.append(y);actual.append(y)
            for key in preds:preds[key].append(estimates[key])
        folds.append({"fold":n,"test_n":len(test),"train_n":len(train),"training_means_ppm":estimates,"test_mean_ppm":statistics.fmean(ys)})
    return {"validation":name,"folds":len(folds),"overall":{k:metrics(actual,v) for k,v in preds.items()},"fold_details":folds}


def main():
    samples=load();origin=(min(s["xyz"][0] for s in samples),min(s["xyz"][1] for s in samples))
    out={"target":"Cu ppm","occupied_tile_width":BLOCK,"buffer":BUFFER,"tile_grid_origin_xy":origin,"interpretation":"Equal-hole and equal-occupied-tile means are declustering baselines, not kriging estimators. Category means are computed from training rows only within each fold.","designs":{}}
    out["designs"]["leave_one_HOLEID_out"]=evaluate(samples,splits_loho(samples),"Leave one HOLEID group out",origin)
    out["designs"]["buffered_XY_tile"]=evaluate(samples,splits_xy(samples),"Leave each occupied 300-unit XY tile out; 50-unit horizontal training buffer",origin)
    OUT.write_text(json.dumps(out,indent=2),encoding="utf-8")
    for design,res in out["designs"].items():print(design,json.dumps(res["overall"],indent=2))


if __name__=="__main__":main()
