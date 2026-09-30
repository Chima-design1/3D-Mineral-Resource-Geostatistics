"""Spatial CV sensitivity for a Cu-Au-Ag-S subcomposition with zero replacement."""
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"/"raw"/"drillholes.csv"
OUT=ROOT/"reports"/"subcomposition_spatial_sensitivity.json"
ELEMENTS=("Cu ppm","Au ppm","Ag ppm","S ppm")
K,POWER,BLOCK,BUFFER=12,2.0,300.0,50.0
ZERO_FRACTIONS=(0.1,0.5)
HIGH_GRADE_CUTOFF=8900.0


def load():
    with DATA.open(encoding="utf-8-sig",newline="") as f: rows=list(csv.DictReader(f))
    return [{"hole":r["HOLEID"],"xyz":tuple(float(r[c]) for c in ("X","Y","Z")),"parts":[float(r[c]) for c in ELEMENTS]} for r in rows]


def loho(samples):
    groups=defaultdict(list)
    for i,s in enumerate(samples):groups[s["hole"]].append(i)
    return [(ids,[i for i,s in enumerate(samples) if s["hole"]!=hole]) for hole,ids in groups.items()]


def spatial(samples):
    x0=min(s["xyz"][0] for s in samples); y0=min(s["xyz"][1] for s in samples); tiles=defaultdict(list)
    for i,s in enumerate(samples):tiles[(math.floor((s["xyz"][0]-x0)/BLOCK),math.floor((s["xyz"][1]-y0)/BLOCK))].append(i)
    out=[]
    for (bx,by),test in sorted(tiles.items()):
        left,bottom=x0+bx*BLOCK,y0+by*BLOCK; right,top=left+BLOCK,bottom+BLOCK; ts=set(test); train=[]
        for i,s in enumerate(samples):
            if i in ts:continue
            x,y,_=s["xyz"]; dx=max(left-x,0,x-right); dy=max(bottom-y,0,y-top)
            if math.hypot(dx,dy)>BUFFER:train.append(i)
        if len(train)>=K:out.append((test,train))
    return out


def transform_training(samples,train,fraction):
    positive_min=[]
    for j in range(4):
        vals=[samples[i]["parts"][j] for i in train if samples[i]["parts"][j]>0]
        positive_min.append(min(vals))
    deltas=[fraction*x for x in positive_min]
    clr={}; logtotal={}
    for i in train:
        adj=[v if v>0 else deltas[j] for j,v in enumerate(samples[i]["parts"])]
        logs=[math.log(v) for v in adj]; center=statistics.fmean(logs)
        clr[i]=[v-center for v in logs]
        logtotal[i]=math.log(sum(adj))
    return deltas,clr,logtotal


def idw(values,neighbors):
    exact=[values[j] for d,j in neighbors if d==0]
    if exact:return exact[0] if not isinstance(exact[0],list) else exact[0][:]
    weights=[1/(d**POWER) for d,_ in neighbors]; sw=sum(weights)
    first=values[neighbors[0][1]]
    if isinstance(first,list):return [sum(w*values[j][k] for w,(_,j) in zip(weights,neighbors))/sw for k in range(len(first))]
    return sum(w*values[j] for w,(_,j) in zip(weights,neighbors))/sw


def predict_composition(clr,logtotal,neighbors):
    c=idw(clr,neighbors); lt=idw(logtotal,neighbors)
    es=[math.exp(v-max(c)) for v in c]; fractions=[v/sum(es) for v in es]
    return fractions[0]*math.exp(lt)


def metrics(actual,pred):
    es=[p-y for y,p in zip(actual,pred)]
    return {"n":len(es),"mae_ppm":statistics.fmean(abs(e) for e in es),"rmse_ppm":math.sqrt(statistics.fmean(e*e for e in es)),"mean_error_ppm":statistics.fmean(es)}


def summarize(actual,preds):
    out={name:metrics(actual,p) for name,p in preds.items()}
    ix=[i for i,y in enumerate(actual) if y>=HIGH_GRADE_CUTOFF]
    out["high_grade_cutoff_ppm"]=HIGH_GRADE_CUTOFF
    out["high_grade_n"]=len(ix)
    out["high_grade_mae_ppm"]={name:statistics.fmean(abs(pred[i]-actual[i]) for i in ix) for name,pred in preds.items()} if ix else {}
    return out


def evaluate(samples,splits,design):
    actual=[]; preds={"raw_cu_idw_p2_k12":[],"training_mean":[]}; preds.update({f"clr_idw_total_idw_zero_{f:g}x_min":[] for f in ZERO_FRACTIONS})
    deltas_by_fraction={f:[] for f in ZERO_FRACTIONS}; fold_ids=[]
    for fold,(test,train) in enumerate(splits,1):
        raw_train={i:samples[i]["parts"][0] for i in train}
        means=statistics.fmean(raw_train.values())
        transforms={f:transform_training(samples,train,f) for f in ZERO_FRACTIONS}
        for f,(deltas,clr,lt) in transforms.items():deltas_by_fraction[f].append(deltas)
        for i in test:
            near=sorted((math.dist(samples[i]["xyz"],samples[j]["xyz"]),j) for j in train)[:K]
            vals=[samples[j]["parts"][0] for _,j in near]
            ds=[d for d,_ in near]
            if ds[0]==0:rawpred=vals[0]
            else:
                ws=[1/(d**POWER) for d in ds];rawpred=sum(w*v for w,v in zip(ws,vals))/sum(ws)
            y=samples[i]["parts"][0];actual.append(y);preds["raw_cu_idw_p2_k12"].append(rawpred);preds["training_mean"].append(means)
            for f,(deltas,clr,lt) in transforms.items():
                preds[f"clr_idw_total_idw_zero_{f:g}x_min"].append(predict_composition(clr,lt,near))
            fold_ids.append(fold)
    delta_summary={}
    for f,fold_deltas in deltas_by_fraction.items():
        delta_summary[str(f)]={ELEMENTS[j]:{"min_across_folds":min(d[j] for d in fold_deltas),"median_across_folds":statistics.median(d[j] for d in fold_deltas),"max_across_folds":max(d[j] for d in fold_deltas)} for j in range(4)}
    return {"n_predictions":len(actual),"folds":len(splits),"overall":summarize(actual,preds),"zero_replacement_deltas_ppm_by_training_fold":delta_summary,"predictions_by_row":[{"row":i+1,"fold":fold_ids[k],"actual_cu_ppm":actual[k],**{name:pred[k] for name,pred in preds.items()}} for k,i in enumerate([i for test,_ in splits for i in test])],"validation_design":design}


def main():
    samples=load()
    out={"purpose":"Test a four-part Cu-Au-Ag-S subcomposition transform as an exploratory Cu prediction alternative; not a whole-rock mass composition.","components":list(ELEMENTS),"method":{"zero_replacement":"Replace each exact zero by fraction times the minimum positive value for that component in the training fold; leave positive values unchanged, then close four parts to sum one.","fractions_of_training_minimum":list(ZERO_FRACTIONS),"transforms":"Interpolate each centered log-ratio coordinate and the natural log of the adjusted four-part total with 3D IDW (p=2, k=12); inverse CLR by softmax and multiply by estimated total to return Cu ppm.","hyperparameters":"Fixed in advance; no tuning on holdout errors.","assay_zero_semantics":"Unknown; scenarios are hypothetical and not detection-limit-based."},"designs":{}}
    out["designs"]["leave_one_HOLEID_out"]=evaluate(samples,loho(samples),"Leave each complete HOLEID out; all samples from a held-out hole excluded.")
    out["designs"]["buffered_XY_tile"]=evaluate(samples,spatial(samples),"Leave each occupied 300-unit XY tile out; remove training points within 50 horizontal units; 3D prediction distances.")
    OUT.write_text(json.dumps(out,indent=2),encoding="utf-8")
    for name,result in out["designs"].items():print(name,json.dumps(result["overall"],indent=2))


if __name__=="__main__":main()
