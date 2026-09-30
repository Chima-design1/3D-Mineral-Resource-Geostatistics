"""Score continuous Cu predictions as screening flags at descriptive thresholds."""
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/"data"/"raw"/"drillholes.csv"
REPORTS=ROOT/"reports"
OUT=REPORTS/"grade_threshold_classification.json"
BLOCK,BUFFER,K,POWER=300.0,50.0,12,2.0
THRESHOLDS={"top_quartile_descriptive":8900.0,"top_decile_descriptive":18100.0}


def read_csv():
    with RAW.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))


def buffered_splits(rows):
    pts=[tuple(float(r[c]) for c in ("X","Y","Z")) for r in rows]
    x0=min(p[0] for p in pts);y0=min(p[1] for p in pts);groups=defaultdict(list)
    for i,p in enumerate(pts):groups[(math.floor((p[0]-x0)/BLOCK),math.floor((p[1]-y0)/BLOCK))].append(i)
    out=[]
    for (bx,by),test in sorted(groups.items()):
        left,bottom=x0+bx*BLOCK,y0+by*BLOCK;right,top=left+BLOCK,bottom+BLOCK;ts=set(test);train=[]
        for i,p in enumerate(pts):
            if i in ts:continue
            dx=max(left-p[0],0,p[0]-right);dy=max(bottom-p[1],0,p[1]-top)
            if math.hypot(dx,dy)>BUFFER:train.append(i)
        out.append((test,train))
    return pts,out


def raw_idw_predictions(rows,pts,splits):
    predictions={};means={}
    for test,train in splits:
        if not train:continue
        mean=statistics.fmean(float(rows[i]["Cu ppm"]) for i in train)
        for i in test:
            near=sorted((math.dist(pts[i],pts[j]),j) for j in train)[:K]
            if near[0][0]==0:pred=float(rows[near[0][1]]["Cu ppm"])
            else:
                ws=[1/(d**POWER) for d,_ in near]
                pred=sum(w*float(rows[j]["Cu ppm"]) for w,(_,j) in zip(ws,near))/sum(ws)
            predictions[i+1]=pred;means[i+1]=mean
    return predictions,means


def map_report_predictions(filename,pred_key):
    data=json.loads((REPORTS/filename).read_text(encoding="utf-8"));out={}
    for row in data["outer_predictions"]:out[int(row["row"])]=float(row[pred_key])
    return out


def confusion(actual,scores,threshold):
    observed=[y>=threshold for y in actual];predicted=[p>=threshold for p in scores]
    tp=sum(a and p for a,p in zip(observed,predicted));fp=sum((not a) and p for a,p in zip(observed,predicted));fn=sum(a and (not p) for a,p in zip(observed,predicted));tn=len(actual)-tp-fp-fn
    precision=tp/(tp+fp) if tp+fp else 0.0;recall=tp/(tp+fn) if tp+fn else 0.0;specificity=tn/(tn+fp) if tn+fp else 0.0
    f1=2*precision*recall/(precision+recall) if precision+recall else 0.0
    return {"threshold_ppm":threshold,"n":len(actual),"actual_positive_n":tp+fn,"predicted_positive_n":tp+fp,"tp":tp,"fp":fp,"fn":fn,"tn":tn,"precision":precision,"recall":recall,"specificity":specificity,"f1":f1,"balanced_accuracy":(recall+specificity)/2}


def average_precision(actual,scores,threshold):
    labels=[y>=threshold for y in actual];positives=sum(labels)
    if not positives:return None
    ranks=sorted(zip(scores,labels),reverse=True);tp=fp=0;prev_recall=0.0;ap=0.0;i=0
    while i<len(ranks):
        score=ranks[i][0];j=i
        while j<len(ranks) and ranks[j][0]==score:
            if ranks[j][1]:tp+=1
            else:fp+=1
            j+=1
        recall=tp/positives;precision=tp/(tp+fp)
        ap+=(recall-prev_recall)*precision;prev_recall=recall;i=j
    return ap


def main():
    rows=read_csv();actual=[float(r["Cu ppm"]) for r in rows];pts,splits=buffered_splits(rows)
    pred,mean=raw_idw_predictions(rows,pts,splits)
    candidates={"raw_Cu_IDW_p2_k12":pred,"sample_weighted_training_mean":mean}
    candidates["nested_ordinary_kriging"]=map_report_predictions("nested_kriging_sensitivity.json","predicted_cu_ppm")
    candidates["conditional_all_assay_ridge"]=map_report_predictions("geochem_covariate_cv.json","predicted_cu_ppm")
    feat=json.loads((REPORTS/"geochem_feature_sensitivity.json").read_text(encoding="utf-8"))["models"]["all_17"]["outer_predictions"]
    candidates["conditional_all_assay_ridge_alpha_100_to_10000"]={int(r["row"]):float(r["predicted_cu_ppm"]) for r in feat}
    sub=json.loads((REPORTS/"subcomposition_spatial_sensitivity.json").read_text(encoding="utf-8"))["designs"]["buffered_XY_tile"]["predictions_by_row"]
    candidates["Cu_Au_Ag_S_subcomposition_IDW"]={int(r["row"]):float(r["clr_idw_total_idw_zero_0.1x_min"]) for r in sub}
    expected=set(range(1,len(rows)+1)); coverage={name:len(set(p)&expected) for name,p in candidates.items()}
    out={"validation":"All methods scored on the same 22 buffered 300-unit XY tile holdouts with a 50-unit horizontal training buffer.","descriptive_grade_thresholds_ppm":THRESHOLDS,"threshold_warning":"These thresholds are distribution-based (top quartile/top decile), not economic cutoffs, ore/waste boundaries, or resource classifications.","model_coverage_rows":coverage,"models":{},"interpretation":"Classification metrics assess whether a predicted Cu concentration crosses a descriptive sample threshold. They do not assess economic value or geological continuity. The conditional assay ridge assumes all non-Cu assays are available at held-out locations."}
    for name,mapping in candidates.items():
        scores=[mapping[i+1] for i in range(len(actual))]
        out["models"][name]={key:{**confusion(actual,scores,thr),"average_precision":average_precision(actual,scores,thr)} for key,thr in THRESHOLDS.items()}
    OUT.write_text(json.dumps(out,indent=2),encoding="utf-8")
    for name,res in out["models"].items():print(name,json.dumps(res,indent=2))


if __name__=="__main__":main()
