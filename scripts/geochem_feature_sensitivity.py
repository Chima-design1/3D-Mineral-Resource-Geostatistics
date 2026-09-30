"""Nested spatial sensitivity: all chemistry versus S/Ag/Au covariates."""
import csv, json, math, statistics
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import fold_kriging as fk
import nested_kriging_sensitivity as nk

ALL=["Ag ppm","Al ppm","Au ppm","C ppm","Ca ppm","Cl ppm","F ppm","Fe ppm","K ppm","Mg ppm","Mn ppm","Na ppm","P ppm","Pb ppm","S ppm","Th ppm","U ppm"]
SETS={"all_17":ALL,"S_Ag_Au":["S ppm","Ag ppm","Au ppm"]}
ALPHAS=(100.0,300.0,1000.0,3000.0,10000.0)

def features(samples, ids, cols, means=None, scales=None):
    x0=[[math.log1p(samples[i]["chem"][c]) for c in cols] for i in ids]
    if means is None:
        means=[statistics.fmean(row[j] for row in x0) for j in range(len(cols))]
        scales=[statistics.pstdev(row[j] for row in x0) or 1.0 for j in range(len(cols))]
    return [[(r[j]-means[j])/scales[j] for j in range(len(cols))] for r in x0],means,scales

def fit(samples, train, cols, alpha):
    x,means,scales=features(samples,train,cols)
    y=[samples[i]["cu"] for i in train]; ym=statistics.fmean(y); yc=[v-ym for v in y]
    p=len(cols)
    xtx=[[sum(row[j]*row[k] for row in x)+(alpha if j==k else 0.0) for k in range(p)] for j in range(p)]
    xty=[sum(row[j]*v for row,v in zip(x,yc)) for j in range(p)]
    return {"means":means,"scales":scales,"ym":ym,"beta":fk.solve(xtx,xty)}

def predict(samples, ids, cols, model):
    x,_,_=features(samples,ids,cols,model["means"],model["scales"])
    return [max(0.0,model["ym"]+sum(a*b for a,b in zip(row,model["beta"]))) for row in x]

def metrics(y,p):
    e=[b-a for a,b in zip(y,p)]
    return {"n":len(y),"mae_ppm":statistics.fmean(abs(v) for v in e),"rmse_ppm":math.sqrt(statistics.fmean(v*v for v in e)),"mean_error_ppm":statistics.fmean(e)}

def main():
    samples=fk.load()
    with fk.DATA.open("r",encoding="utf-8-sig",newline="") as f: rows=list(csv.DictReader(f))
    for i,r in enumerate(samples): r["chem"]={c:float(rows[i][c]) for c in ALL}
    outer=fk.grid_folds(samples); x0=min(s["xyz"][0] for s in samples); y0=min(s["xyz"][1] for s in samples)
    results={name:{"y":[],"p":[],"folds":[],"predictions":[]} for name in SETS}
    means_y=[]; means_p=[]
    for ofold,(test,train) in enumerate(outer):
        splits=nk.inner_splits(train,samples,x0,y0)
        for name,cols in SETS.items():
            inner_mae={a:[] for a in ALPHAS}
            for valid,itrain in splits:
                for alpha in ALPHAS:
                    m=fit(samples,itrain,cols,alpha); pp=predict(samples,valid,cols,m)
                    inner_mae[alpha].extend(abs(p-samples[i]["cu"]) for p,i in zip(pp,valid))
            scores={a:statistics.fmean(es) for a,es in inner_mae.items()}
            selected=min(scores,key=scores.get)
            m=fit(samples,train,cols,selected); pp=predict(samples,test,cols,m); yy=[samples[i]["cu"] for i in test]
            results[name]["y"].extend(yy); results[name]["p"].extend(pp)
            results[name]["folds"].append({"outer_fold":ofold+1,"selected_alpha":selected,"inner_mae_by_alpha_ppm":scores,"outer_metrics":metrics(yy,pp)})
            results[name]["predictions"].extend({"row":i+1,"observed_cu_ppm":y,"predicted_cu_ppm":p,"residual_ppm":p-y,"outer_fold":ofold+1} for i,y,p in zip(test,yy,pp))
        avg=statistics.fmean(samples[i]["cu"] for i in train)
        yy=[samples[i]["cu"] for i in test];means_y.extend(yy);means_p.extend([avg]*len(yy))
        print(f"completed sensitivity fold {ofold+1}/{len(outer)}",flush=True)
    all_y=results["all_17"]["y"]; q75=sorted(all_y)[math.ceil(.75*len(all_y))-1]; high=[y>=q75 for y in all_y]
    summary={}
    for name,res in results.items():
        summary[name]={"features":SETS[name],"overall":metrics(res["y"],res["p"]),"top_quartile_threshold_ppm":q75,"top_quartile":metrics([y for y,h in zip(res["y"],high) if h],[p for p,h in zip(res["p"],high) if h]),"selected_alpha_counts":{str(a):sum(f["selected_alpha"]==a for f in res["folds"]) for a in ALPHAS},"folds":res["folds"],"outer_predictions":res["predictions"]}
    summary["same_split_training_mean"]={"overall":metrics(means_y,means_p),"top_quartile":metrics([y for y,h in zip(means_y,high) if h],[p for p,h in zip(means_p,high) if h])}
    report={"design":"Same 22 outer 300-unit XY tiles with 50-unit buffer. Each feature set tunes ridge alpha inside two distributed inner buffered tiles; transforms and scaling are fitted within each training fold.","alpha_candidates":list(ALPHAS),"models":summary,"limitations":["The S/Ag/Au shortlist was motivated by prior descriptive whole-file correlations; treat that model comparison as exploratory feature sensitivity.","Evaluation assumes companion assays are available at outer held-out locations while Cu is withheld; this is conditional prediction or assay-imputation, not estimation where all chemistry is absent.","Inner tuning uses only two inner spatial tiles per outer fold.","Neither covariate model defines geological domains or a resource estimate."]}
    out=ROOT/"reports"/"geochem_feature_sensitivity.json";out.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps({"all_17":summary["all_17"]["overall"],"S_Ag_Au":summary["S_Ag_Au"]["overall"],"mean":summary["same_split_training_mean"]["overall"],"high_quartile":{"all_17":summary["all_17"]["top_quartile"],"S_Ag_Au":summary["S_Ag_Au"]["top_quartile"],"mean":summary["same_split_training_mean"]["top_quartile"]},"alpha_counts":{"all_17":summary["all_17"]["selected_alpha_counts"],"S_Ag_Au":summary["S_Ag_Au"]["selected_alpha_counts"]},"report":out.name},indent=2))

if __name__=="__main__": main()
