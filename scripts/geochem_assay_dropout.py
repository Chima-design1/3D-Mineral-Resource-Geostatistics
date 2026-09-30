"""Nested spatial stress test when secondary analytes are absent at targets."""
import csv,json,math,statistics,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import fold_kriging as fk
import nested_kriging_sensitivity as nk
import geochem_feature_sensitivity as gf

SCENARIOS={"none":[],"S_missing":["S ppm"],"Ag_missing":["Ag ppm"],"Au_missing":["Au ppm"],"S_Ag_Au_missing":["S ppm","Ag ppm","Au ppm"]}
ALPHAS=gf.ALPHAS

def predict_masked(samples,ids,cols,model,missing):
    x=[]
    for i in ids:
        row=[]
        for j,c in enumerate(cols):
            v=math.log1p(samples[i]["chem"][c])
            z=(v-model["means"][j])/model["scales"][j]
            row.append(0.0 if c in missing else z)
        x.append(row)
    return [max(0.0,model["ym"]+sum(a*b for a,b in zip(r,model["beta"]))) for r in x]

def metrics(y,p):
    e=[b-a for a,b in zip(y,p)]
    return {"n":len(y),"mae_ppm":statistics.fmean(abs(v) for v in e),"rmse_ppm":math.sqrt(statistics.fmean(v*v for v in e)),"mean_error_ppm":statistics.fmean(e)}

def main():
    samples=fk.load()
    with fk.DATA.open("r",encoding="utf-8-sig",newline="") as f: rows=list(csv.DictReader(f))
    for i,r in enumerate(samples): samples[i]["chem"]={c:float(rows[i][c]) for c in gf.ALL}
    outer=fk.grid_folds(samples);x0=min(s["xyz"][0] for s in samples);y0=min(s["xyz"][1] for s in samples)
    out={name:{"y":[],"p":[],"folds":[]} for name in SCENARIOS}
    mean_y=[];mean_p=[]
    for ofold,(test,train) in enumerate(outer):
        splits=nk.inner_splits(train,samples,x0,y0)
        for name,missing in SCENARIOS.items():
            losses={a:[] for a in ALPHAS}
            for valid,itrain in splits:
                for a in ALPHAS:
                    m=gf.fit(samples,itrain,gf.ALL,a)
                    pp=predict_masked(samples,valid,gf.ALL,m,missing)
                    losses[a].extend(abs(p-samples[i]["cu"]) for p,i in zip(pp,valid))
            cv={a:statistics.fmean(v) for a,v in losses.items()}
            selected=min(cv,key=cv.get)
            model=gf.fit(samples,train,gf.ALL,selected)
            pp=predict_masked(samples,test,gf.ALL,model,missing);yy=[samples[i]["cu"] for i in test]
            out[name]["y"].extend(yy);out[name]["p"].extend(pp)
            out[name]["folds"].append({"outer_fold":ofold+1,"missing_at_target":missing,"selected_alpha":selected,"inner_mae_by_alpha_ppm":cv,"outer_metrics":metrics(yy,pp)})
        mu=statistics.fmean(samples[i]["cu"] for i in train);yy=[samples[i]["cu"] for i in test]
        mean_y.extend(yy);mean_p.extend([mu]*len(yy))
        print(f"completed dropout outer fold {ofold+1}/{len(outer)}",flush=True)
    q75=sorted(out["none"]["y"])[math.ceil(.75*len(out["none"]["y"]))-1]
    results={}
    for name,r in out.items():
        top=[(y,p) for y,p in zip(r["y"],r["p"]) if y>=q75]
        results[name]={"missing_at_target":SCENARIOS[name],"overall":metrics(r["y"],r["p"]),"top_quartile_threshold_ppm":q75,"top_quartile":metrics([y for y,p in top],[p for y,p in top]),"selected_alpha_counts":{str(a):sum(f["selected_alpha"]==a for f in r["folds"]) for a in ALPHAS},"folds":r["folds"]}
    report={"design":"Outer buffered XY folds; ridge fitted on complete training assays, and specified secondary analytes masked to training-fold means at validation locations. Alpha selected within two inner spatial tiles with the same target-location masking.","data_missingness":"Hypothetical stress scenarios; source CSV contains no missing assays.","results":results,"same_split_training_mean":{"overall":metrics(mean_y,mean_p),"top_quartile":metrics([y for y in mean_y if y>=q75],[p for y,p in zip(mean_y,mean_p) if y>=q75])},"caveats":["This tests one analyte or the S/Ag/Au group absent at prediction locations; it does not model assay campaign, detection-limit, or geology-driven missingness.","Training locations retain their complete companion assays.","All target locations are the same 22 spatially held-out tiles.","Top-quartile membership uses observed Cu only for post-hoc evaluation."]}
    path=ROOT/"reports"/"geochem_assay_dropout.json";path.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps({"scenarios":{k:{"overall":v["overall"],"top_quartile":v["top_quartile"],"alpha":v["selected_alpha_counts"]} for k,v in results.items()},"mean":report["same_split_training_mean"],"report":path.name},indent=2))

if __name__=="__main__":main()
