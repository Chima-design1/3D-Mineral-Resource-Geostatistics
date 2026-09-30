"""Nested buffered spatial CV for ridge regression using observed chemistry."""
import json
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import fold_kriging as fk
import nested_kriging_sensitivity as nk

FEATURES = ["Ag ppm","Al ppm","Au ppm","C ppm","Ca ppm","Cl ppm","F ppm","Fe ppm","K ppm","Mg ppm","Mn ppm","Na ppm","P ppm","Pb ppm","S ppm","Th ppm","U ppm"]
ALPHAS = (0.1, 1.0, 10.0, 100.0)


def design(samples, indices, means=None, scales=None):
    raw = [[math.log1p(samples[i]["chem"][c]) for c in FEATURES] for i in indices]
    if means is None:
        means = [statistics.fmean(r[j] for r in raw) for j in range(len(FEATURES))]
        scales = [statistics.pstdev(r[j] for r in raw) or 1.0 for j in range(len(FEATURES))]
    return [[(r[j]-means[j])/scales[j] for j in range(len(FEATURES))] for r in raw], means, scales


def fit(samples, train, alpha):
    x, means, scales = design(samples, train)
    y = [samples[i]["cu"] for i in train]
    ym = statistics.fmean(y)
    yc = [v-ym for v in y]
    p = len(FEATURES)
    xtx = [[sum(row[j]*row[k] for row in x)+(alpha if j==k else 0.0) for k in range(p)] for j in range(p)]
    xty = [sum(row[j]*v for row,v in zip(x,yc)) for j in range(p)]
    beta = fk.solve(xtx, xty)
    return {"means":means,"scales":scales,"y_mean":ym,"beta":beta}


def predict(samples, ids, model):
    x, _, _ = design(samples, ids, model["means"], model["scales"])
    return [max(0.0, model["y_mean"]+sum(a*b for a,b in zip(row,model["beta"]))) for row in x]


def metrics(y,p):
    e=[a-b for a,b in zip(p,y)]
    return {"n":len(y),"mae_ppm":statistics.fmean(abs(v) for v in e),"rmse_ppm":math.sqrt(statistics.fmean(v*v for v in e)),"mean_error_ppm":statistics.fmean(e)}


def main():
    raw=fk.load()
    with fk.DATA.open("r",encoding="utf-8-sig",newline="") as f:
        import csv
        chemistry=list(csv.DictReader(f))
    for i,r in enumerate(raw):
        r["chem"]={c:float(chemistry[i][c]) for c in FEATURES}
    outer=fk.grid_folds(raw)
    x0=min(s["xyz"][0] for s in raw); y0=min(s["xyz"][1] for s in raw)
    ys_all=[]; pred_all=[]; fold_results=[]; outer_rows=[]
    for ofold,(test,train) in enumerate(outer):
        splits=nk.inner_splits(train,raw,x0,y0)
        cv_errors={a:[] for a in ALPHAS}
        for inner_valid, inner_train in splits:
            for alpha in ALPHAS:
                model=fit(raw,inner_train,alpha)
                pp=predict(raw,inner_valid,model)
                cv_errors[alpha].extend(p-raw[i]["cu"] for p,i in zip(pp,inner_valid))
        inner_mae={a:statistics.fmean(abs(e) for e in es) for a,es in cv_errors.items() if es}
        selected=min(inner_mae,key=inner_mae.get)
        model=fit(raw,train,selected)
        pp=predict(raw,test,model); yy=[raw[i]["cu"] for i in test]
        ys_all.extend(yy);pred_all.extend(pp)
        fold_results.append({"outer_fold":ofold+1,"inner_fold_count":len(splits),"selected_alpha":selected,"inner_mae_by_alpha_ppm":inner_mae,"outer_metrics":metrics(yy,pp)})
        for i,y,p in zip(test,yy,pp):
            outer_rows.append({"row":i+1,"outer_fold":ofold+1,"observed_cu_ppm":y,"predicted_cu_ppm":p,"residual_ppm":p-y})
        print(f"completed geochemical outer fold {ofold+1}/{len(outer)} alpha={selected}",flush=True)
    q75=sorted(ys_all)[math.ceil(.75*len(ys_all))-1]
    high=[(y,p) for y,p in zip(ys_all,pred_all) if y>=q75]
    mean_preds=[]
    for test,train in outer:
        mean=statistics.fmean(raw[i]["cu"] for i in train)
        mean_preds.extend([mean]*len(test))
    is_high=[y>=q75 for y in ys_all]
    report={"design":"Outer 300-unit XY tile CV with 50-unit buffer. Ridge alpha selected within two geographically distributed buffered inner tiles. All non-Cu assay concentrations are log1p transformed and standardized using training-only statistics.","target":"Cu ppm","features":FEATURES,"alphas":list(ALPHAS),"outer_metrics":metrics(ys_all,pred_all),"same_split_mean_baseline":metrics(ys_all,mean_preds),"top_quartile_threshold_ppm":q75,"top_quartile_metrics":metrics([y for y,p in high],[p for y,p in high]),"top_quartile_mean_baseline":metrics([y for y,h in zip(ys_all,is_high) if h],[p for p,h in zip(mean_preds,is_high) if h]),"selected_alpha_counts":{str(a):sum(f["selected_alpha"]==a for f in fold_results) for a in ALPHAS},"folds":fold_results,"outer_predictions":outer_rows,"caveats":["Validation assumes companion-element assays are known at every held-out location while Cu is the unavailable target; this is a covariate-assisted prediction experiment, not conventional estimation into unassayed blocks.","All chemistry columns are exploratory covariates; no geological interpretations are assigned.","Nested validation uses two inner tiles per outer fold, a limited tuning sample.","This test is predictive and does not establish cause, ore domains, or economic/resource validity."]}
    out=ROOT/"reports"/"geochem_covariate_cv.json"
    out.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps({"outer":report["outer_metrics"],"mean":report["same_split_mean_baseline"],"top_quartile":report["top_quartile_metrics"],"alpha_counts":report["selected_alpha_counts"],"report":out.name},indent=2))


if __name__=="__main__": main()
