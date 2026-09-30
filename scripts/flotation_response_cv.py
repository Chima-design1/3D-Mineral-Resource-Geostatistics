"""Nested leave-one-HOLEID-out prediction of flotation LCT from chemistry."""
import csv,json,math,statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"/"raw"/"flotation.csv"
OUT=ROOT/"reports"/"flotation_response_cv.json"
FEATURES=["Ag ppm","Al ppm","Au ppm","C ppm","Ca ppm","Cl ppm","Cu ppm","F ppm","Fe ppm","K ppm","Mg ppm","Mn ppm","Na ppm","P ppm","S ppm","Th ppm","Ti ppm","U ppm"]
ALPHAS=(0.1,1.0,10.0,100.0,1000.0)

def make_folds(records, groups, n_splits):
    buckets=[[] for _ in range(n_splits)]
    sizes=[0]*n_splits
    group_sizes=defaultdict(int)
    for i in groups: group_sizes[records[i]["HOLEID"]]+=1
    ordered=sorted(group_sizes,key=lambda h:(-group_sizes[h],str(h)))
    assignment={}
    for h in ordered:
        k=min(range(n_splits),key=lambda j:sizes[j])
        assignment[h]=k;sizes[k]+=group_sizes[h]
    for i in groups: buckets[assignment[records[i]["HOLEID"]]].append(i)
    return [(b,[i for i in groups if i not in set(b)]) for b in buckets if b]

def matrices(rows, ids, means=None, scales=None):
    raw=[]
    for i in ids:
        vals=[]
        for c in FEATURES:
            s=rows[i][c].strip()
            vals.append(math.log1p(float(s)) if s else None)
        raw.append(vals)
    if means is None:
        means=[];scales=[]
        for j in range(len(FEATURES)):
            valid=[r[j] for r in raw if r[j] is not None]
            m=statistics.fmean(valid) if valid else 0.0
            means.append(m);scales.append(statistics.pstdev(valid) or 1.0)
    x=[[0.0 if v is None else (v-means[j])/scales[j] for j,v in enumerate(r)] for r in raw]
    return x,means,scales

def solve(a,b):
    n=len(b);m=[a[i][:]+[b[i]] for i in range(n)]
    for c in range(n):
        p=max(range(c,n),key=lambda i:abs(m[i][c]))
        if abs(m[p][c])<1e-12: raise ArithmeticError("singular ridge system")
        m[c],m[p]=m[p],m[c];d=m[c][c]
        for j in range(c,n+1):m[c][j]/=d
        for i in range(n):
            if i==c:continue
            f=m[i][c]
            for j in range(c,n+1):m[i][j]-=f*m[c][j]
    return [m[i][n] for i in range(n)]

def fit(rows,train,alpha):
    x,means,scales=matrices(rows,train)
    y=[float(rows[i]["LCT"]) for i in train];ym=statistics.fmean(y);yc=[v-ym for v in y];p=len(FEATURES)
    xtx=[[sum(r[j]*r[k] for r in x)+(alpha if j==k else 0) for k in range(p)] for j in range(p)]
    xty=[sum(r[j]*v for r,v in zip(x,yc)) for j in range(p)]
    return {"means":means,"scales":scales,"ym":ym,"beta":solve(xtx,xty)}

def predict(rows,ids,model):
    x,_,_=matrices(rows,ids,model["means"],model["scales"])
    return [model["ym"]+sum(a*b for a,b in zip(r,model["beta"])) for r in x]

def metrics(y,p):
    e=[a-b for a,b in zip(p,y)]
    return {"n":len(y),"mae":statistics.fmean(abs(x) for x in e),"rmse":math.sqrt(statistics.fmean(x*x for x in e)),"mean_error":statistics.fmean(e)}

def main():
    with DATA.open("r",encoding="utf-8-sig",newline="") as f: allrows=list(csv.DictReader(f))
    rows=[r for r in allrows if r["LCT"].strip()]
    groups=defaultdict(list)
    for i,r in enumerate(rows):groups[r["HOLEID"]].append(i)
    folds=[];actual=[];pred=[];meanpred=[]
    for n,(hole,test) in enumerate(sorted(groups.items(),key=lambda x:str(x[0])),1):
        train=[i for i in range(len(rows)) if i not in set(test)]
        inner=make_folds(rows,train,min(4,len(set(rows[i]["HOLEID"] for i in train))))
        scores={a:[] for a in ALPHAS}
        for valid,itrain in inner:
            if not valid or not itrain:continue
            for alpha in ALPHAS:
                model=fit(rows,itrain,alpha)
                pv=predict(rows,valid,model)
                scores[alpha].extend(abs(p-float(rows[i]["LCT"])) for p,i in zip(pv,valid))
        cv={a:statistics.fmean(v) for a,v in scores.items() if v}
        alpha=min(cv,key=cv.get)
        model=fit(rows,train,alpha);p=predict(rows,test,model);y=[float(rows[i]["LCT"]) for i in test]
        baseline=statistics.fmean(float(rows[i]["LCT"]) for i in train)
        actual.extend(y);pred.extend(p);meanpred.extend([baseline]*len(test))
        folds.append({"held_out_HOLEID":hole,"n_test":len(test),"n_train":len(train),"selected_alpha":alpha,"inner_mae_by_alpha":cv,"ridge_metrics":metrics(y,p),"training_mean_metrics":metrics(y,[baseline]*len(test)),"predictions":[{"row":i+2,"HOLEID":hole,"observed_LCT":float(rows[i]["LCT"]),"predicted_LCT":v,"residual":v-float(rows[i]["LCT"])} for i,v in zip(test,p)]})
    report={"target":"LCT (flotation table field)","n_records_all":len(allrows),"n_records_with_target":len(rows),"distinct_HOLEID_with_target":len(groups),"features":FEATURES,"validation":"Nested leave-one-HOLEID-out outer CV; inner folds group by HOLEID. Ridge alpha selected by inner pooled MAE.","outer_metrics":{"ridge":metrics(actual,pred),"training_mean":metrics(actual,meanpred)},"selected_alpha_counts":{str(a):sum(f["selected_alpha"]==a for f in folds) for a in ALPHAS},"folds":folds,"caveats":["fr and xr were excluded because they are flotation test intermediates closely related to the LCT response; only chemistry was used as predictors.","Predictor assays are available on each held-out flotation test record; this is conditional laboratory-response prediction, not geospatial interpolation.","Only 52 target records across a small set of holes are available; scores are exploratory and may be unstable.","No test responses were joined to the 2,000 drillhole table because a direct sample crosswalk and support definition are absent."]}
    OUT.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps({"target_records":len(rows),"holes":len(groups),"metrics":report["outer_metrics"],"alpha_counts":report["selected_alpha_counts"],"report":OUT.name},indent=2))

if __name__=="__main__":main()
