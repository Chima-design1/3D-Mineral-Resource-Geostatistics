"""Descriptive paired bootstrap over buffered spatial CV folds."""
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REPORTS=ROOT/"reports"
OUT=REPORTS/"spatial_fold_bootstrap.json"
N_BOOT=10000
SEED=20260930
HIGH_Q75,HIGH_Q90,REVIEW=8900.0,18100.0,0.20


def prediction_map(filename,key):
    data=json.loads((REPORTS/filename).read_text(encoding="utf-8"))
    return {int(r["row"]):float(r[key]) for r in data["outer_predictions"]}


def interval(values):
    a=sorted(values)
    return {"q025":a[int(.025*(len(a)-1))],"median":statistics.median(a),"q975":a[int(.975*(len(a)-1))]}


def main():
    base=json.loads((REPORTS/"grid_support_validation.json").read_text(encoding="utf-8"))["records"]
    krig=prediction_map("nested_kriging_sensitivity.json","predicted_cu_ppm")
    chem=prediction_map("geochem_covariate_cv.json","predicted_cu_ppm")
    feature=json.loads((REPORTS/"geochem_feature_sensitivity.json").read_text(encoding="utf-8"))["models"]["all_17"]["outer_predictions"]
    chem_high={int(r["row"]):float(r["predicted_cu_ppm"]) for r in feature}
    subrows=json.loads((REPORTS/"subcomposition_spatial_sensitivity.json").read_text(encoding="utf-8"))["designs"]["buffered_XY_tile"]["predictions_by_row"]
    sub={int(r["row"]):float(r["clr_idw_total_idw_zero_0.1x_min"]) for r in subrows}
    predictions={"raw_Cu_IDW_p2_k12":{int(r["row"]):float(r["idw"]) for r in base},"nested_ordinary_kriging":krig,"conditional_all_assay_ridge_alpha_0.1_to_100":chem,"conditional_all_assay_ridge_alpha_100_to_10000":chem_high,"Cu_Au_Ag_S_subcomposition_IDW":sub}
    groups=defaultdict(list)
    for r in base:groups[tuple(r["fold"])].append(int(r["row"]))
    folds=[]
    for fold_id,rows in sorted(groups.items()):
        actual={int(r["row"]):float(r["actual"]) for r in base if tuple(r["fold"])==fold_id}
        mean={int(r["row"]):float(r["mean"]) for r in base if tuple(r["fold"])==fold_id}
        aggregate={}
        for name,pmap in predictions.items():
            abs_err=[abs(pmap[i]-actual[i]) for i in rows];sq_err=[(pmap[i]-actual[i])**2 for i in rows]
            high=[i for i in rows if actual[i]>=HIGH_Q75]
            rank=sorted(rows,key=lambda i:(-pmap[i],i));selected=rank[:max(1,math.ceil(len(rank)*REVIEW))]
            aggregate[name]={"n":len(rows),"abs_error_sum":sum(abs_err),"sq_error_sum":sum(sq_err),"high_n":len(high),"high_abs_error_sum":sum(abs(pmap[i]-actual[i]) for i in high),"top_q90_hits":sum(actual[i]>=HIGH_Q90 for i in selected),"selected_n":len(selected),"q90_positive_n":sum(actual[i]>=HIGH_Q90 for i in rows)}
        aggregate["training_mean"]={"n":len(rows),"abs_error_sum":sum(abs(mean[i]-actual[i]) for i in rows),"sq_error_sum":sum((mean[i]-actual[i])**2 for i in rows),"high_n":len([i for i in rows if actual[i]>=HIGH_Q75]),"high_abs_error_sum":sum(abs(mean[i]-actual[i]) for i in rows if actual[i]>=HIGH_Q75)}
        folds.append({"id":list(fold_id),"stats":aggregate})
    rng=random.Random(SEED);names=list(predictions);diffs={name:{"mae":[],"rmse":[],"high_grade_mae":[],"top_q90_recall":[]} for name in names}
    per_method={name:{"mae":[],"rmse":[],"top_q90_recall":[],"top_q90_precision":[]} for name in names}
    for _ in range(N_BOOT):
        draw=[rng.randrange(len(folds)) for _ in folds]
        summed={name:defaultdict(float) for name in [*names,"training_mean"]}
        for ix in draw:
            for name,stats in folds[ix]["stats"].items():
                for key,value in stats.items():summed[name][key]+=value
        b=summed["training_mean"];base_mae=b["abs_error_sum"]/b["n"];base_rmse=math.sqrt(b["sq_error_sum"]/b["n"]);base_high=b["high_abs_error_sum"]/b["high_n"]
        for name in names:
            a=summed[name];mae=a["abs_error_sum"]/a["n"];rmse=math.sqrt(a["sq_error_sum"]/a["n"]);high=a["high_abs_error_sum"]/a["high_n"]
            recall=a["top_q90_hits"]/a["q90_positive_n"] if a["q90_positive_n"] else None
            precision=a["top_q90_hits"]/a["selected_n"]
            per_method[name]["mae"].append(mae);per_method[name]["rmse"].append(rmse)
            per_method[name]["top_q90_recall"].append(recall);per_method[name]["top_q90_precision"].append(precision)
            diffs[name]["mae"].append(mae-base_mae);diffs[name]["rmse"].append(rmse-base_rmse);diffs[name]["high_grade_mae"].append(high-base_high)
            diffs[name]["top_q90_recall"].append(recall)
    result={}
    for name in names:
        result[name]={"model_metric_intervals":{metric:interval(values) for metric,values in per_method[name].items()},"paired_difference_vs_training_mean":{metric:{**interval(values),"bootstrap_probability_difference_below_zero":sum(v<0 for v in values)/len(values)} for metric,values in diffs[name].items() if metric!="top_q90_recall"},"top20_percent_within_fold_q90_screen":{"recall_interval":interval(per_method[name]["top_q90_recall"]),"precision_interval":interval(per_method[name]["top_q90_precision"])}}
    out={"purpose":"Assess spatial-fold stability with a paired cluster bootstrap over the 22 buffered XY holdout tiles.","folds":len(folds),"bootstrap_replicates":N_BOOT,"seed":SEED,"comparison":"All model errors are compared with the sample-weighted training mean on identical held-out rows; high-grade errors use observed Cu >= 8,900 ppm; top-rank screen uses top 20% within each fold and observed Cu >= 18,100 ppm.","results":result,"limitations":["Only 22 spatial folds are available; percentile intervals are descriptive stability intervals, not formal population confidence intervals.","Spatial tiles may not be independent and some folds contain few high-grade observations.","Conditional all-assay ridge assumes non-Cu assays are available at held-out locations.","Top-rank review fractions are operational examples, not economic cutoffs."]}
    OUT.write_text(json.dumps(out,indent=2),encoding="utf-8")
    for name,res in result.items():print(name,"MAE delta vs mean",res["paired_difference_vs_training_mean"]["mae"],"top20% Q90 recall",res["top20_percent_within_fold_q90_screen"]["recall_interval"])


if __name__=="__main__":main()
