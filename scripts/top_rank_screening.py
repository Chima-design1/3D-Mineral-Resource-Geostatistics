"""Evaluate capture of high-Cu samples among top-ranked buffered holdouts."""
import json
import math
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REPORTS=ROOT/"reports"
OUT=REPORTS/"top_rank_screening.json"
Q75,Q90=8900.0,18100.0
FRACTIONS=(0.05,0.10,0.20)


def load_prediction_map(filename,key):
    data=json.loads((REPORTS/filename).read_text(encoding="utf-8"));return {int(r["row"]):float(r[key]) for r in data["outer_predictions"]}


def main():
    krig=load_prediction_map("nested_kriging_sensitivity.json","predicted_cu_ppm")
    chem=load_prediction_map("geochem_covariate_cv.json","predicted_cu_ppm")
    feature=json.loads((REPORTS/"geochem_feature_sensitivity.json").read_text(encoding="utf-8"))["models"]["all_17"]["outer_predictions"]
    chem_high={int(r["row"]):float(r["predicted_cu_ppm"]) for r in feature}
    sub=json.loads((REPORTS/"subcomposition_spatial_sensitivity.json").read_text(encoding="utf-8"))["designs"]["buffered_XY_tile"]["predictions_by_row"]
    submap={int(r["row"]):float(r["clr_idw_total_idw_zero_0.1x_min"]) for r in sub}
    idw=json.loads((REPORTS/"grid_support_validation.json").read_text(encoding="utf-8"))["records"]
    maps={"raw_Cu_IDW_p2_k12":{int(r["row"]):float(r["idw"]) for r in idw},"nested_ordinary_kriging":krig,"conditional_all_assay_ridge_alpha_0.1_to_100":chem,"conditional_all_assay_ridge_alpha_100_to_10000":chem_high,"Cu_Au_Ag_S_subcomposition_IDW":submap}
    actual={int(r["row"]):float(r["actual"]) for r in idw}; rows=sorted(actual);n=len(rows)
    fold_groups={}
    for r in idw:fold_groups.setdefault(tuple(r["fold"]),[]).append(int(r["row"]))
    prevalence={"top_quartile":sum(actual[i]>=Q75 for i in rows)/n,"top_decile":sum(actual[i]>=Q90 for i in rows)/n}
    total_cu=sum(actual.values());results={}
    for name,scores in maps.items():
        ranked=sorted(rows,key=lambda i:(-scores[i],i)); summaries=[]
        for fraction in FRACTIONS:
            take=max(1,math.ceil(n*fraction)); selected=ranked[:take]
            selected_mean=sum(actual[i] for i in selected)/take
            item={"review_fraction":fraction,"selected_n":take,"mean_observed_cu_ppm":selected_mean,"share_of_total_observed_Cu_in_selected":sum(actual[i] for i in selected)/total_cu}
            for label,threshold,prev in (("top_quartile",Q75,prevalence["top_quartile"]),("top_decile",Q90,prevalence["top_decile"])):
                hits=sum(actual[i]>=threshold for i in selected)
                item[label]={"threshold_ppm":threshold,"hits":hits,"precision":hits/take,"recall":hits/sum(y>=threshold for y in actual.values()),"lift_over_random":(hits/take)/prev if prev else None}
            summaries.append(item)
        results[name]=summaries
    within_fold={}
    for name,scores in maps.items():
        bands=[]
        for fraction in FRACTIONS:
            selected_all=[];per_fold=[]
            for fold_id,fold_rows in sorted(fold_groups.items()):
                ordered=sorted(fold_rows,key=lambda i:(-scores[i],i));take=max(1,math.ceil(len(ordered)*fraction));selected=ordered[:take];selected_all.extend(selected)
                stats={"fold":list(fold_id),"n_holdout":len(fold_rows),"selected_n":take}
                for label,threshold in (("top_quartile",Q75),("top_decile",Q90)):
                    hits=sum(actual[i]>=threshold for i in selected);positives=sum(actual[i]>=threshold for i in fold_rows)
                    stats[label]={"hits":hits,"positive_n":positives,"precision":hits/take,"recall":hits/positives if positives else None}
                per_fold.append(stats)
            item={"review_fraction_per_fold":fraction,"selected_total":len(selected_all),"folds":len(per_fold)}
            for label,threshold,prev in (("top_quartile",Q75,prevalence["top_quartile"]),("top_decile",Q90,prevalence["top_decile"])):
                hits=sum(actual[i]>=threshold for i in selected_all);positives=sum(y>=threshold for y in actual.values());prec=hits/len(selected_all)
                precision_folds=[r[label]["precision"] for r in per_fold]
                recall_folds=[r[label]["recall"] for r in per_fold if r[label]["recall"] is not None]
                item[label]={"threshold_ppm":threshold,"hits":hits,"actual_positive_n":positives,"precision":prec,"recall":hits/positives,"lift_over_random":prec/prev,"fold_precision_q1_median_q3":[statistics.quantiles(precision_folds,n=4)[0],statistics.median(precision_folds),statistics.quantiles(precision_folds,n=4)[2]],"fold_recall_q1_median_q3_positive_folds":[statistics.quantiles(recall_folds,n=4)[0],statistics.median(recall_folds),statistics.quantiles(recall_folds,n=4)[2]],"folds_with_at_least_one_actual_positive":len(recall_folds)}
            bands.append(item)
        within_fold[name]=bands
    random_baseline={"top_quartile_precision":prevalence["top_quartile"],"top_decile_precision":prevalence["top_decile"],"recall_at_fraction":{str(f):f for f in FRACTIONS},"mean_observed_cu_ppm":sum(actual.values())/n,"share_of_total_Cu_at_fraction":{str(f):f for f in FRACTIONS}}
    out={"validation":"Pooled predictions from the same buffered 300-unit XY tiles with a 50-unit horizontal training buffer.","n":n,"fold_count":len(fold_groups),"thresholds":{"top_quartile_descriptive_ppm":Q75,"top_decile_descriptive_ppm":Q90},"observed_prevalence":prevalence,"review_fractions":list(FRACTIONS),"random_selection_reference":random_baseline,"pooled_global_score_ranking":results,"within_fold_percentile_ranking":within_fold,"caveats":["Selection fractions are workload scenarios, not economic cutoffs.","Within-fold percentile ranking avoids comparing absolute score scales across folds; pooled global score ranking may be affected by fold calibration differences.","All metrics remain descriptive and do not establish deposit-scale prospectivity.","The conditional all-assay ridge ranking assumes 17 non-Cu assays are available at held-out locations.","Spatial-only methods use limited/experimental interpolation and do not define geological domains."]}
    OUT.write_text(json.dumps(out,indent=2),encoding="utf-8")
    for name,vals in within_fold.items():
        print(name)
        for r in vals:print(r["review_fraction_per_fold"],"Q75 precision/recall",round(r["top_quartile"]["precision"],3),round(r["top_quartile"]["recall"],3),"Q90",round(r["top_decile"]["precision"],3),round(r["top_decile"]["recall"],3))


if __name__=="__main__":main()
