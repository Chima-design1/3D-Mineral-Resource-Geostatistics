"""Compare predictions under explicit co-located assay availability scenarios."""
import csv, json, math, statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"/"raw"/"drillholes.csv"

def metrics(rows):
    e=[r["predicted_cu_ppm"]-r["observed_cu_ppm"] for r in rows]
    return {"n":len(rows),"mae_ppm":statistics.fmean(abs(x) for x in e),"rmse_ppm":math.sqrt(statistics.fmean(x*x for x in e)),"mean_error_ppm":statistics.fmean(e)}

def main():
    with DATA.open("r",encoding="utf-8-sig",newline="") as f: data=list(csv.DictReader(f))
    missing={c:sum(not r[c].strip() for r in data) for c in data[0]}
    feature=json.loads((ROOT/"reports"/"geochem_feature_sensitivity.json").read_text(encoding="utf-8"))
    spatial=json.loads((ROOT/"reports"/"nested_kriging_sensitivity.json").read_text(encoding="utf-8"))
    sets=feature["models"]
    scenarios={
        "Cu_missing_other_17_assays_available":{"prediction_source":"nested all-element ridge","availability":"All non-Cu chemical assays available at each target location; Cu unavailable.","rows":sets["all_17"]["outer_predictions"]},
        "Cu_missing_S_Ag_Au_available":{"prediction_source":"nested S/Ag/Au ridge","availability":"Only S, Ag, and Au available at each target location; Cu unavailable.","rows":sets["S_Ag_Au"]["outer_predictions"]},
        "all_assays_missing_at_target_location":{"prediction_source":"nested raw-Cu ordinary kriging","availability":"No target-location assays available; neighboring training locations retain Cu measurements.","rows":[{"row":r["row"],"outer_fold":r["outer_fold"],"observed_cu_ppm":r["observed_cu_ppm"],"predicted_cu_ppm":r["predicted_cu_ppm"]} for r in spatial["outer_predictions"]]}
    }
    all_rows=scenarios["Cu_missing_other_17_assays_available"]["rows"]
    q75=sorted(r["observed_cu_ppm"] for r in all_rows)[math.ceil(.75*len(all_rows))-1]
    for s in scenarios.values():
        s["overall"]=metrics(s.pop("rows"))
        # Reconstruct rows to compute consistent observed-grade subgroup.
        key=next(k for k,v in scenarios.items() if v is s)
        if key.startswith("Cu_missing_other"):
            source=sets["all_17"]["outer_predictions"]
        elif key.startswith("Cu_missing_S"):
            source=sets["S_Ag_Au"]["outer_predictions"]
        else:
            source=spatial["outer_predictions"]
        top=[r for r in source if r["observed_cu_ppm"]>=q75]
        s["top_quartile_threshold_ppm"]=q75
        s["top_quartile"]=metrics(top)
        s["source_rows_match_outer_validation"]=len(source)==2000 and len({r["row"] for r in source})==2000
    report={"dataset_assay_missing_values_by_field":missing,"source_file_complete_case":all(v==0 for v in missing.values()),"availability_scenarios":scenarios,"design_notes":["The source has no naturally missing chemical observations, so the scenarios are availability assumptions, not measured missingness patterns.","All three predictions are out-of-fold on the same 22 buffered XY tiles.","The all-assays-missing scenario means no co-located target assays at the prediction site, while nearby training sites still have Cu values.","Covariate models assume the listed auxiliary assays are genuinely available at every held-out target site; they cannot estimate into locations with no assays.","The S/Ag/Au subset was selected after descriptive inspection of full-data correlations, so that model comparison is exploratory."]}
    out=ROOT/"reports"/"assay_availability_scenarios.json";out.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps({"source_complete_case":report["source_file_complete_case"],"scenarios":{k:{"overall":v["overall"],"top_quartile":v["top_quartile"],"row_check":v["source_rows_match_outer_validation"]} for k,v in scenarios.items()},"report":out.name},indent=2))

if __name__=="__main__": main()
