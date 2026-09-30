"""Check whether the supplied assay columns can be treated as a composition."""
import csv
import json
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"/"raw"/"drillholes.csv"
OUT=ROOT/"reports"/"composition_readiness.json"
FOUR=("Cu ppm","Au ppm","Ag ppm","S ppm")


def summarize(rows, columns):
    totals=[sum(float(r[c]) for c in columns) for r in rows]
    return {"components":list(columns),"n_components":len(columns),"row_sum_ppm":{"min":min(totals),"q1":statistics.quantiles(totals,n=4)[0],"median":statistics.median(totals),"mean":statistics.fmean(totals),"q3":statistics.quantiles(totals,n=4)[2],"max":max(totals)},"rows_with_sum_between_900000_and_1100000_ppm":sum(900000<=x<=1100000 for x in totals),"zeros_by_component":{c:sum(float(r[c])==0 for r in rows) for c in columns},"negative_by_component":{c:sum(float(r[c])<0 for r in rows) for c in columns}}


def main():
    with DATA.open(encoding="utf-8-sig",newline="") as f: rows=list(csv.DictReader(f))
    assays=[c for c in rows[0] if c.endswith(" ppm")]
    out={"dataset":"GeoMet drillholes.csv","rows":len(rows),"assay_columns":len(assays),"all_18_assays":summarize(rows,assays),"selected_Cu_Au_Ag_S_subcomposition":summarize(rows,FOUR),"zero_rows_in_selected_subcomposition":sum(any(float(r[c])==0 for c in FOUR) for r in rows),"source_metadata":{"assay_sum_closure_definition":None,"below_detection_or_zero_semantics":None,"full_major_oxide_or_elemental_composition":False},"decision":"Do not apply CLR as though the 18 reported analytes form a closed whole-rock composition. A Cu-Au-Ag-S CLR could only be described as a four-part subcomposition; it contains zeros whose analytical meaning and replacement policy are undocumented, so any transformed model requires an explicit sensitivity analysis and must not be interpreted as whole-rock mass fractions."}
    OUT.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__":main()
