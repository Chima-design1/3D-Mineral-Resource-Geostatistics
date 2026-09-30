"""Describe high-Cu samples using only fields observed in drillholes.csv."""
import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "drillholes.csv"
OUT = ROOT / "reports" / "high_grade_context.json"
ANALYTES = ["Ag ppm", "Al ppm", "Au ppm", "C ppm", "Ca ppm", "Cl ppm", "Cu ppm", "F ppm", "Fe ppm", "K ppm", "Mg ppm", "Mn ppm", "Na ppm", "P ppm", "Pb ppm", "S ppm", "Th ppm", "U ppm"]


def pearson(a, b):
    ma, mb = statistics.fmean(a), statistics.fmean(b)
    da = [x-ma for x in a]; db = [x-mb for x in b]
    den = math.sqrt(sum(x*x for x in da)*sum(x*x for x in db))
    return sum(x*y for x,y in zip(da,db))/den if den else None


def main():
    with DATA.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    q = {c: [float(r[c]) for r in rows] for c in ANALYTES}
    cu = q["Cu ppm"]
    ordered = sorted(cu)
    q90, q75 = ordered[math.ceil(.9*len(cu))-1], ordered[math.ceil(.75*len(cu))-1]
    high = [i for i,v in enumerate(cu) if v >= q90]
    coords = [tuple(float(rows[i][c]) for c in ("X","Y","Z")) for i in range(len(rows))]
    x0, y0 = min(x for x,_,_ in coords), min(y for _,y,_ in coords)
    tile_counts = Counter((math.floor((coords[i][0]-x0)/300), math.floor((coords[i][1]-y0)/300)) for i in high)
    hole_counts = Counter(rows[i]["HOLEID"] for i in high)
    correlations = []
    logcu = [math.log1p(v) for v in cu]
    for c in ANALYTES:
        if c == "Cu ppm": continue
        r = pearson(logcu, [math.log1p(v) for v in q[c]])
        correlations.append({"analyte":c,"pearson_log1p":r,"median_in_top_decile_ppm":statistics.median(q[c][i] for i in high),"median_all_ppm":statistics.median(q[c])})
    correlations.sort(key=lambda r: abs(r["pearson_log1p"] or 0), reverse=True)
    # Count high-grade pairs close in supplied 3D coordinates; descriptive only.
    near_pairs_50 = near_pairs_100 = 0
    for ii, i in enumerate(high):
        for j in high[ii+1:]:
            d = math.dist(coords[i], coords[j])
            near_pairs_50 += d <= 50
            near_pairs_100 += d <= 100
    report = {
        "source_file":"drillholes.csv", "n":len(rows), "target":"Cu ppm", "top_decile_cutoff_ppm":q90, "top_decile_n":len(high),
        "top_decile_hole_counts":hole_counts.most_common(20),
        "top_decile_300_unit_xy_tile_counts":[{"tile":list(k),"n":v} for k,v in tile_counts.most_common()],
        "top_decile_close_pairs_xyz_units":{"within_50":near_pairs_50,"within_100":near_pairs_100,"interpretation":"Pair counts among top-decile samples; not independent cluster counts."},
        "log1p_pearson_correlations_with_cu":correlations,
        "quartile_high_support_summary":{"threshold_ppm":q75,"n":sum(v>=q75 for v in cu),"median_cu_ppm":statistics.median(v for v in cu if v>=q75),"max_cu_ppm":max(v for v in cu if v>=q75)},
        "limitations":["Only chemical concentrations, HOLEID, and XYZ sample-centroid coordinates are observed in drillholes.csv.","No lithology, alteration, mineralization, assay quality, interval from/to, sample length, recovery, density, collar/survey, or domain field is present.","Chemical associations and spatial concentrations are exploratory covariation, not geological domains or causal explanations.","Pearson correlations use log1p concentrations; compositional closure and censoring/detection limits are not characterized."]
    }
    OUT.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps({"q90_cutoff":q90,"top_decile_n":len(high),"top_holes":hole_counts.most_common(8),"top_tiles":tile_counts.most_common(8),"near_pairs":report["top_decile_close_pairs_xyz_units"],"top_correlations":correlations[:8],"report":OUT.name},indent=2))


if __name__ == "__main__": main()
