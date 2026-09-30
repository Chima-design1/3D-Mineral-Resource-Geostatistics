"""Screen high Cu observations and large consecutive within-hole gaps."""
import csv
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "drillholes.csv"
OUT = ROOT / "reports" / "sample_influence.json"


def main():
    with DATA.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    vals = [float(r["Cu ppm"]) for r in rows]
    ordered = sorted(vals)
    q99 = ordered[math.ceil(.99 * len(ordered)) - 1]
    q1, q3 = statistics.quantiles(vals, n=4, method="inclusive")[:1][0], statistics.quantiles(vals, n=4, method="inclusive")[2]
    fence = q3 + 1.5 * (q3 - q1)
    top = sorted(enumerate(rows, 1), key=lambda ir: float(ir[1]["Cu ppm"]), reverse=True)[:20]
    by_hole = {}
    for rowno, r in enumerate(rows, 1):
        by_hole.setdefault(r["HOLEID"], []).append((rowno, tuple(float(r[c]) for c in ("X", "Y", "Z"))))
    gaps = []
    for hole, pts in by_hole.items():
        for (r0, p0), (r1, p1) in zip(pts, pts[1:]):
            d = math.dist(p0, p1)
            if d > 100:
                gaps.append({"hole": hole, "row_before": r0, "row_after": r1, "distance": d, "xyz_before": p0, "xyz_after": p1})
    report = {
        "target": "Cu ppm", "n": len(rows), "zero_count": sum(v == 0 for v in vals),
        "q99_ppm": q99, "upper_iqr_fence_ppm": fence, "above_fence_count": sum(v > fence for v in vals),
        "top_20": [{"row": i, "hole": r["HOLEID"], "cu_ppm": float(r["Cu ppm"]), "xyz": [float(r[c]) for c in ("X", "Y", "Z")]} for i, r in top],
        "consecutive_row_gaps_over_100_units": sorted(gaps, key=lambda g: g["distance"], reverse=True),
        "interpretation": "Screening only: high values and large row-order gaps are not presumed erroneous. Source does not document assay detection limits or establish row-order interval support. Verify against primary sample records before exclusion or compositing."
    }
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"zeros": report["zero_count"], "q99": q99, "upper_fence": fence, "above_fence": report["above_fence_count"], "gaps_over_100": len(gaps), "top_cu": report["top_20"][0], "report": OUT.name}, indent=2))


if __name__ == "__main__":
    main()
