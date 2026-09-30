"""Create a strictly support-limited exploratory 3D Cu IDW grid."""
import csv
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "drillholes.csv"
SUPPORT = ROOT / "reports" / "spatial_support_audit.json"
OUT = ROOT / "outputs" / "copper_idw_grid.csv"
META = ROOT / "reports" / "copper_idw_grid.json"
STEP = 50.0
POWER = 2.0
K = 12


def hull(points):
    pts = sorted(set(points))
    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def inside(point, polygon):
    x, y = point
    hit = False
    j = len(polygon)-1
    for i, (xi, yi) in enumerate(polygon):
        xj, yj = polygon[j]
        if (yi > y) != (yj > y) and x < (xj-xi)*(y-yi)/(yj-yi)+xi:
            hit = not hit
        j = i
    return hit


def main():
    with DATA.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    pts = [(tuple(float(r[c]) for c in ("X","Y","Z")), float(r["Cu ppm"])) for r in rows]
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    max_nearest = support["nearest_sample_spacing_units"]["3d"]["q3"]
    max_kth = 2 * support["nearest_sample_spacing_units"]["nearest_sample_from_different_hole_3d"]["median"]
    polygon = hull([(p[0][0],p[0][1]) for p in pts])
    mins = [min(p[0][i] for p in pts) for i in range(3)]
    maxs = [max(p[0][i] for p in pts) for i in range(3)]
    axes = [[mins[i] + j*STEP for j in range(math.ceil((maxs[i]-mins[i])/STEP)+1)] for i in range(3)]
    candidate_count = accepted = 0
    nearest_accepted, kth_accepted = [], []
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=["X","Y","Z","Cu_IDW_ppm","nearest_sample_distance","12th_neighbor_distance","support_status"])
        writer.writeheader()
        for x in axes[0]:
            for y in axes[1]:
                if not inside((x,y),polygon):
                    continue
                for z in axes[2]:
                    candidate_count += 1
                    nearest=sorted((math.dist((x,y,z),xyz),cu) for xyz,cu in pts)[:K]
                    d1=nearest[0][0]; dk=nearest[-1][0]
                    if d1 > max_nearest or dk > max_kth:
                        continue
                    if d1 == 0:
                        estimate=nearest[0][1]
                    else:
                        weights=[1/(d**POWER) for d,_ in nearest]
                        estimate=sum(w*v for w,(_,v) in zip(weights,nearest))/sum(weights)
                    writer.writerow({"X":f"{x:.3f}","Y":f"{y:.3f}","Z":f"{z:.3f}","Cu_IDW_ppm":f"{estimate:.3f}","nearest_sample_distance":f"{d1:.3f}","12th_neighbor_distance":f"{dk:.3f}","support_status":"within provisional thresholds"})
                    accepted += 1; nearest_accepted.append(d1); kth_accepted.append(dk)
    meta={"purpose":"Exploratory Cu concentration interpolation, not a mineral resource estimate.","method":{"idw_power":POWER,"neighbors":K,"grid_spacing_coordinate_units":STEP,"domain":"XY convex hull of supplied sample centroids, with Z limited to observed min/max","support_limits":{"nearest_3d_distance_max":max_nearest,"basis":"Q3 of observed nearest-sample 3D spacing","12th_neighbor_3d_distance_max":max_kth,"basis_2":"2x median observed distance to nearest sample from a different hole"}},"coordinates":{"units":"Undocumented source-local units","CRS":"Not provided"},"grid_counts":{"candidate_nodes_inside_xy_hull_and_z_range":candidate_count,"nodes_passing_support_limits":accepted,"nodes_rejected_for_weak_support":candidate_count-accepted,"support_fraction":accepted/candidate_count if candidate_count else None},"accepted_node_support":{"nearest_sample_distance_median":statistics.median(nearest_accepted) if accepted else None,"12th_neighbor_distance_median":statistics.median(kth_accepted) if accepted else None},"output_csv":str(OUT.relative_to(ROOT)),"limitations":["Convex hull is geometric and may span unsampled geological gaps.","IDW is a simple baseline and does not use geological domains or directional anisotropy.","Support thresholds are provisional screening rules based on this dataset's spacing, not regulatory or resource-classification criteria.","No block support, compositing, density, topography, economics, or mineralized-domain boundary is available.","Predictions outside the accepted nodes are omitted, not extrapolated."]}
    META.write_text(json.dumps(meta,indent=2),encoding="utf-8")
    print(json.dumps(meta,indent=2))


if __name__ == "__main__":
    main()
