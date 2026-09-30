"""Assess how provisional distance cutoffs change the candidate grid footprint."""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "drillholes.csv"
SUPPORT = ROOT / "reports" / "spatial_support_audit.json"
OUT = ROOT / "reports" / "grid_support_sensitivity.json"
STEP = 50.0


def hull(points):
    pts = sorted(set(points))
    def cross(o, a, b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lo=[]
    for p in pts:
        while len(lo)>=2 and cross(lo[-2],lo[-1],p)<=0: lo.pop()
        lo.append(p)
    hi=[]
    for p in reversed(pts):
        while len(hi)>=2 and cross(hi[-2],hi[-1],p)<=0: hi.pop()
        hi.append(p)
    return lo[:-1]+hi[:-1]


def inside(point, polygon):
    x,y=point; hit=False; j=len(polygon)-1
    for i,(xi,yi) in enumerate(polygon):
        xj,yj=polygon[j]
        if (yi>y)!=(yj>y) and x<(xj-xi)*(y-yi)/(yj-yi)+xi: hit=not hit
        j=i
    return hit


def main():
    with DATA.open(encoding="utf-8-sig",newline="") as f: rows=list(csv.DictReader(f))
    pts=[tuple(float(r[c]) for c in ("X","Y","Z")) for r in rows]
    support=json.loads(SUPPORT.read_text(encoding="utf-8"))
    nearest_base=support["nearest_sample_spacing_units"]["3d"]["q3"]
    kth_base=support["nearest_sample_spacing_units"]["nearest_sample_from_different_hole_3d"]["median"]
    mins=[min(p[i] for p in pts) for i in range(3)]; maxs=[max(p[i] for p in pts) for i in range(3)]
    axes=[[mins[i]+j*STEP for j in range(math.ceil((maxs[i]-mins[i])/STEP)+1)] for i in range(3)]
    poly=hull([(p[0],p[1]) for p in pts]); candidates=[]
    for x in axes[0]:
        for y in axes[1]:
            if inside((x,y),poly):
                for z in axes[2]:
                    nearest=sorted(math.dist((x,y,z),p) for p in pts)
                    candidates.append((nearest[0],nearest[11]))
    nearest_limits=[nearest_base,2*nearest_base, max(support["nearest_sample_spacing_units"]["3d"]["max"],2*nearest_base)]
    kth_limits=[kth_base,2*kth_base,3*kth_base]
    combinations=[]
    for a in nearest_limits:
        for b in kth_limits:
            count=sum(d1<=a and d12<=b for d1,d12 in candidates)
            combinations.append({"max_nearest_distance":a,"max_12th_neighbor_distance":b,"accepted_nodes":count,"coverage_percent":100*count/len(candidates)})
    out={"candidate_nodes":len(candidates),"grid_spacing":STEP,"domain":"same XY convex hull and observed Z extent used by copper_idw_grid.py","sensitivity_only":True,"baseline_cutoffs":{"nearest":nearest_base,"12th_neighbor":2*kth_base},"cutoff_basis":{"nearest_levels":["observed nearest-spacing Q3","2x observed nearest-spacing Q3","max observed nearest-spacing (or 2x Q3 if larger)"],"12th_neighbor_levels":["1x","2x","3x median distance to nearest sample from another hole"]},"results":combinations,"caveat":"Changing a distance cutoff changes geometric coverage only. It does not validate geology, interval support, interpolation accuracy, or resource confidence."}
    OUT.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__": main()
