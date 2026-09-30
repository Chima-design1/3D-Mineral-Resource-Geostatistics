"""Clip the exploratory Cu IDW points below a proxy surface from first rows."""
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"/"raw"/"drillholes.csv"
GRID=ROOT/"outputs"/"copper_idw_grid.csv"
OUT=ROOT/"outputs"/"copper_idw_grid_below_proxy_surface.csv"
META=ROOT/"reports"/"surface_clip_audit.json"
K,POWER=12,2.0


def read(path):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))


def main():
    rows=read(DATA); groups=defaultdict(list)
    for r in rows: groups[r["HOLEID"]].append(tuple(float(r[c]) for c in ("X","Y","Z")))
    collars=[(hole,ps[0],ps[-1]) for hole,ps in groups.items()]
    first_above=sum(first[2]>last[2] for _,first,last in collars)
    first_below=sum(first[2]<last[2] for _,first,last in collars)
    same=sum(first[2]==last[2] for _,first,last in collars)
    zsource=[r["Z"] for r in read(GRID)]
    kept=[]; discarded=0; surfaces=[]
    for r in read(GRID):
        x,y,z=(float(r[c]) for c in ("X","Y","Z"))
        near=sorted((math.hypot(x-first[0],y-first[1]),first[2]) for _,first,_ in collars)[:K]
        if near[0][0]==0: surface=near[0][1]
        else:
            ws=[1/(d**POWER) for d,_ in near]
            surface=sum(w*v for w,(_,v) in zip(ws,near))/sum(ws)
        surfaces.append(surface)
        if z < surface:
            row=dict(r); row["proxy_surface_Z"]=f"{surface:.3f}"; row["depth_below_proxy_surface"]=f"{surface-z:.3f}"; kept.append(row)
        else: discarded+=1
    with OUT.open("w",encoding="utf-8",newline="") as f:
        fields=list(kept[0].keys()) if kept else list(read(GRID)[0].keys())+["proxy_surface_Z","depth_below_proxy_surface"]
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(kept)
    meta={"input_grid_nodes":len(zsource),"output_nodes_below_proxy_surface":len(kept),"excluded_as_above_or_on_proxy_surface":discarded,"holes":len(collars),"first_sample_vs_last_sample_z_order":{"first_above_last":first_above,"first_below_last":first_below,"equal":same,"interpretation":"First sample point has greater Z than the last for most hole groups, supporting a top-down ordering proxy; this does not establish surveyed collar elevations or terrain."},"proxy_surface":{"source_points":"first row per HOLEID in drillholes.csv","interpolation":"2D IDW on XY, power 2, 12 nearest HOLEID first-row centroids","output_columns":["proxy_surface_Z","depth_below_proxy_surface"],"CRS_and_vertical_datum":"Not supplied"},"output_csv":str(OUT.relative_to(ROOT)),"limitations":["This is a proxy surface from sample centroids and row ordering, not a surveyed digital terrain model.","The public worked example for GeoMet also derives a simple terrain surface from first trajectory points; it does not make those centroids equivalent to official collars.","The clipped interpolation domain is geometric and topographic only; it is not a geological or mineralized domain.","Existing IDW support limits remain unvalidated; surface clipping does not resolve spatial validation or resource-estimation limitations."]}
    META.write_text(json.dumps(meta,indent=2),encoding="utf-8")
    print(json.dumps(meta,indent=2))


if __name__=="__main__":main()
