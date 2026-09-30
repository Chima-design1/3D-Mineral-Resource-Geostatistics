"""Verify GeoMet companion tables and assess (but do not assert) linkage."""
import csv, hashlib, json, math, statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/"data"/"raw"
FILES={"comminution.csv":"1a33df8ba77f5d49c281ba3fff16b20b","flotation.csv":"f2e90da6bfa81de1261177ee85a91570"}
DRILL=RAW/"drillholes.csv"

def read(path):
    with path.open("r",encoding="utf-8-sig",newline="") as f: return list(csv.DictReader(f))

def coord(r): return tuple(float(r[c]) for c in ("X","Y","Z"))

def main():
    drill=read(DRILL)
    drill_by_hole={}
    for n,r in enumerate(drill,1): drill_by_hole.setdefault(r["HOLEID"],[]).append((n,coord(r)))
    out={"drillholes":{"n":len(drill),"columns":list(drill[0]),"holes":len(drill_by_hole)}}
    candidates={}
    for filename,expected in FILES.items():
        path=RAW/filename; digest=hashlib.md5(path.read_bytes()).hexdigest(); rows=read(path)
        missing={c:sum(not r[c].strip() for r in rows) for c in rows[0] if any(not r[c].strip() for r in rows)}
        table={"n":len(rows),"columns":list(rows[0]),"distinct_hole_ids":len(set(r["HOLEID"] for r in rows)),"missing_values":missing,"source_md5_expected":expected,"local_md5":digest,"md5_verified":digest==expected,"duplicate_xyz_count":len(rows)-len({coord(r) for r in rows})}
        same=[]; any_near=[]; same_id=0; exact=0; unmatched=[]; perrow=[]
        for i,r in enumerate(rows,1):
            q=coord(r); matching=drill_by_hole.get(r["HOLEID"],[])
            any_d, any_j, any_dr=min((math.dist(q,coord(dr)),j,dr) for j,dr in enumerate(drill,1))
            any_near.append(any_d)
            if matching:
                same_id+=1
                d,j=min((math.dist(q,p),j) for j,p in matching)
                same.append(d); exact+=d<1e-9
                same_candidate={"row":j,"distance_xyz":d}
            else:
                unmatched.append({"test_row":i,"HOLEID":r["HOLEID"]})
                same_candidate=None
                perrow.append({"test_row":i,"test_HOLEID":r["HOLEID"],"nearest_same_HOLEID_candidate":same_candidate,"nearest_any_HOLEID_candidate":{"row":any_j,"HOLEID":any_dr["HOLEID"],"distance_xyz":any_d},"warning":"Candidate only; no direct sample identifier or interval key."})
        table.update({"rows_with_HOLEID_present_in_drillholes":same_id,"same_HOLEID_exact_xyz_matches":exact,"nearest_same_HOLEID_distance_units":{"n":len(same),"median":statistics.median(same) if same else None,"max":max(same) if same else None,"within_5":sum(d<=5 for d in same),"within_10":sum(d<=10 for d in same),"within_20":sum(d<=20 for d in same)},"nearest_any_drillhole_distance_units":{"median":statistics.median(any_near),"max":max(any_near),"within_10":sum(d<=10 for d in any_near),"within_20":sum(d<=20 for d in any_near)},"test_HOLEIDs_absent_from_drillholes":sorted(set(r["HOLEID"] for r in rows)-set(drill_by_hole)),"interpretation":"Spatially close points and matching HOLEID values are candidate linkage evidence only. Coordinates do not match exactly, and neither direct sample IDs nor sample intervals are supplied; do not automatically join or transfer target labels."})
        out[filename]=table; candidates[filename]=perrow
    (ROOT/"reports"/"companion_table_audit.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    (ROOT/"reports"/"companion_link_candidates.json").write_text(json.dumps(candidates,indent=2),encoding="utf-8")
    print(json.dumps({k:{"n":v["n"],"columns":len(v["columns"]),"missing":v["missing_values"],"hash_ok":v["md5_verified"],"same_id_rows":v["rows_with_HOLEID_present_in_drillholes"],"exact":v["same_HOLEID_exact_xyz_matches"],"same_id_distance":v["nearest_same_HOLEID_distance_units"],"unmatched_ids":v["test_HOLEIDs_absent_from_drillholes"]} for k,v in out.items() if k.endswith(".csv")},indent=2))

if __name__=="__main__": main()
