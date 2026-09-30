"""Standard-library QA/QC for GeoMet drillholes.csv."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "drillholes.csv"
OUT = ROOT / "reports" / "initial_qaqc.json"
COORDS = ("X", "Y", "Z")
TARGETS = ("Cu ppm", "Au ppm", "Ag ppm", "S ppm")


def quantile(sorted_values: list[float], p: float) -> float:
    pos = (len(sorted_values) - 1) * p
    lo, hi = math.floor(pos), math.ceil(pos)
    return sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * (pos - lo)


def main() -> None:
    if not DATA.is_file():
        raise SystemExit(f"Place the source CSV at {DATA} and rerun.")
    with DATA.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    with DATA.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames or []
        rows = list(reader)
    if not columns:
        raise SystemExit("CSV has no header row.")

    missing = {c: sum(not (r.get(c) or "").strip() for r in rows) for c in columns}
    exact_duplicates = len(rows) - len({tuple((r.get(c) or "") for c in columns) for r in rows})
    report: dict[str, object] = {
        "source_file": DATA.name,
        "source_url": "https://zenodo.org/records/7051975/files/drillholes.csv?download=1",
        "record_doi": "10.5281/zenodo.7051975",
        "sha256": digest,
        "published_md5": "bf48a6b1135113b6e6cf7b32f95315ad",
        "dimensions": {"rows": len(rows), "columns": len(columns)},
        "columns": columns,
        "missing_by_column": missing,
        "exact_duplicate_rows": exact_duplicates,
        "coordinate_ranges": {},
        "coordinate_parse_failures": {},
        "duplicate_coordinate_rows": None,
        "unique_hole_ids": None,
        "samples_per_hole": None,
        "within_hole_centroid_spacing": None,
        "candidate_targets_present": [c for c in TARGETS if c in columns],
        "numeric_descriptives": {},
    }

    numeric: dict[str, list[float]] = {}
    for c in columns:
        parsed = []
        failures = 0
        for row in rows:
            raw = (row.get(c) or "").strip()
            if not raw:
                continue
            try:
                value = float(raw)
                if math.isfinite(value):
                    parsed.append(value)
                else:
                    failures += 1
            except ValueError:
                failures += 1
        if parsed or c in COORDS or c in TARGETS:
            numeric[c] = parsed
            if c in COORDS:
                report["coordinate_parse_failures"][c] = failures + missing[c]
                report["coordinate_ranges"][c] = {
                    "min": min(parsed) if parsed else None,
                    "max": max(parsed) if parsed else None,
                }
        if c != "HOLEID" and parsed:
            vals = sorted(parsed)
            q1, q3 = quantile(vals, .25), quantile(vals, .75)
            iqr = q3 - q1
            report["numeric_descriptives"][c] = {
                "n": len(vals), "missing": missing[c], "min": vals[0],
                "q1": q1, "median": statistics.median(vals),
                "mean": statistics.fmean(vals), "q3": q3,
                "max": vals[-1],
                "sample_sd": statistics.stdev(vals) if len(vals) > 1 else None,
                "outside_1_5_iqr_fences": sum(x < q1 - 1.5 * iqr or x > q3 + 1.5 * iqr for x in vals),
            }

    if all(c in columns for c in COORDS):
        coords = []
        for r in rows:
            try:
                coords.append(tuple(float(r[c]) for c in COORDS))
            except (ValueError, TypeError):
                continue
        report["duplicate_coordinate_rows"] = len(coords) - len(set(coords))
    if "HOLEID" in columns:
        groups: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
        for r in rows:
            hole = (r.get("HOLEID") or "").strip()
            try:
                xyz = tuple(float(r[c]) for c in COORDS)
            except (ValueError, TypeError, KeyError):
                continue
            if hole:
                groups[hole].append(xyz)
        counts = [len(v) for v in groups.values()]
        report["unique_hole_ids"] = len(groups)
        if counts:
            report["samples_per_hole"] = {
                "min": min(counts), "median": statistics.median(counts),
                "mean": statistics.fmean(counts), "max": max(counts),
            }
        # Row order is retained within each hole because the file provides centroid
        # samples in trajectory order; treat this as a diagnostic, not surveyed depth.
        distances = []
        for points in groups.values():
            for a, b in zip(points, points[1:]):
                distances.append(math.dist(a, b))
        if distances:
            distances.sort()
            report["within_hole_centroid_spacing"] = {
                "n_adjacent_pairs": len(distances), "min": distances[0],
                "q1": quantile(distances, .25), "median": statistics.median(distances),
                "mean": statistics.fmean(distances), "q3": quantile(distances, .75),
                "max": distances[-1],
            }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"\nReport saved: {OUT}")


if __name__ == "__main__":
    main()
