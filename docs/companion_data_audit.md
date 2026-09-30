# GeoMet companion-table audit

## Source and integrity

The public [Zenodo GeoMet record](https://zenodo.org/records/7051975) lists `comminution.csv`, `flotation.csv`, and `drillholes.csv`. The companion CSV files are now in `data/raw/`; their local MD5 hashes match the record metadata. The associated paper describes the comminution table as DWT/BWI laboratory data and the flotation table as locked-cycle flotation data.

## Observed tables

| File | Rows | Columns | Response/test fields observed | Missing values |
|---|---:|---:|---|---|
| `comminution.csv` | 60 | 28 | `th1`, `th2`, `th3`, `F80`, `P80`, `M`, `A` plus chemistry | Carbono Grafite 12; Cl 1; Th 34; U 2 |
| `flotation.csv` | 53 | 26 | `fr`, `xr`, `LCT` plus chemistry | LCT 1; Cl 4 |

Both include `HOLEID`, X/Y/Z, and chemical concentrations, so the test observations are spatially described. Neither table provides interval from/to, sample length, an explicit shared sample identifier, or geological domain. Do not infer meanings of the experimental fields beyond the published record without consulting its definitions.

## Candidate correspondence to drillholes

The file audit found 43/60 comminution records and 51/53 flotation records whose `HOLEID` occurs in `drillholes.csv`. None of those records has an exact XYZ match to a drillhole row. For matching IDs, nearest same-ID centroid distance has median 3.03 coordinate units for comminution (42/43 within 20) and 4.36 units for flotation (50/51 within 20). The unmatched test IDs are 120–133 for comminution and 132–133 for flotation.

These close coordinates make spatial correspondence plausible, but are not enough to prove that a test response belongs to the nearest drillhole assay interval. In addition, some IDs do not occur in the main drillhole file, and the tables contain no sample-length or interval keys. `reports/companion_link_candidates.json` stores nearest-coordinate candidates for manual review. The analysis deliberately does not merge test responses onto drillhole rows.

## Candidate ambiguity check

The follow-up in `reports/companion_link_sensitivity.json` ranks both the nearest and second-nearest drillhole centroid within the same `HOLEID` for every companion record. For comminution, 43/60 rows have a matching hole ID (42 are within 20 coordinate units); for flotation, 51/53 do (50 are within 20). Among rows with at least two same-hole candidates, the closest centroid is at least twice as close as the runner-up for 32/43 comminution and 41/51 flotation records. Two comminution records and one flotation record have a runner-up within one coordinate unit of the nearest.

This quantifies geometric ranking, not identity. Even an apparently isolated nearest centroid may represent a different sample interval or support, and coordinate units/CRS and sampling procedures are undocumented. No response-label transfer or combined model is justified by this sensitivity check.

## Next data gate

Before training a model to predict LCT/DWT/BWI responses on the 2,000 drillhole samples, obtain the source's sample-ID/crosswalk, interval support or compositing definition, coordinate reference and survey method, and meanings of `fr`, `xr`, `th1`–`th3`, `F80`, `P80`, `M`, and `A`. If the source has no crosswalk, any spatial nearest-neighbor transfer must be treated as an explicit sensitivity study, not ground-truth pairing. Also request assay QA metadata (detection limits, below-detection encoding, lab method) and geological/lithological domains if the intended endpoint is resource-oriented estimation.

## Next modeling decision

The immediate defensible step is to keep the two analyses separate: use drillhole Cu for the already completed spatial-estimation benchmark, and use the flotation table only for an explicitly limited same-table LCT prediction study. The LCT pilot did not beat its mean baseline, so do not expand that model until response definitions and a stronger validation rationale are available. For a response-on-drillholes study, the next required input is an authoritative crosswalk plus support definitions. For a geology-aware Cu model, the next required input is lithology/domain data and documented sample intervals. See `docs/modeling_plan.md` for the project-wide gates.

An initial chemistry-only LCT pilot was run directly within `flotation.csv`, without joining to drillholes. Nested leave-one-HOLEID-out ridge regression on the 52 non-missing LCT records produced MAE 0.0475 and RMSE 0.0576, slightly worse than the training-mean baseline (0.0434/0.0544). This small test set does not yet support the chemistry model; see `reports/flotation_response_cv.json`.
