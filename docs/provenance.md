# Data provenance

- Dataset: GeoMet dataset, version 4.
- Publisher: Zenodo.
- DOI: [10.5281/zenodo.7051975](https://doi.org/10.5281/zenodo.7051975).
- Record: [https://zenodo.org/records/7051975](https://zenodo.org/records/7051975).
- File of interest: `drillholes.csv`, record listing size 362.2 kB and MD5 `bf48a6b1135113b6e6cf7b32f95315ad`.
- Record description: curated geometallurgical dataset; this file contains chemical analyses and sample geospatial coordinates.
- License: the live record displays a blank License field (checked 2026-09-30). This project does not assume redistribution permission; the raw CSVs are excluded from version control. Reuse and redistribution terms should be clarified with the record authors before republishing source files.
- Retrieved 2026-09-30 from the public Zenodo file URL to `data/raw/drillholes.csv` (362,226 bytes). The source-published MD5 was verified against the downloaded file: `bf48a6b1135113b6e6cf7b32f95315ad`. A reproducible local report now records SHA-256 and the full observed schema and QA/QC in `reports/initial_qaqc.json`.
- Schema/context source: [Geospatial Data Science with Julia, Mineral deposits chapter](https://juliaearth.github.io/geospatial-data-science-with-julia/12-mining.html), which documents centroid X/Y/Z, HOLEID, chemical concentrations in ppm, 2,000 samples, and examples including Cu, Au, Ag, S. These are secondary documentation observations and must be reconciled against the CSV itself.
- The same public chapter defines an illustrative interpolation footprint using the XY convex hull and estimates a simple terrain surface from the first sample point per hole. We replicated only the first-point surface as an explicitly labeled screening proxy in `scripts/clip_copper_grid_to_surface.py`; it is not an official collar or surveyed digital terrain model.
- The Zenodo record also lists `comminution.csv` (25,239 bytes; MD5 `1a33df8ba77f5d49c281ba3fff16b20b`) and `flotation.csv` (17,531 bytes; MD5 `f2e90da6bfa81de1261177ee85a91570`). They were retrieved from the public Zenodo file previews on 2026-09-30 and the published MD5 values were verified locally. The record describes the tables as DWT/BWI and copper locked-cycle flotation measurements, respectively.
- The associated paper is Hoffimann et al. (2022), [Modeling Geospatial Uncertainty of Geometallurgical Variables with Bayesian Models and Hilbert-Kriging](https://doi.org/10.1007/s11004-022-10013-1). It explains the use of assay chemistry to predict test responses and discusses sample support assumptions. The published drillhole table itself does not contain geological domains or interval-length fields, so those cannot be reconstructed from the available CSV alone.
- A structural and candidate-link audit is in `reports/companion_table_audit.json`, with per-row nearest-coordinate candidates in `reports/companion_link_candidates.json`. The companion files have not yet been merged with drillhole assays; no modeling result depends on them.

When the file is downloaded, keep the original unchanged at `data/raw/drillholes.csv`; run `scripts/initial_qaqc.py` to record SHA-256 and observed structure in `reports/initial_qaqc.json`.
