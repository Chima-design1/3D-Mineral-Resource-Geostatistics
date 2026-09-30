# 3D Geomodeling and Geostatistical Mineral Resource Estimation Using Public Drillhole Data

**Portfolio technical report — 30 September 2026**  
**Objective:** demonstrate reproducible geospatial data management, 3D sample QA/QC, spatial validation, interpolation, and uncertainty-aware decision making relevant to MSc Geomatics for Mineral Resource Management preparation.

## Executive summary

The project uses the public GeoMet drillhole dataset from Zenodo. Its observed drillhole table has 2,000 sample-centroid records, 119 `HOLEID` trajectories, XYZ coordinates, and 18 assay variables in ppm. QA/QC found no missing values, exact duplicate records, or duplicate XYZ triplets. The source does not supply collar/survey/interval tables, geological domains, density, a coordinate reference system, or assay detection-limit metadata.

Spatial models were assessed using whole-hole holdouts and buffered XY-tile holdouts. On the selected 300-unit tile / 50-unit buffer design, fixed IDW did not beat the training-mean baseline, and nested kriging's median MAE improvement over the mean was small and uncertain across 22 tiles. High-Cu screening performance from spatial models was limited. A ridge model using the 17 other assays performed much better, but it answers a conditional question—predicting Cu when co-located chemistry is already available—and cannot be treated as estimation into fully unassayed blocks.

An exploratory 3D IDW grid and a proxy surface clip were created to make spatial support limits visible. Only 366 of 9,270 candidate nodes passed the provisional distance screen; 355 remained below a surface proxy derived from first-row samples per hole. The grid's support thresholds were not validated by spatial holdouts, so the grid is a visualization prototype rather than a resource estimate.

## Data and provenance

- **Source:** [GeoMet dataset, Zenodo](https://doi.org/10.5281/zenodo.7051975), current project copy checksum-verified against the source record.
- **Related paper:** Hoffimann et al. (2022), [Modeling Geospatial Uncertainty of Geometallurgical Variables with Bayesian Models and Hilbert–Kriging](https://doi.org/10.1007/s11004-022-10013-1).
- **Observed drillhole table:** 2,000 rows × 22 columns; `HOLEID`, X/Y/Z, and 18 analytes in ppm. The CSV contains sample-centroid coordinates, not interval lengths or surveyed collars.
- **Companion tests:** `comminution.csv` (60 rows) and `flotation.csv` (53 rows) were checksum-verified. Their tables have no explicit shared sample/interval key and no exact XYZ matches to assay rows. Nearest-coordinate links remain candidates only; no test-response labels were transferred.
- **Schema and retrieval details:** [data provenance](../docs/provenance.md) and [data dictionary](../docs/data_dictionary.md).

- **Source-gap review:** the current v4 record lists only the three analyzed CSVs. The accompanying worked example describes XYZ as sample-centroid coordinates. A separate DOI cited by the paper for obfuscated copies remains uninspected; see the [source-gap audit](geomet_source_gap_audit.md).

The associated [GeoMet worked example](https://juliaearth.github.io/geospatial-data-science-with-julia/12-mining.html) documents a geometric interpolation footprint and a simple terrain proxy from first sample points. This project reproduced that surface idea as a screening proxy only; it is not official topography or a geological boundary.

## QA/QC and spatial characterization

- No missing values, exact duplicate rows, or duplicate XYZ triplets in the 2,000-row drillhole file.
- Cu concentration ranges from 0 to 60,200 ppm; Q1/median/Q3 are approximately 2,800/5,000/8,900 ppm.
- 119 hole IDs, with a median 16 samples per hole.
- 22 of 48 possible 300-unit XY tiles contain samples. Median nearest 3D sample spacing is 15.9 supplied coordinate units; coordinate units and CRS are undocumented in the CSV.
- The 18 reported analytes sum to a median 361,988 ppm, with no rows near 1,000,000 ppm. They are not a closed whole-rock composition. A Cu–Au–Ag–S subcomposition has zeros in seven records, and the source does not state whether these represent below-detection values.

Detailed reports: `reports/initial_qaqc.json`, `reports/spatial_support_audit.json`, and `reports/composition_readiness.json`.

## Modeling and validation

### Validation designs

1. **Leave-one-HOLEID-out:** holds out all samples from one trajectory.
2. **Buffered XY tiles:** holds out occupied XY grid tiles and removes nearby training points within a horizontal buffer. The principal comparison uses 300-unit tiles with a 50-unit buffer.
3. **Sensitivity:** fixed IDW and mean baselines were compared across tile widths 200/300/500 and buffers 0/50/100. These are relative source-coordinate units, not confirmed metres.

All model selection and transformations were intended to use training-fold data. Metrics include MAE, RMSE, bias, high-grade errors, and threshold-screening performance. The 8,900 and 18,100 ppm thresholds are descriptive top-quartile/top-decile cutoffs, not economic cutoffs.

### Main results

| Task and validation | Method | MAE (ppm) | RMSE (ppm) | Interpretation |
|---|---|---:|---:|---|
| Leave-one-hole-out | Training mean | 5,314 | 7,734 | Reference baseline |
| Leave-one-hole-out | Raw Cu IDW, p=2, k=12 | 5,174 | 7,562 | Slight improvement on this split |
| Buffered XY tiles, 300/50 | Training mean | 5,353 | 7,804 | Spatial reference baseline |
| Buffered XY tiles, 300/50 | Raw Cu IDW, p=2, k=12 | 5,479 | 8,046 | Worse than mean on this split |
| Buffered XY tiles, 300/50 | Nested ordinary kriging | 5,200 | 7,624 | Small pooled improvement; tile-bootstrap MAE difference vs mean median −149 ppm, interval −549 to +171 ppm |
| Buffered XY tiles, 300/50 | Conditional all-assay ridge, alpha 0.1–100 | 1,795 | 2,897 | Strong conditional chemistry prediction; assumes the other 17 assays are known at each target location |
| Buffered XY tiles, 300/50 | Conditional all-assay ridge, alpha 100–10,000 | 1,898 | 3,248 | Same predictors and folds, different regularization search range |
| Buffered XY tiles, 300/50 | Cu–Au–Ag–S CLR/total IDW sensitivity | 5,146 | 8,192 | Lower MAE, but worse RMSE, high-grade MAE, and bias than the mean |

Both conditional all-assay variants use the same 17 companion assay predictors and buffered folds, but search different ridge regularization ranges. They are reported separately to preserve reproducibility; neither is conventional spatial estimation into locations with no assays. Target-location dropout scenarios show that performance depends on which assays remain available. Companion LCT modeling on 52 rows did not beat its mean baseline.

The assay-availability sensitivity found that hypothetical target-location masking of Ag raised all-17 ridge MAE to 2,734 ppm, while masking S/Ag/Au raised it to 3,866 ppm; masking S alone did not degrade the score. These are availability stress tests, not observed missing-data patterns. See the [assay-availability decision note](assay_availability_decision.md).

![Observed-versus-predicted Cu validation diagnostics for the spatial and conditional models](model_validation_figure.svg)

The conditional model panels should be read as assay-assisted prediction: the held-out Cu value is predicted with the other 17 measured assays supplied. Their apparent advantage does not transfer to fully unassayed locations.

### High-Cu screening

On buffered holdouts, threshold classifiers and top-ranked review workloads were evaluated using the descriptive 18,100 ppm top-decile threshold. Ranking is done **within each held-out tile** to avoid comparing absolute scores from separately fitted folds.

- Reviewing the top 20% within each tile captures 22.3% of observed top-decile samples with nested kriging, 19.3% with raw IDW, and 23.3% with the subcomposition model.
- The alpha 0.1–100 conditional ridge captures 62.9% within its top 10% and 88.1% within its top 20%; the alpha 100–10,000 variant captures 62.4% and 86.1%. Both require the other assays at target locations.
- Pooled cross-fold rankings gave higher capture rates because model score scales differed by fold; pooled rates are retained only as a diagnostic.
- Paired bootstrap intervals across 22 spatial folds indicate that spatial-only MAE differences versus the mean are uncertain. The intervals are descriptive and tiles may not be independent.

See `reports/grade_threshold_classification.json`, `reports/top_rank_screening.json`, and `reports/spatial_fold_bootstrap.json`.

### Residual spatial structure

A 3D k-nearest-neighbour Moran's I diagnostic on the same 2,000 out-of-fold predictions found positive autocorrelation in both signed and absolute residuals across k=4, 8, and 12. The pattern remains after permutations are restricted within outer folds. Conditional chemistry reduces residual clustering relative to spatial-only kriging in the signed errors, but does not eliminate it; absolute-error clustering also remains. The test is descriptive because the coordinate system and units are undocumented and geological domains are unavailable. See the [residual spatial-autocorrelation note](spatial_residual_autocorrelation.md) and `reports/residual_spatial_autocorrelation.json`.

![Paired difference in absolute Cu prediction error between nested kriging and conditional all-17 assay ridge at the same held-out samples](paired_error_comparison.svg)

On the same 2,000 held-out samples, the low-alpha conditional ridge had lower absolute error than nested kriging on 81.0% overall and all 202 samples above the descriptive top-decile cutoff. Its mean absolute-error reduction was 13,100 ppm for those top-decile samples. This comparison assumes the other 17 assays are available at the target and is summarized by grade range in the [paired error note](paired_model_error_comparison.md); it is not spatial-only performance.

## Exploratory grid and domain limits

![Observed Cu assays and the support-screened grid in three Z slices](portfolio_spatial_figure.svg)

The map colors are log-scaled for visibility; each view is an XY projection, with grid estimates separated into three Z slices.

The 3D IDW grid uses 50-unit spacing, p=2, k=12, an XY convex hull, the observed Z range, and provisional distance limits derived from sample spacing. Only **366 of 9,270** candidate nodes passed those limits (**3.95%**). Buffered holdout validation had **zero** samples within the current support band, so the cutoff is not empirically validated. Across tested cutoff combinations, apparent grid coverage ranged from 0.5% to 34.1%.

A surface proxy based on the first sample row per hole (118/119 groups are top-down in Z) removed 11 additional nodes, leaving **355**. This proxy is not a surveyed digital terrain model. Neither the convex hull nor the surface proxy defines a geological or mineralized domain. No uncertainty estimate, block support, tonnage, resource class, mine plan, or economic value is reported.

Exploratory outputs: `outputs/copper_idw_grid.csv` and `outputs/copper_idw_grid_below_proxy_surface.csv`; method records: `reports/copper_idw_grid.json` and `reports/surface_clip_audit.json`.

## Project contribution

The project demonstrates geospatial data provenance, 3D sample QA/QC, whole-hole and spatially buffered validation, variogram/kriging experimentation, interpolation support diagnostics, geochemical prediction, and careful separation of conditional assay imputation from resource estimation. Its methodological contribution is an end-to-end analysis with explicit validation limits, not a claim that public data supports a classified mineral resource.

## Reproduction

- Open `notebooks/01_geomet_project_walkthrough.ipynb` for the guided results review.
- Run individual analysis scripts from the project root; the command list is in `README.md`.
- Raw source files are kept under `data/raw/` and excluded from version control; provenance/checksums are in `docs/provenance.md`.
- Reports are saved under `reports/`; exploratory CSV outputs are saved under `outputs/`.

## Priority data required for a defensible next model

1. Collar and downhole survey tables, with CRS and vertical datum.
2. Sample from/to depths, support lengths, and compositing definitions.
3. Lithology, alteration, mineralization, and interpreted estimation domains.
4. Assay QA metadata, analytical methods, and detection limits/zero semantics.
5. Density and a defensible block size/support definition.
6. An authoritative sample crosswalk for the comminution and flotation test tables.

Until these are available, results remain exploratory assay-concentration analyses; they do not support resource classification, tonnage, economic decisions, or mine planning.
