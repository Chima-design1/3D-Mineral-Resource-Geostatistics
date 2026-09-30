# Initial GeoMet drillhole QA/QC and baseline

## Retrieval and integrity

`data/raw/drillholes.csv` was downloaded from the public [Zenodo record](https://zenodo.org/records/7051975), DOI 10.5281/zenodo.7051975. The file is 362,226 bytes; its MD5 `bf48a6b1135113b6e6cf7b32f95315ad` matches the source record. SHA-256: `e0636576be0a0f9d2d270f34871ae16be8a42ff35457ae2c0179304de1abdc10`.

## Structure and integrity checks

- Dimensions: 2,000 rows × 22 columns.
- Fields: `HOLEID`, `X`, `Y`, `Z`, `Ag ppm`, `Al ppm`, `Au ppm`, `C ppm`, `Ca ppm`, `Cl ppm`, `Cu ppm`, `F ppm`, `Fe ppm`, `K ppm`, `Mg ppm`, `Mn ppm`, `Na ppm`, `P ppm`, `Pb ppm`, `S ppm`, `Th ppm`, `U ppm`.
- Missing values: zero in every field; coordinate parse failures: zero.
- Exact duplicate rows: 0; duplicate numeric coordinate triplets: 0.
- Distinct hole IDs: 119. Samples per hole: minimum 1, median 16, mean 16.81, maximum 41.
- X range: −1,150.00 to 954.58; Y: −882.90 to 756.51; Z: −399.92 to 419.30. Dataset documentation calls these local Cartesian sample-centroid coordinates; a geodetic CRS/datum is not supplied in the reviewed metadata.
- Consecutive row-order centroid spacing within holes (diagnostic): 1,881 adjacent pairs; median 18.05 coordinate units, IQR 13.25–25.00, maximum 190.96. Verify any large gaps and the row ordering before treating this as sample support or downhole interval length.

## Candidate target distributions

All values below are ppm. The IQR fences are a screening flag only: flagged values are not automatically errors and were not removed.

| Variable | Min | Q1 | Median | Mean | Q3 | Max | Outside 1.5×IQR fences |
|---|---:|---:|---:|---:|---:|---:|---:|
| Cu | 0 | 2,800 | 5,000 | 7,540.4 | 8,900 | 60,200 | 202 |
| Au | 0 | 0.09 | 0.21 | 0.471 | 0.51 | 5.89 | 210 |
| Ag | 0.01 | 0.60 | 1.14 | 1.862 | 2.25 | 15.38 | 184 |
| S | 100 | 900 | 1,700 | 2,712.95 | 3,125 | 21,500 | 209 |

The candidate grade-like distributions are strongly right-tailed (means exceed medians), so inspect histograms, assay reporting limits, and high values before selecting transformations or caps. Copper remains the provisional first target because it aligns with the research example, not because QA/QC establishes it as an economic resource grade.

## First baseline: leave-one-hole-out

All samples from each `HOLEID` were excluded together, yielding 2,000 out-of-fold predictions. The local mean, 3D nearest neighbour, and isotropic IDW (power 2, 12 neighbours) were compared on the same folds.

| Method | MAE (ppm) | RMSE (ppm) | Mean error (ppm) |
|---|---:|---:|---:|
| Training mean | 5,314.2 | 7,734.0 | −4.8 |
| 3D nearest neighbour | 6,299.4 | 9,785.3 | −119.5 |
| IDW, power 2, 12 neighbours | 5,173.6 | 7,561.9 | 264.6 |

IDW is the best of these three simple baselines on this grouped split, with modest improvement over the training mean. This is not yet evidence for kriging: whole-hole folds may be spatially interspersed, and isotropic 3D distance is an unvalidated baseline assumption.

## Spatial exploration and spatial-block baseline

A plan view and X-Z section are in `spatial_eda.svg`, with Cu color shown on a log scale. Experimental variograms were computed using 20-unit lag bins, 500-unit maximum lag, and four horizontal azimuths plus a near-vertical direction (22.5° tolerance). Raw-Cu semivariances at the 110-unit bin ranged from about 51 to 61 million ppm² across directions; `log1p(Cu)` semivariances ranged from 0.67 to 1.18. The Cressie-Hawkins robust raw-Cu values ranged from about 46.8 to 65.3 million ppm². Directional variation is visible, but estimates remain noisy; transformation choice and geometry affect apparent continuity. These exploratory curves do not yet support fitting a stable anisotropic model.

An exploratory spatial-block validation leaves one occupied 300-unit XY grid tile out at a time (22 folds). Without a buffer, IDW (MAE 5,324.2 ppm; RMSE 7,699.6 ppm) was essentially tied with training mean (MAE 5,343.2 ppm; RMSE 7,778.7 ppm). Adding a 50-unit horizontal buffer around each held-out tile changed the result: IDW MAE/RMSE were 5,478.6/8,046.1 ppm versus 5,353.4/7,804.5 ppm for the training mean. Here IDW is worse than the mean, suggesting it does not reliably generalize to the spatial gaps created by this buffer. This does not prove kriging will improve results.

## Buffered ordinary-kriging experiment

Candidate isotropic variograms (spherical, exponential, Gaussian) were fitted independently from each training fold's empirical semivariogram (40,000 sampled pairs per fold; 20-unit bins to 500 units), then used for local ordinary kriging with 12 nearest training samples. The 22 folds held out occupied 300-unit XY tiles and excluded training samples within 50 horizontal units of tile boundaries. Holdout grades were used only for scoring.

| Method | MAE (ppm) | RMSE (ppm) | Mean error (ppm) |
|---|---:|---:|---:|
| Buffered training mean | 5,353.4 | 7,804.5 | −44.0 |
| Buffered IDW, p=2, k=12 | 5,478.6 | 8,046.1 | −794.3 |
| Local ordinary kriging, raw Cu | 5,427.2 | 8,001.3 | −800.3 |
| Local ordinary kriging, log1p with lognormal mean correction | 5,862.7 | 8,266.7 | 553.8 |

Neither kriging variant beats the buffered training-mean baseline; raw-Cu kriging is marginally better than IDW but remains worse than the mean. The log1p back-transform performs worst on these pooled scores. Gaussian variograms were selected in 21/22 raw-Cu folds and 20/22 log1p folds by the training-only weighted fit criterion. This is exploratory: parameters were selected by in-sample variogram fit rather than nested spatial CV; folds share samples and are not independent replications; horizontal buffering is combined with full 3D distances; the coordinate CRS, geology/domains, and sample support are not documented.

Influence screening found 4 zero Cu values, q99=36,000 ppm, and 202 values above the 1.5×IQR upper fence of 18,050 ppm. The maximum is 60,200 ppm at row 512, hole 33, XYZ (444.60, −663.83, 274.25). Five consecutive row-order within-hole gaps exceed 100 coordinate units. These values and gaps are flagged for source-record review only; none were removed or capped. See `sample_influence.json`.

Fold diagnostics show raw-Cu kriging errors are heterogeneous (fold MAE 761–8,839 ppm). The correlation between median nearest-training distance and fold MAE is weakly positive (r≈0.27), so distance alone does not explain difficult tiles. Results and per-fold variograms are in `fold_diagnostics.json`.

A first nested spatial sensitivity selected neighborhood count from 4, 8, 12, and 20; all 22 folds selected k=20, so the search was expanded. The expanded nested search compared k=8/16/24/32, spherical/Gaussian families, and ranges 80/200/400 units using two geographically distributed buffered inner tiles per outer training set. On untouched outer tiles, pooled MAE/RMSE were 5,200.2/7,624.2 ppm, compared with 5,353.4/7,804.5 ppm for the same-split training-mean baseline and 5,478.6/8,046.1 ppm for IDW. Improvement over the mean is modest (about 2.9% MAE), and mean error is −503 ppm versus −44 ppm for the mean baseline. Twenty folds selected k=32, again at the upper search boundary; 14 selected spherical and 14 selected range 200. Therefore neighborhood stability remains unresolved and family/range choices vary by fold. The search used 10,000 sampled pairs per fitted variogram. Full fold selections and out-of-fold predictions are in `nested_kriging_sensitivity.json`.

Residuals expose substantial grade shrinkage: the lowest observed-Cu quartile is overpredicted by a mean 4,510 ppm, while the highest quartile (9,000–60,200 ppm) is underpredicted by a mean 10,267 ppm and has MAE 10,506 ppm. This is consistent with local smoothing and warns against interpreting the pooled improvement as reliable high-grade estimation. See `residual_diagnostics.json` and the plan/section view `nested_residuals.svg`.

## High-grade context and geochemical covariates

The observed `drillholes.csv` fields are `HOLEID`, X/Y/Z centroids, and 18 chemical concentrations. It contains no lithology, alteration, geological domain, interval from/to or length, recovery, density, collar/survey, or assay QA/QC fields. The Zenodo record separately lists comminution and flotation test tables, but they are not part of this Cu analysis. The source paper frames those tables as geometallurgical laboratory responses, and notes that support comparability between test samples and drillhole samples is an assumption; it does not provide geological domain labels for this CSV.

Top-decile Cu (≥18,100 ppm; 202 samples) occurs across multiple holes and XY tiles rather than a single hole. Among log1p concentrations, the strongest Pearson associations with log1p(Cu) are S (r=0.857), Ag (0.830), and Au (0.712). These are data-driven chemical associations, not validated geological domains, and do not establish causality.

A nested buffered spatial test used all 17 non-Cu assays as co-located covariates, with log1p transforms, training-only standardization, and ridge strength chosen from inner spatial folds. When the alpha grid was expanded to 100–10,000, the all-element model achieved outer MAE/RMSE 1,898/3,248 ppm versus 5,353/7,804 ppm for the same-split mean. In the top observed Cu quartile (≥8,900 ppm), MAE was 3,787 ppm versus 10,543 ppm for the mean prediction. A pre-specified parsimonious S/Ag/Au model scored 2,075/3,593 ppm overall and 4,294 ppm MAE in that high quartile. All 17 elements performed better in this comparison. Inner folds selected alpha=1,000 in 18/22 all-element folds and 21/22 S/Ag/Au folds, showing that the earlier alpha=100 boundary result was sensitive to the tested tuning grid. This remains a covariate-assisted prediction result and assumes the other element assays are available where Cu is withheld; it does not demonstrate estimation into unassayed blocks. See `geochem_feature_sensitivity.json`.

## Next modeling step

An availability-scenario comparison confirms the distinction: the drillhole CSV is complete case (zero missing entries in every observed column), so the scenarios are counterfactual. On the same spatial outer folds, all-element ridge with non-Cu assays available scored 1,898/3,248 ppm MAE/RMSE; S/Ag/Au-only scored 2,075/3,593; and raw-Cu kriging with no target-location assays scored 5,200/7,624. For the top observed Cu quartile, corresponding MAEs were 3,787, 4,294, and 10,407 ppm. Here “no assays at target” still allows Cu data from neighboring training samples. This comparison supports a clear use distinction but does not establish realistic missingness, because no assay-availability pattern is present in the source. See `assay_availability_scenarios.json`.

A further target-location dropout stress test masks companion assays to their training-fold means while retaining complete training assays, and retunes ridge strength in inner spatial folds. Removing Ag only at prediction locations raised overall MAE from 1,898 to 2,734 ppm and top-quartile MAE from 3,787 to 5,571 ppm. Masking Au raised those to 2,219 and 4,482 ppm. Masking S alone did not worsen scores (MAE 1,822 ppm), indicating S's strong marginal association is not the same as unique predictive value when other elements are present. Masking all three S/Ag/Au increased MAE to 3,866 ppm overall and 7,650 ppm in the top Cu quartile. These synthetic dropout patterns quantify sensitivity only; the source contains no observed missingness process. See `geochem_assay_dropout.json`.

The next concrete step is to obtain a documented assay-coverage pattern (or validation data with intentional assay omissions) and rerun the nested evaluation against that pattern. The current complete-case CSV cannot establish operational missingness; the stress scenarios are hypothetical. To move toward a geological model, obtain lithology/domain and interval-support data; the current CSV does not support domain-aware estimation. Additional geological domains, density, and resource reporting inputs are absent, so do not report tonnage, classification, or economic value.

## Companion geometallurgical test tables

The Zenodo record's `comminution.csv` and `flotation.csv` have now been retrieved and their source MD5 hashes verified. They contain 60 comminution rows (28 columns; `th1`–`th3`, F80, P80, M, A plus assays) and 53 flotation rows (26 columns; `fr`, `xr`, `LCT` plus assays). Missing test-table values include 34 missing Th and 12 missing graphite-carbon values in comminution, and one missing LCT plus four missing Cl values in flotation.

Both test tables contain `HOLEID` and X/Y/Z. However, only 43/60 comminution and 51/53 flotation rows have an ID present in the drillhole file, and none of those has an exact XYZ match. For records with a shared ID, nearest centroid distances have medians 3.03 and 4.36 coordinate units, respectively; most are within 20 units. This suggests plausible spatial correspondence but does not prove that the laboratory response belongs to the nearest drillhole assay support. There are no sample-length, from/to, or shared sample-key fields in the observed tables. No test labels have been merged onto drillhole rows. See `docs/companion_data_audit.md`, `companion_table_audit.json`, and `companion_link_candidates.json`.

Before modeling DWT/BWI/LCT on the 2,000 drillhole points, obtain a sample crosswalk and support/compositing definition, plus field definitions for the lab variables. If those are unavailable, nearest-coordinate pairing can only be explored as a sensitivity analysis, never treated as ground truth.

Separately, a chemistry-only pilot on the 52 flotation rows with observed LCT used nested leave-one-HOLEID-out validation. Ridge regression scored MAE/RMSE 0.0475/0.0576; the training-mean baseline scored 0.0434/0.0544. Chemistry alone did not improve held-out-hole prediction in this small sample. This model does not require a crosswalk because its response and chemistry are recorded on the same flotation row, but it is conditional lab-response prediction, not block interpolation. See `flotation_response_cv.json`.

Full machine-readable results are in `initial_qaqc.json`, `baseline_cv.json`, `spatial_eda.json`, `fold_kriging.json`, `sample_influence.json`, `fold_diagnostics.json`, `nested_kriging_sensitivity.json`, `residual_diagnostics.json`, `high_grade_context.json`, `geochem_covariate_cv.json`, `geochem_feature_sensitivity.json`, `assay_availability_scenarios.json`, `geochem_assay_dropout.json`, `companion_table_audit.json`, `companion_link_candidates.json`, and `flotation_response_cv.json`. Reproduce companion audits and the LCT pilot with `python scripts/companion_table_audit.py` and `python scripts/flotation_response_cv.py`.
