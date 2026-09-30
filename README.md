# Exploratory 3D geostatistical analysis of public GeoMet drillhole data

This project audits the public GeoMet v4 tables and evaluates spatial prediction of copper concentration from sample-centroid data. It demonstrates reproducible data QA/QC, spatial validation, geostatistical sensitivity analysis, and support-aware visualization.

**Scope:** exploratory sample-level concentration prediction only. The available tables do not provide geology/domains, interval lengths, collars and surveys, density, coordinate reference system, or assay QA/QC. Results are not a block model, mineral resource estimate, classification, or economic assessment.

## Data and provenance

The source is the [GeoMet v4 Zenodo record](https://doi.org/10.5281/zenodo.7051975). The local drillhole table contains 2,000 rows and 22 columns; the two associated lab tables contain 60 and 53 rows. Checksums, retrieval notes, and schema are in [data provenance](docs/provenance.md); the observed drillhole fields are documented in the [data dictionary](docs/data_dictionary.md).

The Zenodo record currently displays no license value. Raw source files are therefore excluded from this repository. Obtain them from the source record and place the original CSVs under `data/raw/`; do not redistribute them from this project unless reuse rights are clarified.

## Reproduction

Requires Python 3.11 or newer. The analysis scripts use only the Python standard library; no third-party packages are required.

From the project root:

```powershell
python scripts/run_all.py
```

The pipeline verifies the downloaded drillhole structure, runs analyses in dependency order, then generates JSON reports under `reports/`, grid tables and figures under `outputs/`, and the [walkthrough notebook](notebooks/01_geomet_project_walkthrough.ipynb). The row-level reports and generated visual/data files are excluded from publication while GeoMet reuse terms are unspecified; they remain reproducible locally. The full run can take several minutes because nested spatial validation repeats variogram fitting.

For the full analysis inventory and interpretation, see [the modeling plan](docs/modeling_plan.md), [companion-table audit](docs/companion_data_audit.md), [source gap audit](outputs/geomet_source_gap_audit.md), and [technical report](outputs/portfolio_technical_report.md).

The search for richer collar/survey/interval data is documented in [public drillhole dataset screening](docs/dataset_screening.md). NTGS DIP001 is a promising licensed candidate, but its 1.04 GB archive and actual table schemas have not yet been inspected, so it is not integrated into the analysis.

## Key findings and limitations

- The drillhole file has no missing values, exact duplicate rows, or duplicate XYZ triplets; it has 119 hole identifiers.
- Spatial-only Cu prediction is sensitive to validation design and does not consistently outperform simple mean baselines. High-Cu predictions show shrinkage and residuals retain spatial structure.
- Shifting the origin of the fixed buffered spatial-CV grid changed IDW MAE from 5,248 to 5,763 ppm; IDW beat the training mean on MAE in 2 of 9 tested origins. See [the origin-sensitivity summary](outputs/spatial_origin_sensitivity.md).
- Conditional all-assay models predict Cu more accurately when the other assays are known at each withheld location. This is not a method for estimating into unsampled or unassayed blocks.
- The support-limited exploratory IDW grid covers only a small portion of candidate nodes. Its distance cutoffs are not validated as uncertainty thresholds.
- No geology-aware resource estimation or classification is supported by the available data.

## Repository contents

- `data/raw/`: locally obtained source tables; excluded from version control.
- `docs/`: provenance, field dictionary, companion-table audit, and methodological plan.
- `scripts/`: standard-library analysis and figure-generation scripts.
- `reports/`: generated QA/QC and analysis results (not published; regenerate locally).
- `notebooks/`: reproducible results walkthrough.
- `outputs/`: exploratory figures and grids (generated locally), plus aggregate technical summaries.

## Citation

Hoffimann, J. et al. (2022). *Modeling Geospatial Uncertainty of Geometallurgical Variables with Bayesian Models and Hilbert-Kriging*. DOI: [10.1007/s11004-022-10013-1](https://doi.org/10.1007/s11004-022-10013-1). Dataset: [GeoMet v4](https://doi.org/10.5281/zenodo.7051975).

