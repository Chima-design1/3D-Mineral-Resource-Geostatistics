"""Build the standard-library-only GeoMet portfolio walkthrough notebook."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "notebooks" / "01_geomet_project_walkthrough.ipynb"


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(keepends=True)}


def code(text):
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": text.splitlines(keepends=True)}


cells = [
    md("""# GeoMet 3D Geomodeling Portfolio Walkthrough

This notebook presents the verified GeoMet drillhole analysis as a reproducible walkthrough. It reads the original source tables and saved analysis reports; it does not claim a mineral resource estimate.

**Scope:** data QA/QC, Cu concentration prediction validation, spatial support screening, and an explicitly exploratory grid. The dataset supplies sample-centroid coordinates and assays, but no surveyed collar/interval tables, geological domains, density, CRS, or vertical datum.

Run from the project root after obtaining the GeoMet files described in `docs/provenance.md`. Python standard library only.
"""),
    code("""import csv, json, statistics
from pathlib import Path

ROOT = Path.cwd()
if not (ROOT / 'data' / 'raw').exists() and (ROOT.parent / 'data' / 'raw').exists():
    ROOT = ROOT.parent
DATA = ROOT / 'data' / 'raw' / 'drillholes.csv'
REPORTS = ROOT / 'reports'
with DATA.open(encoding='utf-8-sig', newline='') as f:
    rows = list(csv.DictReader(f))
print(f'Loaded {len(rows):,} rows and {len(rows[0]):,} columns')
"""),
    md("""## 1. Source structure and QA/QC

The original CSV is kept unchanged. Check dimensions, missing values, hole count, duplicate coordinates, and the observed coordinate ranges before modeling.
"""),
    code("""columns = list(rows[0])
missing = {c: sum(not r[c].strip() for r in rows) for c in columns}
xyz = [tuple(float(r[c]) for c in ('X', 'Y', 'Z')) for r in rows]
print('dimensions:', (len(rows), len(columns)))
print('columns:', columns)
print('nonzero missing:', {k:v for k,v in missing.items() if v})
print('holes:', len({r['HOLEID'] for r in rows}))
print('duplicate XYZ triplets:', len(xyz) - len(set(xyz)))
print('coordinate ranges:', {axis: (min(p[i] for p in xyz), max(p[i] for p in xyz))
      for i, axis in enumerate(('X', 'Y', 'Z'))})
"""),
    md("""## 2. Copper assay distribution

Cu is analyzed as an assay concentration in ppm, not a declared mineral resource grade or economic value.
"""),
    code("""cu = [float(r['Cu ppm']) for r in rows]
{'n': len(cu), 'min_ppm': min(cu), 'q1_ppm': statistics.quantiles(cu, n=4)[0],
 'median_ppm': statistics.median(cu), 'mean_ppm': statistics.fmean(cu),
 'q3_ppm': statistics.quantiles(cu, n=4)[2], 'max_ppm': max(cu)}
"""),
    md("""## 3. Check compositional-data assumptions

Some published workflows for this dataset apply log-ratio transforms to Cu, Au, Ag, and S. Before doing that, check whether the assay columns form a closed whole-rock composition and how many zeros require a documented treatment.
"""),
    code("""composition = json.loads((REPORTS / 'composition_readiness.json').read_text(encoding='utf-8'))
print('18-assay row-sum summary:', composition['all_18_assays']['row_sum_ppm'])
print('rows near 1,000,000 ppm closure:',
      composition['all_18_assays']['rows_with_sum_between_900000_and_1100000_ppm'])
print('Cu-Au-Ag-S subcomposition zeros:', composition['zero_rows_in_selected_subcomposition'])
print('decision:', composition['decision'])
"""),
    md("""## 4. Spatial subcomposition sensitivity

The next experiment treats Cu, Au, Ag, and S only as a four-part subcomposition. Zeros are replaced using training-fold minima under two hypothetical factors, CLR coordinates and log total are spatially interpolated, and Cu ppm is reconstructed. This is compared with raw Cu IDW and a training-mean baseline under both validation designs.
"""),
    code("""subcomp = json.loads((REPORTS / 'subcomposition_spatial_sensitivity.json').read_text(encoding='utf-8'))
for design, result in subcomp['designs'].items():
    print('\\n', design)
    for name, metrics in result['overall'].items():
        if isinstance(metrics, dict) and 'mae_ppm' in metrics:
            print(name, metrics)
"""),
    md("""## 5. Account for clustered sampling in the mean baseline

An equal-hole or equal-occupied-tile mean changes the weighting of the sampled population. Compare them as alternative reference baselines; they are not geostatistical estimates and do not establish the correct resource-support weighting.
"""),
    code("""declustered = json.loads((REPORTS / 'declustering_baselines.json').read_text(encoding='utf-8'))
for design, result in declustered['designs'].items():
    print('\\n', design)
    for name, metrics in result['overall'].items():
        print(name, metrics)
"""),
    md("""## 6. Grade-threshold screening

Treat top-quartile/top-decile flags as descriptive geochemical screens only. They are not economic ore/waste labels. The conditional chemistry model is only applicable where the companion assay covariates are known.
"""),
    code("""screen = json.loads((REPORTS / 'grade_threshold_classification.json').read_text(encoding='utf-8'))
print(screen['descriptive_grade_thresholds_ppm'])
for model, thresholds in screen['models'].items():
    print('\\n', model)
    for name, result in thresholds.items():
        print(name, {k:result[k] for k in ('precision','recall','f1','balanced_accuracy','average_precision')})
"""),
    md("""## 7. Top-ranked review workload

Ranking metrics show the fraction of high-Cu assays captured when only a chosen share of predicted samples is reviewed. These shares describe review workloads, not mine-planning cutoffs.
"""),
    code("""ranked = json.loads((REPORTS / 'top_rank_screening.json').read_text(encoding='utf-8'))
for model, bands in ranked['within_fold_percentile_ranking'].items():
    print('\\n', model)
    for band in bands:
        print(band['review_fraction_per_fold'], 'top-decile precision/recall:',
              band['top_decile']['precision'], band['top_decile']['recall'])
"""),
    md("""## 8. Spatial-fold stability

A paired tile bootstrap summarizes how model-versus-mean error differences vary across the 22 buffered spatial folds. These are descriptive stability intervals; folds may be dependent and the number of spatial groups is small.
"""),
    code("""stability = json.loads((REPORTS / 'spatial_fold_bootstrap.json').read_text(encoding='utf-8'))
for model, result in stability['results'].items():
    print(model, result['paired_difference_vs_training_mean']['mae'])
"""),
    md("""### Residual spatial structure

Moran's I is calculated for signed and absolute out-of-fold residuals over 3D k-nearest-neighbour graphs. Permutations are restricted within each outer fold to retain fold-specific residual distributions. This remains a descriptive diagnostic: the CRS and coordinate units are undocumented, and geology domains are unavailable.
"""),
    code("""residual_space = json.loads((REPORTS / 'residual_spatial_autocorrelation.json').read_text(encoding='utf-8'))
for model, by_k in residual_space['results'].items():
    print('\\n', model)
    for k, stats in by_k.items():
        print('k=', k,
              'signed I=', round(stats['signed_residual_ppm']['moran_i'], 3),
              'absolute I=', round(stats['absolute_residual_ppm']['moran_i'], 3),
              'permutation p=', stats['signed_residual_ppm']['two_sided_permutation_p'])
"""),
    md("""### Paired sample-level error comparison

Compare the conditional all-assay ridge with kriging on identical held-out samples. This performance case is conditional on all 17 secondary assays being available at target locations; the source-coordinate map is not a geographic map.
"""),
    code("""paired = json.loads((REPORTS / 'paired_error_comparison.json').read_text(encoding='utf-8'))
for group, result in paired['groups'].items():
    print(group, 'n=', result['n'],
          'kriging MAE=', round(result['kriging_mae_ppm']),
          'conditional ridge MAE=', round(result['ridge_mae_ppm']),
          'ridge wins=', round(result['ridge_lower_absolute_error_fraction'], 3))
"""),
    md("""## 9. Validation-design sensitivity

The fixed IDW baseline is also compared across alternative XY tile widths and buffers. This shows how conclusions change with the spatial holdout design; values are in undocumented local coordinate units.
"""),
    code("""design_sensitivity = json.loads((REPORTS / 'spatial_cv_design_sensitivity.json').read_text(encoding='utf-8'))
for r in design_sensitivity['designs']:
    print(r['tile_width_coordinate_units'], r['buffer_coordinate_units'],
          r['idw']['mae_ppm'], r['training_mean']['mae_ppm'],
          r['idw']['rmse_ppm'], r['training_mean']['rmse_ppm'])
"""),
    md("""## 10. Validation and support

Leave-one-hole-out tests transfer to withheld trajectories. Buffered XY-tile validation withholds occupied 300-unit tiles and excludes training points within a 50-unit horizontal buffer. The support cutoffs used for the exploratory grid were based on dense full-data sample spacing; test their relationship to buffered holdouts before treating them as meaningful.
"""),
    code("""baseline = json.loads((REPORTS / 'baseline_cv.json').read_text(encoding='utf-8'))
spatial = json.loads((REPORTS / 'spatial_eda.json').read_text(encoding='utf-8'))
support_cv = json.loads((REPORTS / 'grid_support_validation.json').read_text(encoding='utf-8'))
print('Leave-one-hole-out methods:', baseline['methods'])
print('Buffered spatial tiles:', spatial['xy_blocked_validation_buffered'])
print('Buffered holdout IDW:', support_cv['all_holdout_metrics']['idw'])
print('Buffered holdout mean:', support_cv['all_holdout_metrics']['training_mean'])
print('Current cutoff band cases:', support_cv['support_band_metrics'][1]['idw']['n'])
"""),
    md("""## 11. Exploratory interpolation grid

The grid uses 3D IDW (power 2, 12 neighbors), a sample XY convex hull, observed Z bounds, and explicit support cutoffs. A subsequent screen estimates a proxy surface from each hole's first sample row. This is not a surveyed terrain model, geological boundary, block model, or resource estimate.
"""),
    code("""grid_meta = json.loads((REPORTS / 'copper_idw_grid.json').read_text(encoding='utf-8'))
surface_meta = json.loads((REPORTS / 'surface_clip_audit.json').read_text(encoding='utf-8'))
print('Grid counts:', grid_meta['grid_counts'])
print('Surface clip:', surface_meta['output_nodes_below_proxy_surface'],
      'remaining of', surface_meta['input_grid_nodes'])
print('First sample above last sample Z:',
      surface_meta['first_sample_vs_last_sample_z_order'])
"""),
    md("""## 12. Assay-availability sensitivity

The source contains no missing chemical assays. These explicit stress scenarios mask specified assays only at held-out locations, while training assays remain complete; they are not measurements of natural missingness. Compare scenarios before selecting a prediction workflow.
"""),
    code("""dropout = json.loads((REPORTS / 'geochem_assay_dropout.json').read_text(encoding='utf-8'))
for scenario, result in dropout['results'].items():
    print(scenario,
          'MAE/RMSE:', round(result['overall']['mae_ppm']), round(result['overall']['rmse_ppm']),
          'top-quartile MAE:', round(result['top_quartile']['mae_ppm']))
print('Scenario caveats:', dropout['caveats'])
"""),
    md("""## 13. Companion table linkage and decision

The comminution and flotation tables have no shared sample/interval key with the drillhole assay table. Spatially close coordinates are candidates only; do not transfer lab response values without an authoritative crosswalk. The small same-row LCT chemistry pilot did not beat its mean baseline.

The project demonstrates reproducible QA/QC, spatial validation, geostatistical experimentation, and explicit uncertainty/data limitations. The next useful inputs are collar and downhole surveys, interval/compositing definitions, lithology or mineralization domains, coordinate reference and vertical datum, density and assay QA metadata, and a companion-test sample crosswalk. Do not report tonnage, classification, or economic value from the current tables.

Reproduction commands for all analyses are listed in `README.md`. Reports are in `reports/`; exploratory CSVs are in `outputs/`; raw GeoMet source files remain unchanged.
"""),
]

notebook = {
    "cells": cells,
    "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                 "language_info": {"name": "python", "version": "3"}},
    "nbformat": 4,
    "nbformat_minor": 5,
}
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
print(OUT)
