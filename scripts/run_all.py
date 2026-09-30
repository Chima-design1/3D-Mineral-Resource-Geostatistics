"""Reproduce the complete GeoMet analysis in dependency order.

Run from any working directory with Python 3.11+: python scripts/run_all.py
The source CSVs must be present under data/raw/.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = (
    "initial_qaqc.py",
    "spatial_eda.py",
    "baseline_cv.py",
    "fold_kriging.py",
    "fold_diagnostics.py",
    "nested_kriging_sensitivity.py",
    "residual_diagnostics.py",
    "geochem_covariate_cv.py",
    "geochem_feature_sensitivity.py",
    "geochem_assay_dropout.py",
    "assay_availability_scenarios.py",
    "companion_table_audit.py",
    "companion_link_sensitivity.py",
    "sample_influence.py",
    "high_grade_context.py",
    "spatial_support_audit.py",
    "copper_idw_grid.py",
    "grid_support_sensitivity.py",
    "grid_support_validation.py",
    "support_error_calibration.py",
    "composition_readiness.py",
    "subcomposition_spatial_sensitivity.py",
    "declustering_baselines.py",
    "grade_threshold_classification.py",
    "top_rank_screening.py",
    "spatial_fold_bootstrap.py",
    "spatial_cv_design_sensitivity.py",
    "spatial_origin_sensitivity.py",
    "residual_spatial_autocorrelation.py",
    "clip_copper_grid_to_surface.py",
    "flotation_response_cv.py",
    "paired_error_comparison_map.py",
    "create_validation_figure.py",
    "create_portfolio_figure.py",
    "create_portfolio_notebook.py",
)


def main() -> None:
    for script in SCRIPTS:
        path = ROOT / "scripts" / script
        print(f"Running {script}", flush=True)
        subprocess.run([sys.executable, str(path)], cwd=ROOT, check=True)
    print(f"Completed {len(SCRIPTS)} analysis steps.", flush=True)


if __name__ == "__main__":
    main()
