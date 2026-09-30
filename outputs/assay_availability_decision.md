# Assay-availability sensitivity

**Project:** 3D Geomodeling and Geostatistical Mineral Resource Estimation Using Public Drillhole Data  
**Purpose:** determine how conditional Cu prediction changes when selected co-located assays are unavailable at held-out locations.

## Design

All scenarios use the same 22 buffered XY holdout tiles and 2,000 out-of-fold sample predictions. The source table is complete case; every missing-assay scenario is therefore hypothetical. Ridge models are trained on complete training-fold chemistry. For a specified missing assay at a held-out location, its log-transformed standardized value is set to the training-fold mean (zero after standardization). Alpha is selected in the inner spatial splits under that same target-location masking scenario. Nearby training locations retain their observed Cu values in every scenario.

The spatial-only comparator is nested ordinary kriging and the chemistry comparators are the nested all-17 assay and S/Ag/Au ridge variants. These estimate different information-availability cases and should not be conflated.

## Out-of-fold results

| Information available at Cu target | MAE (ppm) | RMSE (ppm) | Top-quartile Cu MAE (ppm) | Interpretation |
|---|---:|---:|---:|---|
| No target assays; neighboring training Cu available (nested kriging) | 5,200 | 7,624 | 10,407 | Spatial-only baseline; still interpolates within the sampled area |
| Other 17 assays available | 1,898 | 3,248 | 3,787 | Best overall conditional scenario in this comparison |
| S, Ag, Au available | 2,075 | 3,593 | 4,294 | Three-assay conditional scenario |
| All 17 except S | 1,822 | 3,143 | 3,699 | No degradation in this hypothetical single-assay mask |
| All 17 except Ag | 2,734 | 4,449 | 5,571 | Largest single-assay effect; performance worsens substantially |
| All 17 except Au | 2,219 | 3,861 | 4,482 | Moderate degradation |
| All 17 except S, Ag, and Au | 3,866 | 5,906 | 7,650 | Material degradation; remains conditional on the other 14 assays |

Top-quartile membership is based on the descriptive 8,900 ppm Cu cutoff. It is not an economic or resource-class threshold. The no-target-assay kriging comparator is from the nested buffered-kriging predictions; its top-quartile MAE is calculated on the same held-out observations and cutoff.

## Decision and limitation

If the practical workflow can guarantee Cu plus co-located secondary assays at the prediction site, the ridge results support further conditional prediction experiments. Ag availability appears particularly consequential in the tested masks. Where no target chemistry is available, the relevant comparator is spatial-only prediction, which has substantially higher error here. These figures do not establish that any assay is causally or universally necessary: the source has no naturally missing assays, and masks were imposed only at validation locations.

This is not resource estimation. The source lacks documented CRS/units, assay detection-limit handling, sample intervals, geological domains, density, and verified links from these drillhole samples to the companion metallurgical tables.

## Next concrete step

Request or locate the GeoMet sample crosswalk and assay metadata, plus collar/survey, interval, lithology/domain, and coordinate-reference documentation. Until those records are available, do not transfer flotation/comminution responses to drillholes or interpret interpolated concentrations as block grades. With those inputs, the next defensible modeling stage is domain-aware compositing followed by spatially validated estimation on a documented support.

Reproduction: `python scripts/geochem_assay_dropout.py` and `python scripts/assay_availability_scenarios.py`. Source results: `reports/geochem_assay_dropout.json` and `reports/assay_availability_scenarios.json`.
