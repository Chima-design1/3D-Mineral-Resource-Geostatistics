# Spatial validation sensitivity to grid origin

## Question

Does the fixed buffered XY-tile validation result depend on where the arbitrary tile grid begins? The earlier design sensitivity varied tile width and buffer but retained a single grid origin. This analysis tests nine origins while holding tile width (300 source-coordinate units), buffer (50 units), and prediction methods fixed.

## Method

For each X/Y origin shift in {0, 100, 200}, the 2,000 Cu samples were assigned to XY tiles. Each occupied tile was held out in turn; training samples within the 50-unit horizontal buffer were excluded. The test predictions were 3D Euclidean IDW (power 2, 12 neighbors) and a training-sample mean recomputed for each fold. All 2,000 samples were evaluated in every origin design. The shifts are relative sensitivity settings; the source CSV does not document coordinate units or CRS.

## Results

| X/Y origin shift | Folds | IDW MAE (ppm) | Training-mean MAE (ppm) | Lower MAE |
|---|---:|---:|---:|---|
| 0 / 0 | 22 | 5,478.6 | 5,353.4 | Mean |
| 0 / 100 | 23 | 5,376.3 | 5,358.3 | Mean |
| 0 / 200 | 22 | 5,762.8 | 5,370.9 | Mean |
| 100 / 0 | 22 | 5,514.1 | 5,350.5 | Mean |
| 100 / 100 | 23 | 5,522.6 | 5,348.8 | Mean |
| 100 / 200 | 25 | 5,674.0 | 5,370.9 | Mean |
| 200 / 0 | 23 | 5,325.7 | 5,352.1 | IDW |
| 200 / 100 | 23 | 5,247.6 | 5,343.6 | IDW |
| 200 / 200 | 23 | 5,526.5 | 5,372.2 | Mean |

The IDW MAE ranges from 5,247.6 to 5,762.8 ppm across origins; the training-mean MAE ranges from 5,343.6 to 5,372.2 ppm. IDW beats the mean on MAE in only two of nine origins. This adds evidence that the apparent performance of a spatial interpolation baseline depends on an arbitrary partition choice. It does not select a best origin: the grid should not be tuned to minimize validation error, and the nine partitions are not independent replications.

The result supports reporting validation-design sensitivity and simple baselines together. It does not validate any unsampled grid node, estimate block support, or support a mineral-resource claim.

## Reproduction and files

Run `python scripts/spatial_origin_sensitivity.py` or include it in `python scripts/run_all.py`. The script writes `reports/spatial_origin_sensitivity.json` and [the comparison figure](spatial_origin_sensitivity.svg). It uses only the source CSV; the derived JSON and SVG are generated locally and excluded from publication.
