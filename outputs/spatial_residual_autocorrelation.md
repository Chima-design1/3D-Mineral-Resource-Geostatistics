# Spatial structure in out-of-fold residuals

**Target:** Cu ppm. **Predictions:** buffered 300-unit XY tile / 50-unit buffer out-of-fold predictions across 2,000 samples. **Neighbour graph:** 3D Euclidean k-nearest neighbours, tested at k=4, 8, and 12. **Test:** 499 label permutations restricted within each outer fold, preserving fold-specific residual distributions. The coordinate CRS and units are undocumented, so results are descriptive.

## Results

| Model | Signed-residual Moran's I, range across k | Absolute-residual Moran's I, range across k | Permutation result |
|---|---:|---:|---|
| Nested ordinary kriging | 0.165–0.237 | 0.125–0.142 | Exceeded the 95% within-fold permutation range at k=4, 8, 12; empirical two-sided p=0.002 (the minimum attainable with 499 permutations) |
| Conditional all-17 ridge, alpha 0.1–100 | 0.141–0.165 | 0.097–0.130 | Same |
| Conditional all-17 ridge, alpha 100–10,000 | 0.090–0.114 | 0.094–0.119 | Same |
| Conditional S/Ag/Au ridge | 0.076–0.085 | 0.101–0.112 | Same |

Positive signed-residual Moran's I means nearby samples tend to share the sign of their prediction error. Positive absolute-residual Moran's I indicates nearby error magnitudes also tend to cluster. The pattern persists across the three neighbour counts and after preserving each outer fold's residual distribution. Conditional chemistry models reduce overall error, but they do not remove all local structure from their errors.

## Interpretation limits

This is a diagnostic of out-of-fold prediction residuals, not a formal population inference. The empirical p-value is resolution-limited by 499 permutations. The nearest-neighbour graph assumes comparable coordinate scaling across X, Y, and Z; the source provides sample centroids but no coordinate reference, vertical datum, or geological domains. The result cannot identify whether the remaining structure arises from geological continuity, domain boundaries, sampling, measurement, or fold design. Do not use it as an uncertainty model or as justification for a residual-correction surface.

Reproduce with `python scripts/residual_spatial_autocorrelation.py`. Full values, graph details, permutations, and caveats are in `reports/residual_spatial_autocorrelation.json`.
