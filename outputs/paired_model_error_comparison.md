# Paired spatial and chemistry model error comparison

The comparison uses 2,000 identical buffered XY holdout rows. Each sample has one nested ordinary-kriging Cu prediction and three assay-assisted ridge predictions, all out of fold. The map shows the per-row difference `|kriging error| − |all-17 ridge error|`: green means the conditional ridge had lower absolute error at that held-out sample; orange means kriging had lower absolute error.

![Paired absolute-error comparison map in the source XY coordinate frame](paired_error_comparison.svg)

## Summary by observed Cu range

| Observed Cu group | n | Kriging MAE | All-17 ridge alpha 0.1–100 MAE | All-17 ridge alpha 100–10,000 MAE | S/Ag/Au ridge MAE | Low-alpha ridge had lower absolute error |
|---|---:|---:|---:|---:|---:|---:|
| Below Q75 (<8,900 ppm) | 1,498 | 3,455 | 1,260 | 1,265 | 1,332 | 78.7% |
| Q75 to Q90 (8,900–18,100 ppm) | 300 | 5,458 | 2,540 | 2,182 | 2,118 | 79.3% |
| At or above Q90 (≥18,100 ppm) | 202 | 17,758 | 4,658 | 6,171 | 7,527 | 100.0% |
| All held-out samples | 2,000 | 5,200 | 1,795 | 1,898 | 2,075 | 81.0% |

On this split, the low-alpha all-17 ridge reduced average absolute error relative to kriging by 3,405 ppm across all rows and 13,100 ppm among top-decile Cu samples. It had lower absolute error than kriging for 81.0% of all samples and all 202 top-decile samples.

## Interpretation

This is a paired sample-level result under a specific information scenario: all 17 non-Cu assays are measured and available at each held-out target. It does not establish performance at unassayed locations or validate resource estimates. The XY display uses undocumented source coordinates without CRS information; patterns must not be interpreted as geographic priorities or geological domains. Full row-level metrics and fold summaries are in `reports/paired_error_comparison.json`; reproduce the figure and report data with `python scripts/paired_error_comparison_map.py`.
