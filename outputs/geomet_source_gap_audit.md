# GeoMet source-gap audit

**Checked:** 30 September 2026  
**Purpose:** determine whether the public GeoMet files support progression from sample-point interpolation to drillhole desurveying, interval compositing, and domain-aware estimation.

## Verified package

The current GeoMet version 4 record [10.5281/zenodo.7051975](https://zenodo.org/records/7051975) lists exactly three files:

| File | Record size | MD5 shown by Zenodo |
|---|---:|---|
| `drillholes.csv` | 362.2 kB | `bf48a6b1135113b6e6cf7b32f95315ad` |
| `comminution.csv` | 25.2 kB | `1a33df8ba77f5d49c281ba3fff16b20b` |
| `flotation.csv` | 17.5 kB | `f2e90da6bfa81de1261177ee85a91570` |

The repository copy is the same checksum-verified v4 source already analyzed in this project. The record does not list separate collar, downhole survey, assay interval, lithology/domain, density, or sample-crosswalk files.

The GeoMet worked-example chapter describes the drillhole table as processed sample chemistry with Cartesian X/Y/Z **centroids** of sample cylinders. It separately notes that typical drillhole projects have `COLLAR`, `SURVEY`, and `INTERVAL` tables for desurveying and compositing. That supports interpreting the provided coordinates as sample locations, not documented collar coordinates or a recoverable survey trace. [GeoMet worked example, Data section](https://juliaearth.github.io/geospatial-data-science-with-julia/12-mining.html#121-data)

## Related paper data record

The associated 2022 paper's data-availability section points to [Zenodo DOI 10.5281/zenodo.6336137](https://doi.org/10.5281/zenodo.6336137) for obfuscated versions of the comminution, flotation, and drillhole tables, and states that obfuscation does not affect the presented results. The paper describes a generated predictions table in its methods, but does not identify that as a fourth public source table in the data-availability statement. [Paper manuscript and data-availability section](https://doi.org/10.1007/s11004-022-10013-1)

The DOI 6336137 record could not be opened or independently inspected through the available source viewer during this audit; a direct read-only request also failed at the connection layer. Its file list and columns therefore remain **unverified**; it must not be treated as a source of collar, survey, interval, domain, or linkage fields until inspected directly. The paper's phrasing is evidence that a second copy of the three tables exists, not proof that it contains the missing drillhole database components.

## Project decision

The current GeoMet files support geospatial QA/QC of sample-centroid points, spatial cross-validation, and exploratory sample-point interpolation. They do not, on the evidence currently inspectable, support reconstructing true drillhole trajectories, compositing interval assays, geological domaining, or a defensible block/resource estimate. Existing outputs should remain labelled exploratory or conditional prediction.

### Next acquisition action

1. Inspect the 6336137 Zenodo record and compare its actual files and schemas with the v4 package.
2. If it contains only the three named tables, request from the dataset authors/source provider the collar/survey and interval tables, geological logs/domains, coordinate reference and vertical datum, assay QA/detection-limit metadata, and authoritative links from metallurgical samples to drillhole samples.
3. If those inputs are not obtainable, continue this portfolio as a spatial-geostatistics and data-limitations study using GeoMet; do not present it as a full mineral-resource estimation exercise.

## Source evidence

- [Zenodo GeoMet v4 record](https://zenodo.org/records/7051975)
- [Hoffimann et al. (2022), paper DOI](https://doi.org/10.1007/s11004-022-10013-1)
- [Paper's cited obfuscated-table record](https://doi.org/10.5281/zenodo.6336137)
- [GeoMet worked example](https://juliaearth.github.io/geospatial-data-science-with-julia/12-mining.html#121-data)
