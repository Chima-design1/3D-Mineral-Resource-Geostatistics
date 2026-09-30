# Public drillhole dataset screening

**Assessment date:** 2026-09-30  
**Purpose:** identify data that could extend the GeoMet sample-level analysis into a defensible 3D drillhole workflow. This screening does not claim that any dataset supports a compliant mineral resource estimate.

## Requirements from the current modeling plan

The existing [modeling plan](modeling_plan.md) identifies the main blockers in GeoMet v4: no collar/survey tables, interval support, geological domains, density, documented CRS/vertical datum, or assay QA metadata. A candidate should provide linked collars, downhole surveys, from/to sample assays, usable geology, and enough spatial/analytical metadata to audit the observations. Dataset availability and license must be independently checked before incorporation.

## Candidate comparison

| Candidate | Verified strengths | Gaps or qualification | Decision |
|---|---|---|---|
| **Northern Territory Geological Survey (NTGS), DIP001, Northern Territory geochemical and drilling datasets** | The current GEMIS record lists drillhole collars, downhole surveys, drillhole sample assays, raw geochemistry including repeats, and associated laboratory metadata. Its stated access condition is CC BY 4.0. September 2025 release notes report 313,533 collar records, 149,659 survey records, 2,210,334 drillhole assay records, 1,923 raw assay CSVs, and 43,371 lab-metadata records. The April 2026 update reports additions of 7,679 holes and 444,810 analyses. | The current archive is 1.04 GB. Its actual archive and CSV headers could not be retrieved in this environment, so exact join keys, interval-depth fields, CRS fields, units, and per-record completeness are not yet verified. The record does not list a system-wide lithology interval table. NTGS metadata describes collar positional accuracy ranging from map-derived tens to hundreds of metres in some records and notes that original source/lab methods must be checked. The territory-wide collection will need filtering to a coherent project/area and CRS/datum before modeling. | **Best real-data candidate; do not integrate yet.** Proceed only after obtaining an inspectable subset and verifying actual headers, joins, intervals, geology, and coordinate metadata. |
| **NTGS Geochemistry – Drillhole Samples and Drillhole Collars open data layers** | Official metadata describes open-file drillhole sample geochemistry, publishes collar data through CSV/SHP/TAB, and states CC BY 4.0. The collar metadata documents MGA94/map-grid coordinates and decimal latitude/longitude, while noting mixed positional accuracy. | The public drillhole-sample layer directs users to DIP001 for data; NTGS says assay/lab-method detail and downhole survey are not all exposed in the flat-file web layers. The collar layer alone does not supply assays, intervals, or lithology. Exact field schema has not been locally verified. | **Useful access route for collar reconnaissance, insufficient alone.** |
| **Ontario Drill Hole Database (ODHD)** | The Ontario Data Catalogue says it contains more than 164,000 holes and records location, company/hole number, orientation, depth, overburden depth, and whether assay results exist for selected commodities. | The catalogue description does not establish availability of linked from/to assay intervals, downhole survey stations, geology logs, analytical QA/QC, or a project-specific consistent CRS/schema in this screening. | **Not sufficient on the verified description.** |
| **GeoMet synthetic porphyry datasets** | The public repository documents synthetic geometallurgical datasets, including one described as a full dataset, with repository documentation and CC BY-NC-SA 4.0 terms. | Synthetic data cannot validate geological realism or support claims about a real deposit or resource. Dataset schemas and individual content have not been inspected here. | **Potentially useful only as a separate software/workflow test fixture, not as replacement real evidence.** |

## NTGS DIP001 verification details

### Provenance and documentation

- Publisher: Northern Territory Geological Survey, Northern Territory Government.
- Current record: [GEMIS DIP001 record](https://geoscience.nt.gov.au/gemis/ntgsjspui/handle/1/94622); record metadata identifies the edition as 31 August 2026 and lists a 1.04 GB ZIP.
- Record contents: collar, downhole-survey, drillhole sample-assay, maximum-assay, surface geochemistry, raw assay, and lab-metadata products.
- Summary of the package and historical record counts: [NTGS DIP001 September 2025 release notes](https://geoscience.nt.gov.au/gemis/ntgsjspui/bitstream/1/81743/4/DIP001%20-%20Release%20Notes.pdf).
- Later additions: [NTGS April 2026 DIP001 update notice](https://resourcingtheterritory.nt.gov.au/news-and-events/news/2026/new-data-dip001-update).
- Reuse condition: the current GEMIS record states **Creative Commons Attribution 4.0 International (CC BY 4.0)**; retain NTGS attribution and notices.
- Collar metadata: [NTGS drillhole-collar metadata](https://www.ntlis.nt.gov.au/metadata/export_data?type=html&metadata_id=10E907B8F993795DE050CD9B21447B20). It describes positions sourced from map capture and GPS, MGA94 or AGD66-era map grids, decimal lat/long availability, and estimated positional accuracy of roughly 20–150 m for the documented layer. Do not assume that the entire territory-wide package has one CRS or uniform accuracy.
- Drillhole geochemistry metadata: [NTGS drillhole-sample metadata](https://www.ntlis.nt.gov.au/metadata/export_data?type=html&metadata_id=6740424DE611B24CE050CD9B2144744A). It distinguishes open-file/company and NTGS data, describes legacy digitization, and says laboratory methods for flat-file layers may require checking original source reports. It also describes QA/QC practices for NTGS-collected samples; those practices must not be generalized to all company records.
- NTGS describes the package's current file groups as CSV, ESRI, and MapInfo products and says raw assay repeats and lab metadata are included. These are source-record/release-note descriptions; the project has not yet verified their exact local file headers.

### Requirement-by-requirement status

| Requirement | Evidence available from the official documentation | Status before using data |
|---|---|---|
| Drillhole collars | Explicit NTGS collar table; large historical count; separate public collar layer | Promising; retrieve subset and inspect keys, coordinates, duplicate IDs, and project coverage |
| Downhole surveys | Explicit NTGS survey table; release notes report a large survey table | Promising; verify stations, azimuth convention, dip sign, units, collar joins, and survey coverage |
| Sample intervals | Product is described as drillhole sample assays and drill samples, but actual current CSV headers/interval fields are not verified locally | Unconfirmed; require actual FROM/TO or documented downhole depth/support fields |
| Assay values | Drillhole assay and raw geochemistry tables are explicitly included | Promising; inspect analyte, result, unit, qualifier, detection-limit, and repeat fields |
| Geology/lithology | No dedicated drillhole lithology-interval table appears in the current DIP001 product listing reviewed | Missing/unconfirmed; obtain project-specific open report logs or select an additional documented table before domain modeling |
| CRS / vertical datum | Metadata documents MGA94/AGD66-era map-grid sources and lat/long availability for collar layers | Mixed/record-specific; identify datum, projection zone, coordinate accuracy, and elevation datum for the selected subset |
| Assay QA/QC | Raw repeat assays and lab metadata are included; NTGS describes QA/QC for NTGS-collected samples | Partial; separate company vs NTGS provenance, inspect standards/blanks/duplicates and method/detection-limit coverage |
| Density / block support | No density or block-support table is identified in the DIP001 product list reviewed | Not supplied by the verified product description; would require an additional documented source |
| License | Current GEMIS record says CC BY 4.0 | Verified at package-record level; preserve attribution and check notices within the downloaded package |
| Local schema and completeness | Archive advertised at 1.04 GB; no files have been downloaded or locally inspected | Not verified; do not model or publish claims based only on catalogue descriptions |

## Screening result and next data gate

DIP001 is a substantially stronger **candidate** than the current GeoMet table because official documentation explicitly lists collar, downhole-survey, drillhole-sample geochemistry, raw repeats, laboratory metadata, and a reuse license. It is not yet a verified modeling dataset for this project: the package's actual headers and joins have not been inspected, the product listing does not establish sample interval fields or a lithology-interval table, and its large territory-wide contents may mix coordinate systems, survey quality, source practices, and geological settings.

The next data gate is to obtain an inspectable DIP001 subset for one named open-file project or coherent area, together with its release dictionary and source report. Before adding any rows to the existing pipeline:

1. Inventory actual files and columns; preserve source checksums and the exact edition.
2. Verify that collar, survey, assay, and any lithology tables join on documented identifiers.
3. Confirm actual sample FROM/TO or downhole depth fields, assay units/qualifiers, and sample support.
4. Identify coordinate datum, projection/zone, vertical datum, positional accuracy, azimuth reference, and dip convention for every selected hole.
5. Audit duplicate identifiers, survey coverage, interval overlap/gaps, missing assays, repeated assays, and QA/QC sample types.
6. Confirm a usable lithology/domain source before any domain-based interpolation or resource-style block modeling.
7. Keep the NTGS data under its source attribution/license and outside this repository until the release audit confirms what can be redistributed.

If this gate cannot be passed, the scientifically defensible stage remains the existing GeoMet sample-level study: test validation and geostatistical assumptions, report support and uncertainty limitations, and stop short of block/resource estimation. A synthetic dataset may be used only to test data plumbing and desurvey/compositing code, with all results clearly labeled synthetic and separate from empirical findings.
