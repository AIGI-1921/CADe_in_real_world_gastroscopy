# Data Dictionary

This dictionary includes only input variables directly referenced by the shared analysis scripts and variables derived within those scripts. The cohort CSV files contain additional fields that are not required by this workflow and are therefore not listed. Binary variables are coded `0 = no` and `1 = yes` unless otherwise stated.

## Analysis cohort files

Files: `cohort_1_analysis.csv` and `cohort_2_analysis.csv`

### Identifiers, exposure, and timing

| Variable | Type | Definition / coding |
| --- | --- | --- |
| `Barcode` | string | De-identified examination identifier used to link authorized auxiliary files. |
| `cade` | category | Study exposure: `A` = CADe-assisted; `M` = control without CADe assistance. |
| `报告医师` | string | Reporting endoscopist identifier; used only for physician-related sensitivity analyses. |
| `platform_code` | category | Endoscopy platform identifier: `1` or `2`. |
| `platform_group` | category | Labeled platform category used as an adjustment covariate. |
| `report_date` | date | Calendar date of the examination. |
| `platform_session_cluster` | string | Platform-session identifier combining platform, calendar date, and half-day session; used as the clustering variable in cluster-robust models. |
| `study_phase` | category | Prespecified broad study phase (`Phase 1` or `Phase 2`). |
| `fatigue_time` | boolean | Indicator for an examination in the post hoc fatigue-time period. |

### Baseline covariates

| Variable | Type | Definition / coding |
| --- | --- | --- |
| `age_years` | numeric | Age at examination, in years. |
| `sex_group` | category | `Male` or `Female`. |
| `indication` | category | Examination indication: `diagnosis`, `screening`, or `surveillance`. |
| `patient_source_group` | category | `Gastroenterology` or `Non-gastroenterology`. |
| `hp_detail` | category | Detailed HP-RUT result used for descriptive reporting. |
| `hp_model` | category | Model category: `Negative_or_indeterminate`, `Positive`, or `Missing`. |
| `endoscopist_volume_group` | category | Prior procedure-volume category: `<1000` or `>=1000` examinations. |

### Pathology-confirmed yield outcomes

| Variable | Type | Definition |
| --- | --- | --- |
| `outcome_neoplasm` | binary | At least 1 pathology-confirmed gastric neoplastic lesion in the examination. |
| `outcome_atrophy_im` | binary | At least 1 pathology-confirmed finding of intestinal metaplasia and/or gastric atrophy. |
| `outcome_polyp` | binary | At least 1 pathology-confirmed gastric polyp. |
| `outcome_gastritis` | binary | At least 1 pathology-confirmed gastritis finding. |

### Biopsy-rate outcomes

| Variable | Type | Definition |
| --- | --- | --- |
| `biopsy_any` | binary | At least 1 gastric biopsy was performed. |
| `biopsy_rate_cardia_fundus` | binary | At least 1 biopsy mapped to the cardia-fundus region. |
| `biopsy_rate_body` | binary | At least 1 biopsy mapped to the gastric body. |
| `biopsy_rate_angulus` | binary | At least 1 biopsy mapped to the gastric angulus. |
| `biopsy_rate_antrum_pylorus` | binary | At least 1 biopsy mapped to the antrum-pylorus region. |
| `biopsy_rate_remnant` | binary | At least 1 biopsy mapped to the gastric remnant or thoracic stomach. |

### Biopsy-count and procedure-time outcomes

| Variable | Type | Definition |
| --- | --- | --- |
| `biopsy_count_any` | numeric | Total number of mapped gastric biopsy specimens/lesions per examination. |
| `biopsy_count_cardia_fundus` | numeric | Number mapped to the cardia-fundus region. |
| `biopsy_count_body` | numeric | Number mapped to the gastric body. |
| `biopsy_count_angulus` | numeric | Number mapped to the gastric angulus. |
| `biopsy_count_antrum_pylorus` | numeric | Number mapped to the antrum-pylorus region. |
| `biopsy_count_remnant` | numeric | Number mapped to the gastric remnant or thoracic stomach. |
| `turnover_time` | numeric | Examination turnover time in minutes. |
| `turnover_valid` | binary | Indicator that `turnover_time` met the prespecified validity rules. |

### Derived variables created by the scripts

| Variable | Definition |
| --- | --- |
| `g` | Numeric CADe indicator: `1` for `cade == "A"`, otherwise `0`. |
| `cluster` | Copy of `platform_session_cluster`, representing a half-day session on a given platform and calendar date. |
| `age_group` | `<40` or `>=40` years for subgroup analysis. |
| `fatigue_subgroup` | `Yes` or `No`, derived from `fatigue_time`. |
| `platform_subgroup` | `Platform 1` or `Platform 2`, derived from `platform_code`. |
| `time_block` | One of 10 contiguous blocks formed from chronologically ordered active examination dates. |
| `phys_group` | Each of the 10 highest-volume reporting endoscopists as a separate category; all others pooled as `Others`. |

## Auxiliary file: first-examination mapping

File: `first_examination_mapping.xlsx`; sheet: `Sheet1`

| Source column | Definition |
| --- | --- |
| `Barcode`, `Barcode.1` | Examination identifiers from the 2 source blocks. |
| `病人号`, `病人号.1` | De-identified patient identifiers used to identify repeated examinations. |
| `报告时间`, `报告时间.1` | Examination/report timestamps used to select the first eligible examination. |

## Auxiliary file: lesion documentation

File: `lesion_documentation.xlsx`

Every sheet listed below must contain a `Barcode` column. Membership in a sheet indicates that the corresponding lesion category was documented in the routine endoscopy report.

| Sheet | Derived category |
| --- | --- |
| `S3_screen_by_type` | Any documented gastric lesion |
| `benmen_weidi` | Lesion involving cardia-fundus |
| `weiti` | Lesion involving gastric body |
| `weijiao` | Lesion involving gastric angulus |
| `weidou_youmen` | Lesion involving antrum-pylorus |
| `canwei` | Lesion involving gastric remnant or thoracic stomach |
| `kuiyang_aoxian` | Ulcerative-depressed lesion |
| `milan_fahong` | Erosive-erythematous lesion |
| `xirou_longqi` | Polypoid-elevated lesion |
| `qita` | Other documented lesion |
| `weiti_baxiangbingbian` | Target lesion involving the gastric body |

## Auxiliary file: gastric neoplastic lesions

File: `gastric_neoplastic_lesions.xlsx`

| Source column | Definition / transformation |
| --- | --- |
| `Barcode` | Links the lesion record to the cohort exposure. |
| `病变1病理类型` | Principal pathological subtype, mapped to adenocarcinoma, HGIN, LGIN, adenoma, NET, lymphoma, GIST, or metastasis. |
| `病变1受累部位` | Narrative anatomical involvement; individual sites are identified by keyword matching. |
| `病变1大小` | Lesion size, categorized as `<=5`, `>5 to <=10`, `>10 to <=20`, or `>20` mm; `大` is treated as `>20 mm`. |
| `病变1形态` | Elevated (`隆起`) or flat (`平坦`) morphology. |
| `病变1是否有溃疡` | Surface ulceration (`有` or `无`). |

## Auxiliary file: gastric polypoid lesions

File: `gastric_polypoid_lesions.xlsx`; sheet: `Sheet1`

Up to 3 lesions are stored per examination. The first lesion uses the full source labels shown below; lesions 2 and 3 use the corresponding abbreviated labels `息肉病变{n}病理`, `息肉病变位置{n}`, `息肉病变大小{n}`, and `息肉病变形态{n}`.

| Source field | Coding |
| --- | --- |
| `Barcode` | Links the lesion record to the cohort exposure. |
| Pathological subtype | `1` fundic gland; `2` hyperplastic; `3` hamartomatous; `4` inflammatory; `5` Peutz-Jeghers/hamartomatous. |
| Site | `1` cardia-fundus; `2` body; `3` angulus; `4` antrum-pylorus; `5` gastric remnant. |
| Size | Millimeters; categorized as `<=5`, `>5 to <=10`, `>10 to <=20`, or `>20` mm. |
| Morphology | `1` non-flat; `2` flat; `3` multifocal non-flat; `4` multifocal flat. |

## Missing data handling

- Model functions construct complete-case datasets from the outcome, exposure, cluster, and covariates required for that model.
- `hp_model == "Missing"` is retained as a category in the main adjusted models and excluded from the HP-RUT interaction comparison.
- Turnover-time models are restricted to `turnover_valid == 1`.
- Missing lesion characteristics are retained as `NA` where the corresponding lesion-level script specifies an explicit missing category.
