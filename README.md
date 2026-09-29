# Analysis Code for CADe in Real-World Gastroscopy

This repository contains the statistical analysis scripts and executable workflow used in this study.

## Files

- `run_workflow.py`: Runs the analysis scripts in sequence and saves their outputs.
- `config.example.json`: Template for configuring input file paths.
- `DATA_DICTIONARY.md`: Definitions of the analysis variables and auxiliary data fields.
- `requirements.txt`: Python package dependencies.

## Input Data

The workflow requires the following authorized, de-identified data files:

| Configuration key | File |
| --- | --- |
| `cohort1` | `cohort_1_analysis.csv` |
| `cohort2` | `cohort_2_analysis.csv` |
| `first_study` | `first_examination_mapping.xlsx` |
| `screening_result` | `lesion_documentation.xlsx` |
| `cancer` | `gastric_neoplastic_lesions.xlsx` |
| `polyp` | `gastric_polypoid_lesions.xlsx` |

## Environment Setup

Python 3.10 or later is recommended.

```bash
python -m pip install -r requirements.txt
```

## Configuration and Execution

Copy the configuration template:

```powershell
Copy-Item config.example.json config.json
```

Enter the authorized data-file paths in `config.json`, and then run:

```bash
python run_workflow.py --config config.json --cohort both --output-dir results
```

Analysis outputs are saved in `results/cohort1/` and `results/cohort2/`. Runtime copies of the analysis scripts are saved in `results/staged_scripts/`.

## Data Availability

Patient-level data are accessible only subject to institutional approval and applicable data-use conditions.
