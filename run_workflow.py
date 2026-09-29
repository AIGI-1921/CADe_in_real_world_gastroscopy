import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


CORE_SCRIPTS = [
    "baseline_report.py",
    "outcome_report.py",
    "covariate_report.py",
    "cluster_report.py",
    "risk_diff_report.py",
    "fdr_report.py",
    "subgroup_ana.py",
    "subgroup_others.py",
]

COHORT1_ONLY_SCRIPTS = [
    "first_exam_sensitivity.py",
    "biopsied_sensitivity.py",
    "physician_detail.py",
    "time_block_detail.py",
    "neoplasm_report.py",
    "lesion_documentation.py",
    "polyp_report.py",
    "body_biopsy_subgroup.py",
]

PATH_PLACEHOLDERS = {
    'r"__COHORT_FILE__"': "cohort",
    'r"__FIRST_EXAMINATION_MAPPING_FILE__"': "first_study",
    'r"__LESION_DOCUMENTATION_FILE__"': "screening_result",
    'r"__GASTRIC_NEOPLASTIC_LESIONS_FILE__"': "cancer",
    'r"__GASTRIC_POLYPOID_LESIONS_FILE__"': "polyp",
}

PANDAS_DTYPE_ORIGINAL = 'if str(work[col].dtype) in ("object", "category", "bool"):'
PANDAS_DTYPE_PORTABLE = "if not pd.api.types.is_numeric_dtype(work[col]):"


def parse_args():
    parser = argparse.ArgumentParser(description="Run the CADe gastroscopy analysis workflow.")
    parser.add_argument("--config", default="config.json", help="Path to the JSON input configuration.")
    parser.add_argument("--output-dir", default="results", help="Directory for text output and staged scripts.")
    parser.add_argument(
        "--cohort",
        choices=["cohort1", "cohort2", "both"],
        default="both",
        help="Run the core scripts for one or both cohorts.",
    )
    parser.add_argument("--core-only", action="store_true", help="Skip Cohort 1-only sensitivity and lesion scripts.")
    parser.add_argument("--continue-on-error", action="store_true", help="Continue after a script fails.")
    return parser.parse_args()


def load_config(path):
    config_path = Path(path).expanduser().resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    required = {"cohort1", "cohort2", "first_study", "screening_result", "cancer", "polyp"}
    missing = sorted(required.difference(config.get("data", {})))
    if missing:
        raise KeyError("Missing configuration entries: " + ", ".join(missing))
    paths = {}
    for key, value in config["data"].items():
        candidate = Path(value).expanduser()
        if not candidate.is_absolute():
            candidate = config_path.parent / candidate
        paths[key] = candidate.resolve()
    return paths


def stage_script(source, destination, inputs, cohort_path):
    text = source.read_text(encoding="utf-8-sig")
    replacements = dict(inputs)
    replacements["cohort"] = cohort_path
    for literal, key in PATH_PLACEHOLDERS.items():
        text = text.replace(literal, repr(str(replacements[key])))
    text = text.replace(PANDAS_DTYPE_ORIGINAL, PANDAS_DTYPE_PORTABLE)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")


def run_script(staged_script, output_path, continue_on_error):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    print(staged_script.parent.name, staged_script.name)
    completed = subprocess.run(
        [sys.executable, str(staged_script)],
        cwd=staged_script.parent,
        env=env,
        text=True,
        capture_output=True,
        encoding="utf-8",
    )
    output_path.write_text(completed.stdout + completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        message = f"{staged_script.name} failed with exit code {completed.returncode}; see {output_path}"
        if not continue_on_error:
            raise RuntimeError(message)
        print("WARNING:", message)


def prepare_and_run(script_dir, script_names, cohort_name, cohort_path, inputs, output_root, continue_on_error):
    stage_dir = output_root / "staged_scripts" / cohort_name
    result_dir = output_root / cohort_name
    for script_name in script_names:
        staged_script = stage_dir / script_name
        stage_script(script_dir / script_name, staged_script, inputs, cohort_path)
        run_script(staged_script, result_dir / f"{Path(script_name).stem}.txt", continue_on_error)
    manifest = (
        "Staged scripts generated for portable workflow execution.\n"
        "Portability operations:\n"
        "1. Source-data locations were resolved from the workflow configuration.\n"
        "2. The selected cohort file was bound to each analysis script.\n"
        "3. Categorical dtype detection was normalized across supported pandas versions.\n"
        "Statistical specifications were preserved during staging.\n"
    )
    (stage_dir / "STAGING_MANIFEST.txt").write_text(manifest, encoding="utf-8")


def main():
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    inputs = load_config(args.config)
    output_root = Path(args.output_dir).expanduser().resolve()
    requested = ["cohort1", "cohort2"] if args.cohort == "both" else [args.cohort]

    for cohort_name in requested:
        cohort_path = inputs[cohort_name]
        if not cohort_path.exists():
            raise FileNotFoundError(f"Configured {cohort_name} file does not exist: {cohort_path}")
        prepare_and_run(
            script_dir,
            CORE_SCRIPTS,
            cohort_name,
            cohort_path,
            inputs,
            output_root,
            args.continue_on_error,
        )

    if "cohort1" in requested and not args.core_only:
        for key in ["first_study", "screening_result", "cancer", "polyp"]:
            if not inputs[key].exists():
                raise FileNotFoundError(f"Configured input does not exist: {inputs[key]}")
        prepare_and_run(
            script_dir,
            COHORT1_ONLY_SCRIPTS,
            "cohort1",
            inputs["cohort1"],
            inputs,
            output_root,
            args.continue_on_error,
        )

    print(f"Workflow completed. Outputs: {output_root}")


if __name__ == "__main__":
    main()
