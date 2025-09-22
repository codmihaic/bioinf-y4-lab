# labs/00_smoke/mlflow_smoke.py
"""CI smoke test for MLflow logging (file-based store)."""

import argparse
import os
import sys
import time

import mlflow


def main() -> None:
    # Ensure repo root on sys.path so 'mlops' can be imported when run from labs/*
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)

    # Import after sys.path fix
    from mlops.mlflow_utils import start_run, log_params, log_metrics

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--experiment",
        default="BIOINF-Y4 Demo",
        help="MLflow experiment name",
    )
    args = parser.parse_args()

    # Local file store; CI typically sets MLFLOW_TRACKING_DIR=/tmp/mlruns
    tracking_dir = os.environ.get("MLFLOW_TRACKING_DIR", "mlruns")
    os.environ["MLFLOW_TRACKING_URI"] = f"file:{tracking_dir}"

    with start_run(args.experiment, run_name="mlflow_smoke"):
        log_params({"seed": 42, "model": "toy"})
        for step, val in enumerate([0.65, 0.72, 0.78, 0.81]):
            time.sleep(0.05)  # keep CI snappy
            log_metrics({"auc": val}, step=step)

        # Log a tiny artifact
        artifact_path = os.path.join(os.getcwd(), "artifact.txt")
        with open(artifact_path, "w", encoding="utf-8") as handle:
            handle.write("hello mlflow")
        mlflow.log_artifact(artifact_path)

    print(
        "OK: logged metrics/params/artifact to local MLflow file store at "
        f"{tracking_dir}"
    )


if __name__ == "__main__":
    main()
