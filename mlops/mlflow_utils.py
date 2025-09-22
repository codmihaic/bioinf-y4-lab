# labs/00_smoke/mlflow_smoke.py
import argparse
import time

import mlflow

from mlops.mlflow_utils import start_run, log_params, log_metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", default="BIOINF-Y4 Demo")
    args = parser.parse_args()

    with start_run(args.experiment, run_name="mlflow_smoke"):
        log_params({"seed": 42, "model": "toy"})
        for step, val in enumerate([0.65, 0.72, 0.78, 0.81]):
            time.sleep(0.05)
            log_metrics({"auc": val}, step=step)
        with open("artifact.txt", "w") as f:
            f.write("hello mlflow")
        mlflow.log_artifact("artifact.txt")

    print("OK: metrics/params/artifact logged to ./mlruns (file store)")


if __name__ == "__main__":
    main()
