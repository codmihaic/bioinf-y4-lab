# labs/00_smoke/mlflow_smoke.py
import argparse, time, os, sys, mlflow
# ensure repo root on sys.path for 'mlops'
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)
from mlops.mlflow_utils import start_run, log_params, log_metrics

parser = argparse.ArgumentParser()
parser.add_argument("--experiment", default="BIOINF-Y4 Demo", help="MLflow experiment name")
args = parser.parse_args()

with start_run(args.experiment, run_name="mlflow_smoke") as _:
    log_params({"seed": 42, "model": "toy"})
    for step, val in enumerate([0.65, 0.72, 0.78, 0.81]):
        time.sleep(0.05)  # keep it snappy
        log_metrics({"auc": val}, step=step)
    # artifact
    with open("artifact.txt", "w") as f:
        f.write("hello mlflow")
    mlflow.log_artifact("artifact.txt")
print("OK: logged metrics/params/artifact to local mlruns/ (file store)")
