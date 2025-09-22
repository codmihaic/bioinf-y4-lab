from mlops.mlflow_utils import start_run, log_metrics, log_params
with start_run(lab="00_smoke", pair="demo"):
    log_params({"model": "dummy", "version": "0.1"})
    log_metrics({"accuracy": 0.73})
print("mlflow demo ok")
