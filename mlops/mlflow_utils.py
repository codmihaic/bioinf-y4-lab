from contextlib import contextmanager
from typing import Optional, Dict, Any
import mlflow
from pathlib import Path

MLRUNS_PATH = Path("./mlruns").resolve()

@contextmanager
def start_run(lab: str, pair: str, tags: Optional[Dict[str, Any]] = None):
    mlflow.set_tracking_uri(f"file:{MLRUNS_PATH.as_posix()}")
    with mlflow.start_run(run_name=f"{lab}::{pair}"):
        mlflow.set_tags({"lab": lab, "pair": pair, **(tags or {})})
        yield

def log_params(d: Dict[str, Any]):
    mlflow.log_params(d)

def log_metrics(d: Dict[str, float], step: Optional[int] = None):
    for k, v in d.items():
        mlflow.log_metric(k, float(v), step=step)

def log_artifact(path: str):
    mlflow.log_artifact(path)
