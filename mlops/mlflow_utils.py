"""Tiny MLflow helpers used across labs and CI."""

from __future__ import annotations

import os
from contextlib import contextmanager
from typing import Any, Dict, Optional

import mlflow

__all__ = ["start_run", "log_params", "log_metrics"]


def _clean_params(params: Dict[str, Any]) -> Dict[str, str]:
    """Coerce params to strings; keep keys as strings."""
    out: Dict[str, str] = {}
    for k, v in params.items():
        key = str(k)
        if isinstance(v, (str, int, float, bool)) or v is None:
            out[key] = str(v)
        else:
            out[key] = str(v)
    return out


def _clean_metrics(metrics: Dict[str, Any]) -> Dict[str, float]:
    """Coerce metrics to floats; drop values that cannot be parsed."""
    out: Dict[str, float] = {}
    for k, v in metrics.items():
        key = str(k)
        if isinstance(v, (int, float)):
            out[key] = float(v)
            continue
        try:
            out[key] = float(v)  # type: ignore[arg-type]
        except Exception:
            # skip non-numeric metrics
            continue
    return out


@contextmanager
def start_run(
    experiment: str,
    run_name: Optional[str] = None,
    tags: Optional[Dict[str, str]] = None,
    nested: bool = False,
):
    """
    Context manager to open an MLflow run under a named experiment.

    If MLFLOW_TRACKING_DIR is set (e.g., /tmp/mlruns in CI) and MLFLOW_TRACKING_URI
    is not, we set a local file-based store: file:${MLFLOW_TRACKING_DIR}.
    """
    tracking_dir = os.environ.get("MLFLOW_TRACKING_DIR")
    if tracking_dir and "MLFLOW_TRACKING_URI" not in os.environ:
        os.environ["MLFLOW_TRACKING_URI"] = f"file:{tracking_dir}"

    mlflow.set_experiment(experiment)
    with mlflow.start_run(run_name=run_name, tags=tags, nested=nested):
        yield


def log_params(params: Dict[str, Any]) -> None:
    """Log a flat dict of parameters (values coerced to strings)."""
    mlflow.log_params(_clean_params(params))


def log_metrics(metrics: Dict[str, Any], step: Optional[int] = None) -> None:
    """Log numeric metrics (values coerced to floats)."""
    cleaned = _clean_metrics(metrics)
    if cleaned:
        mlflow.log_metrics(cleaned, step=step)
