ML_IMAGE    ?= ghcr.io/bozdogalex/bioinf-y4-lab:base
ML_STORE    ?= /work/mlruns
ML_EXPER    ?= BIOINF-Y4 Demo
PORT        ?= 5000
JUPYTER_PORT  ?= 8890

.PHONY: jupyter
jupyter:
	# Start Jupyter Lab in the project container (no auth token)
	docker run --rm -it -p $(JUPYTER_PORT):8888 -v "$(CURDIR):/work" -w /work \
		$(ML_IMAGE) \
		bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"

.PHONY: mlflow-smoke
mlflow-smoke:
	docker run --rm -v "${PWD}:/work" -w /work \
		-e MLFLOW_TRACKING_DIR=$(ML_STORE) \
		$(ML_IMAGE) \
		python labs/00_smoke/mlflow_smoke.py --experiment "$(ML_EXPER)"

.PHONY: mlflow-ui
mlflow-ui:
	docker run --rm -p $(PORT):5000 -v "${PWD}:/work" -w /work \
		$(ML_IMAGE) \
		mlflow ui --backend-store-uri file:$(ML_STORE) --host 0.0.0.0 --port 5000

.PHONY: mlflow-clean
mlflow-clean:
	-rm -rf mlruns
