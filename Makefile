# --- Common vars ---
ML_IMAGE       ?= ghcr.io/bozdogalex/bioinf-y4-lab:base
ML_EXPER       ?= BIOINF-Y4 Demo
ML_STORE       ?= mlruns           # relative path in repo (persisted)
PORT           ?= 5000
JUPYTER_PORT   ?= 8890

# Helper: 'HAS_DOCKER' is non-empty if docker is available
HAS_DOCKER := $(shell command -v docker 2>/dev/null)

.PHONY: jupyter
jupyter:
ifeq ($(HAS_DOCKER),)
	# No docker (Codespaces/devcontainer) → run Jupyter directly
	jupyter lab --ip=0.0.0.0 --port=$(JUPYTER_PORT) --no-browser --IdentityProvider.token='' --allow-root
else
	# Local host with docker → run in canonical image
	docker run --rm -it -p $(JUPYTER_PORT):8888 -v "$(CURDIR):/work" -w /work \
		$(ML_IMAGE) \
		bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"
endif

.PHONY: mlflow-smoke
mlflow-smoke:
ifeq ($(HAS_DOCKER),)
	# No docker → run smoke directly (file store in ./mlruns)
	MLFLOW_TRACKING_DIR=$(ML_STORE) \
	python labs/00_smoke/mlflow_smoke.py --experiment "$(ML_EXPER)"
else
	# With docker → run smoke in image (store /work/mlruns -> ./mlruns)
	docker run --rm -v "$(CURDIR):/work" -w /work \
		-e MLFLOW_TRACKING_DIR=/work/$(ML_STORE) \
		$(ML_IMAGE) \
		python labs/00_smoke/mlflow_smoke.py --experiment "$(ML_EXPER)"
endif

.PHONY: mlflow-ui
mlflow-ui:
ifeq ($(HAS_DOCKER),)
	# No docker → serve UI directly from ./mlruns
	mlflow ui --backend-store-uri file:$(ML_STORE) --host 0.0.0.0 --port $(PORT)
else
	# With docker → serve UI from /work/mlruns
	docker run --rm -p $(PORT):5000 -v "$(CURDIR):/work" -w /work \
		$(ML_IMAGE) \
		mlflow ui --backend-store-uri file:/work/$(ML_STORE) --host 0.0.0.0 --port 5000
endif
