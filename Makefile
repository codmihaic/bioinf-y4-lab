setup: ; pip install -r requirements.txt
test: ; echo "smoke"
jupyter:
	docker run -it --rm -p 8890:8888 -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/ghcr.io/bozdogalex/bioinf-y4-lab:base:base \
	bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"

smoke:
	docker run --rm -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bioinf-y4-lab:base python labs/00_smoke/smoke.py

lint:
	docker run --rm -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bioinf-y4-lab:base flake8

mlflow-smoke:
\tpython labs/00_smoke/mlflow_smoke.py --experiment "BIOINF-Y4 Demo"

mlflow-ui:
\tmlflow ui --backend-store-uri file://$$(pwd)/mlruns --host 0.0.0.0 --port 5000