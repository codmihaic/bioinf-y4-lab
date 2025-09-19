setup: ; pip install -r requirements.txt
test: ; echo "smoke"
jupyter:
	docker run -it --rm -p 8890:8888 -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bdhb:base \
	bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"

smoke:
	docker run --rm -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bdhb:base python labs/00_smoke/smoke.py

lint:
	docker run --rm -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bdhb:base flake8
