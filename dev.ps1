param(
  [ValidateSet("jupyter","smoke","lint")] [string]$cmd = "jupyter"
)

if ($cmd -eq "jupyter") {
  docker run -it --rm -p 8890:8888 -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bioinf-y4-lab:base `
    bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"
}
elseif ($cmd -eq "smoke") {
  docker run --rm -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bioinf-y4-lab:base python labs/00_smoke/smoke.py
}
elseif ($cmd -eq "lint") {
  docker run --rm -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bioinf-y4-lab:base flake8
}
