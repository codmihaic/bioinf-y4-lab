FROM python:3.11-slim

ENV PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential git curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /work
COPY requirements.txt /work/requirements.txt

RUN pip install --upgrade pip \
 && pip install -r requirements.txt

# Default shell; Jupyter launched by you (README command) or devcontainer UI
CMD ["/bin/bash"]
