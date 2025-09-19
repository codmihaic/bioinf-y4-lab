# BIOINF-Y4 Onboarding — Environment & First Run

This guide shows two supported ways to run the labs **identically** everywhere:

- **GitHub Codespaces** (fastest; pulls our prebuilt image from GHCR)
- **Local Docker** (Windows/macOS/Linux)

> The canonical environment is the prebuilt image: `ghcr.io/bozdogalex/BIOINF-Y4:base`.

---

## Option A — GitHub Codespaces (recommended)

1. On the repository page, click **Code → Create codespace on main**.
2. Wait while it **pulls** the image `ghcr.io/bozdogalex/BIOINF-Y4:base`. (No build needed.)
3. In VS Code (web), open the **Terminal** and run the smoke test:
   ```bash
   python labs/00_smoke/smoke.py
   ```
   Expected output:
   ```
   ok
   ```
4. Open Jupyter (optional):
   - **View → Command Palette** → “Jupyter: Create New Notebook”
   - Or open any `.ipynb` and run a cell.

### Codespaces tips
- If a notebook asks for a **kernel**, pick `Python 3.11` (container).
- If you see a “rebuild container” prompt, accept; it will re-pull the same GHCR image.

---

## Option B — Local Docker

### Windows PowerShell
```powershell
docker pull ghcr.io/bozdogalex/BIOINF-Y4:base
docker run -it --rm -p 8890:8888 -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/BIOINF-Y4:base `
  bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"
```
Open: **http://localhost:8890/lab**  
New cell → `print("ok")`

### macOS / Linux bash
```bash
docker pull ghcr.io/bozdogalex/BIOINF-Y4:base
docker run -it --rm -p 8888:8888 -v "$PWD:/work" -w /work ghcr.io/bozdogalex/BIOINF-Y4:base \
  bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"
```
Open: **http://localhost:8888/lab**

---

## Verification checklist

- [ ] **Smoke**: `python labs/00_smoke/smoke.py` → prints `ok`
- [ ] **Jupyter opens** in browser; a cell runs `print("ok")`
- [ ] **Same behavior** in Codespaces and local Docker

---

## Troubleshooting

- **Port already in use** → change the left side of `-p`, e.g. `-p 8891:8888` and open that port.  
- **Windows bind mount fails** → restart Docker Desktop; ensure WSL integration + drive sharing are enabled.  
- **Jupyter asks for token** → we run with no token locally; Codespaces adds its own auth.  
- **“ERR_EMPTY_RESPONSE”** → be sure Jupyter is started with `--allow-root` and open the correct mapped port.  
- **CI differs from local** → CI runs *inside* the same GHCR image. If it passes locally in that image, CI will pass too.
