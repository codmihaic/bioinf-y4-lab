#!/usr/bin/env python3
"""
Generate a third-party license report for the CURRENT Python environment.

Writes a Markdown file (default: LICENSES-THIRD-PARTY.md) listing:
- package name
- version
- license (best-effort)
- homepage / project URL (if available)

Usage:
  python scripts/generate_licenses.py
  python scripts/generate_licenses.py --output docs/THIRD_PARTY_LICENSES.md
"""
from __future__ import annotations
import argparse, re
from importlib import metadata
from typing import Optional, Tuple, List

LICENSE_CLASSIFIER_RE = re.compile(r"License ::(?:.*)::\s*(.+)$")

def best_effort_license(meta: metadata.PackageMetadata) -> Optional[str]:
    lic = meta.get("License")
    if lic:
        txt = lic.strip()
        if txt and txt.lower() not in {"unknown"}:
            return txt
    classifiers = meta.get_all("Classifier") or []
    for c in classifiers:
        m = LICENSE_CLASSIFIER_RE.search(c)
        if m:
            candidate = m.group(1).strip()
            if candidate and candidate.lower() not in {"other/proprietary license", "unknown"}:
                return candidate
    return None

def best_effort_homepage(meta: metadata.PackageMetadata) -> Optional[str]:
    hp = meta.get("Home-page")
    if hp:
        return hp.strip()
    project_urls = meta.get_all("Project-URL") or []
    preferred_labels = ("Homepage", "homepage", "Repository", "Source", "source")
    first_url = None
    for pu in project_urls:
        parts = re.split(r"[,|:]\s*", pu, maxsplit=1)
        if len(parts) == 2:
            label, url = parts[0].strip(), parts[1].strip()
            if not first_url:
                first_url = url
            if label in preferred_labels:
                return url
        elif len(parts) == 1 and parts[0].startswith("http"):
            if not first_url:
                first_url = parts[0].strip()
    return first_url

def safe_get_meta(dist: metadata.Distribution) -> metadata.PackageMetadata:
    try:
        return dist.metadata
    except Exception:
        class _Stub(dict):
            def get(self, k, default=None): return default
            def get_all(self, k): return []
        return _Stub()  # type: ignore

def collect_packages() -> List[Tuple[str, str, str, str]]:
    rows: List[Tuple[str, str, str, str]] = []
    for dist in sorted(metadata.distributions(), key=lambda d: d.metadata.get("Name","").lower()):
        meta = safe_get_meta(dist)
        name = meta.get("Name") or "unknown"
        version = meta.get("Version") or "unknown"
        lic = best_effort_license(meta) or "Unknown"
        homepage = best_effort_homepage(meta) or ""
        rows.append((name, version, lic, homepage))
    return rows

def render_markdown(rows: List[Tuple[str, str, str, str]]) -> str:
    out = []
    out.append("# 📜 Third-Party Licenses (Python)\n")
    out.append("> Generated from the active Python environment. Each package retains its original license.\n")
    out.append("\n| Package | Version | License | Homepage |")
    out.append("|---|---:|---|---|")
    for name, version, lic, homepage in rows:
        hp = f"[link]({homepage})" if homepage else ""
        lic_clean = lic.replace("|", "\\|")
        out.append(f"| {name} | {version} | {lic_clean} | {hp} |")
    out.append("\n\n_Last updated automatically; rerun the generator after dependency changes._\n")
    return "\n".join(out)

def main():
    ap = argparse.ArgumentParser(description="Generate third-party license inventory (Markdown).")
    ap.add_argument("--output", "-o", default="LICENSES-THIRD-PARTY.md",
                    help="Output Markdown path (default: LICENSES-THIRD-PARTY.md)")
    args = ap.parse_args()
    rows = collect_packages()
    md = render_markdown(rows)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write(md)
    print(f"Wrote {args.output} with {len(rows)} packages.")

if __name__ == "__main__":
    main()
