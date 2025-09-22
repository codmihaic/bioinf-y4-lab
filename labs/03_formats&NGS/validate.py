import json, sys, pathlib
p = pathlib.Path("results")
qc = p/"qc.json"
if not qc.exists():
    print("Missing results/qc.json"); sys.exit(1)
data = json.loads(qc.read_text())
mr = float(data.get("mapping_rate", -1))
if not (0.5 <= mr <= 1.0):
    print(f"Bad mapping_rate: {mr}"); sys.exit(1)
print("NGS validation: ok")
