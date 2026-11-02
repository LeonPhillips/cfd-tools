"""python grid_lines_json.py <grid_lines.npz> <out.json>: the line-profile data as JSON (4 significant figures)."""
import json
import sys
import numpy as np
z = np.load(sys.argv[1], allow_pickle=True)
r4 = lambda a: [None if not np.isfinite(v) else float(f"{v:.4g}") for v in a]
out = {"labels": [str(s) for s in z["labels"]], "lines": []}
for i, (x0, r, x1, n) in enumerate(z["lines_mm"]):
    x = np.linspace(x0, x1, z[f"{out['labels'][0]}_Mach"].shape[1])
    out["lines"].append({"r": float(r), "x": r4(x), **{f"{l}_{v}": r4(z[f"{l}_{v}"][i]) for l in out["labels"]
                                                         for v in ("Mach", "Density")}})
json.dump(out, open(sys.argv[2], "w"), separators=(",", ":"))
print("saved", sys.argv[2])
