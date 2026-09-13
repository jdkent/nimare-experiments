"""Run CoordinateEffectSize on the shared simulated studysets."""
import glob, json, sys, warnings
import numpy as np
warnings.simplefilter("ignore")
from nimare.studyset import Studyset
from nimare.meta.cbma.effectsize import CoordinateEffectSize

kind = sys.argv[1]
out = {}
for path in sorted(glob.glob(f"/tmp/claude-0/cmp/{kind}_*.json")):
    ss = Studyset(json.load(open(path)))
    est = CoordinateEffectSize(radius=10.0)
    res = est.fit(ss)
    table = res.tables["clusters"]
    if len(table) == 0:
        out[path] = dict(n_clusters=0, any_p05=False, any_fdr05=False, max_g=None)
        continue
    p = table["p"].to_numpy(dtype=float)
    pf = table["p_corr_fdr"].to_numpy(dtype=float)
    g = table["g"].to_numpy(dtype=float)
    out[path] = dict(
        n_clusters=int(len(table)),
        any_p05=bool(np.nanmin(p) < 0.05),
        any_fdr05=bool(np.nanmin(pf) < 0.05),
        max_g=float(np.nanmax(np.abs(g))) if np.isfinite(g).any() else None,
        min_p=float(np.nanmin(p)),
    )
json.dump(out, open(f"/tmp/claude-0/cmp/theirs_{kind}.json", "w"), indent=1)
rows = list(out.values())
print(f"{kind}: {len(rows)} sims, mean clusters {np.mean([r['n_clusters'] for r in rows]):.1f}")
print(f"  any cluster p<.05      : {np.mean([r['any_p05'] for r in rows]):.2f}")
print(f"  any cluster FDR q<.05  : {np.mean([r['any_fdr05'] for r in rows]):.2f}")
gs = [r['max_g'] for r in rows if r['max_g'] is not None]
print(f"  max |g| across clusters: {np.mean(gs):.3f}")
