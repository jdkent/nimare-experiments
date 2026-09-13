"""Does rho_k track the per-study inflation it is supposed to remove?

The claim behind peak_bias="per-study" is falsifiable on real data. For each study we have
both the reported peaks (at a threshold we impose) and the *image* those peaks came from, so
the true inflation of study k is measurable:

    inflation_k = mean |g| study k reports  /  mean |g| the pooled images say at those voxels

If the model is right, inflation_k should be proportional to null_peak_mean_g(u_k, N_k), and
dividing by it should collapse the spread across studies. If it is wrong, the spread will not
move. Each study gets its own randomly assigned threshold, which is the heterogeneous
literature the correction exists for.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma.effectsize import null_peak_mean_g, infer_threshold_from_minimum
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.utils import get_masker, mm2vox

THRESHOLDS = np.array([2.3263, 3.0902, 3.2905, 4.2649])

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
ids = list(ss.ids)
rng = np.random.default_rng(0)
assigned = {str(sid): float(rng.choice(THRESHOLDS)) for sid in ids}

# Per-study g images, and the leave-one-out pooled reference (inverse-variance, no selection).
images = ss.images
g_maps, v_maps, sizes = {}, {}, {}
sample_sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))
for sid, gp, vp in zip(images["id"].astype(str), images["g"], images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g_maps[sid], v_maps[sid] = np.where(ok, g, 0.0), np.where(ok, v, np.inf)
    sizes[sid] = float(sample_sizes[sid])

stack_g = np.array([g_maps[s] for s in g_maps])
stack_w = np.array([1.0 / v_maps[s] for s in g_maps])
keys = list(g_maps)

# Coordinates at each study's own threshold.
rows = {}
for u in THRESHOLDS:
    cs = ImagesToCoordinates(
        merge_strategy="demolish", z_threshold=float(u), two_sided=True, remove_subpeaks=True
    ).transform(ss)
    df = cs.coordinates
    for sid, sub in df.groupby("id"):
        sid = str(sid)
        if assigned.get(sid) == float(u):
            rows[sid] = sub

print(f"{'study':>12s} {'u_k':>6s} {'N_k':>5s} {'m':>4s} {'|g| rep':>8s} {'|g| ref':>8s} "
      f"{'inflation':>10s} {'null_g':>7s} {'u_hat':>6s}")
records = []
for sid, sub in sorted(rows.items()):
    u = assigned[sid]
    idx = mm2vox(sub[["x", "y", "z"]].values, masker.mask_img.affine)
    flat = np.ravel_multi_index(idx.T, masker.mask_img.shape[:3], mode="clip")
    inmask = np.flatnonzero(masker.mask_img.get_fdata().ravel() > 0)
    lookup = {v: i for i, v in enumerate(inmask)}
    cols = np.array([lookup.get(f, -1) for f in flat])
    keep = cols >= 0
    cols = cols[keep]
    if not len(cols):
        continue

    # Leave-one-out reference at exactly the voxels this study reported.
    others = np.array([k != sid for k in keys])
    w = stack_w[others][:, cols]
    ref = np.abs((stack_g[others][:, cols] * w).sum(0) / w.sum(0))

    z = np.abs(sub["z_stat"].astype(float).to_numpy()[keep])
    from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
    g_rep, _ = peak_stat_to_hedges_g(z, np.full(len(z), sizes[sid]), stat_type="z")
    g_rep = np.abs(g_rep)

    null_g = null_peak_mean_g(u, sizes[sid])
    u_hat = infer_threshold_from_minimum(z.min(), len(z))
    records.append((sid, u, sizes[sid], len(z), g_rep.mean(), ref.mean(),
                    g_rep.mean() / ref.mean(), null_g, u_hat))
    print(f"{sid[:12]:>12s} {u:6.2f} {sizes[sid]:5.0f} {len(z):4d} {g_rep.mean():8.3f} "
          f"{ref.mean():8.3f} {g_rep.mean()/ref.mean():10.2f} {null_g:7.3f} {u_hat:6.2f}")

infl = np.array([r[6] for r in records])
null_g = np.array([r[7] for r in records])
u_true = np.array([r[1] for r in records])
u_hat = np.array([r[8] for r in records])

def cv(x):
    return float(np.std(x) / np.mean(x))

rho_oracle = np.median(null_g) / null_g                     # threshold known
null_g_hat = np.array([null_peak_mean_g(uh, r[2]) for uh, r in zip(u_hat, records)])
rho_inferred = np.median(null_g_hat) / null_g_hat           # threshold guessed from the minimum

print(f"\nstudies: {len(records)}")
print(f"corr(inflation, null_peak_mean_g)      = {np.corrcoef(infl, null_g)[0,1]:+.3f}")
print(f"threshold recovery: mean |u_hat - u|   = {np.mean(np.abs(u_hat - u_true)):.3f} "
      f"(bias {np.mean(u_hat - u_true):+.3f})")
print(f"\nspread of per-study inflation (coefficient of variation, lower is better)")
print(f"  uncorrected                          = {cv(infl):.3f}")
print(f"  rho from the true threshold          = {cv(infl * rho_oracle):.3f}")
print(f"  rho from the guessed threshold       = {cv(infl * rho_inferred):.3f}")
for name, rho in (("true", rho_oracle), ("guessed", rho_inferred)):
    lenient = infl[u_true < 3.0].mean() * rho[u_true < 3.0].mean()
    strict = infl[u_true > 3.5].mean() * rho[u_true > 3.5].mean()
    print(f"  strict/lenient gap ({name:>7s})          = {strict/lenient:.2f}")
raw_lenient, raw_strict = infl[u_true < 3.0].mean(), infl[u_true > 3.5].mean()
print(f"  strict/lenient gap (uncorrected)     = {raw_strict/raw_lenient:.2f}")
