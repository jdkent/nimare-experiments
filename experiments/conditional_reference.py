"""Validate CBES's g against a reference for mu, not for pi*mu.

Every real-data comparison so far scored CBES against an inverse-variance pooled image map.
That map averages every study, including those with nothing at the voxel, so it estimates the
marginal pi(v)*mu(v). CBES's ``g`` is mu(v), the effect among the studies that have an effect
there. Comparing them measures 1/pi and says nothing about whether the magnitude is right.

The reference built here conditions on the same event CBES conditions on: at each voxel,
average g over only the studies whose own map is active there. ``prevalence`` gets the same
treatment -- the fraction of studies active at the voxel is a direct empirical check on a map
that has so far only been validated against a simulator.

One caveat drives the design. Selecting studies on their own observed value at v induces the
same winner's curse in the reference, biasing it up, and the size of that bias depends on how
strict "active" is. So the activity threshold is swept rather than fixed, and CBES's estimate
is located against the resulting range instead of against a single number: a lenient threshold
gives a nearly unselected reference that is diluted by studies with little effect, a strict one
gives a heavily selected reference biased high, and the truth is bracketed between them.
"""
import json, os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
MAX_PEAKS = 10
ACTIVITY = (1.0, 2.0, 3.2905)

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))
coords = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                             remove_subpeaks=True).transform(ss).coordinates

# Per-study g and its z, so "active here" can be judged on the study's own evidence.
g_rows, z_rows = [], []
for row in ss.images.itertuples():
    if row.g is None or row.g_var is None:
        continue
    g = masker.transform(str(row.g)).ravel()
    v = masker.transform(str(row.g_var)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g = np.where(ok, g, 0.0)
    se = np.sqrt(np.where(ok, v, np.inf))
    g_rows.append(g)
    z_rows.append(np.where(ok, g / np.maximum(se, 1e-9), 0.0))
G, Z = np.array(g_rows), np.array(z_rows)
W = np.where(np.isfinite(G), 1.0, 0.0)
marginal = G.mean(axis=0)
print(f"{len(G)} studies, {G.shape[1]} voxels\n")

ids = sorted(set(str(i) for i in coords["id"]))
studies = []
for sid in ids:
    sub = coords[coords["id"].astype(str) == sid].copy()
    sub["_a"] = np.abs(sub["z_stat"].astype(float))
    sub = sub.sort_values("_a", ascending=False).head(MAX_PEAKS)
    meta = {"sample_sizes": [int(sizes[sid])], "reporting_threshold": U}
    studies.append({"id": sid, "name": sid, "metadata": meta, "analyses": [
        {"id": sid, "name": "1", "metadata": meta,
         "points": [{"space": "MNI", "coordinates": [float(r.x), float(r.y), float(r.z)],
                     "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                    for r in sub.itertuples()]}]})
studyset = Studyset({"id": "pain", "name": "pain", "studies": studies},
                    target=None, mask=masker.mask_img)
est = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
           threshold="reporting_threshold")
result = est.fit(studyset)
g_cbes = np.abs(result.get_map("g", return_type="array").ravel())
pi_cbes = result.get_map("prevalence", return_type="array").ravel()
covered = result.get_map("n_studies", return_type="array").ravel() > 0
print(f"CBES covers {covered.sum()} voxels; mean |g| there {g_cbes[covered].mean():.3f}, "
      f"mean prevalence {pi_cbes[covered].mean():.3f}")
print(f"marginal image mean |g| over covered: {np.abs(marginal[covered]).mean():.3f}\n")

print(f"{'active if |z|>':>14} {'voxels n>=2':>12} {'mu_hat':>8} {'CBES g':>8} {'ratio':>7} "
      f"{'r(g,mu)':>8} {'pi_hat':>7} {'CBES pi':>8} {'r(pi)':>7}")
for thresh in ACTIVITY:
    active = np.abs(Z) >= thresh
    n_active = active.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu_hat = np.where(n_active > 0, (np.abs(G) * active).sum(axis=0) / np.maximum(n_active, 1), np.nan)
    pi_hat = n_active / G.shape[0]
    use = covered & (n_active >= 2) & np.isfinite(mu_hat) & (g_cbes > 0)
    if use.sum() < 100:
        print(f"{thresh:14.2f} {use.sum():12d}  too few voxels")
        continue
    r_g = stats.pearsonr(g_cbes[use], mu_hat[use])[0]
    r_pi = stats.pearsonr(pi_cbes[use], pi_hat[use])[0]
    print(f"{thresh:14.2f} {use.sum():12d} {mu_hat[use].mean():8.3f} {g_cbes[use].mean():8.3f} "
          f"{g_cbes[use].mean()/mu_hat[use].mean():7.3f} {r_g:8.3f} "
          f"{pi_hat[use].mean():7.3f} {pi_cbes[use].mean():8.3f} {r_pi:7.3f}")
