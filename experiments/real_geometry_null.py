"""The mismatch on real anatomy, with the signal taken out.

On the toy mask the relocation null had *fewer* studies per covered voxel than the observed
foci (1.95 against 3.81) and the test was anticonservative. On the real MNI mask with real pain
peaks the mean runs the other way -- 8.86 relocated against 7.06 observed -- because coverage
saturates there: when every voxel is reached by someone, spreading the foci uniformly maximises
the number of study-voxel pairs and clustered foci overlap within study instead.

But a mean over a saturated brain hides the thing that matters. Inference happens at the voxels
with the *most* studies on them, and the observed map's multiplicity is concentrated where the
literature converged while the null's is flat. So the comparison has to be made in the upper
tail of the multiplicity distribution, not at its mean.

Two measurements, both on the real geometry:

    multiplicity quantiles      observed against relocated, through the upper tail. If the
                                observed tail sits above the null's, the null is too tight
                                where it counts, whatever the means do.

    false positive rate         real focus positions, so the real clumping and the real mask,
                                but each peak's magnitude redrawn from the null distribution of
                                a peak that just cleared that study's threshold. No effect
                                anywhere, real anatomy everywhere.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats as st
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.cbma.effectsize import _p_from_histogram, _null_bin_edges, _NULL_MAX_Z

ITERS = int(sys.argv[1]) if len(sys.argv) > 1 else 30
SIMS = int(sys.argv[2]) if len(sys.argv) > 2 else 10
U = 3.2905

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
ss = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                         remove_subpeaks=True).transform(ss)

est = CBES(fwhm=10.0, null_method="none", threshold="study-min", peak_bias="per-study",
           use_images=False, seed=0)
est.fit(ss)
table = est._focus_table_
ijk = est._in_mask_ijk()
print(f"NIDM pain peaks on the real MNI mask: {len(table)} foci, "
      f"{int(est._mask_bool().sum())} voxels, {ITERS} relocations.\n", flush=True)


def stats(t):
    fit, z = est._statistic(t, est._sample_sizes_, est._thresholds_, est._image_studies_)
    return fit["n_studies"][fit["covered"]], np.abs(z)


QS = [0.50, 0.90, 0.99, 0.999, 1.0]
n_obs, _ = stats(table)
rng = np.random.default_rng(11)
n_null = []
for _ in range(ITERS):
    t = table.copy()
    t[["i", "j", "k"]] = ijk[rng.integers(0, len(ijk), size=len(t))]
    n_null.append(stats(t)[0])
n_null = np.concatenate(n_null)

print("  studies contributing per covered voxel", flush=True)
print(f"    {'q':>7s} {'observed':>9s} {'relocated':>10s} {'ratio':>7s}", flush=True)
for q in QS:
    o, r = float(np.quantile(n_obs, q)), float(np.quantile(n_null, q))
    print(f"    {q:7.3f} {o:9.2f} {r:10.2f} {o / r:7.2f}", flush=True)

# ---- false positive rate on the real geometry, magnitudes redrawn under the null
print(f"\n  false positive rate, {SIMS} redraws (nominal 0.050 / 0.0010)", flush=True)
print(f"    {'selection':>14s} {'unc p<.05':>10s} {'unc p<.001':>11s}", flush=True)
# the z-scale cutoffs, not est._thresholds_, which _apply_peak_bias puts on the g scale
cut = {str(sid): abs(float(v)) for sid, v in est._cutoffs_z_.items()}
for sel in ("zero-inflated", "none"):
    u5, u1 = [], []
    for sim in range(SIMS):
        rng = np.random.default_rng(500 + sim)
        e = CBES(fwhm=10.0, null_method="none", selection_model=sel, threshold="study-min",
                 peak_bias="per-study", use_images=False, seed=sim)
        e.fit(ss)
        t = e._focus_table_.copy()
        # a peak of pure noise that just cleared this study's threshold: |z| from the tail of a
        # standard normal above the cutoff, converted back onto the study's g scale
        thr = np.array([cut.get(str(s), U) for s in t["id"]], dtype=float)
        z_null = st.norm.isf(st.norm.sf(thr) * rng.random(len(t)))
        stat = t["stat"].values.astype(float)
        scale = np.where(stat != 0, t["g"].values / np.where(stat != 0, stat, 1.0),
                         np.sqrt(t["var_g"].values))
        t["g"] = z_null * scale * rng.choice([-1.0, 1.0], size=len(t))
        t["stat"] = z_null
        e._focus_table_ = t
        _, z = stats(t)
        hist = np.zeros(len(_null_bin_edges()) - 1, dtype=float)
        for _ in range(ITERS):
            tt = t.copy()
            tt[["i", "j", "k"]] = ijk[rng.integers(0, len(ijk), size=len(tt))]
            _, zz = stats(tt)
            counts, _ = np.histogram(np.clip(zz, 0, _NULL_MAX_Z), bins=_null_bin_edges())
            hist += counts
        p = _p_from_histogram(z, hist)
        u5.append(float(np.mean(p < 0.05)))
        u1.append(float(np.mean(p < 0.001)))
    print(f"    {sel:>14s} {np.mean(u5):10.4f} {np.mean(u1):11.5f}", flush=True)
