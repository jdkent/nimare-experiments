"""Conditional reference again, with the reporting threshold treated consistently.

Two flaws in the first version. It capped each study at its ten strongest peaks -- which is
itself a stringent threshold, the effective cut being the tenth largest peak -- while telling
CBES the threshold was 3.29, so the estimator was handed a misspecified bound. And it defined
"active" for the reference at |z| > 3.29 uncorrected, whereas a real paper reports peaks that
survived FWE or FDR, a far higher bar; the effect among studies that actually report is
conditioned on that, not on an uncorrected cut.

Both are fixed here. CBES infers each study's threshold from its own smallest reported peak
(threshold="study-min") rather than being told one, and the reference's activity threshold is
swept across the range real corrections imply, up to the whole-brain FWE region. The question
is no longer "is CBES above the reference" but "at what reporting threshold does the reference
agree with CBES, and is that the threshold the coordinates actually came from".
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates

EXTRACT_U = 3.2905          # the cut used to pull peaks out of the images
ACTIVITY = (2.0, 3.29, 4.0, 4.5, 5.0, 5.5, 6.0)
CAPS = (10, 40, None)

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))
coords = ImagesToCoordinates(merge_strategy="demolish", z_threshold=EXTRACT_U, two_sided=True,
                             remove_subpeaks=True).transform(ss).coordinates

g_rows, z_rows = [], []
for row in ss.images.itertuples():
    if row.g is None or row.g_var is None:
        continue
    g = masker.transform(str(row.g)).ravel()
    v = masker.transform(str(row.g_var)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g_rows.append(np.where(ok, g, 0.0))
    se = np.sqrt(np.where(ok, v, np.inf))
    z_rows.append(np.where(ok, g_rows[-1] / np.maximum(se, 1e-9), 0.0))
G, Z = np.array(g_rows), np.array(z_rows)
print(f"{len(G)} studies, {G.shape[1]} voxels, peaks extracted at |z| > {EXTRACT_U}\n")

ids = sorted(set(str(i) for i in coords["id"]))


def build(cap):
    studies, effective = [], []
    for sid in ids:
        sub = coords[coords["id"].astype(str) == sid].copy()
        sub["_a"] = np.abs(sub["z_stat"].astype(float))
        sub = sub.sort_values("_a", ascending=False)
        if cap:
            sub = sub.head(cap)
        effective.append(float(sub["_a"].min()))
        # No reporting_threshold in the metadata: the estimator infers it from the peaks.
        meta = {"sample_sizes": [int(sizes[sid])]}
        studies.append({"id": sid, "name": sid, "metadata": meta, "analyses": [
            {"id": sid, "name": "1", "metadata": meta,
             "points": [{"space": "MNI", "coordinates": [float(r.x), float(r.y), float(r.z)],
                         "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                        for r in sub.itertuples()]}]})
    return Studyset({"id": "pain", "name": "pain", "studies": studies},
                    target=None, mask=masker.mask_img), np.median(effective)


def matched_reference(cap):
    """mu_hat where each study is judged against *its own* effective threshold.

    A paper reporting five coordinates is not choosing five; five is what survived its
    correction, and it reports all of them. So the event "this study has an effect here" is
    "this study's own map clears the cut that its own reported set implies" -- a per-study
    bound, not a common one. The cut is taken as the study's smallest retained peak, which is
    exactly what threshold="study-min" gives the estimator.
    """
    cuts = []
    for sid in ids:
        sub = coords[coords["id"].astype(str) == sid].copy()
        sub["_a"] = np.abs(sub["z_stat"].astype(float))
        sub = sub.sort_values("_a", ascending=False)
        if cap:
            sub = sub.head(cap)
        cuts.append(float(sub["_a"].min()))
    cuts = np.asarray(cuts)[:, None]
    active = np.abs(Z) >= cuts
    n_active = active.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu_hat = np.where(n_active > 0,
                          (np.abs(G) * active).sum(axis=0) / np.maximum(n_active, 1), np.nan)
    return mu_hat, n_active / G.shape[0], n_active, cuts.ravel()


# The reference, as a function of how stringently "active" is defined.
print("reference mu_hat, by the threshold a study's own map must clear at the voxel:")
print(f"{'|z| >':>7} {'voxels n>=2':>12} {'mu_hat':>8} {'pi_hat':>8}")
ref = {}
for thresh in ACTIVITY:
    active = np.abs(Z) >= thresh
    n_active = active.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu_hat = np.where(n_active > 0,
                          (np.abs(G) * active).sum(axis=0) / np.maximum(n_active, 1), np.nan)
    ref[thresh] = (mu_hat, n_active / G.shape[0], n_active)
    keep = (n_active >= 2) & np.isfinite(mu_hat)
    print(f"{thresh:7.2f} {int(keep.sum()):12d} {mu_hat[keep].mean():8.3f} "
          f"{(n_active / G.shape[0])[keep].mean():8.3f}")

print(f"\n{'cap':>6} {'eff. cut':>9} {'CBES g':>8} {'CBES pi':>8} "
      f"{'matches mu_hat at |z|>':>23} {'best r(g,mu)':>13}")
for cap in CAPS:
    studyset, effective = build(cap)
    est = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
               threshold="study-min")
    result = est.fit(studyset)
    g = np.abs(result.get_map("g", return_type="array").ravel())
    pi = result.get_map("prevalence", return_type="array").ravel()
    covered = result.get_map("n_studies", return_type="array").ravel() > 0
    inferred = np.median(np.abs(est._cutoffs_z_.values))
    # Which activity threshold makes the reference agree with CBES?
    best, best_r = None, None
    for thresh in ACTIVITY:
        mu_hat, pi_hat, n_active = ref[thresh]
        use = covered & (n_active >= 2) & np.isfinite(mu_hat) & (g > 0)
        if use.sum() < 500:
            continue
        ratio = g[use].mean() / mu_hat[use].mean()
        r = stats.pearsonr(g[use], mu_hat[use])[0]
        if best is None or abs(np.log(ratio)) < abs(np.log(best[1])):
            best = (thresh, ratio)
        if best_r is None or r > best_r[1]:
            best_r = (thresh, r)
    use0 = covered & (g > 0)
    print(f"{str(cap):>6} {effective:9.2f} {g[use0].mean():8.3f} {pi[use0].mean():8.3f} "
          f"{f'{best[0]:.2f} (ratio {best[1]:.2f})':>23} "
          f"{f'{best_r[1]:.3f} at {best_r[0]:.2f}':>13}   [inferred cut {inferred:.2f}]",
          flush=True)
    # The matched comparison: each study judged against its own effective cut.
    mu_m, pi_m, n_m, cuts = matched_reference(cap)
    use = covered & (n_m >= 2) & np.isfinite(mu_m) & (g > 0)
    if use.sum() >= 500:
        print(f"       matched per-study cuts (median {np.median(cuts):.2f}): "
              f"mu_hat {mu_m[use].mean():.3f}  CBES g {g[use].mean():.3f}  "
              f"ratio {g[use].mean() / mu_m[use].mean():.3f}  "
              f"r {stats.pearsonr(g[use], mu_m[use])[0]:.3f}  |  "
              f"pi_hat {pi_m[use].mean():.3f}  CBES pi {pi[use].mean():.3f}  "
              f"r {stats.pearsonr(pi[use], pi_m[use])[0]:.3f}", flush=True)
    else:
        print(f"       matched per-study cuts: only {int(use.sum())} voxels, skipped",
              flush=True)
