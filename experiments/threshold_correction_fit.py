"""Is the residual magnitude bias one function of the reporting threshold, across collections?

On the pain collection the ratio of CBES's g to a matched conditional reference moves smoothly
with the inferred cut -- 0.825 at 3.28, 0.955 at 3.63, 1.143 at 4.27. If that is one function
of an observable, then dividing by it makes g absolute without any unidentified constant, since
the estimator already recovers the cut to 0.01 z. But those three points came from capping the
same collection's peaks to a fixed count, which is not what a correction does: it fixes the
count and lets the cut float, and that cut is then each study's Nth peak -- high for a study
with strong signal, low for a weak one. The x-axis was therefore a function of the very effect
being estimated, not a reporting convention.

The separation comes from collections that differ in threshold for their own reasons, each with
its own smoothness, sample sizes and peak count. Pain plus the four NeuroVault paradigm groups
give five such points, and each is additionally run at several thresholds, so the design has
both within-collection and between-collection variation in the cut.

Reads: one function fitting every collection means the correction is estimable from coordinates.
Each collection needing its own intercept means the residual is carried by something unobserved
-- smoothness being the candidate -- and images stay necessary.
"""
import json, os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
import subprocess
from scipy import stats
from scipy.ndimage import maximum_filter
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates, t_to_z

U = 3.2905
# Vary the *threshold* and keep everything that survives, which is what a correction does.
# Capping to the top N instead would fix the count and let the cut float -- and that cut is
# then the value of each study's Nth peak, high for a strong study and low for a weak one, so
# it is a function of the very effect being estimated rather than a reporting convention.
THRESHOLDS = (3.2905, 3.75, 4.25, 4.75)
MIN_VOXELS = 100
NV_MAPS = "/tmp/claude-0/paradigms/maps.json"
NV_CACHE = "/tmp/claude-0/paradigms/nii"
EXCLUDE = ("none", "other", "null", "rest eyes open", "rest eyes closed", "none / other")

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
shape = mask_img.shape
mask_bool = np.asarray(mask_img.get_fdata() > 0)
affine = mask_img.affine


def peaks_from(z_masked, threshold):
    """Every local maximum of |z| clearing ``threshold`` -- all of them, as a paper reports."""
    volume = np.zeros(shape, dtype=float)
    volume[mask_bool] = z_masked
    magnitude = np.abs(volume)
    is_peak = (magnitude == maximum_filter(magnitude, size=3)) & (magnitude >= threshold) \
        & mask_bool
    idx = np.argwhere(is_peak)
    if not len(idx):
        return []
    values = volume[tuple(idx.T)]
    order = np.argsort(-np.abs(values))
    return [(idx[i], float(values[i])) for i in order]


def evaluate(name, per_study, threshold):
    """Fit CBES on everything above ``threshold`` and score against the matched reference.

    The reference judges a study active where its own map clears the same threshold the peaks
    were taken at -- a real reporting bound, common to every study, rather than a per-study
    value that moves with how much signal the study happens to have.
    """
    studies, cuts, rows_g, rows_z = [], [], [], []
    dropped = 0
    for sid, (z, g, n) in per_study.items():
        found = peaks_from(z, threshold)
        if len(found) < 2:
            dropped += 1
            continue
        cuts.append(threshold)
        rows_g.append(g)
        rows_z.append(z)
        meta = {"sample_sizes": [int(n)]}
        points = [{"space": "MNI",
                   "coordinates": [float(c) for c in nib.affines.apply_affine(
                       affine, np.asarray(ijk, dtype=float))],
                   "values": [{"kind": "Z", "value": value}]} for ijk, value in found]
        studies.append({"id": sid, "name": sid, "metadata": meta,
                        "analyses": [{"id": sid, "name": "1", "metadata": meta,
                                      "points": points}]})
    if len(studies) < 5:
        return None
    G, Z = np.array(rows_g), np.array(rows_z)
    bound = np.asarray(cuts)[:, None]
    active = np.abs(Z) >= bound
    n_active = active.sum(axis=0)
    with np.errstate(invalid="ignore"):
        mu = np.where(n_active > 0,
                      (np.abs(G) * active).sum(axis=0) / np.maximum(n_active, 1), np.nan)
    estimator = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
                     threshold="study-min")
    try:
        result = estimator.fit(Studyset({"id": name, "name": name, "studies": studies},
                                        target=None, mask=mask_img))
    except Exception as exc:
        print(f"  {name} thr={threshold}: fit failed ({exc})", flush=True)
        return None
    g = np.abs(result.get_map("g", return_type="array").ravel())
    covered = result.get_map("n_studies", return_type="array").ravel() > 0
    use = covered & (n_active >= 2) & np.isfinite(mu) & (g > 0)
    if use.sum() < MIN_VOXELS:
        return None
    inferred = float(np.median(np.abs(estimator._cutoffs_z_.values)))
    kept = [len(peaks_from(z, threshold)) for z, _, _ in per_study.values()]
    n_peaks = float(np.mean([k for k in kept if k >= 2])) if any(k >= 2 for k in kept) else 0.0
    return {
        "collection": name, "threshold": threshold, "studies": len(studies),
        "dropped": dropped, "voxels": int(use.sum()),
        "inferred_cut": inferred, "actual_cut": threshold, "peaks": n_peaks,
        "mu": float(mu[use].mean()), "g": float(g[use].mean()),
        "ratio": float(g[use].mean() / mu[use].mean()),
        "r": float(stats.pearsonr(g[use], mu[use])[0]),
    }


# ---- pain -------------------------------------------------------------------------------
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain := __import__("load_pain").load_pain())
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))
pain = {}
for row in ss.images.itertuples():
    if row.g is None or row.g_var is None:
        continue
    g_img = resample_to_img(nib.load(str(row.g)), mask_img, interpolation="continuous",
                            force_resample=True, copy_header=True)
    v_img = resample_to_img(nib.load(str(row.g_var)), mask_img, interpolation="continuous",
                            force_resample=True, copy_header=True)
    g = masker.transform(g_img).ravel()
    v = masker.transform(v_img).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g = np.where(ok, g, 0.0)
    z = np.where(ok, g / np.sqrt(np.where(ok, v, np.inf)), 0.0)
    pain[str(row.id)] = (z, g, sizes.get(str(row.id), 30))
print(f"pain: {len(pain)} studies with images", flush=True)

# ---- NeuroVault paradigm groups ---------------------------------------------------------
groups = {}
by = {}
for m in json.load(open(NV_MAPS)):
    if (m["paradigm"] or "").strip().lower() in EXCLUDE:
        continue
    by.setdefault(m["paradigm"], {}).setdefault(m["collection"], m)
for paradigm, members in sorted(by.items(), key=lambda kv: -len(kv[1]))[:4]:
    if len(members) < 8:
        continue
    per_study = {}
    for cid, entry in members.items():
        path = os.path.join(NV_CACHE, f"{entry['image']}.nii.gz")
        if not os.path.exists(path):
            subprocess.run(["curl", "-sL", "--max-time", "120", "-o", path, entry["url"]],
                           capture_output=True)
        if not os.path.exists(path) or os.path.getsize(path) < 2000:
            continue
        try:
            img = nib.load(path)
            if img.ndim > 3:
                img = nib.Nifti1Image(np.asarray(img.dataobj)[..., 0], img.affine, img.header)
            data = masker.transform(resample_to_img(
                img, mask_img, interpolation="continuous", force_resample=True,
                copy_header=True)).ravel().astype(float)
        except Exception:
            continue
        n = float(entry["n"])
        z = np.nan_to_num(t_to_z(data, dof=n - 1) if entry["map_type"] == "t" else data)
        g = z / np.sqrt(n)
        per_study[str(cid)] = (z, g, n)
    if len(per_study) >= 8:
        groups[paradigm] = per_study
        print(f"{paradigm[:44]}: {len(per_study)} collections", flush=True)

rows = []
print(f"\n{'collection':>26} {'thr':>5} {'st':>3} {'drop':>4} {'peaks':>6} {'infer':>6} "
      f"{'mu':>6} {'g':>6} {'ratio':>6} {'r':>6} {'vox':>7}")
for name, per_study in [("pain", pain)] + list(groups.items()):
    for threshold in THRESHOLDS:
        got = evaluate(name, per_study, threshold)
        if not got:
            continue
        rows.append(got)
        print(f"{got['collection'][:26]:>26} {got['threshold']:5.2f} {got['studies']:>3} "
              f"{got['dropped']:>4} {got['peaks']:6.1f} {got['inferred_cut']:6.2f} "
              f"{got['mu']:6.3f} {got['g']:6.3f} {got['ratio']:6.3f} {got['r']:6.3f} "
              f"{got['voxels']:7d}", flush=True)

json.dump(rows, open("/tmp/claude-0/paradigms/threshold_fit.json", "w"))
if len(rows) >= 4:
    cut = np.array([r["inferred_cut"] for r in rows])
    ratio = np.array([r["ratio"] for r in rows])
    peaks = np.array([r["peaks"] for r in rows])
    slope, intercept, r_val, p_val, se = stats.linregress(cut, ratio)
    print(f"\none function of the inferred cut, pooled over collections:")
    print(f"  ratio = {intercept:.3f} + {slope:.3f} * cut   r={r_val:.3f} p={p_val:.4f} "
          f"slope se {se:.3f}")
    resid = ratio - (intercept + slope * cut)
    print(f"  residual sd {resid.std(ddof=1):.3f}")
    for name in sorted({r["collection"] for r in rows}):
        sel = np.array([r["collection"] == name for r in rows])
        if sel.sum() >= 2:
            print(f"    {name[:36]:>36}: mean residual {resid[sel].mean():+.3f} "
                  f"over {int(sel.sum())} thresholds")
    s2, i2, r2, p2, _ = stats.linregress(np.log(peaks), ratio)
    print(f"  against log peaks instead: r={r2:.3f} p={p2:.4f}")
