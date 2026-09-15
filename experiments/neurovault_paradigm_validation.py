"""Validate CBES on independent NeuroVault collections that share a cognitive paradigm.

Every real-data number so far comes from the NIDM pain collection: 21 studies, one paradigm,
one curation effort. This is the independent replication. Studies are separate NeuroVault
collections annotated with the same ``cognitive_paradigm_cogatlas``, each contributing one
group-level t or z map with a declared sample size.

The protocol matches the pain check so the two are comparable:

  * truth is the inverse-variance pooled Hedges' g across the collections' own maps;
  * CBES sees only coordinates -- local maxima of |z| clearing a fixed threshold, capped
    strongest-first at a realistic number per study, which is what a paper tabulates;
  * scored on the magnitude ratio where the truth is largest, and on the spatial correlation
    across covered voxels.

Peaks are extracted here rather than through ImagesToCoordinates so that the threshold, the
cap and the t-to-z conversion are explicit and identical for every map.
"""
import json, os, subprocess, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from scipy.ndimage import maximum_filter
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import d_to_g, t_to_d, t_to_z

MAPS = "/tmp/claude-0/paradigms/maps.json"
CACHE = "/tmp/claude-0/paradigms/nii"
os.makedirs(CACHE, exist_ok=True)
U = 3.2905
MIN_COLLECTIONS = int(sys.argv[1]) if len(sys.argv) > 1 else 8
MAX_PEAKS = 10

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
affine = mask_img.affine
shape = mask_img.shape
mask_bool = np.asarray(mask_img.get_fdata() > 0)


def fetch(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 2000:
        return True
    subprocess.run(["curl", "-sL", "--max-time", "120", "-o", path, url], capture_output=True)
    return os.path.exists(path) and os.path.getsize(path) > 2000


def load_as_z(entry):
    """Return the map on the z scale, masked, or None if it cannot be used."""
    path = os.path.join(CACHE, f"{entry['image']}.nii.gz")
    if not entry.get("url") or not fetch(entry["url"], path):
        return None
    try:
        img = nib.load(path)
        if img.ndim > 3:
            img = nib.Nifti1Image(np.asarray(img.dataobj)[..., 0], img.affine, img.header)
        img = resample_to_img(img, mask_img, interpolation="continuous", force_resample=True,
                              copy_header=True)
        data = masker.transform(img).ravel().astype(float)
    except Exception:
        return None
    if not np.isfinite(data).any() or np.allclose(np.nan_to_num(data), 0):
        return None
    n = float(entry["n"])
    if entry["map_type"] == "t":
        data = t_to_z(data, dof=n - 1)
    return np.nan_to_num(data)


def peaks_of(z_masked, n_subjects, max_peaks):
    """Local maxima of |z| above the threshold, strongest first, as (ijk, z)."""
    volume = np.zeros(shape, dtype=float)
    volume[mask_bool] = z_masked
    magnitude = np.abs(volume)
    is_peak = (magnitude == maximum_filter(magnitude, size=3)) & (magnitude >= U) & mask_bool
    idx = np.argwhere(is_peak)
    if not len(idx):
        return []
    values = volume[tuple(idx.T)]
    order = np.argsort(-np.abs(values))[:max_peaks]
    return [(idx[i], float(values[i])) for i in order]


def to_g(z_masked, n_subjects):
    """Whole map on the Hedges' g scale, for the truth."""
    t = np.sign(z_masked) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z_masked)), n_subjects - 1))
    t = np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0)
    d = t_to_d(t, n_subjects)
    return d_to_g(d, n_subjects)


# Resting state has no task contrast, so a "effect size" there is not the quantity this
# estimator is about; the unlabelled bucket is not a paradigm at all.
EXCLUDE = ("none", "other", "null", "rest eyes open", "rest eyes closed", "none / other",
           "resting state")
entries = json.load(open(MAPS))
by = {}
for m in entries:
    if (m["paradigm"] or "").strip().lower() in EXCLUDE:
        continue
    # One map per collection per paradigm, so no single collection dominates a group.
    by.setdefault(m["paradigm"], {}).setdefault(m["collection"], m)
groups = sorted(((p, list(c.values())) for p, c in by.items()),
                key=lambda kv: -len(kv[1]))
groups = [(p, ms) for p, ms in groups if len(ms) >= MIN_COLLECTIONS]
print(f"{len(groups)} paradigms with >= {MIN_COLLECTIONS} collections\n", flush=True)

for paradigm, members in groups[:4]:
    loaded = []
    for entry in members:
        z = load_as_z(entry)
        if z is not None:
            loaded.append((entry, z))
    if len(loaded) < MIN_COLLECTIONS:
        print(f"{paradigm[:44]}: only {len(loaded)} maps loaded, skipping", flush=True)
        continue

    # Truth: inverse-variance pooled g over the same maps.
    g_stack, w_stack = [], []
    for entry, z in loaded:
        n = float(entry["n"])
        g = to_g(z, n)
        var = 1.0 / n + g**2 / (2.0 * n)
        g_stack.append(g)
        w_stack.append(1.0 / np.maximum(var, 1e-9))
    G, W = np.array(g_stack), np.array(w_stack)
    truth = np.sum(G * W, axis=0) / np.maximum(np.sum(W, axis=0), 1e-12)

    studies = []
    total_peaks = 0
    for entry, z in loaded:
        found = peaks_of(z, float(entry["n"]), MAX_PEAKS)
        if not found:
            continue
        total_peaks += len(found)
        meta = {"sample_sizes": [int(entry["n"])], "reporting_threshold": U}
        points = []
        for ijk, value in found:
            xyz = nib.affines.apply_affine(affine, np.asarray(ijk, dtype=float))
            points.append({"space": "MNI", "coordinates": [float(c) for c in xyz],
                           "values": [{"kind": "Z", "value": value}]})
        sid = str(entry["collection"])
        studies.append({"id": sid, "name": sid, "metadata": meta,
                        "analyses": [{"id": sid, "name": "1", "metadata": meta,
                                      "points": points}]})
    if len(studies) < MIN_COLLECTIONS:
        print(f"{paradigm[:44]}: only {len(studies)} report peaks, skipping", flush=True)
        continue
    studyset = Studyset({"id": "nv", "name": "nv", "studies": studies}, target=None,
                        mask=mask_img)

    hot = np.abs(truth) >= np.percentile(np.abs(truth), 75)
    print(f"--- {paradigm[:60]} ---", flush=True)
    print(f"    {len(studies)} collections, {total_peaks} peaks "
          f"({total_peaks / len(studies):.1f}/study), N "
          f"{min(int(e['n']) for e, _ in loaded)}-{max(int(e['n']) for e, _ in loaded)}, "
          f"truth top-quartile mean {np.abs(truth)[hot].mean():.3f}", flush=True)
    for peak_bias in (None, "per-study"):
        estimator = CBES(fwhm=10.0, mask=masker, peak_bias=peak_bias, null_method="none",
                         threshold="reporting_threshold")
        try:
            result = estimator.fit(studyset)
        except Exception as exc:
            print(f"    peak_bias={peak_bias}: fit failed ({exc})", flush=True)
            continue
        g = result.get_map("g", return_type="array").ravel()
        covered = result.get_map("n_studies", return_type="array").ravel() > 0
        both = covered & hot & np.isfinite(g)
        r = stats.pearsonr(np.abs(g[covered]), np.abs(truth[covered]))[0]
        print(f"    peak_bias={str(peak_bias):>9}: |g| hot {np.abs(g[both]).mean():6.3f}  "
              f"ratio {np.abs(g[both]).mean() / np.abs(truth[both]).mean():6.3f}  "
              f"r {r:6.3f}  covered {covered.sum()}", flush=True)
