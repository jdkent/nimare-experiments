"""The convergence-versus-magnitude comparison on a second, independent collection.

On the 21-study NIDM pain collection, CBES's `g` was statistically indistinguishable from MKDA
density at localising a held-out image reference (AUC difference +0.001, paired p 0.980) while
`g_marginal` beat it by +0.111 at p 0.001. One collection cannot separate a property of the
method from a property of that literature: pain is a strong, spatially consistent, high-prevalence
effect, the friendliest case there is. And a single-source result is exactly what the retracted
height-flattening claim rested on before a real-data check reversed it.

So the same procedure on the NeuroVault "animal" set: 11 studies, movie-watching, as-Animal
contrast, with unthresholded z maps to extract from and g/g_var maps for the reference. Split in
half, one half's tables against the other half's inverse-variance pooling, every estimator reading
exactly the same tables, scored on localisation only.

Eleven studies means five per half, which is thin -- so the question this answers is whether the
*ordering* of the arms survives a different literature, not whether the margins are the same size.
"""
import glob, logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES, ALE, MKDADensity, KDA
from nimare.studyset import Studyset
from reporting import report_peaks

SCHEME, FOCUS = "cluster", "max"
N_SPLITS = int(os.environ.get("NSPLITS", 12))
HOME = os.path.expanduser("~/.nimare")
#: The animal collection's per-study N, which its metadata does not carry alongside the maps.
#: Movie-watching samples of this vintage are small; 20 is used for every study, which makes the
#: z-to-g conversion uniform across the roster rather than accidentally informative.
ASSUMED_N = 20

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


def load(pattern):
    out = {}
    for path in sorted(glob.glob(f"{HOME}/{pattern}")):
        study = os.path.basename(path).split("-")[1]
        img = resample_to_img(nib.load(path), mask_img, interpolation="continuous",
                              force_resample=True, copy_header=True)
        out[study] = np.nan_to_num(masker.transform(img).ravel().astype(float))
    return out


def pooled(members, g_maps, var_maps):
    stack = np.array([g_maps[i] for i in members])
    weights = np.array([1.0 / np.maximum(var_maps[i], 1e-9) for i in members])
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def build(members, z_maps):
    studies = []
    for i in members:
        foci, height = report_peaks(z_maps[i], mask_bool, shape, zooms,
                                    scheme=SCHEME, focus=FOCUS)
        if not foci:
            continue
        meta = {"sample_sizes": [ASSUMED_N], "reporting_threshold": float(height)}
        studies.append({"id": f"c{i}", "name": f"c{i}", "metadata": meta, "analyses": [
            {"id": f"c{i}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(c) for c in nib.affines.apply_affine(
                     affine, np.asarray(ijk, dtype=float))],
                 "values": [{"kind": "Z", "value": float(v)}]} for ijk, v in foci]}]})
    if len(studies) < 2:
        return None
    return Studyset({"id": "a", "name": "a", "studies": studies}, target=None, mask=mask_img)


def localisation(est, truth, use):
    a, t = est[use], truth[use]
    if not np.isfinite(a).all() or np.allclose(a, a[0]):
        return np.nan, np.nan
    top = t >= np.percentile(t, 90)
    order = stats.rankdata(a)
    n_pos, n_neg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
           if n_pos and n_neg else np.nan)
    return float(stats.spearmanr(a, t)[0]), float(auc)


if __name__ == "__main__":
    z_maps = load("study-*-animal_*_z.nii.gz")
    g_maps = load("study-*-animal_*[0-9]_g.nii.gz")
    var_maps = load("study-*-animal_*_g_var.nii.gz")
    ids = sorted(set(z_maps) & set(g_maps) & set(var_maps))
    print(f"NeuroVault animal: {len(ids)} studies, {SCHEME}/{FOCUS} reporting, "
          f"{len(ids) // 2} per half over {N_SPLITS} splits, assumed N = {ASSUMED_N}")
    print("Scored on localisation only; the question is whether the ordering survives.\n")

    names = ("CBES g", "CBES g_marginal", "CBES prevalence", "ALE", "MKDA density", "KDA")
    rows = {k: [] for k in names}
    rng = np.random.default_rng(0)
    for _ in range(N_SPLITS):
        order = rng.permutation(len(ids))
        work = [ids[j] for j in order[: len(ids) // 2]]
        hold = [ids[j] for j in order[len(ids) // 2:]]
        truth = np.abs(pooled(hold, g_maps, var_maps))
        studyset = build(work, z_maps)
        if studyset is None:
            continue
        estimates = {}
        cbes = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
                    threshold="reporting_threshold")
        res = cbes.fit(studyset)
        g = np.abs(res.get_map("g", return_type="array").ravel())
        estimates["CBES g"] = g
        if "prevalence" in res.maps:
            pi = res.get_map("prevalence", return_type="array").ravel()
            estimates["CBES prevalence"] = pi
            estimates["CBES g_marginal"] = g * pi
        covered = res.get_map("n_studies", return_type="array").ravel() > 0
        for label, cls in (("ALE", ALE), ("MKDA density", MKDADensity), ("KDA", KDA)):
            out = cls(null_method="approximate", mask=masker).fit(studyset)
            estimates[label] = out.get_map("stat", return_type="array").ravel()
        use = covered & np.isfinite(truth)
        for label in names:
            if label not in estimates:
                continue
            rho, auc = localisation(estimates[label], truth, use)
            if np.isfinite(rho):
                rows[label].append((rho, auc))

    print(f"{'estimate':20s} {'rank r':>8s} {'(sd)':>7s} {'AUC':>8s} {'(sd)':>7s} {'splits':>7s}")
    for label in names:
        a = np.array(rows[label])
        if not a.size:
            print(f"{label:20s} {'--':>8s} {'--':>7s} {'--':>8s} {'--':>7s} {0:7d}")
            continue
        print(f"{label:20s} {a[:,0].mean():8.3f} {a[:,0].std(ddof=1):7.3f} "
              f"{a[:,1].mean():8.3f} {a[:,1].std(ddof=1):7.3f} {len(a):7d}")
    base = np.array(rows["MKDA density"])
    if base.size:
        print("\nPaired against MKDA density, the strongest convergence arm on pain:")
        for label in ("CBES g", "CBES g_marginal", "CBES prevalence", "ALE"):
            a = np.array(rows[label])
            if a.shape != base.shape or not a.size:
                continue
            for j, what in ((0, "rank r"), (1, "AUC   ")):
                d = a[:, j] - base[:, j]
                t = stats.ttest_rel(a[:, j], base[:, j])
                print(f"  {label:18s} {what}  {d.mean():+.3f}  (sd {d.std(ddof=1):.3f}, "
                      f"paired p {t.pvalue:.3f})")
    print("\nOn pain, g was level with MKDA (+0.001 AUC, p 0.980) and g_marginal beat it")
    print("(+0.111, p 0.001). The question here is whether that ordering is a property of")
    print("the method or of the pain literature.")
