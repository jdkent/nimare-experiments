"""Does CBES's magnitude map localise better than a plain convergence map? (task #38)

The question a reviewer will ask. CBES carries a censored mixture likelihood, a selection model,
a prevalence parameter and an effect-size conversion; ALE and MKDA carry a smoothing kernel and a
null. If the extra machinery does not localise the effect better, it is being paid for with
complexity and getting nothing back on the use a meta-analytic map is actually put to.

Note this is *not* the same question as the height-flattening test. Flattening gives every focus
the same magnitude and lets the likelihood run, which produces a nearly constant map -- the
variation had nowhere to go. A convergence statistic is a kernel-weighted density of foci, which
varies spatially by construction. So the two tests ask different things and the first cannot
answer this one.

Design follows the protocol: the 21-study NIDM pain collection split in half, one half's maps
extracted the way a paper would print them (multiplicity-corrected, whole surviving clusters, one
focus each, nothing capped, the real threshold handed in), the other half's inverse-variance
pooling as the reference. So no study is on both sides, and every estimator sees exactly the same
tables.

Scored on localisation, which is the comparable quantity: a convergence statistic has no
effect-size scale, so correlating its magnitude with the truth's would be a category error. Rank
correlation and the AUC for the truth's top decile are scale-free and are what a reader uses the
map for.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
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
from nimare.transforms import d_to_g, t_to_d, ImageTransformer
from load_pain import load_pain
from reporting import report_peaks

SCHEME, FOCUS = "cluster", "max"
N_SPLITS = int(os.environ.get("NSPLITS", 10))

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def pooled(members, maps, sizes):
    stack, weights = [], []
    for i in members:
        g = to_g(maps[i], sizes[i])
        stack.append(g)
        weights.append(1.0 / np.maximum(1.0 / sizes[i] + g**2 / (2.0 * sizes[i]), 1e-9))
    stack, weights = np.array(stack), np.array(weights)
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def build(members, maps, sizes):
    """One studyset the way a paper would have published it, shared by every estimator."""
    studies = []
    for i in members:
        foci, height = report_peaks(maps[i], mask_bool, shape, zooms,
                                    scheme=SCHEME, focus=FOCUS)
        if not foci:
            continue
        meta = {"sample_sizes": [int(sizes[i])], "reporting_threshold": float(height)}
        studies.append({"id": f"c{i}", "name": f"c{i}", "metadata": meta, "analyses": [
            {"id": f"c{i}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(c) for c in nib.affines.apply_affine(
                     affine, np.asarray(ijk, dtype=float))],
                 "values": [{"kind": "Z", "value": float(v)}]} for ijk, v in foci]}]})
    if len(studies) < 2:
        return None
    return Studyset({"id": "v", "name": "v", "studies": studies}, target=None, mask=mask_img)


def localisation(est, truth, use):
    """Rank correlation and top-decile AUC -- the scale-free scores the comparison allows."""
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
    ss = ImageTransformer(target="z").transform(load_pain())
    maps, sizes = [], []
    for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
        path = getattr(row, "z", None)
        if path is None or not os.path.isfile(str(path)) or not np.isfinite(float(n)):
            continue
        img = resample_to_img(nib.load(str(path)), mask_img, interpolation="continuous",
                              force_resample=True, copy_header=True)
        maps.append(np.nan_to_num(masker.transform(img).ravel().astype(float)))
        sizes.append(float(n))
    sizes = np.asarray(sizes, dtype=float)
    total = len(maps)
    print(f"NIDM pain: {total} studies, {SCHEME}/{FOCUS} reporting, "
          f"{total // 2} per half over {N_SPLITS} splits")
    print("Scored on localisation only: a convergence statistic has no effect-size scale.\n")

    names = ("CBES g", "CBES g_marginal", "CBES prevalence", "ALE", "MKDA density", "KDA")
    rows = {k: [] for k in names}
    rng = np.random.default_rng(0)
    for _ in range(N_SPLITS):
        order = rng.permutation(total)
        work, hold = order[: total // 2], order[total // 2:]
        truth = np.abs(pooled(hold, maps, sizes))
        studyset = build(list(work), maps, sizes)
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

    # The splits share studies, so the arms are paired and a paired test is the right one --
    # the between-split variance is common to all of them and would swamp an unpaired
    # comparison. The reference for the comparison is the best convergence statistic.
    print("\nPaired against MKDA density, the strongest convergence arm, across splits:")
    base = np.array(rows["MKDA density"])
    for label in ("CBES g", "CBES g_marginal", "CBES prevalence", "ALE"):
        a = np.array(rows[label])
        if a.shape != base.shape or not a.size:
            continue
        for j, what in ((0, "rank r"), (1, "AUC   ")):
            d = a[:, j] - base[:, j]
            t = stats.ttest_rel(a[:, j], base[:, j])
            print(f"  {label:18s} {what}  {d.mean():+.3f}  (sd {d.std(ddof=1):.3f}, "
                  f"paired p {t.pvalue:.3f})")
    print("\nIf a convergence statistic localises as well as the magnitude map, the censored")
    print("likelihood is not earning its complexity on the use a reader puts the map to.")
