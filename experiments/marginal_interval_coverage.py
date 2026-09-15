"""Does the interval on `g_marginal` cover pi*mu, and does the interval on `g` cover mu?

`se_marginal` was added to the estimator this session and validated against a numerical Hessian
of the observed log-likelihood -- so the arithmetic of the delta method through the 2x2 inverse is
right. What was never measured is whether the resulting interval *covers*. Validating a standard
error against a Hessian and validating an interval against a truth are different claims, and only
the second is what a user experiences.

The existing coverage bed cannot ask this. It sets prevalence to 1 at every site deliberately, so
that `g` and `g_marginal` coincide and an interval could not miss for reasons of estimand rather
than width. That made it the right bed for `g` and a useless one for `g_marginal`.

Here every study independently *has* the effect at a site with probability pi_v, and has no effect
there otherwise. Two truths follow and they are far apart at the read-out site:

    mu      = 0.80   the effect among studies that have it            <- what `g` estimates
    pi      = 0.50   the fraction of studies that have it             <- what `prevalence` does
    pi * mu = 0.40   the average effect over all studies              <- what `g_marginal` does

A study that does not have the effect at a site is not silent about it by choice -- it genuinely
has nothing there -- which is exactly the zero-inflation the model posits, rather than the
censoring. Both are present: a study that does have the effect still only reports a peak if its
own realisation survives its own threshold.

Nothing caps the number of peaks. Every study reports whatever survives a cluster-forming cut,
which is the only reporting rule these beds use.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
from nimare.studyset import Studyset
import interval_coverage as bed
import reporting

N = int(os.environ.get("NSIMS", 80))
RADIUS_VOX = 2.5

#: Site magnitudes and prevalences. The read-out site is the first: a large effect present in
#: only half the studies, which is where mu and pi*mu are furthest apart.
SITE_IJK = [(10, 15, 15), (20, 15, 15), (15, 9, 18), (15, 21, 12), (15, 15, 22)]
SITE_MU = [0.80, 0.80, 0.60, 0.60, 0.80]
SITE_PI = [0.50, 1.00, 0.75, 0.50, 1.00]
READ_SITE = 0
TRUE_MU = SITE_MU[READ_SITE]
TRUE_PI = SITE_PI[READ_SITE]
TRUE_MARGINAL = TRUE_MU * TRUE_PI
READ_AT = SITE_IJK[READ_SITE]

_GRID = np.indices(bed.SHAPE).astype(float)


def _blob(centre, height):
    i, j, k = centre
    d2 = (_GRID[0] - i) ** 2 + (_GRID[1] - j) ** 2 + (_GRID[2] - k) ** 2
    return height * np.exp(-d2 / (2 * RADIUS_VOX**2))


def study_truth(rng):
    """One study's own effect map: each site present independently with its own prevalence."""
    out = np.zeros(bed.SHAPE)
    present = []
    for centre, mu, pi in zip(SITE_IJK, SITE_MU, SITE_PI):
        has = rng.random() < pi
        present.append(has)
        if has:
            out = np.maximum(out, _blob(centre, mu))
    return out, present


def one(seed, n_studies, n_image):
    rng = np.random.default_rng(seed)
    studies, had = [], 0
    for k in range(n_studies):
        n = int(rng.integers(20, 41))
        truth, present = study_truth(rng)
        had += bool(present[READ_SITE])
        t_map = reporting.study_t_field(truth, n, bed.SMOOTH_VOX, rng, shape=bed.SHAPE)
        foci, _ = reporting.report_peaks(t_map[bed.MASK_BOOL], bed.MASK_BOOL, bed.SHAPE,
                                         bed.ZOOMS, scheme="cluster", focus="max")
        meta = {"sample_sizes": [n]}
        points = [{"space": "MNI",
                   "coordinates": [float(v) for v in nib.affines.apply_affine(bed.AFF, ijk)],
                   "values": [{"kind": "T", "value": float(zv)}]} for ijk, zv in foci]
        analysis = {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}
        if k < n_image:
            flat = t_map.ravel()
            g_flat, var_flat = peak_stat_to_hedges_g(
                flat, np.full(flat.size, float(n)), stat_type="t", design="one-sample")
            tag = f"mg{seed}_{n_studies}_{n_image}_{k}"
            gp, vp = bed.OUT / f"{tag}_g.nii.gz", bed.OUT / f"{tag}_v.nii.gz"
            nib.save(nib.Nifti1Image(g_flat.reshape(bed.SHAPE).astype(np.float32), bed.AFF), gp)
            nib.save(nib.Nifti1Image(var_flat.reshape(bed.SHAPE).astype(np.float32), bed.AFF), vp)
            analysis["images"] = [
                {"url": str(gp), "filename": "g.nii.gz", "space": "MNI", "value_type": "g"},
                {"url": str(vp), "filename": "v.nii.gz", "space": "MNI", "value_type": "g_var"}]
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                        "analyses": [analysis]})
    ss = Studyset({"id": "mg", "name": "mg", "studies": studies}, target=None, mask=bed.MASK)
    est = CBES(fwhm=10.0, mask=bed.MASK, null_method="none", use_images=n_image > 0,
               peak_bias="per-study", peak_bias_scale="images" if n_image else 1.0)
    res = est.fit(ss)
    pos = int(np.ravel_multi_index(READ_AT, bed.SHAPE))

    def at(name):
        if name not in res.maps:
            return np.nan
        return float(res.get_map(name, return_type="array").ravel()[pos])

    return (at("g"), at("se"), at("g_marginal"), at("se_marginal"), at("prevalence"),
            had / n_studies)


def _summarise(label, rows, idx_est, idx_se, truth):
    est = np.array([r[idx_est] for r in rows], dtype=float)
    se = np.array([r[idx_se] for r in rows], dtype=float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    if not ok.any():
        print(f"  {label:34s}  no usable fits")
        return
    cover = np.mean((est[ok] - 1.96*se[ok] <= truth) & (truth <= est[ok] + 1.96*se[ok]))
    sd = est[ok].std(ddof=1)
    print(f"  {label:34s} {truth:7.3f} {est[ok].mean():8.3f} {est[ok].mean()-truth:+8.3f} "
          f"{se[ok].mean():8.3f} {sd:7.3f} {se[ok].mean()/sd:6.2f} {cover:6.2f}  (n={ok.sum()})",
          flush=True)


if __name__ == "__main__":
    _probe = reporting.study_t_field(np.zeros(bed.SHAPE), 30, bed.SMOOTH_VOX,
                                     np.random.default_rng(2024), shape=bed.SHAPE)
    reporting.assert_statistic_convention(_probe, 30, "T")
    print("statistic convention check passed: studies report a t on n - 1 df")
    print(f"read-out site: mu {TRUE_MU:.2f}, pi {TRUE_PI:.2f}, pi*mu {TRUE_MARGINAL:.2f}; "
          f"{N} replications\n")
    for n_studies, n_image in ((12, 2), (24, 2), (12, 0)):
        rows = [r for r in Parallel(n_jobs=4)(
            delayed(one)(s, n_studies, n_image) for s in range(N)) if r is not None]
        share = np.mean([r[5] for r in rows])
        print(f"{n_studies} studies, {n_image} images "
              f"(effect actually present in {share:.2f} of studies)")
        print(f"  {'estimand':34s} {'truth':>7s} {'mean':>8s} {'bias':>8s} {'mean se':>8s} "
              f"{'sd':>7s} {'se/sd':>6s} {'cover':>6s}")
        _summarise("g vs mu", rows, 0, 1, TRUE_MU)
        _summarise("g_marginal vs pi*mu", rows, 2, 3, TRUE_MARGINAL)
        # The wrong pairing, reported on purpose: an estimand mismatch and a width problem look
        # identical in a coverage number, and this is what the difference costs.
        _summarise("g vs pi*mu (mismatched)", rows, 0, 1, TRUE_MARGINAL)
        _summarise("g_marginal vs mu (mismatched)", rows, 2, 3, TRUE_MU)
        prev = np.array([r[4] for r in rows], dtype=float)
        print(f"  prevalence: mean {np.nanmean(prev):.3f} against a true {TRUE_PI:.2f} "
              f"and a realised {share:.2f}\n", flush=True)
