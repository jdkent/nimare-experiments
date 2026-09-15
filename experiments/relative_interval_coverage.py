"""Does an interval on `g_relative` cover, in the regime where the interval on `g` does not?

The coverage table's verdict on `g` is that coordinates alone leave a bias of about 30% of the
effect that no interval width absorbs -- 0.75 coverage at twelve studies, 0.35 at twenty-four --
and that only image donors pinning the peak-height scale repair it. But the docstring already
says `g_relative` is the map to read, and `g_relative` is built precisely to be immune to the
thing that breaks `g`:

    g_relative = g / P95(|g| over covered voxels)

If the peak-height inflation were exactly one common multiplicative factor `c`, then g = c*mu
everywhere and the factor cancels *exactly* in that ratio. So the estimand `g_relative` targets is
`mu / P95(|mu|)`, and the bias that destroys the interval on `g` should not appear in it. Nobody
has checked, and no standard error for it is emitted.

To first order the error is `se / P95(|g|)`: the denominator is a map-wide quantile over thousands
of voxels, estimated far more precisely than any single voxel's `g`, so treating it as fixed
should be adequate. That is the candidate `se_relative`, and it is one line -- but only worth
proposing if the interval it builds actually covers.

Three things are measured, because there are three ways this can fail:

  bias      does `g_relative` hit `mu / P95(|mu|)`? If the inflation is not a pure common factor
            the cancellation is incomplete and this is non-zero.
  se/sd     is `se / P95(|g|)` the right width for the spread `g_relative` actually has? The
            delta-method term for the randomness of the denominator is being dropped, and if that
            term matters this comes out below 1.
  coverage  the thing a user experiences.

The comparison against `g` in the same fits is the point: same collections, same foci, same
replications, so any difference is the estimand and not the bed.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
import interval_coverage as bed
import reporting

N = int(os.environ.get("NSIMS", 100))
PCT = 95


def true_relative():
    """`mu / P95(|mu|)` at the read-out voxel, over the voxels a fit can cover.

    The truth field is mostly zero, so the 95th percentile of its magnitude is taken over the
    same whole-mask voxel set the estimator normalizes over. Using a different voxel set here
    than the estimator uses would make the comparison meaningless -- that is the mistake task
    #38 recorded as "it depends entirely on the voxel set".
    """
    magnitude = np.abs(bed.TRUTH.ravel())
    reference = float(np.percentile(magnitude, PCT))
    return float(bed.TRUTH[bed.READ_AT] / reference), reference


def one(seed, n_studies, n_image):
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(n_studies):
        n = int(rng.integers(20, 41))
        gmap, var_map, t_map = bed.study_fields(rng, n, 0.0)
        foci, _ = reporting.report_peaks(t_map[bed.MASK_BOOL], bed.MASK_BOOL, bed.SHAPE,
                                         bed.ZOOMS, scheme="cluster", focus="max")
        meta = {"sample_sizes": [n]}
        points = [{"space": "MNI",
                   "coordinates": [float(v) for v in nib.affines.apply_affine(bed.AFF, ijk)],
                   "values": [{"kind": "T", "value": float(zv)}]} for ijk, zv in foci]
        analysis = {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}
        if k < n_image:
            tag = f"rl{seed}_{n_studies}_{n_image}_{k}"
            gp, vp = bed.OUT / f"{tag}_g.nii.gz", bed.OUT / f"{tag}_v.nii.gz"
            nib.save(nib.Nifti1Image(gmap.astype(np.float32), bed.AFF), gp)
            nib.save(nib.Nifti1Image(var_map.astype(np.float32), bed.AFF), vp)
            analysis["images"] = [
                {"url": str(gp), "filename": "g.nii.gz", "space": "MNI", "value_type": "g"},
                {"url": str(vp), "filename": "v.nii.gz", "space": "MNI", "value_type": "g_var"}]
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                        "analyses": [analysis]})
    ss = Studyset({"id": "rl", "name": "rl", "studies": studies}, target=None, mask=bed.MASK)
    est = CBES(fwhm=10.0, mask=bed.MASK, null_method="none", use_images=n_image > 0,
               peak_bias="per-study", peak_bias_scale="images" if n_image else 1.0)
    res = est.fit(ss)
    pos = int(np.ravel_multi_index(bed.READ_AT, bed.SHAPE))
    g = res.get_map("g", return_type="array").ravel()
    se = res.get_map("se", return_type="array").ravel()
    rel = res.get_map("g_relative", return_type="array").ravel()
    covered = res.get_map("n_studies", return_type="array").ravel() > 0
    magnitude = np.abs(g[covered])
    magnitude = magnitude[np.isfinite(magnitude)]
    if not magnitude.size:
        return None
    reference = float(np.percentile(magnitude, PCT))
    if not np.isfinite(reference) or reference <= 0:
        return None
    # The candidate se_relative: the reported error divided by the same normalizer the map used.
    return float(g[pos]), float(se[pos]), float(rel[pos]), float(se[pos]) / reference, reference


def _row(label, est, se, truth):
    est, se = np.asarray(est, float), np.asarray(se, float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    if not ok.any():
        print(f"  {label:30s}  no usable fits")
        return
    sd = est[ok].std(ddof=1)
    cover = np.mean((est[ok] - 1.96*se[ok] <= truth) & (truth <= est[ok] + 1.96*se[ok]))
    print(f"  {label:30s} {truth:7.3f} {est[ok].mean():8.3f} {est[ok].mean()-truth:+8.3f} "
          f"{se[ok].mean():8.3f} {sd:7.3f} {se[ok].mean()/sd:6.2f} {cover:6.2f}  (n={ok.sum()})",
          flush=True)


if __name__ == "__main__":
    truth_rel, truth_ref = true_relative()
    print(f"truth: mu {bed.TRUE_G:.3f}, P95(|mu|) {truth_ref:.3f}, "
          f"mu/P95 {truth_rel:.3f}; {N} replications\n")
    for n_studies, n_image in ((12, 0), (24, 0), (12, 2)):
        rows = [r for r in Parallel(n_jobs=4)(
            delayed(one)(s, n_studies, n_image) for s in range(N)) if r is not None]
        if not rows:
            print(f"{n_studies} studies, {n_image} images: no usable fits\n"); continue
        print(f"{n_studies} studies, {n_image} images  "
              f"(fitted P95 |g| = {np.mean([r[4] for r in rows]):.3f})")
        print(f"  {'estimand':30s} {'truth':>7s} {'mean':>8s} {'bias':>8s} {'se':>8s} "
              f"{'sd':>7s} {'se/sd':>6s} {'cover':>6s}")
        _row("g vs mu", [r[0] for r in rows], [r[1] for r in rows], bed.TRUE_G)
        _row("g_relative vs mu/P95(|mu|)", [r[2] for r in rows], [r[3] for r in rows], truth_rel)
        print()
    print("If the relative rows cover near nominal where the g rows do not, CBES has a validated")
    print("interval in the coordinates-only regime and `se_relative` is worth emitting. If they")
    print("do not, the inflation is not a common factor and the relative map inherits the bias.")
