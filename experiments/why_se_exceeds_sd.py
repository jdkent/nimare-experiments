"""Why is the reported `se` larger than the estimator's own spread in almost every arm?

`se/sd` runs from 1.10 to 2.14 across the coverage table at tau = 0, and the one arm with real
heterogeneity is the only one near 1:

    12 studies, 6 images, tau 0     se/sd 1.26
    12 studies, 6 images, tau 0.3   se/sd 1.06

An interval wider than the estimator's variability is the *safe* direction, but "conservative" is
a description rather than an explanation, and it has been used as one several times in these
notes. `tau2` is not held fixed -- it is estimated per voxel by DerSimonian-Laird and truncated at
zero. A truncated estimator of a quantity whose true value is zero has a positive mean, so at
tau = 0 the fit should be charging the interval for heterogeneity that is not there, and the
charge should vanish once the truth moves away from the boundary.

That predicts two things, both read off the `tau2` map the estimator already emits:

  1. mean fitted tau2 > 0 when the true tau is 0, by enough to explain the excess width;
  2. at tau = 0.3 the fitted value should be near 0.09 and the excess width should mostly go.

And it predicts the remedy is not to widen or narrow anything, but that `tau2_method="none"`
should bring `se/sd` toward 1 at tau = 0 -- at the cost of being wrong whenever heterogeneity is
real, which is the honest trade rather than a fix.
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

N = int(os.environ.get("NSIMS", 60))


def one(seed, n_studies, n_image, tau, tau2_method):
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(n_studies):
        n = int(rng.integers(20, 41))
        gmap, var_map, t_map = bed.study_fields(rng, n, tau)
        foci, _ = reporting.report_peaks(t_map[bed.MASK_BOOL], bed.MASK_BOOL, bed.SHAPE,
                                         bed.ZOOMS, scheme="cluster", focus="max")
        meta = {"sample_sizes": [n]}
        points = [{"space": "MNI",
                   "coordinates": [float(v) for v in nib.affines.apply_affine(bed.AFF, ijk)],
                   "values": [{"kind": "T", "value": float(zv)}]} for ijk, zv in foci]
        analysis = {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}
        if k < n_image:
            tag = f"tw{seed}_{n_studies}_{n_image}_{tau:.2f}_{tau2_method}_{k}"
            gp, vp = bed.OUT / f"{tag}_g.nii.gz", bed.OUT / f"{tag}_v.nii.gz"
            nib.save(nib.Nifti1Image(gmap.astype(np.float32), bed.AFF), gp)
            nib.save(nib.Nifti1Image(var_map.astype(np.float32), bed.AFF), vp)
            analysis["images"] = [
                {"url": str(gp), "filename": "g.nii.gz", "space": "MNI", "value_type": "g"},
                {"url": str(vp), "filename": "v.nii.gz", "space": "MNI", "value_type": "g_var"}]
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                        "analyses": [analysis]})
    ss = Studyset({"id": "tw", "name": "tw", "studies": studies}, target=None, mask=bed.MASK)
    est = CBES(fwhm=10.0, mask=bed.MASK, null_method="none", use_images=n_image > 0,
               peak_bias="per-study", peak_bias_scale="images" if n_image else 1.0,
               tau2_method=tau2_method)
    res = est.fit(ss)
    pos = int(np.ravel_multi_index(bed.READ_AT, bed.SHAPE))
    tau2 = (float(res.get_map("tau2", return_type="array").ravel()[pos])
            if "tau2" in res.maps else np.nan)
    return (float(res.get_map("g", return_type="array").ravel()[pos]),
            float(res.get_map("se", return_type="array").ravel()[pos]), tau2)


if __name__ == "__main__":
    print(f"truth {bed.TRUE_G:.3f}; {N} replications; 12 studies, 6 image donors")
    print(f"{'true tau':>9s} {'tau2 est':>9s} {'fitted tau2':>12s} {'mean se':>8s} {'sd':>7s} "
          f"{'se/sd':>6s} {'bias':>8s} {'cover':>6s}")
    for tau in (0.0, 0.3):
        for method in ("dl", "none"):
            rows = [r for r in Parallel(n_jobs=4)(
                delayed(one)(s, 12, 6, tau, method) for s in range(N)) if r is not None]
            g = np.array([r[0] for r in rows]); se = np.array([r[1] for r in rows])
            t2 = np.array([r[2] for r in rows], dtype=float)
            ok = np.isfinite(g) & np.isfinite(se) & (se > 0)
            sd = g[ok].std(ddof=1)
            cover = np.mean((g[ok] - 1.96*se[ok] <= bed.TRUE_G)
                            & (bed.TRUE_G <= g[ok] + 1.96*se[ok]))
            print(f"{tau:9.2f} {method:>9s} {np.nanmean(t2):12.4f} {se[ok].mean():8.3f} "
                  f"{sd:7.3f} {se[ok].mean()/sd:6.2f} {g[ok].mean()-bed.TRUE_G:+8.3f} "
                  f"{cover:6.2f}  (n={ok.sum()}, true tau2 {tau**2:.4f})", flush=True)
    print("\nThe prediction is a positive fitted tau2 at true tau = 0 large enough to account for")
    print("the excess width, and se/sd near 1 in the tau2_method='none' row at tau = 0.")
