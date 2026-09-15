"""What makes an all-image CBES fit land 2% low when the images themselves are unbiased?

The coverage table's all-image arms sit at -0.018 (12 studies) and -0.022 (24 studies) on a
truth of 0.800, and every calibrated row is read against that floor. Two candidates have been
eliminated:

  * the bed's own conversion -- donor images average +0.003 at the read-out voxel, measured
    across 400 draws at each of n = 20/30/40/60 without CBES in the loop at all;
  * the spatial kernel -- swept over 4, 6, 10 and 16 mm the arm returns *identical* numbers to
    four decimals, because with every study donating an image the kernel has no foci to spread.
    The parameter did not move the thing, so the peak-attenuation story was untestable that way
    rather than merely wrong.

What is left is the weighting. Hedges' variance is ``1/n + g^2 / (2(n - 1))``, a function of the
observed effect, so a voxel that drew high gets a larger variance and less inverse-variance
weight than one that drew low. That is a downward bias by construction, and 2% is the right
order at n ~ 30 for g ~ 0.8: the g^2 term is 0.011 against 1/n = 0.033, a third of the weight
varying with the square of a noisy quantity.

The test substitutes a variance that does *not* depend on the draw -- the same ``1/n`` for every
voxel of a donor -- and changes nothing else. If the floor goes away the mechanism is identified,
and it is not a property of this bed: real ``g_var`` maps are built with the g^2 term, so any
image-based pooling that uses them inherits the same bias.
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


def one(seed, n_studies, variance):
    """All studies donate images. ``variance`` is "hedges" (as reported) or "flat" (1/n)."""
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(n_studies):
        n = int(rng.integers(20, 41))
        gmap, var_map, t_map = bed.study_fields(rng, n, 0.0)
        if variance == "flat":
            # The part of Hedges' variance that does not depend on the draw. Everything else
            # about the study -- its effect map, its peaks, its sample size -- is untouched.
            var_map = np.full_like(var_map, 1.0 / n)
        foci, _ = reporting.report_peaks(t_map[bed.MASK_BOOL], bed.MASK_BOOL, bed.SHAPE,
                                         bed.ZOOMS, scheme="cluster", focus="max")
        meta = {"sample_sizes": [n]}
        points = [{"space": "MNI",
                   "coordinates": [float(v) for v in nib.affines.apply_affine(bed.AFF, ijk)],
                   "values": [{"kind": "T", "value": float(zv)}]} for ijk, zv in foci]
        tag = f"fl{seed}_{n_studies}_{variance}_{k}"
        gp, vp = bed.OUT / f"{tag}_g.nii.gz", bed.OUT / f"{tag}_v.nii.gz"
        nib.save(nib.Nifti1Image(gmap.astype(np.float32), bed.AFF), gp)
        nib.save(nib.Nifti1Image(var_map.astype(np.float32), bed.AFF), vp)
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points, "images": [
                {"url": str(gp), "filename": "g.nii.gz", "space": "MNI", "value_type": "g"},
                {"url": str(vp), "filename": "v.nii.gz", "space": "MNI",
                 "value_type": "g_var"}]}]})
    ss = Studyset({"id": "fl", "name": "fl", "studies": studies}, target=None, mask=bed.MASK)
    est = CBES(fwhm=10.0, mask=bed.MASK, null_method="none", use_images=True,
               peak_bias="per-study", peak_bias_scale="images")
    res = est.fit(ss)
    pos = int(np.ravel_multi_index(bed.READ_AT, bed.SHAPE))
    return (float(res.get_map("g", return_type="array").ravel()[pos]),
            float(res.get_map("se", return_type="array").ravel()[pos]))


if __name__ == "__main__":
    print(f"truth {bed.TRUE_G:.3f}; every study donates an image; {N} replications")
    print(f"{'studies':>7s} {'variance':>9s} {'mean g':>8s} {'bias':>8s} {'sem':>7s} "
          f"{'mean se':>8s} {'sd':>7s} {'cover':>6s}")
    for n_studies in (12, 24):
        for variance in ("hedges", "flat"):
            rows = [r for r in Parallel(n_jobs=4)(
                delayed(one)(s, n_studies, variance) for s in range(N)) if r is not None]
            g = np.array([r[0] for r in rows]); se = np.array([r[1] for r in rows])
            ok = np.isfinite(g) & np.isfinite(se) & (se > 0)
            cover = np.mean((g[ok] - 1.96*se[ok] <= bed.TRUE_G)
                            & (bed.TRUE_G <= g[ok] + 1.96*se[ok]))
            sd = g[ok].std(ddof=1)
            print(f"{n_studies:7d} {variance:>9s} {g[ok].mean():8.3f} "
                  f"{g[ok].mean()-bed.TRUE_G:+8.3f} {sd/np.sqrt(ok.sum()):7.4f} "
                  f"{se[ok].mean():8.3f} {sd:7.3f} {cover:6.2f}  (n={ok.sum()})", flush=True)
    print("\nPaired on seeds, so the two variance rows see the same studies. The floor is")
    print("identified only if 'flat' moves it by more than a couple of standard errors.")
