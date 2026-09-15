"""Does the reported interval cover a known truth, and is it narrow enough to mean anything?

The estimator reports `g` and an `se`, and the docstring invites `g +/- 1.96*se`. Nothing had
ever checked whether that interval covers. Every earlier accuracy measurement was a correlation
or an RMSE against a reference, and neither says anything about calibration.

Two confounds are removed by construction.

*Estimand mismatch.* `g` is a conditional magnitude -- the effect among studies that have it --
while a truth built from all studies is marginal, so an interval could miss for a reason that
has nothing to do with its width. Every study here has the effect at every site, so prevalence
is 1 and the two coincide. What a study *reports* is still selected by the threshold, which is
the point; what it *has* is not.

*Statistic convention.* Studies report a genuine **t**, built by ``reporting.study_t_field``
and declared as ``"T"``. An earlier version of this bed reported ``(truth + noise/sqrt(n)) *
sqrt(n)``, a normal statistic with known variance, while the estimator -- correctly for real
data -- reads a reported statistic as a t on ``n - 1`` degrees of freedom. In the far tail where
every reported peak lives that mismatch inflated the recovered effect size by about a third, and
the inflation was read as a property of the estimator. The bed now asserts its own convention
before measuring anything.

*Regime.* The bed is calibrated first (`calibrate_coverage_bed`) to the regime coordinate
meta-analysis draws from: ~100% of studies report something and a typical table lists 3-4
clusters. A first attempt at this run, set up by guess, sat at 18% reporting -- a nominally
12-study analysis was really a 2-study one, a third of replications could not be fitted, and the
survivors were selected for high signal, biasing the very quantity being measured.

Four numbers are reported, because coverage alone cannot distinguish the failures:

  bias         mean(g) - truth. A point-estimate problem.
  se / sd      the reported width against the estimator's own variability across replications.
               Under 1 the interval is too narrow; well over 1 it overstates uncertainty.
  coverage     of g +/- 1.96*se. What a user experiences.
  half / truth the interval's half-width as a fraction of the effect. An interval wide enough
               to cover everything covers, and says nothing; coverage without this is a metric
               that lies.

Heterogeneity is varied because the estimator holds tau2 fixed at a value earlier work found
low-biased: if that matters, coverage should fall when the studies genuinely disagree.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from pathlib import Path
from scipy import ndimage
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
from nimare.studyset import Studyset
import reporting

SHAPE, VOXEL_MM = (30, 30, 30), 4.0
SMOOTH_VOX = 0.8           # calibrated: 7.6 mm FWHM, a realistic reported smoothness
PEAK_G = 0.8               # calibrated: 100% of studies report, 3.7 foci each
RADIUS_VOX = 2.5
N_SIMS = int(os.environ.get("NSIMS", 100))
#: Per-process by default, because two concurrent runs of this script generated identical image
#: filenames from (seed, studies, images, tau, k) and overwrote each other's files. That did not
#: fail loudly -- it silently dropped arms from the output table, which is how it was noticed.
OUT = Path(os.environ.get("COVOUT", f"/tmp/claude-0/cov_imgs_{os.getpid()}"))
OUT.mkdir(parents=True, exist_ok=True)

AFF = np.eye(4)
AFF[:3, :3] *= VOXEL_MM
AFF[:3, 3] = -VOXEL_MM * (np.array(SHAPE) - 1) / 2.0
MASK = nib.Nifti1Image(np.ones(SHAPE, dtype=np.int32), AFF)
MASK_BOOL = np.ones(SHAPE, dtype=bool)
ZOOMS = np.full(3, VOXEL_MM)

SITE_IJK = [(10, 15, 15), (20, 15, 15), (15, 9, 18), (15, 21, 12), (15, 15, 22)]
SITE_WEIGHT = [1.0, 1.0, 0.75, 0.75, 0.5]
#: Coverage is read at the strongest site, where the truth is PEAK_G exactly.
READ_AT = SITE_IJK[0]


def truth_field():
    out = np.zeros(SHAPE)
    grid = np.indices(SHAPE).astype(float)
    for (i, j, k), w in zip(SITE_IJK, SITE_WEIGHT):
        d2 = (grid[0] - i) ** 2 + (grid[1] - j) ** 2 + (grid[2] - k) ** 2
        out = np.maximum(out, PEAK_G * w * np.exp(-d2 / (2 * RADIUS_VOX**2)))
    return out


TRUTH = truth_field()
TRUE_G = float(TRUTH[READ_AT])


def study_fields(rng, n, tau):
    """One study's ``(g map, var map, t map)``: the image it would share and the table it prints.

    The image is on the effect-size scale, which is what a shared ``g``/``g_var`` pair carries.
    The statistic is a genuine t on ``n - 1`` degrees of freedom, which is what a paper prints
    and what the estimator's conversion assumes. Both come from the same study-level truth, so
    a donor's image and its table describe the same study.
    """
    shift = rng.normal(0.0, tau) if tau > 0 else 0.0
    field = TRUTH * (1.0 + shift / PEAK_G)      # the whole pattern scales, keeping zeros zero
    t_map = reporting.study_t_field(field, n, SMOOTH_VOX, rng, shape=SHAPE)
    # The image the same study would share is the effect-size map implied by its own t, through
    # the estimator's own conversion -- so the image arm and the coordinate arm are on exactly
    # one convention and the two are not independent draws of the same study. Building it as
    # `t / sqrt(n)` instead is Cohen's d, which is high by the Hedges factor (2.6% at n = 30)
    # and showed up as a residual bias in the all-donor reference arm.
    flat = t_map.ravel()
    g_flat, var_flat = peak_stat_to_hedges_g(
        flat, np.full(flat.size, float(n)), stat_type="t", design="one-sample"
    )
    return g_flat.reshape(SHAPE), var_flat.reshape(SHAPE), t_map


def one(seed, n_studies, n_image, tau, peak_bias=None, fwhm=10.0,
        peak_bias_scale=1.0):
    """Every study publishes a coordinate table; the first `n_image` also supply images.

    CBES refuses a collection with images and no coordinates, so an images-only arm does not
    exist to be measured. The refusal is right -- that would be an IBMA -- but it means the
    contrast here is "how much do images narrow the interval", not "images versus coordinates".
    A study that shares its map also published its table, so both are given.
    """
    rng = np.random.default_rng(seed)
    studies, reported = [], 0
    for k in range(n_studies):
        n = int(rng.integers(20, 41))
        gmap, var_map, t_map = study_fields(rng, n, tau)
        foci, _ = reporting.report_peaks(t_map[MASK_BOOL], MASK_BOOL, SHAPE,
                                         ZOOMS, scheme="cluster", focus="max")
        reported += bool(foci)
        meta = {"sample_sizes": [n]}
        points = [{"space": "MNI",
                   "coordinates": [float(v) for v in nib.affines.apply_affine(AFF, ijk)],
                   "values": [{"kind": "T", "value": float(zv)}]} for ijk, zv in foci]
        analysis = {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}
        if k < n_image:
            tag = f"{seed}_{n_studies}_{n_image}_{tau:.2f}_{k}"
            gp, vp = OUT / f"{tag}_g.nii.gz", OUT / f"{tag}_v.nii.gz"
            nib.save(nib.Nifti1Image(gmap.astype(np.float32), AFF), gp)
            nib.save(nib.Nifti1Image(var_map.astype(np.float32), AFF), vp)
            analysis["images"] = [
                {"url": str(gp), "filename": "g.nii.gz", "space": "MNI", "value_type": "g"},
                {"url": str(vp), "filename": "v.nii.gz", "space": "MNI", "value_type": "g_var"}]
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                        "analyses": [analysis]})

    if reported < 2:
        return None
    ss = Studyset({"id": "cov", "name": "cov", "studies": studies}, target=None, mask=MASK)
    est = CBES(fwhm=fwhm, mask=MASK, null_method="none", use_images=n_image > 0,
               peak_bias=peak_bias, peak_bias_scale=peak_bias_scale)
    res = est.fit(ss)
    pos = int(np.ravel_multi_index(READ_AT, SHAPE))
    pi = (float(res.get_map("prevalence", return_type="array").ravel()[pos])
          if "prevalence" in res.maps else np.nan)
    return (float(res.get_map("g", return_type="array").ravel()[pos]),
            float(res.get_map("se", return_type="array").ravel()[pos]), pi, reported)


#: Kernel widths to sweep in the coordinates-only arm. Widening the kernel was measured on real
#: data to improve both the accuracy of the estimate and its spatial extent (see the kernel-width
#: note in cbes-open-program); whether it also repairs the *interval* is a separate question,
#: since the coverage failure here is entirely bias rather than width.
FWHM_SWEEP = (10.0, 16.0, 24.0)

#: The configuration the docstring recommends for a mixed collection: the per-study correction
#: with its overall scale read off the image studies. Every arm below with `peak_bias=None`
#: measures only the *dilution* effect of images -- them contributing unbiased values alongside
#: the coordinates -- and not the *calibration* effect, where they pin the coordinate arm's own
#: scale. Those are different mechanisms and the weight-share model describes only the first.
CALIBRATED = ("per-study", "images")

#: (label, studies, of which supplying images, tau, peak_bias)
#: The `per-study` rows test the remedy the docstring itself recommends when no images are
#: available, so that the coordinates-only failure is measured against its own best defence
#: rather than against a configuration nobody is advised to use.
ARMS = [
    ("12 studies,  0 images, tau 0.0", 12, 0, 0.0, None),
    ("12 studies,  0 images, per-study", 12, 0, 0.0, "per-study"),
    ("12 studies,  2 images, tau 0.0", 12, 2, 0.0, None),
    ("12 studies,  6 images, tau 0.0", 12, 6, 0.0, None),
    ("12 studies, 12 images, tau 0.0", 12, 12, 0.0, None),
    ("12 studies,  0 images, tau 0.3", 12, 0, 0.3, None),
    ("12 studies,  6 images, tau 0.3", 12, 6, 0.3, None),
    ("12 studies, 12 images, tau 0.3", 12, 12, 0.3, None),
    ("24 studies,  0 images, tau 0.0", 24, 0, 0.0, None),
    ("24 studies,  0 images, per-study", 24, 0, 0.0, "per-study"),
    ("24 studies,  6 images, tau 0.0", 24, 6, 0.0, None),
    ("24 studies, 24 images, tau 0.0", 24, 24, 0.0, None),
]

if __name__ == "__main__":
    # The convention check this bed exists to respect, run on a null field before anything is
    # measured. A t reported as a z, or the reverse, shifts every magnitude by about a third.
    _probe = reporting.study_t_field(np.zeros(SHAPE), 30, SMOOTH_VOX,
                                     np.random.default_rng(12345), shape=SHAPE)
    reporting.assert_statistic_convention(_probe, 30, "T")
    print("statistic convention check passed: studies report a t on n - 1 degrees of freedom")
    print(f"truth at the read-out voxel: g = {TRUE_G:.3f}; prevalence 1 at every site")
    print(f"{N_SIMS} replications per arm; interval is g +/- 1.96*se\n")
    print(f"{'arm':34s} {'mean g':>7s} {'bias':>7s} {'mean se':>8s} {'sd of g':>8s} "
          f"{'se/sd':>6s} {'cover':>6s} {'half/truth':>10s} {'mean pi':>8s} {'report':>7s}")
    arms = [(f"{lab} @ fwhm {w:.0f}" if w != 10.0 else lab, ns, ni, tau, pb, w, 1.0)
            for lab, ns, ni, tau, pb in ARMS
            for w in (FWHM_SWEEP if (ni == 0 and tau == 0.0 and pb is None) else (10.0,))]
    # The recommended mixed configuration, on the arms where it can do anything: it needs both
    # images to read the scale off and coordinates for that scale to apply to.
    arms += [(f"{ns} studies, {ni:2d} images, calibrated", ns, ni, 0.0,
              CALIBRATED[0], 10.0, CALIBRATED[1])
             for ns, ni in ((12, 2), (12, 6), (24, 6))]
    for label, ns, ni, tau, pb, width, pbs in arms:
        rows = Parallel(n_jobs=8)(delayed(one)(s, ns, ni, tau, pb, width, pbs)
                                  for s in range(N_SIMS))
        refused = sum(r is None for r in rows)
        rows = [r for r in rows if r is not None]
        if not rows:
            print(f"{label:34s}   no usable fits ({refused} refused)")
            continue
        g = np.array([r[0] for r in rows]); se = np.array([r[1] for r in rows])
        pi = np.array([r[2] for r in rows]); rep = np.array([r[3] for r in rows])
        ok = np.isfinite(g) & np.isfinite(se) & (se > 0)
        cover = np.mean((g[ok] - 1.96 * se[ok] <= TRUE_G) & (TRUE_G <= g[ok] + 1.96 * se[ok]))
        sd = g[ok].std(ddof=1)
        print(f"{label:34s} {g[ok].mean():7.3f} {g[ok].mean()-TRUE_G:+7.3f} {se[ok].mean():8.3f} "
              f"{sd:8.3f} {se[ok].mean()/max(sd,1e-9):6.2f} {cover:6.2f} "
              f"{1.96*se[ok].mean()/TRUE_G:10.2f} {np.nanmean(pi):8.3f} "
              f"{rep.mean()/ns:7.2f}   (n={ok.sum()}, {refused} refused)", flush=True)
    print("\nse/sd near 1 means the width matches the estimator's real variability.")
    print("coverage near 0.95 with half/truth well under 1 is the only combination that works:")
    print("a wide interval covers by covering everything, and a biased one misses however wide.")
