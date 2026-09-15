"""Is the `se` gap arithmetic, or is it the local-maximum rule the model does not contain?

The model's censoring rule is `|g| < c` at a voxel: silence there means the effect is below the
cutoff. The simulator's reporting rule -- and a paper's -- is `local maximum of |z| AND above
threshold`. Those differ, and the difference is not small: a voxel on the shoulder of a blob can
sit far above the cutoff and still be silent, because it is not a local maximum. The model reads
that silence as evidence against the effect.

So the reported `se` may be exactly right about its own model and wrong about the data. This
distinguishes the two by generating data from the model's rule and from the simulator's, changing
nothing else:

  * **model-matched** -- report every voxel whose |g| clears the threshold. No local-maximum
    condition, so silence really does mean `|g| < c`. If the observed information is sound,
    se/sd must come out near 1 here.
  * **as published** -- the shipped rule, local maxima only.

This is not a claim that papers report whole blobs. It is an oracle: the only way to know whether
a se is mis-calibrated or merely mis-specified is to feed the estimator its own assumptions once.

No cap is applied in either arm -- the count is whatever clears the threshold.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tempfile
import numpy as np
import nibabel as nib
from nilearn.maskers import NiftiMasker
import nimare.generate as gen
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma.effectsize import CBES

ZOOMS, EXTENT, BLOB = 4.0, 48.0, 10.0
HALF = int(np.ceil(EXTENT / ZOOMS))
SHAPE = (2 * HALF + 1,) * 3
AFFINE = np.array([[ZOOMS, 0, 0, -ZOOMS * HALF], [0, ZOOMS, 0, -ZOOMS * HALF],
                   [0, 0, ZOOMS, -ZOOMS * HALF], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones(SHAPE, np.int32), AFFINE)
MASKER = NiftiMasker(MASK).fit()
FOCI = [(-28, -28, 0), (28, -28, 0), (-28, 28, 0), (28, 28, 0)]
EFFECTS = [0.2, 0.4, 0.6, 0.8]
N_REPS = int(os.environ.get("NREPS", 24))

grid = np.stack(np.indices(SHAPE), -1) * ZOOMS - ZOOMS * HALF
sd = BLOB / (2 * np.sqrt(2 * np.log(2)))
volume = np.zeros(SHAPE)
for focus, effect in zip(FOCI, EFFECTS):
    volume += effect * np.exp(
        -((grid - np.asarray(focus, float)) ** 2).sum(-1) / (2 * sd ** 2))
TRUTH = MASKER.transform(nib.Nifti1Image(volume.astype(np.float32), AFFINE)).ravel()

STRATA = [("null", TRUTH < 0.02)]
for focus, effect in zip(FOCI, EFFECTS):
    ijk = np.round(np.linalg.inv(AFFINE) @ np.array(list(focus) + [1.0]))[:3].astype(int)
    hit = np.zeros(SHAPE, bool)
    hit[tuple(ijk)] = True
    STRATA.append((f"g={effect}", MASKER.transform(
        nib.Nifti1Image(hit.astype(np.int16), AFFINE)).ravel() > 0))

# The shipped rule, captured so the patch can restore it.
_REAL_MAXIMUM_FILTER = None


def _identity_maximum_filter(magnitude, size=3):
    """Make the local-maximum test vacuous, so `is_peak` reduces to `magnitude >= threshold`."""
    return magnitude


def run(model_matched):
    import scipy.ndimage
    gs, ses, counts = [], [], []
    for rep in range(N_REPS):
        if model_matched:
            saved = scipy.ndimage.maximum_filter
            scipy.ndimage.maximum_filter = _identity_maximum_filter
        try:
            ss = create_effect_size_coordinate_studyset(
                FOCI, effect_sizes=EFFECTS, n_studies=20, sample_size=(20, 40), tau=0.1,
                seed=5000 + rep, simulate_field=True, n_image_studies=2,
                image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
                blob_fwhm=BLOB)
        finally:
            if model_matched:
                scipy.ndimage.maximum_filter = saved
        n_foci = sum(len(a.points) for s in ss.studies for a in s.analyses)
        counts.append(n_foci)
        try:
            res = CBES(mask=MASK, null_method="none",
                       threshold="reporting_threshold").fit(ss)
        except ValueError:
            continue
        gs.append(res.get_map("g", return_type="array").ravel())
        ses.append(res.get_map("se", return_type="array").ravel())
    return np.array(gs), np.array(ses), counts


if __name__ == "__main__":
    print(f"four foci at {EFFECTS}, {N_REPS} replications, 20 studies, 2 images, tau=0.1\n")
    head = (f"{'arm':>26} {'stratum':>8} {'bias':>7} {'sd(g)':>7} {'se':>7} "
            f"{'se/sd':>6} {'se/rmse':>8}")
    print(head)
    for label, matched in [("model-matched (|g|>=c)", True), ("as published (maxima)", False)]:
        G, SE, counts = run(matched)
        if G.shape[0] < min(5, N_REPS):
            print(f"{label:>26} {'--':>8}  too few usable replications")
            continue
        for name, sel in STRATA:
            g, se = G[:, sel], SE[:, sel]
            truth = TRUTH[sel]
            ok = np.isfinite(se) & (se > 0)
            bias = float(np.mean(g - truth))
            spread = float(np.mean(np.std(g, axis=0, ddof=1)))
            reported = float(np.mean(se[ok])) if ok.any() else np.nan
            rmse = float(np.sqrt(np.mean((g - truth) ** 2)))
            print(f"{label:>26} {name:>8} {bias:+7.3f} {spread:7.4f} {reported:7.4f} "
                  f"{reported / spread:6.2f} {reported / rmse:8.2f}")
        print(f"{'':>26} {'foci/studyset':>8} median {int(np.median(counts))}\n")
