"""Does the redesigned estimator still compress the dynamic range?

The old design's headline failure: judged against references the coordinates never touched, a
reference effect spanning elevenfold across its strata came back spanning about 1.2-fold. An
unknown overall scale would leave that ratio alone; it did not, so the compression was real and
not a units problem. It was never measured for the redesign.

Three-bin stratification cannot answer this -- it mixes the slope with the floor, and the top
bin here holds 7 voxels. So regress the estimate on the truth instead:

    estimate ~ a + b * truth

`b` is the slope the map recovers. b = 1 is calibrated, b < 1 is compression, and `a` is the
floor that survives where the truth is zero. Both are needed: an estimator can have b = 1 and
still be useless if a is large, and the absolute-value floor guarantees a > 0 for any map read
as |g|.

The truth range comes from four well-separated foci at 0.2, 0.4, 0.6 and 0.8, so the range is
inside one map and one fit rather than across collections -- which is how a reader would use it.
Scored on the voxels the simulator actually gave an effect, since a regression dominated by
9000 near-zero voxels estimates the floor and nothing else.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
import tempfile
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.maskers import NiftiMasker
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
N_SEEDS = int(os.environ.get("NSEEDS", 8))

grid = np.stack(np.indices(SHAPE), -1) * ZOOMS - ZOOMS * HALF
sd = BLOB / (2 * np.sqrt(2 * np.log(2)))
volume = np.zeros(SHAPE)
for focus, effect in zip(FOCI, EFFECTS):
    volume += effect * np.exp(
        -((grid - np.asarray(focus, float)) ** 2).sum(-1) / (2 * sd ** 2)
    )
TRUTH = MASKER.transform(nib.Nifti1Image(volume.astype(np.float32), AFFINE)).ravel()
SIGNAL = TRUTH >= 0.05


def fit_line(est):
    """Slope and intercept of estimate on truth, over the voxels carrying signal."""
    result = stats.linregress(TRUTH[SIGNAL], est[SIGNAL])
    return float(result.slope), float(result.intercept)


def peak_values(est):
    """The estimate at each focus's own voxel, so the four truths can be read off."""
    out = []
    for focus in FOCI:
        ijk = np.rint(np.linalg.solve(AFFINE, np.array([*focus, 1.0]))[:3]).astype(int)
        flat = np.ravel_multi_index(tuple(ijk), SHAPE)
        lookup = MASKER.transform(
            nib.Nifti1Image(np.arange(np.prod(SHAPE)).reshape(SHAPE).astype(np.float32), AFFINE)
        ).ravel()
        out.append(float(est[int(np.flatnonzero(lookup == flat)[0])]))
    return out


if __name__ == "__main__":
    arms = {"images only": [], "CBES g": [], "CBES g_marginal": []}
    for seed in range(N_SEEDS):
        ss = create_effect_size_coordinate_studyset(
            FOCI, effect_sizes=EFFECTS, n_studies=20, sample_size=(20, 40), tau=0.1,
            seed=seed, simulate_field=True, n_image_studies=2,
            image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
            blob_fwhm=BLOB)
        plain = CBES(mask=MASK, null_method="none", selection_model="none").fit(ss)
        arms["images only"].append(
            np.abs(plain.get_map("g", return_type="array").ravel()))
        full = CBES(mask=MASK, null_method="none", threshold="reporting_threshold").fit(ss)
        arms["CBES g"].append(np.abs(full.get_map("g", return_type="array").ravel()))
        arms["CBES g_marginal"].append(
            np.abs(full.get_map("g_marginal", return_type="array").ravel()))

    print(f"four foci at {EFFECTS}, {N_SEEDS} collections, 20 studies, 2 images")
    print(f"truth range over the signal voxels: {TRUTH[SIGNAL].min():.3f} to "
          f"{TRUTH[SIGNAL].max():.3f} ({TRUTH[SIGNAL].max() / TRUTH[SIGNAL].min():.1f}-fold), "
          f"{int(SIGNAL.sum())} voxels\n")
    print(f"{'estimate':>16} {'slope':>7} {'intercept':>10}    " +
          "  ".join(f"{'g@' + str(e):>18}" for e in EFFECTS))
    for label, maps in arms.items():
        slopes, intercepts = zip(*[fit_line(m) for m in maps])
        peaks = np.array([peak_values(m) for m in maps])
        cells = []
        for j, e in enumerate(EFFECTS):
            mean, sem = peaks[:, j].mean(), peaks[:, j].std(ddof=1) / np.sqrt(len(peaks))
            cells.append(f"{mean:.3f}+-{sem:.3f} ({100 * (mean - e) / e:+.0f}%)")
        print(f"{label:>16} {np.mean(slopes):7.3f} {np.mean(intercepts):+10.3f}    " +
              "  ".join(f"{c:>18}" for c in cells))
    print("\nrecovered range at the four foci, and mean |relative error|:")
    for label, maps in arms.items():
        peaks = np.array([peak_values(m) for m in maps]).mean(axis=0)
        err = np.mean(np.abs(100 * (peaks - np.array(EFFECTS)) / np.array(EFFECTS)))
        print(f"  {label:>16} {peaks[-1] / max(peaks[0], 1e-9):5.2f}-fold "
              f"(truth {EFFECTS[-1] / EFFECTS[0]:.0f}-fold)   mean |err| {err:5.1f}%")
