"""Validate CBES's magnitude against a truth measured on subjects it never saw.

Every reference used before this one was built from the same maps the coordinates came from,
so "this study has an effect at this voxel" was never observable independently of the effect's
size -- conditioning a study on clearing its own threshold selects the noise it then measures.
That is structural with one map per study, and no stratification fixes it.

HCP collection 4337 fixes it, because it has per-subject maps. Subjects are split in two:

  * one part is cut into synthetic studies of ``n_per_study``. Each study's own one-sample map
    is computed, thresholded, and its local maxima extracted -- the coordinates a paper would
    tabulate. CBES sees nothing else.
  * the other part, never used to make a coordinate, gives the truth at every voxel.

The selection that produces the peaks is then statistically independent of the quantity being
compared against, which is the whole point.

Units are consistent by construction. A synthetic study's t is ``mean / (sd / sqrt(n))`` over
its own subjects and its z is that mapped through the t distribution on ``n - 1`` df, which is
exactly the reporting model CBES inverts; its effect size and the held-out truth are both
Hedges' g on the subject-level scale.

Prevalence is 1 here: every synthetic study draws from one population, so mu is the marginal
and the pi/mu split is not under test. That is the regime where the corrected simulator put the
required scale at 1.05 +- 0.02, so this is a direct real-data test of that claim.

Coordinates come from :mod:`reporting`, which corrects each synthetic study's map for
multiplicity and keeps every surviving local maximum 8 mm apart without capping the count. The
earlier version of this script fixed an uncorrected height and took the strongest peaks, which
is not a reporting practice and which ties the effective cut to the signal.
"""
import os, sys, glob, warnings; warnings.simplefilter("ignore")
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
from nimare.transforms import t_to_z
from reporting import report_peaks

CONTRAST = sys.argv[1] if len(sys.argv) > 1 else "MOTOR_LH"
# Reporting schemes rather than bare cut-offs. Earlier runs fixed an uncorrected height and
# then capped the table at the strongest few peaks, which is not what papers do and which lets
# the effective cut float with the signal. :mod:`reporting` corrects for multiplicity, keeps
# every surviving local maximum 8 mm apart, and caps nothing.
SCHEMES = sys.argv[2:] or ["fdr", "fwe", "cluster"]
DATA = f"/tmp/claude-0/hcp/{CONTRAST}"
CACHE = f"/tmp/claude-0/hcp/{CONTRAST}_masked.npy"
DESIGNS = ((20, 16), (30, 12), (40, 9))   # (subjects per study, number of studies)
rng = np.random.default_rng(0)

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine


def load_all():
    if os.path.exists(CACHE):
        return np.load(CACHE)
    paths = sorted(glob.glob(f"{DATA}/*.nii.gz"))
    rows = []
    for i, path in enumerate(paths):
        try:
            img = nib.load(path)
            if img.ndim > 3:
                img = nib.Nifti1Image(np.asarray(img.dataobj)[..., 0], img.affine, img.header)
            rows.append(masker.transform(resample_to_img(
                img, mask_img, interpolation="continuous", force_resample=True,
                copy_header=True)).ravel().astype(np.float32))
        except Exception:
            continue
        if (i + 1) % 100 == 0:
            print(f"  loaded {i + 1}/{len(paths)}", flush=True)
    out = np.array(rows, dtype=np.float32)
    np.save(CACHE, out)
    return out


def hedges(mean, sd, n):
    """Hedges' g from a one-sample mean and sd over subjects."""
    correction = 1.0 - 3.0 / (4.0 * (n - 1) - 1.0)
    return correction * mean / np.maximum(sd, 1e-9)


ZOOMS = mask_img.header.get_zooms()[:3]


S = load_all()
print(f"\n{CONTRAST}: {S.shape[0]} subjects, {S.shape[1]} voxels at 4mm\n", flush=True)

for SCHEME in SCHEMES:
  for n_per_study, n_studies in DESIGNS:
      need = n_per_study * n_studies
      if need >= S.shape[0] - 100:
          print(f"{n_per_study}x{n_studies}: not enough subjects, skipped")
          continue
      order = rng.permutation(S.shape[0])
      used, held = order[:need], order[need:]
      truth_g = hedges(S[held].mean(axis=0), S[held].std(axis=0, ddof=1), len(held))

      studies, reported = [], 0
      for k in range(n_studies):
          block = S[used[k * n_per_study:(k + 1) * n_per_study]]
          mean, sd = block.mean(axis=0), block.std(axis=0, ddof=1)
          t = mean / np.maximum(sd / np.sqrt(n_per_study), 1e-9)
          z = np.nan_to_num(t_to_z(t, dof=n_per_study - 1))
          found, height = report_peaks(z, mask_bool, shape, ZOOMS, scheme=SCHEME)
          if len(found) < 2:
              continue
          reported += len(found)
          meta = {"sample_sizes": [n_per_study], "reporting_threshold": float(height)}
          studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
              {"id": f"s{k}", "name": "1", "metadata": meta,
               "points": [{"space": "MNI",
                           "coordinates": [float(c) for c in nib.affines.apply_affine(
                               affine, np.asarray(ijk, dtype=float))],
                           "values": [{"kind": "Z", "value": value}]} for ijk, value in found]}]})
      if len(studies) < 5:
          print(f"{n_per_study}x{n_studies}: only {len(studies)} studies reported, skipped")
          continue

      est = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
                 threshold="reporting_threshold")
      result = est.fit(Studyset({"id": "hcp", "name": "hcp", "studies": studies},
                                target=None, mask=mask_img))
      g = np.abs(result.get_map("g", return_type="array").ravel())
      pi = (result.get_map("prevalence", return_type="array").ravel()
            if "prevalence" in result.maps else np.ones_like(g))
      marg = pi * g
      covered = result.get_map("n_studies", return_type="array").ravel() > 0
      truth = np.abs(truth_g)
      used_cut = float(np.median(np.abs(est._cutoffs_z_.values)))

      print(f"--- {SCHEME}: {n_per_study} subjects x {n_studies} studies "
            f"({len(studies)} reported, {reported / max(len(studies),1):.0f} peaks each, "
            f"{len(held)} held out, height z = {used_cut:.2f}) ---")
      print(f"  {'truth stratum':>18} {'vox':>7} {'truth g':>8} {'CBES g':>8} {'ratio':>7} "
            f"{'pi*g':>8} {'ratio':>7} {'r':>7}")
      for lo, hi in ((0, 50), (50, 75), (75, 90), (90, 99), (99, 100)):
          band = (truth >= np.percentile(truth, lo)) & (
              truth < np.percentile(truth, hi) if hi < 100 else np.ones_like(truth, bool))
          use = covered & band & np.isfinite(g) & (g > 0)
          if use.sum() < 50:
              print(f"  {f'{lo}-{hi}%':>18} {int(use.sum()):>7}   too few")
              continue
          t = truth[use].mean()
          print(f"  {f'{lo}-{hi}%':>18} {int(use.sum()):>7} {t:8.3f} "
                f"{g[use].mean():8.3f} {g[use].mean()/max(t,1e-9):7.3f} "
                f"{marg[use].mean():8.3f} {marg[use].mean()/max(t,1e-9):7.3f} "
                f"{stats.pearsonr(g[use], truth[use])[0]:7.3f}")
      use = covered & np.isfinite(g) & (g > 0)
      t = truth[use].mean()
      print(f"  {'all covered':>18} {int(use.sum()):>7} {t:8.3f} "
            f"{g[use].mean():8.3f} {g[use].mean()/max(t,1e-9):7.3f} "
            f"{marg[use].mean():8.3f} {marg[use].mean()/max(t,1e-9):7.3f} "
            f"{stats.pearsonr(g[use], truth[use])[0]:7.3f}   "
            f"(pi*g correlates {stats.pearsonr(marg[use], truth[use])[0]:+.3f})\n", flush=True)
