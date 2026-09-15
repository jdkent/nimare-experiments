"""Test the central claim: does peak *count* beat peak *height* at predicting the truth?

The first-principles note argues that the information about the effect field is in the intensity
of the reported point process, not in its marks. The evidence so far is a simulation sweep
(`peak_height_curve`) plus a per-focus regression (`information_ceiling`). Neither compares the
two signals head to head on real data against an independent reference, which is what would
actually settle it.

For each voxel, three quantities are built from exactly the same coordinate tables:

  * **height** -- the kernel-weighted mean of the reported effect sizes reaching the voxel. This
    is the quantity CBES's `g` is built from.
  * **study count** -- how many distinct studies placed a focus within a radius. This is the
    intensity, and it is what `prevalence` approximates.
  * **foci count** -- how many foci in total, so the study-level and focus-level versions of the
    intensity can be told apart.

Each is scored against a truth measured on subjects that produced no coordinate. If the counts
win, the reframing holds and the marks really are close to uninformative. If the height wins, the
first-principles note is wrong and the simulation misled.

Scored two ways, because they answer different questions: correlation over all covered voxels,
and correlation restricted to the voxels where the truth is in its top decile, which asks whether
the signal is only "is there an effect here at all".
"""
import glob, logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from scipy.ndimage import gaussian_filter
from nilearn.datasets import load_mni152_brain_mask
from nilearn.maskers import NiftiMasker
from nimare.transforms import t_to_z, d_to_g, t_to_d
from reporting import report_peaks

CONTRAST = "MOTOR_LH"
CACHE = f"/tmp/claude-0/hcp/{CONTRAST}_masked.npy"
DESIGNS = ((20, 16), (30, 12))
SCHEMES = ("cluster", "fdr")
RADII_MM = (10.0, 20.0)
KERNEL_FWHM = 10.0

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape = mask_img.shape
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)
sigma_vox = (KERNEL_FWHM / 2.3548) / zooms

S = np.load(CACHE)
print(f"{CONTRAST}: {S.shape[0]} subjects, {S.shape[1]} voxels at 4mm\n", flush=True)


def hedges(mean, sd, n):
    return (1.0 - 3.0 / (4.0 * (n - 1) - 1.0)) * mean / np.maximum(sd, 1e-9)


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def ball(radius):
    r = int(np.ceil(radius / zooms.min()))
    grid = np.stack(np.meshgrid(*[np.arange(-r, r + 1)] * 3, indexing="ij"), -1)
    return grid[np.linalg.norm(grid * zooms, axis=-1) <= radius]


rng = np.random.default_rng(0)
print(f"  {'design':>10} {'scheme':>8} {'signal':>22} {'r, covered':>11} {'r, top decile':>14}")
for n_per, n_studies in DESIGNS:
    need = n_per * n_studies
    order = rng.permutation(S.shape[0])
    used, held = order[:need], order[need:]
    truth = np.abs(hedges(S[held].mean(0), S[held].std(0, ddof=1), len(held)))

    for scheme in SCHEMES:
        height_num = np.zeros(shape)
        height_den = np.zeros(shape)
        study_hit = np.zeros(shape, dtype=np.int32)
        foci_hit = {r: np.zeros(shape, dtype=np.int32) for r in RADII_MM}
        study_cnt = {r: np.zeros(shape, dtype=np.int32) for r in RADII_MM}
        offsets = {r: ball(r) for r in RADII_MM}
        reported = 0

        for k in range(n_studies):
            block = S[used[k * n_per:(k + 1) * n_per]]
            t = block.mean(0) / np.maximum(block.std(0, ddof=1) / np.sqrt(n_per), 1e-9)
            z = np.nan_to_num(t_to_z(t, dof=n_per - 1))
            foci, _ = report_peaks(z, mask_bool, shape, zooms, scheme=scheme, focus="max")
            if not foci:
                continue
            reported += 1
            # Height: each focus's g smeared by the pooling kernel, inverse-variance weighted
            # exactly as the estimator does it.
            spike = np.zeros(shape)
            for ijk, value in foci:
                spike[tuple(ijk)] += abs(float(to_g(np.array([value]), n_per)[0]))
            smeared = gaussian_filter(spike, sigma_vox, mode="constant")
            reach = gaussian_filter((spike > 0).astype(float), sigma_vox, mode="constant")
            height_num += smeared
            height_den += reach
            # Counts, at each radius.
            for r in RADII_MM:
                seen = np.zeros(shape, dtype=bool)
                for ijk, _ in foci:
                    pts = np.asarray(ijk) + offsets[r]
                    ok = np.all((pts >= 0) & (pts < np.array(shape)), axis=1)
                    pts = pts[ok]
                    idx = (pts[:, 0], pts[:, 1], pts[:, 2])
                    foci_hit[r][idx] += 1
                    seen[idx] = True
                study_cnt[r] += seen

        height = np.divide(height_num, height_den,
                           out=np.zeros(shape), where=height_den > 1e-9)
        covered = (height_den > 1e-9) & mask_bool
        t_flat = np.zeros(shape); t_flat[mask_bool] = truth
        use = covered & np.isfinite(height)
        top = use & (t_flat >= np.percentile(t_flat[use], 90))

        signals = [("height (what g uses)", height)]
        for r in RADII_MM:
            signals.append((f"study count r={r:.0f}mm", study_cnt[r].astype(float)))
            signals.append((f"foci count  r={r:.0f}mm", foci_hit[r].astype(float)))
        for name, sig in signals:
            r_all = stats.pearsonr(sig[use], t_flat[use])[0]
            r_top = stats.pearsonr(sig[top], t_flat[top])[0] if top.sum() > 50 else np.nan
            print(f"  {f'{n_per}x{n_studies}':>10} {scheme:>8} {name:>22} "
                  f"{r_all:+11.3f} {r_top:+14.3f}", flush=True)
        print(f"  {'':>10} {'':>8} {f'({reported} studies reported)':>22}\n", flush=True)
