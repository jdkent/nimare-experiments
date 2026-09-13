"""Split-half calibration of peak_bias, with no reference to the answer.

rho is estimated only from studies that supply BOTH an image and coordinates: fit those
studies from their images, fit the same studies from their coordinates, and take the ratio.
It is then applied to a disjoint set of studies and checked against their images. Nothing
in the calibration sees the held-out studies or the 21-study reference.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.cbma import CBES

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
cs = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
ids = np.array(list(cs.ids))

def restrict(studyset, keep_ids):
    """A studyset view holding only these studies."""
    return studyset.slice(analyses=[a for a in studyset.ids if a in set(keep_ids)])

def fit_map(studyset, use_images, peak_bias=None):
    est = CBES(fwhm=10.0, null_method="parametric", use_images=use_images,
               peak_bias=peak_bias)
    return est.fit(studyset).get_map("g", return_type="array").ravel()

def summarise(a, b):
    """Agreement between two maps over the voxels both cover."""
    ok = (a != 0) & (b != 0) & np.isfinite(a) & np.isfinite(b)
    return stats.pearsonr(a[ok], b[ok])[0], a[ok].mean(), b[ok].mean()

rng = np.random.default_rng(1)
rhos, agree_before, agree_after = [], [], []
for rep in range(4):
    perm = rng.permutation(len(ids))
    calib, test = ids[perm[:10]], ids[perm[10:]]

    cal_set = restrict(cs, calib)
    from_images = fit_map(cal_set, True)
    from_coords = fit_map(cal_set, False)
    ok = (from_images != 0) & (from_coords != 0)
    # The fit scales exactly linearly in rho, so the calibration target is the ratio of the
    # summaries a reader would compare -- not a regression slope, which is attenuated by the
    # many voxels where one study peaked and the images say nothing.
    rho = float(from_images[ok].mean() / from_coords[ok].mean())
    rhos.append(rho)

    test_set = restrict(cs, test)
    truth = fit_map(test_set, True)                       # held-out studies' own images
    raw = fit_map(test_set, False)                        # their coordinates, uncorrected
    fixed = fit_map(test_set, False, peak_bias=rho)       # corrected with calibration rho
    r0, m_raw, m_truth = summarise(raw, truth)
    r1, m_fixed, _ = summarise(fixed, truth)
    agree_before.append(m_raw / m_truth)
    agree_after.append(m_fixed / m_truth)
    print(f"rep {rep}: rho={rho:.3f} | held-out images {m_truth:.3f} | "
          f"coords raw {m_raw:.3f} ({m_raw/m_truth:.2f}x) | "
          f"corrected {m_fixed:.3f} ({m_fixed/m_truth:.2f}x) | r={r1:.3f}", flush=True)

print(f"\ncalibrated rho: {np.mean(rhos):.3f} +- {np.std(rhos):.3f}")
print(f"coordinate/image ratio before correction: {np.mean(agree_before):.2f}x")
print(f"coordinate/image ratio after  correction: {np.mean(agree_after):.2f}x   (1.00 = coherent)")
