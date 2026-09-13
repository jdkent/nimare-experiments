"""Simulate studies as smooth random fields, so reported peaks are genuine local maxima.

The existing generator draws one value at the ground-truth location and thresholds it. That
models the reporting threshold but not the *selection of a location*, so it cannot produce --
and therefore cannot be used to test a correction for -- the peak-height bias that dominates
real data.
"""
import numpy as np
from scipy import ndimage

SHAPE = (30, 30, 30)
VOXEL_MM = 3.0

def true_field(peak_g=0.6, radius_vox=4.0):
    """A single Gaussian blob of true effect at the centre of the volume."""
    grid = np.indices(SHAPE).astype(float)
    centre = (np.array(SHAPE) - 1) / 2.0
    d2 = sum((grid[i] - centre[i]) ** 2 for i in range(3))
    return peak_g * np.exp(-d2 / (2 * radius_vox**2))

def study_image(truth, n_subjects, smooth_vox, rng):
    """One study's observed Hedges' g map: truth plus smooth noise of the right scale."""
    noise = rng.standard_normal(SHAPE)
    noise = ndimage.gaussian_filter(noise, smooth_vox)
    noise *= 1.0 / (np.std(noise) + 1e-12)      # restore unit variance after smoothing
    noise *= 1.0 / np.sqrt(n_subjects)          # sampling SD of g for a one-sample design
    return truth + noise

def report_peaks(image, n_subjects, threshold_z):
    """Local maxima clearing the study's reporting threshold, as a paper would list them."""
    cutoff = threshold_z / np.sqrt(n_subjects)
    local_max = ndimage.maximum_filter(image, size=3) == image
    hits = np.argwhere(local_max & (image > cutoff))
    return hits, image[tuple(hits.T)] if len(hits) else np.array([])

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    truth = true_field()
    for smooth in (1.0, 2.0, 3.0):
        at_peak, at_truth, counts = [], [], []
        for _ in range(40):
            n = int(rng.integers(15, 40))
            img = study_image(truth, n, smooth, rng)
            hits, values = report_peaks(img, n, 3.2905)
            if not len(hits):
                continue
            counts.append(len(hits))
            at_peak.extend(np.abs(values))
            at_truth.extend(truth[tuple(hits.T)])
        print(f"smoothing {smooth:.1f} vox ({smooth*VOXEL_MM:.0f} mm): "
              f"{np.mean(counts):5.1f} peaks/study, reported |g| {np.mean(at_peak):.3f}, "
              f"true g there {np.mean(at_truth):.3f}, ratio {np.mean(at_truth)/np.mean(at_peak):.3f}")
