"""Calibrate the coverage testbed before trusting anything measured in it.

The first coverage run was set up by guess and turned out to sit in a regime where only 18% of
studies reported a single focus, so a nominally 12-study meta-analysis was really a 2-study one
and 30% of replications could not be fitted at all -- with the survivors selected for high
signal, which biases the very number the run was measuring.

Realistic is the regime coordinate meta-analysis actually draws from: a published contrast
almost always reports something (that is partly why it is published) and typically lists
several foci, not one. So the bed is tuned until a typical study reports a handful of clusters
and nearly every study reports at least one, and only then used.

Sweeps the true peak, the smoothing, and the number of true blobs, and prints the measured
field smoothness and the extent threshold random field theory demands alongside -- because the
extent demand is what silently dominated the first attempt.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import ndimage
import reporting

SHAPE, VOXEL_MM = (30, 30, 30), 4.0
ZOOMS = np.full(3, VOXEL_MM)
MASK_BOOL = np.ones(SHAPE, dtype=bool)
N_STUDIES, RADIUS_VOX = 30, 2.5

#: Several true sites of differing strength, as a real contrast has.
SITE_IJK = [(10, 15, 15), (20, 15, 15), (15, 9, 18), (15, 21, 12), (15, 15, 22)]
SITE_WEIGHT = [1.0, 1.0, 0.75, 0.75, 0.5]


def truth_field(peak, n_sites):
    out = np.zeros(SHAPE)
    grid = np.indices(SHAPE).astype(float)
    for (i, j, k), w in list(zip(SITE_IJK, SITE_WEIGHT))[:n_sites]:
        d2 = (grid[0] - i) ** 2 + (grid[1] - j) ** 2 + (grid[2] - k) ** 2
        out = np.maximum(out, peak * w * np.exp(-d2 / (2 * RADIUS_VOX**2)))
    return out


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    print(f"{SHAPE} at {VOXEL_MM:.0f} mm, radius {RADIUS_VOX} vox, n in [20,40]\n")
    print(f"{'sigma':>5s} {'FWHM mm':>8s} {'min ext':>8s} {'peak g':>6s} {'sites':>5s} "
          f"{'peak z@30':>9s} {'% reporting':>11s} {'foci/study':>10s} {'foci if any':>11s}")
    for sigma in (0.8, 1.0, 1.5, 2.0):
        probe = ndimage.gaussian_filter(rng.standard_normal(SHAPE), sigma)
        probe *= 1.0 / (probe.std() + 1e-12)
        fwhm = reporting.smoothness_fwhm(probe, MASK_BOOL, ZOOMS)
        ext = reporting.cluster_extent_threshold(probe, MASK_BOOL, ZOOMS)
        for peak in (0.5, 0.8, 1.2):
            for n_sites in (5,):
                truth = truth_field(peak, n_sites)
                counts = []
                for _ in range(N_STUDIES):
                    n = int(rng.integers(20, 41))
                    nz = ndimage.gaussian_filter(rng.standard_normal(SHAPE), sigma)
                    nz *= 1.0 / (nz.std() + 1e-12)
                    z = (truth + nz / np.sqrt(n)) * np.sqrt(n)
                    foci, _ = reporting.report_peaks(z[MASK_BOOL], MASK_BOOL, SHAPE, ZOOMS,
                                                     "cluster", "max")
                    counts.append(len(foci))
                c = np.array(counts)
                nonzero = c[c > 0]
                print(f"{sigma:5.1f} {fwhm.mean():8.1f} {ext:8.1f} {peak:6.2f} {n_sites:5d} "
                      f"{truth.max()*np.sqrt(30):9.2f} {np.mean(c > 0):11.2f} "
                      f"{c.mean():10.2f} "
                      f"{(nonzero.mean() if nonzero.size else 0):11.2f}", flush=True)
    print("\nTarget: '% reporting' above ~0.9 and 'foci if any' in the 3-10 range that real")
    print("coordinate tables show. Anything far below that is a different experiment.")
