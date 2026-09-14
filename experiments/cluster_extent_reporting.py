"""Items 6 and 9: what cluster-extent reporting does to threshold inference and to the fit.

`infer_threshold_from_minimum` inverts the order-statistic bias of the *smallest* reported peak
under the assumption that peaks were selected by height, one local maximum per cluster. Modern
papers overwhelmingly threshold by cluster extent or TFCE instead, which keeps a whole cluster
whose peak may sit only just above the height cut while discarding isolated peaks well above it.
So the modal real input violates the assumption, and rho_k depends on the inferred threshold,
which means the magnitude correction inherits whatever error that introduces.

This also breaks the circularity that makes the pain validation flattering: here the reporting
model that generates the peaks is deliberately *not* the one the estimator assumes.

Simulated directly from smooth fields so that "cluster" means something: a field is thresholded
at a height cut, connected components are found, and components are kept or dropped by extent.
Under height reporting every suprathreshold peak is reported; under extent reporting only peaks
in surviving clusters are, and the surviving set is biased toward broad, lower peaks.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from scipy import ndimage
from nimare.meta.cbma.effectsize import infer_threshold_from_minimum

RNG = np.random.default_rng(0)
SHAPE = (40, 40, 40)
TRUE_U = 3.2905


def smooth_field(fwhm_voxels, effect=0.0, centre=None):
    """A smooth Gaussian random field, optionally with a blob of true signal in it."""
    noise = RNG.normal(size=SHAPE)
    sigma = fwhm_voxels / np.sqrt(8 * np.log(2))
    field = ndimage.gaussian_filter(noise, sigma)
    field /= field.std()
    if effect and centre is not None:
        grid = np.indices(SHAPE).astype(float)
        blob = np.exp(-sum((grid[i] - centre[i]) ** 2 for i in range(3)) / (2 * 6.0**2))
        field = field + effect * blob
    return field


def peaks_from(field, u, min_extent):
    """Local maxima above ``u``, keeping only those in clusters of at least ``min_extent``."""
    above = field > u
    if not above.any():
        return np.array([])
    labels, n = ndimage.label(above)
    if not n:
        return np.array([])
    sizes = ndimage.sum(above, labels, index=np.arange(1, n + 1))
    keep = {i + 1 for i, s in enumerate(sizes) if s >= min_extent}
    # One local maximum per surviving cluster, which is what a paper reports.
    out = []
    for lab in keep:
        where = labels == lab
        out.append(float(field[where].max()))
    return np.array(out)


print(f"true height threshold z = {TRUE_U}; 200 fields per row\n")
print(f"{'reporting':>22s} {'min extent':>11s} {'peaks/field':>12s} {'min peak':>9s} "
      f"{'inferred u':>11s} {'error':>8s}")
for label, min_extent in (("height (as assumed)", 1),
                          ("extent >= 10 voxels", 10),
                          ("extent >= 50 voxels", 50),
                          ("extent >= 200 voxels", 200)):
    minima, inferred, counts = [], [], []
    for _ in range(200):
        peaks = peaks_from(smooth_field(6.0), TRUE_U, min_extent)
        if peaks.size < 2:
            continue
        counts.append(peaks.size)
        minima.append(peaks.min())
        inferred.append(infer_threshold_from_minimum(peaks.min(), peaks.size))
    if not inferred:
        print(f"{label:>22s} {min_extent:11d} {'no clusters survived':>12s}")
        continue
    got = float(np.mean(inferred))
    print(f"{label:>22s} {min_extent:11d} {np.mean(counts):12.1f} {np.mean(minima):9.3f} "
          f"{got:11.3f} {got - TRUE_U:+8.3f}", flush=True)

print()
print("The same, with real signal present, since a paper's clusters contain the effect:")
print(f"{'reporting':>22s} {'min extent':>11s} {'peaks/field':>12s} {'inferred u':>11s} "
      f"{'error':>8s}")
centre = tuple((np.array(SHAPE) - 1) / 2.0)
for label, min_extent in (("height (as assumed)", 1), ("extent >= 50 voxels", 50)):
    inferred, counts = [], []
    for _ in range(200):
        peaks = peaks_from(smooth_field(6.0, effect=3.0, centre=centre), TRUE_U, min_extent)
        if peaks.size < 2:
            continue
        counts.append(peaks.size)
        inferred.append(infer_threshold_from_minimum(peaks.min(), peaks.size))
    got = float(np.mean(inferred))
    print(f"{label:>22s} {min_extent:11d} {np.mean(counts):12.1f} {got:11.3f} "
          f"{got - TRUE_U:+8.3f}", flush=True)
