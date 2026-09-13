"""How well does the smallest reported statistic recover a known threshold?

Also tests an order-statistic correction: the minimum of m peaks drawn from the RFT null
peak-height distribution above u sits above u by a predictable amount, so it can be removed.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from scipy.optimize import brentq
from load_pain import load_pain
from nimare.transforms import ImageTransformer, ImagesToCoordinates

SQRT3 = np.sqrt(3.0)

def survival(z, u):
    """P(peak height > z | peak above u) for a smooth 3D Gaussian field."""
    z = np.asarray(z, float)
    return np.clip((z**2 - 1) * np.exp(-0.5 * z**2) / ((u**2 - 1) * np.exp(-0.5 * u**2)), 0, 1)

def expected_min(u, m, span=8.0, n=600):
    """E[min of m peaks above u] = u + integral of S(z|u)^m."""
    grid = np.linspace(u, u + span, n)
    return u + np.trapezoid(survival(grid, u) ** m, grid)

def infer_threshold(z_min, m):
    """Invert expected_min: the u whose m-peak minimum is z_min."""
    if m <= 0 or not np.isfinite(z_min):
        return np.nan
    lo, hi = 0.5, float(z_min)
    try:
        return brentq(lambda u: expected_min(u, m) - z_min, lo, hi, xtol=1e-4)
    except ValueError:
        return float(z_min)

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
print(f"{'imposed z':>10s} {'studies':>8s} {'peaks/study':>12s} "
      f"{'study-min (raw)':>16s} {'order-stat corrected':>21s}")
for true_u in (2.3263, 3.0902, 3.2905, 4.2649):
    cs = ImagesToCoordinates(merge_strategy="demolish", z_threshold=true_u, two_sided=True,
                             remove_subpeaks=True).transform(ss)
    df = cs.coordinates
    if not len(df):
        print(f"{true_u:10.3f}   (no study reported anything)")
        continue
    raw, fixed, counts = [], [], []
    for sid, sub in df.groupby("id"):
        z = np.abs(sub["z_stat"].astype(float).to_numpy())
        raw.append(z.min()); counts.append(len(z))
        fixed.append(infer_threshold(z.min(), len(z)))
    print(f"{true_u:10.3f} {len(raw):8d} {np.mean(counts):12.1f} "
          f"{np.mean(raw):16.3f} {np.nanmean(fixed):21.3f}")
