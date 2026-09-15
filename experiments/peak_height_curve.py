"""How much does a reported peak's height actually say about the effect underneath it?

``peak_selection_event`` showed the truncated normal is the wrong selection model for a local
maximum: at a true mean of 0.5 its MLE returns 0.26, and at 2.0 it returns 3.6. So the question
is what the right one is, and whether it carries any usable information. Both are answered by
one curve -- the mean of ``|peak|`` against the mean of the field it sits in, at each threshold.

Read two things off it:

  * the *floor*: what ``|peak|`` is when the field is pure noise. That is the threshold plus the
    overshoot, and it is the quantity that gets mistaken for an effect.
  * the *slope*: how fast the mean peak rises with the true mean. A slope near zero says peak
    heights are uninformative at that threshold and no correction recovers the effect; a slope
    near one says the height is the effect plus a constant, and subtracting the floor is enough.

The field is smoothed white noise at unit variance, so the mean of the field *is* Hedges' g on
the subject scale times sqrt(n) -- the same z scale the reporting threshold lives on.
"""
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter

rng = np.random.default_rng(1)
SHAPE = (48, 48, 48)
N_FIELDS = 200
MUS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)
CUTS = (2.5, 3.29, 4.0, 5.0)


def peaks_abs(volume, cut):
    mag = np.abs(volume)
    hit = (mag == maximum_filter(mag, size=3)) & (mag >= cut)
    return mag[hit]


for smooth in (1.5, 2.5):
    print(f"\nsmoothing sigma {smooth} voxels")
    header = "".join(f"{f'cut {c}':>22}" for c in CUTS)
    print(f"{'true mean':>10}{header}")
    print(f"{'':>10}" + "".join(f"{'mean|peak|':>12}{'per field':>10}" for _ in CUTS))
    rows = {}
    for mu in MUS:
        cells = []
        for cut in CUTS:
            found, counts = [], []
            for _ in range(N_FIELDS):
                raw = gaussian_filter(rng.standard_normal(SHAPE), smooth)
                raw /= raw.std()
                got = peaks_abs(mu + raw, cut)
                counts.append(got.size)
                if got.size:
                    found.append(got)
            n = sum(counts)
            value = np.concatenate(found).mean() if found else np.nan
            cells.append((value, n / N_FIELDS))
        rows[mu] = cells
        print(f"{mu:10.2f}" + "".join(f"{v:12.3f}{c:10.1f}" for v, c in cells))

    print(f"\n  {'cut':>5} {'floor at mu=0':>14} {'slope 0->1':>11} {'slope 1->3':>11}")
    for j, cut in enumerate(CUTS):
        floor = rows[0.0][j][0]
        s01 = rows[1.0][j][0] - rows[0.0][j][0]
        s13 = (rows[3.0][j][0] - rows[1.0][j][0]) / 2.0
        print(f"  {cut:5.2f} {floor:14.3f} {s01:11.3f} {s13:11.3f}")
