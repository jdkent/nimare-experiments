"""Does bias(f) = b0 (1-f)/((1-f) + r f) transfer from 12 studies to 24?

Fitted on the twelve-study rows only, then evaluated at the twenty-four-study image fractions,
which the fit never saw. The all-image arm is the floor both channels sit on, so the model is
fitted to the *excess* over it rather than to the raw bias.
"""
import numpy as np
from scipy.optimize import least_squares

FLOOR = {12: -0.018, 24: -0.022}          # all-image arms, 100 reps each
NONE = {12: {2: +0.127, 6: +0.027, 12: -0.018}, 24: {2: +0.163, 6: +0.073}}
CAL = {12: {2: -0.038, 6: -0.026, 12: -0.018}, 24: {2: -0.064, 6: -0.049, 24: -0.022}}
COORDS_ONLY = {12: +0.255, 24: +0.246}


def model(f, b0, r):
    return b0 * (1 - f) / ((1 - f) + r * f)


def fit(rows, ns, b0):
    """One free parameter, r: b0 is the coordinates-only excess, measured not fitted."""
    f = np.array([ni / ns for ni in rows])
    y = np.array([rows[ni] - FLOOR[ns] for ni in rows])
    out = least_squares(lambda p: model(f, b0, p[0]) - y, [4.0], bounds=([0.1], [50.0]))
    return float(out.x[0])


for label, data in (("peak_bias=None", NONE), ("calibrated", CAL)):
    b0 = COORDS_ONLY[12] - FLOOR[12] if label == "peak_bias=None" else None
    if b0 is None:
        # No coordinates-only arm exists for the calibrated configuration (no donors to read a
        # scale from), so b0 is fitted too -- it is the corrected channel's own residual.
        f = np.array([ni / 12 for ni in data[12]])
        y = np.array([data[12][ni] - FLOOR[12] for ni in data[12]])
        o = least_squares(lambda p: model(f, p[0], p[1]) - y, [-0.04, 4.0],
                          bounds=([-1.0, 0.1], [1.0, 50.0]))
        b0, r = float(o.x[0]), float(o.x[1])
    else:
        r = fit(data[12], 12, b0)
    print(f"\n{label}: fitted on 12 studies -> b0 = {b0:+.4f}, r = {r:.2f}")
    print(f"  {'arm':>18s} {'f':>6s} {'measured':>9s} {'predicted':>10s} {'error':>7s}")
    for ns in (12, 24):
        for ni, meas in data.get(ns, {}).items():
            f = ni / ns
            pred = model(f, b0, r) + FLOOR[ns]
            tag = "(fitted)" if ns == 12 else "(HELD OUT)"
            print(f"  {ni:2d} of {ns:2d} {tag:>9s} {f:6.3f} {meas:+9.3f} {pred:+10.3f} "
                  f"{meas - pred:+7.3f}")
print("\nr is how much pooling weight one image study carries relative to one coordinate study.")
