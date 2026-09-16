"""Does the implementation behave the way the algebra says it must?

`proofs/marginal_control_variate.py` establishes four things about the control-variate estimator.
This bed checks that the code in `nimare/meta/cbma/marginal.py` exhibits all four against a known
truth, and that the two documented failure modes actually fail. The algebra is settled; what is
being tested here is the implementation and the finite-sample behaviour the algebra is silent on.

**Kill conditions, stated before running.**

  1. If the estimator is biased with a correct predictor at any of the three regimes, the
     implementation does not match the proof and nothing else here matters.
  2. If the reported standard error is out by more than 15% of the sampling spread, the variance
     formula is not being computed as derived.
  3. If the achieved variance ratio beats the proven floor 1 - rho^2 anywhere, either the floor
     is wrong or the ratio is being computed against the wrong baseline.
  4. If a non-exchangeable coordinate cohort does *not* bias the estimate, the failure mode has
     been coded away rather than left visible, which would be worse than the bias.

Regimes are the three the design document asks for: 1 image with 20 tables, 1 with 500, and 8
with 100. The one-image rows are expected to report no precision at all -- there is no covariance
to read from a single study -- and that is a result rather than a gap.

Environment knobs: REPS, VOXELS, RHO, TRUTH, TAU.
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/home/user/NiMARE")
from nimare.meta.cbma.marginal import (  # noqa: E402
    achievable_ratio,
    control_variate_mean,
    optimal_coefficient,
)

REPS = int(os.environ.get("REPS", "2000"))
VOXELS = int(os.environ.get("VOXELS", "1"))
RHO = float(os.environ.get("RHO", "0.6"))
TRUTH = float(os.environ.get("TRUTH", "0.4"))
TAU = float(os.environ.get("TAU", "1.0"))

REGIMES = ((1, 20), (1, 500), (8, 100))


def draw(rng, n_images, n_coordinates, rho, shift=0.0):
    """Per-study effects and a predictor correlating ``rho`` with them.

    The predictor is standardised, and the effect is ``rho`` of it plus independent noise, so
    Corr(Y, f) = rho exactly by construction rather than approximately. ``shift`` moves the
    coordinate cohort's predictor only, which is the non-exchangeability failure.
    """
    total = n_images + n_coordinates
    common = rng.normal(size=(total, VOXELS))
    independent = rng.normal(size=(n_images, VOXELS))
    predictions = common.copy()
    predictions[n_images:] += shift
    effects = TRUTH + TAU * (rho * common[:n_images] + np.sqrt(1 - rho**2) * independent)
    return effects, predictions[:n_images], predictions[n_images:]


def run(rng, n_images, n_coordinates, rho, shift=0.0, coefficient="oracle"):
    """One regime. Returns per-replication estimates, reported standard errors and ratios."""
    estimates, reported, ratios, images_only = [], [], [], []
    # The oracle coefficient from the proof, using the covariance the simulator actually used
    # rather than one estimated from the same data -- that is what makes it an oracle.
    oracle = rho * TAU / ((1.0 + n_images / n_coordinates) * 1.0)
    for _ in range(REPS):
        effects, f_image, f_coord = draw(rng, n_images, n_coordinates, rho, shift=shift)
        lam = oracle if coefficient == "oracle" else optimal_coefficient(
            effects, f_image, n_images, n_coordinates, mode="pooled"
        )
        out = control_variate_mean(effects, f_image, f_coord, lam)
        estimates.append(out["estimate"])
        reported.append(out["se"])
        ratios.append(out["variance_ratio"])
        images_only.append(effects.mean(axis=0))
    return (
        np.concatenate(estimates),
        np.concatenate(reported),
        np.concatenate(ratios),
        np.concatenate(images_only),
    )


def main():
    print(f"true marginal mean {TRUTH}, per-study sd {TAU}, predictor correlation {RHO}")
    print(f"{REPS} replications per regime, {VOXELS} voxel(s) each")
    print()
    print(f"{'regime':>12} {'bias':>9} {'sd':>8} {'mean se':>9} {'se/sd':>7} "
          f"{'emp ratio':>10} {'proven':>8} {'floor':>7}")

    failures = []
    for n_images, n_coordinates in REGIMES:
        rng = np.random.default_rng(0)
        est, se, _ratio, img = run(rng, n_images, n_coordinates, RHO)
        bias = est.mean() - TRUTH
        sd = est.std(ddof=1)
        empirical_ratio = est.var(ddof=1) / img.var(ddof=1)
        proven = float(achievable_ratio(RHO, n_images, n_coordinates))
        floor = 1.0 - RHO**2
        finite = np.isfinite(se)
        mean_se = se[finite].mean() if finite.any() else float("nan")
        shown_se = f"{mean_se:9.4f}" if finite.any() else f"{'none':>9}"
        shown_ratio = f"{mean_se / sd:7.2f}" if finite.any() else f"{'--':>7}"
        print(f"{f'{n_images}+{n_coordinates}':>12} {bias:+9.4f} {sd:8.4f} {shown_se} "
              f"{shown_ratio} {empirical_ratio:10.3f} {proven:8.3f} {floor:7.3f}")

        # Kill condition 1: unbiased, to within three standard errors of the replication mean.
        if abs(bias) > 3.0 * sd / np.sqrt(REPS):
            failures.append(f"biased at {n_images}+{n_coordinates}: {bias:+.4f}")
        # Kill condition 2: the reported se describes the spread, where it is reported at all.
        if finite.any() and abs(mean_se / sd - 1.0) > 0.15:
            failures.append(f"se off at {n_images}+{n_coordinates}: {mean_se / sd:.2f}")
        # Kill condition 3: the floor holds.
        if empirical_ratio < floor - 3.0 * np.sqrt(2.0 / REPS):
            failures.append(
                f"ratio beat the floor at {n_images}+{n_coordinates}: "
                f"{empirical_ratio:.3f} < {floor:.3f}"
            )

    print()
    print("A '--' standard error is the one-image case: a single study carries no covariance, so")
    print("the estimate is produced and its precision is reported as unavailable.")

    print()
    print("non-exchangeable cohorts, shift of 1.0 in the coordinate predictor:")
    print(f"{'regime':>12} {'bias':>9} {'expected':>10}")
    for n_images, n_coordinates in REGIMES:
        rng = np.random.default_rng(1)
        oracle = RHO * TAU / (1.0 + n_images / n_coordinates)
        est, _se, _r, _i = run(rng, n_images, n_coordinates, RHO, shift=1.0)
        bias = est.mean() - TRUTH
        print(f"{f'{n_images}+{n_coordinates}':>12} {bias:+9.4f} {oracle:+10.4f}")
        # Kill condition 4: the failure mode must be present, and equal to lambda times the shift.
        if abs(bias - oracle) > max(0.05, 4.0 * est.std(ddof=1) / np.sqrt(REPS)):
            failures.append(
                f"non-exchangeability did not bias as predicted at {n_images}+{n_coordinates}"
            )

    print()
    print("estimated rather than oracle coefficient, 8+100 only (1 image cannot estimate one):")
    rng = np.random.default_rng(2)
    est, se, _r, img = run(rng, 8, 100, RHO, coefficient="estimated")
    print(f"  bias {est.mean() - TRUTH:+.4f}, sd {est.std(ddof=1):.4f}, "
          f"empirical ratio {est.var(ddof=1) / img.var(ddof=1):.3f} against "
          f"proven {float(achievable_ratio(RHO, 8, 100)):.3f}")
    print("  Estimating the coefficient on the same eight studies costs some of the reduction,")
    print("  which is the small-sample price the proof is silent about.")

    print()
    if failures:
        print("KILL CONDITIONS MET:")
        for line in failures:
            print(f"  {line}")
        return 1
    print("All four kill conditions passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
