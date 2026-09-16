"""Does the scalar censored reference reproduce the design document's own section 8.1 table?

The document reports a 4,000-replication scalar experiment computed independently of anything
here. That makes it the calibration target for ``nimare/meta/cbma/censored.py``: not a bed of my
own construction, which could agree with my implementation because both share my mistakes.

Its setup, quoted: true mean .4; within-study SD .2; between-study SD .15; total SD .25; known
variance; directional threshold .6; all studies remain on the roster; half of threshold
exceedances are independently retained. Intervals are likelihood-ratio intervals with known
variance, and the image-only baseline is likewise given the true heterogeneity, "deliberately
[removing] heterogeneity estimation as a confound". Monte Carlo SE near 95% coverage is 0.34
percentage points.

    regime    image-only   correct retention   ignoring retention
    1 + 20    .256         .129  (93.95%)      .157  (77.60%)
    1 + 500   .251         .023  (94.93%)      .114  ( 0.00%)
    8 + 100   .090         .047  (94.75%)      .103  (23.45%)

All three columns are checked. The module's plain interval likelihood *is* the "ignoring
retention" model -- a reported study contributes Y > .6 and an unreported one Y < .6, with no
retention probability anywhere -- and its retention-aware path, derived in
``proofs/retention_in_the_reporting_model.py`` and shown there to be bit-identical to the plain
one at rho = 1, is the "correct retention" model with rho supplied at its true .5.

The algebra also predicts the ignoring-retention column in closed form before any of this runs:
the misspecified population root satisfies S(m_hat) = rho * S(m_0), giving an asymptotic bias of
-.1121 at these settings against the document's measured .114 at 1+500. So that column has two
independent predictions to match, not one.

**Stated kill condition, before running.** The image-only and ignoring-retention RMSEs must
agree with the table to within about three Monte Carlo standard errors, and the
ignoring-retention coverages must reproduce the collapse (77.6%, 0%, 23.45%) rather than merely
being "somewhat low". A near miss on the coverages is not a pass: their whole point is that a
misspecified reporting model fails catastrophically and visibly, and an implementation that made
them merely mediocre would be a different model from the document's.

Environment knobs: REPS, REGIMES.
"""
import os
import sys

import numpy as np
from scipy.stats import chi2, norm

sys.path.insert(0, "/home/user/NiMARE")
from nimare.meta.cbma.censored import (  # noqa: E402
    ObservationState,
    bounds_from_states,
    censored_loglik,
    retention_roles,
)

TRUE_MEAN = 0.4
WITHIN_SD = 0.2
BETWEEN_SD = 0.15
THRESHOLD = 0.6
RETENTION = 0.5

WITHIN_VAR = WITHIN_SD**2
BETWEEN_VAR = BETWEEN_SD**2
TOTAL_SD = np.sqrt(WITHIN_VAR + BETWEEN_VAR)
CUTOFF = chi2.ppf(0.95, 1) / 2.0

#: Grid for the likelihood-ratio interval, wide enough that a bound running off it is a
#: recorded event rather than a silent truncation. The document records one such case.
GRID = np.linspace(-1.5, 2.0, 1401)


def draw(n_images, n_coordinates, rng):
    """One replication: image estimates, and the reporting outcome of each coordinate study."""
    images = rng.normal(TRUE_MEAN, TOTAL_SD, size=n_images)
    coordinate_values = rng.normal(TRUE_MEAN, TOTAL_SD, size=n_coordinates)
    exceeded = coordinate_values > THRESHOLD
    reported = exceeded & (rng.random(n_coordinates) < RETENTION)
    return images, reported


def image_only(images):
    """The baseline: the mean of the images, with the true total variance."""
    mean = float(images.mean())
    se = TOTAL_SD / np.sqrt(images.size)
    half = norm.ppf(0.975) * se
    return mean, mean - half, mean + half


def _bounds(images, reported):
    """Interval bounds and retention roles for one replication.

    A reported study is read as Y > c and an unreported one as Y < c, which is what the
    implementation's ``DIRECTION_ONLY`` and one-sided ``ABSENT_COMPLETE_TABLE`` states mean. The
    roles are what the retention-aware arm needs; the ignoring-retention arm never looks at them.
    """
    count = images.size + reported.size
    states = [ObservationState.IMAGE] * images.size + [
        ObservationState.DIRECTION_ONLY if flag else ObservationState.ABSENT_COMPLETE_TABLE
        for flag in reported
    ]
    values = np.concatenate([images, np.full(reported.size, np.nan)])
    thresholds = np.concatenate([np.full(images.size, np.nan), np.full(reported.size, THRESHOLD)])
    signs = np.concatenate([np.full(images.size, np.nan), np.ones(reported.size)])
    lower, upper = bounds_from_states(
        states, values=values, thresholds=thresholds, signs=signs
    )
    return lower, upper, np.full(count, WITHIN_VAR), retention_roles(states)


def profile_on_grid(lower, upper, variances, retention=None, roles=None):
    """Likelihood-ratio interval with the heterogeneity held at its true value.

    The grid is used rather than the module's own search so that the comparison isolates the
    likelihood from the optimiser: the document's intervals are likelihood-ratio intervals with
    known variance, and this computes exactly that.
    """
    profile = np.array(
        [
            censored_loglik(
                mean, BETWEEN_VAR, lower, upper, variances, retention=retention, roles=roles
            )
            for mean in GRID
        ]
    )
    best = int(np.argmax(profile))
    inside = profile >= profile[best] - CUTOFF
    indices = np.flatnonzero(inside)
    touched = bool(indices[0] == 0 or indices[-1] == GRID.size - 1)
    return float(GRID[best]), float(GRID[indices[0]]), float(GRID[indices[-1]]), touched


def run(n_images, n_coordinates, reps, seed):
    rng = np.random.default_rng(seed)
    arms = ("images", "ignoring", "correct")
    rows = {name: {"error": [], "covered": [], "touched": 0} for name in arms}
    for _ in range(reps):
        images, reported = draw(n_images, n_coordinates, rng)

        mean, low, high = image_only(images)
        rows["images"]["error"].append(mean - TRUE_MEAN)
        rows["images"]["covered"].append(low <= TRUE_MEAN <= high)

        lower, upper, variances, roles = _bounds(images, reported)

        mean, low, high, touched = profile_on_grid(lower, upper, variances)
        rows["ignoring"]["error"].append(mean - TRUE_MEAN)
        rows["ignoring"]["covered"].append(low <= TRUE_MEAN <= high)
        rows["ignoring"]["touched"] += int(touched)

        mean, low, high, touched = profile_on_grid(
            lower, upper, variances, retention=RETENTION, roles=roles
        )
        rows["correct"]["error"].append(mean - TRUE_MEAN)
        rows["correct"]["covered"].append(low <= TRUE_MEAN <= high)
        rows["correct"]["touched"] += int(touched)

    out = {}
    for name, row in rows.items():
        errors = np.asarray(row["error"])
        covered = np.asarray(row["covered"], dtype=float)
        out[name] = {
            "rmse": float(np.sqrt(np.mean(errors**2))),
            "rmse_se": float(np.std(errors**2) / (2 * np.sqrt(np.mean(errors**2)) * np.sqrt(reps))),
            "coverage": float(covered.mean()),
            "coverage_se": float(np.sqrt(covered.mean() * (1 - covered.mean()) / reps)),
            "touched": row["touched"],
        }
    return out


TARGETS = {
    (1, 20): {
        "images": (0.256, None),
        "ignoring": (0.157, 0.7760),
        "correct": (0.129, 0.9395),
    },
    (1, 500): {
        "images": (0.251, None),
        "ignoring": (0.114, 0.0000),
        "correct": (0.023, 0.9493),
    },
    (8, 100): {
        "images": (0.090, None),
        "ignoring": (0.103, 0.2345),
        "correct": (0.047, 0.9475),
    },
}


if __name__ == "__main__":
    reps = int(os.environ.get("REPS", 4000))
    wanted = os.environ.get("REGIMES")
    regimes = [(1, 20), (1, 500), (8, 100)]
    if wanted:
        regimes = [regimes[int(index)] for index in wanted.split(",")]

    print(f"Calibration against the design document's section 8.1, {reps} replications each.")
    print(f"{'regime':>10} {'arm':>10} {'rmse':>16} {'target':>8} "
          f"{'coverage':>18} {'target':>8}  verdict")
    verdicts = []
    for index, (n_images, n_coordinates) in enumerate(regimes):
        result = run(n_images, n_coordinates, reps, seed=1000 + index)
        for arm in ("images", "ignoring", "correct"):
            row = result[arm]
            rmse_target, coverage_target = TARGETS[(n_images, n_coordinates)][arm]
            rmse_gap = abs(row["rmse"] - rmse_target) / max(row["rmse_se"], 1e-9)
            verdict = "ok" if rmse_gap < 3 else f"RMSE OFF BY {rmse_gap:.1f} se"
            coverage_shown = "--"
            if coverage_target is not None:
                # Compare against the standard error the *target* implies, not the one the
                # sample happens to show. An all-or-nothing sample gives a binomial se of
                # exactly zero, which is not a standard error, and dividing by it turns a
                # one-se miss into an infinite one. The target proportion is clipped away
                # from 0 and 1 by half an observation for the same reason.
                clipped = min(max(coverage_target, 0.5 / reps), 1 - 0.5 / reps)
                null_se = np.sqrt(clipped * (1 - clipped) / reps)
                coverage_gap = abs(row["coverage"] - coverage_target) / null_se
                coverage_shown = f"{coverage_target:.4f}"
                if coverage_gap >= 3 and verdict == "ok":
                    verdict = f"COVERAGE OFF BY {coverage_gap:.1f} se"
            verdicts.append(verdict)
            print(
                f"{n_images:4d}+{n_coordinates:<5d} {arm:>10} "
                f"{row['rmse']:.4f} +- {row['rmse_se']:.4f} {rmse_target:>8.3f} "
                f"{row['coverage']:.4f} +- {row['coverage_se']:.4f} {coverage_shown:>8}  "
                f"{verdict}"
                + (f"  [{row['touched']} touched the grid]" if row["touched"] else "")
            )
    print()
    if all(verdict == "ok" for verdict in verdicts):
        print("PASS: the implementation reproduces all three columns of the independently")
        print("computed table, so the interval likelihood, its one-sided absence state and its")
        print("retention-aware path are the models the document describes.")
    else:
        print("FAIL: at least one arm disagrees with the document. Per the stated kill")
        print("condition this stops everything built on this module until it is resolved --")
        print("the defect is in the implementation, in my reading of the setup, or in the")
        print("table, and which one is not something to guess at.")
