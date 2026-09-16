"""Does the block model reproduce the design document's section 8.2 table?

Its second experiment, quoted: 1,000 replications; each study has a nine-element region
``Y_ij = m + U_i + eps_ij`` with ``U ~ N(0, .15^2)``, ``eps ~ N(0, .2^2)``, ``m = .4``; the
reporting rule retains the regional maximum if it exceeds .75; no-report studies remain
included; known covariance.

    regime    image-only   correct peak likelihood   correct indicator-only
    1 + 20    .165         .050                      .053
    1 + 500   .167         .010                      .011
    8 + 100   .057         .022                      .024

**The one-maximum cap is this measurement's explicit subject**, as the document says of its own
experiment -- "this experiment intentionally studies a one-maximum cap; it is not claimed to
reproduce whole-brain peak extraction". Every study here gets exactly one block, so the cap is
the design rather than a convenience, and nothing measured here is evidence about a corpus whose
studies report however many peaks survive.

The naive "treating peaks as ordinary values" column is not attempted. Four candidate
constructions give asymptotic biases of +0.060, +0.130, +0.490 and +0.000 against the +0.256 its
RMSE implies, so its construction is not pinned down by the description and guessing at it would
produce a number that means nothing.

**Stated kill condition, before running.** Both modelled columns must agree with the table to
within three standard errors of the difference between two 1,000-replication estimates quoted to
three decimals, and the ratio between them -- the thing the document actually concludes from,
that heights add little beyond the indicator -- must land in 0.90 to 0.95, the range its own
three regimes give. The analytic prediction from the Fisher informations is 0.929.

A first check on the reading, before any of it: an image study observes the whole region, so its
estimate has variance tau^2 + sigma^2 / M and standard deviation 0.164. The document's image-only
column is .165, .167 and .057, and 0.164/sqrt(8) is 0.058. That the baseline follows
arithmetically from the stated setup is what makes the rest worth running.

Environment knobs: REPS, REGIMES.
"""
import os
import sys

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.stats import chi2, norm

sys.path.insert(0, "/home/user/NiMARE")
from nimare.meta.cbma.blocks import block_loglik  # noqa: E402

TRUE_MEAN = 0.4
BETWEEN_SD = 0.15
ELEMENT_SD = 0.2
ELEMENTS = 9
THRESHOLD = 0.75

BETWEEN_VAR = BETWEEN_SD**2
ELEMENT_VAR = ELEMENT_SD**2
#: An image study observes every element, so it estimates the regional mean this precisely.
IMAGE_SD = np.sqrt(BETWEEN_VAR + ELEMENT_VAR / ELEMENTS)
CUTOFF = chi2.ppf(0.95, 1) / 2.0

GRID = np.linspace(-0.6, 1.4, 2001)
#: Offsets at which the reported maximum's log-density is tabulated. The height is a location
#: family in the mean -- proved in proofs/block_peak_selection.py -- so one curve in h - m serves
#: every (height, candidate mean) pair, which is what makes this affordable.
#:
#: The curve is then read by a cubic spline rather than linear interpolation. Linear
#: interpolation on this spacing carries an error of about 2e-6 per point, which accumulates to
#: 3e-5 over a handful of heights -- enough to fail the exactness check below for a reason that
#: has nothing to do with the location-family identity it is testing. A cubic spline drops that
#: to the point where the check is testing the identity again.
OFFSETS = np.linspace(THRESHOLD - TRUE_MEAN - 2.0, 4.0, 4001)


def tabulate_height_density():
    """Tabulate the reported maximum's log-density against ``h - m``, from the module itself."""
    values = np.array(
        [
            block_loglik(
                0.0,
                BETWEEN_VAR,
                np.array([ELEMENT_VAR]),
                heights=np.array([offset]),
                thresholds=np.array([-np.inf]),
                elements=np.array([ELEMENTS]),
            )
            for offset in OFFSETS
        ]
    )
    if not np.all(np.isfinite(values)):
        raise AssertionError(
            "the tabulated height density underflows inside the offsets the grid can reach; "
            "widen the tabulation rather than interpolating through an infinity"
        )
    return CubicSpline(OFFSETS, values)


def tabulate_block_curves():
    """Log-probability of one silent block, and of one report without its height, over the grid."""
    silent, indicator = [], []
    for mean in GRID:
        silent.append(
            block_loglik(
                mean,
                BETWEEN_VAR,
                np.array([ELEMENT_VAR]),
                heights=np.array([np.nan]),
                thresholds=np.array([THRESHOLD]),
                elements=np.array([ELEMENTS]),
            )
        )
        indicator.append(
            block_loglik(
                mean,
                BETWEEN_VAR,
                np.array([ELEMENT_VAR]),
                heights=np.array([THRESHOLD]),
                thresholds=np.array([THRESHOLD]),
                elements=np.array([ELEMENTS]),
                use_heights=False,
            )
        )
    return np.asarray(silent), np.asarray(indicator)


def draw(n_images, n_coordinates, rng):
    """One replication: image regional means, and each coordinate study's reporting outcome."""
    image_means = rng.normal(TRUE_MEAN, IMAGE_SD, size=n_images)
    shared = rng.normal(0.0, BETWEEN_SD, size=n_coordinates)
    region = TRUE_MEAN + shared[:, None] + rng.normal(
        0.0, ELEMENT_SD, size=(n_coordinates, ELEMENTS)
    )
    peaks = region.max(axis=1)
    return image_means, peaks[peaks > THRESHOLD], int(np.sum(peaks <= THRESHOLD))


def interval_from_profile(profile):
    best = int(np.argmax(profile))
    indices = np.flatnonzero(profile >= profile[best] - CUTOFF)
    touched = bool(indices[0] == 0 or indices[-1] == GRID.size - 1)
    return float(GRID[best]), float(GRID[indices[0]]), float(GRID[indices[-1]]), touched


def assert_interpolation_is_exact(heights, density_table):
    """Check the tabulated height density against direct evaluation, rather than trust it."""
    if heights.size == 0:
        return
    sample = heights[: min(8, heights.size)]
    for mean in (0.2, 0.4, 0.6):
        direct = block_loglik(
            mean,
            BETWEEN_VAR,
            np.full(sample.size, ELEMENT_VAR),
            heights=sample,
            thresholds=np.full(sample.size, -np.inf),
            elements=np.full(sample.size, ELEMENTS),
        )
        interpolated = float(density_table(sample - mean).sum())
        gap = abs(direct - interpolated)
        if gap > 1e-5 * max(abs(direct), 1.0):
            raise AssertionError(
                f"the tabulated height density disagrees with direct evaluation by {gap:.2e}; "
                "the location-family shortcut is wrong, not the table"
            )


def run(n_images, n_coordinates, reps, seed, density_table, silent_curve, indicator_curve):
    rng = np.random.default_rng(seed)
    arms = ("images", "heights", "indicator")
    rows = {name: {"error": [], "covered": [], "touched": 0} for name in arms}
    checked = False

    for _ in range(reps):
        image_means, heights, n_silent = draw(n_images, n_coordinates, rng)
        if not checked:
            assert_interpolation_is_exact(heights, density_table)
            checked = True

        mean = float(image_means.mean())
        half = norm.ppf(0.975) * IMAGE_SD / np.sqrt(n_images)
        rows["images"]["error"].append(mean - TRUE_MEAN)
        rows["images"]["covered"].append(mean - half <= TRUE_MEAN <= mean + half)

        image_profile = norm.logpdf(
            image_means[:, None], loc=GRID[None, :], scale=IMAGE_SD
        ).sum(axis=0)
        base = image_profile + n_silent * silent_curve

        reported_term = density_table(heights[None, :] - GRID[:, None]).sum(axis=1)
        for arm, profile in (
            ("heights", base + reported_term),
            ("indicator", base + heights.size * indicator_curve),
        ):
            estimate, low, high, touched = interval_from_profile(profile)
            rows[arm]["error"].append(estimate - TRUE_MEAN)
            rows[arm]["covered"].append(low <= TRUE_MEAN <= high)
            rows[arm]["touched"] += int(touched)

    out = {}
    for name, row in rows.items():
        errors = np.asarray(row["error"])
        covered = np.asarray(row["covered"], dtype=float)
        rmse = float(np.sqrt(np.mean(errors**2)))
        out[name] = {
            "rmse": rmse,
            "rmse_se": float(np.std(errors**2) / (2 * max(rmse, 1e-12) * np.sqrt(reps))),
            "coverage": float(covered.mean()),
            "touched": row["touched"],
            "squared": errors**2,
        }

    # The ratio of the two modelled arms is what the document actually concludes from, and it is
    # a quotient of two *correlated* random quantities -- both arms see the same replication. It
    # therefore needs its own standard error, from a bootstrap that resamples replications and
    # so keeps the pairing, rather than a point value compared against a band.
    bootstrap = np.random.default_rng(seed + 90000)
    heights_squared = out["heights"]["squared"]
    indicator_squared = out["indicator"]["squared"]
    draws = np.empty(2000)
    for index in range(draws.size):
        picked = bootstrap.integers(0, reps, reps)
        draws[index] = np.sqrt(heights_squared[picked].mean()) / np.sqrt(
            indicator_squared[picked].mean()
        )
    out["ratio"] = {
        "value": float(np.sqrt(heights_squared.mean()) / np.sqrt(indicator_squared.mean())),
        "se": float(draws.std()),
    }
    return out


TARGETS = {
    (1, 20): {"images": 0.165, "heights": 0.050, "indicator": 0.053},
    (1, 500): {"images": 0.167, "heights": 0.010, "indicator": 0.011},
    (8, 100): {"images": 0.057, "heights": 0.022, "indicator": 0.024},
}


if __name__ == "__main__":
    reps = int(os.environ.get("REPS", 1000))
    regimes = [(1, 20), (1, 500), (8, 100)]
    wanted = os.environ.get("REGIMES")
    if wanted:
        regimes = [regimes[int(index)] for index in wanted.split(",")]

    print(f"Image-only standard deviation implied by the stated setup: {IMAGE_SD:.4f}")
    print("Tabulating the block curves from the module ...")
    density_table = tabulate_height_density()
    silent_curve, indicator_curve = tabulate_block_curves()

    print(f"Calibration against the design document's section 8.2, {reps} replications each.")
    print(f"{'regime':>10} {'arm':>10} {'rmse':>16} {'target':>8} {'coverage':>9}  verdict")
    verdicts = []
    for index, (n_images, n_coordinates) in enumerate(regimes):
        result = run(
            n_images, n_coordinates, reps, 2000 + index,
            density_table, silent_curve, indicator_curve,
        )
        for arm in ("images", "heights", "indicator"):
            row = result[arm]
            target = TARGETS[(n_images, n_coordinates)][arm]
            rounding = 0.001 / np.sqrt(12)
            document_se = row["rmse_se"] * np.sqrt(1000 / reps)
            difference_se = np.sqrt(row["rmse_se"] ** 2 + document_se**2 + rounding**2)
            gap = abs(row["rmse"] - target) / max(difference_se, 1e-9)
            verdict = "ok" if gap < 3 else f"OFF BY {gap:.1f} se"
            verdicts.append(verdict)
            print(
                f"{n_images:4d}+{n_coordinates:<5d} {arm:>10} "
                f"{row['rmse']:.4f} +- {row['rmse_se']:.4f} {target:>8.3f} "
                f"{row['coverage']:>9.4f}  {verdict}"
                + (f"  [{row['touched']} touched the grid]" if row["touched"] else "")
            )
        ratio, ratio_se = result["ratio"]["value"], result["ratio"]["se"]
        # The analytic prediction from the Fisher informations, and the document's own three
        # regimes, both sit in 0.90-0.95. A miss is a miss only if it is outside that band by
        # more than the ratio's own bootstrap standard error allows.
        distance = max(0.90 - ratio, ratio - 0.95, 0.0)
        band = "ok" if distance <= 3 * ratio_se else f"OUTSIDE THE BAND BY {distance / ratio_se:.1f} se"
        print(
            f"{'':>10} {'ratio':>10} {ratio:.4f} +- {ratio_se:.4f} {'0.929':>8} {'':>9}  {band}"
        )
        verdicts.append(band)

    print()
    if all(verdict == "ok" for verdict in verdicts):
        print("PASS: the block peak-selection likelihood reproduces both modelled columns, and")
        print("the height-to-indicator ratio lands where the Fisher informations said it would.")
    else:
        print("FAIL: per the stated kill condition this stops anything built on the block")
        print("model until it is resolved. The defect is in the implementation, in my reading")
        print("of the setup, or in the table, and which one is not something to guess at.")
