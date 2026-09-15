"""A one-dimensional testbed for intuitions about thresholded, selectively reported data.

Every claim this project has argued over -- whether peak height or peak count carries the signal,
whether a truncated normal is the right selection model, whether the effect scale and the
reporting threshold are separately identifiable -- is a claim about the statistics of a smooth
random field that has been thresholded and reduced to a few reported maxima. None of it is about
brains. In one dimension, on a few thousand grid points, each of those questions answers in
milliseconds instead of the tens of minutes a whole-brain fit costs, and the truth is known
exactly rather than estimated from held-out subjects.

The trade is real and worth stating: a 1-D field has different topology from a 3-D one, so the
*numbers* here will not match a brain. Resel counts, peak densities and the Euler characteristic
all differ. What transfers is the *direction and mechanism* of an effect. So the rule for using
this file is: form the intuition here, then confirm the one that survives on real data. To keep
that honest, `calibration()` reproduces a result already established in 3-D, and anything this
testbed says should be discounted if that check drifts.

Run with no arguments to execute every test.
"""
import sys
import numpy as np
from scipy import optimize, stats
from scipy.ndimage import gaussian_filter1d

LENGTH = 4096
SMOOTH = 4.0          # grid points; the field's FWHM is about 2.35 times this
rng = np.random.default_rng(0)


#: Standard deviation of the smoothing filter's output for unit-variance white input. Computed
#: once from the filter's own coefficients rather than from each realisation, because dividing a
#: realisation by its *own* empirical standard deviation gives the field a Student-t marginal on
#: roughly the effective degrees of freedom. That is a real defect -- the tail is what every
#: question here depends on -- though it turned out not to be the one that broke the calibration.
_impulse = np.zeros(LENGTH)
_impulse[LENGTH // 2] = 1.0
NOISE_SD = float(np.sqrt((gaussian_filter1d(_impulse, SMOOTH, mode="wrap") ** 2).sum()))


def noise(n=1):
    """Unit-variance smooth Gaussian noise, `n` independent realisations."""
    raw = gaussian_filter1d(rng.standard_normal((n, LENGTH)), SMOOTH, axis=-1, mode="wrap")
    return raw / NOISE_SD


def truth_field(n_bumps=12, width=30.0):
    """A latent effect field: bumps of widely varying height, so strata exist to compare."""
    x = np.arange(LENGTH)
    field = np.zeros(LENGTH)
    centres = rng.choice(LENGTH, n_bumps, replace=False)
    heights = np.exp(rng.uniform(np.log(0.1), np.log(1.2), n_bumps))
    for c, h in zip(centres, heights):
        d = np.minimum(np.abs(x - c), LENGTH - np.abs(x - c))
        field += h * np.exp(-0.5 * (d / width) ** 2)
    return field


def peaks_above(z, cut):
    """Indices of local maxima of |z| at or above `cut`."""
    mag = np.abs(z)
    interior = (mag[1:-1] >= mag[:-2]) & (mag[1:-1] >= mag[2:]) & (mag[1:-1] >= cut)
    return np.flatnonzero(interior) + 1


def report(z, cut, per_cluster=False):
    """What a paper prints: every local maximum, or one representative per supra-threshold run.

    The distinction matters more than the threshold does. Reporting every maximum thins the
    pattern as the cut rises but leaves each reported value free to vary. Reporting one
    representative per cluster pins every value near the top of whatever cluster it came from,
    which is what destroys the magnitude channel -- and it is what cluster-extent reporting does.
    """
    if not per_cluster:
        return peaks_above(z, cut)
    above = np.abs(z) >= cut
    if not above.any():
        return np.array([], dtype=int)
    edges = np.flatnonzero(np.diff(above.astype(np.int8)))
    starts = np.r_[0 if above[0] else [], edges[::2] + 1 if above[0] else edges[0::2] + 1]
    bounds = np.flatnonzero(np.diff(np.r_[False, above, False].astype(np.int8)))
    out = []
    for lo, hi in zip(bounds[0::2], bounds[1::2]):
        out.append(lo + int(np.argmax(np.abs(z[lo:hi]))))
    return np.array(out, dtype=int)


# --------------------------------------------------------------- calibration against 3-D
def calibration():
    """Reproduce, in 1-D, a result already established in 3-D: the truncated-normal MLE

    applied to *peaks* under-recovers a small mean and over-recovers a large one, while the same
    MLE applied to ordinary supra-threshold *draws* is roughly right. In 3-D at a cut of 3.29 the
    peak MLE returned 0.26 for a true 0.5. If this testbed cannot show the same direction, nothing
    else it says should be believed.
    """
    print("calibration: truncated-normal MLE, peaks versus draws (3-D gave 0.26 for a true 0.5)")
    print(f"  {'true mean':>10} {'cut':>5} {'mean peak':>10} {'MLE on peaks':>13} "
          f"{'MLE on draws':>13} {'eff. draws':>11}")
    for mu in (0.5, 2.0):
        for cut in (2.5, 3.29):
            found, drawn = [], []
            for _ in range(2500):
                z = mu + noise()[0]
                # Signed, not absolute. The likelihood below is for values drawn from
                # N(m, 1) and kept when |x| >= cut, so folding the negative exceedances onto
                # the positive side inflates the mean -- by about 0.19 at a true 0.5 and a cut
                # of 3.29, where 2.8% of exceedances are negative, and by nothing at a true 2.0
                # where essentially none are. That asymmetry was the whole discrepancy.
                found.append(z[peaks_above(z, cut)])
                drawn.append(z[np.abs(z) >= cut])
            found = np.concatenate(found)
            drawn = np.concatenate(drawn)

            def mle(values):
                def neg(m):
                    surv = max(stats.norm.sf(cut - m) + stats.norm.cdf(-cut - m), 1e-12)
                    return -(stats.norm.logpdf(values, m, 1.0).sum()
                             - values.size * np.log(surv))
                return optimize.minimize_scalar(neg, bounds=(-4, 6), method="bounded").x

            sub = lambda a: a if a.size < 8000 else rng.choice(a, 8000, replace=False)
            effective = drawn.size / (SMOOTH * np.sqrt(4 * np.pi))
            print(f"  {mu:10.1f} {cut:5.2f} {found.mean():10.3f} "
                  f"{mle(sub(found)):13.3f} {mle(sub(drawn)):13.3f} {effective:11.0f}")
    print()


# --------------------------------------------------------- T1: which channel carries the signal
def channel_sweep():
    """Height versus count as the reporting threshold slides, so density varies continuously.

    The 3-D runs sampled two regimes and found the answer flipped between them. Here the cut moves
    smoothly, so the crossover itself is visible and can be located rather than bracketed.
    """
    print("T1: which channel predicts the truth, by reporting rule and density")
    print(f"  {'rule':>14} {'cut':>5} {'foci/study':>11} {'r, height':>10} {'r, count':>9} "
          f"{'winner':>8}")
    mu = truth_field()
    n_studies, n_subj = 20, 25
    smooth_kernel = 12.0
    for per_cluster in (False, True):
      for cut in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0):
        height_num = np.zeros(LENGTH)
        height_den = np.zeros(LENGTH)
        count = np.zeros(LENGTH)
        total = 0
        for _ in range(n_studies):
            z = mu * np.sqrt(n_subj) + noise()[0]
            idx = report(z, cut, per_cluster)
            total += idx.size
            if not idx.size:
                continue
            spike = np.zeros(LENGTH)
            spike[idx] = np.abs(z[idx]) / np.sqrt(n_subj)       # reported effect size
            reach = np.zeros(LENGTH)
            reach[idx] = 1.0
            height_num += gaussian_filter1d(spike, smooth_kernel, mode="wrap")
            height_den += gaussian_filter1d(reach, smooth_kernel, mode="wrap")
            near = np.zeros(LENGTH)
            near[idx] = 1.0
            count += (gaussian_filter1d(near, smooth_kernel, mode="wrap") > 1e-6)
        label = "one per cluster" if per_cluster else "every maximum"
        use = height_den > 1e-9
        if use.sum() < 100:
            print(f"  {label:>14} {cut:5.2f} {total / n_studies:11.1f}   too few reported")
            continue
        height = np.zeros(LENGTH)
        height[use] = height_num[use] / height_den[use]
        r_h = stats.pearsonr(height[use], mu[use])[0]
        r_c = stats.pearsonr(count[use], mu[use])[0]
        print(f"  {label:>14} {cut:5.2f} {total / n_studies:11.1f} {r_h:10.3f} {r_c:9.3f} "
              f"{('count' if r_c > r_h else 'height'):>8}")
      print()


# ------------------------------------------- T3: is the effect scale identifiable from the cut?
def identifiability():
    """Can the effect scale and the reporting threshold be told apart from coordinates alone?

    The claim borrowed from ecology was that they cannot -- that intensity and detectability are
    confounded, and unbiased data is needed to separate them. That claim was overstated when
    repeated from secondary sources, so it is tested here directly.

    The shape of the effect is treated as known and only its scale `a` is unknown, alongside the
    cut `u`. If the profile of the log-likelihood in `a` is flat, the two are confounded; if it
    has a well-separated minimum, they are not. Then one "image" study -- a direct noisy
    observation of the field, which is what an unthresholded map is -- is added, and the profile
    re-examined.
    """
    print("T3: is the effect scale identifiable from coordinates alone?")
    shape = truth_field()
    shape = shape / shape.max()
    a_true, u_true, n_subj, n_studies = 0.8, 3.0, 25, 30

    tables = []
    for _ in range(n_studies):
        z = a_true * shape * np.sqrt(n_subj) + noise()[0]
        tables.append(peaks_above(z, u_true))

    def neg_loglik(a, u):
        """Poisson likelihood for the peak pattern, intensity from the expected-maxima rate.

        The rate of local maxima above `u` in a smooth Gaussian field with local mean `m` is
        proportional to the zero-mean rate evaluated at `u - m`, which is the approximation the
        random field literature uses at high thresholds.
        """
        m = a * shape * np.sqrt(n_subj)
        rate = np.exp(-0.5 * np.maximum(u - m, 0.0) ** 2) + 1e-12
        rate = rate / rate.sum()
        total = sum(t.size for t in tables)
        ll = 0.0
        for t in tables:
            ll += np.log(rate[t]).sum()
        return -(ll - total * np.log(1.0))

    print(f"  {'':>16} {'a = 0.4':>9} {'a = 0.6':>9} {'a = 0.8':>9} {'a = 1.0':>9} "
          f"{'a = 1.5':>9}")
    grid_a = (0.4, 0.6, 0.8, 1.0, 1.5)
    row = []
    for a in grid_a:
        # Profile over the cut: the best the model can do at this scale.
        best = min(neg_loglik(a, u) for u in np.linspace(1.5, 5.0, 36))
        row.append(best)
    base = min(row)
    print(f"  {'coordinates only':>16} " + " ".join(f"{v - base:9.1f}" for v in row))

    # One image study: a direct noisy observation of the same field.
    image = a_true * shape + noise()[0] / np.sqrt(n_subj)
    sigma = 1.0 / np.sqrt(n_subj)
    row_img = []
    for a in grid_a:
        best = min(neg_loglik(a, u) for u in np.linspace(1.5, 5.0, 36))
        gauss = 0.5 * (((image - a * shape) / sigma) ** 2).sum()
        row_img.append(best + gauss)
    base_img = min(row_img)
    print(f"  {'plus one image':>16} " + " ".join(f"{v - base_img:9.1f}" for v in row_img))
    print("\n  Values are the increase in negative log-likelihood from the best scale; the true"
          "\n  scale is 0.8. A flat row means the scale is not identifiable from that data.\n")




# ------------------------------- T4: a fitted detection curve versus an assumed hard threshold
def detection_curve():
    """Does jointly estimating a soft detection curve beat inferring a hard cut?

    The estimator infers each study's reporting threshold from its smallest reported value, which
    is both a hard cut and a function of that study's signal. Seismology does not do this: Ogata
    and Katsura write the observed magnitude distribution as the true law times a smooth detection
    probability and estimate both together, so completeness is a derived quantity rather than an
    input.

    Here each study gets its own reporting behaviour drawn from one of three conventions, which is
    what a real literature search returns. Three estimators of the common effect scale are
    compared: one told the true per-study cut, one inferring it from the smallest reported value
    as the estimator currently does, and one fitting a shared logistic detection curve jointly
    with the scale. The question is which tracks the truth when the conventions are mixed.
    """
    print("T4: recovering the effect scale under mixed reporting conventions")
    shape = truth_field()
    shape = shape / shape.max()
    a_true, n_subj, n_studies = 0.8, 25, 30

    tables, true_cuts = [], []
    for k in range(n_studies):
        # Three conventions: a permissive cut, a strict one, and a cluster-extent-like rule that
        # keeps only the top of each supra-threshold run at a low forming cut.
        convention = k % 3
        z = a_true * shape * np.sqrt(n_subj) + noise()[0]
        if convention == 0:
            cut, idx = 2.5, report(z, 2.5, per_cluster=False)
        elif convention == 1:
            cut, idx = 4.5, report(z, 4.5, per_cluster=False)
        else:
            cut, idx = 3.0, report(z, 3.0, per_cluster=True)
        if idx.size:
            tables.append((idx, np.abs(z[idx])))
            true_cuts.append(cut)

    def rate_of(m, u, soft):
        """Density of reported maxima: the zero-mean maxima rate, shifted by the local mean.

        A local maximum above `u` in a field whose mean is `m` there is as likely as one above
        `u - m` in a zero-mean field, which is the high-threshold approximation. `soft` replaces
        the hard cut with a logistic reporting probability of that width.
        """
        if soft is None:
            return np.exp(-0.5 * np.maximum(u - m, 0.0) ** 2)
        return 1.0 / (1.0 + np.exp(-(m - u) / soft))

    def fit_scale(cuts, soft=None, shared_intensity=True):
        """Poisson likelihood in the scale, with the overall maxima density shared or per-study.

        This distinction is the whole test. An inhomogeneous Poisson log-likelihood is
        ``sum_i log lambda(x_i) - integral lambda``. Writing ``lambda_k = C_k * r_k(x)`` and
        letting each study have its own ``C_k`` profiles out to a multinomial over locations,
        which discards every count and leaves only the *shape* of the intensity -- deleting
        exactly the channel that was argued to carry the signal.

        With one ``C`` shared across studies the counts are back in play, and the scale is
        identified by how the count changes with each study's threshold and sample size. That is
        the difference between "how many foci did this study report" being a free parameter and
        being a prediction.
        """
        def neg(a):
            m = a * shape * np.sqrt(n_subj)
            rates = [rate_of(m, u, soft) for u in cuts]
            counts = [idx.size for idx, _ in tables]
            if not shared_intensity:
                total = 0.0
                for (idx, _), r in zip(tables, rates):
                    dens = r / max(r.sum(), 1e-12) + 1e-12
                    total += np.log(dens[idx]).sum()
                return -total
            # One shared C, profiled out analytically: C_hat = sum(counts) / sum_k sum_x r_k(x).
            integral = sum(float(r.sum()) for r in rates)
            n_total = sum(counts)
            c_hat = n_total / max(integral, 1e-12)
            total = 0.0
            for (idx, _), r in zip(tables, rates):
                total += np.log(c_hat * r[idx] + 1e-300).sum()
            total -= c_hat * integral
            return -total
        return optimize.minimize_scalar(neg, bounds=(0.1, 2.0), method="bounded").x

    inferred = [float(np.min(v)) for _, v in tables]
    print(f"  {'estimator':>34} {'intensity':>10} {'scale':>7} {'error':>7}")
    for label, cuts, soft, shared in (
        ("ORACLE true cut, hard", true_cuts, None, False),
        ("ORACLE true cut, hard", true_cuts, None, True),
        ("cut from smallest reported, hard", inferred, None, True),
        ("cut from smallest reported, soft", inferred, 0.75, True),
        ("ORACLE true cut, soft", true_cuts, 0.75, True),
    ):
        a = fit_scale(cuts, soft, shared)
        kind = "shared" if shared else "per-study"
        print(f"  {label:>34} {kind:>10} {a:7.3f} {a - a_true:+7.3f}")
    print(f"\n  true scale {a_true}. Mean true cut {np.mean(true_cuts):.2f}, "
          f"mean inferred cut {np.mean(inferred):.2f}.\n")




# ------------------------------------------------------- the rate function, measured not derived
_PEAK_CACHE = {}


def peak_height_law(n_fields=3000):
    """The zero-mean field's own peak-height law: maxima per point, and their height survival.

    T4 stalled on writing an expected-maxima rate by hand. That was never the blocker -- the
    generator is right here, so the law can be *measured* rather than approximated. Draw
    zero-mean fields, take every interior local maximum, and keep the sorted heights. Then the
    density of reported maxima where the field's mean is ``m`` and the cut is ``u`` is

        rho_max * [ Sbar(u - m) + Sbar(u + m) ]

    with ``Sbar`` the survival of that height law -- the upper term for maxima clearing ``+u``
    and the lower for minima clearing ``-u``, which two-sided reporting also keeps. This is exact
    for this field instead of being a high-threshold approximation, so a failure downstream is
    the architecture's and not the rate function's.

    In three dimensions the same law comes from random field theory (Cheng and Schwartzman, with
    the non-zero-mean extension of Zhao, Cheng and Schwartzman) or from the data's own smoothness,
    so measuring it here is a shortcut for the testbed, not a cheat that hides a hard step.
    """
    if n_fields in _PEAK_CACHE:
        return _PEAK_CACHE[n_fields]
    heights = []
    for _ in range(n_fields):
        z = noise()[0]
        interior = (z[1:-1] >= z[:-2]) & (z[1:-1] >= z[2:])
        heights.append(z[1:-1][interior])
    heights = np.sort(np.concatenate(heights))
    rho_max = heights.size / float(n_fields * LENGTH)
    _PEAK_CACHE[n_fields] = (heights, rho_max)
    return heights, rho_max


def survival(heights, value):
    """Fraction of peak heights above `value`, with a floor so the log-likelihood stays finite."""
    idx = np.searchsorted(heights, value, side="right")
    return np.maximum((heights.size - idx) / heights.size, 1e-9)


def joint_model():
    """Option C: one latent field, images as Gaussian observations, coordinates as a point process.

    The two failed schemes both treated a coordinate as a noisy measurement of the same number the
    image measures, and lost to the images alone. This does not: a coordinate study contributes
    only through *where and how often* it reported, via a rate that knows its threshold and its
    sample size, and never through the value it printed.

    Four estimators of the common scale, all given the same studies:

      * images only -- the benchmark the other schemes could not beat;
      * coordinates only, through the point process;
      * joint, images plus the point process;
      * ORACLE, the joint fit handed each study's true cut rather than inferring it.

    The intensity constant is shared across studies, without which the Poisson likelihood
    profiles out to a multinomial over locations and the counts -- the whole signal -- disappear.
    """
    print("T5: one latent field, images as Gaussians and coordinates as a point process")
    heights, rho_max = peak_height_law()
    shape_field = truth_field()
    shape_field = shape_field / shape_field.max()
    a_true = 0.8

    def draw(n_images, n_coords, seed):
        local = np.random.default_rng(seed)
        images, tables = [], []
        for _ in range(n_images):
            n_subj = int(local.integers(15, 40))
            obs = a_true * shape_field + noise()[0] / np.sqrt(n_subj)
            images.append((obs, n_subj))
        for j in range(n_coords):
            n_subj = int(local.integers(15, 40))
            z = a_true * shape_field * np.sqrt(n_subj) + noise()[0]
            cut, per_cluster = ((2.5, False), (4.5, False), (3.0, True))[j % 3]
            idx = report(z, cut, per_cluster)
            if idx.size:
                tables.append((idx, cut, n_subj, float(np.abs(z[idx]).min())))
        return images, tables

    def fit(images, tables, use_images, use_coords, oracle_cut):
        def neg(a):
            total = 0.0
            if use_images:
                for obs, n_subj in images:
                    sd = 1.0 / np.sqrt(n_subj)
                    total += 0.5 * (((obs - a * shape_field) / sd) ** 2).sum()
            if use_coords and tables:
                rates, counts = [], []
                for idx, cut, n_subj, inferred in tables:
                    u = cut if oracle_cut else inferred / np.sqrt(1.0)
                    m = a * shape_field * np.sqrt(n_subj)
                    r = rho_max * (survival(heights, u - m) + survival(heights, u + m))
                    rates.append(r)
                    counts.append(idx.size)
                integral = sum(float(r.sum()) for r in rates)
                n_total = sum(counts)
                c_hat = n_total / max(integral, 1e-12)      # one shared intensity, profiled out
                for (idx, _, _, _), r in zip(tables, rates):
                    total -= np.log(c_hat * r[idx] + 1e-300).sum()
                total += c_hat * integral
            return total
        return optimize.minimize_scalar(neg, bounds=(0.1, 2.0), method="bounded").x

    print(f"  {'images':>7} {'coords':>7} {'estimator':>26} {'scale':>7} {'error':>7}")
    for n_images, n_coords in ((2, 12), (2, 30), (5, 12)):
        rows = {}
        for seed in range(12):
            images, tables = draw(n_images, n_coords, seed)
            rows.setdefault("images only", []).append(fit(images, tables, True, False, False))
            rows.setdefault("coordinates only", []).append(fit(images, tables, False, True, False))
            rows.setdefault("joint", []).append(fit(images, tables, True, True, False))
            rows.setdefault("ORACLE joint, true cuts", []).append(
                fit(images, tables, True, True, True))
        for name in ("images only", "coordinates only", "joint", "ORACLE joint, true cuts"):
            a = float(np.mean(rows[name]))
            print(f"  {n_images:>7} {n_coords:>7} {name:>26} {a:7.3f} {a - a_true:+7.3f}")
        print()
    print(f"  true scale {a_true}. The joint beating images alone is what the two failed"
          f"\n  combination schemes could not do.\n")




def joint_field():
    """The real problem: recover the whole field, not one scalar with a known shape.

    T5 showed the *scale* is identifiable from coordinates once the intensity is modelled
    correctly. That is a much easier problem than the one that matters, because it assumes the
    spatial shape is already known and only its amplitude is in question. Here the shape is
    unknown: the field is expanded on a fixed set of evenly spaced bumps and every coefficient is
    fitted, from images alone, from coordinates alone, and from both.

    This is where joining should finally pay. Images are unbiased but few, so the field they
    recover is noisy; coordinates are many and say where effects repeatedly appear. If the point
    process contributes real spatial information, the joint fit beats the images by themselves --
    which is exactly what pooling, rescaled pooling and gating all failed to do.
    """
    print("T6: recovering the whole field, images and coordinates jointly")
    heights, rho_max = peak_height_law()
    x = np.arange(LENGTH)
    n_basis = 24
    centres = np.linspace(0, LENGTH, n_basis, endpoint=False)
    width = LENGTH / n_basis / 1.5
    basis = np.stack([
        np.exp(-0.5 * (np.minimum(np.abs(x - c), LENGTH - np.abs(x - c)) / width) ** 2)
        for c in centres
    ])
    truth = truth_field()

    def field_of(coef):
        return coef @ basis

    def draw(n_images, n_coords, seed):
        local = np.random.default_rng(seed)
        images, tables = [], []
        for _ in range(n_images):
            n_subj = int(local.integers(15, 40))
            images.append((truth + noise()[0] / np.sqrt(n_subj), n_subj))
        for j in range(n_coords):
            n_subj = int(local.integers(15, 40))
            z = truth * np.sqrt(n_subj) + noise()[0]
            cut, per_cluster = ((2.5, False), (4.5, False), (3.0, True))[j % 3]
            idx = report(z, cut, per_cluster)
            if idx.size:
                tables.append((idx, float(np.abs(z[idx]).min()), n_subj))
        return images, tables

    def fit(images, tables, use_images, use_coords, ridge=1e-2):
        def neg(coef):
            field = field_of(coef)
            total = ridge * float(coef @ coef)
            if use_images:
                for obs, n_subj in images:
                    total += 0.5 * (((obs - field) * np.sqrt(n_subj)) ** 2).sum()
            if use_coords and tables:
                rates = []
                for idx, u, n_subj in tables:
                    m = field * np.sqrt(n_subj)
                    rates.append(rho_max * (survival(heights, u - m)
                                            + survival(heights, u + m)))
                integral = sum(float(r.sum()) for r in rates)
                n_total = sum(idx.size for idx, _, _ in tables)
                c_hat = n_total / max(integral, 1e-12)
                for (idx, _, _), r in zip(tables, rates):
                    total -= np.log(c_hat * r[idx] + 1e-300).sum()
                total += c_hat * integral
            return total
        start = np.full(n_basis, 0.3)
        out = optimize.minimize(neg, start, method="Powell",
                                options={"maxiter": 4000, "xtol": 1e-3, "ftol": 1e-3})
        return field_of(out.x)

    print(f"  {'images':>7} {'coords':>7} {'estimator':>20} {'r':>7} {'rmse':>7} {'ratio':>7}")
    for n_images, n_coords in ((2, 30), (5, 30)):
        rows = {}
        for seed in range(6):
            images, tables = draw(n_images, n_coords, seed)
            for name, ui, uc in (("images only", True, False),
                                 ("coordinates only", False, True),
                                 ("joint", True, True)):
                est = np.abs(fit(images, tables, ui, uc))
                rows.setdefault(name, []).append((
                    stats.pearsonr(est, truth)[0],
                    float(np.sqrt(np.mean((est - truth) ** 2))),
                    float(est.mean() / max(truth.mean(), 1e-9)),
                ))
        for name in ("images only", "coordinates only", "joint"):
            r, rmse, ratio = np.mean(np.array(rows[name]), axis=0)
            print(f"  {n_images:>7} {n_coords:>7} {name:>20} {r:+7.3f} {rmse:7.3f} {ratio:7.2f}")
        print()
    print("  The joint beating images alone on r or rmse is what every scheme so far has failed"
          "\n  to do. A ratio near 1.00 means the scale came out right as well.\n")


def smoothness_sensitivity():
    """How wrong can the assumed smoothness be before the recovered scale suffers?

    The intensity model needs a peak-height law, and that law depends on the field's smoothness.
    In the testbed the law is measured from the generator, which is exactly right and therefore
    says nothing about robustness. On real data it would come from random field theory given an
    estimated smoothness, or from the images -- and no method exists for estimating it from
    coordinates alone, which is a real gap.

    So the question that decides whether the gap matters: mis-specify the smoothness deliberately
    and see what it costs. Data are generated at one smoothness and fitted with a peak-height law
    measured at another, from half to double. If the recovered scale barely moves, smoothness is
    a detail to estimate roughly; if it tracks the error, it is a blocker and needs its own method.
    """
    print("T7: what mis-specified smoothness costs the recovered scale")
    x = np.arange(LENGTH)
    a_true, n_studies = 0.8, 24
    shape_field = truth_field()
    shape_field = shape_field / shape_field.max()

    def noise_at(sigma, n=1):
        raw = gaussian_filter1d(rng.standard_normal((n, LENGTH)), sigma, axis=-1, mode="wrap")
        impulse = np.zeros(LENGTH)
        impulse[LENGTH // 2] = 1.0
        sd = np.sqrt((gaussian_filter1d(impulse, sigma, mode="wrap") ** 2).sum())
        return raw / sd

    def law_at(sigma, n_fields=1500):
        heights = []
        for _ in range(n_fields):
            z = noise_at(sigma)[0]
            interior = (z[1:-1] >= z[:-2]) & (z[1:-1] >= z[2:])
            heights.append(z[1:-1][interior])
        heights = np.sort(np.concatenate(heights))
        return heights, heights.size / float(n_fields * LENGTH)

    true_sigma = SMOOTH
    tables = []
    for j in range(n_studies):
        n_subj = int(rng.integers(15, 40))
        z = a_true * shape_field * np.sqrt(n_subj) + noise_at(true_sigma)[0]
        cut, per_cluster = ((2.5, False), (4.5, False), (3.0, True))[j % 3]
        idx = report(z, cut, per_cluster)
        if idx.size:
            tables.append((idx, float(np.abs(z[idx]).min()), n_subj))

    print(f"  {'assumed sigma':>14} {'vs truth':>9} {'maxima/point':>13} {'scale':>7} {'error':>7}")
    for factor in (0.5, 0.75, 1.0, 1.5, 2.0):
        heights, rho_max = law_at(true_sigma * factor)

        def neg(a):
            field = a * shape_field
            rates = []
            for idx, u, n_subj in tables:
                m = field * np.sqrt(n_subj)
                rates.append(rho_max * (survival(heights, u - m) + survival(heights, u + m)))
            integral = sum(float(r.sum()) for r in rates)
            n_total = sum(idx.size for idx, _, _ in tables)
            c_hat = n_total / max(integral, 1e-12)
            total = 0.0
            for (idx, _, _), r in zip(tables, rates):
                total -= np.log(c_hat * r[idx] + 1e-300).sum()
            return total + c_hat * integral

        a = optimize.minimize_scalar(neg, bounds=(0.1, 2.0), method="bounded").x
        print(f"  {true_sigma * factor:14.2f} {f'x{factor:g}':>9} {rho_max:13.4f} "
              f"{a:7.3f} {a - a_true:+7.3f}")
    print(f"\n  true scale {a_true}, true sigma {true_sigma}. Flat errors mean smoothness is a"
          f"\n  detail; errors tracking the factor mean it needs its own estimator.\n")


def separability():
    """Is the *variation* of intensity shape across studies what identifies the effect scale?

    CBMR's predictor is separable: one spatial log-intensity times a scalar per experiment, so
    every study shares a spatial shape and differs only in overall rate. The model proposed here
    is not separable -- the threshold and sample size enter inside the nonlinearity, so a small
    study with a lenient cut has a broad intensity while a large one with a strict cut has an
    intensity concentrated on the peaks of the same field.

    The claim is that this shape variation is the mechanism: it is what pins the effect scale,
    the way varying detection across observers pins the Gutenberg-Richter law. If so, removing the
    variation should destroy identification even though the model and the data-generating process
    are otherwise unchanged.

    Three collections, all with the same true field and the same number of studies. In the
    homogeneous ones every study shares one threshold and one sample size, so the intensity shape
    is identical across studies and the model is effectively separable. In the heterogeneous one
    they vary as a literature does.
    """
    print("T8: does intensity-shape variation across studies identify the scale?")
    heights, rho_max = peak_height_law()
    shape_field = truth_field()
    shape_field = shape_field / shape_field.max()
    a_true, n_studies = 0.8, 24

    def collect(kind, seed):
        local = np.random.default_rng(seed)
        tables = []
        for j in range(n_studies):
            if kind == "homogeneous":
                n_subj, cut, per_cluster = 25, 3.0, False
            elif kind == "varying N only":
                n_subj, cut, per_cluster = int(local.integers(12, 60)), 3.0, False
            else:
                n_subj = int(local.integers(12, 60))
                cut, per_cluster = ((2.5, False), (4.5, False), (3.0, True))[j % 3]
            z = a_true * shape_field * np.sqrt(n_subj) + noise()[0]
            idx = report(z, cut, per_cluster)
            if idx.size:
                tables.append((idx, cut, n_subj))
        return tables

    def fit(tables):
        def neg(a):
            rates = []
            for idx, u, n_subj in tables:
                m = a * shape_field * np.sqrt(n_subj)
                rates.append(rho_max * (survival(heights, u - m) + survival(heights, u + m)))
            integral = sum(float(r.sum()) for r in rates)
            n_total = sum(idx.size for idx, _, _ in tables)
            c_hat = n_total / max(integral, 1e-12)
            total = 0.0
            for (idx, _, _), r in zip(tables, rates):
                total -= np.log(c_hat * r[idx] + 1e-300).sum()
            return total + c_hat * integral
        # Curvature of the profile at the optimum: how sharply the scale is pinned.
        a = optimize.minimize_scalar(neg, bounds=(0.1, 2.0), method="bounded").x
        step = 0.05
        curvature = (neg(a + step) - 2 * neg(a) + neg(a - step)) / step ** 2
        return a, curvature

    print(f"  {'collection':>22} {'scale':>7} {'error':>7} {'curvature':>11} {'spread':>8}")
    for kind in ("homogeneous", "varying N only", "varying N and threshold"):
        scales, curves = [], []
        for seed in range(8):
            tables = collect(kind, seed)
            a, curvature = fit(tables)
            scales.append(a)
            curves.append(curvature)
        print(f"  {kind:>22} {np.mean(scales):7.3f} {np.mean(scales) - a_true:+7.3f} "
              f"{np.mean(curves):11.0f} {np.std(scales):8.3f}")
    print(f"\n  true scale {a_true}. Higher curvature is a more sharply identified scale. If the"
          f"\n  homogeneous collection is flat and the varying one is not, shape variation across"
          f"\n  studies is the mechanism, and a separable model cannot have it.\n")


if __name__ == "__main__":
    which = sys.argv[1:] or ["calibration", "channel_sweep", "identifiability", "detection_curve", "joint_model", "joint_field", "smoothness_sensitivity", "separability"]
    for name in which:
        globals()[name]()
