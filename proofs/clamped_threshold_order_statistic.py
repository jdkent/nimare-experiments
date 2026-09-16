r"""What does ``clamp_threshold`` do to the fit, given that it reads an order statistic?

``CBES`` defaults to ``clamp_threshold=True``, which lowers each study's *assumed* reporting
cutoff to that study's own smallest reported :math:`|z|`. The docstring justifies it as a bound:
anything reported cleared the cut, so :math:`c_k \le \min_j |z_{kj}|`.

The worry that prompted this proof: :math:`\min_j |z_{kj}|` is the **last order statistic** of a
study's reported heights, and it falls as the peak count rises, and the peak count rises with how
much signal the study had. So the *effective* assumed threshold floats with per-study signal --
the same mechanism the never-cap rule describes, reached through an order statistic rather than a
cap. Nothing is capped and the contamination path is still open.

The algebra says the worry is real in form and **backwards in consequence**, which is the whole
reason to write it down before touching code.

  * the overshoot :math:`\mathbb{E}[\hat c_k] - c_k = 1/(M_k\theta)` is strictly positive and
    strictly decreasing in the peak count, so the clamp does float with signal (claims 1-3);
  * an overstated cutoff **attenuates** both the score and the information a silence carries
    about :math:`\mu` (claims 4-7); every silence's score is negative, so attenuating it raises
    the root of the estimating equation (claim 8);
  * therefore the *unclamped* default, which overstates the cutoff by the full
    :math:`c_{\text{default}} - c_k`, biases :math:`\hat\mu` **upward more** than the clamp does.
    The clamp is a partial correction to that bias, not a new source of one (claim 9).

What survives of the worry is not a bias in :math:`\hat\mu` but a **differential weighting**: the
correction is larger for peak-rich studies, so silences are weighted by signal. That is harmless
for the direction of :math:`\hat\mu` and not harmless for anything that compares studies or reads
a spread across them. Declaring the threshold per study removes the overshoot *and* the
differential; the clamp removes part of the first and creates the second.

**Regime this cannot speak to.** Claim 1 models reported heights as independent exceedances over
the cut, which is the high-threshold tail limit for a *voxelwise* cut with no further selection.
A cluster-forming analysis reports the maxima of clusters that also passed an extent criterion,
which is a different order statistic under a selection this model does not contain. I do not know
the direction or size of the difference: a first attempt to check it numerically had no extent
selection in it and so measured nothing, and it is removed rather than left reading as evidence.
Until that is done properly, no claim here should be applied to a cluster-forming corpus.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

t, x, mu, sigma, c, alpha = sp.symbols("t x mu sigma c alpha", real=True)
theta = sp.Symbol("theta", positive=True)
M = sp.Symbol("M", positive=True)

SQRT2 = sp.sqrt(2)
phi = sp.exp(-t**2 / 2) / sp.sqrt(2 * sp.pi)
Phi = (1 + sp.erf(t / SQRT2)) / 2

proof = Proof(
    "clamped_threshold_order_statistic",
    "What clamping the assumed cutoff to the smallest reported height does to the fit",
    __doc__,
)

# ------------------------------------------------------------------ the overshoot

# Exceedances over a high cut are exponential in the tail limit. The minimum of M independent
# exceedances has the survival function of a single exponential with rate M*theta.
survival_of_one = sp.exp(-theta * x)
survival_of_min = sp.exp(-M * theta * x)
proof.define(r"P(X_j - c_k > x)", survival_of_one)
proof.claim(
    "the-minimum-exceedance-is-exponential-with-the-count-in-its-rate",
    sp.simplify(survival_of_one**M - survival_of_min),
    r"P(\min_j (X_j - c_k) > x) = e^{-M\theta x}",
)

overshoot = sp.integrate(survival_of_min, (x, 0, sp.oo))
proof.claim(
    "so-the-clamp-overshoots-the-true-cut-by-the-reciprocal-of-count-times-rate",
    sp.simplify(overshoot - 1 / (M * theta)),
    r"\mathbb{E}[\hat c_k] - c_k = \frac{1}{M_k\theta} > 0",
)
proof.claim(
    "and-the-overshoot-shrinks-as-the-study-reports-more-peaks",
    sp.simplify(sp.diff(1 / (M * theta), M) + 1 / (M**2 * theta)),
    r"\partial_M \frac{1}{M\theta} = -\frac{1}{M^2\theta} < 0",
)

# ------------------------------------------------------- what a silence carries about mu

# A silence at a voxel is the event that the study's effect fell below its cutoff, so on the
# standardised scale t = (c - mu)/sigma the silence probability is Phi(t).
t_of_mu = (c - mu) / sigma
silence_probability = Phi.subs(t, t_of_mu)
proof.define(r"P(\text{silent}) = \Phi(t)", silence_probability)

# Score of a silence: the inverse Mills ratio, and it is negative for every t, so every silence
# pushes mu down.
silence_score = sp.diff(sp.log(silence_probability), mu)
proof.claim(
    "a-silence-scores-minus-the-inverse-mills-ratio-so-it-always-pushes-mu-down",
    sp.simplify(silence_score + (phi / Phi).subs(t, t_of_mu) / sigma),
    r"\partial_\mu \log\Phi(t) = -\frac{1}{\sigma}\frac{\varphi(t)}{\Phi(t)} < 0",
)

# The push weakens monotonically as the assumed cutoff is raised: the Mills ratio's derivative is
# negative for every t, because t*Phi + phi is the (positive) truncated mean times Phi.
mills = phi / Phi
proof.claim(
    "raising-the-assumed-cutoff-weakens-that-push-at-every-t",
    sp.simplify(sp.diff(mills, t) + phi * (t * Phi + phi) / Phi**2),
    r"\partial_t \frac{\varphi}{\Phi} = -\frac{\varphi\,(t\Phi + \varphi)}{\Phi^2} < 0",
)

# Information a single silence indicator carries about mu: the probit information.
information = sp.diff(silence_probability, mu) ** 2 / (
    silence_probability * (1 - silence_probability)
)
target = (phi**2 / (Phi * (1 - Phi))).subs(t, t_of_mu) / sigma**2
proof.claim(
    "the-information-in-a-silence-is-the-probit-information",
    sp.simplify(information - target),
    r"I_\mu(\text{silence}) = \frac{1}{\sigma^2}\frac{\varphi(t)^2}{\Phi(t)(1-\Phi(t))}",
)

# That information is even in t, hence maximised at t = 0 and falling away from it in both
# directions: a cutoff pushed away from mu costs information whichever way it is pushed.
weight = phi**2 / (Phi * (1 - Phi))
proof.claim(
    "and-it-is-even-in-t-so-it-peaks-where-the-cutoff-sits-at-mu",
    sp.simplify(weight - weight.subs(t, -t)),
    r"g(t) = g(-t), \quad g(t) = \frac{\varphi^2}{\Phi(1-\Phi)}",
)
proof.claim(
    "with-a-log-derivative-that-is-exactly-the-mills-difference-minus-2t",
    sp.simplify(sp.diff(sp.log(weight), t) - (-2 * t + phi * (2 * Phi - 1) / (Phi * (1 - Phi)))),
    r"\partial_t \log g = -2t + \frac{\varphi\,(2\Phi-1)}{\Phi(1-\Phi)}",
)

# --------------------------------------------- attenuating the silences raises the fitted mu

# Caricature of the estimating equation: reported studies contribute A(mu), silent ones contribute
# -B(mu) with B > 0, and an overstated cutoff multiplies the silent part by alpha < 1.
A = sp.Function("A")(mu)
B = sp.Function("B")(mu)
score_total = A - alpha * B
implicit = sp.diff(score_total, mu)
root_shift = -sp.diff(score_total, alpha) / implicit
proof.claim(
    "attenuating-the-silences-moves-the-root-by-B-over-the-total-curvature",
    sp.simplify(root_shift - B / implicit),
    r"\frac{d\mu^*}{d\alpha} = \frac{B(\mu^*)}{\partial_\mu S} < 0",
)

# B > 0 and the total score decreases in mu at a maximum, so the derivative is negative: a smaller
# alpha -- a more overstated cutoff -- gives a larger mu*. The unclamped default overstates by the
# full c_default - c_k, the clamp by only 1/(M*theta), so the clamp's alpha is closer to one.
clamped_gap = 1 / (M * theta)
unclamped_gap = sp.Symbol("Delta", positive=True)
# The clamp's overshoot is min(Delta, 1/(M*theta)) and the unclamped default's is Delta, so the
# clamp's is never the larger. That is a property of a minimum, not a claim worth dressing as
# one, so it is recorded as a definition and the comparison is left to the reader.
gap_used = sp.Min(unclamped_gap, clamped_gap)
proof.define(r"\hat c_k - c_k \ \text{(clamped)}", gap_used)
proof.define(r"\hat c_k - c_k \ \text{(unclamped)}", unclamped_gap)

# --------------------------------------------------------------- numeric checks and calibration

if __name__ == "__main__":
    import numpy as np
    from scipy.special import ndtr
    from scipy.stats import norm

    proof.report()

    # The evenness claim fixes the peak at t = 0; monotonicity away from it is checked
    # numerically rather than asserted, because sympy will not order erf expressions.
    grid = np.linspace(0.01, 6.0, 600)
    g = norm.pdf(grid) ** 2 / (ndtr(grid) * (1 - ndtr(grid)))
    if not np.all(np.diff(g) < 0):
        raise AssertionError("the probit information is not monotone decreasing for t > 0")
    mills_grid = norm.pdf(grid) / ndtr(grid)
    if not np.all(np.diff(mills_grid) < 0):
        raise AssertionError("the inverse Mills ratio is not monotone decreasing")
    print("  [ok] numeric: probit information and the Mills ratio both fall monotonically for t>0")

    # Calibration against the corpus already measured at full scale. The pain collection carries
    # 267 published peaks across 21 studies, so about 12.7 peaks per study. The Gaussian tail
    # hazard at the default cut is theta = phi(u)/(1-Phi(u)).
    u = 3.2905267314919255
    hazard = norm.pdf(u) / (1 - ndtr(u))
    peaks_per_study = 267.0 / 21.0
    predicted = 1.0 / (peaks_per_study * hazard)
    print(f"  calibration: tail hazard at the default cut is {hazard:.3f}; at "
          f"{peaks_per_study:.1f} peaks per study the clamp overshoots by {predicted:.4f} z")
    if not predicted < 0.05:
        raise AssertionError(
            "the predicted overshoot is not small at the measured pain peak density, which "
            "contradicts the shipped claim that the clamp is bit-identical where no table "
            "contradicts the assumption"
        )
    print("  [ok] calibration: consistent with the shipped 'bit-identical where no table")
    print("       contradicts the assumption' behaviour for a voxelwise cut")

    print()
    print("Outcome: the flagged mechanism is real -- the clamp's assumed cutoff does float with")
    print("per-study signal -- but its effect on mu-hat's bias is the helpful direction, because")
    print("the unclamped default overstates the cutoff by more and an overstated cutoff")
    print("attenuates the only terms that push mu down. What remains is a differential weighting")
    print("of silences by signal, which no claim here shows to be benign. Declaring the")
    print("threshold per study removes the overshoot and the differential together.")
