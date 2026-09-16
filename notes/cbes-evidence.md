# CBES: the evidence behind the estimator

This is the measurement record that used to live in the `CBES` class docstring in
`nimare/meta/cbma/effectsize.py`. It was moved here because a 903-line docstring -- a quarter of
the module -- is a research log in the wrong place: a user opening `help(CBES)` needs to know
what the maps mean and which ones to trust, not the autopsy of every route that failed on the way.

Nothing is abridged. Every table, every retraction and every failed route is reproduced verbatim
below, and the docstring now carries the conclusions with a pointer here. The running log with
dates and reasoning is `cbes-open-program.md`; this file is the curated version, organised the way
the docstring organised it.

The scripts that produced each number are in `experiments/`, and the symbolic proofs behind the
identifiability and asymptote claims are in `proofs/`.

---

    Notes
    -----
    Where ALE and (M)KDA ask *where do studies agree something happened*, this asks *how big is
    the effect there*, and answers it from two channels that are deliberately unlike each other:

    **Images give the magnitude.** Each voxel pools the studies supplying ``g``/``g_var`` maps
    by inverse variance, with a local DerSimonian-Laird :math:`\tau^2`:

    .. math::

        \hat{g}(v) = \frac{\sum_k g_k(v) / (s^2_k(v) + \tau^2(v))}
                          {\sum_k 1 / (s^2_k(v) + \tau^2(v))}

    which with one study is that study's map and with several is a textbook random-effects
    meta-analysis. No spatial kernel: an image already says what it says at every voxel.

    **Coordinates give the selection.** Every study on the collection's roster that reported no
    focus within ``coverage_radius`` of :math:`v` enters a zero-inflated censored likelihood
    there, contributing the probability that a study of its sample size, applying its reporting
    threshold, would have stayed silent. Two quantities come out of that fit, and keeping them
    apart is the point of it: :math:`\pi(v)`, the fraction of studies with a non-null effect
    here, and :math:`\mu(v)`, the effect size *among the studies that have one*. A study silent
    in a region either has no effect there or has one that failed to clear its threshold, and
    the mixture lets the data decide, so silence need not be explained as a small-but-real
    common effect -- which is what drags a plain Tobit fit below the truth. Nothing is imputed
    :footcite:p:`tench2017coordinate`.

    Available maps:

    ============== ===============================================================
    "g"            Pooled Hedges' g among the studies with an effect, on the
                   effect-size scale the images arrive on.
    "prevalence"   Fitted fraction of studies with a non-null effect here. Added
                   under the zero-inflated selection model. Read ordinally; see
                   Warnings.
    "g_marginal"   ``g`` times ``prevalence``: the effect averaged over *all*
                   studies rather than over those that have one, which is the
                   estimand an image-based meta-analysis reports. Added alongside
                   ``prevalence``.
    "se_marginal"  Standard error of ``g_marginal``, by the delta method on the
                   same observed information. Zero where there is none.
    "coordinate\_  Share of the Fisher information about ``g`` that came from the
     share"        coordinate tables rather than from the images' values, in [0, 1].
                   **Read this before trusting the coordinate channel anywhere.**
                   Added under the zero-inflated selection model.
    "se"           Standard error of the pooled estimate. See ``se_method``.
    "z"            ``g / se``. Two-tailed.
    "p", "logp"    p-value for ``z``, and its ``-log10``.
    "tau2"         Local between-study variance.
    "n_studies"    Number of image studies contributing a value here.
    "n_eff"        Kish effective number of image studies, ``(sum w)^2 / sum w^2``.
    "dof"          Degrees of freedom to refer ``se`` to; see below.
    ============== ===============================================================

    Why the indicator and not the heights
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    An earlier design pooled the reported peak heights as effect sizes alongside the images,
    with a spatial kernel and a per-study peak-height correction. Stratifying its error by how
    many studies reported near a voxel showed the two channels pulling opposite ways: against a
    held-out reference the bias reduction from adding a coordinate corpus to a two-image
    collection was **54% where no focus reached the voxel and 6% where two studies reported**,
    shrinking monotonically as more reported. That is the signature of the reporting *pattern*
    helping and the magnitudes hurting, not of a magnitude channel working.

    On the 21-study NIDM pain collection, split in half so the reference comes from studies the
    coordinates never touched, eight splits scored paired, and run twice: once on the tables
    **the papers actually printed** -- the collection's own 267 transcribed peaks -- and once on
    tables re-extracted from the same studies' maps the way a paper would produce them
    (cluster-forming cut, whole clusters kept, one focus per cluster).

    ==========================  ======  ========  ==========  =====  ======
    estimate                      bias  at top      rmse      rank r  AUC
    ==========================  ======  ========  ==========  =====  ======
    images only, pooled         +0.136    +0.137       0.269   0.484   0.893
    CBES with the silence off   +0.142    +0.155       0.272   0.499   0.899
    CBES ``g``, published        +0.069    +0.005       0.230   0.480   0.886
    CBES ``g``, extracted        +0.075    -0.065       0.210   0.486   0.889
    CBES ``g_marginal``, publ.   -0.064    -0.291       0.184   0.483   0.887
    ==========================  ======  ========  ==========  =====  ======

    **The coordinate channel corrects the level and leaves the pattern alone, and this holds on
    real tables rather than only on a proxy for them.** On the published coordinates ``rmse``
    falls 14% against pooling the images alone (paired p = 0.0024) and the mean bias 49%
    (p = 0.0001), while the +0.137 overestimate at the strongest voxels becomes **+0.005** --
    the best-centred top stratum of any arm measured here. The ordering barely moves: rank
    correlation -0.004 (p = 0.84) and AUC -0.006 (p = 0.30). Pearson ``r`` costs 0.030 and that
    cost is no longer significant (p = 0.075). Read it as a correction to the magnitude, not as
    a better map.

    Running both table sources matters because they are not interchangeable, which was measured
    rather than assumed: the cluster scheme recovers only 23% of pain's published peaks within
    8 mm and 42% within 20 mm, and finds 119 peaks where the papers printed 267. So the
    extracted arm is a different and sparser input, not a faithful copy -- it reads more voxels
    as silent than the literature does. That it gives a *larger* rmse gain (22%, p = 0.0003)
    and a worse-centred top stratum (-0.065) is consistent with that: more silence, more
    shrinkage. The published-table row is the one to quote.

    The reason the heights were never going to work is in the input rather than the fit.
    Regressing a held-out truth at a focus on the effect size that focus's own table reports
    gives a slope of 0.08 to 0.18 with most of the value in the intercept: **one tabulated
    coordinate explains 5% to 9% of the variance in the effect at its own location.** And the
    level is not even a property of the studies -- holding the studies fixed and changing only
    how a paper would have tabulated them, the ratio of the old ``g`` to the held-out truth ran
    from 0.82 to 2.29 across FDR, voxelwise FWE and cluster-extent thresholding and across
    tabulating a cluster by its maximum or its centre of mass. Measured and rejected as
    remedies before the heights were dropped: the truncated-normal selection correction, which
    is the wrong event for a local maximum and returns 0.26 for a true 0.5; a per-study
    peak-height rescaling; subtracting the censoring floor; reporting by centre of mass; and
    widening the assumed cluster.

    Removing the heights removes what came with them. There is no longer an unidentified
    overall scale -- the images arrive on the effect-size scale, so ``g`` is in the units it
    claims and the ``g_relative``/``g_absolute`` distinction is gone. There is no kernel width
    to choose, no peak-height correction to calibrate, and no threshold to infer.

    Both values of the indicator, over different extents
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A silence and a report are complementary events -- :math:`|g| < c` and :math:`|g| \ge c`
    for the study's own cut :math:`c` -- and **both are evidence**. Dropping the report limb
    leaves the silences as the only evidence about the indicator, so the model reads the
    observed silence fraction against a denominator that excludes every study that reported,
    and the magnitude shrinks past the threshold toward zero. On a field simulator with the
    truth known exactly, where the effect is largest that cost an rmse of 0.091 against 0.074
    and a bias of -0.060 against -0.042; at the focus itself, 0.352 for a true 0.500 against
    0.444.

    But the two assert over **different extents**, and getting that wrong is far worse than
    dropping the limb altogether. A silence is a statement about a neighbourhood: nothing
    within ``coverage_radius`` cleared this study's cut. A report is a statement about one
    voxel, because a reported peak is a local maximum selected for being large and displaced
    from wherever the effect is -- so "someone reported 18 mm away" is not evidence that the
    effect *here* cleared anything. Asserting the report across the same sphere as the silence
    gave an rmse of 0.457 against 0.114 for the named voxel alone, and turned a -0.042 bias
    where the truth is largest into +0.097. So a report asserts its indicator at the voxel it
    names; a voxel a study reached but did not name carries no indicator at all.

    That asymmetry is invisible without a spatial dimension. A one-voxel likelihood, with no
    displacement and no radius, says "restore the report limb" unconditionally -- 0.351 to
    0.524 for a true 0.500 -- and cannot see the radius problem at any signal level.

    **The floor is the threshold, not zero.** Because :math:`\mu` is bounded below the cut only
    as far as the silences push it, a quiet region does not read as an effect of zero; it reads
    as an effect the reporting studies could not detect. The pull toward zero comes through
    :math:`\pi` instead -- many studies, all silent, means few of them have an effect -- which
    is why ``g_marginal`` is the better map where nothing was reported (bias +0.056 against
    +0.095 for the images alone) and the worse one where an effect exists (-0.096).

    The interval
    ~~~~~~~~~~~~
    Refer ``se`` to a *t*, not to a normal, and **mask on coverage before you do**. ``dof`` is
    the number of studies on the censoring roster minus one, zeroed where no study speaks about
    the voxel. The roster rather than a Kish count over the pooling weights, for two reasons,
    and the second is decisive: the likelihood uses every study on the roster -- the images
    through their values, the rest through their silence -- so crediting only the weighted
    contributors denies the model information it demonstrably used; and with the coordinate
    magnitudes gone the only weighted contributions are images at weight 1, so a Kish count is
    just the image count, making ``dof`` **zero with one image** and the recommended interval
    ``nan``. Whether a censored-likelihood observed information really carries the roster's
    degrees of freedom is a genuine question rather than a settled one; a profile-likelihood
    interval would need no ``dof`` at all.

    Under the selection model ``se`` is the observed information of the censored mixture
    likelihood at the fitted point, with the prevalence profiled out by a Schur complement, so
    it carries both the uncertainty about which component an observation came from and the cost
    of not knowing the prevalence. On the estimator's own censored mixture with known variances
    it covers 94.5% to 98.4% of nominal-95% intervals across prevalences, cutoffs and study
    counts, erring conservative.

    **End to end the interval is conservative rather than wrong, and coverage will not tell you
    that either way.** Measured on the field simulator with the truth known exactly, 40
    replications per arm, stratified by the truth because 9204 of 9261 voxels sit near zero and
    a whole-map figure rewards any estimator that shrinks. ``width`` is the half-width of the
    documented interval as a fraction of the truth, so the quiet stratum's is meaningless by
    construction and omitted:

    =========================  ======  =========  ======  ========  =========  ========
    arm                          bias  ``se/sd``  cov(t)  bias      ``se/sd``  width
                                       (quiet)    (quiet) (effect)  (effect)   (effect)
    =========================  ======  =========  ======  ========  =========  ========
    20 studies, 1 image        -0.000       1.51    0.98    -0.089       1.39      0.99
    20 studies, 2 images       +0.001       1.80    0.98    -0.086       1.66      0.89
    20 studies, 5 images       -0.000       2.06    0.99    -0.038       1.83      0.59
    20 studies, 20 images      +0.000       1.16    0.98    -0.023       1.19      0.26
    2 images, silence off      +0.001       1.23    1.00    -0.040       1.22      5.58
    2 images, tau = 0.3        +0.001       1.80    0.98    -0.097       1.78      1.54
    =========================  ======  =========  ======  ========  =========  ========

    **An earlier version of this table was wrong in a way worth recording, because the error
    was in the measurement and not the model.** It reported a quiet-stratum bias of +0.045 to
    +0.114 and ``se/sd`` of 2.05 to 3.74. Both came from taking ``numpy.abs`` of ``g`` before
    the spread across replications was computed. ``g`` is a *signed* inverse-variance mean, so
    at a voxel whose truth is zero the absolute value shrinks the spread to about 0.6 of the
    real one and puts the mean about 0.8 spreads above zero: it halved the denominator of
    ``se/sd`` and manufactured the bias out of nothing, at once. At the foci, where the estimate
    sits far from zero, the absolute value is nearly a no-op -- which is exactly why the old
    table's *quiet* column ran 2.71 to 3.74 while its *effect* column ran 1.55 to 1.97, and
    nobody looked twice. **The estimator is not biased upward where nothing was reported.**

    The defect was confined to this: a spread of :math:`|g|` read as the estimator's sampling
    spread, or a difference of :math:`|g|` from a near-zero truth read as a bias. Comparisons
    between arms that take the absolute value of the *reference* as well -- the collection
    comparisons below, and every ranking against SDM-PSI or an images-only pool -- apply the
    same transformation to every arm and are unaffected.

    The check that would have caught it on the first run, and is now the top row of the arm
    list: give every study an image and switch the selection model off, and the fit is a
    textbook local inverse-variance random-effects meta-analysis whose ``se/sd`` must come out
    near 1. It reads 1.00 to 1.24. A little of the residue above 1 is the bed rather than the
    estimator -- the simulator's per-study noise measures 4% below the ``1/n + g^2/(2n)`` the
    model assumes, so every ``se`` is 4% generous by construction.

    Two things follow.

    **The inflation was the prevalence, not the censoring term.** The old reading -- "the excess
    is in the censoring term, since switching the silence off halves it" -- compared two
    different estimators and was never a localisation. The configuration that settles it is the
    all-image one: with every study carrying an image the reporting indicator is *structurally*
    empty, ``coordinate_share`` is identically zero, and the coordinate channel cannot be
    responsible for anything. That arm was nonetheless the worst measured, ``se/sd`` of 2.21 at
    quiet voxels against 1.14 for the same images fitted without the mixture, in the same bed on
    the same seeds. What it was paying for was a prevalence that nothing identified: only the
    indicator separates "no effect in this study" from "a small effect plus noise", and left
    free against 20 Gaussian values the mixture returned :math:`\pi = 0.577` against a true
    1.000, whose uncertainty was then profiled out of the information about :math:`\mu`. The
    inflation was concentrated at quiet and weak voxels (2.21 and 2.30) and absent at strong
    ones (0.99 and 1.06), which is the signature of that trade-off: where the effect is weak,
    "a small effect in every study" and "a large effect in a few" fit equally well.
    :math:`\pi` is now held at 1 wherever no study contributes an indicator, and the all-image
    arm returns the non-mixture fit exactly -- 1.16 in the table above, and asserted as an
    identity in the tests rather than as a tolerance.

    **What remains is a conservative interval where the coordinates do act, 1.28 to 2.03.** It
    is anomalous in one direction only: the likelihood conditions on where the foci fell while
    the replication spread is marginal over that, so a calibrated *conditional* error should sit
    *below* the marginal spread, not above it. Two candidate causes are now ruled out. It is not
    the reporting rule: the model censors on :math:`|g| < c` while a paper reports local maxima,
    and feeding the estimator its own rule instead -- every supra-threshold voxel reported, 386
    foci per collection against 172 -- moves ``se/sd`` from 1.83 to 1.80 at quiet voxels and
    1.90 to 1.90 at the strongest focus, and the bias not at all. It is not the baseline either,
    at 1.00 to 1.24. The residue tracks how much the indicator says about :math:`\pi`:
    ``se/sd`` falls monotonically with ``coordinate_share``, 1.87 at under 0.05 to 1.33 above
    0.50. So the same channel that identifies the prevalence is what makes its cost bearable,
    and where the tables are thin the interval on :math:`\mu` is wider than the estimate
    deserves.

    **Changing the estimand does not fix it, and this was tested rather than assumed.** The
    ridge runs along curves of roughly constant :math:`\pi\mu`, so the product ought to be the
    determined combination and ``g_marginal`` the sound thing to bracket. Half of that is true
    and the useful half is not. The product really is the steadier *estimate* -- its spread
    across replications is 0.076 at quiet voxels against 0.112 for ``g`` -- but its reported
    error is worse, 0.349 against 0.204, and ``se_marginal/sd`` comes out above ``se/sd`` in
    every stratum but the strongest: 4.58 against 1.83 at quiet voxels, then 2.04, 2.48, 2.43
    and 1.58. The reason is mechanical: :math:`\operatorname{Var}(\mu\pi)` needs the whole
    2x2 inverse rather than a Schur complement, and a near-singular information matrix amplifies
    there instead of cancelling. **So ``g_marginal`` is the more stable estimate and the less
    trustworthy interval.**

    **And the width is not the problem. There is no finite width to be had.** The obvious
    remedy for both -- a profile likelihood, which inverts nothing and needs no ``dof`` -- was
    implemented (``interval="profile"``) and it answers the question in a way that disqualifies
    the question. As :math:`\pi \to 0` the active component explains nothing, the mixture
    density tends to the null one at every observation whose probability depends on
    :math:`\mu`, and the profile log-likelihood approaches a *horizontal asymptote*,

    .. math::

        \ell_\infty = \max_\pi \left[
            \sum_{i \,\notin\, \mathrm{reported}} \log\left((1 - \pi) b_i\right)
            + \sum_{i \,\in\, \mathrm{reported}} \log\left(\pi + (1 - \pi) b_i\right)
        \right],

    where :math:`b_i` is observation :math:`i`'s density or probability under an effect of
    exactly zero. Image values and silences drop their :math:`\mu` dependence because both
    probabilities vanish as :math:`|\mu| \to \infty`; a *report* does not, its probability
    tending to 1 instead, which is the only thing keeping :math:`\pi` in the expression. So the
    interval is bounded exactly when :math:`2(\hat\ell - \ell_\infty)` exceeds the critical
    value -- and at a voxel where no study reported, where the maximum is attained as
    :math:`\pi \to 0`, that is precisely the likelihood-ratio test of :math:`\pi = 0`. With a
    handful of image studies it almost never fires. Measured on the field bed, the fraction of
    voxels where the interval is bounded at all:

    ==============  ================  =============  =============
    image studies   bounded overall   within 10 mm   beyond 30 mm
    ==============  ================  =============  =============
    2                         0.025          0.086          0.029
    5                         0.055          0.148          0.065
    10                        0.072          0.284          0.083
    ==============  ================  =============  =============

    Reparametrising does not escape it: :math:`\pi\mu` has the same asymptote, reached along
    :math:`\pi \to 0` with :math:`\mu \to \infty`.

    **What does escape it is a literature whose studies differ in size, and the reason is
    algebraic.** A silence constrains :math:`(\pi, \mu)` only through the probability of the
    event observed, :math:`P = \pi S(\mu) + (1 - \pi) S(0)` with
    :math:`S(\mu) = P(|g| < c \mid \mu)` -- one equation in two unknowns, which is the ridge.
    Studies with *different* :math:`(\sigma, c)` supply different equations, but not all
    heterogeneity helps, because the reporting cutoff measured in sampling standard deviations
    is just the reported statistic back again: :math:`c / \sigma \approx z`. So

    .. math::

        S(0) = 2\Phi(z) - 1, \qquad
        S(\mu) = \Phi(z - \mu\sqrt{n}) - \Phi(-z - \mu\sqrt{n}).

    The null component's silence probability depends on the *threshold alone*. Varying the
    threshold moves both components together through the same tail and leaves the equations
    nearly collinear; varying the sample size moves the active component through
    :math:`\mu\sqrt{n}` while leaving :math:`S(0)` exactly fixed, which is the contrast that
    separates the two parameters. Measured, 20 studies throughout, with the bounded fraction of
    the profile interval near the signal as the identifiability probe:

    ==========================  ==============  ===================
    studies                     bounded near    fitted :math:`\pi`
    ==========================  ==============  ===================
    alike, n 28-32, one cut              0.422                0.823
    sample size spread 12-120            0.691                0.913
    threshold spread 2.3-4.5             0.414                0.815
    both spread                          0.676                0.901
    ==========================  ==============  ===================

    True :math:`\pi` is 1.0. Spreading the sample size identifies it; spreading the threshold
    does nothing whatever, exactly as the cancellation above says. **So a magnitude is
    recoverable from a literature of widely differing sample sizes and not from a literature of
    uniformly sized studies, however many of them there are** -- which is the opposite of
    treating heterogeneity as a nuisance, and is the one piece of advice here that bears on
    whether to run this estimator on a given collection at all.

    **Necessary, not sufficient, and the boundary was measured rather than assumed.** On the HCP
    held-out bed -- a strong motor contrast, about fourteen peaks per table -- the interval is
    already bounded at 97% of the scored voxels whether the table sizes are uniform or spread
    from 9 to 70, and spreading them changes the recovered magnitude not at all (0.680 against
    0.671, with ``r`` and AUC identical to three decimals). The derivation predicts exactly
    that: where :math:`\pi = 0` is decisively rejected there is no ridge for heterogeneity to
    break. What it also shows is that the shortfall there is *not* a flat likelihood -- with
    :math:`\mu` identified almost everywhere, ``g`` still recovers 0.68 of the truth against
    the images' 0.92, and :math:`\pi` still fits 0.66 against a true 1.000. Those are
    properties of what the censoring term assumes a silence to mean, and they are the subject
    of the warnings below.

    **So ``se`` reports a curvature at the point the EM selected, and outside that regime the
    likelihood does not support that precision about :math:`\mu`.** The estimate's stability
    across replications -- a spread of 0.11 where ``se`` says 0.20 -- comes from the pooled
    image mean it starts from and the ``max_iter`` it stops at, not from the data pinning it
    down. That is worth
    knowing before reading ``g`` as a magnitude, and it is the strongest statement available
    about why the ``se/sd`` question never resolved: it was asking whether an interval was
    calibrated for a parameter the data do not bound.

    None of this touches the p-values, which come from the permutation null and need no
    calibrated ``se``, nor the *ordering* of ``g``, which every collection comparison here
    scores and which holds up. It bears on reading a single voxel's ``g`` as a number with an
    error bar.

    **The ``silence off`` row's width is not comparable.** Without the selection model there is
    no censoring roster, so ``dof`` falls back to the Kish count over the image weights, which
    at two images is 1 -- and a *t* on one degree of freedom has a critical value of 12.71. That
    row is a correct statement about two studies, not a wider interval for the same
    information.

    **Read ``se/sd`` and the width, never coverage.** Every arm covers 0.95 to 1.00, including
    the ones whose interval admits almost any magnitude: at two images the half-width is 0.89 of
    the effect, and coverage alone cannot distinguish that from the all-image row's 0.26.

    **P-values are unaffected by any of this.** They come from the permutation null, which is
    valid for whatever statistic it is computed on and does not require a calibrated ``se``.

    An earlier version reported the curvature of the EM's *Q function* instead of the observed
    information, which holds the responsibilities fixed and so overstates the information; that
    covered 62.5% to 89.8% and did not improve with more studies.

    :meth:`correct_fwe_montecarlo` adds ``logp_level-voxel``,
    ``logp_desc-size_level-cluster`` and ``logp_desc-mass_level-cluster`` (each with a signed
    ``z_*`` companion), matching the names :class:`~nimare.meta.cbma.ale.ALE` uses.
    :class:`~nimare.correct.FDRCorrector` and ``FWECorrector(method="bonferroni")`` work off
    the uncorrected ``"p"`` map instead, and are only meaningful when that map came from the
    permutation null.

    Warnings
    --------
    This estimator is new and has not been validated against a reference implementation.

    **The correction can make ``g`` worse than doing nothing, and the regime where it does is
    not exotic.** Judged against a reference built from subjects used to make no coordinate at
    all -- 786 HCP subjects, 480 cut into 16 synthetic studies of 30 with coordinates extracted
    the way papers produce them, 306 held out to give the truth:

    ======================  ======  ========  =====  ===========
    estimate                     r  rank r     AUC   magnitude
    ======================  ======  ========  =====  ===========
    two images, pooled      +0.845    +0.576  0.973         0.85
    ``g``                   +0.830    +0.575  0.967         0.63
    ``g_marginal``          +0.785    +0.564  0.943         0.54
    ======================  ======  ========  =====  ===========

    Pooling the two images alone wins on every metric, magnitude included -- **and so does an
    independent method given the same data.** SDM-PSI, which takes the same mixture by design,
    was run on the same studies with the same two supplied as maps: it returns r +0.722, AUC
    0.927 and 0.43 of the reference, behind ``g`` on every column and behind the two images
    alone on every column. Two unrelated methods both do worse with the fourteen coordinate
    tables than without them, which makes this a property of thresholded coordinate tables in
    this regime rather than a quirk of this estimator's censoring term. **That design has a
    prevalence of exactly 1**: every synthetic study draws from the same population, so there is
    no between-study absence for the mixture to find, and an unthresholded map of 30 of those
    subjects is already nearly unbiased. There is nothing for a selection correction to correct,
    and it does harm anyway -- fitted prevalence comes back at 0.664 against a true 1.0 (0.929 at
    the strongest decile), because "failed to clear its threshold" and "has no effect" both
    explain a silence and the mixture splits the difference.

    **But the prevalence is the symptom, not the cause.** Refitting the same collection with the
    prevalence *fixed* at 1 -- which is the truth here -- makes the magnitude **worse**, not
    better: 0.60 of the reference against 0.63, and 0.411 against 0.422 for a true 0.5 on the
    simulator. Removing the "this study has no effect" escape forces every silence to be
    explained by a small :math:`\mu`, so :math:`\mu` falls further. So the **censoring term
    over-shrinks :math:`\mu` whatever the prevalence does**, and a fitted :math:`\pi` below 1 is
    the model partly *absorbing* that over-shrinkage rather than adding to it. ``g_marginal``
    comes back worst of all (0.54) because it multiplies the two together.

    Which makes this the same defect as the overstated reporting probability above, from another
    direction: too high a :math:`P(\text{report})` against an under-observed count drags
    :math:`\mu` down, and here there is no genuine absence for the prevalence to absorb it
    into. Pinning the prevalence is therefore not the middle option it sounds like.

    Contrast the 21-study NIDM pain collection above, where the same correction cuts rmse 23%
    and bias 47%. **Why the two disagree is not settled, and the obvious explanations have been
    tested and rejected.**

    It is not the prevalence. Building the regime as a dial -- MOTOR_LH studies that carry the
    effect mixed with EMOTION_FACES studies that do not, so the true :math:`\pi` is designed
    rather than assumed -- ``g`` recovers 0.68, 0.63, 0.56 and 0.53 of :math:`\mu` at true
    prevalences of 1.00, 0.75, 0.50 and 0.25. It degrades *monotonically* as prevalence falls,
    rather than improving.

    Nor is it the statistic. On the pain bed's own whole-map rmse the two arms tie on that dial
    (0.130 against 0.129 at :math:`\pi = 1`) or the images win (0.281 against 0.241 at 0.50).

    Nor is it between-study heterogeneity, the other obvious candidate: adding a relative
    :math:`\tau` of 0.0, 0.3 and 0.6 to the dial's effect studies leaves rmse at 0.126, 0.135
    and 0.159 against the images' 0.126, 0.123 and 0.140 -- no crossing anywhere.

    **What does move it is how the reference is built, and that is a caution about the pain
    number rather than an endorsement of it.** The pain reference was an inverse-variance mean
    of 19 study-level ``g`` maps. Hedges' variance is a function of the *observed* effect, so a
    study that drew high gets less weight and such a reference is itself pulled downward --
    and ``g`` is pulled downward too, so an estimator biased low scores better against a
    reference biased low. Rebuilding the dial's reference the same way, from synthetic
    reference studies rather than pooled subjects, takes ``g`` from tied (0.126 against 0.126)
    to winning (0.122 against 0.125). **But that accounts for about 2 points of the pain gap's
    22**, so the direction of the artefact is demonstrated and its magnitude is not. The one
    remaining difference -- real studies against synthetic, with their different scanners,
    paradigms and sample sizes -- cannot be dialled.

    So the advice here is empirical rather than principled: **fit it both ways.**
    ``selection_model="none"`` reduces to an inverse-variance meta-analysis of the images and
    costs about a tenth of the runtime, so the comparison is cheap. Where the two agree, little
    turns on the choice; where they disagree sharply, the correction is doing something
    load-bearing that nothing measured here can yet vouch for.

    **``coordinate_share`` is the diagnostic every caveat here needs, and it is the one thing a
    reader could not otherwise get.** The warnings above are all about *when* to trust the
    coordinate channel -- the prevalence-1 regime where it does harm, the unexplained
    disagreement between collections, the over-shrinkage at the window. None of them can be
    checked on a given collection. What can be checked is whether the channel is even acting at
    a voxel: ``coordinate_share`` is the fraction of the information about ``g`` contributed by
    the indicators rather than the images' values, and it comes free because the two are
    accumulated separately inside the likelihood. At 0 the images carry the estimate alone and
    the tables changed nothing here, so none of the coordinate caveats apply; at 1 the
    indicators carry it and all of them do.

    On a 20-study collection with two image donors it runs from 0.02 to 1.00, with a **median of
    0.10 and 0.79 at the focus** -- so on a typical map the images carry the estimate almost
    everywhere and the coordinates take over exactly where studies reported. That is the
    stratification the design rests on, now readable per voxel rather than only in aggregate.

    **It also says where the interval is trustworthy, in the direction opposite to the obvious
    guess.** If the ``se``'s conservatism were produced by the censoring term, ``se/sd`` would
    be worst where the share is high. It is the reverse, and monotonically so: across bands of
    the share from below 0.05 to above 0.50 it runs 1.87, 1.82, 1.87, 1.63 and 1.33, with the
    top decile at 1.54 against the bottom's 1.86. That is not a signal effect masquerading as a
    share effect -- only a few hundred of these 15625 voxels carry any truth, so the band from
    0.25 to 0.50 is almost entirely quiet, and the gradient holds inside it. The reason is the
    one the interval section gives: what the ``se`` is paying for is an imprecisely known
    prevalence, the indicator is the only thing that identifies it, and the share is precisely
    how much indicator reached this voxel. So a high share is the regime where the coordinate
    caveats bite *and* where the interval is soundest, and the two readings do not conflict.

    **Read ``prevalence`` ordinally, not as a fraction, and not within one map.** On the
    designed-prevalence dial just described -- real subjects, a true :math:`\pi` set by how many
    studies carry the effect -- it reads 0.918, 0.714, 0.590 and 0.540 against true values of
    1.00, 0.75, 0.50 and 0.25, measured at the strongest decile. So it tracks well down to about
    0.5 and then **floors near 0.54**, which is the compression in its least flattering place:
    a rare effect and a common one come back nearly the same. Against a simulator the same
    compression is worse still -- a true 0.25 comes back as 0.49 to 0.60 depending on
    ``coverage_radius``, a true 0.50 as 0.65 to 0.81, a true 1.00 as 0.74 to 0.94.
    Its map-wide median sits near 0.4 whatever the truth,
    so a map cannot be summarised by it.

    Worse, the ordinal reading holds on average over many maps and rarely within any one of
    them. Tested as the claim is made -- four sites in a single fit at true prevalences 0.25,
    0.50, 0.75 and 1.00 -- the rank correlation against the truth averages +0.76 for a strong
    effect, but the four-site ranking is exactly right in only **19%** of maps; for a weak
    effect it averages +0.33 with a standard deviation of 0.60 and is exactly right in **6%**.
    The per-voxel scatter is largest at the low end, so rare sites are both biased upward and
    noisier -- the worst combination for the use this invites, picking out which region is the
    least consistent.

    That is structural rather than a calibration that could be fixed. :math:`\pi` and
    :math:`\mu` are separably estimable only in a window of detectability: where a study's
    effect lands near its own reporting threshold, so that the chance of reporting responds to
    the magnitude. Below that window nothing is detected and the prevalence is not identified
    at all; above it detection saturates, the magnitude stops being constrained from above and
    the prevalence absorbs the level instead -- which is why a strongly reported site returns a
    prevalence near 1 whatever its truth. A map spans magnitudes and therefore spans the
    window, so comparing two voxels compares quantities identified to different degrees. A
    spread of sample sizes and reporting thresholds across the collection widens the window; a
    roster of identically powered studies narrows it.

    **``g`` and ``g_marginal`` are different estimands, and the second is the one an image-based
    meta-analysis reports.** Every IBMA -- DerSimonian-Laird, Hedges, weighted least squares,
    the likelihood estimators -- pools per-study effect maps around a single mean, so a study
    with no effect at a voxel enters that average as a zero and the quantity estimated is
    :math:`\pi(v)\,\mu(v)`. ``g`` is :math:`\mu(v)`, the effect over the studies that have one.
    The two differ by a factor of :math:`1/\pi`, which on the NIDM pain collection is 1.4 to
    1.8. So ``g`` cannot be checked against an image-based reference even in principle;
    ``g_marginal`` can, and is the map to compare. Conversely :math:`\mu(v)` may not be
    identifiable from images at all: computing it requires classifying every study as having an
    effect at every voxel or not, which is a thresholding decision and reintroduces the
    selection this estimator exists to correct.

    **``g_marginal`` does not work for the reason its name gives.** In a held-out-subject design
    where every synthetic study is drawn from one population, the true prevalence is exactly 1
    and ``g_marginal`` should equal ``g``; instead ``prevalence`` comes back near 0.68 and the
    product is the better estimate. Multiplying by it is shrinking a magnitude by a data-driven
    factor, not averaging over studies that have no effect. That the two errors cancel is why it
    is worth reporting and also why it should not be trusted outside the regimes it has been
    measured in. It degrades when studies report few foci, because ``prevalence`` then falls
    toward its floor: at six foci per study it came back at 0.70 times the held-out truth.

    **The dynamic range is recovered, which the earlier design's was not.** Compression was that
    design's headline failure: a reference effect spanning elevenfold across its strata came
    back spanning about 1.2-fold, and an unknown overall scale would have left that ratio alone,
    so it was a real defect and not a units problem. Re-measured here on four well-separated
    foci at true ``g`` of 0.2, 0.4, 0.6 and 0.8 inside one map, 8 collections, regressing the
    estimate on the truth over the voxels carrying signal:

    Values are the estimate at each focus's own voxel, with its standard error across the 24
    collections and its relative error:

    ==============  =====  =========  ==================  ==================
    estimate        slope  intercept  @ 0.2               @ 0.4
    ==============  =====  =========  ==================  ==================
    images only     0.911     +0.043  0.230 +- .025 (+15%)  0.419 +- .029 (+5%)
    ``g``           0.759     +0.049  0.194 +- .019 (-3%)   0.329 +- .022 (-18%)
    ``g_marginal``  0.769     +0.018  0.151 +- .018 (-25%)  0.267 +- .024 (-33%)
    ==============  =====  =========  ==================  ==================

    ==============  ==================  ==================  =========  ==========
    estimate        @ 0.6               @ 0.8               range      mean |err|
    ==============  ==================  ==================  =========  ==========
    images only     0.635 +- .027 (+6%)   0.792 +- .025 (-1%)  3.44-fold        6.7%
    ``g``           0.634 +- .016 (+6%)   0.825 +- .017 (+3%)  4.26-fold        7.5%
    ``g_marginal``  0.576 +- .015 (-4%)   0.780 +- .019 (-3%)  5.18-fold       16.1%
    ==============  ==================  ==================  =========  ==========

    **The range is recovered and the level is not uniformly better.** ``g`` returns 4.26-fold
    for a true 4-fold against 3.44-fold for pooling the images alone, so the compression that
    was the earlier design's headline failure is gone. But its *mean* absolute relative error
    over the four foci is 7.5% against the images' 6.7% -- a wash, slightly the wrong way --
    because the two arms trade errors focus by focus rather than one dominating.

    Where they trade is informative, and it maps onto the window of detectability. At the
    weakest focus, where **no** study reports (0 of 18 in a representative collection), the
    images are inflated 15% and ``g`` is 3% off: the correction is working, though note that the
    coordinate channel contributes nothing usable there -- with nothing reported and silence
    nearly flat in :math:`\mu` below the cut, ``g`` is the two-image estimate and its accuracy
    is the images' sampling noise. At the strongest focus, where 12 of 18 report, both are
    within a few percent.

    **The 0.4 focus is the clear remaining defect: 18% low, about three standard errors, where
    the images are 5% high. It is located, and it is the reporting model's functional form.**
    The censored likelihood gives a reported voxel the probability :math:`P(|g| \ge c)`, but a
    paper reports a voxel only if it cleared :math:`c` **and** was a local maximum -- a strictly
    smaller event. Instrumenting the same bed to compare the observed reporting rate against the
    rate the model computes at the true :math:`\mu`:

    =======  ========  ==============  =============  =====================
    truth    reports   observed rate   model's rate   ratio (95% CI)
    =======  ========  ==============  =============  =====================
    0.2             3           0.007          0.007  0.97 (unmeasured)
    0.4            21           0.050          0.086  **1.72 [1.19, 3.03]**
    0.6           163           0.392          0.392  1.00 [0.87, 1.18]
    0.8           299           0.702          0.789  1.12 [1.01, 1.27]
    =======  ========  ==============  =============  =====================

    Counts are totals over 24 collections and the interval is Poisson on them, because this is a
    ratio of two small rates and a point estimate would not be a measurement. **The
    over-statement is real in the middle of the window** -- the 0.4 interval excludes 1 -- and
    it is what pulls :math:`\mu` down there, which is the whole of that focus's deficit. Above
    the window it is absent or mild, and the two intervals overlap, so the shape is "well below
    1 at the window, at or just past 1 above it" rather than a peak with two sides. At 0.2
    nothing is reported, so nothing is measured.

    **The mechanism is confirmed by intervention, and the shape of the fix is known.** The
    correction is not a reweighting -- that was tried twice and made every focus worse. It is
    that both limbs must be complementary probabilities of the same event: writing
    :math:`E = P(|g| \ge c \mid \mu)` and :math:`q` for the chance that a study exceeding here
    actually *names* this voxel, a reported pair carries :math:`qE` and a silent one
    :math:`1 - qE`. That is still a proper likelihood, and the arithmetic says where it acts --
    the report limb's score is :math:`(qE)'/(qE) = E'/E`, so :math:`q` cancels there, while the
    silent limb's becomes :math:`-qE'/(1 - qE)`, weakened by roughly :math:`q`. Exactly the
    over-shrinkage above.

    Measured with :math:`q` fixed at 0.56, the reciprocal of the 1.78 at the 0.4 focus, that
    focus closes precisely: -11% to +1%. **But no constant :math:`q` helps overall** -- mean
    absolute error over the four foci runs 8.9%, 9.4%, 10.3%, 11.0% and 12.2% at
    :math:`q` of 1.00, 0.90, 0.80, 0.70 and 0.56, because a scalar lifts the weak foci and
    overshoots the strong ones. The shipped :math:`q = 1` is the best constant.

    So the fix needs :math:`q(\mu)`, the probability that an exceeding voxel is the one actually
    named. **A random-field expected-maxima density is not that function, and gets its sign
    backwards.** The usual clump argument gives :math:`q \sim 1/\text{clump size} \sim u^3`
    for a standardised threshold :math:`u = (c - \mu)/\sigma`, which *falls* as :math:`\mu`
    rises -- :math:`u` runs +1.10, 0.00, -1.10 across the three foci above -- while the measured
    :math:`q` *rises* (0.58, 1.00, 0.89). That density describes a zero-mean field, and these are
    signal peaks: at a strong focus the blob's own curvature makes that voxel the local maximum,
    so :math:`q` approaches 1, while at a marginal focus the noise decides which of several
    exceeding voxels is the maximum. The governing quantity is the signal's curvature against the
    noise smoothness, which a coordinate table does not carry and a reported FWHM does not
    supply.

    **That is now derived rather than inferred, which bounds it.** Write the observed field as
    :math:`Z = m + e` with :math:`e` smooth, stationary and mean-zero. A local maximum needs
    :math:`Z'(0) = 0` and :math:`Z''(0) < 0`. At a signal peak the first condition reduces to
    the null one, because :math:`m'(0) = 0` there; the second does not, since
    :math:`Z''(0) = -\kappa + e''(0)` with :math:`\kappa = -m''(0) > 0`, so

    .. math::

        P\big(Z''(0) < 0\big) = \Phi(\kappa / \sigma_2), \qquad
        \sigma_2 = \operatorname{sd}\big(e''(0)\big).

    A zero-mean density is the case :math:`\kappa = 0`, where this is exactly one half, and
    :math:`\Phi(\kappa/\sigma_2)` increases strictly in :math:`\kappa`. So the standard
    density is not merely wrong at a signal peak, it is a **lower bound**, understating the
    chance of a maximum by a factor of :math:`2\Phi(\kappa/\sigma_2)` -- one where the mean is
    flat, tending to two as the peak sharpens. Exactly one of the two conditions defining a
    maximum carries the signal, and it is the one that density pins at a half.

    Which settles the status of the limb rather than only its sign: the error is understood, it
    is bounded by two, and **it is not fixable from tables alone**, because the correction needs
    a per-study peak sharpness no paper reports. Open, and now open for a stated reason.

    Worth stating alongside, because it is easy to assume otherwise: **the indicator channel does
    not dominate the fit.** At the 0.4 focus the observed count alone implies
    :math:`\mu = 0.161` and the images imply about 0.4; the fit lands at 0.337. Below the window
    the count implies a nonsensical :math:`\mu = -0.144` and the fit is within 6% of the truth.
    The images carry the magnitude; the indicators move it.

    The slope of 0.759 is not the foci: it comes from the blob skirts, where the truth runs 0.05
    to 0.2 and both arms are dominated by the floor that reading a map as ``|g|`` imposes. Read
    the per-focus columns, not the slope.

    ``g_marginal`` should be read narrowly: within 4% at the two strong foci and 25% to 33% low
    at the two weak ones, because ``prevalence`` falls toward its floor exactly where few studies
    reported. Prefer it to ``g`` only when comparing against an image-based reference, where it
    is the matching estimand.

    A corollary for reading any three-bin summary of this estimator, including the ones above
    under "The interval": a top bin spanning 0.25 to 0.50 of truth averages voxels whose
    estimate is slightly high with voxels whose estimate is low, and reports the mixture as a
    bias. The slope, the intercept and the per-focus values are the honest summary.

    **What the null tests is not what a reader may expect.** The null is that *within a study,
    effect size is unrelated to location*. A voxel is significant when the image studies'
    effects near it are large relative to what the same studies show elsewhere -- not when the
    pooled effect differs from zero, and not when studies converge there. A collection with a
    genuine effect of the same size everywhere has nothing for this null to find. The
    zero-effect null is not available: a reported peak exists only because it cleared a
    threshold, so "no effect anywhere" predicts no coordinates at all and the observed table
    falsifies it before any voxel is examined; testing it needs subject-level images, which is
    what :footcite:t:`albajes2019meta` imputes in order to permute. Sign-flipping the image
    studies while shuffling the coordinates would test it for part of the collection only, and
    the two hypotheses then combine into a rejection either can cause -- with 20 coordinate and
    5 image analyses the sign flips alone floored the p-value at 1/32 whatever the locations
    said.

    **The collection must supply an image, and one image is a thin basis for a magnitude.**
    ``g`` at a voxel reached by a single image is that image's value corrected by the others'
    silence, with no between-study spread to estimate and no replication to average over. That
    configuration is supported, and measured better than two images on real collections, but the
    magnitude it reports rests on one study's map.



---

## The parameter justifications, as the docstring carried them

Condensed in the code to what changes a user's choice. The measurements behind each default are reproduced here verbatim.

```
    Parameters
    ----------
    design : {"one-sample", "two-sample"}, default="one-sample"
        Design behind the collection, used for the sampling variance a silent study is judged
        against.
    tau2_method : {"dl", "none"}, default="dl"
        ``"dl"`` estimates a local between-study variance with a DerSimonian-Laird moment
        estimator over the image studies; ``"none"`` fits a fixed-effects model
        (:math:`\tau^2 \equiv 0`).

        ``"dl"`` is estimated once, about the naive weighted mean, and then held fixed while the
        selection model fits :math:`\mu` -- which keeps each EM iteration one-dimensional and
        concave, at a known cost. Because the weighted mean is by construction the centre that
        *minimises* the moment estimator's ``Q``, taking ``Q`` about the value finally reported
        can only raise :math:`\tau^2`, and the shipped estimate is therefore biased low.
        Measured against a known :math:`\tau = 0.35`, alternating the two recovers
        :math:`\tau^2` of 0.067, 0.083, 0.093, 0.100 over three extra rounds against a true
        0.1225, and moves ``g`` at the focus from 0.883 to 0.822 against a true 0.8. No spurious
        heterogeneity appears where there is none. Not done, because each round is a full refit
        and a principled joint estimate is a larger change than alternation.
    selection_model : {"zero-inflated", "none"}, default="zero-inflated"
        Whether the silence of the coordinate studies is read at all.

        ``"zero-inflated"``
            Each coordinate study's reporting indicator enters a zero-inflated censored
            (Tobit) likelihood, under a mixture in which the study either has a real effect or
            none at all: the probability of staying silent where it has no focus within
            ``coverage_radius``, and the probability of clearing its cut at a voxel it named.
            This is the entire reason the coordinates are in the model. Nothing is imputed.
        ``"none"``
            Silence is not read, so the fit reduces to an inverse-variance random-effects
            meta-analysis of the images alone and the coordinate tables have no effect
            whatever. Useful as the control arm -- it is what the coordinates are being
            credited against -- and roughly ten times faster.
    se_method : {"model", "hksj"}, default="model"
        Standard error of the pooled estimate. ``"model"`` is the inverse-variance expression,
        which treats the estimated :math:`\tau^2` as known; ``"hksj"`` is the
        Hartung-Knapp-Sidik-Jonkman residual-variance form, which does not and covers better
        with few studies. **Requires** ``selection_model="none"``, the zero-inflated model
        reporting the censored likelihood's curvature instead. Changes ``se`` and so ``z``;
        p-values come from the permutation null either way.
    interval : {"wald", "profile"}, default="wald"
        How to bracket ``g``. ``"wald"`` reports ``se`` alone, referred to a *t* on ``dof``.
        ``"profile"`` additionally emits ``g_lower`` and ``g_upper`` from the profile
        likelihood -- the set of :math:`\mu` whose log-likelihood, with the prevalence
        maximised out at each point, sits within half a chi-square critical value of the
        maximum. It inverts no matrix and needs no ``dof``, which is why it is worth having
        here specifically: both of the quantities it replaces degrade as the prevalence becomes
        weakly identified, the Schur complement and the delta method alike. The bounds are
        asymmetric, as a likelihood region generally is, so they are reported rather than
        summarised as a half-width. Costs roughly a second fit, and **requires**
        ``selection_model="zero-inflated"``. Provisional: ``se`` is still what ``z`` is built
        from, and the interval's calibration against the arm table above is not yet measured.
    analysis_mask : :obj:`str` or None, optional
        ``value_type`` of a per-study image marking the voxels that study examined, nonzero
        meaning examined. Studies without one are taken to have examined the whole analysis
        volume, which is what the censoring term assumes of every study by default.

        This is what an ROI or partial-coverage study needs: its silence outside the region it
        analysed is not evidence that nothing is there, and the censoring term would otherwise
        read it as evidence against an effect. Voxels a study did not examine contribute
        neither a value nor a silence for it. Without this the only remedy was
        ``selection_model="none"``, which discards the whole coordinate channel.
    threshold : :obj:`float`, :obj:`str`, or None, optional
        Reporting threshold each study applied, on the z scale. This is the one number a
        silence cannot do without: "study k reported nothing here" is evidence about the
        effect only against how large an effect k would have needed to report it.

        A string names a metadata field holding per-study values, which is the right answer
        when the papers state their thresholds; studies missing the field fall back to the
        median of those that have it. A float applies one cut to every study. None assumes
        two-tailed p < .001, the most common screening threshold in the literature.

        **It is no longer inferred from the reported heights, because the heights are no longer
        read.** The removed rules undid the order statistic on a study's smallest reported
        value, which cannot distinguish a voxelwise height cut from a cluster-forming one and
        overshot the second by about 1 z: against a true forming cut of z = 3.1 the inference
        returned 4.0, and the prevalences it produced were inflated at every site (a true 0.25
        reading 0.48, a true 0.50 reading 0.86) where a plausible constant recovered 0.21 and
        0.47. Assuming a constant is both simpler and more accurate than inferring one.

        What the choice costs divides sharply. ``g`` barely notices -- across a +1.1 z error it
        stays within 10% of its value at the correct threshold, non-monotonically. Supplying
        the real thresholds matters far more for ``prevalence``, which does not survive it:
        0.73, 0.96, 0.99, 1.00 as the threshold given is inflated by 0, 0.4, 0.8 and 1.1 z,
        against a true 0.60.

        This threshold is also what stops a quiet region reading as an effect of zero. It is
        converted onto the effect-size scale by :func:`reporting_cutoff_to_g` before the
        likelihood sees it -- a z of 3.29 is about 0.65 g at ``n = 30`` -- and the silences
        push :math:`\mu` down only as far as they can carry it, which is toward that cut
        rather than toward nothing. Leaving the cut on the z scale would put it eighteen
        sampling standard deviations out, make every silence certain whatever the effect, and
        take the whole coordinate channel inert.
    clamp_threshold : :obj:`bool`, default=True
        Lower each study's assumed threshold to its own smallest reported statistic, where the
        table carries one. Anything a study reported cleared its cut, so this is a hard
        inequality and not an inference: it can only move a cutoff *down*, only for a study
        whose table contradicts the assumption, and never below the truth. That is what
        distinguishes it from the retired ``"study-min"`` rule, which tried to recover the cut
        by undoing an order statistic it could not identify.

        This is the one remaining use of the reported statistics, and it is a bound on the
        threshold rather than a magnitude -- nothing here reaches the pooled effect size.

        On a simulator whose studies applied 2.4, 2.8, 3.29 and 3.8 z against an assumption of
        3.29, the clamp moved 9 of 20 studies' cutoffs and improved the rmse against a known
        truth by 0.014 where the truth is near zero (paired p = 0.0001) and 0.008 in the middle
        stratum (p = 0.015), with no change where the effect is largest (p = 0.79). It was
        **bit-identical** in both regimes where it should do nothing: where every study really
        applied the assumed cut, and where every study thresholded above it so the bound is
        vacuous -- which is also the thin-table case, a paper reporting only its strongest
        peaks. Safe to leave on; turn it off to hold an assumed threshold exactly.

        One caveat, recorded rather than relied on: in that last regime the *true* thresholds
        were worse than the too-low assumption (rmse 0.128 against 0.106 near zero), so a
        cutoff slightly below the truth is compensating for something. The likely cause is that
        this model treats a report as :math:`|g| \ge c` while a reported peak is
        :math:`|g| \ge c` **and** a local maximum, a strictly smaller event -- so the
        probability of reporting is overstated and a lower cut offsets it. Do not read the
        clamp as more accurate than a stated threshold; read it as a bound that cannot hurt.
    coverage_radius : :obj:`float`, default=20.0
        Radius, in mm, within which a reported focus counts as this study having said
        *something* about a voxel; a study with no focus inside it is silent there and
        contributes a censoring term. This is the only geometry left in the model: papers do
        not report cluster extent reliably, so the extent a focus stands in for has to be
        assumed rather than read, and assuming it is not the same as treating it as silence.

        Used only when ``selection_model="zero-inflated"``. ``g`` is insensitive to it at
        realistic focus counts; ``prevalence`` is not. Against a known prevalence the estimate
        rises monotonically with this radius at every true value (a true 0.50 reads 0.65, 0.73,
        0.76, 0.81 at 8, 14, 20 and 28 mm) and no radius recovers the truth: mean absolute
        error runs 0.17 to 0.21 over that range, 14 mm marginally best and 20 mm close behind.
        Left at 20 mm because the differences are small beside the bias itself. On dense focus
        tables ``prevalence`` saturates at 1.0 here.
    max_iter : :obj:`int`, default=25
        Maximum Newton iterations for the censored likelihood. Voxels reached by a single image
        do not converge at any value of this, and raising it does not help: their likelihood is
        flat in the magnitude over the whole plausible range, so the iteration count only
        decides which point on a plateau is reported.
    null_method : {"permute-images", "none"}, default="permute-images"
        How uncorrected p-values are obtained. ``g / se`` is not null-referenced -- the
        standard error treats :math:`\tau^2` as known and ignores the selection the censoring
        term is modelling -- so p comes from a randomization null instead.

        ``"permute-images"`` reassigns each image study's effect sizes among **its own** voxels
        and refits, holding the coordinate tables exactly as they are. The hypothesis is that
        within a study, effect size is unrelated to location. Because nothing moves between
        studies and the silence pattern is identical in the observed fit and in every
        permutation, each voxel keeps its own studies and its own censoring roster throughout
        and is referred to a null of its own -- the exchangeability the test needs
        (:footcite:t:`winkler2014permutation`). Whatever the censoring term contributes cancels
        between observed and null, which is why a null over the coordinates is neither needed
        nor available.

        The p-value is ``(1 + #{null >= observed}) / (1 + n_iters)`` and so cannot fall below
        ``1 / (1 + n_iters)``, which is where the default ``cluster_threshold`` of .001 sits
        unless ``n_iters`` is raised past 1000; familywise correction has no such floor. This
        is deliberately **not** a test of spatial convergence, which is what a null that
        relocates the foci -- ALE's and MKDA's -- would give instead, nor a test of whether the
        effect is zero, which is what sign-flipping the images would give.

        ``"none"`` returns ``p = 1`` everywhere, for inspecting the estimates at no cost.
    cluster_threshold : :obj:`float` or None, default=0.001
        Cluster-forming threshold, as an uncorrected p-value, for the cluster-level FWE null
        that :meth:`fit` builds alongside the voxel-level one. Set to None to skip it, which
        makes :meth:`correct_fwe_montecarlo` pay for a second pass over the permutations if
        cluster correction is then requested.
    n_iters : :obj:`int`, default=1000
        Permutations for the null. Each is a full refit, which makes this the dominant cost of
        the estimator. It also sets the resolution of the uncorrected p, which cannot fall
        below ``1 / (1 + n_iters)``.
    n_cores : :obj:`int`, default=1
        Processes used for the permutation null, which is where nearly all the time goes.
        ``-1`` uses every available core and is close to linear. The iterations run one block
        per core rather than one task per iteration, because the null is accumulated per voxel
        and shipping each iteration's whole map back would cost more than the refits.
    seed : :obj:`int`, default=0
        Seed for the permutation draws.
    memory, memory_level, generate_description
        As in every other :class:`~nimare.estimator.Estimator`.

```


---

## `_censoring_terms`, as its docstring carried it

The measurements behind the report limb, the restored `-1` pairs, and the profiling that decided against rewriting it.

```
def _censoring_terms(mu, cutoff_scaled, twice_cutoff_scaled, inv_sigma, inv_sigma_sq, sign):
    r"""Probability of each observed *reporting indicator*, and the pieces of its derivatives.

    A coordinate table carries one bit per study per voxel: the study reported something near
    here, or it did not. Both values of that bit are informative, and the two are complementary
    probabilities of the same event, so they are computed together and told apart by ``sign``:
    ``+1`` for a silent pair, whose probability is :math:`P(|g| < c \mid \mu)`, and ``-1`` for
    a pair that reported, whose probability is :math:`1 - P(|g| < c \mid \mu)` with its height
    discarded.

    **The ``-1`` limb's probability is overstated, and the error peaks mid-window.** A paper
    reports a voxel only if it cleared ``c`` *and* the value there was a local maximum, which
    this :math:`P(|g| \\ge c)` does not require. Measured against a known truth, the model's
    reporting rate exceeded the observed one by 1.00, 1.78, 1.08 and 1.13 at true ``g`` of 0.2,
    0.4, 0.6 and 0.8 -- non-monotone, so no uniform reweighting of the limb can absorb it, and
    raising its probability to a power was measured making every focus worse. Correcting it needs
    the survival of a suprathreshold local maximum, and that needs a field smoothness this model
    does not carry. See ``CBES`` under "Why the indicator and not the heights".

    **Dropping the ``-1`` pairs biases the magnitude down, and hard.** They used to contribute
    nothing at all, on the reasoning that a coordinate carries no usable height -- but omitting
    them leaves the *silent* pairs as the only evidence about the indicator, so the model reads
    the observed silence fraction against a denominator that excludes every study that reported.
    On a one-voxel likelihood with the truth known exactly, 20 studies of which 2 supply images
    and a cutoff of 0.60 g, that returned a mean :math:`\hat\mu` of 0.351 for a true 0.500
    (rmse 0.192); with the indicator restored, 0.524 (rmse 0.129).

    Returned as one dict because the E step and the M step both need these at the same ``mu``.
    **This is the most expensive single function in the estimator** -- 59% of a whole-brain fit
    with a permutation null, which at the default ``n_iters`` is most of the wall clock, since
    each permutation runs a full EM.

    Profiled rather than assumed, and the assumption was wrong: the cost is *not* spread evenly
    over memory-bound kernels. On 600,000 pairs the two ``ndtr`` calls take 10.7 ms and 8.1 ms
    against 3.0 ms each for the two densities and 0.4 ms for an arithmetic pass, so **45% of
    the time is two normal CDFs** and no rearrangement of the surrounding algebra reaches it.
    Three were measured -- one reciprocal in place of two divisions, ``second`` obtained from
    ``first`` through
    :math:`u\phi(u) - l\phi(l) = u(\phi(u) - \phi(l)) + k\phi(l)`, and both together -- and
    they came out at 1.00x, 1.09x and 1.05x with up to 5e-14 of drift. Not worth the churn.

    Nor can the lower tail be dropped to save its CDF: its median contribution is 2e-5 of the
    silent probability, which sounds negligible, but its maximum is 0.30 and it exceeds 1% of
    the score's numerator for 38% of pairs. The two levers that do work are ``n_cores``, the
    null being a thousand independent fits, and the compaction in :meth:`CBES._fit_chunk`,
    which shrinks the array as voxels retire.

    Everything that does not move between EM iterations is passed in already divided: ``mu`` is
    the only argument that changes, so ``cutoffs / sigma`` and the reciprocals are hoisted to
    the caller. The lower tail looks negligible and is not -- at a typical cutoff it is a third
    of the score's numerator -- so it is kept. The arithmetic writes into its own temporaries
    wherever numpy allows it, each avoided temporary being hundreds of megabytes of traffic.
    """
    # upper = (c - mu) / sigma;  lower = (-c - mu) / sigma = upper - 2c/sigma
    upper = mu * -inv_sigma
```


---

## reporting_cutoff_to_g: the degrees-of-freedom sensitivity table

```
    **The assumed degrees of freedom are load-bearing, because a threshold sits far into the
    tail where that map is steep.** Holding ``n`` at 30 and varying only the assumed residual
    degrees of freedom:

    ============  =======  =======  ========  =========  ======
    cutoff z       df=29    df=60    df=120    df=1000   spread
    ============  =======  =======  ========  =========  ======
    3.30           0.653    0.626     0.614      0.604    1.08x
    4.00           0.830    0.776     0.752      0.733    1.13x
    5.00           1.133    1.009     0.959      0.918    1.23x
    6.00           1.522    1.273     1.178      1.105    1.38x
    ============  =======  =======  ========  =========  ======

    The effective degrees of freedom of a published map are frequently *above* ``n - 1`` --
    variance smoothing raises them, and some mixed-effects tools do that deliberately -- and
    papers seldom state them, so a threshold read off a published z is likely to be placed a
    little too high, which makes a silence look less surprising than it was.
```


---

## _study_cutoffs_z: what the retired inference cost

```
        It can no longer be inferred. The old rules read it off the smallest reported statistic,
        undoing the order statistic for the number of peaks; with the heights no longer read,
        there is nothing to read it from. So it is **supplied or assumed**, which is also the
        honest position: the earlier inference was measured at 0.201 of prevalence error against
        0.008 for a fixed constant on cluster-extent tables, because it cannot tell a
        cluster-forming cut from a voxelwise one and overshoots the first by about 1 z.
```


---

## _indicator_entries: the extent and limb measurements

```
        Measured on a known truth, asserting the report across the sphere gave an rmse of
        0.457 against 0.114 for the named voxel alone.

        Omitting the ``-1`` limb entirely is worse than including it at the named voxel: the
        silent pairs are then the only evidence about the indicator, so the model reads the
        observed silence fraction against a denominator that excludes every study that
        reported, and over-shrinks. Where the truth is largest that cost 0.091 of rmse against
        0.074 and -0.060 of bias against -0.042.
```


---

## `_hartung_knapp_se`, removed

Reachable only with `selection_model="none"`, which makes the coordinate tables inert -- so it only ever corrected the standard error of a plain image-based random-effects fit, a correction `nimare.meta.ibma` already offers as `small_sample_correction='knapp-hartung'`. Kept here in case the derivation is wanted again.

```python
def _hartung_knapp_se(*, g_hat, sum_a, sum_a_g2, n_eff, covered, fallback):
    r"""Hartung-Knapp-Sidik-Jonkman standard error of a kernel-weighted pooled estimate.

    The model-based SE treats :math:`\hat{\tau}^2` as if it were the true heterogeneity, so its
    intervals are too short exactly when heterogeneity is large and the studies are few. HKSJ
    replaces it with the weighted spread of the studies about the pooled value,

    .. math::

        \mathrm{SE}^2 = \frac{\sum_k a_k (g_k - \hat{g})^2}{(k_{\mathrm{eff}} - 1)\sum_k a_k},
        \qquad a_k = \frac{w_k}{s^2_k + \tau^2},

    on :math:`k_{\mathrm{eff}} - 1` degrees of freedom, which gives much better interval
    coverage than the model SE when heterogeneity is large and the studies are few.

    :math:`k_{\mathrm{eff}}` is Kish's :math:`(\sum w)^2 / \sum w^2` -- the ``n_eff`` map --
    and not :math:`\sum w`. The two agree when every weight is one, but only Kish's form is
    invariant to rescaling the weights: at a voxel reached only by distant foci the weights sum
    to less than one, and using that as a study count sends the degrees of freedom to zero.
    Voxels with no effective spread to measure keep the model-based value.
    """
    se = np.array(fallback, dtype=float, copy=True)
    usable = covered & (n_eff > 1.0) & (sum_a > 0)
    if not np.any(usable):
        return se
    # sum a (g - ghat)^2, from the identity noted at the call site. Clipped at zero: the two
    # terms are close where the studies agree, so rounding can make the difference negative.
    residual = np.clip(sum_a_g2[usable] - g_hat[usable] ** 2 * sum_a[usable], 0.0, None)
    se[usable] = np.sqrt(residual / ((n_eff[usable] - 1.0) * sum_a[usable]))
    return se
```


---

## `_hartung_knapp_se`, removed

Reachable only with `selection_model="none"`, which makes the coordinate tables inert -- so it only ever corrected the standard error of a plain image-based random-effects fit, a correction `nimare.meta.ibma` already offers as `small_sample_correction='knapp-hartung'`.

```python
def _hartung_knapp_se(*, g_hat, sum_a, sum_a_g2, n_eff, covered, fallback):
    r"""Hartung-Knapp-Sidik-Jonkman standard error of a kernel-weighted pooled estimate.

    The model-based SE treats :math:`\hat{\tau}^2` as if it were the true heterogeneity, so its
    intervals are too short exactly when heterogeneity is large and the studies are few. HKSJ
    replaces it with the weighted spread of the studies about the pooled value,

    .. math::

        \mathrm{SE}^2 = \frac{\sum_k a_k (g_k - \hat{g})^2}{(k_{\mathrm{eff}} - 1)\sum_k a_k},
        \qquad a_k = \frac{w_k}{s^2_k + \tau^2},

    on :math:`k_{\mathrm{eff}} - 1` degrees of freedom, which gives much better interval
    coverage than the model SE when heterogeneity is large and the studies are few.

    :math:`k_{\mathrm{eff}}` is Kish's :math:`(\sum w)^2 / \sum w^2` -- the ``n_eff`` map --
    and not :math:`\sum w`. The two agree when every weight is one, but only Kish's form is
    invariant to rescaling the weights: at a voxel reached only by distant foci the weights sum
    to less than one, and using that as a study count sends the degrees of freedom to zero.
    Voxels with no effective spread to measure keep the model-based value.
    """
    se = np.array(fallback, dtype=float, copy=True)
    usable = covered & (n_eff > 1.0) & (sum_a > 0)
    if not np.any(usable):
        return se
    # sum a (g - ghat)^2, from the identity noted at the call site. Clipped at zero: the two
    # terms are close where the studies agree, so rounding can make the difference negative.
    residual = np.clip(sum_a_g2[usable] - g_hat[usable] ** 2 * sum_a[usable], 0.0, None)
    se[usable] = np.sqrt(residual / ((n_eff[usable] - 1.0) * sum_a[usable]))
    return se
```

## When do the coordinates stop helping?

`when_do_coordinates_start_hurting.py`, pain, published tables, 8 splits, paired: the same
permutations and the same held-out reference at every image count, so only the division of the
working half into image donors and table donors changes.

| images | tables | rmse CBES | rmse pool | diff | p | err at top CBES | pool | diff | p |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 9 | 0.271 | 0.363 | −0.092 | 0.000 | −0.047 | +0.141 | −0.188 | 0.000 |
| 2 | 8 | 0.230 | 0.269 | −0.040 | 0.002 | +0.005 | +0.137 | −0.132 | 0.000 |
| 3 | 7 | 0.210 | 0.217 | −0.007 | 0.432 | −0.038 | +0.054 | −0.092 | 0.003 |
| 4 | 6 | 0.207 | 0.196 | +0.011 | 0.152 | −0.043 | +0.020 | −0.064 | 0.017 |
| 5 | 5 | 0.198 | 0.171 | +0.027 | 0.001 | −0.023 | +0.020 | −0.044 | 0.037 |
| 6 | 4 | 0.200 | 0.161 | +0.039 | 0.001 | −0.020 | −0.014 | −0.006 | 0.524 |
| 7 | 3 | 0.203 | 0.147 | +0.056 | 0.000 | −0.017 | −0.027 | +0.011 | 0.290 |
| 8 | 2 | 0.203 | 0.143 | +0.060 | 0.000 | −0.014 | −0.041 | +0.027 | 0.022 |
| 9 | 1 | 0.196 | 0.142 | +0.054 | 0.000 | −0.015 | −0.058 | +0.043 | 0.002 |

**The crossing on rmse is 3 images; on bias at the top it is 6.** The image pool's rmse keeps
falling, 0.363 to 0.142, because its error is all variance. CBES flattens at about 0.20 — that
floor is the coordinate channel's bias, and it does not shrink with anything. Which crossing
matters depends on the question: the pool swings from +0.141 to −0.058 at the top decile as
images accumulate, while CBES holds −0.015 to −0.047 throughout, so for a calibrated magnitude at
the peaks the coordinates pay for twice as long as they do for whole-map error.

### The crossing condition, in closed form

`micro_when_coordinates_hurt.py`. Two channels: images unbiased with variance `1/(k n)`, falling
with image count; coordinates with variance `1/(m I)` that falls with table count, plus a bias
`b` from the assumptions they are read through, which does not fall with anything. Inverse-
variance weighting puts `w = R/(1+R)` on the coordinates, `R = m I/(k n)`. Writing out the mean
squared error of the combination and asking when it beats `1/(k n)`, every term in `R` cancels:

    coordinates help   <=>   b^2 < 1/(m I) + 1/(k n)

Squared channel bias below the **sum** of the two channels' variances, which is the variance of
their difference. So the crossing is a Hausman statistic, `Q = (g_c - g_i)^2/(var_i + var_c) < 1`,
and `k* = 1/(n (b^2 - 1/(m I)))` — **infinite whenever the coordinate channel's own variance
already exceeds its squared bias**, so a well-specified threshold never hurts at any image count.

Verified in the scalar bed, where the reports are generated at a true cut and the estimator reads
them at a wrong one. The `Q > 1` crossing matched the measured MSE crossing in all five
scenarios, and predicted `k*` 1.7/0.4/0.1 landed on observed 2/1/1. Median Q is 0.455 when the
bias is exactly zero, which is the median of a chi-square on 1 df, as it must be.

### Q on real data: read it at the strong voxels, and do not trust the absolute threshold yet

In terms of what CBES reports, with `w = coordinate_share` and `T = 1/se^2`:

    Q = (g - g_selection_off)^2 / se^2 * (1 - w)/w

Two things went wrong on first contact with real data, both in the estimator of Q rather than
the condition. **Algebraically `Q -> 0` as `w -> 0`, because `delta = w (g_c - g_i)`; but delta is
a difference of two separately fitted maps and carries float noise that does not shrink with
`w`, so dividing it by `w^2` manufactures values in the thousands.** Restricting to `w >= 0.10`
fixes most of it, and is not enough at 9 images where only 6% of voxels act and the reading is
2e17. And the **map median does not locate the crossing** — it sits at 0.09 to 0.38 throughout,
far below 1 — while **Q at the strongest decile does**: 0.966, 0.661, 0.844, 0.996, 1.067, 1.136,
1.031, 1.051, 3.790 for 1 to 9 images, crossing 1 between 4 and 5, where the rmse penalty turns
significant.

That agreement is provisional, and the reason turned out to be decisive. Under no channel bias Q
has median 0.455, and an se inflated by `f` deflates Q by `f^2`. The scalar bed was given an
`INFLATE` knob to measure how much that matters, at the marginal scenario whose true crossing is
two images:

| se inflation | median Q at k = 1, 2, 4, 8, 16, 32 | Q > 1 crosses | truth |
|---|---|---|---|
| 1.0 | 0.809, 1.077, 2.061, 2.951, 4.388, 5.403 | **2** | 2 |
| 1.5 | 0.359, 0.479, 0.916, 1.311, 1.950, 2.401 | **8** | 2 |
| 2.2 | 0.167, 0.223, 0.426, 0.610, 0.907, 1.116 | **32** | 2 |

**A 1.5x inflation moves the recommended image count from two to eight, and a 2.2x from two to
thirty-two.** Across the three bias scenarios, the image count at which `Q > 1` fires:

| se inflation | cut 0.70 (true 2) | cut 0.85 (true 1) | cut 1.10 (true 1) |
|---|---|---|---|
| 1.0 | **2** | **1** | **1** |
| 1.5 | 8 | **1** | **1** |
| 2.2 | 32 | 4 | **1** |

At the bottom of the measured `se/sd` range the criterion is exact; at the top only the most
extreme bias still locates the crossing. It is not that the statistic is fragile only in the
marginal regime -- across the measured range it is unusable absolutely anywhere.

The correction is exact in the algebra: reported `Q = Q_true / f^2`, so the threshold is
`Q > 1/f^2`, and at `f = 1.5` the threshold 0.444 recovers `k = 2` from the inflated series. It
still cannot be applied to pain, because `f = 1.5` there would put the crossing at one image
against a directly measured tie at three. Either the inflation at those voxels is well below 1.5
or `Q@top` is not deflated the way a uniform model says; `se/sd` is not known per voxel, so the
question is not decidable from what is measured.

**So Q is ordinal, not absolute.** Its rise with image count, and `coordinate_share`'s fall, flag
the direction reliably. The decision itself should be measured directly, which is cheap and needs
no calibrated `se`: fit with `selection_model="none"` and with it on at the real image count and
score both against a held-out subset. A calibrated `se` is the blocker for any absolute form of
the criterion, which is the practical reason the `se` work matters.

### The radius is an exclusion zone, and a wider one means *less* silence

First, the semantics, because getting them backwards inverts every prediction. From
`effectsize.py`, `sign = where(at_focus, -1, where(reached, 0, 1))`: the named voxel is reported,
a voxel within `coverage_radius` of one of that study's peaks but not named gets **no indicator
at all**, and everything else the study examined is silent. So the radius is an *exclusion zone*
around reported peaks, not the reach of the silence. **A wider radius means less silence**, hence
less coordinate information and less misassignment bias -- the opposite of the reading that comes
naturally from the parameter's name.

The measured `coordinate_share` confirms which way it runs: 0.472 at 4 mm falling to 0.443 at
20 mm with one image, 0.229 to 0.201 with five. It could only fall if wider means less silence.

`does_shrinking_the_omission_help.py`, pain, published tables, 6 splits, `coverage_radius` swept
4 to 20 mm:

| radius | rmse, 1 image | err at top | rmse, 5 images | err at top |
|---|---|---|---|---|
| 4 mm | 0.282 | −0.075 | **0.199** | −0.038 |
| 8 mm | 0.282 | −0.072 | 0.199 | −0.034 |
| 14 mm | **0.281** | −0.062 | 0.201 | −0.019 |
| 20 mm | 0.281 | −0.049 | 0.204 | −0.002 |

The nominal best radius is 14 mm at one image and 4 mm at five, but **that shift is not a
finding**: the rmse spread across radii is 0.001 to 0.005 over six splits, which is noise, and
0.4% to 2.5% of rmse against the 35% the image count moves. Claiming it confirmed a predicted
inward shift would have been reading a direction into noise, and with the semantics the right way
round the prediction pointed outward anyway.

What is clean is the **top-decile bias, monotone in the radius in both regimes**: −0.075 to
−0.049 at one image, −0.038 to −0.002 at five, improving as the radius widens. That is exactly
what less misassigned silence should do. So the defensible statement is the narrow one: widening
the radius trades coordinate information for centring at the peaks, and rmse is nearly flat in
it. There is no single best radius, and the knob is worth far less than the image count.

The asymmetry between silence and reporting is not a knob at all, and should not be softened. A
report contributes `P(|g| > c)`, which saturates: once `mu sqrt(n)` is past `c sqrt(n)` the factor
is nearly 1 and its derivative in mu nearly 0. A report genuinely says only that the effect is
not small, so the confidence at a reported voxel is carried by the images. Down-weighting the
omission would add nothing where coordinates are reported, remove information where they are
silent, and cost the observed-information SE its meaning as the curvature of a real likelihood.

## A working-memory studyset from the NeuroStore base-study endpoint

`fetch_working_memory.py` and `wm_sdm_vs_cbes.py`. Pain is one hand-curated collection and HCP is
one task cut into synthetic studies; neither says what happens on a studyset assembled the way a
meta-analyst assembles one, by asking an index for the studies on a topic and taking whatever
maps they shared. Eighteen of NeuroStore's 684 base studies name a working-memory paradigm in
their title, abstract or keywords and carry group t or z maps.

Three things about that studyset are worth more than the estimator comparison run on it.

**The corpus is barely coherent.** Median pairwise spatial correlation between the 18 study maps
is +0.081, 71% positive; curating to the 13 that are activation contrasts rather than searchlight
decoding accuracy, white-matter FA association or resting fALFF moves it only to +0.087. A
held-out reference pooled over such maps is mostly noise, which caps what any comparison here can
resolve — measured r against the reference is 0.15 uncurated and 0.20 to 0.26 curated, against
0.5-plus on pain.

**Sample sizes span 21 to 1369**, so one study carries more inverse-variance weight than all the
others together. In 17 of 24 splits the largest held-out study holds over 80% of the reference
weight, making the reference very nearly that one study. Any summary averaged over that and the
pooled-reference splits describes neither, so the bed reports them separately.

**Most studies do not survive their own correction.** Only 2 to 5 of 5 table donors report
anything under cluster/max, and the two maps whose |z| never exceeds 1 report nothing at all,
correctly dropping out.

The consequence for CBES is a clean null result. Curated, 24 splits, 1 image:

| estimate | r | rank r | AUC | magnitude |
|---|---|---|---|---|
| images only | +0.201 | +0.233 | 0.622 | 0.901 |
| CBES g | +0.196 | +0.231 | 0.623 | 0.779 |

Identical on pattern, in both reference regimes separately as well as pooled. The estimator's own
diagnostic says why: `coordinate_share` at the strongest decile is 0.000 in three of six splits.
With only two to five studies contributing an indicator at all, the coordinate channel is roughly
ten times thinner than pain's 267 transcribed peaks over ten studies, and there is nothing for it
to add. The one place CBES still differs is calibration — magnitude 0.78 against the pool's 0.90.

### Two harness defects this bed uncovered

**NeuroVault changed its `map_type` encoding** from short codes to `T map`/`Z map`. Both
`fetch_corpus.py` and the first version of `fetch_working_memory.py` filtered on `("t", "z")`,
matched nothing, and exited zero with an empty corpus — the worst way for a fetcher to fail. Both
now normalise.

**`sdm_parse` cannot find a silent study whose name carries two or more digits.** `c16` and `c99`
are reported missing while `c16.no_peaks.txt` sits next to the table; `c3`, `c7` and any
digit-free name are found. Studies that report peaks are unaffected, so the bug bites only the
silent studies, which is the arm that matters here, and it is why the pain bed never hit it: its
published tables always had peaks. Study labels are now letters only. Worth noting that the first
three diagnoses were all wrong — a zero-byte file, then row order, then a nine-study limit — and
one run passed on input that fails deterministically, which was stale state in the install
directory, not evidence.

## The familywise defect was the bed, not the estimator

`fpr_after_redesign.py` on the whole simulated cube: uncorrected 0.0533, **familywise 0.375**
against a nominal 0.05, binomial SE 0.034. A defect that spares the marginal and wrecks the
maximum lives in the dependence structure, and for this estimator there is only one candidate.

The argument, which is sound and worth keeping: the fit is voxelwise, so permuting an image's
values among its own voxels changes no value, only which voxel each is paired with. Were the
indicator spatially constant, the permuted `z` map would be a permutation of the observed one and
the maxima would agree exactly. They can differ only because the indicator is held fixed while
the image moves — so the observed map is scored on the *actual* pairing of image magnitude to
indicator, and every permutation on a random one. The null is therefore exchangeable if and only
if, across voxels, a study's image magnitude is independent of how many studies were silent there.

`is_the_image_permutation_exchangeable.py` measures that in seconds, with no CBES, no permutations
and no correction — just the simulator's studyset and a spatial query. Mean Spearman rho −0.115,
negative in all eight realisations, per-realisation p from 4e-09 to 1e-29. The competing
explanation fails on direction: scrambling destroys smoothness, which gives the permuted field
*more* effectively independent tests and a *larger* maximum, making a test conservative rather
than liberal.

**And then the interior-shell check overturned the diagnosis.**

| voxels used | mean rho(&#124;g&#124;, silent count) |
|---|---|
| all of the simulated cube | −0.125 |
| beyond 2 voxels (8 mm) from the edge | **−0.0005** |
| beyond 4 voxels (16 mm) | +0.032 |

with `corr(|g|, edge distance) = −0.212` and `corr(silent count, edge distance) = +0.469`. Both
are what smoothing a field inside a finite cube does: the kernel runs off the edge and leaves the
outer shell with inflated variance and inflated peak density, which couples magnitude to peak
density there and nowhere else.

So `permute-images` is non-exchangeable **in the bed that measured the false-positive rate**, for
reasons belonging to the bed. `fpr_after_redesign.py` now takes `ERODE`, and the rate is being
re-measured with 8 mm off every face. What the erosion cannot settle either way: a real brain mask
also has a boundary, and smoothing near it inflates variance the same way, so whether a real
analysis's familywise correction is affected is a question about data preparation that needs
measuring on a real collection.

This is the third harness defect in this project mistaken for a model defect. The check that
overturned it took ninety seconds and should have run before any mechanism was reported at all.

## SDM-PSI on the working-memory studyset, and HCP emotion as a second task

Both run at one image, the configuration under question, with the same parity: the same studies
supply images, the rest supply coordinates, and the reference is held out from both.

**Working memory** (curated to 13 activation contrasts, 4 splits, cluster/max extraction):

| estimate | r | rank r | AUC |
|---|---|---|---|
| images only | +0.337 | +0.369 | 0.733 |
| CBES g | +0.328 | +0.363 | 0.732 |
| SDM-PSI | +0.075 | +0.139 | 0.580 |

CBES tracks the donated image and SDM-PSI is far behind both. SDM's magnitude column is not
comparable, its coefficient not being on a `g` scale. This needed a harness fix first: `sdm_parse`
cannot find a silent study whose name carries two or more digits, so every `no_peaks` study was
being dropped.

**HCP EMOTION_FACES**, 30 subjects x 16 studies, 3 splits, scored against 306 held-out subjects:

| estimate | r | rank r | AUC | mag ratio |
|---|---|---|---|---|
| SDM-PSI coeff | +0.482 | +0.463 | 0.737 | 0.15 |
| images only | +0.879 | +0.676 | 0.984 | **0.96** |
| CBES g | **+0.882** | **+0.680** | 0.984 | 0.68 |
| CBES g_marginal | +0.312 | +0.439 | 0.552 | 0.28 |

This reproduces the prevalence-1 failure on a second HCP task, independently of MOTOR_LH. CBES
matches the image pool on pattern and beats SDM-PSI on every metric, but the pool recovers 0.96 of
the true magnitude where `g` recovers 0.68 -- on MOTOR_LH, 0.85 against 0.63. The cause is the same
and it is not subtle: every HCP study genuinely carries the effect, so there is no absence to find
and the correction can only shrink. The fitted prevalence is 0.626 at the median against a true
1.0 and **0.077 at the strongest decile**, so it is worst exactly where the estimate matters most.

`g_marginal` is unusable here for the same reason -- it multiplies by that biased prevalence, and
its AUC of 0.552 is barely above chance on a task where the image pool reaches 0.984.

## Global-null error rates on the redesigned null, all three arms

`fpr_after_redesign.py`, 40 simulations, 200 permutations, nominal 0.05, binomial SE 0.034:

| arm | uncorrected | familywise |
|---|---|---|
| 20 studies, 2 images | 0.0533 | **0.375** |
| 20 studies, 1 image | 0.0496 | **0.200** |
| 12 studies, 2 images | 0.0534 | **0.350** |

The voxelwise rate is calibrated in every arm. The familywise rate is four to seven times
nominal, and scales with the image count, which is what the exchangeability argument predicts --
the non-exchangeable pairing enters only through the image channel.

**Eroding the boundary removes most of it**, confirming that the coupling was the cube's edge.
`ERODE=2` drops 8 mm from every face, leaving 3,375 of 6,859 voxels, and changes nothing about how
the data are generated -- the simulator works from `noise_extent`, so the same studies report the
same peaks and the edge voxels are simply not read:

| arm | uncorrected | familywise, whole cube | familywise, interior |
|---|---|---|---|
| 20 studies, 2 images | 0.0508 | 0.375 | **0.075** |
| 20 studies, 1 image | 0.0498 | 0.200 | **0.025** |
| 12 studies, 2 images | 0.0501 | 0.350 | **0.150** |

With a binomial standard error of 0.034, the two 20-study arms are within noise of nominal: 0.075
is 0.7 SE above 0.05, and 0.025 is conservative. **These are the numbers to quote**, and they say
the estimator's null is sound where the field is stationary.

**A residual remains in the 12-study arm at 0.150**, which is 2.9 SE above nominal and not the
boundary. Adding two arms discriminates what drives it:

| arm | image share of the roster | familywise |
|---|---|---|
| 20 studies, 1 image | 5% | **0.025** |
| 20 studies, 2 images | 10% | 0.075 |
| 12 studies, 1 image | 8.3% | **0.100** |
| 12 studies, 2 images | 16.7% | **0.150** |

Not purely the image share -- 8.3% gives 0.100 where 10% gives 0.075, and those two are within
noise of each other. It reads as roughly additive: dropping 20 studies to 12 adds about 0.075, and
a second image adds about 0.05. **The study count is the larger driver**, with the image count
secondary, which is consistent with the non-exchangeable pairing entering through the image channel
while the max-statistic null gets coarser with fewer studies.

So the shipping statement has to be that the familywise correction is calibrated at twenty studies
and anti-conservative below that, with these numbers, rather than calibrated. Whether the
degenerate-null guard is firing at all in the 12-study arm is a separate question: `refused` reads
0 for all 40 simulations, and the first attempt to record the distinct-maxima count returned `nan`
because `correct_fwe_montecarlo` writes its diagnostics onto the result's own estimator rather than
the one the bed holds. Re-measuring with that fixed.

## One image against a wall of tables: where the estimator earns its keep, and its bias floor

jdkent's configuration, and the one that actually arises: the literature is fixed and large, and
what a meta-analyst gains over time is images. `build_pain_table_corpus.py` takes the NeuroStore
2026-09 release (32,444 studies, 871,671 coordinates), keeps the 1,522 studies naming pain, and
`one_image_many_tables.py` holds that wall constant while adding NIDM pain images one at a time.
Ten of the 21 NIDM images are held out as a reference that is fixed within a draw, so every row is
paired.

**The contamination guard earned its place.** The NIDM study ids are anonymised, so overlap cannot
be matched by name, but the NIDM collection carries each study's published peaks -- so a NeuroStore
study reporting the same peaks in the same places is the same study. Twelve were found and dropped.
Without that the reference would have been inside the coordinate channel.

Three draws, rmse against the held-out pool, paired:

| tables | 1 img | 2 | 3 | 5 | 8 | 11 | crossing |
|---|---|---|---|---|---|---|---|
| 25 | **−0.118** (p=.020) | −0.010 | −0.005 | −0.004 | +0.001 | +0.003 | ~8 images |
| 200 | **−0.124** (p=.051) | +0.007 | +0.011 | +0.011 | +0.021 | +0.025 | 2 images |
| 1443 | **−0.131** (p=.057) | +0.020 | +0.026 | +0.050 | +0.078 | **+0.085** (p=.022) | 2 images |

**The crossing moves to fewer images as the wall grows** -- eight at 25 tables, two at 200 and at
1,443 -- which is what `b^2 < 1/(m I_c) + 1/(k n)` predicts once the first term vanishes. And the
estimate stops responding to images at all: at 1,443 tables rmse goes 0.236 to 0.224 across one to
eleven images while the pool falls 0.366 to 0.138. With `coordinate_share` at 0.98, eleven images
cannot outvote the wall. **A large table wall does not average away its own bias; it protects it.**

At eleven images the ranking inverts by wall size -- 0.141 at 25 tables, 0.163 at 200, 0.224 at
1,443 -- so more coordinate studies is strictly worse once images are available, which is the
opposite of the usual meta-analytic intuition.

### The rmse win at one image is shrinkage, not centring

| tables | CBES err at top, 1 image | pool |
|---|---|---|
| 25 | −0.193 | +0.065 |
| 200 | −0.289 | +0.065 |
| 1443 | −0.326 | +0.065 |

reaching −0.469 at 1,443 tables and eleven images. So at one image CBES beats the pool on rmse by
shrinking an estimate that is mostly noise, while being badly biased low at the strong voxels --
where the pool is nearly unbiased. **The over-shrinkage grows with the wall.** This is the opposite
of the pain-collection result, where CBES was *better* centred at the top (+0.005 against +0.137),
and the difference is the coordinate set: there, nine tables of the same studies' own published
peaks; here, 1,443 foreign studies.

One explanation was checked and is wrong in direction. 56% of the corpus's studies report a peak
below the assumed z of 3.09, so the assumed cut is too high -- but a cut that is too high makes the
bound `|g| <= c` *weaker*, which would shrink less, not more. Fixing it would deepen the problem.

### It is not a population mismatch, and the cause is an assumption I made

A population mismatch was the live explanation: `mu` is the effect among studies whose effect is
non-null, and 1,443 heterogeneous pain studies (analgesia, chronic-pain contrasts, modulation)
genuinely do not carry this contrast's effect, so their silence is truthful for them while the
reference is a different population. That predicts **small `pi` with `mu` staying large**. The
prevalence diagnostic falsifies it. At the top decile, against a truth of 0.622:

| tables | images | pi | pi at top | &#124;g&#124; at top |
|---|---|---|---|---|
| 25 | 11 | 0.658 | **0.933** | **0.472** |
| 200 | 11 | 0.487 | 0.506 | 0.312 |
| 1443 | 1 | 0.490 | 0.324 | 0.232 |
| 1443 | 11 | 0.425 | **0.295** | **0.159** |

The wall drags *both* down together. At 1,443 tables `pi` at the top is 0.295 **and** `g` is 0.159,
a quarter of the truth; at 25 tables with eleven images the fit recovers `pi` at 0.933 and `g` at
0.472, close to right.

The cause is almost certainly the bed rather than the estimator, and it is the identifiability
result already on record. Every one of the 1,443 coordinate studies was given the same assumed
`n = 20`, because the release reports a sample size for 0.5% of its analyses. A silence constrains
`(pi, mu)` through one scalar per distinct *(sampling error, cutoff)* pair, and only sample-size
spread separates the components -- so 1,443 studies at one shared `n` and one clamped threshold
supply **one constraint repeated 1,443 times**. That is a tight ridge in `(pi, mu)` with nothing
fixing the position along it, and the fit slides down it to low `pi` and low `mu` while
`coordinate_share` at 0.98 stops the images pulling back.

**That prediction was wrong too, and in an informative direction.** `NSPREAD=1` draws each
coordinate study's `n` from 12 to 120, and it makes everything worse:

| 1,443 tables | images | rmse CBES | pool | diff | pi at top | &#124;g&#124; at top | share |
|---|---|---|---|---|---|---|---|
| n = 20 fixed | 1 | 0.236 | 0.366 | −0.131 | 0.324 | 0.232 | 0.982 |
| n drawn 12–120 | 1 | 0.306 | 0.408 | −0.102 | **0.100** | 0.156 | 0.998 |
| n = 20 fixed | 11 | 0.224 | 0.138 | +0.085 | 0.295 | 0.159 | 0.937 |
| n drawn 12–120 | 11 | 0.291 | 0.120 | **+0.172** | **0.019** | 0.236 | **1.000** |

I had the bound's direction backwards. A silence from a study with `n = 120` says
`|g| < 3.09/sqrt(120) = 0.28`, where `n = 20` says `|g| < 0.69`: **a larger assumed sample size
makes each silence a tighter bound.** Spreading `n` upward therefore adds hundreds of very tight
bounds, the likelihood is told the effect is small nearly everywhere, and `pi` collapses to 0.019
while `coordinate_share` reaches 1.000.

### The structural reason: the wall is silence, not information

Counting what a 1,443-study wall actually contributes per voxel, at `coverage_radius` 20 mm:

| | studies saying something | silent |
|---|---|---|
| median voxel | 316 (21.9%) | 1,127 |
| 99th percentile voxel | 651 (45.1%) | 792 |
| best voxel | 724 (50.2%) | 719 |

and only **8.8% of voxels are named by any study at all, with at most 4 studies naming the same
voxel.** The report limb -- the only limb that says the effect here is at least this big -- is four
studies against a thousand silences, a ratio near 1:280. On the pain collection, where the
estimator works well, a named voxel has one report against about eight silences.

So the count of coordinate studies is not the currency. What decides whether a coordinate corpus
can carry a magnitude is the ratio of studies reporting *at a voxel* to studies silent there, and
adding topically-related-but-different studies makes that ratio worse rather than better. Silence
can shrink a noisy estimate in the right direction, which is exactly what the one-image win is; it
cannot supply a magnitude. That is the honest statement of what the extreme configuration does.

## Widening the report limb fixes the wall's over-shrinkage

jdkent's suggestion: assert the report over a small radius, 4 to 8 mm, so a peak says the effect is
at least this large *somewhere here* rather than at one named voxel, with the ring out to
`coverage_radius` still carrying no indicator. Three zones rather than two.

The sweep on record said no -- rmse 0.070 at the named voxel against 0.113 at 4 mm -- but that was
measured on a small collection where the report limb was already well fed, about one report per
eight silences. It does not transfer to a corpus running at 1089:1, and the bias widening adds is
*positive*, which is the direction the wall's error needs.

Shipped as `report_radius`, default `None` (named voxel, bit-identical). On the 1,443-study pain
corpus against one image, three draws:

| report radius | rmse | against pool | error at top | &#124;g&#124; at top (truth 0.622) | pi at top |
|---|---|---|---|---|---|
| named voxel | 0.234 | −0.137 (p=.037) | **−0.303** | 0.250 | 0.289 |
| 4 mm | **0.216** | **−0.155** (p=.021) | −0.122 | 0.450 | 0.306 |
| 6 mm | 0.256 | −0.115 (p=.044) | −0.031 | 0.559 | 0.330 |
| 8 mm | 0.289 | −0.082 (p=.094) | **+0.022** | **0.622** | 0.351 |

**The over-shrinkage is cured.** Error at the strongest voxels goes −0.303 to +0.022, and `g` there
goes 0.250 to 0.622 against a truth of 0.622. Whole-map rmse is lowest at 4 mm, so 4 mm minimises
error while 8 mm centres the peaks, and there is no single best radius -- the same shape as the
`coverage_radius` trade-off but an order of magnitude larger.

It also moves the crossing later. At three images the named voxel has already lost (+0.023) while
4 mm still wins (−0.013); at eleven images the deficit halves, +0.087 to +0.050. And `g_marginal`
improves monotonically with the radius, 0.269 to 0.205 at one image.

So the radius is not a free parameter but a function of the corpus: small where the coordinate
studies are few and their peaks trustworthy, larger where the silence limb would otherwise swamp
the reports. The quantity that decides it is the silence-to-report ratio, which is measurable
before fitting.

## Calibrating the adaptive report radius

The report extent is not a free parameter and not a constant either: which setting wins depends on
how badly the reports are outnumbered, and that spans two orders of magnitude across real
collections. Measured as the median silences per report over the voxels some study named:

| collection | ratio |
|---|---|
| NIDM pain, 9 tables | 5 |
| NeuroStore pain, 25 tables | 18 |
| NeuroStore pain, 200 tables | 140 |
| NeuroStore pain, 1,443 tables | 604 |

Sweeping the radius at each, as rmse for the named voxel against 4 mm:

| ratio | 1 image | 3 images | error at top, 1 image |
|---|---|---|---|
| 5 | named by 0.006 | -- | named by 0.008 |
| 18 | named by 0.006 | named by 0.003 | **4 mm by 0.022** |
| 140 | tie | **4 mm by 0.017** | **4 mm by 0.123** |
| 604 | **4 mm by 0.018** | **4 mm by 0.035** | **4 mm by 0.181** |

**The rmse crossing is bracketed between 18 and 140**, whose geometric midpoint is 50.2, and
`ADAPTIVE_REPORT_RATIO` is 50. Centring at the strongest voxels crosses earlier -- 4 mm is better
centred from 18 upward, and at 140 it takes `|g|` at the top from 0.271 to 0.433 against a truth of
0.640 -- so a user who cares about the peaks more than the whole map should lower the threshold.

Shipped as `report_radius="adaptive"`, which measures the ratio from the tables and picks, records
the choice on `report_radius_`, and logs it.

### Two harness defects this uncovered

The first mid-scale sweep returned **byte-identical results for 0 mm and 4 mm**, which is
impossible. It had caught the file between adding `report_radius` as a parameter of
`_indicator_entries` and updating the call site to pass it, so both arms silently ran at the named
voxel. Identical numbers across a swept parameter are the cheapest possible signal that the
parameter is not connected, and worth checking for deliberately.

Then `pkill -f hcp_sdm_vs_cbes`, meant to clear a contaminated run, matched its own shell's command
line -- the pattern appears in it -- and killed the replacement runs along with it.

### The adaptive rule across all four datasets, including where it misses

| collection | ratio | adaptive picks | better option | cost of the pick |
|---|---|---|---|---|
| NIDM pain, 9 tables | 5 | named voxel | named voxel | — |
| HCP MOTOR_LH / EMOTION_FACES | 10 | named voxel | **4 mm** | 0.02 of magnitude |
| working-memory studyset | thin | named voxel | tie | — |
| NeuroStore pain wall, 1,443 tables | 604 | 4 mm | 4 mm | — |

HCP, named voxel against 4 mm:

| task | r | AUC | magnitude | `g_marginal` magnitude |
|---|---|---|---|---|
| MOTOR_LH | +0.772 / +0.774 | 0.939 / 0.941 | 0.62 / **0.64** | 0.42 / 0.44 |
| EMOTION_FACES | +0.882 / +0.883 | 0.984 / 0.983 | 0.68 / **0.70** | 0.28 / 0.29 |

4 mm is consistently a little better on HCP, across both tasks and both maps, in the direction the
mechanism predicts: HCP under-recovers and a wider report adds positive bias. But it is 0.02 of
magnitude against a shortfall of 0.28, so it does not touch the actual problem, which is that
every HCP study genuinely has the effect and any fitted prevalence below 1 is error.

The threshold stays at 50. Moving it to catch HCP would chase 0.02 on a bed of synthetic studies
at the cost of the bracketing evidence from four real ratios, and the win it protects at the other
end is an order of magnitude larger. Recorded as a known miss rather than presented as
four-for-four.

### A third knob that was defined but never connected

`hcp_sdm_vs_cbes.py` reads `SKIP_SDM`, and the `SDM=0` switch added to match the other beds was
never wired to the call site, so four arms of a CBES-only comparison each ran the full SDM-PSI
pipeline. The symptom was duration -- one arm still going after nineteen minutes -- with no error.
That is the third instance in this session of a parameter that existed and did nothing: the
others were `report_radius` unplumbed through `_indicator_entries`, caught by 0 mm and 4 mm
returning byte-identical results, and this same import leaving `DEFAULT_REPORT_RADIUS_MM`
undefined on the default path only, because the ternary short-circuits for every value except the
default. **Sweeping a parameter and checking the results actually differ is the cheapest guard
available, and none of these three would have been caught by a test suite or a linter.**

## Prior art: does this estimator already exist?

Searched PubMed. The premise is not new; the estimator is a different treatment of it.

**MetaNSUE** (Albajes-Eizagirre, Solanes, Radua, Stat Methods Med Res 2018,
doi:10.1177/0962280218811349) states the exact premise: a study that reports "not significant"
without an effect size cannot be dropped (biased) and cannot be entered as zero (also biased).
Indexed under "interval censoring" and "multiple imputation". Reaches neuroimaging as SDM-PSI.

The difference is what happens to the bound. MetaNSUE draws plausible values inside it and pools
with Rubin's rules. CBES writes it into the likelihood as a censored term and maximises by EM.
Nothing is drawn. Consequences: no proposal distribution to misspecify, and the observed
information falls out directly, which is where the CBES standard errors come from. The PR
docstring's "Non-reporting is censoring, not imputation" survives contact with the literature.

Other ingredients, all previously published:

  * images + coordinates in one model -- ES-SDM, Radua et al. Eur Psychiatry 2011,
    doi:10.1016/j.eurpsy.2011.04.001
  * voxelwise effect size from sparse peaks -- Salimi-Khorshidi, Nichols, Smith, Woolrich,
    IEEE TMI 2011, doi:10.1109/TMI.2011.2122341 (GPR; interpolates from peak *height*, which
    CBES deliberately does not read)
  * filling in around a peak -- anisotropic kernels, Radua et al. Front Psychiatry 2014,
    doi:10.3389/fpsyt.2014.00013 (same problem as report_radius, solved by spatial correlation
    rather than an inequality)
  * random-effects CBMA using height, and the subsample-vs-high-powered-group validation design
    we are using -- Bossier et al. Front Neurosci 2018, doi:10.3389/fnins.2017.00745
  * point-process CBMA -- Montagna, Wager, Barrett, Johnson, Nichols, Biometrics 2018,
    doi:10.1111/biom.12713

Outside neuroimaging the machinery is routine and older: Tobit regression, survival
right-censoring, astronomical upper limits, environmental non-detects. PubMed does not index
those fields.

What appears new is the combination:

  1. censored likelihood instead of imputation
  2. silence as a spatial pattern read at every voxel, with the three-zone indicator
     (named / ambiguous ring / silent), rather than one per-study "this contrast was null"
  3. zero inflation as a separate parameter, so pi and mu are distinct reportable quantities
     (ES-SDM and MetaNSUE estimate one pooled effect)
  4. the declared correction threshold as the censoring bound, rather than inferring it from
     the smallest reported peak

Caveat, recorded so it is not forgotten: this cannot prove absence. English-language PubMed,
keyword search, and PubMed ANDs every term so the first two attempts returned zero. A
censored-likelihood CBMA in a statistics journal without neuroimaging keywords would not appear.
Before write-up: Google Scholar forward citations of MetaNSUE, and arXiv stat.ME.

Positioning for the PR and any paper: CBES is the censored-likelihood, voxelwise, zero-inflated
version of the idea MetaNSUE introduced. Cite Albajes-Eizagirre 2018 and Radua 2011 as direct
antecedents. Do not present the premise as novel.

## The spatial null does not fix the small-collection familywise rate

Ran the global-null bed with NULL=spatial-images (Fourier phase randomisation, autocorrelation
preserved) against the shipped permute-images, 40 sims, 200 permutations, mask eroded by 2.

  arm                       permute-images    spatial-images
  12 studies, 2 images          0.150             0.150
  20 studies, 2 images          0.075             0.075
  12 studies, 1 image           0.100             0.100

Identical to three decimals on every arm, with uncorrected rates at 0.050 throughout. Binomial se
on the familywise rate is 0.034, so this cannot resolve a difference smaller than about 0.07, but
there is no hint of one.

Conclusion: the residual inflation at 12 studies is not the permutation destroying spatial
structure. spatial-images stays available -- it is the right null to reach for if the question is
ever about spatial specificity -- but it does not become the default, since the randomiser costs
8-10x for no change in calibration. #82 needs a different mechanism.

## Censored likelihood against multiple imputation, and two bugs it found

`experiments/censoring_versus_imputation.py`. Scalar bed, truth known exactly, 25,600 fits per
arm: true mu 0.5, tau 0.15, 2 image studies, 20 coordinate studies, cut 0.691 in g units. Four
arms plus a brute-force grid MLE as an oracle, so a disagreement can be attributed rather than
argued about. The censored arm is not a reimplementation -- it calls CBES._fit_chunk on
hand-built arrays.

  arm                       bias      sd    rmse   CI width  coverage
  grid MLE (oracle)      -0.0112  0.0999  0.1006     0.3047     0.954
  censored (CBES)        -0.0112  0.0999  0.1006     0.3047     0.954
  imputation (MetaNSUE)  -0.0115  0.1001  0.1008     0.2721     0.910
  drop silent            +0.2364  0.0745  0.2478         --        --
  zero fill              -0.2705  0.0750  0.2807         --        --

The point estimates tie. A 20-rep pilot showed the censored arm 3.7% better on rmse; at 400 reps
that is 0.2%, i.e. nothing. Both methods are Monte Carlo and closed-form evaluations of the same
integral, so a tie is what theory predicts, and the pilot's margin was noise.

What survives is the interval: the imputation arm's is 11% too narrow and covers 0.910 against
nominal 0.95. Part of that is finite M -- Rubin's (1 + 1/M) correction is asymptotic and a
Barnard-Rubin df adjustment would recover some of it. So the defensible claim is narrower than
the docstring's "censoring, not imputation": the censored likelihood gets the interval right with
no tuning, where imputation needs M large enough and a df correction. Not "imputation is worse".

Both beat the two estimators MetaNSUE's abstract warns about by a wide margin: +0.24 and -0.27
bias against -0.011.

### Two shipped bugs the oracle exposed

**The EM retired at its start value.** Where a study reported, log(1 - P(silent)) is convex in mu,
so at a voxel with enough reporters the total observed curvature is positive and `-score /
curvature` points uphill. The `curvature < 0` guard then set the step to 0, `settled` fired, and
the voxel kept its start value with an se taken at a non-stationary point. 2.2% of voxels, mean
|error| 0.34 at each -- exactly the gap between the estimator's sd of 0.095 and the MLE's 0.080.
Elsewhere the EM matched the grid to 0.0006, so the EM was never wrong, only stalled. Fixed with a
Fisher-information denominator, non-negative by construction, same fixed point. A scan finds 48
stalling configurations; they all sit at small start values, where the reported limb's convexity
outweighs everything else.

**selection_model="none" profiled out a prevalence it never estimated.** It passed
identified=None, which skips the guard that stops the Schur complement subtracting for an unfitted
parameter. Information understated about eightfold, nominal-95% coverage 0.998. Fixed by always
passing the array: all-False without a mixture, which is the truth. The default zero-inflated path
was never affected.

After both, CBES matches the oracle to four decimals on estimate and interval alike.

### Three of my own errors on the way, all in the bed

  * started the EM at zero, where the same convexity refuses the first step. The shipped code
    starts at the pooled estimate (`start=fit["g"]`). A bed that does not copy the caller's start
    is measuring a different algorithm.
  * passed inverse-variance weights where `_accumulate` supplies unit weights ("Every image
    contributes weight 1"). The point estimate was unaffected -- weights scale the whole score, so
    the root does not move -- but the information was multiplied by the weight.
  * with those weights the bed read se/sd 0.745, and #71 records the old scalar bed reading 0.74.
    I thought I had found #71. I had not: both older scalar beds (censored_se.py,
    observed_information_se.py) already use unit weights. Coincidence, withdrawn.

## How a CBES impute() would differ from SDM-PSI's

SDM-PSI's algorithm, from Albajes-Eizagirre, Solanes, Vieta & Radua, NeuroImage 2019,
doi:10.1016/j.neuroimage.2018.10.077: "a) multiple imputation of study images; b) imputation of
subject images; and c) subject-based permutation test to control the familywise error rate."

Three differences, all of which have to be stated rather than glossed:

  1. Two imputation levels against one. SDM imputes study maps and then synthetic subject images,
     which is what lets it run Freedman-Lane permutation with TFCE. A CBES impute() gives level
     (a) only: downstream meta-analytic tools, not subject-level permutation.
  2. SDM conditions on reported peak heights and recreates the map with the anisotropic kernel.
     CBES reads only the +-1 indicator, so its bounds are |g| < c at a silence and |g| >= c at a
     report, with no point value anywhere. Imputed CBES maps are necessarily more diffuse. That is
     the price of refusing heights, made visible rather than hidden.
  3. CBES fits pi, so a completed value must come from the mixture: with probability 1 - pi the
     study has no effect at this voxel and the draw belongs to the null component. SDM imputes
     from a single distribution. This is a structural difference and probably an improvement -- a
     completed dataset can represent "this study genuinely had nothing here".

Do not smooth after imputing. A smoothed draw can leave its own interval, and the marginal
variance shrinks, so the completed values stop satisfying the censoring they were built from. SDM
avoids this because its kernel is part of the imputation model, not a post-hoc filter. The
construction that works is a Gaussian copula: draw a smooth GRF with the target covariance, map
through Phi to uniforms, then through each voxel's truncated-normal inverse CDF. Exact on the
marginal truncation, approximate on the correlation, and about as cheap as drawing the field.
Gibbs under a GMRF prior is exact and much slower.

Noted for #82: SDM-PSI's own validation found its FWER control "might be too conservative". Ours
is liberal at small collections, so their diagnosis will not transfer.

## The Fisher fix is inert on real data, and the earlier "slight weakening" was mine

Clean A/B on published pain, 8 splits, same seeds, the pre-fix code run from a git worktree at
41f9213 through PYTHONPATH:

                              before        after
  shipped CBES g, r           +0.583       +0.575
  mean err                    +0.072       +0.073
  err at top                  +0.013       +0.014
  rmse                         0.236        0.240

Everything else -- images only, silence switched off, g_marginal -- is bit-identical. So the
Fisher fallback moves pain by 0.004 rmse and nothing else, which matches the stall rate measured
at whole-brain scale: 18 to 46 (voxel, iteration) pairs out of 3.7 million on synthetic fits with
12 to 30 studies, i.e. at most about 0.1% of voxels.

I earlier read the post-fix numbers as weaker than rr_pain.log's. That comparison was against a
differently-configured run, not against the same code. The worktree A/B is the right comparison
and it says the fix is inert.

## _observed_information is exactly right, and max_iter=25 is doing something else

`experiments/zero_inflated_oracle.py`. The first version compared the fit against the argmax of a
2-D grid over (mu, pi) and produced an oracle *worse* than the estimator -- bias +0.19 against
-0.01, sd 0.48 against 0.40. That is the bed's fault. The mixture likelihood has a near-flat
ridge, so a global argmax wanders while an EM started at the pooled estimate settles nearby.
"The MLE" is not a well-defined target here the way it is without the mixture. Rebuilt to ask two
questions that need no unique maximum.

**Is the interval right?** Reported se against the exact Schur standard error from a
finite-difference 2-D Hessian of the analytic log-likelihood, evaluated at the estimator's own
(mu, pi): ratio median 1.0000, 10th-90th percentile 1.0000 to 1.0000, at every max_iter tried.
_observed_information is arithmetically correct. Whatever se/sd is in a given configuration is the
likelihood's, not the code's. Worth noting that in this bed se/sd is 0.78, below 1, where the
Notes say 1.2 to 2.2 -- so that figure is configuration-dependent and not a property of the
estimator.

**Is the fit at a maximum?** No. At max_iter=25, 75% of voxels sit more than 0.01 log-likelihood
below the best a coarse grid can find, median gap 0.096. At 100 it is 21%, at 400 it is 3.6%.

Converging changes the estimator consistently across four configurations (20 tables throughout):

  images  true pi | rmse(mu) @25  @400 | pi hat @25   @400
       6      1.0 |       0.083  0.085 |      0.970  0.999
       6      0.6 |       0.218  0.302 |      0.650  0.666
       2      1.0 |       0.210  0.294 |      0.847  0.998
       2      0.6 |       0.482  0.739 |      0.554  0.611

Early stopping shrinks mu toward the images-only pool and leaves the prevalence short. The
shrinkage is buying real accuracy on mu -- up to 0.26 of rmse -- and costing the prevalence up to
0.15. Three consequences, none of them currently documented:

  * the estimator does not compute the censored MLE at its default, and says nothing about it
  * the reported se is the exact curvature at a point that is not a maximum, so it is the right
    se for the wrong estimator
  * "prevalence is ordinal, not a fraction" may be an artefact of stopping at 25: at 400 the
    fitted prevalence reads 0.999, 0.666, 0.998, 0.611 against truths 1.0, 0.6, 1.0, 0.6

A full sweep over max_iter at four regimes is running, plus pain at 400. Nothing changes until
those land, and max_iter must not be tuned to a collection.

## CORRECTION: converging does not fix the prevalence, it overshoots

The 20-rep pilot above read a fitted prevalence of 0.611 at true 0.6 with max_iter=400, and I
wrote that the "prevalence is ordinal, not a fraction" caveat might be an artefact of stopping at
25. At 200 reps it reads 0.723. The pilot was noise. Full sweep, 6,400 fits per cell, 2 image
studies and 20 tables:

  max_iter | true pi 1.0: pi hat  rmse(mu)  short% | true pi 0.6: pi hat  rmse(mu)  short%
        10 |               0.709     0.183   0.943 |               0.541     0.410   0.885
        25 |               0.854     0.205   0.852 |               0.576     0.459   0.763
        50 |               0.955     0.243   0.412 |               0.615     0.540   0.504
       100 |               0.996     0.281   0.111 |               0.657     0.654   0.208
       200 |               0.998     0.287   0.022 |               0.699     0.695   0.073
       400 |               0.998     0.289   0.007 |               0.723     0.707   0.035

pi hat climbs monotonically with the iteration cap and simply crosses the truth on the way past.
It lands on 1.0 only because 1.0 is a boundary. The ordinality caveat stands.

What survives: max_iter is an undeclared shrinkage parameter, monotone in both directions, and 25
is a point on a bias-against-ordering trade rather than a convergence criterion. On published
pain, 25 to 400 gives mean error 0.073 to 0.069 and error at the top 0.014 to 0.004 (better
level), against r 0.575 to 0.560, AUC 0.886 to 0.882 and rmse 0.240 to 0.244 (worse ordering).
Documented in the max_iter entry; the default is unchanged, because nothing here says a different
value is better overall.

Third time this session a 20-rep pilot has pointed the wrong way (the others: the censored arm's
"3.7% better rmse" than imputation, which was 0.2% at 400 reps; and reading the post-fix pain run
as weaker when the comparison run was differently configured). Pilots size a run. They are not
findings.

## A docstring claim, falsified by my own measurement

The null_method entry said the family-wise rate is liberal in small collections because
permute-images destroys the image's spatial autocorrelation, so "the observed map can lay a
coherent blob of large values over a region where few studies were silent and a scattered map
cannot". It does destroy the autocorrelation. That is not the cause: spatial-images preserves it
and returns 0.150, 0.075 and 0.100, identical to permute-images on all three arms. Replaced with
the measurement and an explicit statement that the mechanism is unknown.

Worth noting the direction was never right either. A rougher null field has more effective resels,
so its maximum is larger, which would make the test conservative rather than liberal.

## The two imputation claims, now measured instead of asserted

`experiments/copula_imputation_check.py`, 1-D, 4096 points, 93% of voxels silent, target field
autocorrelation 0.802 at lag 4.

  construction                      bound violations     sd   autocorr
  independent draw                            0.0000  0.267      0.143
  independent draw, then smoothed             0.0432  0.118      0.848
  Gaussian copula                             0.0000  0.258      0.696

Claim 1 holds: smoothing after imputation puts 4.3% of voxels outside the interval they were
drawn from and collapses the sd by 56%. A completed map built that way no longer satisfies the
censoring it was built from, which is the whole point of building it.

Claim 2 holds, with a number attached: the copula violates nothing and recovers 0.696 of the
target's 0.802, i.e. 87%, not 100%. A Gaussian copula preserves rank correlation, and the
truncation's nonlinear marginal transform attenuates the Pearson correlation on top of that. Quote
87%, not "preserves the correlation".

1-D caveat, per the protocol: a kernel of a given FWHM couples far more neighbours in a volume, so
read the retained fraction as a direction rather than a number until it is checked in 3-D.

## se/sd is prevalence-dependent, and the Notes have the sign wrong for the useful case

12,800 fits per cell against a known truth, 20 tables throughout.

  images  true pi  max_iter   se/sd
       2      1.0        25    1.25
       2      1.0       400    1.84
       2      0.6        25    0.80
       6      1.0        25    1.52
       6      0.6        25    0.80

The estimator's Notes say "se is conservative" with se/sd from 1.2 to 2.2. That range is real and
it is what a prevalence of 1 gives -- which is where it was measured. At a prevalence of 0.6 it is
0.80, i.e. the interval is *anti-conservative*, and 0.6 is the regime the estimator exists for:
studies genuinely differing in whether they carry the effect. Where every study has the effect
there is no absence to find, and the Warnings already say the correction can only do harm there.

This is not an arithmetic error. The zero-inflated oracle verified _observed_information against
an exact finite-difference 2-D Hessian at ratio 1.0000. The observed information at the fitted
point simply does not capture how far mu wanders along the near-flat ridge between simulations
when the mixture is real.

The Notes entry has to change: se/sd is not a property of the estimator, it is a property of the
configuration, and the direction flips.

## The coordinate channel's information about (mu, pi) is rank 1

`proofs/mixture_identification.py`, 10 sympy claims plus a calibration check against #73.

A reporting indicator is a Bernoulli whose silence probability is the mixture integrated over
(-c, c): S = pi*s_a(mu) + (1-pi)*s_0. Its two scores are

    dS/dmu = pi * u      with u = d s_a / d mu
    dS/dpi = v           with v = s_a - s_0

so every block of its information matrix is built from the same two scalars, and

    I_mumu * I_pipi - I_mupi^2 = 0     identically

whatever the number of coordinate studies. The Schur complement from indicators alone is exactly
zero. **Coordinate tables cannot separate mu from pi at all.** Every bit of that separation comes
from the images, whose two scores are not proportional: the mu score carries a factor (g - mu)
that vanishes at g = mu while the pi score does not.

That is the algebraic version of the wall result (#81). A corpus of 1,443 tables adds no ability
to tell a large effect in few studies from a small effect in many, because it adds rank-1 matrices
to a rank-1 matrix.

### Unless the studies differ, and then Lagrange's identity says by how much

For heterogeneous studies the determinant is not zero but

    det I = pi^2 * sum_{j<k} w_j w_k (u_j v_k - u_k v_j)^2

so the entire identifying power of the coordinate channel is the spread of u/v across studies, and
nothing else. Evaluated on the two designs #73 measured:

  sample sizes 12-120 at z=3.09    3.0864e+00
  thresholds 2.3-4.5 at n=20       8.8768e-02
  no spread at all                 0

35x in favour of spreading sample sizes, which is the measured direction: sample-size spread
lifted the bounded fraction from 0.42 to 0.69 and threshold spread changed nothing. The ratios u/v
spread about equally in the two designs (1.23 against 0.93), so the factor is not about the ratio.
Spreading the threshold moves both mixture components together and leaves v = s_a - s_0 small,
which is what the cross-product is built from. That was the verbal explanation in #73; it is now a
number, and the proof carries the check so it cannot drift.

### What this decides

  * #81, algebraically: the wall of tables is rank 1 and cannot be anything else.
  * #73 explained rather than recorded, with a design rule: to identify the magnitude, collect
    studies of differing *sample size*. Differing thresholds are worth about 3% as much.
  * The se/sd problem: the Wald se inverts I_mumu - I_mupi^2 / I_pipi, which goes to zero exactly
    at the ridge. With two images the non-degenerate part of the information is rank 2 from two
    observations, so the Schur complement is near zero and the Wald interval is unstable in both
    directions -- which is what 0.80 at the default and 2.55 at convergence are.
  * A shippable diagnostic follows: the estimator already computes all three blocks in
    _observed_information, so it can report how close I_mupi^2 is to I_mumu * I_pipi and refuse
    the Wald interval where the ridge is flat. interval="profile" is already implemented for
    exactly that case.

## A rougher null cannot make a max-statistic test liberal

`proofs/roughness_cannot_be_liberal.py`, 6 claims. Closes a direction rather than re-measuring it.

Discrete version, airtight: for n independent standard normals, P(max > u) = 1 - Phi(u)^n, whose
derivative in n is -Phi^n log Phi > 0. Differentiating Phi(q)^n = 1 - alpha implicitly gives
dq/dn = -log Phi(q) Phi(q) / (n phi(q)) > 0, so the critical value rises with the number of
locations.

Field version: at a high threshold the exceedance probability is the expected Euler
characteristic, sum over d of R_d rho_d(u), linear in each resel count. In a volume rho_3 factors
as k_3 (u-1)(u+1) exp(-u^2/2) with k_3 > 0, so every EC density is positive above u = 1 and a
max-statistic threshold sits far above 1. Resel counts go as V / f^3, whose derivative in the FWHM
is -3V/f^4 < 0, so roughening adds resels in every dimension at once.

Chain: rougher null -> more resels -> larger expected EC at any threshold -> larger null maximum
-> higher critical value -> fewer rejections. **A permutation that roughens the map is
conservative.** The observed statistic maps are 2.5x smoother than the permuted ones, which pushes
the test the wrong way to explain a family-wise rate of 0.150, and spatial-images already returned
identical rates. The mechanism is doubly dead: wrong empirically and wrong in sign.

#82 still needs a mechanism. What is now excluded: the boundary artefact (#79, fixed by eroding),
a degenerate null (197-200 distinct maxima of 200, cv 0.07-0.11), the Pareto tail (identical with
TAIL=0), the null being too narrow (percentile 0.467 against 0.462), the p-value arithmetic (it is
the exact (1 + #{null >= obs}) / (1 + B), and the histogram is used only for cluster forming), and
now roughness in both directions.

Still open and worth testing: whether the inflation is there at all with no coordinate tables (the
2-study 2-image arm reads 0.080 at n=100, se 0.022, so 1.4 se above nominal -- suggestive, not
decisive), and whether data-dependent early stopping breaks exchangeability, which the MAXITER=400
arm will answer.

## Why the naive count out-ranks the fitted prevalence (#51)

`proofs/count_versus_fitted_prevalence.py`, 4 claims. Not a failure of the model -- a statement
about what each quantity estimates.

The expected report rate collapses to a single product:

    P(report) = 1 - s_0 - pi * (s_a(mu) - s_0)

affine in pi with slope s_0 - s_a and no curvature, and monotone in mu through the same term. So
two voxels with the same pi(s_a - s_0) give the same expected count however different their
prevalences.

Hold mu roughly fixed across voxels and the product is affine in pi, so the count is a strictly
monotone transform of the prevalence: its rank correlation with the truth is 1 up to binomial
noise. Nothing can rank better, and an estimator that also fits mu must rank worse in finite
samples, because it spends information separating two things that did not need separating. Where
mu does vary, the count confounds them and the fit should win.

So #51 is a measurement of the corpus, not a verdict on the estimator. It also makes a prediction
worth testing: the count's ordering advantage should shrink as the spread of g across voxels
grows. Untested.

Two defects in my own proof, caught before committing. One claim duplicated another's expression.
The level-set claim was a tautology -- I subtracted exactly the terms I had added, so it reduced
to zero by construction -- and its stated conclusion was wrong as well, since the rate depends on
pi(s_a - s_0) rather than pi*s_a. Suspect the proof as readily as the test.

## WITHDRAWN before shipping: identified_share does not flag a bad Wald interval

The algebra gives Schur = I_mumu - I_mupi^2 / I_pipi, which goes to zero at the ridge. I inferred
a diagnostic from that -- report identified_share = Schur / I_mumu, the fraction of information
about mu surviving the prevalence, and distrust the Wald interval where it is small -- wrote it
into the estimator, and then simulated it, in that order.

The simulation says the opposite. 150 reps, 32 voxels, prevalence 0.6:

  identified_share   voxels  coverage  mean se  sd of err
           0.0-0.2     1299     0.943   0.5213     0.3791
           0.2-0.4     2229     0.935   0.2354     0.3373
           0.4-0.6      473     0.924   0.2865     0.4865
           0.6-0.8      201     0.955   0.7445     0.9719
           0.8-1.0      219     0.630   0.5860     0.8669

A flat ridge makes the standard error *large*, so coverage there is fine or conservative. The
interval fails at high identified_share, where the prevalence costs nothing and mu wanders anyway.
Reverted, uncommitted, nothing shipped.

The algebra was not wrong. The step from "Schur goes to zero" to "therefore small Schur flags a
bad interval" was an extra inference I never proved, and it happens to be false: a wide interval
is not a wrong one. That step is exactly what the algebra-first rule is for, and the simulation
stage caught it before real data or the repository.

What would be needed instead: a proof about Var(mu_hat) rather than about the information at the
fitted point -- the failure at high identified_share is between-collection wandering that the
observed information does not see, which is a statement about the sampling distribution of the
estimator, not about the curvature of one likelihood.

## Sample-size spread buys prevalence accuracy, and the cross-product predicts how much

`experiments/does_sample_size_spread_identify_prevalence.py`, 9,600 fits per arm, kill condition
stated in the docstring before running (under 10% and the inference from the theorem fails).

  arm                  cross-prod   pi bias   pi rmse   mu rmse
  no spread at all         0.0000   -0.0252    0.2438    0.4586
  sample sizes 12-40       0.8137   -0.0139    0.2325    0.4781
  thresholds 2.3-4.5       0.8451   -0.0073    0.2363    0.4371
  sample sizes 12-120     37.2951   +0.0303    0.2075    0.4405
  both spread             36.4866   +0.0302    0.2161    0.4301

14.9% of prevalence rmse between no spread and the widest, so the kill condition is not met and
the inference from the theorem stands. A 40-rep pilot read 14.8%, which is the first time this
session a pilot and the full run have agreed.

The better result is the correspondence. Accuracy tracks the cross-product, not the nominal
design: the 12-to-40 sample-size arm and the threshold arm have nearly equal cross-products (0.81
against 0.85) and nearly equal prevalence rmse (0.2325 against 0.2363), despite being different
kinds of heterogeneity. That equivalence is not obvious and the identity predicted it. The
no-spread arm's cross-product is exactly 0.0000, as the algebra requires.

mu rmse does not improve (0.44 to 0.48 across all arms, no ordering). The spread buys prevalence,
not magnitude -- which fits, since mu is carried by the images either way.

### What this means for requiring a sample size

Two reasons, of unequal weight.

  1. The silent drop, which is the bigger one. `_all_sample_sizes` raises only when *every*
     sample size is missing; otherwise `series[series.notna()]` removes those studies and
     `study_ids = list(sample_sizes.index)` makes that series the roster. A study without a
     sample size never enters the model, so its silence -- "the single most informative
     observation about how common the effect is", by the estimator's own docstring -- is
     discarded with no warning. Studies missing the field are not a random subset.
  2. Identification, worth about 15% of prevalence rmse here. Real, confirmed, and modest.

Also worth documenting either way: a collection of uniformly-sized studies at one threshold gets
*no* separating information from its tables, exactly. That belongs in the Warnings beside "the
correction can make g worse than doing nothing".

## Two more proofs: the boundary test for pi = 1, and what max_iter really is

### proofs/boundary_test_for_full_prevalence.py (5 claims)

The worst failure mode -- nothing tells a user whether their collection has pi < 1, and the
correction helps in one regime and hurts in the other -- is a hypothesis test that has not been
written down because the usual recipe does not apply. pi = 1 is on the boundary of the parameter
space, so 2 log Lambda is not chi^2_1 but Chernoff's half-and-half mixture of chi^2_0 and chi^2_1
(Self and Liang 1987, doi:10.1080/01621459.1987.10478472; Chernoff 1954,
doi:10.1214/aoms/1177728725). Marked in the proof as asserted rather than derived, because it is a
distributional result about where the unconstrained maximum falls, not an algebraic identity.

Practical content, a factor of two: the level-0.05 critical value is **2.7055**, not 3.8415. Using
the naive chi^2_1 cut runs the test at half its nominal level.

Verified: the score in pi is (f_1 - f_0)/f, which at the boundary is 1 - f_0/f_1 and needs no fit;
the local expansion about the boundary is eps*r - eps^2 r^2 / 2 with r = f_0/f_1 - 1, so the
statistic is quadratic in the departure and its scale is E[r^2] = I_pipi.

**Where the power comes from, and this is the useful part.** The rank-1 theorem says the
coordinate channel cannot separate mu from pi. It does *not* say the channel is silent about pi
alone: one indicator contributes I_pipi = n_t v^2 / D with v = s_a - s_0, increasing in the number
of tables. So unlike the split, this test does draw power from a wall of tables. Power vanishes
only when v -> 0, i.e. when a study with the effect is no more likely to report than one without.

Next per the order of operations: simulate its size and power before writing any of it into the
estimator.

### proofs/early_stopping_is_shrinkage.py (5 claims)

An EM map converges linearly with rate rho = the fraction of information the missing data hides
(Dempster, Laird and Rubin 1977). So stopping at t iterations does not return an approximate MLE
with unspecified error; it returns

    theta_t = (1 - rho^t) theta* + rho^t theta_0

an exact convex combination of the MLE and the starting value. Equating that to the optimum of a
quadratic penalty towards theta_0 gives the penalty the cap is silently imposing:

    lambda_eff = I rho^t / (1 - rho^t)

and lambda_eff is increasing in the missing information. **So a fixed cap shrinks hardest at the
voxels where the censoring hides the most -- which is where the coordinate channel is carrying the
estimate and the images are carrying least.** A prior nobody chose, different at every voxel, and
invisible in the output.

What this does not settle: whether a stated lambda would be better. The measurements say the
shrinkage is buying real accuracy on mu (rmse 0.183 at 10 iterations against 0.289 at 400), so a
replacement has to beat it in simulation before it is worth writing.

## The familywise inflation is probably not real: 0.055 at 200 simulations

Every "0.150 at 12 studies" figure in this file came from a 40-simulation run. The binomial
standard error of a rate estimate at n=40 near 0.05 is 0.034. A 40-simulation run cannot
distinguish 0.05 from 0.15. It was never evidence.

Re-run at 200 simulations, same erosion, same permutation count, same seeds:

  arm                         sims    se    uncorrected  familywise
  12 studies, 2 images          40  0.034       0.0501       0.150
  12 studies, 2 images         200  0.015       0.0498       0.055
  2 studies, 2 images           100  0.022       0.0507       0.080
  6 studies, 6 images           100  0.022       0.0487       0.030

The 200-simulation estimate is 0.055 against a nominal 0.05, i.e. nominal. The no-table arms are
0.080 and 0.030, both within noise of 0.05 at their precision.

What this invalidates, and it is a lot:

  * the family-wise numbers I wrote into the null_method docstring -- 0.150, 0.075, 0.100 --
    stated as evidence that the rate is "too liberal in small collections". They are not evidence.
  * the spatial-null comparison's *numbers*. "Identical on all three arms" is true and useless at
    se 0.034. Its conclusion survives only because roughness_cannot_be_liberal.py is an
    independent argument.
  * the whole #82 investigation, including the arms still running to test mechanisms for an
    effect that may not exist.

What survives: #79, where eroding the boundary moved the rate 0.375 to 0.075. That is nine
standard errors even at n=40.

The lesson is the same one three times over today, now at its most expensive: a 40-replicate run
is a pilot that sizes the real run. Reporting one as a finding, and worse, writing one into a
shipped docstring, is how a session spends hours chasing an artefact. Simulation counts need a
stated standard error before the number is written down anywhere.

## The boundary test: calibrated above 20 tables, and starved of power by small studies

`experiments/boundary_test_size_and_power.py`, 3000 replicates per cell, kill conditions stated
before running. Two images, coordinate studies at n = 20.

  tables   atom at 0   size @2.71   naive/boundary
       5       0.828       0.0343             0.49
      20       0.796       0.0390             0.51
      80       0.786       0.0447             0.56
     400       0.665       0.0523             0.62

Kill condition 1 fires at 5 tables: 0.0343 is outside 0.05 +- 3 se. From 20 tables up the derived
critical value is usable, and the factor-of-two claim holds. The pattern is asymptotic arrival,
and the atom explains it -- Chernoff's mixture needs the boundary atom at 1/2 and it is 0.83 at 5
tables, falling towards 1/2 as tables grow. The score at the boundary is a sum of 1 - f_0/f_1, and
f_0/f_1 is a right-skewed likelihood ratio with mean 1, so its sample mean falls below 1 more than
half the time and the maximum sits on the boundary too often.

Power, though, is the problem:

  prevalence   5 tables   20    80    400
         0.9      0.060  0.095 0.110  0.103
         0.8      0.083  0.118 0.140  0.161
         0.6      0.096  0.137 0.158  0.178
         0.4      0.081  0.094 0.112  0.108

Never above 0.18, and not monotone in the departure: it peaks near 0.6-0.8 and falls at 0.4,
because with two images and few active studies both images are likely null, mu-hat collapses
toward zero and the mixture explains the data at any pi.

### Why, and the answer is sample size, not table count

At pi = 1 the per-table information about pi is v^2 / (s_a (1 - s_a)) with v = s_a - s_0. The null
component's silence probability is free of the sample size entirely -- the cut in effect-size
units is z/sqrt(n) and the null sd is 1/sqrt(n), so their ratio is z whatever n is. The whole gap
is carried by s_a, whose sd tends to tau rather than to zero while the cut tends to zero. So:

  n        12      20      40      80     120     200
  I_pipi  0.124   0.309   1.079   3.873   8.284  21.237

171x from n = 12 to n = 200. Spread helps too, by Jensen on a convex function, but only +32% over
12-120 at fixed mean -- against two orders of magnitude for the level. That is the opposite
balance from the mu/pi split, where the level is irrelevant (a uniform collection gives exactly
zero at any size) and only the spread matters, because separation is a determinant of differences
while this is a sum.

Stated plainly: **you can only tell that some studies lack the effect if studies that have it
would reliably report it.** A small study is silent either way, so its silence says nothing about
which component it belongs to. And 400 tables at n = 20 carry the same information about pi as 15
at n = 120 -- twenty-seven small studies are worth one large one for this question.

Running: power against coordinate-study sample size at n = 20, 40, 80, 160, and against image
count at 6 and 20 images.

## CORRECTION to the retraction: the family-wise rate is unverified, not nominal

I retracted the inflation on the strength of a 200-simulation run reading 0.055. That run was
launched **before** the Fisher-denominator fix, so it does not describe the shipped code. A
100-simulation run after the fix reads 0.100. Three estimates at 12 studies and two images:

  sims   code        familywise   se
    40   pre-fix          0.150   0.034
   200   pre-fix          0.055   0.016
   100   post-fix         0.100   0.022

The 40 and 200 are nested in seeds and consistent (6 rejections in the first 40, 11 in 200, so the
early clustering was chance). The 100-simulation post-fix run is a different revision. 0.055 and
0.100 differ by 1.7 standard errors -- not a separation -- and neither separates from 0.05.

So the honest state is **unverified**, not "nominal" and not "inflated". I swung from one
confident conclusion to the opposite one within an hour, on runs that could not support either.
The docstring now says unverified below about 20 studies. A 500-simulation post-fix run is going.

Also worth recording: the no-tables arms were measured post-fix and read 0.080 (2 studies, 2
images) and 0.030 (6 studies, 6 images) at 100 simulations each, se 0.022. Both within noise of
nominal, and they do not isolate the coordinate channel as the cause of anything, because there is
not yet an established effect to attribute.

## Boundary power against sample size: the prediction holds

`experiments/boundary_test_size_and_power.py`, 2000 replicates per cell, 2 images.

  n of coordinate studies    size@20 tab   power(pi=0.8)   power(pi=0.6)   [80 tables]
                       20         0.0415           0.141           0.169
                       40         0.0525           0.206           0.276
                       80         0.0440           0.259           0.350
                      160         0.0370           0.338           0.443

Size stays calibrated across the range (0.036 to 0.055 against 0.05). Power at prevalence 0.6 goes
from 0.169 to 0.443 as the coordinate studies grow from 20 to 160 subjects, which is what
I_pipi = v^2 / D predicted: 0.309 to about 14 per table over that range. Table count buys much
less (0.128 to 0.169 going from 20 to 80 tables at n = 20).

Still not a usable test at realistic sizes -- 0.44 power to detect prevalence 0.6 needs 80 tables
of 160 subjects each -- but the mechanism is confirmed and the design implication is clear: for
this question, recruit larger studies rather than more of them.

## The interval is fine where the data identify mu, and bad where they do not

`experiments/does_the_profile_interval_cover.py`, 12,800 fits, prevalence 0.6, 2 images, 20
tables. A 40-rep pilot read profile 0.958 against Wald 0.928 and looked like a fix. It was a
selection confound: the profile bound is finite at only 52% of voxels, and those are selected for
being informative, so scoring profile there against Wald everywhere flatters profile.

Scored on the same voxels:

  max_iter   Wald (all)   Wald (bounded)   profile   width ratio
        25        0.922            0.969     0.964          0.98
       400        0.932            0.967     0.964          1.02

Wald and profile are indistinguishable once compared like with like, so kill condition 1 is met in
substance: the profile interval is not the fix, and swapping the default would buy nothing.

The useful result is what falls out of the two Wald columns. Backing out the unbounded half at
max_iter 25: coverage is about **0.87 where the profile bound is infinite and 0.97 where it is
finite**. The anti-conservatism that se/sd 0.80 was pointing at is entirely concentrated in the
voxels where the data do not identify mu -- and the finiteness of the profile bound already flags
exactly those, at no new cost, in code that already ships.

That is the diagnostic I went looking for earlier and got wrong. identified_share failed because a
flat ridge makes the interval *wide*, not wrong. The profile bound's finiteness is a different
quantity -- whether the data reject pi = 0 at all -- and it is the one that tracks coverage.

Actionable and needing no code: read ``g`` and its interval where ``g_lower``/``g_upper`` are
finite. Elsewhere the point estimate is still the MLE but the interval under-covers by about 8
points.

## The regime question is answerable, with about twenty images

`experiments/boundary_test_size_and_power.py`, 2000-3000 replicates per cell. Power at prevalence
0.6, coordinate studies at n = 20:

  images    20 tables   400 tables
       2        0.128        0.178
       6        0.346        0.407
      20        0.624        0.749

Levers at prevalence 0.6, ranked: images 2 to 20 takes power 0.128 to 0.624; coordinate-study
sample size 20 to 160 takes 0.169 to 0.443 (at 80 tables); table count 20 to 400 takes 0.128 to
0.178. Images dominate, which is what the rank-1 theorem requires -- the separating information is
theirs and the indicator's is rank 1.

Size is conservative everywhere measured, 0.030 to 0.044 against a nominal 0.05, and never
liberal. The cause is the boundary atom: Chernoff's mixture needs 1/2 and the measured atom runs
0.66 to 0.83, falling as information grows. So the derived critical value 2.7055 is *safe* and
leaves power on the table; calibrating the atom empirically would recover some.

The practical statement for the worst failure mode: a collection with two images cannot be told
whether its prevalence is below 1 -- power 0.13 at a true 0.6. A collection with twenty can, at
0.62 to 0.75. That is a design answer, and it is the first one this failure mode has had.

## CORRECTION from the SDM-PSI paper itself: it uses a censored likelihood

I claimed, in this file and in the shipped docstring, that MetaNSUE and SDM-PSI resolve the bound
by multiple imputation while CBES writes it into the likelihood. Having read
doi:10.1016/j.neuroimage.2018.10.077, that is wrong. Section 3.2.2 states the likelihood
explicitly:

    L = prod_i [ Phi((y_upper,i - X_i b) / sqrt(v_upper,i + tau^2))
               - Phi((y_lower,i - X_i b) / sqrt(v_lower,i + tau^2)) ]

described as "not the likelihood of a specific effect size but the likelihood that the unreported
effect size lays within the two effect size bounds", citing Tobin (1958), Costafreda (2012) and
Schnedler, "Likelihood estimation for censored random vectors". That is the same Tobit
construction. The multiple imputation runs *after* it -- "SDM-PSI only uses MLE as a starting
point for the subsequent multiple imputation, avoiding the biases associated with single
imputation" -- to propagate uncertainty and to enable the subject-image permutation.

So the earlier head-to-head in censoring_versus_imputation.py measured CBES against a *pure*
multiple-imputation estimator that nobody actually ships. The measurement stands as a comparison
of two estimators; it does not describe SDM-PSI, and the docstring claim built on it was
retracted in NiMARE fd6d1fb.

The impute() design I proposed under #86 is also not novel: SDM-PSI already draws from a truncated
normal with the MLE and the bounds as parameters, and already imposes positive correlation between
adjacent voxels using the AES-SDM correlation templates -- the spatial-coherence problem I flagged
as the hard part. #86 should be recast as "match what SDM-PSI does", not "invent it".

### What actually differs, from the source

  * **Peak heights.** "The lower and upper effect-size bounds are obvious in a peak: both are the
    effect size of the peak." A reported peak enters as an exact observation at its reported
    height, i.e. with zero interval width, so the upward bias of a selected maximum enters
    directly. CBES discards heights; that costs information and avoids this.
  * **No prevalence.** Their beta is one effect size, "the same for all studies". So there is no
    g against g_marginal distinction and no rank-1 identification problem -- there is nothing to
    separate. CBES's ridge is the price of asking a question SDM-PSI does not ask.
  * **Declared shrinkage.** SDM-PSI damps influential studies on purpose, iteratively discarding
    the study that most increases the absolute MLE, "similar to a trimmed mean", and says this is
    "required for a correct control of the FWER". CBES shrinks too, through max_iter, but by
    accident. Theirs is a documented choice.
  * **Inference.** Permutation of imputed subject images against permutation of image values.

### Their FWER is more nuanced than "too conservative"

Table 1, 400 simulated meta-analyses per cell with Clopper-Pearson intervals: conservative in most
scenarios (0-4%), **but 15-23% for cluster-based statistics at high z-thresholds in small
meta-analyses of small studies**. So SDM-PSI has a liberal corner in exactly the regime I was
chasing in CBES: few studies, small studies, cluster statistics. Worth noting their evidence
standard -- 400 per cell -- against the 40 I drew conclusions from.

## Homogenised sample sizes on real pain data: g is unaffected

`HOMOGENISE=1` in validate_redesign.py declares the corpus mean sample size for every study while
the data keep their true sizes -- the scenario where an analyst has no per-study sample sizes and
defaults them. Published tables, 8 splits, same seeds.

  estimate                  r     rank r   AUC    mean err  err at top   rmse
  true sizes,      g   +0.575     +0.481  0.886     +0.073      +0.014  0.240
  homogenised,     g   +0.582     +0.482  0.885     +0.066      -0.011  0.228
  true sizes,      gm  +0.602     +0.485  0.887     -0.062      -0.282  0.184
  homogenised,     gm  +0.583     +0.476  0.879     -0.062      -0.296  0.181

The spread destroyed is real: the corpus runs 9 to 32 subjects, a 3.56x range with cv 0.38. So
this is not a vacuous test, and flattening it left g slightly *better* (rmse 0.240 to 0.228), not
worse. With 8 splits a 0.012 difference is not a result in itself; the result is the absence of
harm.

**This confirms the simulation rather than contradicting it, and corrects what I told the user.**
The spread bed already showed that sample-size spread buys prevalence rmse (0.244 to 0.208, 15%)
and leaves mu rmse alone (0.4586, 0.4781, 0.4405, 0.4371, 0.4301 across arms, no ordering). The
pain bed scores mu against a held-out reference and has no ground truth for prevalence, so "no
harm to g" is exactly the predicted outcome. I had said a defaulted constant "kills the mu/pi
split outright and misreports the regime test" -- true of the mechanism, overstated as a
consequence for g.

Also worth pricing: pain's 3.56x spread is much narrower than the simulation's 12-to-120. The
closest simulated arm, 12-to-40 at 3.3x, bought only 4.6% of prevalence rmse against 15% for the
wide one. So on a corpus like this the identification benefit of real sample sizes is a few
percent, not fifteen.

### Where the requirement now rests

  1. **The silent drop.** A study with no sample size never enters the roster, so its silence is
     discarded with no warning, and studies missing the field are not a random subset. Unchanged
     and still the strongest reason.
  2. **The threshold conversion.** z to g depends on n, so a wrong n biases the censoring bound.
  3. **Prevalence accuracy.** Real, measured at 15% for a tenfold spread and about 5% for a
     threefold one.
  4. **Not g.** Magnitude is carried by the images and survives a 3.56x spread being flattened.

So: require it, drop-with-warning rather than silently, and do not justify the requirement by
claiming g needs it.

## External support for coverage_radius, from Radua et al. 2014

Full text via PMC3919071, doi:10.3389/fpsyt.2014.00013. They recreate effect-size maps from peak
information at a grid of anisotropy degrees and FWHMs, and score against the effect-size map
computed from the raw statistical parametric map.

Their optimal isotropic kernel is **40-45 mm FWHM**, which they note is "substantially larger than
in previous validations" of 20-25 mm. A 40-45 mm FWHM puts the half-maximum at 20-22 mm from the
peak, which is where DEFAULT_COVERAGE_RADIUS_MM = 20 sits. That parameter was chosen here as
"roughly the extent a paper's peak stands in for", with nothing behind it. It now has an external
number that agrees.

Two further things worth having:

  * They reached the adaptivity conclusion from the other direction. "Optimal FWHM might vary
    largely depending on the specific data meta-analyzed", and their remedy is full anisotropy
    *because* it removes the FWHM dependence entirely (their Eq. 3 at alpha = 1). CBES's adaptive
    report_radius is the same admission with a cruder remedy -- pick one of two radii from a
    corpus statistic, rather than remove the dependence.
  * The correlation templates are published and free, covering all 26 neighbours of every voxel
    for grey matter, white matter, CSF and FA. That is exactly the covariance the Gaussian copula
    in #86 needs, so it does not have to be estimated.

**Caveat, and it matters for how much weight the agreement carries.** Their validation recreates
maps from the peaks of the *same* statistical map it then compares against: 120 IXI subjects split
six ways, each split thresholded at p < 0.001 with a 10-voxel minimum, peaks extracted, map
recreated, compared to that split's own effect-size map. The truth is conditioned on the noise
that produced the peaks, which is the thing PROTOCOL.md forbids and the reason this project splits
studies or holds out subjects. So their 40-45 mm optimum is plausibly tuned to how well a kernel
reproduces that map's own noise, and the agreement with 20 mm should be read as encouraging rather
than as confirmation.

Their cluster counts are worth recording for the reporting bed: mean 22 +- 28 clusters per map,
median 16 +- 14, at p < 0.001 with a 10-voxel minimum extent -- the same scheme reporting.py uses.

## Bossier et al. 2018: the matching validation design explicitly excluded censoring

Full text via PMC5778144, doi:10.3389/fnins.2017.00745. Their design is the one PROTOCOL.md asks
for: resample a large dataset (IMAGEN) into pseudo-studies as the test condition, score against an
independent high-powered group analysis as the evaluation condition, across 10, 12, 14, 16, 18,
20, 30 and 35 studies.

Their Discussion states the gap directly:

  "bias due to missing data if peak effect sizes for some studies are not reported (Wager et al.;
  Costafreda). Seed based-mapping, uses imputations to solve this latter missing data problem. As
  we did not have any missing data in our simulations, we did not evaluate the influence of these
  missing data on the performance of the various CBMA methods."

So the benchmark whose design this project copies was run in a regime with **no censoring at
all**, and its favourable result for effect-size random-effects CBMA says nothing about what
happens when studies fail to report. There is no existing benchmark for the thing CBES does. That
is worth stating plainly when this is written up -- not as a claim to novelty, but because it
means the comparison everyone would reach for does not exist yet.

It also names Costafreda as a source on the missing-data bias, alongside Wager. Second
independent pointer to Costafreda 2012 (doi:10.1016/j.jneumeth.2012.07.016) in an hour, after
SDM-PSI cited it for the censored-likelihood MLE. That is now the paper to get.

## The mu bias at prevalence 1 is a boundary effect, not an O(1/k) term

I claimed the bias looked like a correctable O(1/k) mixture bias, on two points: +0.068 at two
images and +0.022 at six, a factor of three for a factor of three. The sd fell by the same factor
and I did not check it. Swept properly, 300 reps per cell, images 2 to 16:

  images   max_iter      bias        sd   bias/sd
       2         25   +0.0721    0.1965     0.367
       2        400   +0.0645    0.2875     0.224
       4         25   +0.0411    0.1104     0.372
       4        400   +0.0311    0.1331     0.234
       8         25   +0.0231    0.0749     0.308
       8        400   +0.0175    0.0778     0.225
      16         25   +0.0151    0.0561     0.269
      16        400   +0.0122    0.0572     0.213

At convergence bias/sd is flat at about 0.22 across an eightfold change in image count. An O(1/k)
bias would give bias/sd falling as 1/sqrt(k), a factor of 2.83 from 2 to 16 images: 0.224 would
have to become 0.079. It reads 0.213.

So the bias is a fixed fraction of the standard deviation, which is a boundary signature. At a
true prevalence of 1 the mixture can always explain a small magnitude by lowering pi instead of
lowering mu, and that asymmetry floors mu-hat at a constant fraction of its own scale.

Two consequences:

  * Not correctable by subtraction. The correction would need the sd, and removing a bias of this
    kind re-inflates the variance.
  * Not worth chasing. A bias of 0.22 sd contributes 0.22^2 / (1 + 0.22^2), about 5%, of mean
    squared error. The lever on it is the same as the lever on the variance: more images.

At max_iter 25 the ratio declines mildly, 0.37 to 0.27, so early stopping adds a component that
does fall with the image count -- consistent with the shrinkage pulling towards the images-only
pool, whose own bias behaves differently.

Item B of the proof queue is therefore closed without a proof: there is no O(1/k) term to derive.
Third time today an inference drawn from a correct observation has failed its measurement, and the
pattern is the same each time -- a ratio looked constant, or a difference looked large, and I did
not check what the denominator was doing.

# ==================== A second estimator: the marginal effect ====================

Design document supplied 16 September. Working through it in the order the rules now require:
algebra, then simulation, then real data. New NiMARE branch
`claude/marginal-effect-combined-estimator`, taken from main at bc361b1.

## The document's own numbers check out

Verified before building on it, since a design document deserves the same suspicion as a testbed.
All four checkable claims reproduce exactly: the peak ratio 8.331 at u = 0.5 with M = 27;
information eigenvalues [0, 10.593] for identical cutoffs against [0.259, 7.110] for 0.35 and 0.8;
the non-identified pair (active mean 0.500, prevalence 0.3181, marginal 0.1591) and (0.800,
0.1160, 0.0928); and the control-variate reductions of 23% and 59%.

## proofs/marginal_control_variate.py -- 8 claims

The image-corrected predictor m_hat = Ybar_I + lambda (fbar_C - fbar_I). Built from actual random
variables at n = 2, N = 3 and handed to sympy.stats, which independently confirms the mean is m
for any lambda and any predictor, and the variance is Var(Y - lambda f)/n + lambda^2 Var(f)/N
including its cross term. My first draft "verified" unbiasedness by subtracting mu_f from itself,
which proves nothing; that is now a real computation.

  lambda* = Cov(Y,f) / ((1 + n/N) Var(f))          a regression slope, shrunk
  V*/V_images = 1 - rho^2 / (1 + n/N)

**The cap is the useful part.** As N grows the ratio tends to 1 - rho^2, so a perfect predictor
with unlimited coordinate studies still cannot beat the images alone by more than the squared
correlation. There is no configuration in which coordinates substitute for images. At 8 images
and 100 tables: 23% reduction at rho = 0.5, 59% at rho = 0.8, against floors of 0.750 and 0.360.

## proofs/reporting_partition_and_variance.py -- 10 claims

Three defects, two of which correct claims of mine.

**The three reporting events are not a two-outcome experiment.** Report inside one radius, silence
beyond a larger one, annulus dropped. The likelihood uses P(|Y| >= c) and P(|Y| < c) as
complementary. They are not: 1 - (p_A + p_C) = p_B, and a correct conditional likelihood carries
1/(p_A + p_C). The error does not cancel between limbs -- it biases the score by
d/dmu log(p_A + p_C), which moves with mu. This is a real defect in the shipped model, not a
documentation problem.

**An image-only tau2 is not the active component's tau2.** Var(theta) = pi tau_a^2 +
pi(1-pi) mu_a^2, so substituting the total for the active part errs by
-(1-pi) tau_a^2 + pi(1-pi) mu_a^2, vanishing only at pi = 1. At mu_a = 0.5 and tau_a = 0.15 the
error is +0.036, +0.051, +0.047 at prevalence 0.8, 0.6, 0.4 -- two to three times tau_a^2 itself,
and positive, so the active spread is over-stated. `_pool` does exactly this substitution.

**No universal factor-of-two bound on the report limb.** The peak probability is
(1 - Phi(u)^M)/M, verified by differentiating rather than by symbolic integration, and the ratio
of exceedance to peak probability is 8.33 at u = 0.5, growing as M 2^(M-1)/(2^M - 1) as the
threshold vanishes. My docstring's "bounded by a factor of two" is false in general. At the
conventional z = 3.09 the ratio is 1.013, so it happens to hold where it is used -- over-general
rather than wrong in practice, and it should say so.

Also caught by the proof rather than by me: I first wrote that u -> 0 limit as M/2. It is
M 2^(M-1)/(2^M - 1), because Phi(0) is 1/2 and not 0. The claim failed and I fixed it.

## Corrections owed to the estimator's Notes

  * "spreading thresholds 2.3 to 4.5 changes nothing" -- false as stated. Threshold spread does
    restore rank; the document's counterexample gives a minimum eigenvalue of 0.259 where
    identical cutoffs give 0. What I measured was a small *effect on prevalence rmse* over one
    range, which is a different statement.
  * "the error is largely bounded by a factor of two" -- not universal.
  * my reading of Bossier. Their methods threshold at voxelwise FDR 0.05 and extract cluster
    maxima, so their test condition is censored summary data. The sentence I quoted concerns
    missing *peak effect sizes*, not absence of censoring. I over-read it and the note saying
    "no existing benchmark for what CBES does" is too strong.

## Real data kills the control variate on pain, and shows why

`experiments/control_variate_on_pain.py`. 21 NIDM pain studies, 267 published peaks, no cap. Eight
studies held out as the reference so the truth never touches a coordinate used in the fit, six as
the image cohort, seven coordinate-only. 40 splits.

**Two correlations, and only one is the one in the formula.** The variance formula's rho is
Corr(Y, f) *across studies at a fixed voxel*, because that is the covariance the estimator
averages over. I first measured the within-study, across-voxel spatial correlation -- whether the
kernel recreates a study's map shape -- which is a different and much larger number.

  fwhm   rho_spatial   rho_between   ratio at 6+7   floor
    10         0.173         0.115          0.993   0.987
    20         0.266         0.094          0.995   0.991
    30         0.287         0.062          0.998   0.996
    45         0.241        -0.002          1.000   1.000

Stated kill condition was rho below 0.3. The right rho is **0.06 to 0.12**, an order of magnitude
short, and the achievable variance reduction is under 1%. Measured performance agrees: rmse 0.1707
against 0.1710 for images alone, a 0.2% difference that is noise at 40 splits.

**The mechanism, and it is the useful part.** A kernel over reported peaks does capture the shape
of a study's own map -- rho_spatial 0.29 at 30 mm, which is roughly what ES-SDM's recreation
validation would predict. It carries almost nothing about how *this* study differs from *another*
at a given voxel, which is the only thing a control variate can exploit. Where a paper reports is
largely common across a domain; the between-study variation is in magnitudes the table does not
record.

Note the widths disagree about which is best: 30 mm maximises the spatial correlation and 10 mm
the between-study one, and at 45 mm the between-study correlation is zero. A wider kernel washes
out exactly the study-to-study differences the estimator needs.

So the design document's 23% and 59% are attainable only with a predictor far better than a peak
kernel. Its own text says as much -- "algebraic targets, not measured Neurostore prediction
quality" -- and this is that measurement. To make this estimator worth running, the predictor has
to track between-study variation, which means it has to use the reported heights or study
covariates, not just peak locations.

Also recorded: the cohort shift runs 1.11 of its own standard errors, so on pain the image and
coordinate cohorts are not detectably non-exchangeable at this size. That is a necessary check,
not a reassuring one -- 21 studies cannot resolve a shift that matters.

### What survives

The estimator itself is correct and validated: 8,000-replication simulation matched the proven
bias, variance, ratio and floor, and the two failure modes reproduce. It is the *predictor* that
fails on this corpus. That is worth having, because the implementation is now a measuring device
for any future predictor: hand it one, read rho_between, and the ceiling follows without fitting
anything.

## What a table of coordinates plus statistics can and cannot fit (2026-09-16)

Asked directly: does CBES fit on coordinates and peak statistics alone?

No, and the failure is by design in one place and structural in another.

  * **Sample size is required, no default.** `_focus_table` and `_all_sample_sizes` both raise.
    A silence is evidence only against the precision of the study that stayed silent.
  * **Threshold is not required**: default two-tailed p < .001 (z = 3.2905), lowered per study by
    `clamp_threshold=True` to that study's `min|z|`.
  * **Coordinates alone fit nothing in either estimator.** `g` comes from the images (#36);
    `optimal_coefficient` raises below two image studies because lambda is a covariance. Both
    estimators are for a mixed corpus, not a coordinate-only one.

### Candidate defect in clamp_threshold: an order statistic that floats with signal

`min_j |z_kj|` is the last order statistic of a study's reported heights. It falls as the peak
count rises, and the peak count rises with the study's signal. So with the threshold assumed
rather than declared, a peak-rich study's assumed cutoff is pulled near its true cut while a
peak-poor study keeps the default, which may sit above its true cut.

Nothing is capped, yet the effective threshold still floats with per-study signal -- the same
mechanism the never-cap rule describes, reached through the order statistic instead of a cap.

Status: **algebra written, and it reverses my headline.** `proofs/clamped_threshold_order_statistic.py`,
9 claims. Every silence scores minus the inverse Mills ratio, so it always pushes mu down; raising
the assumed cutoff attenuates that push and the probit information alike; and attenuating the
silent terms raises the root of the estimating equation (`dmu*/dalpha = B/S' < 0`). The unclamped
default overstates the cutoff by the full `c_default - c_k`, the clamp by only `1/(M*theta)`, so
**the clamp reduces the upward bias rather than creating one**. What survives of the worry is not a
bias in mu-hat but a *differential weighting* of silences by per-study signal, which nothing there
shows to be benign -- it is the quantity that matters for a spread across studies, not for the
direction of mu-hat.

Calibration: the Gaussian tail hazard at the default cut is 3.554, so at the pain corpus's 12.7
peaks per study the overshoot is 0.022 z -- small, consistent with the shipped "bit-identical where
no table contradicts the assumption". Regime it cannot speak to: a cluster-forming corpus, where
the reported peaks are maxima of clusters that also passed an extent criterion. My first numeric
check of that regime had no extent selection in it and measured nothing; it is removed rather than
left reading as evidence, and the direction there is unknown. Calibration target is clean, because the clamp is documented bit-identical
where no table contradicts the assumption: measure `clamp_threshold=True` against `False` on a
corpus that carries statistics.

Independent of the sign: declare the threshold per study from the paper's stated correction
scheme, and let the clamp catch only contradictions. That removes the order statistic from the
path.
