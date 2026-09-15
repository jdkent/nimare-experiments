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
