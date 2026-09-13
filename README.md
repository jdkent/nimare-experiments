# NiMARE experimentation: effect-size coordinate-based meta-analysis

The design record and experimental scratch behind `nimare.meta.cbma.CBES`, kept out of the
NiMARE tree so the pull request carries the estimator and its tests rather than a research diary.

Nothing here is a NiMARE API, nothing is packaged, and the scripts are scratch quality: they hard
code paths, assume a populated `~/.nimare` cache, and several of them exist only to record an
approach that did not work.

## Layout

    notes/effect_size_cbma.md   the design note: what was tried, what the numbers were, why the
                                surviving choices survived. Sections are referenced by number.
    notes/validate_cbes.py      the two reproducible validations quoted in the note: recovery
                                against the 21 NIDM pain images, and false positive rates.
    notes/peak_height_deconvolution.py
                                whether a reported peak height can be deconvolved back to an
                                effect size. It cannot, in the usual underpowered regime.
    experiments/                one-off harnesses: SDM-PSI comparisons, weighting diagnostics,
                                calibration attempts, null approximations, profiling.
    experiments/logs/           their raw output, kept because several conclusions rest on runs
                                that take hours and nobody should have to repeat them to read
                                the number.

## The short version

The estimator localises effect magnitude well and cannot pin its absolute scale from
coordinates alone.

* Against the 21 NIDM pain images, from their peaks only, it reaches r = 0.80 over the whole
  mask and r = 0.81 where the reference has signal. SDM-PSI on the identical coordinates gets
  0.76 and 0.68.
* Magnitude is the weak point: coordinate-only it runs about 2x high, where SDM-PSI runs about
  2.2x low. Neither is calibrated, and they fail in opposite directions.
* Five images out of twenty-one fix it: r = 0.91 where the signal is, magnitude ratio 0.98,
  against 0.76/0.76 for those five images alone. SDM given the same five reaches 0.76 and 0.53.
* The scale is not recoverable from coordinates in principle, not merely in practice. Seven
  attempts are in section 15; the censored likelihood's value term is exactly invariant to it.

## Negative results worth not repeating

Section 15 (seven attempts at the absolute scale), section 16 (spatial similarity explains 0.26%
of magnitude variance), section 17 (a metric dividing by near-zero inflated every quoted
inflation figure), section 18 (an RFT regional censoring term: right in theory, worse in
practice), section 19 (coordinates look 29x under-weighted against images and correcting it
makes the estimator worse, because weight governs bias leakage rather than noise averaging).
