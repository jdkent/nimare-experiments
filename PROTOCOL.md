# Standing rules for every measurement in this repository

These are not suggestions. A result produced against any of them is void and has to be redone.

## Never cap the number of peaks

Do not take "the strongest N peaks per study", `.head(n)`, `[:n]`, or any other limit on how
many coordinates a study contributes. Not 10, not 40, not "a realistic number".

A cap fixes the count and lets the *effective threshold* float: the smallest value still in the
table is that study's Nth strongest peak, which is a function of how much signal the study had.
Every quantity that depends on the reporting threshold -- the censoring term, `prevalence`, the
inferred cut, the peak-height floor -- then moves with the signal, and any relationship measured
between the estimate and the truth is contaminated by that path. It invalidated the conditional
reference, the f(u) threshold fit, and the first cross-dataset run.

The count is not a free parameter. It is whatever survives the correction.

The one exception is a measurement whose subject *is* capping, and it must say so in its
docstring.

## Extract coordinates the way papers produce them

Use `experiments/reporting.py`. It does three things and no others:

  * **corrects for multiplicity** -- voxelwise FDR at q = 0.05, voxelwise family-wise error, or
    a p < 0.001 cluster-forming cut kept only where the cluster survives. Never a bare
    uncorrected height.
  * **keeps whole clusters**, so a study with nothing surviving reports nothing at all and
    drops out of the meta-analysis. That is a real outcome and the selection model needs to see
    it.
  * **separates peaks by 8 mm**, because adjacent-voxel maxima are not separate rows in a table.
    This is spatial de-duplication, not a count limit.

It returns the real height threshold alongside the peaks. Hand that to the estimator through a
metadata field and `threshold="reporting_threshold"` rather than letting it be inferred from the
smallest reported value.

Report all three schemes where the data allow it. Real papers do not agree on one, and the
estimator's magnitude turns out to depend on which was used.

## Keep the truth independent of the coordinates

Split the studies: one half is reduced to coordinates, the other half forms the reference. With
one map per study there is nothing to split within a study, so a same-map reference conditions
the comparison on the noise that produced the peaks. Held-out subjects are better still where
per-subject maps exist.
