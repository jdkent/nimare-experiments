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

## Waiting on a background job: never `pgrep -f <script-name>`

A waiter of the form

    until ! pgrep -f "myjob.py" >/dev/null; do sleep 20; done

matches **its own command line**, because the waiting shell's argv contains the pattern. It
therefore waits on itself and never exits. This cost three queued experiments a full half hour of
silence in which everything looked like it was running -- and the failure is invisible, because
"still waiting" and "waiting forever" print the same thing.

Use one of these instead:

- **Just run the jobs in sequence** in one background script. No polling, no pattern, nothing to
  get wrong. This is almost always the right answer.
- If a wait is genuinely needed, wait on a **PID**, not a pattern: record `$!` when launching and
  use `wait "$pid"` (same shell) or `kill -0 "$pid"` in a loop.
- If a pattern is unavoidable, exclude self: `pgrep -f "myjob.py" | grep -qv "^$$\$"`, and check
  it actually distinguishes before trusting it.

The general rule this belongs to: **a wait that cannot fail loudly will eventually fail
silently.** Before arming any waiter, ask what it would print if the thing it waits for never
started -- and if the answer is "nothing", fix the waiter.

## Every bed declares and checks its statistic convention

A reported statistic means something specific, and the estimator acts on that meaning. CBES reads
a reported value as a **t on n - 1 degrees of freedom** (a reported z is treated as a
p-value-preserving image of one) and maps it back before converting to an effect size. In the far
tail, where every reported peak lives, that map is strongly expansive: at ``z = 5.33`` with
``n = 30`` it turns ``d = 0.98`` into ``g = 1.25``.

A bed that generates ``(truth + noise / sqrt(n)) * sqrt(n)`` is producing a **normal statistic
with known variance** -- exactly ``d * sqrt(n)`` -- and feeding it in inflates every recovered
effect size by about a third. Measured at the same foci with the same reporting pipeline:

```
known-variance z declared as Z    bias at the foci  +0.630
proper t declared as T            bias at the foci  +0.296
```

This has now happened twice in this program: once in the field simulator (fixed), then again in a
new bed, where it produced a large, stable, reproducible bias that I spent most of a session
attributing to the estimator -- inventing a "pooling step" contribution and then a "conversion
convexity" one before checking the input. It is the recurring failure mode, not a subtle one.

So:

- Build statistics with ``reporting.study_t_field`` and report them as ``"T"``.
- Call ``reporting.assert_statistic_convention(null_field, n, "T")`` once, on a **null** field,
  before reading any estimate out of the bed. Print that it passed.
- Where a bed also supplies images, derive them from the same ``t`` through the estimator's own
  ``peak_stat_to_hedges_g``, so the image arm and the coordinate arm share one convention. Using
  ``t / sqrt(n)`` is Cohen's d and is high by the Hedges factor -- 2.6% at ``n = 30``, which is
  enough to show up as a residual bias in an arm meant to be the unbiased reference.
- Check the *tail*, not the variance. A t on 29 df has variance 1.074 against a normal's 1.000,
  a 7% gap that sampling noise hides; the first version of the check duly passed a field whose t
  had collapsed to a z. ``P(|T_29| > 3) = 0.0055`` against ``P(|Z| > 3) = 0.0027`` is a factor of
  two and cannot be missed.
- Do not make a variance field smooth by smoothing chi-square draws: that averages independent
  variates, the denominator goes nearly constant, and the t collapses into a z -- the very
  confusion being avoided. Use a probability integral transform of a smooth Gaussian field, which
  keeps the marginal exact.

## Before explaining a measurement, run the one-minute check that it is wrong

Three of today's errors are the same error. A number arrived, it was surprising, I went looking
for a mechanism, and I found one in the estimator -- because the estimator was what I was
studying. In each case a check costing about a minute would have shown the number was not what I
thought it was.

| the surprising number | the mechanism I proposed | what it actually was | the check I skipped |
| --- | --- | --- | --- |
| 30% of replications unfittable | a degenerate-null defect | 18% of studies reported anything: the bed was in the wrong regime | print the reporting fraction |
| pooling adds +0.324 of bias | kernel up-weights the larger excursions | my stage-2 baseline used `z/sqrt(n)`; the pooling step adds +0.001 | compare the estimator's own conversion on one number |
| a further +0.323 from "conversion convexity" | Jensen on a convex transform of a selected maximum | my bed reported a known-variance z where the estimator expects a t (Jensen is only +0.043) | check the null variance or tail of the generated statistic |

The rule that would have caught all three, stated so it can be executed rather than merely
agreed with:

> **When a measurement surprises you, write down the cheapest check that would show the
> measurement is wrong, and run it before proposing any mechanism.** If you cannot name such a
> check, that is the finding -- the bed is not instrumented well enough to be believed.

"Suspect the test before the theory" was already in my notes and did not stop me, because it is a
disposition rather than a step. Naming the check is the step. All three checks above take about a
minute and each would have saved an hour of confident wrong explanation.

A corollary on where to look: the mechanism you find first will be in whatever you are currently
studying, because that is where your attention is. That is a property of attention, not evidence
about the system. The input path -- what the bed generates, what convention it declares, what
regime it lands in -- is the part nobody is studying and therefore the part that goes unchecked.

## What worked: write the prediction down first

The counterweight to the above, and it earned its keep today. Four predictions were recorded
before their tests reported:

- the bias tracks image *share* rather than donor count, with a numeric range -- **confirmed**
  (+0.176 measured against +0.12 to +0.20 predicted)
- the prevalence/magnitude split needs a spread of study power -- **confirmed**, and by a
  diagnostic the prediction had not named (fitted `g` swinging 53% across a prevalence sweep
  under a fixed roster and going flat once sample sizes varied)
- the `g_marginal` cancellation depends on reporting density -- **confirmed**
- a couple of images would distort the relative map into a U-shaped ratio spread -- **falsified**,
  monotone 1.65x to 1.28x to 1.06x

The falsification was the most immediately useful of the four, because it retired a mechanism I
would otherwise have carried into the write-up as an explanation for an old puzzle. A prediction
written down cannot be quietly reshaped to fit what arrives; one held in mind can, and will.

### Postscript: I broke the pgrep rule within the hour of writing it

Having written the rule above, I then armed a monitor with
`until ! pgrep -f "is_g_a_scale_error.py" | grep -qv "^$$\$"`. The `grep -v` was an attempt to
exclude the waiter's own PID, which does not work: `pgrep -f` also matches the tail and subshell
processes in the pipeline, each with a different PID from `$$`.

Two conclusions, both duller and more useful than "be careful":

- The rule needs to be **never write a `-f` pattern that could match your own command line**, not
  "exclude yourself when you do". Exclusion is one more thing to get wrong and it silently
  degrades to waiting forever, which is the failure mode that looks like success.
- **And it is not only `pgrep`.** I later ran `pkill -f interval_coverage.py` inside a compound
  command whose own text contained that string, and killed the shell running it -- exit 144, the
  whole command lost including a patch and a file write I had queued behind it. Third occurrence
  in one session, in a third tool. Kill by **PID**:
  `PID=$(ps -eo pid,args | grep "[i]nterval_coverage.py" | awk '{print $1}')`, where the bracket
  trick keeps `grep`'s own line out, then `kill "$PID"`.
- A rule written in a document does not fire at the moment of writing code. The thing that would
  have fired is a *habit*: launch with `run_in_background`, capture the PID, wait on `kill -0`.
  Where a habit is available, prefer it to a rule.

## Check what a knob does before sweeping it

I identified `coverage_radius` as the lever on how much of the brain CBES estimates, swept it from
8 mm to 45 mm, and got a covered-share of **0.089 at every radius, to three decimals**. Five
values of a parameter, one answer: the parameter does not do what I said it does.

It governs how far from a reported focus a study counts as having been *silent* rather than
uninformative -- which changes the estimates without changing where they exist. The extent of the
map is set by the *pooling kernel*, so the lever is `fwhm`.

The tell was there in the output and nearly went past me, because the column I was watching (AUC)
*did* move: 0.796, 0.819, 0.850, 0.863, 0.817 across the radii. A knob that moves one column and
not another is not a broken knob; it is a knob doing something other than what you think. The
check is one line:

> **Before reading a sweep, confirm the parameter moved the thing you are attributing to it.**
> If a column is identical across every setting, stop -- either the parameter is not wired in, or
> it does not control that column, and both change the conclusion.

This is the same family as the statistic-convention bug and the `study-min` degeneracy: in all
three the harness was doing something reasonable that I had mislabelled, and the label was what I
reasoned from. Cheap to catch, expensive to miss.


## Two concurrent runs of one script must not share a scratch directory

`interval_coverage.py` wrote its simulated images to a fixed `/tmp/claude-0/cov_imgs` with names
built from `(seed, studies, images, tau, k)`, so two concurrent runs generated *identical
filenames* and could overwrite each other's files mid-fit. Fixed by deriving the path from the
pid: `f"…/cov_imgs_{os.getpid()}"`, which costs nothing and removes the class.

**That fix is hygiene, not a diagnosis.** It was prompted by a run that silently dropped eight of
nineteen arms from its output table -- no error, exit code 0, and the dropped ones scattered
through the middle of the sequence rather than truncated from the end. I wrote the collision up as
the cause and then checked: the dropped set included `12 studies, 0 images, per-study`, which
writes no images at all, and all three representative dropped configurations run correctly when
called standalone. So the collision cannot explain it and **I did not find out what did.**

Recording it unresolved rather than leaving a tidy wrong answer, because the wrong answer is the
more expensive artefact: the next person to see arms vanish would fix the scratch directory and
believe the matter closed. What is actually known is that the arms are individually fine, the
comprehension that builds the list is fine, and the loop printed 11 of 19 with no diagnostic.
The scientific question that run was asking was re-answered by a separate single-purpose script
instead, which is the right move once a harness is under suspicion -- do not debug the harness on
the critical path of a result.
