"""Read the completed coverage log and answer the three questions it was run to answer.

Takes the log path as an argument so it can be pointed at any run rather than hard-coding one.

1. **What does the documented interval actually cover, and how wide is it?** The log reports
   `cov(z)` at 1.96 and `cov(t)` at each fit's own `dof`, plus the median `dof`, so the width of
   the documented interval is recoverable as `t(median dof) * mean se / truth`.

2. **Does the bias-to-width model still predict coverage under the t interval?** It predicted
   `cov(z)` across sixteen arms to a mean absolute error of 0.034. If it also predicts `cov(t)`
   then nothing new is happening and the t is just a wider interval; if it does not, the per-
   replication variation in `dof` matters and the median is not a sufficient summary.

3. **Would referring the se to the censoring roster instead of n_eff fix the calibration?** dof
   counts only the studies that reported near the voxel while the se draws on the silent ones
   too; at a well-covered voxel the roster was 10 against an n_eff of 4.76. Projecting the roster
   critical value onto each arm's own bias and spread says whether option (a) of that decision is
   a fix or just a smaller error.
"""
import re
import sys
import numpy as np
from scipy.stats import norm, t as student_t

PATH = sys.argv[1] if len(sys.argv) > 1 else (
    "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad/"
    "coverage_both_intervals.log")
TRUTH = 0.800
ROW = re.compile(
    r"^(?P<label>(?:12|24) studies,.*?)\s{2,}"
    r"(?P<mean>[-+]?\d+\.\d+)\s+(?P<bias>[-+]\d+\.\d+)\s+(?P<se>\d+\.\d+)\s+"
    r"(?P<sd>\d+\.\d+)\s+(?P<sesd>\d+\.\d+)\s+(?P<covz>\d+\.\d+)\s+(?P<covt>\d+\.\d+)\s+"
    r"(?P<dof>\d+\.\d+)\s+(?P<half>\d+\.\d+)")


def model(half, bias, sd):
    return norm.cdf((half - bias) / sd) - norm.cdf((-half - bias) / sd)


rows = []
for line in open(PATH):
    m = ROW.match(line.strip())
    if m:
        rows.append({k: (v if k == "label" else float(v)) for k, v in m.groupdict().items()})
if not rows:
    raise SystemExit(f"no rows parsed from {PATH}")

print(f"{len(rows)} arms parsed; truth {TRUTH}\n")
print("--- 1. the documented interval: what it covers and what it costs ---")
print(f"{'arm':36s} {'bias':>7s} {'se/sd':>6s} {'dof':>5s} {'cov(z)':>6s} {'cov(t)':>6s} "
      f"{'width_t':>8s}")
for r in rows:
    crit = float(student_t.ppf(0.975, max(r["dof"], 1.0)))
    print(f"{r['label']:36s} {r['bias']:+7.3f} {r['sesd']:6.2f} {r['dof']:5.1f} "
          f"{r['covz']:6.2f} {r['covt']:6.2f} {crit * r['se'] / TRUTH:8.2f}")

print("\n--- 2. does the bias-to-width model still predict, under the t? ---")
err_z, err_t = [], []
for r in rows:
    crit = float(student_t.ppf(0.975, max(r["dof"], 1.0)))
    err_z.append(model(1.96 * r["se"], r["bias"], r["sd"]) - r["covz"])
    err_t.append(model(crit * r["se"], r["bias"], r["sd"]) - r["covt"])
for name, e in (("cov(z) at 1.96", err_z), ("cov(t) at median dof", err_t)):
    e = np.abs(np.array(e))
    print(f"  {name:24s} mean abs error {e.mean():.3f}   max {e.max():.3f}")
print("  A larger error for cov(t) means the median dof does not summarise the per-replication")
print("  variation in the critical value, and the arm needs its own dof distribution.")

print("\n--- 3. would the censoring roster be a better reference? ---")
print("  Roster taken as the study count less the studies that reach no focus near the voxel;")
print("  measured at 10 of 12 on one fit, so approximated here as 0.83 * studies.")
print(f"{'arm':36s} {'dof(n_eff)':>10s} {'cov(t)':>6s} {'dof(roster)':>11s} "
      f"{'cov(roster)':>11s} {'width':>6s}")
for r in rows:
    n_studies = 24 if r["label"].startswith("24") else 12
    roster_dof = max(0.83 * n_studies - 1.0, 1.0)
    crit = float(student_t.ppf(0.975, roster_dof))
    cov = model(crit * r["se"], r["bias"], r["sd"])
    print(f"{r['label']:36s} {r['dof']:10.1f} {r['covt']:6.2f} {roster_dof:11.1f} "
          f"{cov:11.2f} {crit * r['se'] / TRUTH:6.2f}")
print("\n  A column of values near 0.95 would make the roster a fix. Values spread across the")
print("  range mean the dof is not what is wrong, and the bias is.")
