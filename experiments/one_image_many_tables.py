"""The extreme configuration: a fixed wall of pain tables, images added one at a time.

Every earlier sweep here confounded two things. `when_do_coordinates_start_hurting.py` splits a
working half into image donors and table donors, so adding an image *removes* a table and the two
channels move together. The configuration that actually arises in practice is the opposite: the
literature is fixed and large, and what a meta-analyst gains over time is images.

So the coordinate channel is held constant -- the same NeuroStore pain studies at every image
count -- and images are added one by one. That decouples the channels, and it makes the crossing
condition sharp. With `m` in the thousands the coordinate channel's own variance `1/(m I_c)` goes
to zero, so `b^2 < 1/(m I_c) + 1/(k n)` collapses to `b^2 < 1/(k n)`: the crossing lands at
`k* = 1/(n b^2)`, set by the coordinate channel's bias alone, and **a large table wall cannot
rescue a biased one -- it only makes the bias harder to outvote.** A feasibility probe already
shows median `coordinate_share` climbing 0.303 -> 0.947 as the tables go 25 -> 1443, so at the top
end one image is arguing against a channel that carries 95% of the weight.

The reference is fixed within a draw rather than being the complement of the donors: ten of the 21
NIDM pain images are held out and pooled, and the donors are drawn one at a time from the other
eleven. That keeps the truth identical across image counts, so every comparison along a row is
paired, and it keeps the truth independent of the coordinates -- which `build_pain_table_corpus.py`
enforces by dropping the twelve NeuroStore studies whose peaks reproduce an NIDM table.

Two assumptions are unavoidable and load-bearing, both flagged in the corpus builder: only 0.5% of
the release's analyses report a sample size, so a coordinate study's `n` is assumed, and its
reporting threshold is assumed at p < 0.001 and then lowered per study by `clamp_threshold` to its
own smallest reported statistic. `NSPREAD=1` draws the sample sizes instead of fixing them, since
a spread of sample sizes is the only thing that identifies the magnitude.
"""

import os, sys, warnings

warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import logging

logging.getLogger("nimare").setLevel(logging.ERROR)
import numpy as np
import pandas as pd
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nimare.meta.cbma.effectsize import CBES
from nimare.studyset import Studyset
from load_pain import load_pain
from validate_redesign import g_and_var, masker, pooled, score

CORPUS = "/tmp/claude-0/paintables/pain_tables.parquet"
WORK = f"/tmp/claude-0/onemany_{os.getpid()}"
N_REFERENCE = int(os.environ.get("NREF", 10))
N_DRAWS = int(os.environ.get("NDRAWS", 3))
TABLE_COUNTS = [int(v) for v in os.environ.get("TABLES", "25,200,1443").split(",")]
IMAGE_COUNTS = [int(v) for v in os.environ.get("IMAGES", "1,2,3,5,8,11").split(",")]
#: Assumed sample size for a coordinate study, and whether to spread it.
ASSUMED_N = int(os.environ.get("ASSUMEDN", 20))
SPREAD_N = os.environ.get("NSPREAD", "0") == "1"
ASSUMED_Z = 3.0902
#: Whether to hand each peak its own reported statistic, converted to z, instead of the assumed
#: cut. It changes exactly one thing -- `clamp_threshold` lowers a study's assumed cut to its own
#: smallest reported statistic, and with every peak carrying the same assumed value there is
#: nothing to lower, so the clamp is inert and every study is read as having cut at p < 0.001.
#: 62% of the corpus's peaks report a t, so for most studies the real cut is knowable and lower,
#: and assuming the higher one reads a silence as a stronger bound than the paper justifies.
#: Peaks without a statistic keep the assumed cut.
USE_REPORTED_STATS = os.environ.get("STATS", "0") == "1"
#: Radius over which a report asserts its lower bound. The default named voxel leaves the report
#: limb at four studies against a thousand silences on this corpus; a few mm takes it to 80:1.
REPORT_RADII = [
    None if v in ("", "none") else float(v)
    for v in os.environ.get("REPORTR", "none").split(",")
]

mask_img = load_mni152_brain_mask(resolution=4)


def coordinate_study(sid, points, n, stats_z=None):
    meta = {"sample_sizes": [int(n)]}
    values = stats_z if stats_z is not None else np.full(len(points), ASSUMED_Z)
    return {
        "id": sid,
        "name": sid,
        "metadata": meta,
        "analyses": [
            {
                "id": sid,
                "name": "1",
                "metadata": meta,
                "points": [
                    {
                        "space": "MNI",
                        "coordinates": [float(c) for c in xyz],
                        "values": [{"kind": "Z", "value": float(v)}],
                    }
                    for xyz, v in zip(points, values)
                ],
            }
        ],
    }


def image_study(index, g, var, n):
    os.makedirs(WORK, exist_ok=True)
    gp, vp = f"{WORK}/g{index}.nii.gz", f"{WORK}/v{index}.nii.gz"
    if not os.path.exists(gp):
        nib.save(masker.inverse_transform(g), gp)
        nib.save(masker.inverse_transform(var), vp)
    meta = {"sample_sizes": [int(n)]}
    label = f"img{index}"
    return {
        "id": label,
        "name": label,
        "metadata": meta,
        "analyses": [
            {
                "id": label,
                "name": "1",
                "metadata": meta,
                "points": [],
                "images": [
                    {"url": gp, "filename": "g", "value_type": "g", "space": "MNI"},
                    {"url": vp, "filename": "v", "value_type": "g_var", "space": "MNI"},
                ],
            }
        ],
    }


if __name__ == "__main__":
    frame = pd.read_parquet(CORPUS)
    by_study = {
        sid: g[["x", "y", "z"]].to_numpy(dtype=float) for sid, g in frame.groupby("study_id")
    }
    all_sids = sorted(by_study)
    # A reported t on an assumed n - 1 degrees of freedom, put on the z scale the model reads.
    by_stat = None
    if USE_REPORTED_STATS:
        from nimare.transforms import t_to_z

        by_stat = {}
        for sid, group in frame.groupby("study_id"):
            t = group["t_stat"].to_numpy(dtype=float)
            z = np.full(t.shape, ASSUMED_Z)
            ok = np.isfinite(t) & (np.abs(t) > 0)
            if ok.any():
                z[ok] = np.abs(t_to_z(np.abs(t[ok]), ASSUMED_N - 1))
            by_stat[sid] = z

    collection = load_pain()
    maps, sizes = [], []
    for row, n in zip(collection.images.itertuples(), collection.sample_sizes()):
        maps.append(np.nan_to_num(masker.transform(row.z).ravel()))
        sizes.append(int(n if np.ndim(n) == 0 else np.asarray(n).ravel()[0]))
    maps = np.array(maps)
    total = len(maps)
    donors_available = total - N_REFERENCE
    image_counts = [k for k in IMAGE_COUNTS if k <= donors_available]

    print(f"pain tables: {len(all_sids)} studies, {len(frame)} peaks", flush=True)
    print(
        f"pain images: {total} NIDM studies, {N_REFERENCE} held out as the reference, "
        f"{donors_available} available as donors",
        flush=True,
    )
    print(
        f"coordinate sample size {'drawn 12-120' if SPREAD_N else ASSUMED_N}, "
        f"threshold z {ASSUMED_Z:.2f} clamped per study, peak statistics "
        f"{'as reported' if USE_REPORTED_STATS else 'set to the assumed cut (clamp inert)'}",
        flush=True,
    )
    print(
        f"table counts {TABLE_COUNTS}, image counts {image_counts}, " f"{N_DRAWS} draws\n",
        flush=True,
    )

    rng = np.random.default_rng(0)
    results = {}
    for draw in range(N_DRAWS):
        order = rng.permutation(total)
        held, donors = order[:N_REFERENCE], order[N_REFERENCE:]
        truth = np.abs(pooled(held, maps, sizes))  # fixed for every cell of this draw
        for n_tables in TABLE_COUNTS:
            chosen = list(
                rng.choice(all_sids, size=min(n_tables, len(all_sids)), replace=False)
            )
            table_n = (
                rng.integers(12, 121, size=len(chosen))
                if SPREAD_N
                else np.full(len(chosen), ASSUMED_N)
            )
            tables = [
                coordinate_study(
                    sid,
                    by_study[sid],
                    table_n[i],
                    by_stat[sid] if by_stat is not None else None,
                )
                for i, sid in enumerate(chosen)
            ]
            peaks = sum(len(by_study[sid]) for sid in chosen)
            for k, report_radius in ((k, r) for k in image_counts for r in REPORT_RADII):
                members = list(donors[:k])
                studies = list(tables)
                for i in members:
                    g, var = g_and_var(maps[i], sizes[i])
                    studies.append(image_study(int(i), g, var, sizes[i]))
                est = CBES(
                    mask=masker,
                    null_method="none",
                    threshold=ASSUMED_Z,
                    clamp_threshold=True,
                    report_radius=report_radius,
                )
                res = est.fit(
                    Studyset(
                        {"id": "p", "name": "p", "studies": studies},
                        target=None,
                        mask=mask_img,
                    )
                )
                cbes = np.abs(res.get_map("g", return_type="array").ravel())
                covered = res.get_map("n_studies", return_type="array").ravel() > 0
                share = res.get_map("coordinate_share", return_type="array").ravel()
                # The mixture's own account of itself. If the table wall is truthfully silent
                # because most of these studies do not carry this contrast's effect, then pi
                # should be small while mu stays large. `g` coming out low with a small pi
                # instead means the two are not separating.
                have = set(res.maps)
                prev = (
                    res.get_map("prevalence", return_type="array").ravel()
                    if "prevalence" in have
                    else np.full(cbes.shape, np.nan)
                )
                marg = (
                    np.abs(res.get_map("g_marginal", return_type="array").ravel())
                    if "g_marginal" in have
                    else np.full(cbes.shape, np.nan)
                )
                only = np.abs(pooled(members, maps, sizes))
                use = covered & np.isfinite(cbes) & np.isfinite(only) & np.isfinite(truth)
                if use.sum() < 100:
                    continue
                s_cbes, s_pool = score(cbes, truth, use), score(only, truth, use)
                strong = truth >= np.nanpercentile(truth[use], 90)
                s_marg = (
                    score(marg, truth, use) if np.isfinite(marg[use]).all() else (np.nan,) * 6
                )
                results.setdefault((n_tables, k, report_radius), []).append(
                    (
                        s_cbes[5],
                        s_pool[5],
                        s_cbes[4],
                        s_pool[4],
                        float(np.median(share[use])),
                        float(peaks),
                        float(np.nanmedian(prev[use])),
                        float(np.nanmedian(prev[use & strong])),
                        s_marg[5],
                        s_marg[4],
                        float(np.nanmedian(cbes[use & strong])),
                        float(np.nanmedian(truth[use & strong])),
                    )
                )
                print(
                    f"  draw {draw + 1}, {n_tables} tables ({peaks} peaks), {k} image(s): "
                    f"r{report_radius or 0:.0f}mm rmse {s_cbes[5]:.3f} against {s_pool[5]:.3f}, "
                    f"share {np.median(share[use]):.3f}",
                    flush=True,
                )

    print(
        f"\n{'tables':>7} {'peaks':>7} {'images':>7} | {'rmse CBES':>10} {'rmse pool':>10} "
        f"{'diff':>7} {'p':>7} | {'top CBES':>9} {'top pool':>9} | {'share':>7} "
        f"{'prev':>7} {'prev@top':>8} | {'rmse marg':>9} {'|g|@top':>8} {'truth':>8}"
    )
    for n_tables in TABLE_COUNTS:
        for k, report_radius in ((k, r) for k in image_counts for r in REPORT_RADII):
            rows = np.array(results.get((n_tables, k, report_radius), []), dtype=float)
            if not len(rows):
                continue
            m = rows.mean(axis=0)
            d = rows[:, 0] - rows[:, 1]
            p = stats.ttest_rel(rows[:, 0], rows[:, 1]).pvalue if len(rows) > 1 else np.nan
            print(
                f"{n_tables:7d} {int(m[5]):7d} {k:4d}/{report_radius or 0:.0f}mm | {m[0]:10.3f} {m[1]:10.3f} "
                f"{d.mean():+7.3f} {p:7.3f} | {m[2]:+9.3f} {m[3]:+9.3f} | {m[4]:7.3f} "
                f"{m[6]:7.3f} {m[7]:8.3f} | {m[8]:9.3f} {m[10]:8.3f} {m[11]:8.3f}"
            )
    print(
        "\nA negative rmse difference means the tables are paying for themselves. The "
        "prediction is that\nthe crossing moves to fewer images as the table wall grows, "
        "because a fixed bias gets more weight."
    )
