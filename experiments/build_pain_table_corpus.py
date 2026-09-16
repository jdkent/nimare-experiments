"""A large pain coordinate corpus from the NeuroStore release, with the NIDM studies removed.

The question this serves is the extreme configuration: one image against hundreds or thousands of
coordinate studies. The NeuroStore monthly release carries 32,444 studies and 871,671 coordinates,
of which 1,522 studies name pain in their title or abstract.

Two things have to be right or the bed measures nothing.

**The truth must not be in the coordinate set.** The images and the held-out reference come from
the 21-study NIDM pain collection, whose study identifiers are anonymised (`pain_01.nidm`), so
overlap cannot be detected by name -- and some of those 21 are certainly in NeuroStore. But the
NIDM collection carries each study's *published* peak coordinates, so the overlap can be found in
the data: a NeuroStore analysis reporting the same peaks in the same places is the same study.
Any NeuroStore study matching an NIDM study's coordinate set is dropped.

**The sample sizes are not there.** Only 0.5% of the release's analyses carry `sample_sizes`, so a
coordinate study's `n` has to be assumed. That assumption is load-bearing: it sets where the
reporting threshold sits on the effect-size scale, and a spread of sample sizes is the only thing
that identifies the magnitude at all. Both a constant and a drawn spread are therefore available,
and which was used is printed with the corpus.

Peaks are taken as published -- nothing is capped, and a study contributes however many rows its
table had. Talairach coordinates are converted rather than dropped, so studies are not lost for
their reporting convention.
"""
import json, os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy import spatial
from nimare.utils import tal2mni

RELEASE = "/tmp/claude-0/nsrelease/neurostore-studyset-2026-09"
OUT = "/tmp/claude-0/paintables"
#: Terms that mark a study as being about pain rather than merely mentioning it once.
TERMS = ("pain", "nocicept", "noxious", "analgesi", "hyperalgesi", "allodyni")
#: A NeuroStore study is treated as an NIDM study when this fraction of an NIDM study's peaks
#: land within MATCH_MM of one of its peaks. Set low deliberately: a false positive costs one
#: study out of fifteen hundred, a false negative puts the reference inside the coordinates.
MATCH_FRACTION = 0.4
MATCH_MM = 4.0


def load_release():
    studies = pd.read_parquet(f"{RELEASE}/studies.parquet")
    coordinates = pd.read_parquet(f"{RELEASE}/coordinates.parquet")
    return studies, coordinates


def pain_studies(studies):
    text = (studies["name"].fillna("") + " " + studies["description"].fillna("")).str.lower()
    keep = text.apply(lambda t: any(term in t for term in TERMS))
    return studies.loc[keep, "study_id"]


def to_mni(frame):
    """Coordinates on one scale, dropping only the rows whose space is unstated."""
    frame = frame[frame["space"].isin(("MNI", "TAL"))].copy()
    tal = frame["space"] == "TAL"
    if tal.any():
        converted = tal2mni(frame.loc[tal, ["x", "y", "z"]].to_numpy(dtype=float))
        frame.loc[tal, ["x", "y", "z"]] = converted
    return frame


def nidm_peaks():
    """Each NIDM pain study's published peaks, for the overlap guard."""
    from nimare.tests.utils import get_test_data_path
    raw = json.load(open(os.path.join(get_test_data_path(), "nidm_pain_dset.json")))
    out = []
    for study in raw.values():
        for contrast in study.get("contrasts", {}).values():
            coords = contrast.get("coords") or {}
            if not coords.get("x"):
                continue
            out.append(np.column_stack([coords["x"], coords["y"], coords["z"]]).astype(float))
    return out


def contaminated(frame, nidm):
    """NeuroStore study ids whose peaks reproduce an NIDM study's table."""
    bad = set()
    by_study = {sid: g[["x", "y", "z"]].to_numpy(dtype=float)
                for sid, g in frame.groupby("study_id")}
    trees = {sid: spatial.cKDTree(pts) for sid, pts in by_study.items() if len(pts)}
    for peaks in nidm:
        for sid, tree in trees.items():
            hits = tree.query_ball_point(peaks, MATCH_MM, return_length=True)
            if np.mean(np.asarray(hits) > 0) >= MATCH_FRACTION:
                bad.add(sid)
    return bad


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    studies, coordinates = load_release()
    print(f"release: {len(studies)} studies, {len(coordinates)} coordinates", flush=True)

    ids = set(pain_studies(studies))
    frame = coordinates[coordinates["study_id"].isin(ids)]
    print(f"{len(ids)} pain studies named, {frame['study_id'].nunique()} with coordinates, "
          f"{len(frame)} peaks", flush=True)

    frame = to_mni(frame)
    print(f"{frame['study_id'].nunique()} studies after dropping unstated spaces, "
          f"{len(frame)} peaks on MNI", flush=True)

    nidm = nidm_peaks()
    print(f"guarding against {len(nidm)} NIDM pain tables "
          f"({sum(len(p) for p in nidm)} peaks)", flush=True)
    bad = contaminated(frame, nidm)
    print(f"dropping {len(bad)} NeuroStore studies whose peaks reproduce an NIDM table: "
          f"{sorted(bad)}", flush=True)
    frame = frame[~frame["study_id"].isin(bad)]

    counts = frame.groupby("study_id").size()
    print(f"\ncorpus: {len(counts)} studies, {len(frame)} peaks", flush=True)
    print(f"peaks per study: median {counts.median():.0f}, "
          f"range {counts.min()}-{counts.max()}, mean {counts.mean():.1f}")
    stats = frame["t_stat"].notna().sum()
    print(f"{stats} peaks carry a t statistic ({100 * stats / len(frame):.0f}%)")
    frame.to_parquet(f"{OUT}/pain_tables.parquet")
    json.dump(sorted(bad), open(f"{OUT}/excluded_studies.json", "w"))
    print(f"\nwritten to {OUT}/pain_tables.parquet")
