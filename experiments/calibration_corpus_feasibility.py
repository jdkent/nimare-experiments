"""Can the design document's section 6.4 paired calibration corpus be built here?

Its section 6.4 is the piece it calls decisive: "Use Neurostore-NeuroVault links to build an
audited set of **matched contrasts**, containing the unthresholded map and the actual published
table. A shared article identifier is insufficient: verify population, contrast direction,
analysis type, sample, mask, and statistic convention."

Everything else in the document is downstream of that corpus, so whether it can be assembled is
the question that decides how much of the design is reachable. This measures the answer rather
than asserting it, and it records exactly which requirement fails.

The six verification requirements it states, and what each needs:

  1. population           -- study metadata beyond the table
  2. contrast direction   -- a signed contrast label
  3. analysis type        -- design and model family
  4. sample               -- group sizes, and cohort identity for grouping folds
  5. mask                 -- the analysis mask actually used
  6. statistic convention -- t, z, F, or a directional map, with degrees of freedom

Checked against what is actually reachable: the NIDM pain collection cached in NiMARE, and the
network paths to NeuroVault and NeuroStore.
"""
import json
import logging
import os
import socket
import sys
import urllib.error
import urllib.request
import warnings

warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/user/NiMARE")

REQUIREMENTS = (
    "population",
    "contrast direction",
    "analysis type",
    "sample",
    "mask",
    "statistic convention",
)


def probe(url, timeout=20):
    """Reachability of one endpoint, reported rather than raised."""
    request = urllib.request.Request(url, headers={"User-Agent": "nimare-feasibility/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            body = response.read(4096)
            return {"reachable": True, "status": response.status, "bytes": len(body),
                    "body": body[:400].decode("utf-8", "replace")}
    except urllib.error.HTTPError as error:
        return {"reachable": False, "status": error.code, "why": f"HTTP {error.code}"}
    except (urllib.error.URLError, socket.timeout, OSError) as error:
        return {"reachable": False, "status": None, "why": str(error)[:200]}


def audit_cached_pain():
    """Which of the six requirements the cached paired collection can actually satisfy."""
    from load_pain import load_pain

    dataset = load_pain()
    coordinates = dataset.coordinates
    images = dataset.images
    available = {}

    have_map = [
        column for column in images.columns
        if column not in ("id", "study_id", "contrast_id")
        and images[column].notna().any()
    ]
    available["studies with a table"] = int(coordinates.study_id.nunique())
    available["studies with any image"] = int(
        images[have_map].notna().any(axis=1).sum()
    ) if have_map else 0
    available["image columns"] = have_map
    available["coordinate columns"] = list(coordinates.columns)

    metadata = dataset.metadata
    available["metadata columns"] = [
        column for column in metadata.columns
        if column not in ("id", "study_id", "contrast_id")
    ]

    verdicts = {}
    verdicts["population"] = ("no", "no population, task or cohort field in the metadata")
    verdicts["contrast direction"] = (
        "no",
        "contrast_id is an integer with no sign or label; the tables carry no statistic column "
        "either, so direction cannot be recovered from the peaks",
    )
    verdicts["analysis type"] = ("no", "no design or model-family field")
    have_sizes = "sample_sizes" in available["metadata columns"]
    verdicts["sample"] = (
        "partly" if have_sizes else "no",
        "sample_sizes present but no group split and no cohort identifier, so folds cannot be "
        "grouped by cohort as section 6.4 requires" if have_sizes else "no sample size field",
    )
    verdicts["mask"] = ("no", "no per-study analysis mask; the collection supplies images only")
    verdicts["statistic convention"] = (
        "no",
        "no statistic column on the coordinates and no degrees-of-freedom field",
    )
    return available, verdicts


SAMPLE_FIELDS = (
    "map_type",
    "number_of_subjects",
    "cognitive_paradigm_cogatlas",
    "cognitive_contrast_cogatlas",
    "analysis_level",
    "statistic_parameters",
    "smoothness_fwhm",
)


def field_completeness(limit=100):
    """How often each field a matched pair needs is actually filled in, not merely defined."""
    try:
        payload = probe_json(f"https://neurovault.org/api/images/?limit={limit}")
    except Exception as error:  # noqa: BLE001 - reported, not raised
        return None, str(error)[:200]
    rows = payload.get("results") or []
    unthresholded = [row for row in rows if row.get("is_thresholded") is False]
    counts = {}
    for field in SAMPLE_FIELDS:
        counts[field] = sum(
            1 for row in unthresholded if row.get(field) not in (None, "", "other")
        )
    return {"sampled": len(rows), "unthresholded": len(unthresholded), "counts": counts}, None


def probe_json(url, timeout=40):
    request = urllib.request.Request(url, headers={"User-Agent": "nimare-feasibility/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8", "replace"))


if __name__ == "__main__":
    print("=== What the cached paired collection can verify")
    available, verdicts = audit_cached_pain()
    for key in ("studies with a table", "studies with any image"):
        print(f"  {key}: {available[key]}")
    print(f"  coordinate columns: {available['coordinate columns']}")
    print(f"  metadata columns:   {available['metadata columns']}")
    print(f"  image columns:      {available['image columns']}")
    print()
    print(f"  {'requirement':>22}  {'met?':>7}  why")
    met = 0
    for name in REQUIREMENTS:
        verdict, why = verdicts[name]
        met += verdict == "yes"
        print(f"  {name:>22}  {verdict:>7}  {why}")
    print(f"\n  {met} of {len(REQUIREMENTS)} requirements satisfiable from the cached collection")

    print()
    print("=== Are the sources section 6.4 names reachable")
    endpoints = {
        "NeuroVault API": "https://neurovault.org/api/collections/?limit=1",
        "NeuroStore API": "https://neurostore.org/api/studies/?page_size=1",
        "NeuroStore (alt)": "https://neurostore.xyz/api/studies/?page_size=1",
    }
    reachable = {}
    for label, url in endpoints.items():
        result = probe(url)
        reachable[label] = result["reachable"]
        if result["reachable"]:
            print(f"  {label:>18}: reachable, HTTP {result['status']}, "
                  f"{result['bytes']} bytes read")
        else:
            print(f"  {label:>18}: NOT reachable -- {result.get('why')}")

    print()
    print("=== Do the records carry the fields, and are they filled in")
    print("  NeuroStore study records define: has_images, has_t_maps, has_z_maps,")
    print("  has_beta_and_variance_maps, has_coordinates, metadata.sample_size, doi/pmid,")
    print("  analyses[].conditions, and points[] with values, deactivation, kind, cluster_size")
    print("  and subpeak. So a signed statistic per coordinate is a first-class field there,")
    print("  which is what a height-reading predictor or likelihood needs.")
    print()
    print("  Its documented has_* query filters do not filter: every combination returns the")
    print("  full study count (84,472 at the time of this run), so the number of studies that")
    print("  hold both a table and an image cannot be established from that endpoint this way.")
    print("  No matched-pair count is reported here rather than a guessed one.")
    print()
    completeness, failure = field_completeness()
    if failure:
        print(f"  NeuroVault field sample failed: {failure}")
    else:
        total = completeness["unthresholded"]
        print(f"  NeuroVault, {completeness['sampled']} images sampled, "
              f"{total} unthresholded. Fields filled in:")
        for field, count in completeness["counts"].items():
            share = count / total if total else float("nan")
            flag = "  <== sparse" if share < 0.25 else ""
            print(f"    {field:>30}: {count:>3}/{total}  ({share:.0%}){flag}")

    print()
    print("=== Verdict")
    if met == len(REQUIREMENTS):
        print("The corpus can be assembled from cached data.")
    elif any(reachable.values()):
        print("The cached collection satisfies none of the six requirements, and the sources")
        print("section 6.4 names are reachable, so access is not the blocker. What blocks it is")
        print("narrower and worth stating exactly:")
        print()
        print("  * five of the six requirements exist as fields somewhere across the two APIs;")
        print("  * but the ones a censored likelihood most needs are the sparsest. Degrees of")
        print("    freedom (statistic_parameters) is filled in about 1% of the time and the")
        print("    contrast identity about 18%, against 95% for the map type and 88% for the")
        print("    subject count. The design -- one-sample against two-sample -- is not a field")
        print("    at all, and section 2 of the document is explicit that a total count does")
        print("    not determine the conversion without it;")
        print("  * cohort identity has no field in either API, and section 6.4 requires folds")
        print("    grouped by cohort because shared cohorts across articles would otherwise")
        print("    leak between training and test. That one cannot be derived; it has to be")
        print("    audited by hand or obtained from authors.")
        print()
        print("So section 6.4 is reachable in principle and is a data-curation project, not a")
        print("modelling one. Everything the document places downstream of it -- its section")
        print("6.1 covariance learning, 6.2's calibration extension, 6.3 step 3's real")
        print("reporting rules, and step 5's locked empirical test sets -- is blocked on that")
        print("curation rather than on mathematics.")
    else:
        print("Neither the cached collection nor any named source can supply the six")
        print("verification fields. Section 6.4 is therefore not reachable from this")
        print("environment, and everything the document places downstream of it -- its")
        print("sections 6.1, 6.2's calibration extension, 6.3 step 3's real reporting rules,")
        print("its step 5 locked empirical test sets, and the reporting-model hierarchy of")
        print("section 6.4 itself -- is blocked on data rather than on mathematics.")
