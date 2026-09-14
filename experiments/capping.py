"""Cap a simulated studyset to a realistic number of reported peaks, keeping the truth.

``cap_peaks`` rebuilds the studyset from NIMADS points that carry only a Z value, which drops
the simulator's ``value_trueg`` column -- and that column is the only record of what the true
effect at each reported location was. Anything scoring against it has to take the kept rows
directly rather than re-reading them off the rebuilt studyset.
"""
import numpy as np
from nimare.studyset import Studyset
from nimare.transforms import t_to_z


def metadata_lookup(studyset):
    """``{(study_id, contrast_id): metadata}``, since the coordinate table keys on both."""
    return {
        (str(study.id), str(analysis.id)): dict(analysis.metadata)
        for study in studyset.studies
        for analysis in study.analyses
    }


def cap_rows(studyset, max_peaks):
    """Return the strongest ``max_peaks`` foci per analysis, with every simulator column."""
    coords = studyset.coordinates.copy()
    coords["_abs"] = np.abs(coords["z_stat"].astype(float))
    if max_peaks is None:
        return coords
    return (
        coords.sort_values("_abs", ascending=False)
        .groupby(["study_id", "contrast_id"], sort=False)
        .head(max_peaks)
    )


def build(rows, meta, statistic="reported"):
    """Build a studyset from kept rows, reporting either the real z or the truth's z.

    ``statistic="truth"`` keeps each focus where selection put it but replaces its value with
    the statistic the true effect at that location would have produced, so the locations are
    still selected while the magnitudes carry no winner's curse.
    """
    studies = []
    for (study_id, contrast_id), sub in rows.groupby(["study_id", "contrast_id"], sort=False):
        info = meta[(str(study_id), str(contrast_id))]
        n = float(info["sample_sizes"][0])
        correction = 1.0 - 3.0 / (4.0 * (n - 1) - 1.0)
        points = []
        for row in sub.itertuples():
            if statistic == "truth":
                value = float(
                    t_to_z(
                        np.array([float(row.value_trueg) / correction * np.sqrt(n)]),
                        dof=n - 1,
                    )[0]
                )
            else:
                value = float(row.z_stat)
            points.append({
                "space": "MNI",
                "coordinates": [float(row.x), float(row.y), float(row.z)],
                "values": [{"kind": "Z", "value": value}],
            })
        studies.append({
            "id": str(study_id), "name": str(study_id), "metadata": info,
            "analyses": [{"id": str(contrast_id), "name": "1",
                          "metadata": info, "points": points}],
        })
    return Studyset({"id": "capped", "name": "capped", "studies": studies})
