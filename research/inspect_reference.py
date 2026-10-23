"""Verify SXS 2327 files and report reference physics and data-domain checks."""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent
FIELDS = {"position": "CoordCenterInertial.dat", "mass": "ChristodoulouMass.dat", "spin": "chiInertial.dat"}


def main():
    catalog = json.loads((ROOT / "adapter-sources/sxs-2327-catalog-status.json").read_text())
    current = catalog["entries"]["SXS:BBH:2327"]
    assert "deprecated" not in current.get("keywords", [])
    result = {"simulation": "SXS:BBH:2327", "doi": "10.26138/SXS:BBH:2327v3.0",
              "status": "not_marked_deprecated_in_pinned_official_catalog_or_downloaded_metadata",
              "catalog_source": catalog["source"], "catalog_source_sha256": catalog["source_sha256"],
              "license": "not_declared_in_CaltechDATA_record",
              "limitations": ["Coordinate horizon diagnostics are gauge dependent.",
                              "Preserve per-level PBJ metadata. Lev1/2 report transition=0 and different base_lev; Lev3 reports base_lev=Lev3, transition=1418.353. Do not assume all available levels share an identical prefix.",
                              "0003 and 2327 have different initial-data constructions and reference times; their difference is not a resolution error."],
              "levels": []}
    raw = {}
    for level in ("Lev1", "Lev2", "Lev3"):
        directory = ROOT / "sxs-data/2327" / level
        manifest = json.loads((directory / "manifest.json").read_text())
        for item in manifest["files"]:
            content = (directory / item["file"]).read_bytes()
            assert "md5:" + hashlib.md5(content).hexdigest() == item["md5"]
            assert hashlib.sha256(content).hexdigest() == item["sha256"]
        meta = json.loads((directory / "metadata.json").read_text())
        omega = np.asarray(meta["reference_orbital_frequency"])
        lhat = omega / np.linalg.norm(omega)
        spins = [np.asarray(meta[f"reference_dimensionless_spin{i}"]) for i in (1, 2)]
        masses = [meta[f"reference_mass{i}"] for i in (1, 2)]
        entry = {"level": level, "reference_mass_ratio": masses[0] / masses[1],
                 "reference_mass_total": sum(masses), "reference_time": meta["reference_time"],
                 "reference_spin_magnitudes": [float(np.linalg.norm(s)) for s in spins],
                 "reference_transverse_spins": [float(np.linalg.norm(s - np.dot(s, lhat) * lhat)) for s in spins],
                 "reference_chi_eff": float(sum(m * np.dot(s, lhat) for m, s in zip(masses, spins)) / sum(masses)),
                 "reference_eccentricity": meta["reference_eccentricity"],
                 "pbj_info": meta["pbj_info"], "horizons": {}}
        raw[level] = {}
        with h5py.File(directory / "Horizons.h5") as file:
            entry["hdf5_sxs_format"] = str(file.attrs.get("sxs_format", "not_set; detected documented SpEC groups"))
            assert all(f"Ah{h}.dir" in file for h in "ABC")
            for horizon in "ABC":
                raw[level][horizon] = {}
                stats = {}
                for field, name in FIELDS.items():
                    a = np.asarray(file[f"Ah{horizon}.dir/{name}"])
                    assert a.shape[1] == (2 if field == "mass" else 4)
                    assert np.isfinite(a).all() and np.all(np.diff(a[:, 0]) > 0)
                    if field == "mass":
                        assert np.all(a[:, 1] > 0)
                    if field == "spin":
                        assert np.all(np.linalg.norm(a[:, 1:], axis=1) <= 1)
                    raw[level][horizon][field] = a
                    stats[field] = {"samples": len(a), "time_min": float(a[0, 0]), "time_max": float(a[-1, 0])}
                entry["horizons"][horizon] = stats
        common = max(a[0, 0] for a in raw[level]["C"].values())
        assert abs(common - meta["common_horizon_time"]) < 1e-5
        entry["common_horizon_time"] = common
        binary_end = min(a[-1, 0] for h in "AB" for a in raw[level][h].values())
        entry["pre_common_gap_source_units"] = max(0.0, common - binary_end)
        entry["binary_common_overlap_source_units"] = max(0.0, binary_end - common)
        result["levels"].append(entry)
    # Test rather than assume shared prefixes; available levels name different bases.
    result["prefix_equality_checks"] = []
    transition = max(x["pbj_info"]["transition_time"] for x in result["levels"])
    for level in ("Lev1", "Lev2"):
        for horizon in "AB":
            for field in FIELDS:
                a = raw[level][horizon][field]
                b = raw["Lev3"][horizon][field]
                a = a[a[:, 0] < transition]
                b = b[b[:, 0] < transition]
                assert len(a) > 0 and len(b) > 0
                result["prefix_equality_checks"].append({"level": level, "horizon": horizon, "field": field,
                    "transition": transition, "samples": len(a), "exactly_equal_to_Lev3": np.array_equal(a, b)})
    out = ROOT / "adapter-results/reference-2327-validation.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(out), "levels": [{k: v for k, v in e.items() if k != "horizons"} for e in result["levels"]],
                      "prefix_equality_checks": result["prefix_equality_checks"]}, indent=2))


if __name__ == "__main__":
    main()
