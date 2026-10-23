"""Download pinned SXS:BBH:2327 v3.0 horizon/reference files, never waveforms.

CaltechDATA uses an InvenioRDM record, not the legacy Zenodo record schema.
The record does not declare a license: record that fact rather than inferring one.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import requests

ROOT = Path(__file__).resolve().parent
RECORD_URL = "https://data.caltech.edu/api/records/f70bw-gcd50"
RECORD_PATH = ROOT / "adapter-sources/sxs-2327-record.json"


def fetch_level(level):
    record = json.loads(RECORD_PATH.read_text(encoding="utf-8"))
    assert record["pids"]["doi"]["identifier"] == "10.26138/SXS:BBH:2327v3.0"
    out = ROOT / "sxs-data/2327" / level
    out.mkdir(parents=True, exist_ok=True)
    manifest = {
        "simulation": "SXS:BBH:2327", "record_id": record["id"],
        "repository": "CaltechDATA", "record_url": RECORD_URL,
        "doi": record["pids"]["doi"]["identifier"], "version": "3.0", "level": level,
        "license": {"id": "unspecified", "status": "not_declared_in_record",
                    "note": "Public open-access record; no metadata.rights field, and DataCite rightsList is empty. Do not infer CC-BY-4.0 from older Zenodo simulations.",
                    "additional_source": "https://api.datacite.org/dois/10.26138/SXS:BBH:2327v3.0"},
        "record_sha256": hashlib.sha256(RECORD_PATH.read_bytes()).hexdigest(),
        "files": [],
    }
    for name in ("Horizons.h5", "metadata.json"):
        source = record["files"]["entries"][f"{level}:{name}"]
        dest = out / name
        url = source["links"]["content"]
        if not dest.exists():
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            content = response.content
            assert len(content) == source["size"], (name, "size mismatch")
            assert "md5:" + hashlib.md5(content).hexdigest() == source["checksum"]
            dest.write_bytes(content)
        content = dest.read_bytes()
        assert len(content) == source["size"]
        assert "md5:" + hashlib.md5(content).hexdigest() == source["checksum"]
        manifest["files"].append({"file": name, "url": url, "bytes": len(content),
                                  "md5": source["checksum"],
                                  "sha256": hashlib.sha256(content).hexdigest()})
    metadata = json.loads((out / "metadata.json").read_text())
    manifest["status_checks"] = {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "keywords": metadata.get("keywords", []),
        "superseded_by": metadata.get("superseded_by"),
        "record_is_latest_at_retrieval": record["versions"]["is_latest"],
    }
    catalog_path = ROOT / "adapter-sources/sxs-2327-catalog-status.json"
    if catalog_path.exists():
        catalog = json.loads(catalog_path.read_text())
        current = catalog["entries"]["SXS:BBH:2327"]
        manifest["status_checks"]["catalog"] = {
            "source": catalog["source"], "source_sha256": catalog["source_sha256"],
            "keywords": current.get("keywords", []), "DOI_versions": current["DOI_versions"],
        }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"level": level, "bytes": sum(f["bytes"] for f in manifest["files"]),
                      "status_checks": manifest["status_checks"]}), flush=True)


if __name__ == "__main__":
    if not RECORD_PATH.exists():
        response = requests.get(RECORD_URL, timeout=60)
        response.raise_for_status()
        record = response.json()
        assert record["pids"]["doi"]["identifier"] == "10.26138/SXS:BBH:2327v3.0"
        RECORD_PATH.write_text(json.dumps(record, indent=2), encoding="utf-8")
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(fetch_level, ["Lev1", "Lev2", "Lev3"]))
