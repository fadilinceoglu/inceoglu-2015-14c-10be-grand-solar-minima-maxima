"""Download immutable public upstream inputs and record their provenance."""

from __future__ import annotations

import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from .io import sha256, write_csv, write_json


PUBLIC_FILES = {
    "intcal13.14c": {
        "url": "https://www.intcal.org/curves/intcal13.14c",
        "sha256": "0af281d2b559143ea9230c45e3c709137d868b8dbc09abf1ad16739fdb424cbd",
        "citation_doi": "10.2458/azu_js_rc.55.16947",
        "role": "Published IntCal13 atmospheric calibration curve named in the paper.",
        "rights": "License not stated on the per-file download; scholarly citation required.",
    },
    "intcal04.14c": {
        "url": "https://www.intcal.org/curves/intcal04.14c",
        "sha256": "e119d6235e11a08414aa75aae66db5d2ccfa6dfd75c77ce874f6af1d1c7838f8",
        "citation_doi": "10.1017/S0033822200032999",
        "role": "IntCal04 comparison used by the cited Knudsen et al. processing chain.",
        "rights": "License not stated on the per-file download; scholarly citation required.",
    },
    "grip_10be.txt": {
        "url": "https://www.ncei.noaa.gov/pub/data/paleo/icecore/greenland/summit/grip/cosmoiso/grip_10be.txt",
        "sha256": "ad70f6be5c11adee1db1484c51969dc7497c538580f69d56fc3db744ae521ddb",
        "citation_doi": "10.25921/hj8n-xc03",
        "role": "Authoritative upstream GRIP measurements; not the processed modulation potential.",
        "rights": "Externally contributed NCEI holding; dataset-specific license not stated; citation required.",
    },
    "grip_10be.xls": {
        "url": "https://www.ncei.noaa.gov/pub/data/paleo/icecore/greenland/summit/grip/cosmoiso/grip_10be.xls",
        "sha256": "cb0d30b0b70932a91197fba59c27cb481475cfc83ef342afaa7dcac56f08ad2d",
        "citation_doi": "10.25921/hj8n-xc03",
        "role": "Spreadsheet form of the authoritative upstream GRIP archive.",
        "rights": "Externally contributed NCEI holding; dataset-specific license not stated; citation required.",
    },
}


def prepare_public_intcal13(repository_root: Path) -> dict[str, object]:
    """Select the overlap interval and linearly interpolate IntCal13 annually.

    This performs the public-data preparation described by the paper up to
    Delta14C. Conversion from Delta14C to production and Phi is outside this
    function and is not implemented by the repository.
    """

    source = repository_root / "data" / "raw" / "public" / "intcal13.14c"
    values = np.genfromtxt(
        source, delimiter=",", comments="#", dtype=float, encoding="latin-1"
    )
    values = values[np.all(np.isfinite(values), axis=1)]
    order = np.argsort(values[:, 0])
    age, delta14c, sigma = values[order, 0], values[order, 3], values[order, 4]
    annual_age = np.arange(300.0, 8551.0, 1.0)
    annual_delta = np.interp(annual_age, age, delta14c)
    annual_sigma = np.interp(annual_age, age, sigma)
    target = repository_root / "data" / "processed" / "public" / "intcal13_delta14c_annual.csv"
    rows = (
        {"age_bp": int(item_age), "delta14c_per_mil": item_delta, "sigma_per_mil": item_sigma}
        for item_age, item_delta, item_sigma in zip(annual_age, annual_delta, annual_sigma, strict=True)
    )
    write_csv(target, ["age_bp", "delta14c_per_mil", "sigma_per_mil"], rows)
    return {
        "file": str(target.relative_to(repository_root)),
        "sha256": sha256(target),
        "rows": int(annual_age.size),
        "range_bp": [300, 8550],
        "operation": "select IntCal13 overlap interval and linearly interpolate Delta14C annually",
        "causal_boundary": "Delta14C-to-production/Phi conversion unavailable",
    }


def download_public(repository_root: Path, force: bool = False) -> dict[str, object]:
    destination = repository_root / "data" / "raw" / "public"
    destination.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    for filename, metadata in PUBLIC_FILES.items():
        target = destination / filename
        if target.exists() and not force:
            observed_hash = sha256(target)
            if observed_hash != metadata["sha256"]:
                raise ValueError(f"Existing public input fails checksum: {target}")
            records.append({"file": filename, **metadata, "bytes": target.stat().st_size})
            continue
        request = urllib.request.Request(
            str(metadata["url"]), headers={"User-Agent": "Inceoglu-2015-reproduction/1.0"}
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read()
            response_metadata = {
                "final_url": response.geturl(),
                "content_type": response.headers.get("Content-Type"),
                "content_length": response.headers.get("Content-Length"),
                "etag": response.headers.get("ETag"),
                "last_modified": response.headers.get("Last-Modified"),
            }
        target.write_bytes(payload)
        observed_hash = sha256(target)
        if observed_hash != metadata["sha256"]:
            target.unlink()
            raise ValueError(
                f"Downloaded bytes differ from the expected public-source checksum for {filename}: "
                f"expected {metadata['sha256']}, observed {observed_hash}"
            )
        records.append(
            {
                "file": filename,
                **metadata,
                **response_metadata,
                "bytes": len(payload),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            }
        )
    prepared = prepare_public_intcal13(repository_root)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "records": records,
        "prepared": [prepared],
        "causal_boundary": (
            "The GRIP archive is upstream measurement data. It must not be substituted "
            "for the missing Vonmoos/Knudsen processed production and Phi series."
        ),
    }
    write_json(destination / "manifest.json", manifest)
    return manifest
