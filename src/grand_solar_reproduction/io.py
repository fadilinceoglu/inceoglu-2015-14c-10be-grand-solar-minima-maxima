"""Input/output helpers with explicit format and checksum validation."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        writer.writerows(rows)


def load_two_column_series(path: Path) -> tuple[np.ndarray, np.ndarray]:
    values = np.loadtxt(path, dtype=float)
    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError(f"Expected at least two numeric columns in {path}")
    time = values[:, 0]
    signal = values[:, 1]
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(signal)):
        raise ValueError(f"Non-finite observation in {path}")
    if not np.all(np.diff(time) > 0):
        raise ValueError(f"Time coordinate is not strictly increasing in {path}")
    return time, signal

