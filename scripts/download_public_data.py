#!/usr/bin/env python3
"""Download checksum-pinned authoritative public inputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from grand_solar_reproduction.data import download_public  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="Re-fetch existing files")
    args = parser.parse_args()
    manifest = download_public(ROOT, force=args.force)
    print(f"Verified {len(manifest['records'])} public files under data/raw/public")


if __name__ == "__main__":
    main()

