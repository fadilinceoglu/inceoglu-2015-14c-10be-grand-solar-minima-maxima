#!/usr/bin/env python3
"""Run data preparation, all calculations, all tables/figures, and tests."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib-cache"))
sys.path.insert(0, str(ROOT / "src"))

from grand_solar_reproduction.pipeline import reproduce  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reproduce Inceoglu et al. (2015), A&A 577, A20"
    )
    parser.add_argument(
        "--skip-tests", action="store_true", help="Skip the final pytest invocation"
    )
    args = parser.parse_args()
    report = reproduce(ROOT)
    if not args.skip_tests:
        environment = os.environ.copy()
        environment["PYTHONPATH"] = str(ROOT / "src")
        subprocess.run(
            [sys.executable, "-m", "pytest", "-q"],
            cwd=ROOT,
            env=environment,
            check=True,
        )
    print(
        "Reproduction complete: "
        f"{report['computed']['grand_minima']} grand minima, "
        f"{report['computed']['grand_maxima']} grand maxima; validation passed."
    )


if __name__ == "__main__":
    main()
