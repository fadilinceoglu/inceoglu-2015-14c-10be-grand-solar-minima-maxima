"""One-command reproduction pipeline."""

from __future__ import annotations

from pathlib import Path

from .analysis import run_analysis
from .figures import write_all_figures
from .tables import write_tables
from .validation import validate


def reproduce(repository_root: Path) -> dict[str, object]:
    repository_root = repository_root.resolve()
    aligned_data = repository_root / "data" / "derived"
    analysis = run_analysis(
        aligned_data, repository_root / "catalog" / "manual_event_catalog.csv"
    )
    write_tables(analysis, repository_root / "outputs" / "tables")
    write_all_figures(analysis, repository_root / "outputs" / "figures")
    return validate(analysis, repository_root, repository_root / "outputs" / "validation")
