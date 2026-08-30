"""Generate all five published tables from calculated values."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .analysis import StudyAnalysis
from .events import table_event_rows
from .io import write_csv


def _nearest_integer(value: float) -> int:
    return int(np.floor(value + 0.5))


def write_tables(analysis: StudyAnalysis, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    table1_rows = []
    for nuclide, fit in [
        ("10Be", analysis.classification.be10_mixture),
        ("14C", analysis.classification.c14_mixture),
    ]:
        amplitudes = np.array(
            [
                getattr(pair, "be10" if nuclide == "10Be" else "c14").amplitude_sigma
                for pair in analysis.classification.pairs
            ]
        )
        from .models import fit_gaussian

        gaussian = fit_gaussian(amplitudes)
        table1_rows.append(
            {
                "nuclide": nuclide,
                "gaussian_bic": _nearest_integer(gaussian.criteria.bic),
                "gaussian_aic": _nearest_integer(gaussian.criteria.aic),
                "bimodal_gaussian_bic": _nearest_integer(fit.criteria.bic),
                "bimodal_gaussian_aic": _nearest_integer(fit.criteria.aic),
            }
        )
    fields1 = list(table1_rows[0])
    write_csv(output_dir / "table1_information_criteria.csv", fields1, table1_rows)

    for number, name, pairs in [
        (2, "grand_minima", analysis.classification.minima),
        (3, "grand_maxima", analysis.classification.maxima),
    ]:
        rows = table_event_rows(pairs)
        fields = list(rows[0])
        write_csv(output_dir / f"table{number}_{name}.csv", fields, rows)

    table4_rows = []
    for name, result in [
        ("Grand minima", analysis.minimum_waiting),
        ("Grand maxima", analysis.maximum_waiting),
    ]:
        table4_rows.append(
            {
                "event_type": name,
                "power_law_alpha": f"{result.power_law.alpha:.2f}",
                "power_law_bic": _nearest_integer(result.power_law.criteria.bic),
                "power_law_aic": _nearest_integer(result.power_law.criteria.aic),
                "exponential_tau_years": _nearest_integer(result.exponential.tau),
                "exponential_bic": _nearest_integer(result.exponential.criteria.bic),
                "exponential_aic": _nearest_integer(result.exponential.criteria.aic),
            }
        )
    fields4 = list(table4_rows[0])
    write_csv(output_dir / "table4_waiting_time_models.csv", fields4, table4_rows)

    table5_rows = []
    for name, result in [
        ("Grand maxima", analysis.maximum_durations),
        ("Grand minima", analysis.minimum_durations),
    ]:
        table5_rows.append(
            {
                "event_type": name,
                "gaussian_bic": _nearest_integer(result.gaussian.criteria.bic),
                "gaussian_aic": _nearest_integer(result.gaussian.criteria.aic),
                "bimodal_gaussian_bic": _nearest_integer(result.gaussian_mixture.criteria.bic),
                "bimodal_gaussian_aic": _nearest_integer(result.gaussian_mixture.criteria.aic),
                "lognormal_bic": _nearest_integer(result.lognormal.criteria.bic),
                "lognormal_aic": _nearest_integer(result.lognormal.criteria.aic),
            }
        )
    fields5 = list(table5_rows[0])
    write_csv(output_dir / "table5_duration_models.csv", fields5, table5_rows)
