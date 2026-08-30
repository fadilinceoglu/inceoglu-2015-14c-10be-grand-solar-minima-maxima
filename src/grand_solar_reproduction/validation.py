"""Numerical validation targets transcribed from the publisher paper."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import numpy as np

from .analysis import StudyAnalysis
from .events import table_event_rows
from .io import sha256, write_json


PUBLISHED_MINIMA = [
    (1463, 161, 1450, 167), (1283, 75, 1300, 91), (1027, 47, 1050, 65),
    (901, 31, 900, 28), (676, 107, 690, 102), (417, 106, 434, 54),
    (266, 36, 261, 34), (124, 31, 133, 25), (-344, 82, -348, 107),
    (-658, 35, -660, 22), (-790, 135, -750, 154), (-906, 43, -895, 54),
    (-1190, 41, -1189, 38), (-1358, 125, -1370, 121), (-1488, 61, -1491, 53),
    (-2125, 39, -2132, 42), (-2183, 83, -2195, 25), (-2453, 51, -2461, 74),
    (-2901, 93, -2874, 106), (-3088, 30, -3080, 28), (-3344, 140, -3330, 134),
    (-3479, 79, -3492, 73), (-3627, 80, -3627, 98), (-3699, 24, -3695, 30),
    (-4231, 59, -4224, 67), (-4317, 53, -4322, 72), (-5193, 45, -5205, 53),
    (-5298, 62, -5298, 69), (-5459, 50, -5466, 80), (-5606, 76, -5610, 53),
    (-5718, 51, -5713, 41), (-6445, 123, -6425, 152),
]

PUBLISHED_MAXIMA = [
    (1616, 82, 1604, 83), (1373, 58, 1370, 36), (517, 27, 521, 72),
    (314, 72, 305, 67), (224, 37, 200, 88), (-200, 59, -218, 22),
    (-263, 43, -241, 26), (-447, 67, -433, 94), (-1845, 57, -1838, 39),
    (-2052, 71, -2078, 63), (-2509, 46, -2510, 43), (-2764, 72, -2718, 40),
    (-2947, 62, -2948, 75), (-3127, 76, -3150, 70), (-3406, 70, -3394, 83),
    (-3844, 104, -3854, 78), (-4087, 56, -4090, 63), (-4626, 65, -4630, 51),
    (-4852, 48, -4863, 70), (-6130, 49, -6133, 104), (-6309, 86, -6280, 80),
]

EVENT_CATALOG_SHA256 = "a61021bd9588a58f5da65544b7c0bdf0273e50fed12a47586afe3546e960694b"
ALIGNED_SMP_SHA256 = {
    "smp_10be_aligned.dat": "563da1692bfdf6279370cfcf1c57038fa01f364d6a36b785707fd05ab2e74215",
    "smp_14c_aligned.dat": "404219e6723480b7b9166e2185874c33954e63c97c9c9dbc7a9598cae2fe9538",
}

def _event_tuples(analysis: StudyAnalysis, kind: str) -> list[tuple[int, int, int, int]]:
    pairs = analysis.classification.minima if kind == "minimum" else analysis.classification.maxima
    return [tuple(row.values()) for row in table_event_rows(pairs)]


def _rounded(value: float) -> int:
    return int(np.floor(value + 0.5))


def validate(
    analysis: StudyAnalysis,
    repository_root: Path,
    output_dir: Path,
) -> dict[str, object]:
    classification = analysis.classification
    exact_checks = {
        "manual_catalog_has_160_pairs": len(classification.pairs) == 160,
        "manual_catalog_has_75_dips": sum(pair.be10.sign == 0 for pair in classification.pairs) == 75,
        "manual_catalog_has_85_peaks": sum(pair.be10.sign == 1 for pair in classification.pairs) == 85,
        "table2_all_32_rows_exact": _event_tuples(analysis, "minimum") == PUBLISHED_MINIMA,
        "table3_all_21_rows_exact": _event_tuples(analysis, "maximum") == PUBLISHED_MAXIMA,
        "thresholds_round_to_published": (
            np.allclose(classification.be10_mixture.means, [-0.92, 1.35], atol=0.005)
            and np.allclose(classification.c14_mixture.means, [-0.67, 1.41], atol=0.005)
        ),
        "minimum_duration_is_inclusive_22": (
            len(classification.minima) == 32 and len(classification.maxima) == 21
        ),
        "table4_rounds_exactly": [
            (
                f"{result.power_law.alpha:.2f}",
                _rounded(result.power_law.criteria.bic),
                _rounded(result.power_law.criteria.aic),
                _rounded(result.exponential.tau),
                _rounded(result.exponential.criteria.bic),
                _rounded(result.exponential.criteria.aic),
            )
            for result in [analysis.minimum_waiting, analysis.maximum_waiting]
        ] == [("2.30", 648, 646, 255, 815, 813), ("2.45", 443, 441, 395, 562, 560)],
        "time_fraction_numerators_exact": (
            sum(pair.be10.duration_years for pair in classification.minima) == 2254
            and sum(pair.c14.duration_years for pair in classification.minima) == 2312
            and sum(pair.be10.duration_years for pair in classification.maxima) == 1307
            and sum(pair.c14.duration_years for pair in classification.maxima) == 1347
        ),
        "all_published_grand_pairs_overlap": all(
            min(pair.be10.end_bp, pair.c14.end_bp) - max(pair.be10.onset_bp, pair.c14.onset_bp) > 0
            for pair in classification.minima + classification.maxima
        ),
        "all_ks_comparisons_significant_at_99_percent": all(
            pvalue < 0.01
            for nuclide in analysis.ks_tests.values()
            for comparison in nuclide.values()
            for pvalue in asdict(comparison).values()
        ),
        "event_catalog_checksum": sha256(repository_root / "catalog" / "manual_event_catalog.csv")
        == EVENT_CATALOG_SHA256,
        "aligned_smp_checksums": all(
            sha256(repository_root / "data" / "derived" / filename) == expected
            for filename, expected in ALIGNED_SMP_SHA256.items()
        ),
        "aligned_smp_shared_published_grid": (
            np.array_equal(analysis.be10_signal.age_bp, analysis.c14_signal.age_bp)
            and np.array_equal(
                analysis.c14_signal.age_bp, np.arange(305.0, 8553.0, 1.0)
            )
        ),
    }

    documented_differences = [
        {
            "location": "Table 1, 14C one-Gaussian AIC",
            "calculated": _rounded(
                __import__("grand_solar_reproduction.models", fromlist=["fit_gaussian"])
                .fit_gaussian(np.array([pair.c14.amplitude_sigma for pair in classification.pairs]))
                .criteria.aic
            ),
            "published": 540,
            "explanation": "The conventional likelihood calculation gives 539; the paper prints 540.",
        },
        {
            "location": "Table 5, grand-maximum lognormal BIC",
            "calculated": _rounded(analysis.maximum_durations.lognormal.criteria.bic),
            "published": 384,
            "explanation": "The input durations give 384.55, which rounds to 385; the paper prints 384.",
        },
    ]
    failures = [name for name, passed in exact_checks.items() if not passed]
    report: dict[str, object] = {
        "status": "passed" if not failures else "failed",
        "exact_checks": exact_checks,
        "failures": failures,
        "documented_numerical_differences": documented_differences,
        "computed": {
            "thresholds": classification.thresholds,
            "grand_minima": len(classification.minima),
            "grand_maxima": len(classification.maxima),
            "fraction_of_8250_years": {
                "10Be_minimum": 2254 / 8250,
                "14C_minimum": 2312 / 8250,
                "10Be_maximum": 1307 / 8250,
                "14C_maximum": 1347 / 8250,
            },
            "ks_p_values": {
                nuclide: {metric: asdict(comparison) for metric, comparison in metrics.items()}
                for nuclide, metrics in analysis.ks_tests.items()
            },
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "validation_report.json", report)
    markdown = [
        "# Numerical validation",
        "",
        f"Overall status: **{report['status']}**.",
        "",
        "## Exact checks",
        "",
    ]
    markdown.extend(f"- [{'x' if passed else ' '}] {name}" for name, passed in exact_checks.items())
    markdown.extend(["", "## Documented differences", ""])
    markdown.extend(
        f"- {item['location']}: calculated {item['calculated']}, published {item['published']}. "
        f"{item['explanation']}"
        for item in documented_differences
    )
    (output_dir / "validation_report.md").write_text("\n".join(markdown) + "\n", encoding="utf-8")
    if failures:
        raise AssertionError("Validation failed: " + ", ".join(failures))
    return report
