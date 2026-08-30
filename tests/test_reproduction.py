from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from grand_solar_reproduction.analysis import run_analysis
from grand_solar_reproduction.events import table_event_rows
from grand_solar_reproduction.validation import PUBLISHED_MAXIMA, PUBLISHED_MINIMA


ROOT = Path(__file__).resolve().parents[1]
ALIGNED_DATA = ROOT / "data" / "derived"


@pytest.fixture(scope="module")
def analysis():
    return run_analysis(ALIGNED_DATA, ROOT / "catalog" / "manual_event_catalog.csv")


def test_aligned_smp_inputs_share_published_grid():
    be10 = np.loadtxt(ALIGNED_DATA / "smp_10be_aligned.dat")
    c14 = np.loadtxt(ALIGNED_DATA / "smp_14c_aligned.dat")
    assert be10.shape == c14.shape == (8248, 2)
    assert np.array_equal(be10[:, 0], c14[:, 0])
    assert np.array_equal(be10[:, 0], np.arange(305.0, 8553.0))


def test_all_published_event_rows_match(analysis):
    minima = [tuple(row.values()) for row in table_event_rows(analysis.classification.minima)]
    maxima = [tuple(row.values()) for row in table_event_rows(analysis.classification.maxima)]
    assert minima == PUBLISHED_MINIMA
    assert maxima == PUBLISHED_MAXIMA


def test_mixture_thresholds_match_paper(analysis):
    assert np.allclose(analysis.classification.be10_mixture.means, [-0.92, 1.35], atol=0.005)
    assert np.allclose(analysis.classification.c14_mixture.means, [-0.67, 1.41], atol=0.005)


def test_table4_matches_paper(analysis):
    assert analysis.minimum_waiting.power_law.alpha == pytest.approx(2.30, abs=0.005)
    assert analysis.maximum_waiting.power_law.alpha == pytest.approx(2.45, abs=0.005)
    assert analysis.minimum_waiting.exponential.tau == pytest.approx(255, abs=0.5)
    assert analysis.maximum_waiting.exponential.tau == pytest.approx(395, abs=0.5)


def test_activity_fractions_match_paper(analysis):
    assert sum(pair.be10.duration_years for pair in analysis.classification.minima) == 2254
    assert sum(pair.c14.duration_years for pair in analysis.classification.minima) == 2312
    assert sum(pair.be10.duration_years for pair in analysis.classification.maxima) == 1307
    assert sum(pair.c14.duration_years for pair in analysis.classification.maxima) == 1347


def test_ks_conclusions_match_paper(analysis):
    for metrics in analysis.ks_tests.values():
        for comparison in metrics.values():
            assert comparison.minimum_vs_moderate < 0.01
            assert comparison.maximum_vs_moderate < 0.01
            assert comparison.minimum_vs_maximum < 0.01


def test_all_eight_figures_exist():
    figures = sorted((ROOT / "outputs" / "figures").glob("figure*.pdf"))
    assert len(figures) == 8
    assert all(path.stat().st_size > 10_000 for path in figures)


def test_one_output_format_per_artifact_type():
    assert not list((ROOT / "outputs" / "figures").glob("*.png"))
    assert not list((ROOT / "outputs" / "tables").glob("*.md"))
    assert len(list((ROOT / "outputs" / "tables").glob("*.csv"))) == 5
