"""Assemble the paper's downstream calculations from committed aligned SMP inputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.stats import ks_2samp

from .events import Classification, EventPair, classify_events, load_event_catalog
from .io import load_two_column_series
from .models import (
    ExponentialFit,
    GaussianFit,
    GaussianMixtureFit,
    LognormalFit,
    PowerLawFit,
    fit_exponential,
    fit_gaussian,
    fit_gaussian_mixture,
    fit_lognormal,
    fit_power_law,
)
from .signals import PreparedSignal, activity_masks, prepare_signal
from .wavelet import CrossWaveletResult, band_power, cross_wavelet


# The lower-tail bounds are explicit analysis parameters. The minimum bound is
# also the Clauset KS optimum. For maxima, 180 reproduces the published
# alpha/AIC/BIC values, whereas an unconstrained continuous KS search selects
# 147; this known difference is documented in the reproduction notes.
PUBLISHED_POWER_LAW_XMIN = {"minimum": 112.0, "maximum": 180.0}


@dataclass(frozen=True)
class WaitingTimeAnalysis:
    be10: np.ndarray
    c14: np.ndarray
    pooled: np.ndarray
    exponential: ExponentialFit
    power_law: PowerLawFit


@dataclass(frozen=True)
class DurationAnalysis:
    pooled: np.ndarray
    gaussian: GaussianFit
    gaussian_mixture: GaussianMixtureFit
    lognormal: LognormalFit


@dataclass(frozen=True)
class DistributionComparison:
    minimum_vs_moderate: float
    maximum_vs_moderate: float
    minimum_vs_maximum: float


@dataclass(frozen=True)
class StudyAnalysis:
    classification: Classification
    be10_signal: PreparedSignal
    c14_signal: PreparedSignal
    minimum_waiting: WaitingTimeAnalysis
    maximum_waiting: WaitingTimeAnalysis
    minimum_durations: DurationAnalysis
    maximum_durations: DurationAnalysis
    cross_wavelet: CrossWaveletResult
    be10_band_power: np.ndarray
    c14_band_power: np.ndarray
    distribution_samples: dict[str, dict[str, dict[str, np.ndarray]]]
    ks_tests: dict[str, dict[str, DistributionComparison]]


def _waiting_times(pairs: list[EventPair]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    be10_centers = np.sort(np.array([pair.be10.center_bp for pair in pairs]))
    c14_centers = np.sort(np.array([pair.c14.center_bp for pair in pairs]))
    be10 = np.diff(be10_centers)
    c14 = np.diff(c14_centers)
    return be10, c14, np.concatenate([be10, c14])


def _waiting_analysis(pairs: list[EventPair], kind: str) -> WaitingTimeAnalysis:
    be10, c14, pooled = _waiting_times(pairs)
    return WaitingTimeAnalysis(
        be10,
        c14,
        pooled,
        fit_exponential(pooled),
        fit_power_law(pooled, PUBLISHED_POWER_LAW_XMIN[kind]),
    )


def _duration_analysis(pairs: list[EventPair], kind: str) -> DurationAnalysis:
    pooled = np.array(
        [duration for pair in pairs for duration in (pair.be10.duration_years, pair.c14.duration_years)],
        dtype=float,
    )
    # The published grand-maximum AIC/BIC values correspond to the collapsed
    # local EM solution. The minimum-duration fit uses the global solution.
    initialization = "collapsed" if kind == "maximum" else "global"
    return DurationAnalysis(
        pooled,
        fit_gaussian(pooled),
        fit_gaussian_mixture(pooled, initialization=initialization),
        fit_lognormal(pooled),
    )


def _samples_and_ks(
    classification: Classification,
    be10_signal: PreparedSignal,
    c14_signal: PreparedSignal,
    be10_power: np.ndarray,
    c14_power: np.ndarray,
) -> tuple[
    dict[str, dict[str, dict[str, np.ndarray]]],
    dict[str, dict[str, DistributionComparison]],
]:
    samples: dict[str, dict[str, dict[str, np.ndarray]]] = {}
    tests: dict[str, dict[str, DistributionComparison]] = {}
    for nuclide, prepared, power in [
        ("10Be", be10_signal, be10_power),
        ("14C", c14_signal, c14_power),
    ]:
        masks = activity_masks(
            prepared.age_bp, classification.minima, classification.maxima, nuclide
        )
        samples[nuclide] = {"variance": {}, "power": {}}
        tests[nuclide] = {}
        for metric, values in [("variance", prepared.moving_variance), ("power", power)]:
            finite = np.isfinite(values)
            for state, mask in masks.items():
                samples[nuclide][metric][state] = values[mask & finite]
            minimum = samples[nuclide][metric]["minimum"]
            maximum = samples[nuclide][metric]["maximum"]
            moderate = samples[nuclide][metric]["moderate"]
            tests[nuclide][metric] = DistributionComparison(
                float(ks_2samp(minimum, moderate).pvalue),
                float(ks_2samp(maximum, moderate).pvalue),
                float(ks_2samp(minimum, maximum).pvalue),
            )
    return samples, tests


def run_analysis(aligned_data: Path, event_catalog: Path) -> StudyAnalysis:
    if not event_catalog.is_file():
        raise FileNotFoundError(f"Visual-event catalog is absent: {event_catalog}")
    pairs = load_event_catalog(event_catalog)
    classification = classify_events(pairs)

    age_c14, phi_c14 = load_two_column_series(aligned_data / "smp_14c_aligned.dat")
    age_be10, phi_be10 = load_two_column_series(aligned_data / "smp_10be_aligned.dat")
    if not np.array_equal(age_c14, age_be10):
        raise ValueError("The two committed SMP inputs must share an identical age grid")
    c14_signal = prepare_signal(age_c14, phi_c14)
    be10_signal = prepare_signal(age_be10, phi_be10, target_age_bp=c14_signal.age_bp)

    xwt = cross_wavelet(be10_signal.high_pass, c14_signal.high_pass)
    be10_power, _ = band_power(be10_signal.high_pass)
    c14_power, _ = band_power(c14_signal.high_pass)
    samples, tests = _samples_and_ks(
        classification, be10_signal, c14_signal, be10_power, c14_power
    )

    return StudyAnalysis(
        classification,
        be10_signal,
        c14_signal,
        _waiting_analysis(classification.minima, "minimum"),
        _waiting_analysis(classification.maxima, "maximum"),
        _duration_analysis(classification.minima, "minimum"),
        _duration_analysis(classification.maxima, "maximum"),
        xwt,
        be10_power,
        c14_power,
        samples,
        tests,
    )
