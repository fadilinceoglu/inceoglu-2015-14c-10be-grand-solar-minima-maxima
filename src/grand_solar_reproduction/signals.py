"""Documented signal transformations used by Figures 1 and 5--8."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

from .events import EventPair


ANALYSIS_START_BP = 305.0
ANALYSIS_END_BP = 8552.0
POLYNOMIAL_DEGREE = 5
HIGH_PASS_CUTOFF_PER_YEAR = 1.0 / 30.0
HIGH_PASS_ORDER = 10
MOVING_VARIANCE_YEARS = 25


@dataclass(frozen=True)
class PreparedSignal:
    age_bp: np.ndarray
    calendar_year: np.ndarray
    phi: np.ndarray
    trend: np.ndarray
    residual: np.ndarray
    standardized: np.ndarray
    high_pass: np.ndarray
    moving_variance: np.ndarray


def trim_analysis_period(
    age_bp: np.ndarray,
    values: np.ndarray,
    start_bp: float = ANALYSIS_START_BP,
    end_bp: float = ANALYSIS_END_BP,
) -> tuple[np.ndarray, np.ndarray]:
    mask = (age_bp >= start_bp) & (age_bp <= end_bp)
    age = np.asarray(age_bp[mask], dtype=float)
    value = np.asarray(values[mask], dtype=float)
    if age.size == 0:
        raise ValueError("No observations fall in the study analysis period")
    if not np.allclose(np.diff(age), 1.0):
        target = np.arange(np.ceil(age[0]), np.floor(age[-1]) + 1.0, 1.0)
        value = np.interp(target, age, value)
        age = target
    return age, value


def polynomial_detrend(
    age_bp: np.ndarray, values: np.ndarray, degree: int = POLYNOMIAL_DEGREE
) -> tuple[np.ndarray, np.ndarray]:
    # Center and scale time before fitting. The polynomial is mathematically
    # identical to a degree-five fit in age, but the Vandermonde system is far
    # better conditioned and therefore deterministic across modern platforms.
    center = float(np.mean(age_bp))
    scale = float(np.ptp(age_bp) / 2.0)
    coordinate = (age_bp - center) / scale
    coefficients = np.polynomial.Polynomial.fit(coordinate, values, degree).convert().coef
    trend = np.polynomial.polynomial.polyval(coordinate, coefficients)
    return trend, values - trend


def standardize(values: np.ndarray) -> np.ndarray:
    centered = values - np.mean(values)
    sigma = np.sqrt(np.mean(centered**2))
    if sigma <= 0:
        raise ValueError("Cannot standardize a constant signal")
    return centered / sigma


def high_pass_filter(
    values: np.ndarray,
    cutoff_per_year: float = HIGH_PASS_CUTOFF_PER_YEAR,
    order: int = HIGH_PASS_ORDER,
) -> np.ndarray:
    """Zero-phase, tenth-order Butterworth high-pass filter on annual data.

    This implementation declares degree 10. With ``fs=1``, the physical cutoff
    1/30 yr^-1 is equivalent to the Nyquist-normalized value ``Wn=2/30``.
    """

    sos = signal.butter(order, cutoff_per_year, btype="highpass", fs=1.0, output="sos")
    return signal.sosfiltfilt(sos, np.asarray(values, dtype=float))


def centered_moving_variance(values: np.ndarray, window: int = MOVING_VARIANCE_YEARS) -> np.ndarray:
    if window < 2 or window % 2 == 0:
        raise ValueError("The moving-variance window must be an odd integer >= 3")
    values = np.asarray(values, dtype=float)
    kernel = np.ones(window, dtype=float) / window
    mean = np.convolve(values, kernel, mode="valid")
    mean_square = np.convolve(values * values, kernel, mode="valid")
    variance = mean_square - mean * mean
    pad = window // 2
    return np.pad(variance, (pad, pad), constant_values=np.nan)


def prepare_signal(
    age_bp: np.ndarray, phi: np.ndarray, target_age_bp: np.ndarray | None = None
) -> PreparedSignal:
    if target_age_bp is None:
        age, values = trim_analysis_period(age_bp, phi)
    else:
        age = np.asarray(target_age_bp, dtype=float)
        if age[0] < age_bp[0] or age[-1] > age_bp[-1]:
            raise ValueError("Target grid would extrapolate beyond the source chronology")
        values = np.interp(age, age_bp, phi)
    trend, residual = polynomial_detrend(age, values)
    standardized = standardize(residual)
    high_pass = high_pass_filter(standardized)
    moving_variance = centered_moving_variance(high_pass)
    return PreparedSignal(
        age,
        1950.0 - age,
        values,
        trend,
        residual,
        standardized,
        high_pass,
        moving_variance,
    )


def activity_masks(
    age_bp: np.ndarray,
    minima: list[EventPair],
    maxima: list[EventPair],
    nuclide: str,
) -> dict[str, np.ndarray]:
    if nuclide not in {"10Be", "14C"}:
        raise ValueError("nuclide must be '10Be' or '14C'")
    minimum_mask = np.zeros(age_bp.size, dtype=bool)
    maximum_mask = np.zeros(age_bp.size, dtype=bool)
    attribute = "be10" if nuclide == "10Be" else "c14"
    for pair in minima:
        event = getattr(pair, attribute)
        minimum_mask |= (age_bp >= event.onset_bp) & (age_bp <= event.end_bp)
    for pair in maxima:
        event = getattr(pair, attribute)
        maximum_mask |= (age_bp >= event.onset_bp) & (age_bp <= event.end_bp)
    return {
        "minimum": minimum_mask,
        "maximum": maximum_mask,
        "moderate": ~(minimum_mask | maximum_mask),
    }
