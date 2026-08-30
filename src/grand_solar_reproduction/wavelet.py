"""Morlet continuous and cross-wavelet calculations after Grinsted et al. (2004)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import chi2


MORLET_OMEGA0 = 6.0
FOURIER_FACTOR = 4.0 * np.pi / (MORLET_OMEGA0 + np.sqrt(2.0 + MORLET_OMEGA0**2))


@dataclass(frozen=True)
class WaveletResult:
    coefficients: np.ndarray
    power: np.ndarray
    scales: np.ndarray
    periods: np.ndarray
    cone_of_influence: np.ndarray


@dataclass(frozen=True)
class CrossWaveletResult:
    cross: np.ndarray
    power: np.ndarray
    periods: np.ndarray
    cone_of_influence: np.ndarray
    significance_ratio: np.ndarray


def lag_one_autocorrelation(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    result = float(np.corrcoef(x[:-1], x[1:])[0, 1])
    return float(np.clip(result, -0.99, 0.99))


def morlet_cwt(
    values: np.ndarray,
    *,
    dt: float = 1.0,
    dj: float = 1.0 / 12.0,
    smallest_scale: float = 2.0,
) -> WaveletResult:
    """Torrence--Compo Fourier-domain Morlet CWT used by Grinsted's code."""

    x = np.asarray(values, dtype=float)
    x = x - np.mean(x)
    n = x.size
    padded_n = 1 << int(np.ceil(np.log2(n)))
    padded = np.zeros(padded_n, dtype=float)
    padded[:n] = x
    frequencies = 2.0 * np.pi * np.fft.fftfreq(padded_n, d=dt)
    positive = frequencies > 0
    spectrum = np.fft.fft(padded)
    maximum_j = int(np.floor(np.log2(n * dt / smallest_scale) / dj))
    scales = smallest_scale * 2.0 ** (np.arange(maximum_j + 1) * dj)
    coefficients = np.empty((scales.size, n), dtype=complex)
    normalization_frequency = 2.0 * np.pi / (padded_n * dt)
    for index, scale in enumerate(scales):
        daughter = np.zeros(padded_n, dtype=float)
        daughter[positive] = (
            np.pi ** (-0.25)
            * np.sqrt(scale * normalization_frequency)
            * np.sqrt(padded_n)
            * np.exp(-0.5 * (scale * frequencies[positive] - MORLET_OMEGA0) ** 2)
        )
        coefficients[index] = np.fft.ifft(spectrum * daughter)[:n]
    periods = scales * FOURIER_FACTOR
    edge_distance = np.minimum(np.arange(n) + 1, np.arange(n, 0, -1))
    cone = FOURIER_FACTOR / np.sqrt(2.0) * dt * edge_distance
    return WaveletResult(coefficients, np.abs(coefficients) ** 2, scales, periods, cone)


def red_noise_spectrum(periods: np.ndarray, autocorrelation: float, dt: float = 1.0) -> np.ndarray:
    frequency_angle = 2.0 * np.pi * dt / periods
    return (1.0 - autocorrelation**2) / (
        1.0 + autocorrelation**2 - 2.0 * autocorrelation * np.cos(frequency_angle)
    )


def cross_wavelet(
    first: np.ndarray,
    second: np.ndarray,
    *,
    dt: float = 1.0,
    dj: float = 1.0 / 12.0,
) -> CrossWaveletResult:
    left = morlet_cwt(first, dt=dt, dj=dj)
    right = morlet_cwt(second, dt=dt, dj=dj)
    cross = left.coefficients * np.conjugate(right.coefficients)
    power = np.abs(cross)
    background_left = red_noise_spectrum(
        left.periods, lag_one_autocorrelation(first), dt=dt
    )
    background_right = red_noise_spectrum(
        right.periods, lag_one_autocorrelation(second), dt=dt
    )
    significance = (
        np.sqrt(background_left * background_right) * chi2.ppf(0.95, 2) / 2.0
    )
    ratio = power / significance[:, None]
    return CrossWaveletResult(cross, power, left.periods, left.cone_of_influence, ratio)


def band_power(
    values: np.ndarray,
    low_period_years: float = 20.0,
    high_period_years: float = 40.0,
) -> tuple[np.ndarray, WaveletResult]:
    result = morlet_cwt(values)
    selection = (result.periods >= low_period_years) & (result.periods <= high_period_years)
    if not np.any(selection):
        raise ValueError("Requested wavelet period band is absent")
    return np.mean(result.power[selection], axis=0), result

