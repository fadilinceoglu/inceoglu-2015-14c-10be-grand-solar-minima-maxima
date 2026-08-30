"""Small, explicit statistical models used in the published analysis."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


LOG_2PI = float(np.log(2.0 * np.pi))


@dataclass(frozen=True)
class FitCriteria:
    log_likelihood: float
    parameters: int
    sample_size: int

    @property
    def aic(self) -> float:
        return 2.0 * self.parameters - 2.0 * self.log_likelihood

    @property
    def bic(self) -> float:
        return self.parameters * np.log(self.sample_size) - 2.0 * self.log_likelihood


@dataclass(frozen=True)
class GaussianFit:
    mean: float
    sigma: float
    criteria: FitCriteria

    def pdf(self, values: np.ndarray) -> np.ndarray:
        z = (np.asarray(values) - self.mean) / self.sigma
        return np.exp(-0.5 * z * z) / (self.sigma * np.sqrt(2.0 * np.pi))


@dataclass(frozen=True)
class GaussianMixtureFit:
    weights: np.ndarray
    means: np.ndarray
    sigmas: np.ndarray
    criteria: FitCriteria

    def pdf(self, values: np.ndarray) -> np.ndarray:
        x = np.asarray(values, dtype=float)[:, None]
        z = (x - self.means[None, :]) / self.sigmas[None, :]
        components = np.exp(-0.5 * z * z) / (
            self.sigmas[None, :] * np.sqrt(2.0 * np.pi)
        )
        return np.sum(components * self.weights[None, :], axis=1)


@dataclass(frozen=True)
class LognormalFit:
    log_mean: float
    log_sigma: float
    criteria: FitCriteria

    @property
    def mean(self) -> float:
        return float(np.exp(self.log_mean + 0.5 * self.log_sigma**2))


@dataclass(frozen=True)
class ExponentialFit:
    tau: float
    criteria: FitCriteria


@dataclass(frozen=True)
class PowerLawFit:
    alpha: float
    xmin: float
    tail_size: int
    criteria: FitCriteria


def fit_gaussian(values: np.ndarray) -> GaussianFit:
    x = np.asarray(values, dtype=float)
    mean = float(np.mean(x))
    sigma = float(np.sqrt(np.mean((x - mean) ** 2)))
    if sigma <= 0:
        raise ValueError("A Gaussian variance must be positive")
    log_likelihood = float(
        np.sum(-0.5 * LOG_2PI - np.log(sigma) - 0.5 * ((x - mean) / sigma) ** 2)
    )
    return GaussianFit(mean, sigma, FitCriteria(log_likelihood, 2, x.size))


def _normal_log_density(x: np.ndarray, means: np.ndarray, sigmas: np.ndarray) -> np.ndarray:
    z = (x[:, None] - means[None, :]) / sigmas[None, :]
    return -0.5 * LOG_2PI - np.log(sigmas)[None, :] - 0.5 * z * z


def _logsumexp(values: np.ndarray, axis: int) -> np.ndarray:
    maximum = np.max(values, axis=axis, keepdims=True)
    return np.squeeze(maximum, axis=axis) + np.log(
        np.sum(np.exp(values - maximum), axis=axis)
    )


def _em(
    x: np.ndarray,
    weights: np.ndarray,
    means: np.ndarray,
    sigmas: np.ndarray,
    *,
    tolerance: float = 1e-12,
    max_iterations: int = 100_000,
) -> GaussianMixtureFit:
    previous = -np.inf
    variance_floor = max(float(np.var(x)) * 1e-12, 1e-12)
    for _ in range(max_iterations):
        log_weighted = _normal_log_density(x, means, sigmas) + np.log(weights)[None, :]
        log_norm = _logsumexp(log_weighted, axis=1)
        log_likelihood = float(np.sum(log_norm))
        responsibilities = np.exp(log_weighted - log_norm[:, None])
        counts = np.sum(responsibilities, axis=0)
        weights = counts / x.size
        means = np.sum(responsibilities * x[:, None], axis=0) / counts
        variances = (
            np.sum(responsibilities * (x[:, None] - means[None, :]) ** 2, axis=0)
            / counts
        )
        sigmas = np.sqrt(np.maximum(variances, variance_floor))
        if abs(log_likelihood - previous) <= tolerance * (1.0 + abs(log_likelihood)):
            break
        previous = log_likelihood
    log_weighted = _normal_log_density(x, means, sigmas) + np.log(weights)[None, :]
    log_likelihood = float(np.sum(_logsumexp(log_weighted, axis=1)))
    order = np.argsort(means)
    return GaussianMixtureFit(
        weights[order], means[order], sigmas[order], FitCriteria(log_likelihood, 5, x.size)
    )


def fit_gaussian_mixture(
    values: np.ndarray,
    *,
    initialization: str = "global",
) -> GaussianMixtureFit:
    """Fit a two-component univariate Gaussian mixture by maximum likelihood.

    ``global`` evaluates several deterministic starts and keeps the greatest
    likelihood. ``collapsed`` starts two equal components at the one-Gaussian
    solution. The latter reproduces the documented local solution used for
    grand-maximum durations in Table 5; it is not used elsewhere.
    """

    x = np.asarray(values, dtype=float)
    one = fit_gaussian(x)
    if initialization == "collapsed":
        return _em(
            x,
            np.array([0.5, 0.5]),
            np.array([one.mean, one.mean]),
            np.array([one.sigma, one.sigma]),
        )
    if initialization != "global":
        raise ValueError(f"Unknown GMM initialization: {initialization}")

    starts: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for lower, upper in [(15, 85), (25, 75), (35, 65), (5, 95)]:
        means = np.percentile(x, [lower, upper])
        starts.append((np.array([0.5, 0.5]), means, np.array([one.sigma, one.sigma])))
    negative = x < 0
    if np.any(negative) and np.any(~negative):
        counts = np.array([np.sum(negative), np.sum(~negative)], dtype=float)
        means = np.array([np.mean(x[negative]), np.mean(x[~negative])])
        sigmas = np.array([np.std(x[negative]), np.std(x[~negative])])
        starts.append((counts / x.size, means, np.maximum(sigmas, one.sigma * 0.05)))
    fits = [_em(x, *start) for start in starts]
    return max(fits, key=lambda fit: fit.criteria.log_likelihood)


def fit_lognormal(values: np.ndarray) -> LognormalFit:
    x = np.asarray(values, dtype=float)
    if np.any(x <= 0):
        raise ValueError("Lognormal observations must be positive")
    log_x = np.log(x)
    mean = float(np.mean(log_x))
    sigma = float(np.sqrt(np.mean((log_x - mean) ** 2)))
    log_likelihood = float(
        np.sum(
            -np.log(x)
            - np.log(sigma)
            - 0.5 * LOG_2PI
            - 0.5 * ((log_x - mean) / sigma) ** 2
        )
    )
    return LognormalFit(mean, sigma, FitCriteria(log_likelihood, 2, x.size))


def fit_exponential(values: np.ndarray) -> ExponentialFit:
    x = np.asarray(values, dtype=float)
    if np.any(x < 0):
        raise ValueError("Waiting times cannot be negative")
    tau = float(np.mean(x))
    log_likelihood = float(-x.size * np.log(tau) - np.sum(x) / tau)
    return ExponentialFit(tau, FitCriteria(log_likelihood, 1, x.size))


def fit_power_law(values: np.ndarray, xmin: float) -> PowerLawFit:
    """Continuous Clauset MLE with the finite-sample alpha correction.

    The likelihood and BIC sample size are evaluated on the fitted tail.
    """

    full = np.asarray(values, dtype=float)
    tail = full[full >= xmin]
    if tail.size < 2:
        raise ValueError("Power-law tail is too small")
    mle_alpha = 1.0 + tail.size / np.sum(np.log(tail / xmin))
    alpha = 1.0 + (tail.size - 1.0) / tail.size * (mle_alpha - 1.0)
    log_likelihood = float(
        tail.size * np.log(alpha - 1.0)
        + tail.size * (alpha - 1.0) * np.log(xmin)
        - alpha * np.sum(np.log(tail))
    )
    return PowerLawFit(
        float(alpha),
        float(xmin),
        int(tail.size),
        FitCriteria(log_likelihood, 1, tail.size),
    )


def empirical_cdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.asarray(values, dtype=float))
    return x, np.arange(1, x.size + 1, dtype=float) / x.size


def empirical_ccdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.asarray(values, dtype=float))
    return x, np.arange(x.size, 0, -1, dtype=float) / x.size
