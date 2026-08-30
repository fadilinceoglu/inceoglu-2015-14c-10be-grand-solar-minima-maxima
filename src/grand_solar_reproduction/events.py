"""The final visual catalog and the paper's event-selection calculation."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .models import GaussianMixtureFit, fit_gaussian_mixture


@dataclass(frozen=True)
class VisualEvent:
    onset_bp: float
    end_bp: float
    center_bp: float
    duration_years: float
    amplitude_sigma: float
    sign: int


@dataclass(frozen=True)
class EventPair:
    event_index: int
    be10: VisualEvent
    c14: VisualEvent


@dataclass(frozen=True)
class Classification:
    pairs: list[EventPair]
    minima: list[EventPair]
    maxima: list[EventPair]
    be10_mixture: GaussianMixtureFit
    c14_mixture: GaussianMixtureFit

    @property
    def thresholds(self) -> dict[str, tuple[float, float]]:
        return {
            "10Be": tuple(float(value) for value in self.be10_mixture.means),
            "14C": tuple(float(value) for value in self.c14_mixture.means),
        }


def load_event_catalog(path: Path) -> list[EventPair]:
    """Load the released CSV event catalogue."""

    pairs: list[EventPair] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            sign = int(row["sign"])
            be10 = VisualEvent(
                float(row["be10_onset_bp"]),
                float(row["be10_end_bp"]),
                float(row["be10_center_bp"]),
                float(row["be10_duration_years"]),
                float(row["be10_amplitude_sigma"]),
                sign,
            )
            c14 = VisualEvent(
                float(row["c14_onset_bp"]),
                float(row["c14_end_bp"]),
                float(row["c14_center_bp"]),
                float(row["c14_duration_years"]),
                float(row["c14_amplitude_sigma"]),
                sign,
            )
            pairs.append(EventPair(int(row["event_index"]), be10, c14))
    if [pair.event_index for pair in pairs] != list(range(1, 161)):
        raise ValueError("Event catalog must contain event_index 1 through 160 in order")
    return pairs


def classify_events(pairs: list[EventPair], minimum_duration: float = 22.0) -> Classification:
    be10_mixture = fit_gaussian_mixture(
        np.array([pair.be10.amplitude_sigma for pair in pairs]), initialization="global"
    )
    c14_mixture = fit_gaussian_mixture(
        np.array([pair.c14.amplitude_sigma for pair in pairs]), initialization="global"
    )
    be_low, be_high = be10_mixture.means
    c_low, c_high = c14_mixture.means

    def durations_qualify(pair: EventPair) -> bool:
        # The paper says "longer than 22 years", but both published lists
        # contain a 22-year event and therefore require inclusive >= 22.
        return (
            pair.be10.duration_years >= minimum_duration
            and pair.c14.duration_years >= minimum_duration
        )

    minima = [
        pair
        for pair in pairs
        if pair.be10.sign == 0
        and pair.be10.amplitude_sigma < be_low
        and pair.c14.amplitude_sigma < c_low
        and durations_qualify(pair)
    ]
    maxima = [
        pair
        for pair in pairs
        if pair.be10.sign == 1
        and pair.be10.amplitude_sigma > be_high
        and pair.c14.amplitude_sigma > c_high
        and durations_qualify(pair)
    ]
    return Classification(pairs, minima, maxima, be10_mixture, c14_mixture)


def calendar_year(center_bp: float) -> int:
    return int(np.rint(1950.0 - center_bp))


def rounded_duration(duration: float) -> int:
    return int(np.rint(duration))


def table_event_rows(pairs: list[EventPair]) -> list[dict[str, object]]:
    return [
        {
            "be10_center_ad_bc": calendar_year(pair.be10.center_bp),
            "be10_duration_years": rounded_duration(pair.be10.duration_years),
            "c14_center_ad_bc": calendar_year(pair.c14.center_bp),
            "c14_duration_years": rounded_duration(pair.c14.duration_years),
        }
        for pair in pairs
    ]
