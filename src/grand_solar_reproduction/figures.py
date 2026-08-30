"""Generate modern versions of the paper's eight published figures."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "inceoglu2015-matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

from .analysis import StudyAnalysis, WaitingTimeAnalysis
from .events import EventPair
from .models import empirical_ccdf, empirical_cdf


BE10_COLOR = "#c51b8a"
C14_COLOR = "#238b45"
MINIMUM_COLOR = "#2b4c9b"
MAXIMUM_COLOR = "#d73027"
MODERATE_COLOR = "#777777"
GRID_COLOR = "#d9d9d9"


def configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 220,
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": GRID_COLOR,
            "grid.linewidth": 0.6,
            "legend.frameon": False,
            "lines.linewidth": 1.35,
        }
    )


def _save(figure: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output_dir / f"{stem}.pdf",
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(figure)


def figure1(analysis: StudyAnalysis, output_dir: Path) -> None:
    be10, c14 = analysis.be10_signal, analysis.c14_signal
    figure, axes = plt.subplots(2, 1, figsize=(10.5, 5.8), sharex=True)
    axes[0].plot(be10.calendar_year, be10.phi, color=BE10_COLOR, label=r"$\Phi_{10\mathrm{Be}}$")
    axes[0].plot(c14.calendar_year, c14.phi, color=C14_COLOR, label=r"$\Phi_{14\mathrm{C}}$")
    axes[0].plot(be10.calendar_year, be10.trend, color=BE10_COLOR, linestyle="--", alpha=0.75)
    axes[0].plot(c14.calendar_year, c14.trend, color=C14_COLOR, linestyle="--", alpha=0.75)
    axes[0].set_ylabel("Solar modulation potential (MeV)")
    axes[0].set_title("Solar modulation potential and degree-5 long-term trends")
    axes[0].legend(ncol=2, loc="upper right")

    axes[1].plot(be10.calendar_year, be10.standardized, color=BE10_COLOR, label=r"$^{10}$Be")
    axes[1].plot(c14.calendar_year, c14.standardized, color=C14_COLOR, label=r"$^{14}$C")
    axes[1].axhline(0.0, color="#333333", linewidth=0.8)
    axes[1].set_ylabel("Detrended potential (σ)")
    axes[1].set_xlabel("Calendar year (AD positive, BC negative)")
    axes[1].set_xlim(float(be10.calendar_year[0]), float(c14.calendar_year[-1]))
    figure.tight_layout()
    _save(figure, output_dir, "figure1_solar_modulation_potential")


def figure2(analysis: StudyAnalysis, output_dir: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    for axis, nuclide, attribute, fit, color in [
        (axes[0], r"$^{10}$Be", "be10", analysis.classification.be10_mixture, BE10_COLOR),
        (axes[1], r"$^{14}$C", "c14", analysis.classification.c14_mixture, C14_COLOR),
    ]:
        values = np.array(
            [getattr(pair, attribute).amplitude_sigma for pair in analysis.classification.pairs]
        )
        bins = np.histogram_bin_edges(values, bins="fd")
        axis.hist(values, bins=bins, density=True, color="#d9d9d9", edgecolor="white")
        grid = np.linspace(values.min() - 0.4, values.max() + 0.4, 600)
        axis.plot(grid, fit.pdf(grid), color=color, linewidth=2.2, label="Two-Gaussian mixture")
        for mean, label in zip(fit.means, ["minimum threshold", "maximum threshold"], strict=True):
            axis.axvline(mean, color="#222222", linestyle="--", linewidth=1.0)
            axis.text(mean, axis.get_ylim()[1] * 0.94, f"{mean:.2f}σ", ha="center", va="top")
        axis.set_title(f"Visual events in {nuclide}")
        axis.set_xlabel("Event amplitude (σ)")
        axis.legend(loc="lower center", fontsize=8)
    axes[0].set_ylabel("Probability density")
    figure.tight_layout()
    _save(figure, output_dir, "figure2_event_amplitude_mixtures")


def figure3(analysis: StudyAnalysis, output_dir: Path) -> None:
    minima = {pair.event_index for pair in analysis.classification.minima}
    maxima = {pair.event_index for pair in analysis.classification.maxima}
    figure, axis = plt.subplots(figsize=(8.4, 5.1))
    for attribute, marker, nuclide in [("be10", "D", r"$^{10}$Be"), ("c14", "o", r"$^{14}$C")]:
        for state, indexes, color, order in [
            ("Moderate", None, MODERATE_COLOR, 1),
            ("Grand minimum", minima, MINIMUM_COLOR, 2),
            ("Grand maximum", maxima, MAXIMUM_COLOR, 3),
        ]:
            selected = [
                pair
                for pair in analysis.classification.pairs
                if (indexes is None and pair.event_index not in minima | maxima)
                or (indexes is not None and pair.event_index in indexes)
            ]
            axis.scatter(
                [getattr(pair, attribute).amplitude_sigma for pair in selected],
                [getattr(pair, attribute).duration_years for pair in selected],
                s=28,
                marker=marker,
                facecolor=color,
                edgecolor="white",
                linewidth=0.35,
                alpha=0.88,
                zorder=order,
                label=f"{state}, {nuclide}",
            )
    axis.set_xlabel("Event amplitude (σ)")
    axis.set_ylabel("Duration (years)")
    axis.set_title("Duration and amplitude of all 160 paired visual events")
    axis.legend(ncol=2, fontsize=8)
    figure.tight_layout()
    _save(figure, output_dir, "figure3_event_duration_vs_amplitude")


def _waiting_panel(axis: plt.Axes, result: WaitingTimeAnalysis, title: str) -> None:
    for values, color, label in [
        (result.be10, BE10_COLOR, r"$^{10}$Be"),
        (result.c14, C14_COLOR, r"$^{14}$C"),
    ]:
        x, y = empirical_ccdf(values)
        axis.step(x, y, where="post", color=color, alpha=0.85, label=label)
    grid = np.geomspace(max(1.0, result.pooled.min()), result.pooled.max(), 600)
    axis.plot(grid, np.exp(-grid / result.exponential.tau), color="#f39c12", label="Exponential")
    tail = grid[grid >= result.power_law.xmin]
    axis.plot(
        tail,
        (tail / result.power_law.xmin) ** (1.0 - result.power_law.alpha),
        color="#202020",
        linestyle="--",
        label="Power law",
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_ylim(0.02, 1.1)
    axis.set_title(title)
    axis.set_xlabel("Waiting time (years)")
    axis.set_ylabel("Complementary cumulative probability")


def figure4(analysis: StudyAnalysis, output_dir: Path) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(10.8, 7.7))
    _waiting_panel(axes[0, 0], analysis.minimum_waiting, "Grand-minimum waiting times")
    _waiting_panel(axes[0, 1], analysis.maximum_waiting, "Grand-maximum waiting times")
    axes[0, 0].legend(fontsize=8)
    for axis, result, title in [
        (axes[1, 0], analysis.minimum_durations, "Grand-minimum durations"),
        (axes[1, 1], analysis.maximum_durations, "Grand-maximum durations"),
    ]:
        bins = np.arange(0, 201, 20)
        axis.hist(result.pooled, bins=bins, color="#bdbdbd", edgecolor="white")
        axis.set_xlim(0, 200)
        axis.set_xlabel("Duration (years)")
        axis.set_ylabel("Number per bin")
        axis.set_title(title)
    figure.tight_layout()
    _save(figure, output_dir, "figure4_waiting_times_and_durations")


def _overlay_events(
    axis: plt.Axes,
    analysis: StudyAnalysis,
    pairs: list[EventPair],
    color: str,
) -> None:
    for attribute, prepared in [("be10", analysis.be10_signal), ("c14", analysis.c14_signal)]:
        for pair in pairs:
            event = getattr(pair, attribute)
            mask = (prepared.age_bp >= event.onset_bp) & (prepared.age_bp <= event.end_bp)
            axis.plot(prepared.calendar_year[mask], prepared.standardized[mask], color=color, linewidth=2.3)


def _wavelet_focus(
    analysis: StudyAnalysis,
    output_dir: Path,
    *,
    stem: str,
    title: str,
    pairs: list[EventPair],
    event_color: str,
    left_year: float,
    right_year: float,
) -> None:
    xwt = analysis.cross_wavelet
    years = analysis.be10_signal.calendar_year
    selection = (years <= max(left_year, right_year)) & (years >= min(left_year, right_year))
    figure, axes = plt.subplots(
        2, 1, figsize=(10.6, 6.2), sharex=True, gridspec_kw={"height_ratios": [1, 1.65]}
    )
    axes[0].plot(years, analysis.be10_signal.standardized, color=BE10_COLOR, label=r"$\Phi_{10\mathrm{Be}}$")
    axes[0].plot(years, analysis.c14_signal.standardized, color=C14_COLOR, label=r"$\Phi_{14\mathrm{C}}$")
    _overlay_events(axes[0], analysis, pairs, event_color)
    axes[0].set_ylabel("Detrended potential (σ)")
    axes[0].set_title(title)
    axes[0].legend(ncol=2, loc="upper left")

    periods = xwt.periods
    period_selection = (periods >= 8.0) & (periods <= 40.0)
    power = xwt.power[period_selection][:, selection]
    normalized_log_power = np.log2(power / np.nanmedian(power))
    levels = np.linspace(-4, 4, 17)
    contour = axes[1].contourf(
        years[selection],
        periods[period_selection],
        np.clip(normalized_log_power, -4, 4),
        levels=levels,
        cmap="viridis",
        extend="both",
    )
    significance = xwt.significance_ratio[period_selection][:, selection]
    if np.nanmin(significance) <= 1.0 <= np.nanmax(significance):
        axes[1].contour(
            years[selection], periods[period_selection], significance, levels=[1.0], colors="white", linewidths=0.8
        )
    coi = np.clip(xwt.cone_of_influence[selection], periods[period_selection].min(), periods[period_selection].max())
    axes[1].plot(years[selection], coi, color="#202020", linestyle="--", linewidth=0.9, label="Cone of influence")
    axes[1].set_ylim(40, 8)
    axes[1].set_ylabel("Period (years)")
    axes[1].set_xlabel("Calendar year (AD positive, BC negative)")
    axes[1].set_xlim(left_year, right_year)
    colorbar = figure.colorbar(contour, ax=axes[1], pad=0.015)
    colorbar.set_label("Cross-wavelet power (log₂ relative scale)")
    figure.tight_layout()
    _save(figure, output_dir, stem)


def figure5(analysis: StudyAnalysis, output_dir: Path) -> None:
    _wavelet_focus(
        analysis,
        output_dir,
        stem="figure5_cross_wavelet_grand_minima",
        title="Cross-wavelet behavior during selected grand minima",
        pairs=analysis.classification.minima,
        event_color=MINIMUM_COLOR,
        left_year=-4000,
        right_year=-6000,
    )


def figure6(analysis: StudyAnalysis, output_dir: Path) -> None:
    _wavelet_focus(
        analysis,
        output_dir,
        stem="figure6_cross_wavelet_grand_maxima",
        title="Cross-wavelet behavior during selected grand maxima",
        pairs=analysis.classification.maxima,
        event_color=MAXIMUM_COLOR,
        left_year=600,
        right_year=-550,
    )


def _plot_state_segments(
    axis: plt.Axes,
    years: np.ndarray,
    values: np.ndarray,
    analysis: StudyAnalysis,
    nuclide: str,
) -> None:
    attribute = "be10" if nuclide == "10Be" else "c14"
    for pairs, color in [
        (analysis.classification.minima, MINIMUM_COLOR),
        (analysis.classification.maxima, MAXIMUM_COLOR),
    ]:
        for pair in pairs:
            event = getattr(pair, attribute)
            age = 1950.0 - years
            mask = (age >= event.onset_bp) & (age <= event.end_bp)
            axis.plot(years[mask], values[mask], color=color, linewidth=1.7)


def figure7(analysis: StudyAnalysis, output_dir: Path) -> None:
    figure, axes = plt.subplots(4, 1, figsize=(11.2, 8.4), sharex=True)
    panels = [
        (analysis.be10_signal, "high_pass", r"High-pass $\Phi_{10\mathrm{Be}}$", "10Be"),
        (analysis.be10_signal, "moving_variance", r"25-year variance, $\Phi_{10\mathrm{Be}}$", "10Be"),
        (analysis.c14_signal, "high_pass", r"High-pass $\Phi_{14\mathrm{C}}$", "14C"),
        (analysis.c14_signal, "moving_variance", r"25-year variance, $\Phi_{14\mathrm{C}}$", "14C"),
    ]
    for axis, (prepared, field, label, nuclide) in zip(axes, panels, strict=True):
        values = getattr(prepared, field)
        axis.plot(prepared.calendar_year, values, color="#303030", linewidth=0.85)
        _plot_state_segments(axis, prepared.calendar_year, values, analysis, nuclide)
        axis.set_ylabel(label)
    axes[-1].set_xlabel("Calendar year (AD positive, BC negative)")
    axes[-1].set_xlim(
        float(analysis.be10_signal.calendar_year[0]),
        float(analysis.be10_signal.calendar_year[-1]),
    )
    axes[0].set_title("Thirty-year high-pass signals and moving variance")
    figure.tight_layout()
    _save(figure, output_dir, "figure7_high_pass_and_moving_variance")


def _cdf_panel(axis: plt.Axes, samples: dict[str, np.ndarray], title: str, xlabel: str) -> None:
    for state, label, color in [
        ("minimum", "Grand minima", MINIMUM_COLOR),
        ("maximum", "Grand maxima", MAXIMUM_COLOR),
        ("moderate", "Moderate activity", MODERATE_COLOR),
    ]:
        x, y = empirical_cdf(samples[state])
        axis.step(x, y, where="post", color=color, label=label)
    axis.set_title(title)
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Cumulative probability")
    axis.set_ylim(0, 1.01)


def figure8(analysis: StudyAnalysis, output_dir: Path) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(10.8, 7.4))
    _cdf_panel(
        axes[0, 0], analysis.distribution_samples["10Be"]["variance"], r"Variance of HP $\Phi_{10\mathrm{Be}}$", "Variance"
    )
    _cdf_panel(
        axes[0, 1], analysis.distribution_samples["14C"]["variance"], r"Variance of HP $\Phi_{14\mathrm{C}}$", "Variance"
    )
    _cdf_panel(
        axes[1, 0], analysis.distribution_samples["10Be"]["power"], r"20–40 year power of HP $\Phi_{10\mathrm{Be}}$", "Mean wavelet power"
    )
    _cdf_panel(
        axes[1, 1], analysis.distribution_samples["14C"]["power"], r"20–40 year power of HP $\Phi_{14\mathrm{C}}$", "Mean wavelet power"
    )
    axes[0, 0].legend(fontsize=8)
    for axis in axes.ravel():
        axis.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    figure.tight_layout()
    _save(figure, output_dir, "figure8_activity_state_distributions")


def write_all_figures(analysis: StudyAnalysis, output_dir: Path) -> None:
    configure_style()
    for function in [figure1, figure2, figure3, figure4, figure5, figure6, figure7, figure8]:
        function(analysis, output_dir)
