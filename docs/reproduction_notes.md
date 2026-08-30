# Reproduction method and known differences

## Reproduced outputs

The pipeline generates Figures 1–8, Tables 1–5, and a numerical validation
report. It uses the two committed aligned SMP series and the committed
160-event catalogue; no uncommitted input is required.

The implementation reproduces:

- all 32 Table 2 grand minima and all 21 Table 3 grand maxima in order,
  including centers and durations for both nuclides;
- the four printed two-Gaussian amplitude thresholds;
- the Table 4 waiting-time fits and the Table 5 distribution fits, subject to
  the numerical differences below;
- summed event durations of 2254/2312 years for 10Be/14C minima and 1307/1347
  years for maxima; and
- the paper's 99% Kolmogorov–Smirnov conclusions for activity-state variance
  and 20–40 year wavelet power.

## Declared calculation

1. Fit a two-component univariate Gaussian mixture to all 160 amplitudes for
   each nuclide.
2. Use the two component means as the low- and high-activity thresholds. The
   calculated values are `-0.92478` and `+1.34629` for 10Be and `-0.66866` and
   `+1.41096` for 14C.
3. Require both members of a paired event to cross their nuclide-specific
   threshold in the same direction.
4. Require both event durations to be at least 22 years. The inclusive
   boundary is necessary because the published lists contain a 22-year minimum
   and a 22-year maximum, despite prose that says "longer than" 22 years.
5. Convert event centers from BP to calendar year with `round(1950 - center_BP)`.

The signal analysis fits and subtracts a degree-5 polynomial, then uses
population standardization. AIC and BIC use their conventional definitions,
`2k - 2 log L` and `k log n - 2 log L`; the labels of Equations (1a) and (1b)
in the paper are reversed. A single Gaussian has two fitted parameters, and a
two-component one-dimensional Gaussian mixture has five.

Figure 7 uses a zero-phase, tenth-order Butterworth high-pass filter on annual
data with a 30-year cutoff (`Wn = 2/30` in Nyquist-normalized form), followed by
a centered 25-year moving variance. Figures 5 and 6 use the cross-wavelet
implementation in `src/grand_solar_reproduction/wavelet.py`.

## Catalogue and figure conventions

The catalogue preserves the row pairing used by the study. Ten moderate or
weak pairs do not have positive literal interval overlap (one pair only
touches), although all 53 published grand-event pairs overlap. The pipeline
therefore does not reinterpret or reorder the released catalogue.

Figure 2 uses a correct Freedman–Diaconis histogram containing all 160 events.
The published histogram bars contain 159 observations per proxy, while its
mixture and information-criterion calculations use all 160. The statistical
fit is reproduced from all observations; the bar omission is not repeated.

## Known numerical differences

- **Table 1, 14C Gaussian AIC:** the conventional likelihood calculation gives
  `538.73`, which rounds to `539`; the paper prints `540`. Every other Table 1
  entry rounds as printed.
- **Table 5, grand-maximum two-Gaussian fit:** the printed `380/389` corresponds
  to a collapsed local mixture solution. The global maximum gives `378/387`.
  The pipeline calculates the collapsed solution deterministically for this
  published comparison.
- **Table 5, grand-maximum lognormal BIC:** the calculated value is `384.55`,
  which rounds to `385`; the paper prints `384`. Its AIC and all grand-minimum
  entries agree.
- **Table 4, grand-maximum power-law lower bound:** the implementation declares
  `xmin = 180` years, with 34 tail observations and the finite-sample alpha
  correction, matching the printed alpha/AIC/BIC. An unconstrained continuous
  KS search selects 147 years.
- **Figure 4 duration bars:** the published lower-panel counts differ from
  histograms of the final Table 2/3 durations by one or two observations in
  several bins. The reproduction calculates the bars from the released final
  catalogue rather than copying the plotted heights.

## Limits of exact reproduction

- The committed 10Be SMP is already chronology-aligned to the 14C reference.
  The original running-window length, lag range, maxima-selection rule, and
  endpoint handling are not specified in the released method, so the repository
  does not repeat that upstream alignment.
- The stochastic seed and synthetic-sample size for the waiting-time
  two-sample KS checks are not specified.
- The paper does not specify a test statistic, simulation rule, or sample count
  for the "lumping" Monte Carlo test, so that result cannot be recalculated
  uniquely and is not replaced by a different test.
- The paper specifies the 30-year high-pass cutoff but not the filter order.
  This implementation explicitly declares order 10, annual sampling, and
  forward/backward filtering.
- The standard-deviation denominator, interpolation endpoint convention, and
  filter endpoint padding are not fully specified. This implementation declares
  population standardization, the common observed annual grid, and SciPy's
  forward/backward SOS filtering.
- The cross-wavelet calculation follows the Grinsted/Torrence–Compo equations;
  numerical and rendering differences from the original MATLAB figures are
  possible and bitwise equality is not claimed.

These boundaries are intentional: the implementation exposes each calculation
it performs and does not invent unspecified upstream operations.
