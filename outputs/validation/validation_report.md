# Numerical validation

Overall status: **passed**.

## Exact checks

- [x] manual_catalog_has_160_pairs
- [x] manual_catalog_has_75_dips
- [x] manual_catalog_has_85_peaks
- [x] table2_all_32_rows_exact
- [x] table3_all_21_rows_exact
- [x] thresholds_round_to_published
- [x] minimum_duration_is_inclusive_22
- [x] table4_rounds_exactly
- [x] time_fraction_numerators_exact
- [x] all_published_grand_pairs_overlap
- [x] all_ks_comparisons_significant_at_99_percent
- [x] event_catalog_checksum
- [x] aligned_smp_checksums
- [x] aligned_smp_shared_published_grid

## Documented differences

- Table 1, 14C one-Gaussian AIC: calculated 539, published 540. The conventional likelihood calculation gives 539; the paper prints 540.
- Table 5, grand-maximum lognormal BIC: calculated 385, published 384. The input durations give 384.55, which rounds to 385; the paper prints 384.
