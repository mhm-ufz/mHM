# Test Fixtures

This directory contains non-runnable fixture inputs used by `test_nml` smoke and scenario cases.
It is separate from `examples/domain_01` and `examples/domain_02`, which are runnable v6 example setups.

## Fixture groups

- Synthetic fixtures used only by `test_nml`:
  `fdir_simple.asc`, `scc_gauges_xy.txt`, `runoff.nc`, `lc_static.nc`,
  `lai_monthly_cycle.nc`, `lai_monthly_series.nc`, `lai_yearly.nc`,
  `pre_6h.nc`, `pet_6h.nc`, `tavg_6h.nc`, `tmin_6h.nc`, `tmax_6h.nc`,
  `pre_daily_aligned.nc`, `pet_daily_aligned.nc`, `tavg_daily_aligned.nc`,
  `pre_weights.nc`, `pet_weights.nc`, `tavg_weights.nc`, `scc_gauges.nc`
- Copied legacy-domain assets:
  `dem.asc`, `fdir.asc`, `slope.asc`
- Derived/generated NetCDF fixtures:
  `input.nc`, `lc_periods.nc`

These files intentionally stay flat in this pass so existing `test_nml` path usage remains simple.
