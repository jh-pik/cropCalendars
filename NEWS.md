# cropCalendars NEWS

## 0.2.0 (branch fix-alg-vectorize-phu)

### Bug fixes

- **`calcMonthlyClimate`**: floor the monthly P/PET (`mppet`) denominator with
  `pmax(mpet_y, 1e-6)`, mirroring the daily `dppet`. FAO-56 PET can be clamped to 0 in deep
  cold; a month with zero PET (and, with zero precipitation, `0/0`) produced `Inf`/`NaN` in
  `mppet`, which propagated to high-latitude cells and crashed downstream rules
  (`calcHarvestDateVector`: `min(monthly_ppet)`). Now finite everywhere.

- **`calcPET_FAO56`**: preserve the shape of the inputs. The final `pmax(0, pet)` dropped the
  `dim` attribute, so passing matrices (cell-vectorised callers, e.g. all cells × days at once)
  returned a vector. Now clamps at 0 while restoring `dim`; NA-safe; vector inputs unchanged.

- **`date_to_doy(skip_feb29 = TRUE)`**: corrected the leap-year fold point. The condition
  `doy1 > 28` mis-folded at Jan. 29 (collapsing Jan. 29 onto DOY 28 and shifting the late-Jan/
  Feb day-of-year by one in leap years). Changed to `doy1 > 59` ("after Feb. 28", whose
  day-of-year is 59) so Feb. 29 folds onto DOY 59 (shared with Feb. 28) and Mar. 1–Dec. 31 stay
  aligned with non-leap years. This shifts the daily temperature / P-PET climatologies
  (`dtemp`/`dppet` from `calcMonthlyClimate`) in leap years and hence some sowing /
  threshold-crossing dates.

- **`calcDoyWetMonth`**: rewrote to operate on a 365-value daily P/PET climatology
  (DOY 1–365) instead of 12 monthly values. The previous version interpolated monthly
  values to daily and then folded the result onto a doubled (730-day) array, producing
  index-wrapping errors of up to 351 days in the wet-season onset detection, particularly
  in equatorial regions. The new implementation applies a 120-day circular rolling window
  directly to the daily climatology.

- **`calcDoyCrossThreshold`**: rewrote to operate on a 365-value daily temperature
  climatology instead of 12 monthly values. The previous monthly interpolation produced a
  doubled (730-day) array whose modular index reduction was incorrect for late-year
  crossings. The new implementation uses a direct crossing-detection with a Dec→Jan
  circular wrap via `c(x[365], x[1:364])`.

- **`calcPET_FAO56`**: removed unused `tmax` and `tmin` parameters; FAO-56 Penman-
  Monteith is formulated in terms of daily mean temperature only. The previous function
  signature accepted these parameters and the pipeline passed them, but they were never
  used. Aligned surface albedo to FAO-56 standard value (0.23). Changed
  `return(max(0, pet))` to `return(pmax(0, pet))` to support vectorized calls.

- **`calcPET`**: changed `eeq <- max(0, ...)` to `pmax(0, ...)` in the observed-Rn
  branch to support vectorized calls over all days at once.

- **`generatePHUTserie_isimip3`**: fixed four namespace errors introduced when the
  package API was renamed: `wintercrop` → `isWinterCrop`, `calc.vd` → `calcVd`,
  `calc.vrf` → `calcVrf`, `grid$lat` → `grid_df$lat`.

### New features

- **`calcPET` — observed-radiation branch**: when `swdown` (rsds) and `lwdown` (rlds)
  are supplied, net radiation is computed from actual flux observations rather than
  estimated from orbital geometry. This removes the dependence on a fixed sunshine
  fraction and improves PET accuracy over regions with variable cloud cover.

- **`calcPET_FAO56`**: new FAO-56 Penman-Monteith reference ET function. Requires
  daily mean temperature (`tas`), wind speed at 10 m (`sfcwind`), specific humidity
  (`huss`), downwelling shortwave radiation (`rsds`), downwelling longwave radiation
  (`rlds`), and surface pressure (`ps`, default 101325 Pa). Aerodynamic resistance
  follows LPJmL `getpet.c` (grass reference crop, `rs = 70 s/m`).

- **`calcMonthlyClimate`**: added `pet_method` argument (`"pt"` or `"fao56"`). When
  `pet_method = "fao56"`, passes `sfcwind`, `huss`, `swdown`, `lwdown`, `ps` to
  `calcPET_FAO56` in a single vectorized call. Returns two new list elements `dtemp`
  and `dppet`: 365-value daily climatologies (mean temperature and P/PET ratio by DOY,
  averaged across years) consumed by the revised seasonality functions.

- **`calcSowingDate`**: updated signature to accept `daily_ppet` (365-value vector,
  replaces `monthly_ppet`) and `daily_temp` (365-value vector, replaces `monthly_temp`
  for threshold-crossing detection). Monthly temperature is still used internally for
  coldest-month calculations.

- **`01_calc_crop_calendars.R`** (pipeline script): extended to read five additional
  ISIMIP3b climate variables (`rsds`, `rlds`, `huss`, `sfcwind`, `ps`) and pass them
  to `calcMonthlyClimate(pet_method = "fao56")`. Memory use optimised: replaced
  sequential `abind` calls (O(n) copies) with list accumulation and a single
  `do.call(abind, ...)`, pre-allocated output list to avoid O(n²) `rbind` growth,
  reduced progress printing from every pixel to every 500th.

### Performance

- **`generatePHUTserie_isimip3`**: replaced the `for (i in 1:NCELLS)` cell-by-cell
  loop with fully vectorized operations across all land cells simultaneously:
  - `lin_idx` flat-index mapping (`(ilat−1)×720 + ilon`) pre-computed once before the
    time-slice loop — eliminates 67420 `which(lons == ...)` calls per time slice.
  - Climate reading changed from a single multi-year read to a year-by-year accumulation
    (`tas_sum / nyears`) that caps peak RAM regardless of time-slice width.
  - Sdate/hdate extraction via `matrix(sdate, nrow=720*360)[lin_idx, ]` + `rowMeans`
    replaces per-cell `ncvar_get` slices.
  - Seven unexported helper functions added: `.monthly_temps_vec`, `.calc_vd_vec`,
    `.wintercrop_vec`, `.phu_cumsum`, `.calc_phu_thermal_vec`,
    `.calc_phu_vernal_vec`, `.build_vrf_mat`. All verified numerically equivalent to
    their scalar counterparts (`calcPHU`, `calcVd`, `isWinterCrop`).

---

## 0.1.2 and earlier

See git log on `master` branch.
