# cropCalendars NEWS

## 0.2.0 — annual sliding-window pipeline

### Unified daily-climatology smoothing window
- **One windowing primitive, `.circRoll`.** The forward-sum `.circRollSum` and the
  centred-mean `.smoothCycle` are merged into a single `.circRoll(x, w, position, mean)`
  that re-indexes the one rolling sum by a chosen anchor day: `position = "start"` (default,
  the old forward sum — the 120-day wettest window relies on its start indexing),
  `"center"` (circular shift by `w %/% 2`), or `"end"` (window's last day); `mean = TRUE`
  divides by the window length. The manual half-width re-centring in `.doyWarmestWindow`
  (now a plain centred argmax) is gone, and `.dailyPpetDiff` is rewritten as its natural
  form — the window **ending** at the junction minus the window **starting** at it (the
  trailing-vs-leading moisture trend, exactly the monthly `mppet[m] − mppet[m+1]`). At odd
  `w` (the 31 default) this is identical to the previous centred formulation. `.smoothCycle`
  and `.circRollSum` are removed.
- **`smooth_window` default 30 → 31 (odd).** `.circRoll` no longer force-odds the window, so
  the package now has ONE centring convention (`h = w %/% 2`). An odd window is exactly
  symmetric for the centred smoothing and argmax; an even window centres half a day low
  (documented caveat on `.circRoll`). 31 keeps the ~1-month width while making the crossings,
  `.doyWarmestWindow`, and `.dailyPpetDiff` all exactly symmetric. (Previously
  `smooth_window = 30` was silently force-odded to 29 for the crossings only — this removes
  that hidden asymmetry; it changes the effective crossing window 29 → 31, so re-validate.)
- **Removed the dormant build-time smoother.** `finalizeClimate()` (and `ringClimatology()`,
  `calcClimatology()`) had a `smooth_window` argument that applied a centred circular running
  mean (`.circSmooth`) to the daily climatology *at build time*. The pipeline now smooths at
  rule time only (`smooth_window` → `.smoothCycle`), 01a writes the RAW climatology, and no
  caller passed the build-time knob (it defaulted to off everywhere). The argument, the
  `.circSmooth` helper, and the duplicate centred-mean code path are dropped; the rule-time
  `smooth_window` is now the single, unambiguous smoothing knob.
- **`cross_smooth_window` → `smooth_window`: one global window for BOTH crossings and
  reductions** (`calcCropCalendars` and all callees; default **31**, see above). Previously the
  threshold-crossing inputs used a tunable 15-day `.smoothCycle` while the extremum/driest
  reductions used a structural 30-day window — two windows on the same daily climatology.
  They are unified: `smooth_window` now drives the crossing smoothing AND the
  `.warmestWindowMean`/`.coldestWindowMean`/`.doy*Window`/`.dailyPpetDiff` reductions (and the
  always-wet `min(daily_ppet)` floor; see below). At 30 the reductions are unchanged (already
  30); the crossings move 15→30.
- **Always-wet floor and wet-end crossing share one P/PET series.** The always-wet aridity
  test now reads `min_ppet <- min(daily_ppet)` — the minimum of the very curve the wet-end
  existence/level crossing (`doy_wet1`) tests against `ppet_ratio`. Under the unified
  `smooth_window` this is bit-identical to the old separate `.driestWindowPpet` (the per-cycle
  minimum is anchor-invariant and ratio-of-sums == ratio-of-means), so it is an exact no-op for
  the default path but makes `ppet_min` and `ppet_ratio` compare against the SAME curve, closing
  the former `ppet_min`/`ppet_ratio` dead-zone for free. The `.driestWindowPpet` helper and the
  now-vacuous `cc_harmonize_minppet` toggle are removed. (The toggle's earlier +40% regression
  was a 15-vs-30-day window mismatch, gone once the two windows were unified.)
- **`.dailyPpetDiff` normalised to a per-30-day rate** (`× 30/lag`): the `doy_wet2`
  moisture-trend threshold `ppet_ratio_diff` is calibrated as Δ(P/PET) per 30 days, so the
  diff is rescaled to that horizon and stays valid for any `smooth_window` (no-op at 30).
  Only this reduction needed normalising — the others are means/ratios/argmin compared to
  absolute thresholds, whose units don't scale with the window.

### Deduplicate scalar vs. vectorized rule helpers (PR-B)
- **`.monthlyFromDaily` now subsumes `.monthly_temps_vec`.** The two did the same calendar-month
  aggregation over identical (non-leap) month boundaries — one for a single 365-day vector, the
  other row-wise over an `[ncells x 365]` matrix. `.monthlyFromDaily(daily, agg)` now accepts
  either (vector in → length-12 vector out; matrix in → `[ncells x 12]` out) via `rowMeans`/
  `rowSums`, and `.monthly_temps_vec` is removed; the gridded PHU path calls
  `.monthlyFromDaily(tas_mean_day, "mean")`. Vector results are bit-identical; the matrix path
  reproduces the old helper exactly (verified, diff 0).

### Deduplicate scalar vs. vectorized rule helpers (PR-C)
- **`calcVrf` and `calcPHU` are now thin wrappers over their vectorized cores**
  (`.build_vrf_mat`, `.calc_phu_thermal_vec` / `.calc_phu_vernal_vec`), completing the
  scalar/`_vec` consolidation. The PHU pipeline already ran only the vectorized forms.
  `calcVrf` is equivalent to the old scalar to ~1e-14 (the cumsum core reorders the
  summation) and now returns a plain vector instead of a 1-D array. `calcPHU` matches the
  old scalar everywhere except two edge bugs the vectorized core (and therefore the pipeline)
  already fixed: at `sdate == 1` the old code dropped day 1 from the growing period, and at
  `sdate == hdate` it dropped ~2 days; both are handled correctly now. `calcPHU` also errors
  (rather than printing) on an unsupported `phen_model`.
- **Equivalence harness extended** (`tests/testthat/test-vec-equivalence.R`): adds the
  original `calcPHU`/`calcVrf` bodies as ground truth and fuzzes ~12 000 more inputs —
  `calcPHU` asserted bit-identical off the two edges (and asserted to *correct* them on the
  edges), `calcVrf` to numeric tolerance across all `sdate`/`hdate` including the edges.

### Deduplicate scalar vs. vectorized rule helpers (PR-A)
- **`calcVd` and `isWinterCrop` are now thin wrappers over their vectorized cores**
  (`.calc_vd_vec`, `.wintercrop_vec`). The PHU pipeline already ran only the vectorized
  forms; the public scalars carried a second, drift-prone copy of the same vernalization-day
  and winter-crop rules. Each scalar now calls its `_vec` core with a single row, so there is
  one source of truth per rule. Public signatures are unchanged. `isWinterCrop` now returns NA
  (not an error) when `end` is NA, and `calcVd` keeps exactly `max.vern.months` coldest months
  (the scalar always took 5 regardless) — identical at the production default of 5.
- **Equivalence harness** (`tests/testthat/test-vec-equivalence.R`): embeds the original
  scalar bodies as ground truth and fuzzes ~9 000 inputs (plus rule-boundary and NA-guard
  cases), asserting the wrappers reproduce them bit-for-bit. This is also the scaffold for
  PR-C (`calcVrf`/`calcPHU`, whose control flow differs from their `_vec` cores).
- **Vectorized helpers documented.** Each `_vec`/cumsum helper in
  `generatePHUTserie_isimip3.R` gained a header comment explaining its logic against the
  readable scalar form (the matrix `order`/`cumsum`/`ifelse` tricks that replace the loops).

### Daily-climatology rule port + ring slimming
- **All extremum reductions moved from the 12 calendar months to daily 30-day windows.**
  The coldest/warmest-month temperature (`calcSeasonality` min-temp, `calcSowingDate`
  winter-type tests, `calcHarvestRule`/`calcHarvestDate` warmest-month guards), the
  warmest-month mid-day, the coldest-day sowing anchor, the driest-month P/PET
  (`hd_wetseas` fallback) and the P/PET month-over-month difference (`doy_wet2`) now use
  30-day windows / centres of the daily climatology (`R/zz_daily_reductions.R`). A daily
  mean is already a per-DOY 30-year mean, so the 30-day window reproduces the
  calendar-month value continuously, removing the ~30-day flips when the warmest/coldest
  (or wettest/driest) month alternated between adjacent windows. The warm-winter sowing
  date and the spring up-crossing are now anchored at the **coldest day** (centre of the
  coldest 30-day window). The seasonality **CV classifier stays on 12 monthly bins**
  (now reconstructed from the daily cycle) — the calibrated 0.4/0.010 thresholds are
  hard cutoffs, so a daily/overlapping CV would reclassify grazing cells.
- **Crossing-detector smoothing relocated out of the climatology** and renamed
  `clm_smooth_window` → `cross_smooth_window`. The smoothing (centred running mean,
  reject 1-day jitter) is now applied to **each crossing detector's own input**
  (`.smoothCycle`, threaded through `calcSowingDate`/`calcHarvestDateVector`), not to the
  stored climatology — which stays raw, so the extremum reductions and the 120-day
  wettest window are unaffected. Numerically equivalent (the two smoothers are
  floating-point identical; crossings unchanged). In the split pipeline this makes
  `cross_smooth_window` tunable in `01b` without recomputing the `01a` cache.
- **Bug-free monthly fallbacks.** With no daily series a direct caller now falls back to
  the exact legacy monthly rule (`min/max(monthly_*)`, coldest/warmest-month mid-day),
  and the wettest-window uses a 4-month `Σ₄P/Σ₄PET` ratio-of-sums on the monthly totals
  (`.wetDoyMonthly`) — summing P and PET **separately**, never the monthly P/PET *ratio*
  (the legacy `Σ`-of-ratios form let one near-zero-PET month explode and dominate the
  window — the bug that motivated the daily redesign). Threshold crossings interpolate
  monthly→daily in fallback (no pure-monthly crossing exists), matching the reference.
- **Sliding ring carries only daily data.** The streaming engine gained a `daily_only`
  mode: the ring allocates/accumulates/serves only `dtemp`/`dprec`/`dpet` (no monthly
  sums, no `dppet`), trimming the per-slot footprint (~4% smaller ring → less fork
  pressure). `calcCropCalendars` reconstructs the 12 monthly temp/precip values for the
  seasonality CVs via `.monthlyFromDaily` (mean / sum over each month's DOYs — exact
  except a sub-percent February residual from the Feb-29→DOY-59 fold). `01a` now caches
  only the daily climatology. Verified bit-identical on the daily fields; 42/42 cell×crop
  cases identical in seasonality, sowing and maturity. (An existing `01a` cache must be
  regenerated for the new `01b`.)
- **Climate engine renamed** to drop the now-misleading "Monthly":
  `initMonthlyClimate`→`initClimateAccum`, `addYearMonthlyClimate`→`addYearClimate`,
  `finalizeMonthlyClimate`→`finalizeClimate`, `calcMonthlyClimate`→`calcClimatology`
  (files/man/tests to match). After the daily port these primarily build the *daily*
  climatology. Pure rename, no behaviour change.

### Wet-window hysteresis — anchor fix + Gaussian distance kernel
- **Fixed the wettest-window hysteresis anchor bug.** The per-cell state carried between
  years (`prev_wet`) was taken from crop #1's *sowing day* (`calcCropCalendars` set
  `attr(., "wet_doy") <- sowing_day` for PREC/PRECTEMP cells). That equals the
  wettest-window argmax only for crops that take the wet-season sowing branch — but a
  **vernalizing/winter crop (e.g. `Winter_Wheat`, crop #1 in the combined 7-crop run)
  takes the winter branch even in PREC/PRECTEMP cells**, so its sowing day is NOT the
  wet-window DOY. The wrong anchor was then fed back as `prev_doy` for *every* crop,
  flipping the eps-weighted pick at the cold-start → first-sticky-year boundary: on
  GFDL-ESM4 historical Maize this produced a spurious 1850→1851 jump in 21,015 cells
  (then frozen 1851–1880, since the bogus anchor was constant). Maize-only runs never
  exposed it (there crop #1 *is* Maize). Fix: `wet_doy` is now computed
  **crop-independently** from the climate via `calcDoyWetMonth(dprec, dpet, …)`; after the
  fix the 1850→1851 jump is **0 cells**. (Any combined-run output produced before this fix
  is contaminated and must be regenerated.)
- **`calcDoyWetMonth` distance weight is now Gaussian** with an explicit floor, replacing
  the linear ramp: `w(x) = (1−eps) + eps·exp(−(x/decay)²)`, `x = Δ/(365/2)`. Same anchors
  as before — `w(0)=1` (small drifts essentially free) and `w(antipode)=1−eps` (the floor;
  a far peak must be >1/(1−eps)× wetter to win) — but a **faster mid-distance decline**
  (new `decay` param; config `wet_window_decay=0.3`) that closes the soft 1.5–3.5-month
  band where bimodal "chronic flippers" lived: a rival ~3 months out must now be ~1.9×
  wetter vs ~1.3× under the linear kernel. `eps` keeps its meaning (floor depth);
  `eps=0` / `prev_doy=NA` still reproduces the plain argmax. Threaded as `wet_window_decay`
  through `calcSowingDate`/`calcCropCalendars` and the `01b`/`01` drivers.

### Major
- **Pipeline rebuilt from a 10-year-step scheme to an annual 30-yr sliding window**
  (`R/slidingClimate.R` ring buffer; calendars computed every year, smooth and
  rule-consistent, so the output moving average is dropped). Driver
  `utils/ggcmi_ph3/01_compute_annual_calendars.R` (single-pass) — or the split below.
  Future scenarios seed the ring from historical climate; opt-in seed cache (`SAVE_SEED`).
- **Stage 01 split into 01a (climatology) + 01b (calendars)** for RAM safety and fast
  rule iteration (`utils/ggcmi_ph3/01a_climatology_annual.R`, `01b_calendars_annual.R`).
  `01a` streams the raw climate through the ring and writes the per-year raw daily
  climatology to disk (fork-free; bounded ~25 GB). `01b` loads one year's climatology at
  a time and runs `calcCropCalendars` via `mclapply` — so the heavy ring buffer is NOT
  live in the parent during the 64-worker fork, which removes the copy-on-write OOM
  (the single-pass driver held the ~24 GB ring+output across the fork → 340 GB at 64
  cores; see `R_GC_MEM_GROW=0` note in the `.sh`). 01b re-runs in minutes on any
  sowing/harvest rule change without re-reading the raw climate. Both take a `YEARS=`
  env subset (e.g. `"2000:2014"`) for dev iteration. Cost: ~0.6 GB/year climatology
  cache (~100 GB for a 165-yr historical run), regenerable.
- **Stage 02** (`utils/ggcmi_ph3/02_assemble_annual_ncdf.R`) writes the
  publication-ready ISIMIP3b DRS NetCDF in one pass (final variable names,
  "years since 1601" time axis, **ascending latitude**, `_FillValue`/`missing_value`,
  per-timestep chunking, publish path) — eliminating the NCO/CDO stages `04`–`07`
  (and a `cdo -L` segfault). Fixes a **latitude inversion** vs the official product.

### Algorithm
- `calcDoyWetMonth` + `hd_wetseas` use **ΣP/ΣPET** (ratio of summed P and PET) from new
  daily `dprec`/`dpet` climatologies, fixing P/PET blow-ups where PET≈0 (spurious
  ~80-day sowing flips in monsoon cells).
- `calcDoyWetMonth` gains **distance-weighted, max-normalised selection** (`prev_doy`,
  `eps`) for the sliding window: ~57 % of PREC/PRECTEMP cells have a second 120-day
  ΣP/ΣPET peak within 10 % of the best, so the plain argmax flips between far-apart peaks
  between adjacent years (sowing, and the sowing-anchored harvest, oscillate). The window
  is now chosen as `argmax( (ws/max ws) · max(1 − eps·Δ/(365/2), 0) )`, where `Δ` is the
  circular DOY distance to last year's window. Near peaks are essentially free (the same
  peak drifts), a far peak is penalised but — since the weight floors at `1−eps`, never 0 —
  still wins when decisively better (genuine regime shift). Normalising by the year's best
  makes `eps` dimensionless; `wet_window_eps = 0.5` (config) cut the wet-cell mean
  year-to-year jump 1.59 → 0.13 d and cells-ever-flipping 0.257 → 0.077 on GFDL-ESM4.
  `eps = 0` (or no `prev_doy`) reproduces the plain argmax — backward compatible.
  (`drift_gate` is deprecated/ignored; the smooth weight subsumes the old small-move gate.)
- `calcSeasonality` gains **threshold-deadband hysteresis** (`prev_seas`, `seas_eps`,
  `mtemp_margin`): ~18 % of cells flip their seasonality *class* between years by grazing a
  classifier threshold (`CV_prec` vs 0.4, `CV_temp` vs 0.010, `min_temp` vs 10 °C), which
  swaps the entire sowing rule. Each threshold is now relaxed toward keeping last year's
  class (thermostat deadband: a test the previous class was on the high side of uses
  `thr·(1−seas_eps)`, else `thr·(1+seas_eps)`; the `min_temp` test uses an absolute
  `mtemp_margin` °C). `seas_eps = 0.25` (config) cut year-to-year class flips 18.3 % →
  2.8 %. Class is crop-independent, so the resolved class is returned as
  `attr(., "seas_type")` and carried per-cell by the driver. `seas_eps = 0` reproduces the
  plain Waha thresholds — backward compatible.
- **Daily-climatology smoothing + sustained-crossing guard** — the dominant source of
  year-to-year sowing/harvest oscillation was traced to the *temperature* branch, not the
  wet-window: the per-DOY daily means (`dtemp`/`dprec`/`dpet`) carry ~1 °C / spiky
  day-to-day jitter, and the point detectors (`calcDoyCrossThreshold` for the spring/fall
  temperature crossings, and the wet-season-end P/PET crossing in `calcHarvestDateVector`)
  latch onto single-day blips. A 1-day dip-and-recover through `temp_spring` in the autumn
  descent produces a spurious "first" up-crossing ~130 days before the real spring
  crossing, and sub-1 °C differences between 30-yr windows flip which side wins (verified
  on GFDL-ESM4 cell 136.25/−33.25, S. Australia). Two fixes: (a) a centred **circular
  running mean** applied to each crossing detector's input (config
  `cross_smooth_window`, default 15 d, `.smoothCycle`; see the rule-port entry above for
  the later relocation off the stored climatology); (b) `calcDoyCrossThreshold(min_duration=)` accepts a crossing only if
  the excursion **persists** that many days (config `cross_min_duration`, default 5;
  threaded through `calcSowingDate`/`calcHarvestDateVector`/`calcCropCalendars`). Smoothing
  is the decisive lever (collapses the bistable flip to the genuine spring crossing); the
  duration guard is a secondary net. The 120-day wet-window argmax is essentially unchanged
  (already integrates 120 d). `smooth_window ≤ 1` and `min_duration = 1` reproduce the
  prior behaviour exactly — backward compatible. The wet-window hysteresis is now **off by
  default** (`wet_window_eps = 0`); it remains available for the genuine PREC-branch
  bimodal near-tie, which smoothing does not address.
- `calcHarvestDateVector` harvest dates evaluated on the **daily** climatology:
  `hd_temp_base` = centre of the warmest 30-day window (was warmest-month mid-day);
  `hd_wetseas`/`hd_temp_opt` fixed (had computed month indices as DOYs).
- `calcCropCalendars` ~2.5× faster (cumulative-sum rolling windows; `list2env`
  unpacking; optional pre-extracted `crop_parameters`), output bit-identical.
- `generatePHUTserie_isimip3`: reads `planting_day`/`maturity_day` from the DRS file;
  new `smooth_window` (default 1 = PHU matches the growing period exactly).

### Infrastructure
- Climate files **discovered dynamically** across the official ISIMIP roots (3b
  primary+secondary, 3a obsclim/spinclim); ensemble member / scenario / year-range read
  from file names; product ranges derived from the data. Uniform 67420-cell GGCMI mask
  (`ggcmi_landcells.csv`). Tunables centralised in `00_config.R`;
  `enms`/`syears`/`eyears`/`ccal_years` removed.

### Earlier on this branch (vectorisation, engine, I/O — pre-annual-refactor)

#### Changes

- **`generateCropCalTSerie_isimip3` (performance, ~16× on a GFDL-ESM4 test, 9.3 min → 35 s)**:
  vectorise the array-fill loop (one flat-index assignment instead of a per-pixel × per-year
  double loop with `which()` lookups); build `count.dt` once instead of growing it with `rbind`
  inside the pixel loop (was O(n²)); and parallelise the per-pixel smoothing via
  `parallel::mclapply` — new `ncores` argument (default 1 = serial `lapply`, so the package
  stays portable; the pipeline passes `SLURM_CPUS_PER_TASK`). Output is bit-identical
  (all 8 NetCDF variables, 0 diff vs the serial result).

- **Cell-vectorised, streaming monthly-climate engine** (`initClimateAccum` /
  `addYearClimate` / `finalizeClimate`, new `R/climateAccum.R`).
  Years are added one at a time, vectorised over an arbitrary number of grid cells, so the
  gridded pipeline can process the full grid one year at a time without the full multi-year
  series in memory. `calcClimatology` is now a thin per-pixel wrapper over this engine, so
  the aggregation, PET-method switch, P/PET flooring and leap-year DOY handling live in one
  place (previously duplicated in the pipeline). Output is unchanged (verified bit-identical
  to the previous gridded results); monthly-vector names retained.

- **`calcPET`**: `lat` and `day` are now optional — they are only needed for the orbital
  net-radiation estimate. With `swdown`/`lwdown` supplied (observed Rn), PET no longer
  depends on latitude/day-of-year.

#### Bug fixes

- **`generateCropCalTSerie_isimip3` — `seasonality` / `harv-reason` encoding**: these two
  ncdf variables were encoded with `as.numeric(as.factor(...))` applied *per pixel*, so the
  integer code was alphabetical among only the categories present at that pixel — every
  temporally-constant pixel collapsed to code 1 (≈94 % spurious `NO_SEASONALITY`;
  `harv-reason` → `GPmin`). Now mapped through fixed *global* factor levels
  (`seasonality` 1=NoSeas…5=TempPrec; `harv-reason` 1=GPmin…6=Thigh), matching the ncdf
  `long_name` and the standalone pipeline. Only these two variables change; sowing/harvest
  dates are unaffected (jump detection depends on change-points, not the absolute code).

- **`generateCropCalTSerie_isimip3`**: take `FYs`, `LYs`, `ggdir`, `ncdir`, `csvdir`, `pldir`
  as explicit arguments instead of reading them from the caller's global environment, and
  honour the `years_nc` argument — a stray `years_nc <- YEARSnc` had silently overwritten it
  with a global, making the argument dead. No behaviour change for the existing pipeline call.

- **`plotMap_ggplot` / `plotMapCropCalendars`**: qualify all external functions
  (`ggplot2::`, `scales::squish`, `RColorBrewer::brewer.pal`) so plotting works without those
  packages being attached; fix the `fil = landFill` typo (was a silently-ignored fill colour),
  and replace deprecated `aes_string`/`size` with `aes(.data[[…]])`/`linewidth`. `scales` and
  `RColorBrewer` added to Suggests.

- **`calcClimatology`**: floor the monthly P/PET (`mppet`) denominator with
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
  (`dtemp`/`dppet` from `calcClimatology`) in leap years and hence some sowing /
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

#### New features

- **`calcPET` — observed-radiation branch**: when `swdown` (rsds) and `lwdown` (rlds)
  are supplied, net radiation is computed from actual flux observations rather than
  estimated from orbital geometry. This removes the dependence on a fixed sunshine
  fraction and improves PET accuracy over regions with variable cloud cover.

- **`calcPET_FAO56`**: new FAO-56 Penman-Monteith reference ET function. Requires
  daily mean temperature (`tas`), wind speed at 10 m (`sfcwind`), specific humidity
  (`huss`), downwelling shortwave radiation (`rsds`), downwelling longwave radiation
  (`rlds`), and surface pressure (`ps`, default 101325 Pa). Aerodynamic resistance
  follows LPJmL `getpet.c` (grass reference crop, `rs = 70 s/m`).

- **`calcClimatology`**: added `pet_method` argument (`"pt"` or `"fao56"`). When
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
  to `calcClimatology(pet_method = "fao56")`. Memory use optimised: replaced
  sequential `abind` calls (O(n) copies) with list accumulation and a single
  `do.call(abind, ...)`, pre-allocated output list to avoid O(n²) `rbind` growth,
  reduced progress printing from every pixel to every 500th.

#### Performance

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
