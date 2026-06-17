#' @title Calculate vector of possible harvest dates (Minoli et al., 2019)
#'
#' @description Harvest reasons are the events that can trigger the harvest
#' of a crop within an agro-climatic zone:
#' Earliest-maturing cultivar (GPmin);
#' Cultivar with longest grain filling (GPmaxrp);
#' Latest-maturing cultivar (GPmax);
#' Escape terminal water stress (w. lim);
#' Grain filling in warmest period (mid. t.);
#' Escape high temperature (high t.).
#'
#' The temperature- and water-driven harvest dates (\code{hd_wetseas},
#' \code{hd_temp_base}, \code{hd_temp_opt}) are evaluated on the \emph{daily}
#' climatology, consistent with \code{calcSowingDate}. Evaluating them on the
#' 12 monthly values instead quantises these dates to whole months, which makes
#' the resulting harvest dates flip by ~30 days between climate windows (e.g.
#' when the warmest month alternates between July and August).
#'
#' @param croppar data.frame with crop parematers as returned by getCropParam
#' @param sowing_date numeric value as day of the year (DOY). This can be either
#' caculated with calcSowingDate or prescribed.
#' @param sowing_season character value. Can be either "winter" or "spring". See
#' calcSowingDate.
#' @param monthly_temp numeric vector of length 12. Mean monthly air temperature
#' (degree Celsius).
#' @param monthly_ppet numeric vestor of length 12. Mean Potential Evapotranspiration (mm). See calcClimatology.
#' @param monthly_ppet_diff numeric vestor of length 12. Mean difference of
#' Potential Evapotranspiration (mm). See calcClimatology.
#' @param daily_temp numeric vector of length 365. Climatological daily mean
#' temperature (deg C), one value per DOY (the \code{dtemp} element from
#' \code{calcClimatology}). If \code{NULL}, it is interpolated from
#' \code{monthly_temp}.
#' @param daily_prec numeric vector of length 365. Climatological daily
#' precipitation (mm) per DOY (the \code{dprec} element from
#' \code{calcClimatology}). Used with \code{daily_pet} to form the spike-free
#' daily P/PET (ratio of per-DOY means) for the wet-season-end crossing. If
#' \code{NULL}, the daily P/PET is interpolated from \code{monthly_ppet} instead.
#' @param daily_pet numeric vector of length 365. Climatological daily PET (mm)
#' per DOY (the \code{dpet} element from \code{calcClimatology}).
#' @param cross_min_duration integer minimum sustained-excursion length (days)
#' forwarded to \code{calcDoyCrossThreshold} for the wet-season-end and
#' hot-day crossings (default 1 = off). See \code{?calcDoyCrossThreshold}.
#' @param smooth_window Integer day-window for the daily-climatology smoothing -- the single
#' global window applied to BOTH the crossing inputs (the hot-day \code{daily_temp} crossing,
#' and \code{daily_prec} / \code{daily_pet} before the wet-season-end P/PET ratio) AND the
#' reductions (\code{warmest_day}, driest-window P/PET, the \code{daily_ppet_diff} trend).
#' Default 31 (odd). See \code{calcCropCalendars}.
#' @param prev_always_wet Logical hysteresis state from last year (or \code{NA}): was
#' the cell always-wet (\code{min_ppet >= ppet_min})? Used only when
#' \code{harv_ppet_eps > 0}. The resolved value is returned as
#' \code{attr(., "always_wet")}.
#' @param harv_ppet_eps Non-negative relative deadband (default 0 = off) on the
#' always-wet test (\code{min_ppet} vs \code{ppet_min}): it relaxes \code{ppet_min}
#' toward keeping \code{prev_always_wet}, so a grazing P/PET wobble no longer snaps
#' \code{hd_wetseas} between \code{hd_last} and \code{hd_first}.
#' @param prev_rw1 Integer level wet-end EXISTENCE regime from last year (0 = absent-DRY,
#' 1 = found, 2 = absent-WET), or \code{NA}. Used only on the \code{cc_regime} path when
#' \code{harv_exist_eps > 0}.
#' @param harv_exist_eps Non-negative relative directional deadband (default 0 = off) on the
#' wet-end EXISTENCE boundary (the \code{ppet_ratio} level crossing, \code{doy_wet1}). Leaving
#' the absent-DRY regime needs the peak to clear \code{ppet_ratio * (1 + eps)}; leaving
#' absent-WET needs the trough below \code{ppet_ratio * (1 - eps)}. Experimental, gated by
#' \code{options(cc_regime = TRUE)}.
#'
#' @seealso getCropParam, calcClimatology, calcSowingDate, calcCropCalendars
#' @export
calcHarvestDateVector <- function(croppar,
                                  sowing_date,
                                  sowing_season,
                                  monthly_temp,
                                  monthly_ppet,
                                  monthly_ppet_diff,
                                  daily_temp = NULL,
                                  daily_prec = NULL,
                                  daily_pet  = NULL,
                                  cross_min_duration = 1L,
                                  smooth_window = 31L,
                                  prev_always_wet = NA,
                                  harv_ppet_eps = 0,
                                  prev_rw1 = NA_integer_,
                                  harv_exist_eps = 0
                                  ) {

  # Extract individual parameter names and values
  list2env(croppar, environment())  # 1-row data.frame: columns -> scalar params

  ndays_year <- 365

  # Daily climatologies (one value per DOY 1:365). dtemp/dprec/dpet from the cache
  # are already on this grid; without a daily temperature series, interpolate the
  # monthly means onto the grid for the hot-day crossings (which have no pure-monthly
  # form). have_dtemp also selects the legacy mid-day for warmest_day below.
  have_dtemp <- !is.null(daily_temp)
  if (!have_dtemp) daily_temp <- .monthlyToDoy365(monthly_temp)
  # Spike-free daily P/PET for the wet-season-end crossing: the ratio of the
  # per-DOY mean P and mean PET (never blows up — mean PET on a DOY is never ~0),
  # unlike the mean of daily P/PET ratios. This feeds a threshold CROSSING, so P and
  # PET are smoothed (smooth_window) BEFORE the ratio is formed -- the daily
  # climatology itself is kept raw, and only the crossing input is conditioned here
  # (0 = no-op). Fall back to the interpolated monthly ratio if daily P/PET absent.
  if (!is.null(daily_prec) && !is.null(daily_pet)) {
    daily_ppet <- .circRoll(daily_prec, smooth_window, "center", mean = TRUE) /
                  pmax(.circRoll(daily_pet, smooth_window, "center", mean = TRUE), 1e-6)
  } else {
    daily_ppet <- .monthlyToDoy365(monthly_ppet)
  }
  # P/PET month-over-month difference (declining-moisture trend). Use the daily
  # analogue -- the 30-day Sum P / Sum PET centred at each DOY minus the window
  # centred ~30 days later -- when the daily P and PET are supplied; this removes the
  # residual whole-month quantisation of the second wet-season-end candidate
  # (doy_wet2). Fall back to interpolating the 12 monthly differences only when the
  # daily series are absent.
  if (!is.null(daily_prec) && !is.null(daily_pet)) {
    daily_ppet_diff <- .dailyPpetDiff(daily_prec, daily_pet, width = smooth_window)
  } else {
    daily_ppet_diff <- .monthlyToDoy365(monthly_ppet_diff)
  }

  # NB: the wet-season-end crossing EXISTENCE is deliberately NOT deadbanded. A
  # one-directional nudge of ppet_ratio cannot keep a band-membership condition
  # (min(daily_ppet) < ppet_ratio < max) sticky -- raising it to keep a found crossing
  # instead eliminates the crossing whenever the wet peak sits just above ppet_ratio,
  # forcing a found/not-found 2-cycle. Measured on unbiased cells this manufactured
  # far MORE flicker than it removed, so the wet-end crossing is left raw.

  # Shortest cycle: crop lower biological limit
  hd_first <- sowing_date + min_growingseason
  # Medium cycle: best trade-off vegetative and reproductive growth
  hd_maxrp <- sowing_date + maxrp_growingseason
  # Longest cycle: crop upper biological limit
  hd_last <- ifelse(sowing_season == "winter",
                    sowing_date + max_growingseason_wt,
                    sowing_date + max_growingseason_st)

  # End of wet season ----
  doy_wet1 <- calcDoyCrossThreshold(
    daily_ppet,
    ppet_ratio,
    min_duration = cross_min_duration
    )[["doy_cross_down"]]
  doy_wet2 <- calcDoyCrossThreshold(
    daily_ppet_diff,
    ppet_ratio_diff,
    min_duration = cross_min_duration
    )[["doy_cross_down"]]
  # Wet-season-end estimates (P/PET level doy_wet1 + trend doy_wet2); the crop escapes
  # terminal water stress at the first valid one (a sub-minimum wet-end wraps to next year,
  # see the hd_wetseas branch). wet_ends is assembled inside that branch (the cc_regime path
  # may re-detect the crossings at deadbanded thresholds first).
  # If does not find harvest date and it is always high rainfall. Default: the
  # .driestWindowPpet over the smooth_window window (continuous; no month quantisation),
  # else the calendar-month minimum. cc_harmonize_minppet keys the always-wet test instead
  # to min(daily_ppet) -- the SAME variable as the wet-end EXISTENCE crossing -- closing the
  # ppet_min/ppet_ratio dead-zone. That only pays off if daily_ppet is well smoothed:
  # harmonizing at smooth_window=15 regressed +40% on Uruguay+France (the 15-day cycle is
  # too noisy); at smooth_window=30 min(daily_ppet) ~ the driest-window magnitude AND
  # consistent with the crossing (toggle default off; under A/B).
  min_ppet <- if (!is.null(daily_prec) && !is.null(daily_pet)) {
    if (isTRUE(getOption("cc_harmonize_minppet", FALSE))) min(daily_ppet)
    else .driestWindowPpet(daily_prec, daily_pet, width = smooth_window)
  } else min(monthly_ppet)
  # Escape harvest for a found wet season: wrap a wet-end to next year if its escape
  # harvest would be SUB-MINIMUM (wet_end + rphase < hd_first), not only if strictly
  # before sowing. A wet-end resolving to a sub-minimum season is the tail of the PREVIOUS
  # wet season, spuriously detected near the sowing DOY; the long (wrapped) season is the
  # correct one. This pulls the near-sowing oscillation band onto the long side so it no
  # longer flickers down to a short (hd_first) season -- the dominant ~50% of harvest jumps.
  wetseas_escape <- function(wet_ends) {
    wrap <- wet_ends + rphase_duration < hd_first
    min(ifelse(wrap, wet_ends + ndays_year, wet_ends)) + rphase_duration
  }
  # Always-wet test (min_ppet >= ppet_min) deciding hd_last vs hd_first when no wet-end
  # exists; flips (a big jump) when min_ppet grazes ppet_min. Thermostat: relax ppet_min
  # toward keeping last year's verdict by the given relative deadband.
  always_wet_test <- function(eps) {
    ppet_min_eff <- if (eps > 0 && !is.na(prev_always_wet))
      ppet_min * (if (prev_always_wet) 1 - eps else 1 + eps) else ppet_min
    min_ppet >= ppet_min_eff
  }

  if (isTRUE(getOption("cc_regime", FALSE))) {
    # DIRECTIONAL Schmitt-trigger deadband (harv_exist_eps) on the wet-end EXISTENCE
    # boundary only -- the ppet_ratio LEVEL crossing (doy_wet1). Last year's regime shifts
    # the bar asymmetrically so a grazing signal must move decisively to flip the existence
    # of the crossing, WITHOUT a single nudged threshold (which would destroy the very
    # crossing it tries to keep). doy_wet2 (TREND) is left plain (a deadband there was strictly
    # worse), and the always-wet test keeps its OWN deadband (harv_ppet_eps) -- ppet_min is a
    # separate aridity floor (Rice: ppet_ratio=1.0 but ppet_min=0.5), not the same threshold.
    eps <- harv_exist_eps
    # ppet_ratio -- LEVEL wet-end (doy_wet1). 3-state: 0 absent-DRY / 1 found / 2 absent-WET.
    # Leave dry: peak must clear ppet_ratio*(1+eps); leave wet: trough below ppet_ratio*(1-eps).
    thr1 <- if (is.na(prev_rw1) || eps == 0) ppet_ratio
            else if (prev_rw1 == 0L) ppet_ratio * (1 + eps)
            else if (prev_rw1 == 2L) ppet_ratio * (1 - eps)
            else                     ppet_ratio
    doy_wet1 <- calcDoyCrossThreshold(daily_ppet, thr1,
                                      min_duration = cross_min_duration)[["doy_cross_down"]]
    rw1 <- if (max(daily_ppet) < thr1) 0L else if (min(daily_ppet) >= thr1) 2L else 1L
    # always-wet test keeps its own independent deadband (NOT the regime eps).
    always_wet <- always_wet_test(harv_ppet_eps)

    wet_ends <- c(doy_wet1, doy_wet2); wet_ends <- wet_ends[wet_ends != -9999]
    hd_wetseas <- if (doy_wet1 == -9999) (if (always_wet) hd_last else hd_first)
                  else wetseas_escape(wet_ends)
    harv_hi <- rw1 + 6L * as.integer(always_wet)
  } else {
    always_wet <- always_wet_test(harv_ppet_eps)
    wet_ends <- c(doy_wet1, doy_wet2); wet_ends <- wet_ends[wet_ends != -9999]
    hd_wetseas <- if (doy_wet1 == -9999) (if (always_wet) hd_last else hd_first)
                  else wetseas_escape(wet_ends)
    harv_hi <- 6L * as.integer(always_wet)
  }

  # Warmest period of the year ----
  # Centre DOY of the warmest 30-day window of the daily climatology (the daily
  # analogue of the previous "mid-day of the warmest month"); the legacy monthly
  # fallback is exactly that mid-day when the daily series is absent.
  warmest_day <- if (have_dtemp) .doyWarmestWindow(daily_temp, width = smooth_window) else
    c(15, 43, 74, 104, 135, 165, 196, 227, 257, 288, 318, 349)[which.max(monthly_temp)]
  hd_temp_base <- ifelse(
    sowing_season == "winter", warmest_day, warmest_day + rphase_duration
    )

  # Smoothed daily temperature for the hot-day threshold crossings, using the same global
  # smooth_window as warmest_day and all other reductions.
  daily_temp_x <- .circRoll(daily_temp, smooth_window, "center", mean = TRUE)

  # First hot day ----
  doy_exceed_opt_rp <- calcDoyCrossThreshold(
    daily_temp_x, temp_opt_rphase, min_duration = cross_min_duration
    )[["doy_cross_up"]]
  idx <- which(doy_exceed_opt_rp < sowing_date & doy_exceed_opt_rp != -9999)
  doy_exceed_opt_rp[idx] <- doy_exceed_opt_rp[idx] + ndays_year
  doy_exceed_opt_rp      <- sort(doy_exceed_opt_rp)[1]

  # Last hot day ----
  doy_below_opt_rp <- calcDoyCrossThreshold(
    daily_temp_x,
    temp_opt_rphase,
    min_duration = cross_min_duration
    )[["doy_cross_down"]]
  idx <- which(doy_below_opt_rp < sowing_date & doy_below_opt_rp != -9999)
  doy_below_opt_rp[idx] <- doy_below_opt_rp[idx] + ndays_year
  doy_below_opt_rp      <- sort(doy_below_opt_rp)[1]

  # Winter type: First hot day; Spring type: Last hot day
  doy_opt_rp  <- ifelse(
    sowing_season == "winter", doy_exceed_opt_rp, doy_below_opt_rp
    )
  if (doy_opt_rp == -9999) {
    hd_temp_opt <- hd_maxrp
  } else {
    hd_temp_opt <- ifelse(
      sowing_season == "winter", doy_opt_rp, doy_opt_rp+rphase_duration
      )
  }

  # If harvest date < sowing date, it occurs the following year, so add 365 days
  hd_wetseas    <- ifelse(
    hd_wetseas < sowing_date, hd_wetseas + ndays_year, hd_wetseas
    )
  hd_temp_base  <- ifelse(
    hd_temp_base < sowing_date, hd_temp_base + ndays_year, hd_temp_base
    )
  hd_temp_opt   <- ifelse(
    hd_temp_opt < sowing_date, hd_temp_opt + ndays_year, hd_temp_opt
    )

  # hd_vector ----
  hd_vector        <- c(hd_first, hd_maxrp, hd_last,
                        hd_wetseas, hd_temp_base, hd_temp_opt)
  names(hd_vector) <- c("hd_first", "hd_maxrp", "hd_last",
                        "hd_wetseas", "hd_temp_base", "hd_temp_opt")

  # Moisture-state hysteresis carried forward by calcCropCalendars: the high part of the
  # packed harvest state. 0/1 (always-wet flag) on the default path; 0/1/2 (DRY/NORMAL/WET
  # regime) on the cc_regime path.
  attr(hd_vector, "harv_hi") <- harv_hi
  return(hd_vector)
}

# Map 12 monthly values (referring to month mid-days) to a clean 365-element
# vector indexed by DOY 1:365, via linear interpolation with Dec->Jan wrap.
.monthlyToDoy365 <- function(monthly_value) {
  d   <- interpolateMonthlyToDaily(monthly_value)
  doy <- ((d[["x"]] - 1L) %% 365L) + 1L
  # interpolateMonthlyToDaily covers each DOY (often twice, from the two-year
  # replication); average duplicates and order by DOY.
  as.numeric(tapply(d[["y"]], doy, mean)[as.character(1:365)])
}

# Centre DOY of the warmest `width`-day window of a daily (DOY-indexed)
# climatology, evaluated circularly. Daily analogue of "warmest month mid-day".
# The centred rolling sum is already indexed by window centre, so its argmax IS the
# centre DOY (no manual half-width shift).
.doyWarmestWindow <- function(daily_value, width = 30) {
  as.integer(which.max(.circRoll(daily_value, width, "center")))
}
