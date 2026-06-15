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
#' @param monthly_ppet numeric vestor of length 12. Mean Potential Evapotranspiration (mm). See caclMonthlyClimate.
#' @param monthly_ppet_diff numeric vestor of length 12. Mean difference of
#' Potential Evapotranspiration (mm). See caclMonthlyClimate.
#' @param daily_temp numeric vector of length 365. Climatological daily mean
#' temperature (deg C), one value per DOY (the \code{dtemp} element from
#' \code{calcMonthlyClimate}). If \code{NULL}, it is interpolated from
#' \code{monthly_temp}.
#' @param daily_prec numeric vector of length 365. Climatological daily
#' precipitation (mm) per DOY (the \code{dprec} element from
#' \code{calcMonthlyClimate}). Used with \code{daily_pet} to form the spike-free
#' daily P/PET (ratio of per-DOY means) for the wet-season-end crossing. If
#' \code{NULL}, the daily P/PET is interpolated from \code{monthly_ppet} instead.
#' @param daily_pet numeric vector of length 365. Climatological daily PET (mm)
#' per DOY (the \code{dpet} element from \code{calcMonthlyClimate}).
#' @param cross_min_duration integer minimum sustained-excursion length (days)
#' forwarded to \code{calcDoyCrossThreshold} for the wet-season-end and
#' hot-day crossings (default 1 = off). See \code{?calcDoyCrossThreshold}.
#'
#' @seealso getCropParam, calcMonthlyClimate, calcSowingDate, calcCropCalendars
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
                                  cross_min_duration = 1L
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
  # unlike the mean of daily P/PET ratios. Fall back to the interpolated monthly
  # ratio if the daily P and PET are not supplied.
  if (!is.null(daily_prec) && !is.null(daily_pet)) {
    daily_ppet <- daily_prec / pmax(daily_pet, 1e-6)
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
    daily_ppet_diff <- .dailyPpetDiff(daily_prec, daily_pet)
  } else {
    daily_ppet_diff <- .monthlyToDoy365(monthly_ppet_diff)
  }

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
  doy_wet_vec <- ifelse(
    c(doy_wet1, doy_wet2) < sowing_date & c(doy_wet1, doy_wet2) != -9999,
    c(doy_wet1, doy_wet2) + ndays_year,
    c(doy_wet1, doy_wet2)
    )
  # If more than one wet seasons take the first one, else -9999
  doy_wet_first <- ifelse(
    length(doy_wet_vec[doy_wet_vec != -9999]) > 0,
    min(doy_wet_vec[doy_wet_vec != -9999]),
    -9999
    )
  # If does not find harvest date and it is always high rainfall. The "driest month"
  # P/PET uses the driest 30-day window of the daily Sum P / Sum PET when available
  # (continuous; no month quantisation), else the calendar-month minimum.
  min_ppet <- if (!is.null(daily_prec) && !is.null(daily_pet))
    .driestWindowPpet(daily_prec, daily_pet) else min(monthly_ppet)
  if (doy_wet1 == -9999) {
    if (min_ppet >= ppet_min) {
      hd_wetseas <- hd_last
    } else {
      hd_wetseas <- hd_first
    }
  } else {
    hd_wetseas <- doy_wet_first + rphase_duration
  }

  # Warmest period of the year ----
  # Centre DOY of the warmest 30-day window of the daily climatology (the daily
  # analogue of the previous "mid-day of the warmest month"); the legacy monthly
  # fallback is exactly that mid-day when the daily series is absent.
  warmest_day <- if (have_dtemp) .doyWarmestWindow(daily_temp, width = 30) else
    c(15, 43, 74, 104, 135, 165, 196, 227, 257, 288, 318, 349)[which.max(monthly_temp)]
  hd_temp_base <- ifelse(
    sowing_season == "winter", warmest_day, warmest_day + rphase_duration
    )

  # First hot day ----
  doy_exceed_opt_rp <- calcDoyCrossThreshold(
    daily_temp, temp_opt_rphase, min_duration = cross_min_duration
    )[["doy_cross_up"]]
  idx <- which(doy_exceed_opt_rp < sowing_date & doy_exceed_opt_rp != -9999)
  doy_exceed_opt_rp[idx] <- doy_exceed_opt_rp[idx] + ndays_year
  doy_exceed_opt_rp      <- sort(doy_exceed_opt_rp)[1]

  # Last hot day ----
  doy_below_opt_rp <- calcDoyCrossThreshold(
    daily_temp,
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
.doyWarmestWindow <- function(daily_value, width = 30) {
  n     <- length(daily_value)
  start <- which.max(.circRollSum(daily_value, width))  # first DOY of warmest window
  ((start - 1L + width %/% 2L) %% n) + 1L                # centre DOY
}
