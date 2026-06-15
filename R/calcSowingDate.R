#' @title Calculate sowing date (Waha et al., 2012)
#'
#' @param monthly_temp Numeric vector of length 12. Average (e.g. 20-years)
#' monthly mean temperatures (degree Celsius). Retained as a fallback: the
#' coldest-month temperature and the coldest-day anchor are taken from
#' \code{daily_temp} (coldest 30-day window) when supplied, and only derived from
#' these monthly values if \code{daily_temp} is \code{NULL}.
#' @param daily_prec Numeric vector of length 365. Climatological daily
#' precipitation (mm, one value per DOY, averaged across years). The
#' \code{dprec} element from \code{calcClimatology}. Used by
#' \code{calcDoyWetMonth} (with \code{daily_pet}) to find the start of the wet
#' season as the wettest 120-day window (\eqn{\sum P / \sum PET}). Preferred input
#' for the wettest-window start in PREC/PRECTEMP cells; if absent, the monthly
#' fallback (\code{monthly_prec} / \code{monthly_pet}) is used instead.
#' @param daily_pet Numeric vector of length 365. Climatological daily PET (mm,
#' one value per DOY). The \code{dpet} element from \code{calcClimatology}.
#' Paired with \code{daily_prec} (see there).
#' @param daily_temp Numeric vector of length 365. Climatological daily mean
#' temperature (°C), one value per DOY (the \code{dtemp} element from
#' \code{calcClimatology}). Used for the coldest-month reductions (coldest
#' 30-day window mean and centre DOY) and by \code{calcDoyCrossThreshold} for the
#' spring/fall threshold crossings. \code{NULL} falls back to the legacy monthly
#' rule (calendar-month minimum and coldest-month mid-day) and interpolates the
#' monthly means onto the daily grid for the crossings.
#' @param seasonality character value indicating the seasonality type as
#' computed by calcSeasonality
#' @param prev_wet_doy Integer DOY of last year's wettest-window start (or
#' \code{NA}), forwarded to \code{calcDoyWetMonth} for temporal hysteresis in the
#' PREC/PRECTEMP sowing branch. Only used when \code{wet_window_eps > 0}.
#' @param wet_window_eps Non-negative relative hysteresis band forwarded to
#' \code{calcDoyWetMonth} (default 0 = off). See \code{?calcDoyWetMonth}.
#' @param wet_window_decay Positive Gaussian decay scale forwarded to
#' \code{calcDoyWetMonth} (default 0.3; smaller = stickier). See
#' \code{?calcDoyWetMonth}.
#' @param cross_min_duration Integer minimum sustained-excursion length (days)
#' forwarded to \code{calcDoyCrossThreshold} for the spring/fall temperature
#' crossings (default 1 = off). See \code{?calcDoyCrossThreshold}.
#' @param cross_smooth_window Integer odd day-window for smoothing \code{daily_temp}
#' before the spring/fall threshold crossings only (the coldest-month reductions use
#' the raw series). Default \code{0} = off. The daily climatology is kept raw and the
#' crossing input is conditioned here, so the smoothing does not affect the
#' reductions or the 120-day wettest window.
#' @param monthly_prec,monthly_pet Optional numeric vectors of length 12, the
#' monthly precipitation and PET TOTALS (the \code{mprec} / \code{mpet} elements
#' from \code{calcClimatology}). Used only as the bug-free monthly fallback for
#' the wettest-window start in PREC/PRECTEMP cells, when
#' \code{daily_prec}/\code{daily_pet} (and \code{wet_doy}) are not supplied: the
#' wettest 4-month ratio-of-sums (\eqn{\sum P / \sum PET}). Coarser (~1-month
#' resolution) than the daily 120-day rule, but identical in form (sums P and PET
#' separately, no daily interpolation).
#' @param wet_doy Optional pre-computed wettest-window start DOY (the crop-independent
#' \code{calcDoyWetMonth} result). When supplied it is used directly in the
#' PREC/PRECTEMP spring-sowing branch, so the caller (\code{calcCropCalendars}) can
#' compute it once per cell instead of twice. \code{NULL} (default) recomputes it.
#' @export
calcSowingDate <- function(croppar,
                           monthly_temp,
                           daily_prec            = NULL,
                           daily_pet             = NULL,
                           daily_temp            = NULL,
                           seasonality,
                           lat,
                           prev_wet_doy          = NA_integer_,
                           wet_window_eps        = 0,
                           wet_window_decay      = 0.3,
                           cross_min_duration    = 1L,
                           wet_doy               = NULL,
                           monthly_prec          = NULL,
                           monthly_pet           = NULL,
                           cross_smooth_window   = 0L
                           ) {

  # extract individual parameter names and values
  list2env(croppar, environment())  # 1-row data.frame: columns -> scalar params

  # Coldest-month temperature (coldest_t) and coldest-day anchor (coldest_doy).
  # With the daily climatology these are the coldest 30-day window mean and its
  # centre DOY: a daily mean is already a per-DOY 30-year mean (smoothed by
  # cross_smooth_window), so the 30-day window reproduces the calendar-month coldness
  # continuously -- removing the ~30-day quantisation that made the warm-winter
  # sowing date and the spring-crossing anchor jump between adjacent climate windows.
  # When daily_temp is absent, fall back to the exact legacy monthly rule (the
  # calendar-month minimum and the coldest-month mid-day), and interpolate the
  # monthly means onto the daily grid only for the threshold-crossing scans below
  # (the spring/fall crossings have no pure-monthly form).
  if (!is.null(daily_temp)) {
    coldest_t   <- .coldestWindowMean(daily_temp)
    coldest_doy <- .doyColdestWindow(daily_temp)
    warmest_doy <- .doyWarmestWindow(daily_temp)
  } else {
    midday      <- c(15, 43, 74, 104, 135, 165, 196, 227, 257, 288, 318, 349)
    coldest_t   <- min(monthly_temp)
    coldest_doy <- midday[which.min(monthly_temp)]
    warmest_doy <- midday[which.max(monthly_temp)]
    daily_temp  <- .monthlyToDoy365(monthly_temp)
  }

  # Smoothed daily temperature for the threshold crossings ONLY (the coldest-month
  # reductions above use the raw daily_temp). The daily climatology is kept raw; this
  # is where cross_smooth_window conditions the crossing input (0 = no-op).
  daily_temp_x <- .smoothCycle(daily_temp, cross_smooth_window)

  # Constrain first possible date for winter crop sowing
  earliest_sdate  <- ifelse(lat >= 0, initdate.sdatenh, initdate.sdatesh)
  earliest_smonth <- doy2month(earliest_sdate)
  DEFAULT_DOY     <- ifelse(lat >= 0, 1, 182)
  DEFAULT_MONTH   <- 0

  # What type of winter is it?
  if ((coldest_t > basetemp.low) &
      (seasonality %in% c("TEMP", "TEMPPREC", "PRECTEMP", "PREC"))) {
    # "Warm winter" (allowing non-vernalizing winter-sown crops)
    # sowing 2.5 months before the coldest day
    # it seems a good approximation for both India and South US)
    coldestday     <- coldest_doy
    firstwinterdoy <- ifelse(coldestday-75<=0, coldestday-75+365, coldestday-75)

  } else if ((coldest_t < -10) &
             (seasonality %in% c("TEMP", "TEMPPREC", "PRECTEMP", "PREC"))) {
    # "Cold winter" (winter too harsh for winter crops, only spring sowing possible)
    firstwinterdoy <- -9999

  } else {
    # "Mild winter" (allowing vernalizing crops). Anchor the autumn down-crossing scan
    # at the warmest day -- the mirror of the coldest-day anchor on the spring
    # up-crossing below. Scanning forward from the summer peak catches the genuine
    # autumn cooling and skips a non-monotonic spring ascent grazing temp_fall, which
    # from DOY 1 would register a spurious autumn-cooling date ~half a year early
    # (the same failure the spring fix removed, mirrored).
    firstwinterdoy <- calcDoyCrossThreshold(
      daily_temp_x, temp_fall, min_duration = cross_min_duration,
      from = warmest_doy)[["doy_cross_down"]]

  }

  # First day of winter
  firstwintermonth <- ifelse(
    firstwinterdoy == -9999, DEFAULT_MONTH, doy2month(firstwinterdoy)
    )
  firstwinterdoy   <- ifelse(
    firstwinterdoy == -9999, DEFAULT_DOY, firstwinterdoy
    )

  # First day of spring: the first upward temp crossing AFTER the coldest day.
  # Anchoring the scan to the (very stable) winter minimum skips an autumn
  # temperature plateau grazing temp_spring -- which sits before the coldest day --
  # that otherwise produces a spurious ~half-year-early sowing in mild-winter cells
  # (e.g. Uruguay / S. Brazil), without adding any temporal lag. coldest_doy is the
  # centre of the coldest 30-day window of the daily climatology (computed above).
  firstspringdoy   <- calcDoyCrossThreshold(
    daily_temp_x, temp_spring, min_duration = cross_min_duration,
    from = coldest_doy)[["doy_cross_up"]]
  firstspringmonth <- ifelse(
    firstspringdoy == -9999, DEFAULT_MONTH, doy2month(firstspringdoy)
    )

  # If winter type
  if (calcmethod_sdate == "WTYP_CALC_SDATE") {

    if (firstwinterdoy > earliest_sdate &
        firstwintermonth != DEFAULT_MONTH) {

      sowing_month  <- firstwintermonth
      sowing_doy    <- firstwinterdoy
      sowing_season <- "winter"

    } else if (firstwinterdoy <= earliest_sdate &
               coldest_t > temp_fall &
               firstwintermonth != DEFAULT_MONTH) {

      sowing_month  <- earliest_smonth
      sowing_doy    <- earliest_sdate
      sowing_season <- "winter"

    } else {

      sowing_month  <- firstspringmonth
      sowing_doy    <- firstspringdoy
      sowing_season <- "spring"

    }

  } else {

    if (seasonality == "NO_SEASONALITY") {

      sowing_month <- DEFAULT_MONTH
      sowing_doy <- DEFAULT_DOY
      sowing_season <- "spring"

    } else if (seasonality == "PREC" || seasonality == "PRECTEMP") {

      # The wettest-window DOY is crop-independent, so calcCropCalendars computes
      # it once and passes it in; only recompute if not supplied (direct callers).
      # Prefer the daily 120-day SUM P / SUM PET rule (calcDoyWetMonth). Without the
      # daily series, fall back to the monthly 4-month ratio-of-sums on the monthly P
      # and PET TOTALS (.wetDoyMonthly) -- bug-free (it sums P and PET separately,
      # never the monthly P/PET ratios that a near-zero-PET month makes explode) and
      # using no daily interpolation, just coarser (~1-month resolution).
      if (is.null(wet_doy)) {
        if (!is.null(daily_prec) && !is.null(daily_pet)) {
          wet_doy <- calcDoyWetMonth(daily_prec, daily_pet,
                                     prev_doy = prev_wet_doy, eps = wet_window_eps,
                                     decay = wet_window_decay)
        } else if (!is.null(monthly_prec) && !is.null(monthly_pet)) {
          wet_doy <- .wetDoyMonthly(monthly_prec, monthly_pet)
        } else {
          stop("PREC/PRECTEMP sowing needs P and PET: supply daily_prec/daily_pet ",
               "(preferred, 120-day window), monthly_prec/monthly_pet (4-month ",
               "fallback), or a precomputed wet_doy.")
        }
      }
      sowing_doy    <- wet_doy
      sowing_month  <- doy2month(sowing_doy)
      sowing_season <- "spring"

    } else {

      sowing_month <- firstspringmonth
      sowing_doy <- firstspringdoy
      sowing_season <- "spring"

    }
  } # STYP_CALC_SDATE

  sowing_doy    <- ifelse(
    sowing_month == DEFAULT_MONTH, DEFAULT_DOY, sowing_doy
    )
  sowing_month  <- ifelse(
    calcmethod_sdate == "WTYP_CALC_SDATE" & sowing_season == "spring",
    DEFAULT_MONTH, sowing_month
    )

  sd_vector <- list("sowing_month"  = sowing_month,
                    "sowing_doy"    = sowing_doy,
                    "sowing_season" = sowing_season)

  return(sd_vector)
}

