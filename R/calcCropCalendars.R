#' @title Calculate crop calendars
#'
#' @description Wrapper function calling sub-functions to
#' calculate sowing and harvest dates.
#' @param lon Longitude (decimal degrees)
#' @param lat Latitude (decimal degrees)
#' @param mclimate Monthly climate. A list returned by the calcMonthlyClimate()
#' function.
#' @param crop A crop name (chr), among those specified in the croppar_file
#' @param croppar_file Crop parameter file. If not specified, the default one is
#' used.
#' @param crop_parameters Optional pre-extracted crop-parameter row (one row
#' \code{data.frame}, as returned by \code{getCropParam}). When supplied, the
#' parameter file is not read — essential when calling this per cell in a loop
#' (avoids re-reading the CSV on every call). \code{crop}/\code{croppar_file}
#' are then ignored.
#' @param prev_wet_doy Integer DOY of last year's wettest-window start (or
#' \code{NA}), forwarded to \code{calcSowingDate}/\code{calcDoyWetMonth} for
#' temporal hysteresis (sliding-window pipeline). Only used when
#' \code{wet_window_eps > 0}. The wettest-window argmax is crop-independent, so
#' this state is per-cell. The resolved DOY is also returned as
#' \code{attr(., "wet_doy")} for the caller to carry forward.
#' @param wet_window_eps Non-negative relative hysteresis band forwarded to
#' \code{calcDoyWetMonth} (default 0 = off). See \code{?calcDoyWetMonth}.
#' @param wet_window_drift_gate Integer day distance within which wettest-window
#' moves are always followed (forwarded to \code{calcDoyWetMonth}; default 7).
#' @param cross_min_duration Integer minimum sustained-excursion length (days)
#' forwarded to \code{calcSowingDate}/\code{calcHarvestDateVector} ->
#' \code{calcDoyCrossThreshold} (default 1 = off). See
#' \code{?calcDoyCrossThreshold}.
#' @param prev_seas Character seasonality class from last year (or \code{NA}),
#' forwarded to \code{calcSeasonality} for class hysteresis. The class is
#' crop-independent, so this is per-cell state; the resolved class is returned as
#' \code{attr(., "seas_type")} for the caller to carry forward. See
#' \code{?calcSeasonality}.
#' @param seas_eps Non-negative seasonality-threshold deadband forwarded to
#' \code{calcSeasonality} (default 0 = off).
#' @seealso calcMonthlyClimate
#' @export

calcCropCalendars <- function(lon                   = NULL,
                              lat                   = NULL,
                              mclimate              = NULL,
                              crop                  = NULL,
                              croppar_file          = NULL,
                              crop_parameters       = NULL,
                              prev_wet_doy          = NA_integer_,
                              wet_window_eps        = 0,
                              wet_window_drift_gate = 7L,
                              cross_min_duration    = 1L,
                              prev_seas             = NA_character_,
                              seas_eps              = 0
                              ) {

  # Import crop parameters (unless already supplied by the caller).
  if (is.null(crop_parameters)) {
    if (is.null(croppar_file)) {
      croppar_file <- system.file("extdata", "crop_parameters.csv",
                                  package = "cropCalendars", "mustWork" = TRUE)
    }
    crop_parameters <- getCropParam(
      crops          = crop,
      cropparam_file = croppar_file
    )
  }

  # Get weather data of the grid cell
  mtemp      <- mclimate$mtemp
  mprec      <- mclimate$mprec
  mppet      <- mclimate$mppet
  mppet_diff <- mclimate$mppet_diff
  dtemp      <- mclimate$dtemp
  dprec      <- mclimate$dprec
  dpet       <- mclimate$dpet

  # Seasonality type
  seasonality <- calcSeasonality(
    monthly_temp = mtemp,
    monthly_prec = mprec,
    temp_min     = 10,
    prev_seas    = prev_seas,
    seas_eps     = seas_eps
  )

  # Sowing date
  sowing <- calcSowingDate(
    croppar        = crop_parameters,
    monthly_temp   = mtemp,
    daily_prec     = dprec,
    daily_pet      = dpet,
    daily_temp     = dtemp,
    seasonality           = seasonality,
    lat                   = lat,
    prev_wet_doy          = prev_wet_doy,
    wet_window_eps        = wet_window_eps,
    wet_window_drift_gate = wet_window_drift_gate,
    cross_min_duration    = cross_min_duration
  )

  sowing_month  <- sowing[["sowing_month"]]
  sowing_day    <- sowing[["sowing_doy"]]
  sowing_season <- sowing[["sowing_season"]]

  # Resolved wettest-window start to carry forward as hysteresis state. It only
  # drives sowing in PREC/PRECTEMP cells (where sowing_doy IS the wettest-window
  # DOY); elsewhere there is no wet-window state, so report NA.
  wet_doy <- if (seasonality %in% c("PREC", "PRECTEMP")) as.integer(sowing_day) else NA_integer_

  # Harvest date
  harvest_rule  <- calcHarvestRule(
    croppar      = crop_parameters,
    monthly_temp = mtemp,
    monthly_ppet = mppet,
    seasonality  = seasonality
  )

  harvest_vector <- calcHarvestDateVector(
    croppar           = crop_parameters,
    sowing_date       = sowing_day,
    sowing_season     = sowing_season,
    monthly_temp      = mtemp,
    monthly_ppet      = mppet,
    monthly_ppet_diff = mppet_diff,
    daily_temp        = dtemp,
    daily_prec        = dprec,
    daily_pet         = dpet,
    cross_min_duration = cross_min_duration
  )

  harvest <- calcHarvestDate(
    croppar       = crop_parameters,
    monthly_temp  = mtemp,
    sowing_date   = sowing_day,
    sowing_month  = sowing_month,
    sowing_season = sowing_season,
    seasonality   = seasonality,
    harvest_rule  = harvest_rule,
    hd_vector     = harvest_vector
  )

  harvest_day_rf  <- harvest[["hd_rf"]]
  harvest_day_ir  <- harvest[["hd_ir"]]
  harvest_reas_rf <- names(harvest[["harvest_reason_rf"]])
  harvest_reas_ir <- names(harvest[["harvest_reason_ir"]])

  # Growing period length
  growpriod_rf <- calcGrowingPeriod(sowing_day, harvest_day_rf, 365)
  growpriod_ir <- calcGrowingPeriod(sowing_day, harvest_day_ir, 365)

  # Output table
  pixel_df <- data.frame(
    "lon"              = rep(lon, 2),
    "lat"              = rep(lat, 2),
    "crop"             = rep(crop_parameters$crop_name, 2),
    "irrigation"       = c("Rainfed", "Irrigated"),
    "seasonality_type" = rep(seasonality, 2),
    "sowing_season"    = rep(sowing_season, 2),
    "sowing_month"     = rep(sowing_month, 2),
    "sowing_doy"       = rep(sowing_day, 2),
    "harvest_rule"     = rep(names(harvest_rule), 2),
    "harvest_reason"   = c(harvest_reas_rf, harvest_reas_ir),
    "maturity_doy"     = c(harvest_day_rf, harvest_day_ir),
    "growing_period"   = c(growpriod_rf, growpriod_ir)
  )

  attr(pixel_df, "wet_doy")   <- wet_doy
  attr(pixel_df, "seas_type") <- seasonality   # crop-independent; carried as prev_seas state
  return(pixel_df)

}
