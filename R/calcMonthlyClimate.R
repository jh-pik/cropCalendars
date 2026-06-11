#' @title Calculate monthly climate
#'
#' @description Calculate monthly climate variables needed for the computation
#' of the rule-based crop calendars (Waha et al., 2012; Minoli et al., 2019):
#' average monthly mean temperature (mtemp);
#' average monthly cumulative precipitation (mprec);
#' average monthly cumulative potential evapotranspiration (mpet);
#' dryness index 1, mprec-to-mpet ratio (mppet);
#' dryness index 2, difference of mppet of two consecutive months (mppet_diff).
#' P-to-PET (mppet) ratio indicates the water surplus or deficit with respect to
#' the plant water demand; P-to-PET ratio difference (ppet_diff) indicates the
#' monthly trend in moisture conditions, if mppet_diff[m] > 0,
#' the trend is declining, indicating that the following month (m + 1) is
#' dryer than month m.
#'
#' @param lat latitude (decimal value)
#' @param temp daily mean temperature (degree Celsius) for a number of years
#' (syear:eyear). It should be passed in form of a vector.
#' @param prec daily precipitation (mm) for a number of years (syear:eyear).
#' It should be passed in form of a vector.
#' @param syear start year in the climate time series.
#' @param eyear end year in the climate time series.
#' @param incl_feb29 Does the time series include February 29th in leap years?
#' @param pet_method PET method: \code{"pt"} (default) for Priestley-Taylor
#'   equilibrium ET (calcPET), or \code{"fao56"} for FAO-56 Penman-Monteith
#'   reference ET (calcPET_FAO56).
#' @param swdown  surface downwelling shortwave radiation (W/m2, 24 h mean),
#'   vector matching length of \code{temp}. Used by both methods when supplied;
#'   if NULL with \code{pet_method = "pt"} falls back to orbital-geometry Rn.
#' @param lwdown  surface downwelling longwave radiation (W/m2, 24 h mean),
#'   vector matching length of \code{temp}. Required when \code{swdown} is given.
#' @param windspeed near-surface wind speed at 10 m (m/s). Required for
#'   \code{pet_method = "fao56"}.
#' @param humid   near-surface specific humidity (kg/kg, ISIMIP3b variable
#'   \code{huss}). Required for \code{pet_method = "fao56"}.
#' @param ps      surface air pressure (Pa). Used only with
#'   \code{pet_method = "fao56"}; defaults to 101325 Pa (sea level).
#'   \code{lat} and \code{day} (used internally for date handling and the PT
#'   orbital path) are not passed to \code{calcPET_FAO56}, which uses full
#'   24 h fluxes and requires no geometric daylength.
#'
#' @return list with five monthly vectors (length 12) and two daily vectors
#' (length 365):
#' \describe{
#'   \item{mtemp, mprec, mpet, mppet, mppet_diff}{Monthly climate (length 12).}
#'   \item{dtemp}{Climatological daily mean temperature (°C), one value per
#'     DOY 1–365, averaged across years. Used by \code{calcDoyCrossThreshold}.}
#'   \item{dppet}{Climatological daily P/PET ratio, one value per DOY 1–365,
#'     averaged across years. Used by \code{calcDoyWetMonth}.}
#' }
#'
#' @examples
#' d_temp <- matrix(rnorm(365*3, 15), nrow = 3)
#' d_prec <- matrix(rnorm(365*3, 3, 50), nrow = 3)
#' d_prec[d_prec <= 0] <- 0
#' calcMonthlyClimate(lat = 45, temp = d_temp, prec = d_prec,
#'                    syear = 2001, eyear = 2003)
#' @export

calcMonthlyClimate <- function(lat        = NULL,
                               temp       = NULL,
                               prec       = NULL,
                               syear      = NULL,
                               eyear      = NULL,
                               incl_feb29 = TRUE,
                               pet_method = c("pt", "fao56"),
                               swdown     = NULL,
                               lwdown     = NULL,
                               windspeed  = NULL,
                               humid      = NULL,
                               ps         = 101325
                               ) {

  pet_method <- match.arg(pet_method)

  years   <- syear:eyear
  nyears  <- length(years)
  nmonths <- 12

  if (incl_feb29 == FALSE) {
    dates <- createDateSeq(nstep = 365, years = years)

  } else {
    dates <- seqDates(start_date = paste0(syear, "-01-01"),
                      end_date   = paste0(eyear, "-12-31"),
                      step       = "day")
  }
  y_dates <- date_to_year(dates)
  m_dates <- date_to_month(dates)
  d_dates <- date_to_doy(dates, skip_feb29 = TRUE)

  # Compute daily PET
  if (pet_method == "fao56") {

    required <- list(windspeed = windspeed, humid = humid,
                     swdown = swdown, lwdown = lwdown)
    missing_vars <- names(which(sapply(required, is.null)))
    if (length(missing_vars) > 0)
      stop("pet_method = 'fao56' requires: ", paste(missing_vars, collapse = ", "))

    pet <- mapply(calcPET_FAO56,
                  temp      = temp,
                  windspeed = windspeed,
                  humid     = humid,
                  swdown    = swdown,
                  lwdown    = lwdown,
                  ps        = ps)

  } else {

    if (!is.null(swdown) && !is.null(lwdown)) {
      pet <- mapply(calcPET, temp = temp, lat = lat, day = d_dates,
                    swdown = swdown, lwdown = lwdown)
    } else {
      pet <- mapply(calcPET, temp = temp, lat = lat, day = d_dates)
    }

  }

  # Compute monthly climate for each year
  mtemp_y <- array(NA, dim = c(nyears, nmonths))
  mprec_y <- mpet_y <- mppet_y <- mtemp_y

  for (yy in seq_len(nyears)) {
    for (mm in seq_len(nmonths)) {

      idx <- which(y_dates == years[yy] & m_dates == mm)

      mtemp_y[yy, mm] <- mean(temp[idx])
      mprec_y[yy, mm] <- sum(prec[idx])
      mpet_y[yy, mm]  <- sum(pet[idx])
      mppet_y[yy, mm] <- mprec_y[yy, mm] / mpet_y[yy, mm]

    }
  }

  mtemp      <- round(apply(mtemp_y, 2, mean), digits = 5)
  mprec      <- round(apply(mprec_y, 2, mean), digits = 5)
  mpet       <- round(apply(mpet_y,  2, mean), digits = 5)
  mppet      <- round(apply(mppet_y, 2, mean), digits = 5)
  mppet_diff <- mppet - c(mppet[-1], mppet[1])

  names(mtemp) <- names(mprec) <- c("month" = seq_len(12))
  names(mpet)  <- names(mppet) <- c("month" = seq_len(12))
  names(mppet_diff) <- c("month" = seq_len(12))

  # Climatological daily arrays (365 values, one per DOY, averaged across years)
  ppet_daily <- prec / pmax(pet, 1e-6)
  dtemp <- as.vector(tapply(temp,       d_dates, mean))
  dppet <- as.vector(tapply(ppet_daily, d_dates, mean))

  return(list(mtemp      = mtemp,
              mprec      = mprec,
              mpet       = mpet,
              mppet      = mppet,
              mppet_diff = mppet_diff,
              dtemp      = dtemp,
              dppet      = dppet
              )
         )
}
