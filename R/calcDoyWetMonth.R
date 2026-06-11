#' @title Calculate beginning of the rainy season
#'
#' @description Find the start of the 120 wettest days as the DOY where the
#' 120-day cumulative P/PET sum is maximised, evaluated on a climatological
#' daily P/PET cycle (365 values, one per DOY, averaged across years).
#' The rolling window wraps circularly across the year boundary.
#'
#' @param daily_ppet Numeric vector of length 365. Climatological daily
#'   precipitation-to-PET ratio (P/PET), one value per DOY, averaged across
#'   the reference period. Typically the \code{dppet} element returned by
#'   \code{calcMonthlyClimate}.
#'
#' @return Integer DOY (1–365) of the start of the 120-day wettest window.
#' @export

calcDoyWetMonth <- function(daily_ppet) {

  n <- length(daily_ppet)
  x <- vapply(0:(n - 1),
              function(i) sum(daily_ppet[((0:119 + i) %% n) + 1]),
              numeric(1))
  return(which.max(x))

}
