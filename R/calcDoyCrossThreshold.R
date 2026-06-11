#' @title Calculate day of crossing threshold
#'
#' @description Find the DOY where a climatological daily variable crosses a
#' threshold, evaluated on actual daily data (365 values per DOY, averaged
#' across years). The comparison wraps circularly so crossings between
#' December 31 and January 1 are detected.
#'
#' @param daily_value Numeric vector of length 365. Climatological daily
#'   values of the variable (e.g. temperature), one per DOY, averaged across
#'   the reference period. Typically the \code{dtemp} element returned by
#'   \code{calcMonthlyClimate}.
#' @param threshold Numeric. Threshold value for crossing detection.
#'
#' @return Named list with two elements:
#'   \describe{
#'     \item{doy_cross_up}{First DOY where \code{daily_value} crosses upward
#'       through \code{threshold} (i.e. was below on the previous day, at or
#'       above today). \code{-9999} if no upward crossing exists.}
#'     \item{doy_cross_down}{First DOY where \code{daily_value} crosses
#'       downward through \code{threshold}. \code{-9999} if none.}
#'   }
#' @export

calcDoyCrossThreshold <- function(daily_value, threshold) {

  n <- length(daily_value)

  is_above      <- daily_value >= threshold
  is_above_prev <- c(is_above[n], is_above[1:(n - 1)])

  cross <- as.integer(is_above) - as.integer(is_above_prev)

  day_cross_up   <- which(cross ==  1L)[1]
  day_cross_down <- which(cross == -1L)[1]

  if (is.na(day_cross_up))   day_cross_up   <- -9999L
  if (is.na(day_cross_down)) day_cross_down <- -9999L

  return(list("doy_cross_up"   = day_cross_up,
              "doy_cross_down" = day_cross_down))
}
