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
#'   \code{calcClimatology}.
#' @param threshold Numeric. Threshold value for crossing detection.
#' @param min_duration Integer >= 1. Minimum number of consecutive days the
#'   excursion must persist for a crossing to count (a "sustained-excursion"
#'   guard). With \code{min_duration = 1L} (default) the behaviour is unchanged:
#'   the first day-to-day sign change is returned. With larger values, a crossing
#'   is only accepted if the variable stays on the new side of \code{threshold}
#'   for at least \code{min_duration} days (evaluated circularly), so a single-day
#'   noise blip in the daily climatology cannot register as a spurious crossing
#'   (e.g. a 1-day dip-and-recover through \code{temp_spring} in the autumn
#'   descent producing a sowing date ~130 days early).
#' @param from Integer DOY (1..365, default 1) at which the circular scan for the
#'   "first" crossing begins. The first candidate at or after \code{from} (wrapping
#'   across the year boundary) is returned. The default \code{from = 1L} reproduces
#'   the original "first crossing from Jan 1" behaviour. Anchoring the scan to a
#'   stable phenological reference removes the spurious-crossing ambiguity in mild
#'   climates: e.g. for the spring up-crossing, starting from the coldest day skips
#'   an autumn temperature plateau grazing the threshold (which sits \emph{before}
#'   the winter minimum) and lands on the genuine spring warming crossing -- without
#'   introducing any temporal lag.
#'
#' @return Named list with two elements:
#'   \describe{
#'     \item{doy_cross_up}{First DOY at/after \code{from} where \code{daily_value}
#'       crosses upward through \code{threshold} (i.e. was below on the previous day,
#'       at or above today) and stays above for \code{min_duration} days.
#'       \code{-9999} if no such upward crossing exists.}
#'     \item{doy_cross_down}{First DOY at/after \code{from} where \code{daily_value}
#'       crosses downward through \code{threshold} and stays below for
#'       \code{min_duration} days. \code{-9999} if none.}
#'   }
#' @export

calcDoyCrossThreshold <- function(daily_value, threshold, min_duration = 1L,
                                  from = 1L) {

  n <- length(daily_value)

  is_above      <- daily_value >= threshold
  is_above_prev <- c(is_above[n], is_above[1:(n - 1)])

  cross <- as.integer(is_above) - as.integer(is_above_prev)

  up_cand   <- which(cross ==  1L)
  down_cand <- which(cross == -1L)

  # Sustained-excursion guard: keep only crossings that persist on the new side
  # of the threshold for >= min_duration consecutive days (evaluated circularly).
  if (min_duration > 1L) {
    persists <- function(d, want) {
      idx <- ((d - 1L + 0:(min_duration - 1L)) %% n) + 1L
      all(is_above[idx] == want)
    }
    if (length(up_cand))
      up_cand   <- up_cand[vapply(up_cand,   persists, logical(1), want = TRUE)]
    if (length(down_cand))
      down_cand <- down_cand[vapply(down_cand, persists, logical(1), want = FALSE)]
  }

  # First candidate at/after `from`, scanning circularly (from = 1 -> plain first).
  pick_first <- function(cand) {
    if (!length(cand)) return(-9999L)
    cand[which.min((cand - from) %% n)]
  }
  day_cross_up   <- pick_first(up_cand)
  day_cross_down <- pick_first(down_cand)

  return(list("doy_cross_up"   = day_cross_up,
              "doy_cross_down" = day_cross_down))
}
