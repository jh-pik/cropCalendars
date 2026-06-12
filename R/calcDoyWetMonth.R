#' @title Calculate beginning of the rainy season
#'
#' @description Find the start of the wettest 120 days as the DOY that maximises
#' the 120-day accumulated water balance, computed as the ratio of summed
#' precipitation to summed PET over the window (\eqn{\sum P / \sum PET}). This is
#' evaluated on the climatological daily P and PET cycles (365 values each, one
#' per DOY, averaged across years) and the rolling window wraps circularly across
#' the year boundary.
#'
#' Forming the ratio \emph{after} summing — rather than averaging the daily P/PET
#' ratio — is essential: on the rare day where PET ≈ 0 the daily ratio
#' \eqn{P/PET} explodes (by ~10^6), and a single such day in the 30-year
#' climatology dominates the window and makes the detected onset jump by months.
#' \eqn{\sum PET} over 120 days is never near zero, so the ratio-of-sums is stable.
#'
#' @param daily_prec Numeric vector of length 365. Climatological daily
#'   precipitation (mm), one value per DOY averaged across the reference period
#'   (the \code{dprec} element from \code{calcMonthlyClimate}).
#' @param daily_pet Numeric vector of length 365. Climatological daily PET (mm),
#'   one value per DOY (the \code{dpet} element from \code{calcMonthlyClimate}).
#'
#' @return Integer DOY (1–365) of the start of the 120-day wettest window.
#' @export

calcDoyWetMonth <- function(daily_prec, daily_pet) {
  which.max(.circRollSum(daily_prec, 120) / pmax(.circRollSum(daily_pet, 120), 1e-6))
}

# Circular rolling-window sum: returns, for each start position p (1-based), the
# sum of `w` consecutive values of `x` wrapping across the year boundary. O(n) via
# a cumulative sum (the previous per-window vapply was O(n*w) and dominated the
# crop-calendar runtime). Result is identical to summing each window directly.
.circRollSum <- function(x, w) {
  n  <- length(x)
  cs <- cumsum(c(0, x, x[seq_len(w - 1L)]))   # length n + w
  cs[(1:n) + w] - cs[1:n]
}
