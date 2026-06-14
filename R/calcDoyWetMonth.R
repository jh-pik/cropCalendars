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
#' @param prev_doy Integer DOY (1–365) of the wettest-window start chosen for the
#'   \emph{previous} year, or \code{NA} (the default) for no prior state. Used
#'   only when \code{eps > 0} to apply temporal hysteresis (see \code{eps}).
#' @param eps Non-negative numeric (default 0 = off). Hysteresis strength for the
#'   sliding window. In a sliding-window pipeline the 30-yr climatology shifts only
#'   slightly year to year, but in semi-arid / monsoon-fringe cells several 120-day
#'   windows can be near-tied (the \eqn{\sum P/\sum PET} curve has multiple peaks
#'   within a few percent), so the plain argmax flips between far-apart peaks
#'   between adjacent years (spurious oscillation in sowing, and hence harvest).
#'
#'   With \code{eps > 0} and a \code{prev_doy} supplied, the window is chosen by a
#'   \strong{distance-weighted, max-normalised} score rather than the raw argmax:
#'   \deqn{score(d) = \frac{w(d)}{\max_k w(k)} \,\cdot\, \Big[(1-\varepsilon) +
#'         \varepsilon \, e^{-(x/decay)^2}\Big], \quad x = \frac{\Delta(d, prev)}{365/2},}
#'   where \eqn{w} is \eqn{\sum P/\sum PET} and \eqn{\Delta(d, prev)} is the circular
#'   DOY distance to last year's chosen window. The first factor is the window's
#'   goodness as a fraction of the year's best (so \code{eps} is dimensionless and
#'   comparable across cells); the second is a \strong{Gaussian} down-weight by
#'   normalised distance from last year's pick. It equals 1 at \eqn{\Delta = 0}
#'   (last year's window, plus a near-flat plateau for small drifts), decays on a
#'   scale set by \code{decay}, and \strong{floors at \eqn{1-\varepsilon}} (never 0)
#'   so a decisively better far peak still wins — a genuine regime shift is followed,
#'   a near-tie flip is not. Compared with a linear weight the Gaussian declines
#'   faster through the mid-distances where chronic bimodal flippers live, so it is
#'   stickier there for the same floor.
#'
#'   \code{eps} sets the floor (max penalty: far weight = \eqn{1-\varepsilon};
#'   \code{eps = 0.5} -> floor 0.5, i.e. a far peak must be >2x wetter to win).
#'   \code{decay} sets how fast the weight falls (smaller = faster/stickier; ~0.3
#'   reaches the floor by ~3 months out).
#'
#'   With \code{eps = 0} or \code{prev_doy = NA} the plain argmax is returned
#'   (backward compatible).
#'
#' @param decay Positive numeric (default 0.3). Gaussian decay scale for the
#'   distance weight, in units of half a year (so \code{decay = 0.3} ~= 55 days).
#'   Smaller is stickier. Only used when \code{eps > 0}.
#'
#' @return Integer DOY (1–365) of the start of the 120-day wettest window.
#' @export

calcDoyWetMonth <- function(daily_prec, daily_pet,
                            prev_doy = NA_integer_, eps = 0, decay = 0.3) {
  ws <- .circRollSum(daily_prec, 120) / pmax(.circRollSum(daily_pet, 120), 1e-6)
  if (eps <= 0 || is.na(prev_doy) || prev_doy < 1L || prev_doy > length(ws))
    return(as.integer(which.max(ws)))
  # Distance-weighted, max-normalised selection. q in (0,1] is each window's
  # goodness relative to the year's best; the Gaussian weight down-weights windows
  # by circular distance to last year's pick (faster mid-distance decline than a
  # linear ramp) but floors at 1-eps so a decisively better far peak still wins
  # (regime shift). Near DOYs (weight ~1) let the same peak drift freely.
  n     <- length(ws)
  d     <- abs(seq_len(n) - prev_doy); d <- pmin(d, n - d)   # circular distance (days)
  x     <- d / (n / 2)                                       # normalised to [0, 1]
  q     <- ws / max(ws)
  score <- q * ((1 - eps) + eps * exp(-(x / decay)^2))
  as.integer(which.max(score))
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
