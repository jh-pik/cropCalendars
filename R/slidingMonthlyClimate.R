#' @title Sliding-window monthly-climate ring buffer
#'
#' @description Maintains a running \code{window}-year climatology that can be
#' advanced one year at a time, so a calendar can be computed for every year from
#' its preceding \code{window} years (a 30-year normal slid by 1 year), reading
#' each raw year only once.
#'
#' The per-year aggregation reuses the validated streaming engine: each pushed
#' year is reduced to a one-year accumulator via \code{addYearMonthlyClimate}, and
#' the ring keeps the field-wise running sum of the last \code{window} of these
#' (subtracting the year that falls out). \code{ringClimatology()} then calls
#' \code{finalizeMonthlyClimate()} on that running sum, so its output is identical
#' to building the same window from scratch (bit-identical before any slide;
#' identical to ~1e-9 after, due to the add/subtract floating-point path — the
#' resulting integer calendar dates are unaffected).
#'
#' @param ncells number of grid cells.
#' @param window averaging window length in years (e.g. 30).
#' @param pet_method PET method, passed to the engine (\code{"fao56"} / \code{"pt"}).
#' @return `initClimateRing` returns a ring; `ringClimatology` returns the same
#'   list of fields as `finalizeMonthlyClimate` (mtemp, mprec, mpet, mppet,
#'   mppet_diff, dtemp, dppet, dprec, dpet); `dppet` is left at zero (unused).
#' @name slidingMonthlyClimate
#' @export
initClimateRing <- function(ncells, window, pet_method = c("fao56", "pt")) {
  pet_method <- match.arg(pet_method)
  list(
    window     = as.integer(window),
    ncells     = ncells,
    pet_method = pet_method,
    slots      = vector("list", window),  # circular buffer of one-year accumulators
    pos        = 0L,                       # index of the most recently written slot
    count      = 0L,                       # years currently in the window
    # Running field-wise sum of the slots, as a zeroed engine accumulator. dppet
    # (D_psum) is dropped (no rule uses it), so it is never accumulated.
    run        = initMonthlyClimate(ncells, pet_method)
  )
}

# Fields summed across the ring (everything finalizeMonthlyClimate reads except
# D_psum, which stays zero).
.ringFields <- c("M_tas", "M_pr", "M_pet", "M_ppet",
                 "D_tsum", "D_prsum", "D_petsum", "D_cnt")

#' @rdname slidingMonthlyClimate
#' @param ring ring from \code{initClimateRing}.
#' @param temp,prec,dates,swdown,lwdown,windspeed,humid,ps,lat one calendar year
#'   of forcings, exactly as passed to \code{addYearMonthlyClimate}.
#' @export
pushClimateYear <- function(ring, temp, prec, dates,
                            swdown = NULL, lwdown = NULL, windspeed = NULL,
                            humid = NULL, ps = 101325, lat = NULL) {
  # One-year contribution, via the validated engine (single-year accumulator).
  slot <- addYearMonthlyClimate(
    initMonthlyClimate(ring$ncells, ring$pet_method),
    temp = temp, prec = prec, dates = dates,
    swdown = swdown, lwdown = lwdown, windspeed = windspeed,
    humid = humid, ps = ps, lat = lat
  )
  slot$D_psum <- NULL  # not retained (dppet dropped) — saves ~1/4 of the memory

  newpos <- (ring$pos %% ring$window) + 1L
  if (ring$count == ring$window) {                 # window full → drop the oldest
    old <- ring$slots[[newpos]]
    for (f in .ringFields) ring$run[[f]] <- ring$run[[f]] - old[[f]]
    ring$run$nyears <- ring$run$nyears - 1L
  } else {
    ring$count <- ring$count + 1L
  }
  for (f in .ringFields) ring$run[[f]] <- ring$run[[f]] + slot[[f]]
  ring$run$nyears <- ring$run$nyears + 1L

  ring$slots[[newpos]] <- slot
  ring$pos <- newpos
  ring
}

#' @rdname slidingMonthlyClimate
#' @export
ringClimatology <- function(ring) {
  if (ring$count == 0L) stop("Empty ring: push at least one year first.")
  finalizeMonthlyClimate(ring$run)
}
