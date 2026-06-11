#' @title Streaming, cell-vectorised monthly-climate accumulator
#'
#' @description Building blocks for computing the average monthly climate (and the
#' daily DOY climatologies) used by the rule-based crop calendars, without holding
#' the full multi-year daily series in memory at once. Years are added one at a
#' time; each call is vectorised over an arbitrary number of grid cells.
#'
#' `calcMonthlyClimate()` is a thin per-pixel wrapper around these; the gridded
#' pipeline (utils/ggcmi_ph3/01a) feeds the full grid one year at a time. Both share
#' this single implementation so the aggregation, the PET method switch, the P/PET
#' flooring and the leap-year DOY handling live in exactly one place.
#'
#' Aggregation (reproduces the previous `tapply`-based code):
#' \itemize{
#'   \item monthly fields = mean over years of the per-year monthly statistic
#'         (mean temperature; summed precipitation, PET; summed P/PET ratio);
#'   \item daily fields (`dtemp`, `dppet`) = per-DOY sum / per-DOY count over all
#'         days of all years, i.e. `tapply(., DOY, mean)` (the count is not uniform:
#'         `date_to_doy(skip_feb29 = TRUE)` folds Feb 29 onto DOY 59).
#' }
#'
#' @param ncells number of grid cells the accumulator covers.
#' @param pet_method PET method, \code{"pt"} (Priestley-Taylor, calcPET) or
#'   \code{"fao56"} (FAO-56 Penman-Monteith, calcPET_FAO56).
#' @return `initMonthlyClimate` returns an accumulator (a list of zeroed matrices);
#'   `finalizeMonthlyClimate` returns the monthly climate as a list of seven fields
#'   (vectors when `ncells == 1`, otherwise `[ncells x 12]` / `[ncells x 365]`
#'   matrices), matching `calcMonthlyClimate()`'s return.
#' @name monthlyClimateAccum
#' @export
initMonthlyClimate <- function(ncells, pet_method = c("fao56", "pt")) {
  pet_method <- match.arg(pet_method)
  list(
    M_tas  = matrix(0, ncells, 12),  M_pr   = matrix(0, ncells, 12),
    M_pet  = matrix(0, ncells, 12),  M_ppet = matrix(0, ncells, 12),
    D_tsum = matrix(0, ncells, 365), D_psum = matrix(0, ncells, 365),
    D_cnt  = numeric(365), nyears = 0L,
    pet_method = pet_method, ncells = ncells
  )
}

# Daily PET for one year, vectorised over cells. temp/swdown/... are [ncells x ndays].
# fao56 and pt-with-radiation are elementwise; pt without radiation (orbital fallback)
# is scalar per day and only supported for a single cell (the per-pixel path).
.petDaily <- function(temp, swdown, lwdown, windspeed, humid, ps, lat, day, method) {
  if (method == "fao56") {
    return(calcPET_FAO56(temp, windspeed, humid, swdown, lwdown, ps))
  }
  if (!is.null(swdown) && !is.null(lwdown)) {
    return(calcPET(temp, swdown = swdown, lwdown = lwdown)) # observed Rn: no lat/day needed
  }
  if (nrow(temp) > 1L) {
    stop("pet_method = 'pt' without swdown/lwdown is not cell-vectorised; supply radiation.")
  }
  matrix(mapply(calcPET, temp = temp[1, ], lat = lat, day = day), nrow = 1L)
}

#' @rdname monthlyClimateAccum
#' @param acc accumulator from \code{initMonthlyClimate}.
#' @param temp,prec daily mean temperature (deg C) and precipitation (mm) for one
#'   calendar year: a vector (single cell) or an \code{[ncells x ndays]} matrix.
#' @param dates character vector of the year's dates ("YYYY-MM-DD"), length ndays.
#' @param swdown,lwdown,windspeed,humid,ps daily forcings matching \code{temp}'s
#'   shape (see \code{calcMonthlyClimate}); \code{ps} may be scalar.
#' @param lat latitude(s); only used by the Priestley-Taylor orbital fallback.
#' @export
addYearMonthlyClimate <- function(acc, temp, prec, dates,
                                  swdown = NULL, lwdown = NULL, windspeed = NULL,
                                  humid = NULL, ps = 101325, lat = NULL) {
  as_mat <- function(x) if (is.null(x) || !is.null(dim(x))) x else matrix(x, nrow = 1L)
  temp <- as_mat(temp); prec <- as_mat(prec)
  swdown <- as_mat(swdown); lwdown <- as_mat(lwdown)
  windspeed <- as_mat(windspeed); humid <- as_mat(humid); ps <- as_mat(ps)

  mon <- date_to_month(dates)
  doy <- date_to_doy(dates, skip_feb29 = TRUE)

  pet  <- .petDaily(temp, swdown, lwdown, windspeed, humid, ps, lat, doy, acc$pet_method)
  ppet <- prec / pmax(pet, 1e-6)

  for (m in 1:12) {
    dom     <- which(mon == m)
    pr_mon  <- rowSums(prec[, dom, drop = FALSE])
    pet_mon <- rowSums(pet[,  dom, drop = FALSE])
    acc$M_tas[,  m] <- acc$M_tas[,  m] + rowMeans(temp[, dom, drop = FALSE])
    acc$M_pr[,   m] <- acc$M_pr[,   m] + pr_mon
    acc$M_pet[,  m] <- acc$M_pet[,  m] + pet_mon
    acc$M_ppet[, m] <- acc$M_ppet[, m] + pr_mon / pmax(pet_mon, 1e-6)
  }

  # A DOY can recur within a leap year (Feb 29 -> 59), so accumulate sum and count.
  for (k in seq_len(ncol(temp))) {
    d <- doy[k]
    acc$D_tsum[, d] <- acc$D_tsum[, d] + temp[, k]
    acc$D_psum[, d] <- acc$D_psum[, d] + ppet[, k]
  }
  acc$D_cnt  <- acc$D_cnt + tabulate(doy, nbins = 365)
  acc$nyears <- acc$nyears + 1L
  acc
}

#' @rdname monthlyClimateAccum
#' @export
finalizeMonthlyClimate <- function(acc) {
  if (acc$nyears == 0L) stop("No years accumulated.")
  mtemp      <- round(acc$M_tas  / acc$nyears, 5)
  mprec      <- round(acc$M_pr   / acc$nyears, 5)
  mpet       <- round(acc$M_pet  / acc$nyears, 5)
  mppet      <- round(acc$M_ppet / acc$nyears, 5)
  mppet_diff <- mppet - mppet[, c(2:12, 1), drop = FALSE]
  dtemp      <- sweep(acc$D_tsum, 2, acc$D_cnt, "/")
  dppet      <- sweep(acc$D_psum, 2, acc$D_cnt, "/")

  out <- list(mtemp = mtemp, mprec = mprec, mpet = mpet, mppet = mppet,
              mppet_diff = mppet_diff, dtemp = dtemp, dppet = dppet)
  if (acc$ncells == 1L) out <- lapply(out, as.vector) # match per-pixel contract
  out
}
