# Daily-climatology analogues of the monthly extremum reductions used by the
# rule-based crop calendars. The original rules reduced the 12 calendar-month
# means (max/min/argmin of monthly_temp, min of monthly_ppet, the month-over-month
# mppet_diff), which quantises temperatures and DOYs to whole-month resolution:
# the warmest/coldest month -- and hence the derived sowing/harvest dates and the
# harvest-rule class -- flips by ~30 days (or switches branch) between adjacent
# climate windows when two months are near-tied.
#
# A daily climatological value (dtemp/dprec/dpet) is already a per-DOY mean over
# the 30-year window, so the equivalent statistic over a `width`-day window
# reproduces the calendar-month quantity continuously, with no month-boundary
# discretisation. These helpers are the daily replacements; each is the strict
# analogue of the monthly reduction it supersedes. They build on .circRoll (below)
# and .doyWarmestWindow.

# Circular rolling-window reduction over a DOY/month-indexed cycle -- the single
# windowing primitive behind every daily reduction AND the rule-time smoothing
# (also used by the 120-day wettest window in calcDoyWetMonth). For each anchor day
# it sums `w` consecutive values of `x`, wrapping across the year boundary; O(n) via
# a cumulative sum (an explicit per-window vapply is O(n*w) and dominated the
# crop-calendar runtime). Identical to summing each window directly. w <= 1 returns
# `x` unchanged (reduction off); callers guarantee w < n.
# The same windows, re-indexed by the chosen anchor day (a circular shift of the
# start-positioned sum):
#   position = "start"  : output[d] = window [d, d+w-1], anchored at its FIRST day
#                         (default; the 120-day wettest window relies on this start
#                         indexing). shift 0.
#   position = "center" : output[d] = window centred on d. shift h = w %/% 2.
#   position = "end"    : output[d] = window [d-w+1, d], anchored at its LAST day.
#                         shift w - 1. (Used with "start" to form .dailyPpetDiff as
#                         a trailing-minus-leading difference about a junction day.)
#   mean = FALSE / TRUE : window sum (default) or window mean (sum / w).
#
# CAVEAT (centred means / argmax): for an EVEN w there is no exact symmetric
# centre, so it is taken as h = w %/% 2 and the window sits half a day low
# ([d-h, d+h-1] rather than a symmetric [d-h, d+h]); an ODD w (e.g. the default
# smooth_window = 31) is exactly symmetric. "start"/"end" are exact for any w (a
# window's first/last day is unambiguous), and the extremum reductions (max/min
# over ALL windows) are invariant to the anchor -- only the centred-mean smoothing
# and the centred argmax carry the half-day shift.
.circRoll <- function(x, w, position = c("start", "center", "end"), mean = FALSE) {
  position <- match.arg(position)
  w <- as.integer(w)
  if (is.na(w) || w <= 1L) return(x)
  n     <- length(x)
  cs    <- cumsum(c(0, x, x[seq_len(w - 1L)]))   # length n + w
  s     <- cs[(1:n) + w] - cs[1:n]               # start-positioned window sum
  shift <- switch(position, start = 0L, center = w %/% 2L, end = w - 1L)
  if (shift > 0L) s <- s[((seq_len(n) - 1L - shift) %% n) + 1L]
  if (mean) s / w else s
}

# Reconstruct the 12 calendar-month values from a 365-day DOY climatology: the
# monthly mean (temperature) or monthly sum (precipitation) over each month's DOYs.
# Accepts a single 365-day vector (-> length-12 vector) OR an [ncells x 365] matrix
# (-> [ncells x 12] matrix, aggregated row-wise) -- the per-pixel path passes a vector,
# the gridded PHU path (generatePHUTserie_isimip3) a matrix. This lets the sliding ring
# carry only the daily climatology and build the monthly seasonality stats downstream.
# Exact for every month except February: the daily climatology folds Feb 29 onto DOY 59,
# so the per-year leap structure is not recoverable and Feb differs by < ~1 % (immaterial
# to the seasonality CV classifier, which also has the seas_eps deadband).
.monthlyFromDaily <- function(daily, agg = c("mean", "sum")) {
  agg    <- match.arg(agg)
  ndays  <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)   # 365-day calendar
  end    <- cumsum(ndays); start <- end - ndays + 1L
  is_mat <- !is.null(dim(daily))
  m      <- if (is_mat) daily else matrix(daily, nrow = 1L)      # [ncells x 365]
  aggf   <- if (agg == "mean") rowMeans else rowSums
  out    <- vapply(seq_len(12L),
                   function(j) aggf(m[, start[j]:end[j], drop = FALSE]),
                   numeric(nrow(m)))
  out    <- matrix(out, nrow = nrow(m), ncol = 12L)
  if (is_mat) out else out[1, ]                                  # vector in -> vector out
}

# Mean of the warmest `width`-day window (daily analogue of max(monthly_temp)).
.warmestWindowMean <- function(daily_temp, width = 30L) {
  max(.circRoll(daily_temp, width, mean = TRUE))
}

# Mean of the coldest `width`-day window (daily analogue of min(monthly_temp)).
.coldestWindowMean <- function(daily_temp, width = 30L) {
  min(.circRoll(daily_temp, width, mean = TRUE))
}

# Centre DOY of the coldest `width`-day window (daily analogue of the coldest-month
# mid-day). .doyWarmestWindow on the negated series finds the coldest window.
.doyColdestWindow <- function(daily_temp, width = 30L) {
  .doyWarmestWindow(-daily_temp, width)
}

# Daily analogue of min(monthly_ppet) -- the always-wet aridity floor -- is just
# min(daily_ppet) at the call site (calcHarvestDateVector): under the unified
# smooth_window that minimum is anchor-invariant and bit-identical to the driest
# smooth_window-day window's Sum P / Sum PET, so it needs no separate helper and is
# the SAME series the wet-end crossing reads. (Former .driestWindowPpet, removed.)

# Daily analogue of mppet_diff (= mppet[m] - mppet[m+1], the month-over-month
# moisture trend; > 0 means the next month is drier). Built exactly as the monthly
# version: at each junction day j it is the `width`-day Sum P / Sum PET of the window
# ENDING at j (the trailing/"this month" side) minus the window STARTING at j (the
# leading/"next month" side), so a positive value flags declining moisture across j --
# on the daily grid, removing the residual whole-month quantization. The junction value
# is assigned back to the DOY d = j - h that the two windows straddle (h = width %/% 2).
# NB: the trend wet-end candidate (doy_wet2) this fed was RETIRED (see NEWS); this helper
# is retained only for the diagnostic replay scripts under utils/ggcmi_ph3.
# NORMALISED to a per-30-day rate (x 30/lag): the two window centres are `lag = 2*h`
# days apart, but the comparison threshold ppet_ratio_diff is calibrated as a
# delta(P/PET) per 30 days (= the month-over-month mppet_diff), so the diff is rescaled
# to that horizon and the threshold stays valid for ANY `width`. At width 30 and 31
# lag=30 -> factor 1 (identical to the original).
.dailyPpetDiff <- function(daily_prec, daily_pet, width = 30L) {
  ppet <- function(pos) .circRoll(daily_prec, width, pos) /
                        pmax(.circRoll(daily_pet, width, pos), 1e-6)
  jd  <- ppet("end") - ppet("start")   # (window ending at j) - (window starting at j)
  n <- length(jd); h <- width %/% 2L; lag <- 2L * h
  idx <- ((seq_len(n) - 1L + h) %% n) + 1L   # junction j = d + h, assigned back to DOY d
  jd[idx] * (30 / lag)
}
