# Daily-climatology analogues of the monthly extremum reductions used by the
# rule-based crop calendars. The original rules reduced the 12 calendar-month
# means (max/min/argmin of monthly_temp, min of monthly_ppet, the month-over-month
# mppet_diff), which quantises temperatures and DOYs to whole-month resolution:
# the warmest/coldest month -- and hence the derived sowing/harvest dates and the
# harvest-rule class -- flips by ~30 days (or switches branch) between adjacent
# climate windows when two months are near-tied.
#
# A daily climatological value (dtemp/dprec/dpet) is already a per-DOY mean over
# the 30-year window, smoothed by cross_smooth_window, so the equivalent statistic
# over a 30-day window reproduces the calendar-month quantity continuously, with no
# month-boundary discretisation. These helpers are the daily replacements; each is
# the strict analogue of the monthly reduction it supersedes. They reuse the
# package-internal .circRollSum (O(n) circular window sum) and .doyWarmestWindow.

# Mean of the warmest `width`-day window (daily analogue of max(monthly_temp)).
.warmestWindowMean <- function(daily_temp, width = 30L) {
  max(.circRollSum(daily_temp, width)) / width
}

# Mean of the coldest `width`-day window (daily analogue of min(monthly_temp)).
.coldestWindowMean <- function(daily_temp, width = 30L) {
  min(.circRollSum(daily_temp, width)) / width
}

# Centre DOY of the coldest `width`-day window (daily analogue of the coldest-month
# mid-day). .doyWarmestWindow on the negated series finds the coldest window.
.doyColdestWindow <- function(daily_temp, width = 30L) {
  .doyWarmestWindow(-daily_temp, width)
}

# Driest `width`-day-window P/PET (daily analogue of min(monthly_ppet)): the minimum
# over DOY of the spike-free 30-day ratio-of-sums (Sum P / Sum PET).
.driestWindowPpet <- function(daily_prec, daily_pet, width = 30L) {
  min(.circRollSum(daily_prec, width) / pmax(.circRollSum(daily_pet, width), 1e-6))
}

# Daily analogue of mppet_diff (= mppet[m] - mppet[m+1], the month-over-month
# moisture trend; > 0 means the next month is drier). At each DOY it is the 30-day
# Sum P / Sum PET centred at that DOY minus the same window centred ~`width` days
# later, so a positive value flags a declining-moisture DOY exactly as the monthly
# version flagged a declining-moisture month -- but on the daily grid, removing the
# residual whole-month quantization of the second wet-season-end candidate (doy_wet2).
.dailyPpetDiff <- function(daily_prec, daily_pet, width = 30L) {
  r <- .circRollSum(daily_prec, width) / pmax(.circRollSum(daily_pet, width), 1e-6)
  n <- length(r); h <- width %/% 2L
  idx <- function(k) ((seq_len(n) - 1L + k) %% n) + 1L
  # .circRollSum at DOY d is the forward window [d, d+width-1] (centre d + h), so
  # r[idx(-h)] is the window centred at d and r[idx(h)] the window centred ~width
  # days later: their difference is mppet_diff valued at d.
  r[idx(-h)] - r[idx(h)]
}
