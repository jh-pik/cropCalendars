# Equivalence harness for the scalar->vectorised consolidation (PR-A, and the basis for
# PR-C). Each public scalar function (calcVd, isWinterCrop, ...) now delegates to its
# vectorised core; these tests embed the ORIGINAL scalar body verbatim as ground truth and
# fuzz a wide input range to assert the wrapper reproduces it bit-for-bit at the production
# defaults. When PR-C adds calcVrf / calcPHU wrappers, extend this file the same way.

# ---- ground-truth reference implementations (the pre-refactor scalar bodies) ----

ref_calcVd <- function(temp_mean_month, max.vern.days = 70, max.vern.months = 5,
                       tv2 = 3, tv3 = 10) {
  max.vern.days.month <- max.vern.days / max.vern.months
  temp <- array(0, c(12, 2)); temp[, 1] <- temp_mean_month; temp[, 2] <- 1:12
  sort_temp      <- temp[order(temp[, 1]), ]
  coldest_months <- sort_temp[1:5, 2]            # original always takes the 5 coldest
  days_required  <- 0
  for (m in 1:max.vern.months) {
    t <- temp_mean_month[coldest_months[m]]
    if      (t <= tv2) days <- max.vern.days.month
    else if (t >= tv3) days <- 0
    else               days <- max.vern.days.month * (1 - (t - tv2) / (tv3 - tv2))
    days_required <- days_required + days
  }
  round(days_required)
}

ref_isWinterCrop <- function(start, end, tcm, lat) {
  growp <- ifelse(start <= end, end - start, 365 + end - start)
  wc <- 0
  if (!is.na(start) && start > 0 && !is.na(lat) && !is.na(tcm)) {
    if (lat > 0) {
      if (((start + growp > 365) && (growp >= 150)) && (tcm >= -10 && tcm <= 7)) wc <- 1
    } else {
      if (((start < 182) && (start + growp > 182) && (growp >= 150)) &&
          (tcm >= -10 && tcm <= 7)) wc <- 1
    }
  }
  wc
}

# ---- calcVd ----

test_that("calcVd wrapper reproduces the original scalar (production max.vern.months = 5)", {
  set.seed(20260617)
  for (i in 1:2000) {
    tmm <- runif(12, -25, 35)
    mvd <- sample(40:90, 1)
    tv2 <- sample(0:5, 1); tv3 <- tv2 + sample(3:10, 1)
    expect_equal(
      cropCalendars::calcVd(tmm, max.vern.days = mvd, max.vern.months = 5,
                            tv2 = tv2, tv3 = tv3),
      ref_calcVd(tmm, max.vern.days = mvd, max.vern.months = 5, tv2 = tv2, tv3 = tv3)
    )
  }
})

test_that("calcVd wrapper matches the original for max.vern.months in {1,3,4,5}", {
  # The original picks the 5 coldest months but loops max.vern.months, so it is only
  # self-consistent for max.vern.months <= 5; the vec keeps exactly max.vern.months.
  set.seed(11)
  for (mvm in c(1L, 3L, 4L, 5L)) for (i in 1:300) {
    tmm <- runif(12, -25, 35)
    expect_equal(
      cropCalendars::calcVd(tmm, max.vern.days = 70, max.vern.months = mvm, tv2 = 3, tv3 = 10),
      ref_calcVd(tmm, max.vern.days = 70, max.vern.months = mvm, tv2 = 3, tv3 = 10)
    )
  }
})

test_that("calcVd boundary temperatures (exactly tv2 / tv3) match", {
  tmm <- c(rep(3, 5), rep(20, 7))           # five months exactly at tv2
  expect_equal(cropCalendars::calcVd(tmm, 70, 5, 3, 10), ref_calcVd(tmm, 70, 5, 3, 10))
  tmm <- c(rep(10, 5), rep(20, 7))          # five months exactly at tv3
  expect_equal(cropCalendars::calcVd(tmm, 70, 5, 3, 10), ref_calcVd(tmm, 70, 5, 3, 10))
})

# ---- isWinterCrop ----

test_that("isWinterCrop wrapper reproduces the original scalar over valid inputs", {
  set.seed(99)
  for (i in 1:4000) {
    start <- sample(1:365, 1); end <- sample(1:365, 1)
    tcm   <- runif(1, -25, 25); lat <- runif(1, -60, 60)
    expect_equal(
      as.numeric(cropCalendars::isWinterCrop(start, end, tcm, lat)),
      ref_isWinterCrop(start, end, tcm, lat)
    )
  }
})

test_that("isWinterCrop wrapper matches near the rule boundaries", {
  cases <- expand.grid(
    start = c(1, 181, 182, 183, 215, 366),   # SH mid-winter pivot at 182
    growp = c(149, 150, 151),                # length threshold
    tcm   = c(-11, -10, 7, 8),               # cold-but-not-killing band edges
    lat   = c(-1, 1))                        # hemisphere
  for (r in seq_len(nrow(cases))) {
    s <- cases$start[r]; e <- ((s + cases$growp[r] - 1) %% 365) + 1
    expect_equal(
      as.numeric(cropCalendars::isWinterCrop(s, e, cases$tcm[r], cases$lat[r])),
      ref_isWinterCrop(s, e, cases$tcm[r], cases$lat[r]),
      info = sprintf("start=%d growp=%d tcm=%g lat=%g", s, cases$growp[r], cases$tcm[r], cases$lat[r])
    )
  }
})

test_that("isWinterCrop guard cases return 0 (NA/0 start, NA lat, NA tcm)", {
  expect_equal(as.numeric(cropCalendars::isWinterCrop(NA, 200, 0, 45)), 0)
  expect_equal(as.numeric(cropCalendars::isWinterCrop(0,  200, 0, 45)), 0)
  expect_equal(as.numeric(cropCalendars::isWinterCrop(100, 300, NA, 45)), 0)
  expect_equal(as.numeric(cropCalendars::isWinterCrop(100, 300, 0, NA)), 0)
})

test_that("isWinterCrop is graceful on NA end (returns NA, the old scalar errored)", {
  expect_true(is.na(cropCalendars::isWinterCrop(100, NA, 0, 45)))
})
