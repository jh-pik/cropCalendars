#' @title Generate an annual crop Phenological Heat Unit (PHU) time series
#' for the LPJmL model.
#'
#' @export
generatePHUTserie_isimip3 <- function(
    ncdir         = NULL,
    gcm           = NULL,
    scen          = NULL,
    cro           = NULL,
    irri          = NULL,
    SYs           = NULL,
    EYs           = NULL,
    FYnc          = NULL,
    LYnc          = NULL,
    grid_df       = NULL,
    crop_par_file = NULL,
    ncfile        = NULL,
    smooth_window = 1L

) {
  # grid_df (LPJmL grid: data.frame with lon/lat) is now an explicit argument rather
  # than a global read from the caller's environment. years/nyears are derived from
  # FYnc:LYnc here for the same reason (they were the only remaining hidden globals).
  if (is.null(grid_df)) stop("generatePHUTserie_isimip3: grid_df (LPJmL grid) is required.")
  years  <- FYnc:LYnc
  nyears <- length(years)
  # smooth_window (years, odd, default 1): widens ONLY the temperature averaging
  # used for the heat-unit sum, centred on each period, clamped to [FYnc, LYnc].
  # With 1 (default) the PHU matches the growing period of that period exactly (the
  # original behaviour); larger values damp single-year weather noise in the PHU of
  # the annual product without changing the (already smooth) sowing/harvest dates.
  half <- (as.integer(smooth_window) - 1L) %/% 2L

  cr <- which(crop_ls[["ggcmi"]] == cro)
  ir <- which(irri_ls[["ggcmi"]] == irri)

  # Read the DRS crop-calendar file directly if its path is supplied (annual
  # pipeline), otherwise reconstruct the legacy intermediate name under ncdir.
  ncfname <- if (!is.null(ncfile)) ncfile else paste0(ncdir,
                    crop_ls[["ggcmi"]][cr], "_", irri_ls[["ggcmi"]][ir],
                    "_", gcm, "_", scen, "_", FYnc, "-", LYnc,
                    "_ggcmi_ph3_rule_based_crop_calendar.nc4")
  cat("\nreading:", ncfname)


  # ------------------------------------------------------#
  # Get Crop Parameters ----
  if (is.null(crop_par_file)) {
    crop_par_file <- system.file("extdata", "lpjml_crop_parameters_ggcmi_ph3.csv",
                                package = "cropCalendars", mustWork = TRUE)
  }
  croppar    <- subset(read.csv(crop_par_file, header = T, stringsAsFactors = F),
                       crop == cro)

  basetemp      <- croppar$basetemp
  max.vern.days <- croppar$max.vern.days
  tv1           <- croppar$vern.temp.min
  tv2           <- croppar$vern.temp.opt.min
  tv3           <- croppar$vern.temp.opt.max
  tv4           <- croppar$vern.temp.max

  is_vernal_all  <- crop_ls[["vernal"]][cr] == "yes_all"
  is_vernal_cond <- crop_ls[["vernal"]][cr] == "yes"

  NCELLS <- nrow(grid_df)

  # ------------------------------------------------------#
  # Pre-compute 720x360 <-> NCELLS index mapping (once) ----
  nc_tmp  <- nc_open(ncfname)
  lons    <- ncvar_get(nc_tmp, "lon")
  lats    <- ncvar_get(nc_tmp, "lat")
  pmask   <- !is.na(ncvar_get(nc_tmp, "planting_day", start = c(1, 1, 1),
                              count = c(720, 360, 1)))   # product land mask
  nc_close(nc_tmp)

  # Flat index into [lon x lat] = 720 x 360 array stored column-major in R:
  # matrix(sdate, nrow=720*360) row = (ilat-1)*720 + ilon
  lin_idx <- (match(grid_df$lat, lats) - 1L) * 720L + match(grid_df$lon, lons)

  # Reconcile the LPJmL vs GGCMI grid mismatch (a 1-cell sub-antarctic-island
  # artefact: GGCMI has 178.75,-49.25 while LPJmL has 178.75,-49.75). Any LPJmL cell
  # that lands on a cell absent from the product grid is reassigned the nearest
  # product land cell -- here its adjacent twin -- so it inherits that calendar.
  off <- which(!pmask[lin_idx])
  if (length(off) > 0) {
    land     <- which(as.vector(pmask))               # flat (column-major) land indices
    land_lon <- lons[((land - 1L) %% 720L) + 1L]
    land_lat <- lats[((land - 1L) %/% 720L) + 1L]
    for (i in off) {
      d <- (land_lon - grid_df$lon[i])^2 + (land_lat - grid_df$lat[i])^2
      lin_idx[i] <- land[which.min(d)]
    }
    cat("\nRemapped", length(off), "LPJmL cell(s) off the product grid to nearest land cell")
  }

  # ------------------------------------------------------#
  # Output arrays ----
  phu.annual <- array(NA, c(720L, 360L))
  phu.cube   <- array(NA, c(720L, 360L, nyears))

  # ------------------------------------------------------#
  # Loop through time periods ----
  for (yy in seq_len(length(SYs))) {

    cat("\n--- Doing yy", yy, "---")

    # --------------------------------------------------#
    # Climate: accumulate year-by-year to cap peak RAM ----
    # Temperature window = the period, widened by `half` years each side for
    # smoothing (clamped to the product range). half = 0 -> exactly the period.
    w_lo <- max(FYnc, SYs[yy] - half); w_hi <- min(LYnc, EYs[yy] + half)
    tas_sum <- matrix(0.0, NCELLS, 365L)
    for (yr in w_lo:w_hi) {
      tas_yr  <- get.isimip.tas(gcm, scen, yr, yr, ncells = NCELLS)
      tas_sum <- tas_sum + tas_yr[, , 1L]
      rm(tas_yr)
    }
    tas_mean_day <- tas_sum / length(w_lo:w_hi)
    rm(tas_sum)

    # --------------------------------------------------#
    # Sdate / Hdate: vectorized extraction via flat index ----
    nc    <- nc_open(ncfname)
    sy    <- which(years == SYs[yy])
    nyp   <- length(SYs[yy]:EYs[yy])
    sdate <- ncvar_get(nc, "planting_day", start = c(1, 1, sy), count = c(720, 360, nyp))
    hdate <- ncvar_get(nc, "maturity_day", start = c(1, 1, sy), count = c(720, 360, nyp))
    nc_close(nc)

    sdate_cells <- matrix(sdate, nrow = 720L * 360L)[lin_idx, , drop = FALSE]
    hdate_cells <- matrix(hdate, nrow = 720L * 360L)[lin_idx, , drop = FALSE]
    sdate_avg   <- as.integer(round(rowMeans(sdate_cells, na.rm = TRUE)))
    hdate_avg   <- as.integer(round(rowMeans(hdate_cells, na.rm = TRUE)))
    sdate_avg[is.na(sdate_avg)] <- 1L   # ri2 fill: 1-day growing period
    hdate_avg[is.na(hdate_avg)] <- 2L
    rm(sdate, hdate, sdate_cells, hdate_cells)

    # --------------------------------------------------#
    # Monthly temps (NCELLS x 12) ----
    mtemp_mat <- .monthly_temps_vec(tas_mean_day)

    # --------------------------------------------------#
    # PHU computation (vectorized across all cells) ----

    if (is_vernal_all) {
      # wwh: every cell uses vernal-thermal model
      vd_vec  <- .calc_vd_vec(mtemp_mat, max.vern.days, tv2 = tv2, tv3 = tv3)
      vrf_mat <- .build_vrf_mat(sdate_avg, hdate_avg, tas_mean_day, vd_vec,
                                tv1 = tv1, tv2 = tv2, tv3 = tv3, tv4 = tv4)
      phu_vec <- .calc_phu_vernal_vec(sdate_avg, hdate_avg, tas_mean_day, vrf_mat, basetemp)

    } else if (is_vernal_cond) {
      # rap: vernal-thermal only for winter-type cells
      tcm_vec   <- apply(mtemp_mat, 1L, min)
      wcrop_vec <- .wintercrop_vec(sdate_avg, hdate_avg, tcm_vec, grid_df$lat)

      vd_vec     <- .calc_vd_vec(mtemp_mat, max.vern.days, tv2 = tv2, tv3 = tv3)
      vd_for_vrf <- vd_vec * wcrop_vec   # zero out spring-type cells
      vrf_mat    <- .build_vrf_mat(sdate_avg, hdate_avg, tas_mean_day, vd_for_vrf,
                                   tv1 = tv1, tv2 = tv2, tv3 = tv3, tv4 = tv4)

      phu_thermal <- .calc_phu_thermal_vec(sdate_avg, hdate_avg, tas_mean_day, basetemp)
      phu_vernal  <- .calc_phu_vernal_vec( sdate_avg, hdate_avg, tas_mean_day, vrf_mat, basetemp)
      phu_vec     <- ifelse(wcrop_vec == 1L, phu_vernal, phu_thermal)

    } else {
      # All other crops: pure thermal model
      phu_vec <- .calc_phu_thermal_vec(sdate_avg, hdate_avg, tas_mean_day, basetemp)
    }

    # --------------------------------------------------#
    # Store result in 720x360 annual array via flat index ----
    phu_flat          <- rep(NA_integer_, 720L * 360L)
    phu_flat[lin_idx] <- phu_vec
    phu.annual[]      <- phu_flat

    # Repeat same phu for all years in this time slice
    ys <- which(years %in% SYs[yy]:EYs[yy])
    for (j in ys) phu.cube[, , j] <- phu.annual

    cat(" done\n")

  } # yy


  # ------------------------------------------------------#
  # Save intermediate result ----  (under the output dir, not a global work_dir)
  tmp_dir <- paste0(ncdir, "tmp/")
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)

  fn <- paste0(tmp_dir,
               crop_ls[["ggcmi"]][cr], "_",
               irri_ls[["ggcmi"]][ir], "_",
               gcm, "_", scen, "_", FYnc, "-", LYnc,
               "_ggcmi_ph3_rule_based_phu.Rdata")
  save(phu.cube, file = fn)

  # ------------------------------------------------------#
  # Write Crop-specific Output File ----
  # ------------------------------------------------------#

  ncfname <- paste0(ncdir,
                    crop_ls[["ggcmi"]][cr], "_", irri_ls[["ggcmi"]][ir],
                    "_", gcm, "_", scen, "_", FYnc, "-", LYnc,
                    "_ggcmi_ph3_rule_based_phu.nc4")
  cat("\nwriting:", ncfname)

  # Define dimensions
  londim       <- ncdim_def("lon", "degrees_east",  lons)
  latdim       <- ncdim_def("lat", "degrees_north", lats)
  timdim       <- ncdim_def("time", "year",        years)
  nc_dimension <- list(londim, latdim, timdim)

  # Define variables
  phu_def  <- ncvar_def(name = "phu", units = "degree days", dim = nc_dimension,
                        longname = "Phenological Heat Unit Requirements",
                        prec = "single", compression = 6)

  # Create netCDF file and put arrays
  ncout <- nc_create(ncfname, list(phu_def), verbose = F)

  # Put variables
  ncvar_put(ncout, phu_def, phu.cube)

  # Put additional attributes into dimension and data variables
  ncatt_put(ncout,  "lon", "axis", "X")
  ncatt_put(ncout,  "lat", "axis", "Y")
  ncatt_put(ncout, "time", "axis", "T")

  ncatt_put(ncout, 0, "Crop",
            paste0(crop_ls[["ggcmi"]][cr], "_", irri_ls[["ggcmi"]][ir]))
  ncatt_put(ncout, 0, "Institution",
            "Potsdam Institute for Climate Impact Research (PIK), Germany")
  history <- paste("Created by Sara Minoli on", date(), sep = " ")
  ncatt_put(ncout, 0, "history", history)

  # Close the file, writing data to disk
  nc_close(ncout)

  unlink(fn)

}


# ------------------------------------ #
# Unexported helpers: climate I/O

read.climate.input <- function(fname, ncells, ryear,
                               fyear, lyear, header = 43L,
                               nbands = 365L, dtype = "integer", scalar = 0.1) {

  cat(paste("\nReading climate input:\n-----------------------\n", fname))

  nyears <- lyear - fyear + 1
  dsize  <- ifelse(dtype == "integer", 2, 4)

  cat(paste("\nNBANDS =", sprintf("%.10f", ((file.size(fname) - header) / ncells / dsize / (lyear - ryear + 1)))))

  fcon <- file(fname, "rb")
  x    <- array(data = NA, dim = c(ncells, nbands, nyears))

  for (i in 1:nyears) {
    year <- fyear + (i - 1)
    seek(fcon, where = header + ((year - ryear) * ncells * nbands * dsize), origin = "start")
    for (j in 1:ncells) {
      x[j, , i] <- readBin(fcon, dtype, n = nbands, size = dsize) * scalar
    }
  }
  close.connection(fcon)
  cat(str(x))
  return(x)
}

paste.isimip3b.clm.fn <- function(path, gcm, scen, var, syear, eyear) {
  paste0(path,
         paste(scen, gcm,
               paste(var, tolower(gcm), scen,
                     paste(syear, eyear, sep = "-"), sep = "_"), sep = "/"), ".clm")
}

get.isimip.tas <- function(GCM, SC, SY, EY, ncells) {

  # ISIMIP3a observational forcings (GSWP3-W5E5 etc.): read the single dataset file
  # directly, with NO historical/SSP splice. ryear (= the seek-offset reference) is
  # the file's first year, so FY1 must be the dataset's actual first year. (The
  # ESM branch below still maps to the 1850-2014 / 2015-2100 ISIMIP3b windows.)
  obs_ranges <- list(obsclim     = c(1901, 2019),
                     spinclim    = c(1801, 1900),
                     counterclim = c(1901, 2019))

  if (SC %in% names(obs_ranges)) {

    FY1 <- obs_ranges[[SC]][1]
    LY1 <- obs_ranges[[SC]][2]

  } else {

    if (EY <= 2014) {
      FY1 <- 1850
      LY1 <- 2014
      SC  <- "historical"
    } else if (SY <= 2014 & EY > 2014) {
      FY1 <- 1850
      LY1 <- 2014
      SC1 <- "historical"
      FY2 <- 2015
      LY2 <- 2100
      SC2 <- SC
    } else {
      FY1 <- 2015
      LY1 <- 2100
    }

  }

  if (!(SC %in% names(obs_ranges)) && SY <= 2014 & EY > 2014) {

    tas_fn1 <- paste.isimip3b.clm.fn(isimip3b.path, GCM, SC1, "tas", FY1, LY1)
    tas_fn2 <- paste.isimip3b.clm.fn(isimip3b.path, GCM, SC2, "tas", FY2, LY2)

    tas1 <- read.climate.input(tas_fn1, ncells = ncells, ryear = FY1,
                               fyear = SY, lyear = 2014, header = 43,
                               nbands = 365, dtype = "integer", scalar = 0.1)
    tas2 <- read.climate.input(tas_fn2, ncells = ncells, ryear = FY2,
                               fyear = 2015, lyear = EY, header = 43,
                               nbands = 365, dtype = "integer", scalar = 0.1)

    tas <- array(NA, dim = c(ncells, 365, length(SY:EY)))
    tas[1:ncells, 1:365, 1:dim(tas1)[3]] <- tas1
    tas[1:ncells, 1:365, (dim(tas1)[3] + 1):(dim(tas)[3])] <- tas2
    rm(tas1, tas2)

  } else {

    tas_fn <- paste.isimip3b.clm.fn(isimip3b.path, GCM, SC, "tas", FY1, LY1)

    tas <- read.climate.input(tas_fn, ncells = ncells, ryear = FY1,
                              fyear = SY, lyear = EY, header = 43,
                              nbands = 365, dtype = "integer", scalar = 0.1)
  }

  return(tas)

}


# ------------------------------------ #
# Unexported helpers: vectorized PHU computation

# Monthly mean temperatures: NCELLS x 365 matrix -> NCELLS x 12 matrix.
# Fixed (non-leap) calendar-month day boundaries sday..eday; each output column is the
# row-wise mean of that month's day block. Vectorised twin of .monthlyFromDaily(.,"mean").
.monthly_temps_vec <- function(temp_mat) {
  sday <- c(  1, 32, 60,  91, 121, 152, 182, 213, 244, 274, 305, 335)
  eday <- c( 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
  m_mat <- matrix(0.0, nrow(temp_mat), 12L)
  for (m in 1:12) m_mat[, m] <- rowMeans(temp_mat[, sday[m]:eday[m]])
  m_mat
}

# Vernalization days required (per cell): vectorised twin of calcVd, the canonical
# reference. Per row: rank the 12 monthly means coldest-first (order), keep the
# `max.vern.months` coldest, and give each a vernalization-day credit that is piecewise
# in its temperature -- the full max_per_m below tv2 (vern optimum floor), 0 above tv3
# (too warm to vernalize), a linear ramp down between -- then sum the credits and round.
# (Unlike the scalar, which always picks the 5 coldest regardless of max.vern.months, this
# keeps exactly max.vern.months; identical at the default 5 used in production.)
.calc_vd_vec <- function(mtemp_mat, max.vern.days, max.vern.months = 5L, tv2, tv3) {
  nc         <- nrow(mtemp_mat)
  max_per_m  <- max.vern.days / max.vern.months
  sorted_idx <- t(apply(mtemp_mat, 1, order))
  cold_cols  <- as.vector(sorted_idx[, 1:max.vern.months])
  row_idx    <- rep(seq_len(nc), max.vern.months)
  cold_temp  <- matrix(mtemp_mat[cbind(row_idx, cold_cols)], nrow = nc, ncol = max.vern.months)
  days_m     <- ifelse(cold_temp <= tv2, max_per_m,
                ifelse(cold_temp >= tv3, 0,
                       max_per_m * (1 - (cold_temp - tv2) / (tv3 - tv2))))
  round(rowSums(days_m))
}

# Winter-crop classification (per cell): vectorised twin of isWinterCrop (Portmann 2010,
# with tcm <= 7 not 6). growp is the circular sowing->harvest length (wrapping the year).
# A season is winter crop (1) iff it is long (growp >= 150 d), its coldest month is
# cold-but-not-killing (tcm in [-10, 7]), AND it overwinters: in the N hemisphere it
# crosses the year boundary (sdate + growp > 365); in the S hemisphere it straddles
# mid-winter (sows before and harvests after DOY 182). `valid` masks the NA/sdate<=0 rows
# the scalar guards with its outer `if` (returns 0 there). Returns 0/1.
.wintercrop_vec <- function(sdate_v, hdate_v, tcm_v, lat_v) {
  growp    <- ifelse(sdate_v <= hdate_v, hdate_v - sdate_v, 365L + hdate_v - sdate_v)
  valid    <- !is.na(sdate_v) & sdate_v > 0L & !is.na(lat_v) & !is.na(tcm_v)
  temp_ok  <- tcm_v >= -10 & tcm_v <= 7
  nh <- valid & lat_v >  0 & (sdate_v + growp > 365L) & (growp >= 150L) & temp_ok
  sh <- valid & lat_v <= 0 & (sdate_v < 182L) & (sdate_v + growp > 182L) & (growp >= 150L) & temp_ok
  as.integer(nh | sh)
}

# Shared window-sum core for the PHU models. Given a row-wise cumulative sum of daily
# effective thermal units (cumteff[, d] = sum of teff over days 1..d), return the sum over
# the growing window [sdate, hdate) WITHOUT a per-cell loop, as a difference of cumsums:
#   no wrap (sdate < hdate):  cum(hdate-1) - cum(sdate-1)            = sum[sdate .. hdate-1]
#   wrap    (hdate <= sdate): cum(365) - cum(sdate-1)  +  cum(hdate-1)  (tail + head of year)
# The pmin/pmax/<1 guards keep the column index in 1..365 and treat "sdate-1 = 0" as cum 0.
.phu_cumsum <- function(sdate_v, hdate_v, cumteff) {
  nc      <- nrow(cumteff)
  idx     <- seq_len(nc)
  no_wrap <- sdate_v < hdate_v
  hd1     <- hdate_v - 1L
  cum_hd1 <- ifelse(hd1 >= 1L, cumteff[cbind(idx, pmin(pmax(hd1, 1L), 365L))], 0)
  cum_hd1[hd1 < 1L] <- 0
  sd1     <- sdate_v - 1L
  cum_sd1 <- ifelse(sd1 >= 1L, cumteff[cbind(idx, pmin(pmax(sd1, 1L), 365L))], 0)
  cum_sd1[sd1 < 1L] <- 0
  cum_365 <- cumteff[, 365L]
  ifelse(no_wrap, cum_hd1 - cum_sd1, cum_365 - cum_sd1 + cum_hd1)
}

# PHU thermal model (phen_model "t"), vectorised twin of calcPHU: daily effective thermal
# units teff = max(T - basetemp, 0), cumulated, then summed over the growing window via
# .phu_cumsum. Positive heat-unit sum.
.calc_phu_thermal_vec <- function(sdate_v, hdate_v, temp_mat, basetemp) {
  teff    <- pmax(temp_mat - basetemp, 0)
  cumteff <- t(apply(teff, 1, cumsum))
  as.integer(.phu_cumsum(sdate_v, hdate_v, cumteff))
}

# PHU vernal-thermal model (phen_model "tv"), vectorised twin of calcPHU: as the thermal
# model but teff is scaled by the daily vernalization reduction factor (vrf_mat) before
# accumulation. Returned NEGATED -- LPJmL reads a negative PHU as "this crop needs
# vernalization".
.calc_phu_vernal_vec <- function(sdate_v, hdate_v, temp_mat, vrf_mat, basetemp) {
  teff    <- pmax(temp_mat - basetemp, 0) * vrf_mat
  cumteff <- t(apply(teff, 1, cumsum))
  as.integer(-.phu_cumsum(sdate_v, hdate_v, cumteff))
}

# Vernalization reduction factor matrix (NCELLS x 365), vectorised twin of calcVrf.
# Per day, vernalization EFFECTIVENESS veff(T) is a trapezoid in temperature: 0 below tv1,
# ramping up to 1 across [tv1,tv2], a plateau of 1 on [tv2,tv3], ramping back to 0 across
# [tv3,tv4], 0 above. The year is doubled (730 d) so a season wrapping past Dec 31 stays
# contiguous, and veff is cumulated once per row. Then per cell (the loop): accumulate veff
# from sdate, find the day the required vernalization vd is reached, and set the reduction
# factor vrf along the way -- 0 until vd_b (20%) of vd is banked, then a linear ramp to 1 at
# full vd (so growth is throttled until enough cold has accrued). vd <= 0 leaves vrf = 1
# (no requirement). The cumsum + difference (cumveff730[k] - cum0) replaces the scalar's
# day-by-day running sum; results match.
.build_vrf_mat <- function(sdate_v, hdate_v, temp_mat, vd_vec,
                           vd_b = 0.2, tv1, tv2, tv3, tv4) {
  nc       <- nrow(temp_mat)
  temp730  <- cbind(temp_mat, temp_mat)
  veff730  <- pmin(pmax(
    ifelse(temp730 < tv1,  0,
    ifelse(temp730 < tv2,  (temp730 - tv1) / (tv2 - tv1),
    ifelse(temp730 <= tv3, 1,
    ifelse(temp730 < tv4,  (tv4 - temp730) / (tv4 - tv3), 0)))), 0), 1)
  cumveff730 <- t(apply(veff730, 1, cumsum))

  vrf_mat <- matrix(1.0, nc, 365L)

  for (i in seq_len(nc)) {
    vd <- vd_vec[i]
    if (vd <= 0) next

    sd  <- sdate_v[i]
    hd  <- if (sd < hdate_v[i]) hdate_v[i] else hdate_v[i] + 365L
    cum0 <- if (sd > 1L) cumveff730[i, sd - 1L] else 0

    k_seq   <- sd:min(hd, 730L)
    cumvd   <- cumveff730[i, k_seq] - cum0
    end_pos <- which(cumvd >= vd)[1L]
    endday  <- if (!is.na(end_pos)) k_seq[end_pos] else Inf

    last <- as.integer(min(endday, hd))
    if (!is.finite(last)) next

    for (k in sd:last) {
      vdsum <- cumveff730[i, k] - cum0
      k_mod <- if (k > 365L) k - 365L else k
      vrf_mat[i, k_mod] <-
        if (vdsum < vd * vd_b) 0
        else max(0, min(1, (vdsum - vd * vd_b) / (vd * (1 - vd_b))))
    }
  }
  vrf_mat
}
