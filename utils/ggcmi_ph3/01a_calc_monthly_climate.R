# ---------------------------------------------------------------------------- #
# Step 1a: compute and cache the average monthly climate (GGCMI phase3 / ISIMIP3b)
#
# The monthly climate (and the daily DOY climatologies dtemp/dppet) depend only on
# GCM x scenario x climate-window — NOT on the crop. This script computes them once
# per GCM x scenario x year and caches them, so the per-crop step (01b) can reuse the
# result instead of re-reading daily climate and recomputing PET for every crop.
#
# Design (matches the standalone pipeline): stream over years, holding only one year
# of the full grid in memory at a time, and vectorise FAO-56 PET over all cells. This
# bounds memory to ~one year and needs no spatial chunking, so it runs single-threaded.
# (If the 7-variable I/O ever dominates, parallelise over years via a SLURM job array.)
#
# Output is identical by construction to the per-pixel calcMonthlyClimate() path; the
# aggregation below reproduces each cached field exactly (see equivalences inline).
#
# Args: GCM SCENARIO YEAR   (single-threaded; NNODES/NTASKS accepted but ignored)
# ---------------------------------------------------------------------------- #

rm(list = ls(all.names = TRUE))

starttime <- Sys.time()
print(starttime)

# ------------------------------------ #
# General settings
# Run from the pipeline dir: sbatch passes --chdir=$WD, interactive runs cd there.
work_dir <- getwd()
source(file.path(work_dir, "00_config.R"))

# ------------------------------------ #
# Individual-run settings
if (cluster_job == TRUE) {
  options(echo = FALSE)
  args <- commandArgs(trailingOnly = TRUE)
} else {
  args <- c("GFDL-ESM4", "historical", "1991")
}
print(args)

gcm  <- args[1]
scen <- args[2]
year <- as.numeric(args[3])

cat("\n", gcm, scen, year, "\n")

# Cache output directory (parallel to the DT tree, crop-independent)
clm_dir <- paste0(output_dir, "/crop_calendars/monthly_climate/", scen, "/", gcm, "/")
if (!dir.exists(clm_dir)) dir.create(clm_dir, recursive = TRUE)

NCELLS <- nrow(grid_df)

# ------------------------------------ #
# Build the file list only for the requested GCM/scenario (select by name).
vars <- c("tas", "pr", "rsds", "rlds", "huss", "sfcwind", "ps")
clm_file_list <- list()
for (vv in seq_along(vars)) {
  if (length(grep("ssp", scen)) > 0) {
    # If SSP, concatenate also the historical files (needed for the climate window)
    fnames1 <- paste_ggcmi3_clm_fname(
      path1 = climate_dir, path2 = "", clm_scenario = "historical",
      clm_forcing = gcm, ens_member = enms[gcm], bias_adj = "w5e5",
      clm_var = vars[vv], extent = "global", time_step = "daily",
      start_year = syears[["historical"]], end_year = eyears[["historical"]],
      file_ext = ".nc"
    )
  } else {
    fnames1 <- NULL
  }
  fnames2 <- paste_ggcmi3_clm_fname(
    path1 = climate_dir, path2 = "", clm_scenario = scen,
    clm_forcing = gcm, ens_member = enms[gcm], bias_adj = "w5e5",
    clm_var = vars[vv], extent = "global", time_step = "daily",
    start_year = syears[[scen]], end_year = eyears[[scen]],
    file_ext = ".nc"
  )
  clm_file_list[[vars[vv]]] <- c(fnames1, fnames2)
}

# ------------------------------------ #
# Climate window for this year
syear  <- year - clm_avg_years
eyear  <- year - 1
years  <- syear:eyear
nyears <- length(years)
cat("\nClimate window:", syear, "-", eyear, "(", nyears, "years )\n")

# Per file, the [first_year, last_year] it covers (parsed from the filename).
file_year_range <- function(fnames) {
  m  <- regmatches(fnames, regexpr("[0-9]{4}_[0-9]{4}\\.nc$", fnames))
  fy <- as.integer(substr(m, 1, 4))
  ly <- as.integer(substr(m, 6, 9))
  list(fy = fy, ly = ly)
}
yr_tas <- file_year_range(clm_file_list[["tas"]])
if (!any(yr_tas$fy <= syear & yr_tas$ly >= syear)) stop("Climate file for ", syear, " not found.")

# Guard: this script assumes a real, leap-aware calendar. seqDates() builds a
# proleptic-Gregorian day sequence and the per-year day positions / leap handling
# depend on it; a noleap (365_day) or 360_day dataset would silently mis-position
# days. Check the CF time:calendar attribute of the first file and refuse if it is
# not leap-aware (in case a different climate dataset is plugged in).
accepted_cal <- c("proleptic_gregorian", "standard", "gregorian")
nc_chk <- ncdf4::nc_open(clm_file_list[["tas"]][1])
cal_att <- ncdf4::ncatt_get(nc_chk, "time", "calendar")
ncdf4::nc_close(nc_chk)
if (!isTRUE(cal_att$hasatt)) {
  warning("No time:calendar attribute found — assuming a leap-aware calendar.")
} else if (!(tolower(cal_att$value) %in% accepted_cal)) {
  stop("Unsupported time:calendar = '", cal_att$value, "' in\n  ",
       clm_file_list[["tas"]][1],
       "\n  01a assumes a leap-aware calendar (", paste(accepted_cal, collapse = ", "),
       "); a noleap/365_day/360_day dataset needs different day handling.")
} else {
  cat("time:calendar =", cal_att$value, "\n")
}

# Read one calendar year of one variable for the full grid, returning only the land
# cells as a [NCELLS, ndays] matrix. Reads one variable at a time and frees the full
# grid immediately, so memory stays ~one variable-year (fits the default per-cpu RAM).
# Sets cell_lin (lon-major land-cell index) on first use, from the array dimnames.
read_year_cells <- function(fnames, yr, conv = identity) {
  rng <- file_year_range(fnames)
  wi  <- which(rng$fy <= yr & rng$ly >= yr)[1]
  if (is.na(wi)) stop("No climate file covering year ", yr)
  f       <- fnames[wi]
  dates_f <- seqDates(paste0(rng$fy[wi], "-01-01"), paste0(rng$ly[wi], "-12-31"), "day")
  pos     <- which(as.integer(substr(dates_f, 1, 4)) == yr) # 1-based day positions in file
  arr     <- cropCalendars::readNcdf(f, dim_subset = list(time = (pos[1] - 1):(pos[length(pos)] - 1)))
  if (is.null(cell_lin)) {
    dn       <- dimnames(arr)
    arr_lon  <- as.numeric(dn[[1]]); arr_lat <- as.numeric(dn[[2]])
    cell_lin <<- match(grid_df$lon, arr_lon) + (match(grid_df$lat, arr_lat) - 1) * length(arr_lon)
  }
  conv(matrix(arr, nrow = dim(arr)[1] * dim(arr)[2])[cell_lin, , drop = FALSE])
}

# ------------------------------------ #
# Accumulators (NCELLS rows). Monthly: sum over years of the per-year monthly
# statistic (mean temp / sum prec / sum pet / sum-ratio ppet). DOY: per-DOY sum and
# count across all days of all years, so DTEMP/DPPET = sum/count == tapply(.,mean).
M_tas  <- M_pr <- M_pet <- M_ppet <- matrix(0, NCELLS, 12)
D_tsum <- D_psum <- matrix(0, NCELLS, 365)
D_cnt  <- numeric(365)
valid  <- rep(TRUE, NCELLS)
cell_lin <- NULL # lon-major linear index of the land cells into the [720,360] grid

cat("Streaming", nyears, "years (PET vectorised over", NCELLS, "cells)...\n")
for (yr in years) {
  cat(yr, "")

  tas     <- read_year_cells(clm_file_list[["tas"]],     yr, conv = k2deg)
  pr      <- read_year_cells(clm_file_list[["pr"]],      yr, conv = function(x) 86400 * x)
  rsds    <- read_year_cells(clm_file_list[["rsds"]],    yr)
  rlds    <- read_year_cells(clm_file_list[["rlds"]],    yr)
  huss    <- read_year_cells(clm_file_list[["huss"]],    yr)
  sfcwind <- read_year_cells(clm_file_list[["sfcwind"]], yr)
  ps      <- read_year_cells(clm_file_list[["ps"]],      yr)
  nd      <- ncol(tas)

  # FAO-56 PET, vectorised over all cells x days at once (calcPET_FAO56 is elementwise)
  pet  <- calcPET_FAO56(tas, sfcwind, huss, rsds, rlds, ps)
  ppet <- pr / pmax(pet, 1e-6)
  # Only tas/pr/pet/ppet feed the accumulation below; free the PET-only inputs now.
  rm(rsds, rlds, huss, sfcwind, ps)

  # Pixel validity follows the original code: skip cells with any NA in tas or pr.
  valid <- valid & (rowSums(is.na(tas)) == 0) & (rowSums(is.na(pr)) == 0)

  # Month / DOY labels for this year (same helpers calcMonthlyClimate uses)
  dts <- seqDates(paste0(yr, "-01-01"), paste0(yr, "-12-31"), "day")
  mon <- date_to_month(dts)
  doy <- date_to_doy(dts, skip_feb29 = TRUE)

  # Monthly accumulation (per-year statistic, summed across years; /nyears at the end)
  for (m in 1:12) {
    dom     <- which(mon == m)
    pr_mon  <- rowSums(pr[,  dom, drop = FALSE])
    pet_mon <- rowSums(pet[, dom, drop = FALSE])
    M_tas[,  m] <- M_tas[,  m] + rowMeans(tas[, dom, drop = FALSE])
    M_pr[,   m] <- M_pr[,   m] + pr_mon
    M_pet[,  m] <- M_pet[,  m] + pet_mon
    M_ppet[, m] <- M_ppet[, m] + pr_mon / pet_mon
  }

  # DOY accumulation: add each day to its DOY column (a DOY can recur within a leap
  # year, e.g. DOY 28), and count per DOY -> exact tapply(., DOY, mean) replication.
  for (k in seq_len(nd)) {
    d <- doy[k]
    D_tsum[, d] <- D_tsum[, d] + tas[, k]
    D_psum[, d] <- D_psum[, d] + ppet[, k]
  }
  D_cnt <- D_cnt + tabulate(doy, nbins = 365)

  rm(tas, pr, pet, ppet)
  gc(verbose = FALSE)
}
cat("\n")

# ------------------------------------ #
# Multi-year averages (rounding matches calcMonthlyClimate: monthly rounded to 5
# digits, mppet_diff derived from the rounded mppet; daily fields not rounded).
MTEMP      <- round(M_tas  / nyears, 5)
MPREC      <- round(M_pr   / nyears, 5)
MPET       <- round(M_pet  / nyears, 5)
MPPET      <- round(M_ppet / nyears, 5)
MPPET_DIFF <- MPPET - MPPET[, c(2:12, 1)]
DTEMP      <- sweep(D_tsum, 2, D_cnt, "/")
DPPET      <- sweep(D_psum, 2, D_cnt, "/")

# Keep valid pixels only (same set the original per-pixel code would have produced)
keep <- which(valid)
grid_clm   <- grid_df[keep, c("lon", "lat")]
MTEMP      <- MTEMP[keep, , drop = FALSE]
MPREC      <- MPREC[keep, , drop = FALSE]
MPET       <- MPET[keep, , drop = FALSE]
MPPET      <- MPPET[keep, , drop = FALSE]
MPPET_DIFF <- MPPET_DIFF[keep, , drop = FALSE]
DTEMP      <- DTEMP[keep, , drop = FALSE]
DPPET      <- DPPET[keep, , drop = FALSE]

fnout <- paste0(clm_dir, "monthly_climate_", gcm, "_", scen, "_", syear, "_", eyear, ".Rdata")
cat("\nSaving monthly-climate cache:\n", fnout, "\n")
save(grid_clm, MTEMP, MPREC, MPET, MPPET, MPPET_DIFF, DTEMP, DPPET,
     gcm, scen, syear, eyear, file = fnout)
cat("Cached pixels: ", nrow(grid_clm), "\n")

endtime <- Sys.time()
print(endtime)
print(endtime - starttime)
