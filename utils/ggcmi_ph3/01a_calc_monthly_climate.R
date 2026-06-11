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
# cells (cell_lin, set below from the climate mask) as a [NCELLS, ndays] matrix. Reads
# one variable at a time and frees the full grid immediately, so memory stays ~one
# variable-year (fits the default per-cpu RAM).
read_year_cells <- function(fnames, yr, conv = identity) {
  rng <- file_year_range(fnames)
  wi  <- which(rng$fy <= yr & rng$ly >= yr)[1]
  if (is.na(wi)) stop("No climate file covering year ", yr)
  f       <- fnames[wi]
  dates_f <- seqDates(paste0(rng$fy[wi], "-01-01"), paste0(rng$ly[wi], "-12-31"), "day")
  pos     <- which(as.integer(substr(dates_f, 1, 4)) == yr) # 1-based day positions in file
  arr     <- cropCalendars::readNcdf(f, dim_subset = list(time = (pos[1] - 1):(pos[length(pos)] - 1)))
  conv(matrix(arr, nrow = dim(arr)[1] * dim(arr)[2])[cell_lin, , drop = FALSE])
}

# ------------------------------------ #
# Cell set = the CLIMATE land mask (all non-NA cells of the first tas file), so the
# crop-calendar product covers every land cell that has climate data. The LPJmL grid
# is NOT used here; it is applied later (stage 03) when writing the .clm files.
arr0     <- cropCalendars::readNcdf(clm_file_list[["tas"]][1], dim_subset = list(time = 0:0))
lon_axis <- as.numeric(dimnames(arr0)[[1]]); lat_axis <- as.numeric(dimnames(arr0)[[2]])
nlon     <- length(lon_axis)
cell_lin <- which(!is.na(matrix(arr0, nrow = nlon)))   # lon-major land-cell indices
NCELLS   <- length(cell_lin)
land_lon <- lon_axis[((cell_lin - 1L) %% nlon) + 1L]
land_lat <- lat_axis[((cell_lin - 1L) %/% nlon) + 1L]
rm(arr0)
cat("Climate land cells:", NCELLS, "\n")

# Cell-vectorised monthly-climate accumulator, shared with calcMonthlyClimate
# (init / addYear / finalize). Switch the PET method here: "fao56" (Penman-Monteith,
# 7 vars) or "pt" (Priestley-Taylor; uses radiation when supplied). 02/03 assume fao56.
pet_method <- "fao56"
acc        <- initMonthlyClimate(NCELLS, pet_method = pet_method)
valid      <- rep(TRUE, NCELLS)  # cells with any NA in tas/pr are dropped at the end

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
  # Pixel validity follows the original code: skip cells with any NA in tas or pr.
  valid <- valid & (rowSums(is.na(tas)) == 0) & (rowSums(is.na(pr)) == 0)

  # Accumulate this year into the shared engine (computes PET, monthly & DOY stats).
  dts <- seqDates(paste0(yr, "-01-01"), paste0(yr, "-12-31"), "day")
  acc <- addYearMonthlyClimate(
    acc, temp = tas, prec = pr, dates = dts,
    swdown = rsds, lwdown = rlds, windspeed = sfcwind, humid = huss, ps = ps,
    lat = land_lat
  )

  rm(tas, pr, rsds, rlds, huss, sfcwind, ps)
  gc(verbose = FALSE)
}
cat("\n")

# ------------------------------------ #
# Finalise via the shared engine, then keep only valid pixels (same set the original
# per-pixel code would have produced).
mclm <- finalizeMonthlyClimate(acc)
keep <- which(valid)
grid_clm   <- data.frame(lon = land_lon[keep], lat = land_lat[keep])
MTEMP      <- mclm$mtemp[keep, , drop = FALSE]
MPREC      <- mclm$mprec[keep, , drop = FALSE]
MPET       <- mclm$mpet[keep, , drop = FALSE]
MPPET      <- mclm$mppet[keep, , drop = FALSE]
MPPET_DIFF <- mclm$mppet_diff[keep, , drop = FALSE]
DTEMP      <- mclm$dtemp[keep, , drop = FALSE]
DPPET      <- mclm$dppet[keep, , drop = FALSE]

fnout <- paste0(clm_dir, "monthly_climate_", gcm, "_", scen, "_", syear, "_", eyear, ".Rdata")
cat("\nSaving monthly-climate cache:\n", fnout, "\n")
save(grid_clm, MTEMP, MPREC, MPET, MPPET, MPPET_DIFF, DTEMP, DPPET,
     gcm, scen, syear, eyear, file = fnout)
cat("Cached pixels: ", nrow(grid_clm), "\n")

endtime <- Sys.time()
print(endtime)
print(endtime - starttime)
