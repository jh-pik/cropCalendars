# ---------------------------------------------------------------------------- #
# Step 1a: compute and cache the average monthly climate (GGCMI phase3 / ISIMIP3b)
#
# The monthly climate (and the daily DOY climatologies dtemp/dppet) depend only on
# GCM x scenario x climate-window — NOT on the crop. This script computes them once
# per GCM x scenario x year and caches them, so the per-crop step (01b) can reuse the
# result instead of re-reading daily climate and recomputing PET for every crop.
#
# Args: GCM SCENARIO YEAR NNODES NTASKS   (no crop)
# ---------------------------------------------------------------------------- #

rm(list = ls(all.names = TRUE))

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores()))
cat("\nNumber of cores: ", num_cores, "\n")

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
  args <- c("GFDL-ESM4", "historical", "1991", "1", "1")
}
print(args)

gcm    <- args[1]
scen   <- args[2]
year   <- as.numeric(args[3])
nnodes <- as.numeric(args[4])
ntasks <- as.numeric(args[5])

cat("\n", gcm, scen, year, nnodes, ntasks, "\n")

ncpus        <- ntasks * nnodes
n_lon_chunks <- 120

# Cache output directory (parallel to the DT tree, crop-independent)
clm_dir <- paste0(output_dir, "/crop_calendars/monthly_climate/", scen, "/", gcm, "/")
if (!dir.exists(clm_dir)) dir.create(clm_dir, recursive = TRUE)

# ------------------------------------ #
# Register cluster
if (parallel == TRUE) {
  library(foreach)
  library(doParallel)
  if (!exists("ncpus")) ncpus <- 120
  cl <- makeCluster(ncpus)
  registerDoParallel(cl)
}

# ------------------------------------ #
# Build the file list only for the requested GCM/scenario (select by name).
vars <- c("tas", "pr", "rsds", "rlds", "huss", "sfcwind", "ps")
clm_file_list <- list()
for (vv in seq(length(vars))) {

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

fnames_tas     <- clm_file_list[["tas"]]
fnames_pr      <- clm_file_list[["pr"]]
fnames_rsds    <- clm_file_list[["rsds"]]
fnames_rlds    <- clm_file_list[["rlds"]]
fnames_huss    <- clm_file_list[["huss"]]
fnames_sfcwind <- clm_file_list[["sfcwind"]]
fnames_ps      <- clm_file_list[["ps"]]

# ------------------------------------ #
# Climate window for this year
syear <- year - clm_avg_years
eyear <- year - 1
idx_first_file <- grep(syear, fnames_tas)
idx_last_file  <- grep(eyear, fnames_tas)
cat("\n", fnames_tas[idx_first_file:idx_last_file])
if (length(idx_first_file) < 1) stop("Climate file does not exist!")

# ------------------------------------ #
# Split lon dimension in chunks to be run in parallel
lons       <- seq(min(grid_df$lon), max(grid_df$lon), by = 0.5)
lon_matrix <- matrix(lons, nrow = length(lons) / n_lon_chunks)

# ------------------------------------ #
# Loop through lon chunks: read climate, compute monthly climate per pixel.
# Each chunk returns the per-pixel monthly-climate fields as matrices (valid pixels only).
clm_chunks <- foreach(lo        = seq_len(ncol(lon_matrix)),
                      .inorder  = FALSE,
                      .packages = c("ncdf4", "abind", "cropCalendars"),
                      .verbose  = FALSE
                      ) %dopar% {

  log_file <- paste0(clm_dir, "log_", lo, "_", syear, "_", eyear, ".txt")
  unlink(log_file)
  sink(log_file, append = TRUE)
  cat("\nDoing ", lo, " of ", ncol(lon_matrix), "\n")

  grid_sub <- subset(grid_df, lon %in% lon_matrix[, lo])
  nfiles   <- length(idx_first_file:idx_last_file)
  lon_sub  <- list(lon = lon_matrix[, lo])

  tas_list <- pr_list <- rsds_list <- rlds_list <-
    huss_list <- sfcwind_list <- ps_list <- vector("list", nfiles)

  for (i in idx_first_file:idx_last_file) {
    ii <- i - idx_first_file + 1
    tas_list[[ii]]     <- k2deg(cropCalendars::readNcdf(fnames_tas[i],     dim_subset = lon_sub))
    pr_list[[ii]]      <- 86400 * cropCalendars::readNcdf(fnames_pr[i],      dim_subset = lon_sub)
    rsds_list[[ii]]    <- cropCalendars::readNcdf(fnames_rsds[i],    dim_subset = lon_sub)
    rlds_list[[ii]]    <- cropCalendars::readNcdf(fnames_rlds[i],    dim_subset = lon_sub)
    huss_list[[ii]]    <- cropCalendars::readNcdf(fnames_huss[i],    dim_subset = lon_sub)
    sfcwind_list[[ii]] <- cropCalendars::readNcdf(fnames_sfcwind[i], dim_subset = lon_sub)
    ps_list[[ii]]      <- cropCalendars::readNcdf(fnames_ps[i],      dim_subset = lon_sub)
  }

  bind3 <- function(lst) do.call(abind, c(lst, list(along = 3, use.dnns = TRUE)))
  tas_array     <- bind3(tas_list)
  pr_array      <- bind3(pr_list)
  rsds_array    <- bind3(rsds_list)
  rlds_array    <- bind3(rlds_list)
  huss_array    <- bind3(huss_list)
  sfcwind_array <- bind3(sfcwind_list)
  ps_array      <- bind3(ps_list)
  rm(tas_list, pr_list, rsds_list, rlds_list, huss_list, sfcwind_list, ps_list)

  dates <- seqDates(start_date = paste0(syear, "-01-01"),
                    end_date   = paste0(eyear, "-12-31"),
                    step       = "day")

  ng <- nrow(grid_sub)
  mtemp_m <- mprec_m <- mpet_m <- mppet_m <- mppet_diff_m <- matrix(NA_real_, ng, 12)
  dtemp_m <- dppet_m <- matrix(NA_real_, ng, 365)
  valid   <- logical(ng)

  for (j in seq_len(ng)) {
    lon_pix <- grid_sub$lon[j]
    lat_pix <- grid_sub$lat[j]

    tas_pix     <- tas_array[as.character(lon_pix),     as.character(lat_pix), ]
    pr_pix      <- pr_array[as.character(lon_pix),      as.character(lat_pix), ]
    rsds_pix    <- rsds_array[as.character(lon_pix),    as.character(lat_pix), ]
    rlds_pix    <- rlds_array[as.character(lon_pix),    as.character(lat_pix), ]
    huss_pix    <- huss_array[as.character(lon_pix),    as.character(lat_pix), ]
    sfcwind_pix <- sfcwind_array[as.character(lon_pix), as.character(lat_pix), ]
    ps_pix      <- ps_array[as.character(lon_pix),      as.character(lat_pix), ]

    if (any(is.na(tas_pix)) | any(is.na(pr_pix))) {
      cat("\nMissing values, skipping ", lon_pix, lat_pix)
      next
    }

    names(tas_pix) <- names(pr_pix) <- dates

    mclm <- calcMonthlyClimate(
      lat = lat_pix, temp = tas_pix, prec = pr_pix,
      syear = syear, eyear = eyear, incl_feb29 = TRUE,
      pet_method = "fao56",
      swdown = rsds_pix, lwdown = rlds_pix,
      windspeed = sfcwind_pix, humid = huss_pix, ps = ps_pix
    )

    mtemp_m[j, ]      <- mclm$mtemp
    mprec_m[j, ]      <- mclm$mprec
    mpet_m[j, ]       <- mclm$mpet
    mppet_m[j, ]      <- mclm$mppet
    mppet_diff_m[j, ] <- mclm$mppet_diff
    dtemp_m[j, ]      <- mclm$dtemp
    dppet_m[j, ]      <- mclm$dppet
    valid[j]          <- TRUE

    if (j %% 500 == 0) cat(j, "of", ng, "\n")
  }

  rm(tas_array, pr_array, rsds_array, rlds_array, huss_array, sfcwind_array, ps_array)

  keep <- which(valid)
  unlink(log_file)
  list(grid       = grid_sub[keep, c("lon", "lat")],
       mtemp      = mtemp_m[keep, , drop = FALSE],
       mprec      = mprec_m[keep, , drop = FALSE],
       mpet       = mpet_m[keep, , drop = FALSE],
       mppet      = mppet_m[keep, , drop = FALSE],
       mppet_diff = mppet_diff_m[keep, , drop = FALSE],
       dtemp      = dtemp_m[keep, , drop = FALSE],
       dppet      = dppet_m[keep, , drop = FALSE])
}

cat("\nFinished foreach loop\n")
if (parallel == TRUE) stopCluster(cl)

# ------------------------------------ #
# Assemble per-chunk results into one cache (valid pixels only)
grid_clm   <- do.call(rbind, lapply(clm_chunks, `[[`, "grid"))
MTEMP      <- do.call(rbind, lapply(clm_chunks, `[[`, "mtemp"))
MPREC      <- do.call(rbind, lapply(clm_chunks, `[[`, "mprec"))
MPET       <- do.call(rbind, lapply(clm_chunks, `[[`, "mpet"))
MPPET      <- do.call(rbind, lapply(clm_chunks, `[[`, "mppet"))
MPPET_DIFF <- do.call(rbind, lapply(clm_chunks, `[[`, "mppet_diff"))
DTEMP      <- do.call(rbind, lapply(clm_chunks, `[[`, "dtemp"))
DPPET      <- do.call(rbind, lapply(clm_chunks, `[[`, "dppet"))

fnout <- paste0(clm_dir, "monthly_climate_", gcm, "_", scen, "_", syear, "_", eyear, ".Rdata")
cat("\nSaving monthly-climate cache:\n", fnout, "\n")
save(grid_clm, MTEMP, MPREC, MPET, MPPET, MPPET_DIFF, DTEMP, DPPET,
     gcm, scen, syear, eyear, file = fnout)

cat("\nCached pixels: ", nrow(grid_clm), "\n")

endtime <- Sys.time()
print(endtime)
print(endtime - starttime)
