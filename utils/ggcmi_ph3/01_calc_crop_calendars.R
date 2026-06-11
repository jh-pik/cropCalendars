# ---------------------------------------------------------------------------- #
# Calculate crop calendars for GGCMI phase3 (ISIMIP3b climate)
# run via job script (.sh)

# Author:  Sara Minoli
# Email:   sara.minoli@pik-potsdam.de
# ---------------------------------------------------------------------------- #

rm(list = ls(all = TRUE))

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores()))
cat("\nNumber of cores: ", num_cores, "\n")

log_status <- function(stage) {
  cat("==== ", stage, " ====\n")
  conns <- showConnections(all = TRUE)
  cat("Number of open connections: ", if (is.null(conns)) 0 else nrow(conns), "\n")
}

log_status("Start of script")


starttime <- Sys.time() # Track run-time
print(starttime)

# ------------------------------------ #
# General settings
# Run from the pipeline dir: sbatch passes --chdir=$WD, interactive runs cd there.
work_dir <- getwd()
source(file.path(work_dir, "00_config.R"))

log_status("read config")


# ------------------------------------ #
# Individual-run settings
if (cluster_job == TRUE) {
  # import argument from bash script
  options(echo = FALSE) # if you want see commands in output file
  args <- commandArgs(trailingOnly = TRUE)
} else {
  args <- c('GFDL-ESM4', 'historical', 'Millet', '2071')
}
print(args)

log_status("read args")

# ------------------------------------ #
# Select variable, crop, model, year
gcm    <- args[1]
scen   <- args[2]
cro    <- args[3]
year   <- as.numeric(args[4])
nnodes <- as.numeric(args[5])
ntasks <- as.numeric(args[6])

cat("\n", gcm, scen, cro, year, nnodes, ntasks, "\n")
showConnections(all = TRUE)

ncpus        <- ntasks * nnodes
n_lon_chunks <- 120
# --nodes=1 --ntasks-per-node=16 --exclusive
log_status("computed ncups")

# Output directory
dfout_dir <- paste0(output_dir, "/crop_calendars/DT/", scen, "/", gcm, "/")
if (!dir.exists(dfout_dir)) dir.create(dfout_dir, recursive = TRUE)
plot_dir  <- paste0(output_dir, "/crop_calendars/DT/", scen, "/", gcm, "/plots/")
if (!dir.exists(plot_dir) & plot_results) dir.create(plot_dir, recursive = TRUE)

# ------------------------------------ #
# Register cluster
if(parallel == TRUE) {
  library(foreach)
  library(doParallel)
  # The size of our cluster must match the number of CPUs allocated to us
  # by SLURM.
  # By default, R can see all CPUs, including those not allocated to us.
  #ncpus <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
  if (!exists("ncpus")) ncpus <- 120
  cl <- makeCluster(ncpus)
  registerDoParallel(cl)
  getDoParName()
  getDoParWorkers()
}

# ------------------------------------ #
# Loop through all GCMs and scenarios and list all files needed to
#  calculate crop calendars for one sim configuration.
vars <- c("tas", "pr", "rsds", "rlds", "huss", "sfcwind", "ps")
clm_file_list <- list()
# Build the file list only for the requested GCM/scenario (gcm/scen from args).
# NB: select by name, not by loop counter — the previous gcms[gg]/scens[sc] left
# gg/sc at the last config index, so every job read the last GCM's files.
for (vv in seq(length(vars))) {
  cat("\n", vars[vv], "\n-----")

  if (length(grep("ssp", scen)) > 0) {
    # If SSP, concatenate also the historical files (needed for the climate window)
    fnames1 <- paste_ggcmi3_clm_fname(
      path1         = climate_dir,
      path2         = "",
      clm_scenario  = "historical",
      clm_forcing   = gcm,
      ens_member    = enms[gcm],
      bias_adj      = "w5e5",
      clm_var       = vars[vv],
      extent        = "global",
      time_step     = "daily",
      start_year    = syears[["historical"]],
      end_year      = eyears[["historical"]],
      file_ext      = ".nc"
    )
  } else {
    fnames1 <- NULL
  }

  fnames2 <- paste_ggcmi3_clm_fname(
    path1         = climate_dir,
    path2         = "",
    clm_scenario  = scen,
    clm_forcing   = gcm,
    ens_member    = enms[gcm],
    bias_adj      = "w5e5",
    clm_var       = vars[vv],
    extent        = "global",
    time_step     = "daily",
    start_year    = syears[[scen]],
    end_year      = eyears[[scen]],
    file_ext      = ".nc"
  )

  # Check if all files exist
  cat(paste0("\n", gcm, " ", scen, "\t",
      all(file.exists(c(fnames1, fnames2)))))

  clm_file_list[[vars[vv]]][[gcm]][[scen]] <- c(fnames1, fnames2)
}

# Get all climate files needed for this scenario
fnames_tas     <- clm_file_list[["tas"]][[gcm]][[scen]]
fnames_pr      <- clm_file_list[["pr"]][[gcm]][[scen]]
fnames_rsds    <- clm_file_list[["rsds"]][[gcm]][[scen]]
fnames_rlds    <- clm_file_list[["rlds"]][[gcm]][[scen]]
fnames_huss    <- clm_file_list[["huss"]][[gcm]][[scen]]
fnames_sfcwind <- clm_file_list[["sfcwind"]][[gcm]][[scen]]
fnames_ps      <- clm_file_list[["ps"]][[gcm]][[scen]]

# ------------------------------------ #
# First and last year of crop calendar calculation in this scenario
ccal_first_year <- ccal_years[
  which.min(abs(min(syears[[scen]]) - ccal_years))
  ]
ccal_last_year  <- ccal_years[
  which.min(min(eyears[[scen]]) - ccal_years)
  ]
# Index of the years
ccal_y_idx <- (
  which(ccal_years == ccal_first_year):which(ccal_years == ccal_last_year)
  )
cat("\n", gcm, scen, "--- \nCrop calendars: ", ccal_years[ccal_y_idx])

# ------------------------------------ #
# Get first and last year of climate to calculate crop calendar this year
cat("\nNow calculating crop calendars for year: ", year, "\n")
syear <- year - clm_avg_years
eyear <- year - 1
idx_first_file <- grep(syear, fnames_tas)
idx_last_file  <- grep(eyear, fnames_tas)
cat("\n", fnames_tas[idx_first_file:idx_last_file])

if (length(idx_first_file) < 1) {
  stop("Climate file does not exist!")
}

# ------------------------------------ #
#
# In terminal: scontrol show config
# Look for: DefMemPerCPU = 3500 MB
# https://kb.northwestern.edu/page.php?id=81074
# After loading two lon column: print(pryr::mem_used()) 220 MB
# 3500/110

# ------------------------------------ #
# Split lon dimension in chunks to be run in parallel
lons <- seq(min(grid_df$lon), max(grid_df$lon), by = 0.5)
lon_matrix <- matrix(lons, nrow = length(lons)/n_lon_chunks)

#print(rlimit_all())

# ------------------------------------ #
# Loop through lon chunks
output_df <- foreach(lo        = seq_len(ncol(lon_matrix)),
                     .combine  = "rbind",
                     .inorder  = FALSE,
                     .packages = c("ncdf4", "abind", "cropCalendars"),
                     .verbose  = FALSE
                     ) %dopar% {

  # Log file to track foreach parallel loop
  log_file <- paste0(dfout_dir, "log_", lo, "_", syear, "_", eyear, ".txt")
  unlink(log_file)
  sink(log_file, append = TRUE)
  cat("\nDoing ", lo, " of ", ncol(lon_matrix), "\n")

  # Subset grid for lons of this slice
  grid_sub <- subset(grid_df, lon %in% lon_matrix[, lo])

  nfiles  <- length(idx_first_file:idx_last_file)
  lon_sub <- list(lon = lon_matrix[, lo])

  tas_list <- pr_list <- rsds_list <- rlds_list <-
    huss_list <- sfcwind_list <- ps_list <- vector("list", nfiles)

  for (i in idx_first_file:idx_last_file) {
    ii <- i - idx_first_file + 1

    # Read and convert units on the fly (no intermediate copies)
    tas_list[[ii]]     <- k2deg(cropCalendars::readNcdf(fnames_tas[i],     dim_subset = lon_sub))
    pr_list[[ii]]      <- 86400 * cropCalendars::readNcdf(fnames_pr[i],      dim_subset = lon_sub)
    rsds_list[[ii]]    <- cropCalendars::readNcdf(fnames_rsds[i],    dim_subset = lon_sub)
    rlds_list[[ii]]    <- cropCalendars::readNcdf(fnames_rlds[i],    dim_subset = lon_sub)
    huss_list[[ii]]    <- cropCalendars::readNcdf(fnames_huss[i],    dim_subset = lon_sub)
    sfcwind_list[[ii]] <- cropCalendars::readNcdf(fnames_sfcwind[i], dim_subset = lon_sub)
    ps_list[[ii]]      <- cropCalendars::readNcdf(fnames_ps[i],      dim_subset = lon_sub)
  } # i

  # Single allocation: bind all time slices at once (avoids repeated copies)
  bind3 <- function(lst) do.call(abind, c(lst, list(along = 3, use.dnns = TRUE)))
  tas_array     <- bind3(tas_list)
  pr_array      <- bind3(pr_list)
  rsds_array    <- bind3(rsds_list)
  rlds_array    <- bind3(rlds_list)
  huss_array    <- bind3(huss_list)
  sfcwind_array <- bind3(sfcwind_list)
  ps_array      <- bind3(ps_list)
  rm(tas_list, pr_list, rsds_list, rlds_list, huss_list, sfcwind_list, ps_list)

  # Date sequence is the same for every pixel — compute once
  dates <- seqDates(
    start_date = paste0(syear, "-01-01"),
    end_date   = paste0(eyear, "-12-31"),
    step       = "day"
  )

  # Pre-allocate output list (avoids repeated O(n^2) rbind growth)
  ccal_list <- vector("list", nrow(grid_sub))

  for (j in seq_len(nrow(grid_sub))) {
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

    # Calculate monthly climate
    mclm <- calcMonthlyClimate(
      lat        = lat_pix,
      temp       = tas_pix,
      prec       = pr_pix,
      syear      = syear,
      eyear      = eyear,
      incl_feb29 = TRUE,
      pet_method = "fao56",
      swdown     = rsds_pix,
      lwdown     = rlds_pix,
      windspeed  = sfcwind_pix,
      humid      = huss_pix,
      ps         = ps_pix
      )

    # Calculate crop calendar
    ccal_list[[j]] <- calcCropCalendars(
      lon      = lon_pix,
      lat      = lat_pix,
      mclimate = mclm,
      crop     = cro
      )

    if (j %% 500 == 0) cat(j, "of", nrow(grid_sub), "\n")

  } # j

  rm(tas_array, pr_array, rsds_array, rlds_array, huss_array, sfcwind_array, ps_array)

  ccal_df <- do.call(rbind, ccal_list)
  unlink(log_file)
  return(ccal_df)
} # lo

cat("\nFinished foreach loop\n")
# Always close worker nodes, if running in parallel mode
if(parallel==T) {
  stopCluster(cl)
}

# ------------------------------------ #
# Save data table
DT <- data.table(output_df)
fnout <- paste0(dfout_dir, "DT_output_crop_calendars_",
                cro, "_", gcm, "_", scen, "_", syear, "_", eyear, ".Rdata")
cat("\n", fnout)
save(DT, file = fnout)


# Plot maps
if (plot_results) {
  fnout <- paste0(dfout_dir, "DT_output_crop_calendars_",
                  cro, "_", gcm, "_", scen, "_", syear, "_", eyear, ".Rdata")
  fnpdf <- paste0(plot_dir, "map_crop_calendars_",
                  cro, "_", gcm, "_", scen, "_", syear, "_", eyear, ".pdf")
  plotMapCropCalendars(fnDT = fnout, fnPDF = fnpdf)
}

# ------------------------------------ #
# Done!
endtime <- Sys.time()
print(endtime)

print(endtime-starttime)
