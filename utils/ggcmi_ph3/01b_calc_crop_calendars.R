# ---------------------------------------------------------------------------- #
# Step 1b: compute crop calendars from the cached monthly climate (01a output)
#
# Loads the GCM x scenario x window monthly-climate cache produced by
# 01a_calc_monthly_climate.R and runs calcCropCalendars() per pixel for one crop.
# Produces the same DT_output_crop_calendars_*.Rdata that stage 02 consumes, so
# downstream is unchanged.
#
# Args: GCM SCENARIO CROP YEAR NNODES NTASKS
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
  args <- c("GFDL-ESM4", "historical", "Maize", "1991", "1", "1")
}
print(args)

gcm    <- args[1]
scen   <- args[2]
cro    <- args[3]
year   <- as.numeric(args[4])
nnodes <- as.numeric(args[5])
ntasks <- as.numeric(args[6])

cat("\n", gcm, scen, cro, year, nnodes, ntasks, "\n")

ncpus <- ntasks * nnodes

# Climate window for this year (must match 01a)
syear <- year - clm_avg_years
eyear <- year - 1

# Output directory (same layout as the original stage 01)
dfout_dir <- paste0(output_dir, "/crop_calendars/DT/", scen, "/", gcm, "/")
if (!dir.exists(dfout_dir)) dir.create(dfout_dir, recursive = TRUE)
plot_dir  <- paste0(output_dir, "/crop_calendars/DT/", scen, "/", gcm, "/plots/")
if (!dir.exists(plot_dir) & plot_results) dir.create(plot_dir, recursive = TRUE)

# ------------------------------------ #
# Load the monthly-climate cache from 01a
clm_dir <- paste0(output_dir, "/crop_calendars/monthly_climate/", scen, "/", gcm, "/")
clm_fn  <- paste0(clm_dir, "monthly_climate_", gcm, "_", scen, "_", syear, "_", eyear, ".Rdata")
if (!file.exists(clm_fn)) {
  stop("Monthly-climate cache not found - run 01a first:\n  ", clm_fn)
}
cat("\nLoading monthly-climate cache:\n", clm_fn, "\n")
load(clm_fn) # grid_clm, MTEMP, MPREC, MPET, MPPET, MPPET_DIFF, DTEMP, DPPET
npix <- nrow(grid_clm)
cat("Pixels in cache: ", npix, "\n")

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
# Split pixels into one block per worker. Each block carries only its own slice
# of the cache matrices (avoids copying the full cache to every worker).
nblocks <- min(ncpus, npix)
blk_id  <- cut(seq_len(npix), breaks = nblocks, labels = FALSE)
blocks  <- split(seq_len(npix), blk_id)

block_data <- lapply(blocks, function(idx) {
  list(lon        = grid_clm$lon[idx],
       lat        = grid_clm$lat[idx],
       mtemp      = MTEMP[idx, , drop = FALSE],
       mprec      = MPREC[idx, , drop = FALSE],
       mpet       = MPET[idx, , drop = FALSE],
       mppet      = MPPET[idx, , drop = FALSE],
       mppet_diff = MPPET_DIFF[idx, , drop = FALSE],
       dtemp      = DTEMP[idx, , drop = FALSE],
       dppet      = DPPET[idx, , drop = FALSE])
})

# ------------------------------------ #
# Per-pixel crop calendars (light compute, no I/O)
output_df <- foreach(bd        = block_data,
                     .combine  = "rbind",
                     .inorder  = FALSE,
                     .packages = c("cropCalendars"),
                     .verbose  = FALSE
                     ) %dopar% {

  n   <- length(bd$lon)
  res <- vector("list", n)
  for (k in seq_len(n)) {
    mclm <- list(mtemp      = bd$mtemp[k, ],
                 mprec      = bd$mprec[k, ],
                 mpet       = bd$mpet[k, ],
                 mppet      = bd$mppet[k, ],
                 mppet_diff = bd$mppet_diff[k, ],
                 dtemp      = bd$dtemp[k, ],
                 dppet      = bd$dppet[k, ])

    res[[k]] <- calcCropCalendars(
      lon      = bd$lon[k],
      lat      = bd$lat[k],
      mclimate = mclm,
      crop     = cro
    )
  }
  do.call(rbind, res)
}

cat("\nFinished foreach loop\n")
if (parallel == TRUE) stopCluster(cl)

# ------------------------------------ #
# Save data table (same name/format as the original stage 01)
DT <- data.table(output_df)
fnout <- paste0(dfout_dir, "DT_output_crop_calendars_",
                cro, "_", gcm, "_", scen, "_", syear, "_", eyear, ".Rdata")
cat("\n", fnout)
save(DT, file = fnout)

# Plot maps
if (plot_results) {
  fnpdf <- paste0(plot_dir, "map_crop_calendars_",
                  cro, "_", gcm, "_", scen, "_", syear, "_", eyear, ".pdf")
  plotMapCropCalendars(fnDT = fnout, fnPDF = fnpdf)
}

endtime <- Sys.time()
print(endtime)
print(endtime - starttime)
