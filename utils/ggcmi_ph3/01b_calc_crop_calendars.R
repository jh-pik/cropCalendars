# ---------------------------------------------------------------------------- #
# Step 1b: compute crop calendars from the cached monthly climate (01a output)
#
# Loads the GCM x scenario x window monthly-climate cache produced by
# 01a_calc_monthly_climate.R and runs calcCropCalendars() per pixel for one crop.
# Produces the same DT_output_crop_calendars_*.Rdata that stage 02 consumes, so
# downstream is unchanged.
#
# Single-threaded by design: there is no climate I/O and no per-pixel PET here, just
# calcCropCalendars over the cached monthly climate. Parallelism is at the job level
# (one job per crop x year). If a single crop ever needs to be faster, split it with a
# SLURM job array; no in-script chunking is required.
#
# Args: GCM SCENARIO CROP YEAR
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
  args <- c("GFDL-ESM4", "historical", "Maize", "1991")
}
print(args)

gcm  <- args[1]
scen <- args[2]
cro  <- args[3]
year <- as.numeric(args[4])

cat("\n", gcm, scen, cro, year, "\n")

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
# Per-pixel crop calendars (single-threaded loop over the cached monthly climate)
ccal_list <- vector("list", npix)
for (j in seq_len(npix)) {
  mclm <- list(mtemp      = MTEMP[j, ],
               mprec      = MPREC[j, ],
               mpet       = MPET[j, ],
               mppet      = MPPET[j, ],
               mppet_diff = MPPET_DIFF[j, ],
               dtemp      = DTEMP[j, ],
               dppet      = DPPET[j, ],
               dprec      = DPREC[j, ],
               dpet       = DPET[j, ])

  ccal_list[[j]] <- calcCropCalendars(
    lon      = grid_clm$lon[j],
    lat      = grid_clm$lat[j],
    mclimate = mclm,
    crop     = cro
  )

  if (j %% 5000 == 0) cat(j, "of", npix, "\n")
}

output_df <- do.call(rbind, ccal_list)
cat("\nFinished pixel loop\n")

# ------------------------------------ #
# Save data table (same name/format as the original stage 01)
DT <- data.table(output_df)
fnout <- paste0(dfout_dir, "DT_output_crop_calendars_",
                cro, "_", gcm, "_", scen, "_", syear, "_", eyear, ".Rdata")
cat("\n", fnout)
save(DT, file = fnout)

# Plot maps (best-effort: the DT is already saved, so a plotting failure must not
# fail the job).
if (plot_results) {
  fnpdf <- paste0(plot_dir, "map_crop_calendars_",
                  cro, "_", gcm, "_", scen, "_", syear, "_", eyear, ".pdf")
  tryCatch(
    plotMapCropCalendars(fnDT = fnout, fnPDF = fnpdf),
    error = function(e) warning("Plotting failed (DT already saved): ",
                                conditionMessage(e))
  )
}

endtime <- Sys.time()
print(endtime)
print(endtime - starttime)
