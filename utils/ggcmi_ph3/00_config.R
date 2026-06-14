# ---------------------------------------------------------------------------- #
# Configuration for GGCMI phase3 (ISIMIP3b climate)

# Author:  Sara Minoli
# Email:   sara.minoli@pik-potsdam.de
# ---------------------------------------------------------------------------- #

library(ncdf4)
library(abind)
library(data.table)
library(foreach)
library(cropCalendars)
library(zoo)           # for rolling mean
#library(unix)

# ------------------------------------ #
# General Settings

# Deployment paths and SLURM account live in settings.sh, so the .sh job scripts
# and this config share a single source of truth. Parse the KEY=VALUE lines here.
settings_file <- if (exists("work_dir")) file.path(work_dir, "settings.sh") else "settings.sh"
.settings <- local({
  lines <- readLines(settings_file)
  lines <- lines[grepl("^[A-Za-z_][A-Za-z0-9_]*=", lines)]   # keep KEY=VALUE lines
  keys  <- sub("=.*$", "", lines)
  vals  <- gsub("(^[\"']|[\"']$)", "", trimws(sub("^[^=]*=", "", lines)))
  setNames(as.list(vals), keys)
})

# Output directory: where output data are going to be saved
output_dir <- .settings$OUTPUT_DIR

parallel     <- TRUE
cluster_job  <- TRUE
plot_results <- TRUE

# ------------------------------------ #
# Tunable parameters (annual sliding-window pipeline). The env vars EMIT_STEP and
# PHU_SMOOTH_WINDOW override clm_emit_step / phu_smooth_window for ad-hoc runs.
clm_avg_years     <- 30       # climate-averaging window length (years)
clm_emit_step     <- 1        # compute calendars every N years (1 = fully annual)
pet_method        <- "fao56"  # PET: "fao56" (Penman-Monteith) or "pt" (Priestley-Taylor)
phu_smooth_window <- 1        # PHU temperature-averaging window (years; 1 = period-exact)
# Daily-climatology smoothing + sustained-crossing guard. These attack the
# year-to-year sowing/harvest oscillation at its source: the per-DOY daily means
# (dtemp/dprec/dpet) carry ~1 degC / spiky day-to-day jitter, and the point
# detectors (calcDoyCrossThreshold spring/fall + wet-season-end crossings) latch
# onto single-day blips -> spurious crossings ~130 days off (e.g. a 1-day dip
# through temp_spring=14 in the autumn descent). clm_smooth_window applies a
# centred circular running mean to the daily climatologies; cross_min_duration
# requires a crossing to persist that many days before it counts.
clm_smooth_window <- 15L      # daily-climatology smoothing window (days, odd; 0/1 = off)
cross_min_duration <- 5L      # min sustained-excursion days for calcDoyCrossThreshold (1 = off)

# Wettest-window hysteresis (distance-weighted, max-normalised selection). For the
# ~57% of PREC/PRECTEMP cells with a near-tied second 120-day P/PET peak, the plain
# argmax flips between far-apart peaks year to year. The new rule picks
# argmax( (ws/max ws) * (1 - eps*dist_to_last_year/182.5) ): near peaks essentially
# free (tracks drift), far peaks must be decisively better to win. eps=0.5 cut the
# wet-cell mean year-to-year jump 1.59 -> 0.13 and cells-ever-flipping 0.257 -> 0.077.
wet_window_eps    <- 0.5      # wettest-window distance-weighting strength (0 = plain argmax)

# Seasonality-classifier hysteresis (threshold deadband). ~18% of cells flip their
# seasonality CLASS year to year (grazing the CV_prec/CV_temp/min_temp thresholds),
# which swaps the whole sowing rule. seas_eps relaxes the relative CV thresholds toward
# keeping last year's class (thermostat deadband); seas_mtemp_margin is the absolute
# deadband (deg C) on the min-temp threshold. seas_eps=0.25 cut class flips 18.3% -> 2.8%.
seas_eps          <- 0.25     # seasonality CV-threshold deadband (rel.; 0 = off)
seas_mtemp_margin <- 1        # seasonality min-temp threshold deadband (deg C)

climate_dirs  <- sub("/+$", "", strsplit(.settings$CLIMATE_DIR, ":")[[1]])  # search list
climate_dir   <- climate_dirs[1]                                            # legacy (stage 01a)
isimip3b.path <- .settings$ISIMIP3B_PATH # .clm climate
agmip_dir     <- .settings$AGMIP_DIR     # AgMIP reference crop calendars (used in 02)
grid_file     <- .settings$GRID_BIN      # LPJmL grid path; read only in 03 (.clm writing)

# Climate input forcings (ESMs + the GSWP3-W5E5 observational forcing).
gcms <- c(
  "GSWP3-W5E5",
  "GFDL-ESM4",
  "IPSL-CM6A-LR",
  "MPI-ESM1-2-HR",
  "MRI-ESM2-0",
  "UKESM1-0-LL"
)

# Scenarios available per forcing -- the (GCM x scenario) processing matrix.
scenarios <- list(
  "GSWP3-W5E5"    = c("spinclim", "obsclim"),
  "GFDL-ESM4"     = c("historical", "ssp126", "ssp245", "ssp370", "ssp585"),
  "IPSL-CM6A-LR"  = c("historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp460", "ssp585"),
  "MPI-ESM1-2-HR" = c("historical", "ssp126", "ssp245", "ssp370", "ssp585"),
  "MRI-ESM2-0"    = c("historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp460", "ssp585"),
  "UKESM1-0-LL"   = c("historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585")
)
# Union of all scenarios, used to build the file-window lists below.
scens <- sort(unique(unlist(scenarios, use.names = FALSE)))

# NB: ensemble members and the per-scenario climate-file year windows are no longer
# hardcoded. The annual pipeline (01_compute_annual_calendars.R) discovers the
# climate files by globbing <CLIMATE_DIR>/<scenario>/<gcm>/ and reads the ensemble
# member, scenario tag and year range straight from the file names, so new datasets
# or scenarios need no config edits. (The legacy stage 01a -- now superseded by
# 01_compute_annual_calendars.R -- was the only consumer of enms/syears/eyears.)

# NB: the LPJmL grid is NOT read here. 01a/01b/02 work on the climate land mask; only
# stage 03 (writing the .clm files) needs the LPJmL grid, and reads it there.

# ------------------------------------------------------#
# Crop Names: ----
crop_ls <- list(all_low = c("winter_wheat", "spring_wheat", "maize", "rice1", "rice2",
                            "soybean", "millet", "sorghum","peas","sugar_beat",
                            "cassava","rape_seed","sunflower","nuts","sugarcane"),
                rb_cal  = c("Winter_Wheat", "Spring_Wheat", "Maize", "Rice", NA,
                            "Soybean", "Millet", "Sorghum", NA, NA,
                            NA, NA, NA, NA, NA),
                ggcmi   = c("wwh","swh","mai","ri1","ri2",
                            "soy","mil","sor","pea","sgb",
                            "cas","rap","sun","nut","sgc"),
                # vernal: yes_all = vern. forced in all grid cells;
                #         yes = only if conditions are met, see wintercrop()
                vernal  = c("yes_all","no","no","no","no",
                            "no","no","no","no","no",
                            "no","yes","no","no","no"))

irri_ls <- list(all_low = c("rainfed", "irrigated"),
                rb_cal  = c("Rainfed", "Irrigated"),
                ggcmi   = c("rf", "ir"))

# ------------------------------------ #
cat("\nConfigs imported.\n")