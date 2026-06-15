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
  raw   <- sub("^[^=]*=", "", lines)
  raw   <- sub("[[:space:]]+#.*$", "", raw)                  # strip trailing inline comment
  vals  <- gsub("(^[\"']|[\"']$)", "", trimws(raw))
  setNames(as.list(vals), keys)
})

# Fail loudly if a required deployment key is absent: a missing key would otherwise
# be NULL and either crash (strsplit on CLIMATE_DIR) or silently collapse in paste0
# to a root/relative path. Validate here, at the single config source of truth.
.required_settings <- c("OUTPUT_DIR", "CLIMATE_DIR", "ISIMIP3B_PATH",
                        "AGMIP_DIR", "GRID_BIN")
.missing_settings  <- setdiff(.required_settings, names(.settings))
if (length(.missing_settings) > 0)
  stop("settings.sh (", settings_file, ") is missing required key(s): ",
       paste(.missing_settings, collapse = ", "))

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
# Threshold-crossing noise filters (a matched pair). They attack the year-to-year
# sowing/harvest oscillation at its source: the per-DOY daily climatology
# (dtemp/dprec/dpet) carries ~1 degC / spiky day-to-day jitter, and the point
# detectors (calcDoyCrossThreshold: spring/fall temperature + wet-season-end
# crossings) otherwise latch onto a single-day blip -> a spurious crossing ~130 days
# off (e.g. a 1-day dip through temp_spring=14 in the autumn descent). The two knobs:
#   cross_smooth_window : centred circular running mean applied to the daily
#       climatologies BEFORE crossing detection (days, odd; 0/1 = off).
#   cross_min_duration  : a crossing must stay on the new side of the threshold for
#       at least this many days to count (1 = off).
# NB this is NOT a general climatology smoother. The window/extremum reductions
# (warmest/coldest-month temperature, driest-month P/PET, the 120-day wettest
# window) use their OWN structural 30-day (= 1 month) window -- that is what makes
# them monthly-equivalent -- which dominates cross_smooth_window, so those rules are
# ~invariant to it. The two smoothings therefore serve different rules and do not
# meaningfully compound.
cross_smooth_window <- 15L    # crossing-detector smoothing window (days, odd; 0/1 = off)
cross_min_duration  <- 5L     # min sustained-excursion days for calcDoyCrossThreshold (1 = off)

# Wettest-window hysteresis (distance-weighted, max-normalised selection). For the
# ~57% of PREC/PRECTEMP cells with a near-tied second 120-day P/PET peak, the plain
# argmax flips between far-apart peaks year to year. The new rule picks
# argmax( (ws/max ws) * (1 - eps*dist_to_last_year/182.5) ): near peaks essentially
# free (tracks drift), far peaks must be decisively better to win. eps=0.5 cut the
# wet-cell mean year-to-year jump 1.59 -> 0.13 and cells-ever-flipping 0.257 -> 0.077.
wet_window_eps    <- 0.5      # wettest-window weight FLOOR = 1-eps (0 = plain argmax); far peak must be >1/(1-eps)x wetter to win
wet_window_decay  <- 0.3      # Gaussian decay scale of the distance weight (units of half-year; smaller = stickier)

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

# ISIMIP3a observational-climate scenarios (vs ISIMIP3b ESM scenarios). Stage 02/03
# branch the DRS publish path on these: any forcing run under one of these scenarios
# (GSWP3-W5E5, 20CRv3-W5E5, ... -- the gcm stays a variable token) publishes under
# ISIMIP3a/.../crop_calendar/<soc>/ with soc = histsoc (obsclim/spinclim) or
# countersoc (counterclim), rather than the ISIMIP3b <gcm>/<soc> layout.
isimip3a_scenarios <- c("obsclim", "spinclim", "counterclim")

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