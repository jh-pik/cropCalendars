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