# ---------------------------------------------------------------------------- #
# Calculate PHUs for LPJmL and create NCDF files (multi-year and crop-specific)
#  for GGCMI-phase3 adaptation runs

# Author:  Sara Minoli
# Email:   sara.minoli@pik-potsdam.de
# ---------------------------------------------------------------------------- #

rm(list = ls(all.names = TRUE))

starttime <- Sys.time() # Track run-time
print(starttime)

# ------------------------------------ #
# General settings
# Run from the pipeline dir: sbatch passes --chdir=$WD, interactive runs cd there.
work_dir <- getwd()
source(file.path(work_dir, "00_config.R"))

makeplot <- TRUE

# Read the LPJmL grid here (the only stage that needs it, for the .clm output).
# lpjmlkit::read_io auto-detects the header; round to the grid's native 0.01-degree
# resolution. grid_df / NCELLS are read as globals by generatePHUTserie_isimip3().
library(lpjmlkit)
grid_io <- suppressWarnings(read_io(grid_file, silent = TRUE))
grid_df <- data.frame(
  lon = round(as.numeric(grid_io$data[, 1, 1]), 2),
  lat = round(as.numeric(grid_io$data[, 1, 2]), 2)
)
NCELLS  <- nrow(grid_df)

# ------------------------------------ #
# Individual-run settings
if (cluster_job == TRUE) {
  # import argument from bash script
  options(echo = FALSE) # if you want see commands in output file
  args <- commandArgs(trailingOnly = TRUE)
} else {
  args <- c("GFDL-ESM4", "ssp585", "mai", "rf")
}
print(args)

# ------------------------------------ #
# Select variable, crop, model, year
gcm    <- args[1]
scen   <- args[2]
cro    <- args[3]
irri   <- args[4]

# Product year range (annual product), per scenario.
prod_range <- list("historical" = c(1850, 2014), "picontrol" = c(1601, 2100),
                   "ssp119" = c(2015, 2100), "ssp126" = c(2015, 2100),
                   "ssp245" = c(2015, 2100), "ssp370" = c(2015, 2100),
                   "ssp460" = c(2015, 2100), "ssp585" = c(2015, 2100))[[scen]]
FYnc <- prod_range[1]; LYnc <- prod_range[2]

# Annual product: one "period" per year, so the PHU matches each year's growing
# period exactly. PHU_SMOOTH_WINDOW (default 1) optionally widens the temperature
# averaging only. LYs is a global read by generatePHUTserie (= LYnc skips the
# legacy 2015gs branch).
SYs <- FYnc:LYnc; EYs <- FYnc:LYnc; LYs <- LYnc
smooth_window <- as.integer(Sys.getenv("PHU_SMOOTH_WINDOW", as.character(phu_smooth_window)))
cat(sprintf("PHU: %s %s %s_%s | years %d-%d | smooth_window=%d\n",
            gcm, scen, cro, irri, FYnc, LYnc, smooth_window))

years  <- FYnc:LYnc
nyears <- length(years)

ncdir  <- paste0(output_dir, "crop_calendars/ncdf/", gcm, "/", scen, "/")  # PHU .nc4 output
if (!dir.exists(ncdir)) dir.create(ncdir, recursive = TRUE)

# DRS crop-calendar file to read (output of stage 02).
soc_dir  <- if (scen == "historical") "historical" else paste0(scen, "soc-adapt")
soc_file <- if (scen == "historical") "histsoc" else scen
irr_tok  <- if (irri == "ir") "firr" else "noirr"
ncfile   <- paste0(output_dir, "ISIMIP3b/InputData/socioeconomic/crop_calendar/", gcm, "/", soc_dir,
                   "/ggcmi-crop-calendar_", tolower(gcm), "_", soc_file, "_", cro, "-", irr_tok,
                   "_annual_", FYnc, "_", LYnc, ".nc")

# ------------------------------------------------------#
# Compute PHUs and Write ncdfs

generatePHUTserie_isimip3(
    ncdir         = ncdir,
    gcm           = gcm,
    scen          = scen,
    cro           = cro,
    irri          = irri,
    SYs           = SYs,
    EYs           = EYs,
    FYnc          = FYnc,
    LYnc          = LYnc,
    crop_par_file = NULL,
    ncfile        = ncfile,
    smooth_window = smooth_window
)

# ------------------------------------------------------#

cat("\n", paste("Computation PHUs ended!"),
    "-------------------------------------------------------", sep = "\n")

endtime <- Sys.time()
print(endtime)

print(endtime-starttime)
