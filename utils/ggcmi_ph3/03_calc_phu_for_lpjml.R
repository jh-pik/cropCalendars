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
# resolution. grid_df is passed explicitly to generatePHUTserie_isimip3().
library(lpjmlkit)
grid_io <- suppressWarnings(read_io(grid_file, silent = TRUE))
grid_df <- data.frame(
  lon = round(as.numeric(grid_io$data[, 1, 1]), 2),
  lat = round(as.numeric(grid_io$data[, 1, 2]), 2)
)

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

# DRS crop-calendar file to read (output of stage 02). ISIMIP3a observational scenarios
# (obsclim/spinclim -> histsoc, counterclim -> countersoc; any forcing dataset) ->
# ISIMIP3a; ESMs -> ISIMIP3b per gcm x soc. The directory + filename stem MUST mirror
# stage 02 (02_assemble_annual_ncdf.R) exactly.
irr_tok  <- if (irri == "ir") "firr" else "noirr"
if (scen %in% isimip3a_scenarios) {
  soc_file <- if (scen == "counterclim") "countersoc" else "histsoc"
  cal_dir  <- paste0(output_dir, "ISIMIP3a/InputData/socioeconomic/crop_calendar/", soc_file, "/")
} else {
  soc_dir  <- if (scen == "historical") "historical" else paste0(scen, "soc-adapt")
  soc_file <- if (scen == "historical") "histsoc" else scen
  cal_dir  <- paste0(output_dir, "ISIMIP3b/InputData/socioeconomic/crop_calendar/", gcm, "/", soc_dir, "/")
}
cal_stem <- paste0("ggcmi-crop-calendar_", tolower(gcm), "_", soc_file, "_",
                   cro, "-", irr_tok, "_annual_")

# Derive the product year range from the file stage 02 actually wrote, instead of a
# hardcoded per-scenario table: stage 02 names the file from the discovered emit_years
# (min/max), so reading the range back from the file keeps the two stages on a single
# source of truth (any new dataset / partial run / unlisted scenario just works).
hits <- Sys.glob(paste0(cal_dir, cal_stem, "*.nc"))
if (length(hits) == 0)
  stop("Stage-02 crop-calendar file not found (run stage 02 first): ",
       cal_dir, cal_stem, "*.nc")
ncfile <- hits[1]
m <- regmatches(basename(ncfile),
                regexec("_annual_([0-9]{4})_([0-9]{4})\\.nc$", basename(ncfile)))[[1]]
if (length(m) != 3)
  stop("Cannot parse the product year range from: ", basename(ncfile))
FYnc <- as.integer(m[2]); LYnc <- as.integer(m[3])

# Annual product: one "period" per year, so the PHU matches each year's growing
# period exactly. PHU_SMOOTH_WINDOW (default 1) optionally widens the temperature
# averaging only.
SYs <- FYnc:LYnc; EYs <- FYnc:LYnc
smooth_window <- as.integer(Sys.getenv("PHU_SMOOTH_WINDOW", as.character(phu_smooth_window)))
cat(sprintf("PHU: %s %s %s_%s | years %d-%d | smooth_window=%d\n",
            gcm, scen, cro, irri, FYnc, LYnc, smooth_window))

ncdir  <- paste0(output_dir, "crop_calendars/ncdf/", gcm, "/", scen, "/")  # PHU .nc4 output
if (!dir.exists(ncdir)) dir.create(ncdir, recursive = TRUE)

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
    grid_df       = grid_df,
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
