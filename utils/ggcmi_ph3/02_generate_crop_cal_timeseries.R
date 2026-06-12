# ---------------------------------------------------------------------------- #
# Create NCDF files with sdate and hdate (multi-year and crop-specific)
#  for GGCMI-phase3 adaptation runs

# Author:  Sara Minoli
# Email:   sara.minoli@pik-potsdam.de
# ---------------------------------------------------------------------------- #

rm(list = ls(all.names = TRUE))

stime <- Sys.time() # Track run-time
print(stime)

# ------------------------------------ #
# General settings
# Run from the pipeline dir: sbatch passes --chdir=$WD, interactive runs cd there.
work_dir <- getwd()
source(file.path(work_dir, "00_config.R"))

makeplot <- TRUE


# ------------------------------------ #
# Individual-run settings
if (cluster_job == TRUE) {
  # import argument from bash script
  options(echo = FALSE) # if you want see commands in output file
  args <- commandArgs(trailingOnly = TRUE)
} else {
  args <- c("GFDL-ESM4", "ssp126", "swh", "ir")
}
print(args)

# ------------------------------------ #
# Select variable, crop, model, year
gcm    <- args[1]
scen   <- args[2]
cro    <- args[3]
irri   <- args[4]

if (scen == "2015gs") {

  FYnc <- 2015 # first year ncdf output file
  LYnc <- 2100 # last year ncdf output file

  # Sowing and cultivar to change every 10 years
  SYs <- seq(1991, 2091, by = 10) # Start of the period
  EYs <- seq(2000, 2100, by = 10) # End of the period
  # Computing sowing and harvest dates based on preceding 30-years climate
  FYs <- rep(1981, length(SYs)) # First year DT file (output of main.R)
  LYs <- rep(2010, length(SYs)) # Last  year DT file (output of main.R)

  nhist <- length(LYs[LYs < 2015]) # number of historical time slices
  HYs <- c(rep("historical", nhist), rep(scen, length(SYs) - nhist))

} else if (scen == "historical") {

  # first and last year of ncdf files
  FYnc <- 1850
  LYnc <- 2014

  # Application periods and the 30-yr DT window (FYs-LYs) each uses. Matches the
  # standalone pipeline (.../landuse/.../crop_calendars) so the product starts at
  # 1850: the first four periods (SY 1850-1880) reuse the earliest window 1850-1879
  # as proxy; the last three use the reference-period windows. Pure historical.
  SYs <- c(1850, 1860, 1870, 1880, 1890, 1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1991, 2001, 2011)
  EYs <- c(1859, 1869, 1879, 1889, 1899, 1909, 1919, 1929, 1939, 1949, 1959, 1969, 1979, 1990, 2000, 2010, 2014)
  FYs <- c(1850, 1850, 1850, 1850, 1860, 1870, 1880, 1890, 1900, 1910, 1920, 1930, 1940, 1950, 1961, 1971, 1981)
  LYs <- c(1879, 1879, 1879, 1879, 1889, 1899, 1909, 1919, 1929, 1939, 1949, 1959, 1969, 1979, 1990, 2000, 2010)
  HYs <- rep(scen, length(SYs))

} else {

  FYnc <- 2015
  LYnc <- 2100

  # Sowing and cultivar to change every 10 years
  SYs <- seq(1991, 2091, by = 10) # Start of the period
  EYs <- seq(2000, 2100, by = 10) # End of the period
  # Computing sowing and harvest dates based on preceding 30-years climate
  FYs <- seq(1961, 2061, by = 10) # First year DT file (output of main.R)
  LYs <- seq(1990, 2090, by = 10) # Last  year DT file (output of main.R)

  nhist <- length(LYs[LYs < 2015]) # number of historical time slices
  HYs <- c(rep("historical", nhist), rep(scen, length(SYs) - nhist))

}

# time dimension ncdf files
YEARSnc  <- FYnc:LYnc
NYEARSns <- length(YEARSnc)

print(data.frame(SYs, EYs, FYs, LYs, HYs))


# PATHS ----
dtdir  <- paste0(output_dir, "crop_calendars/DT/")
ncdir  <- paste0(output_dir, "crop_calendars/ncdf/", gcm, "/", scen, "/")
csvdir <- paste0(ncdir, "csv/")
pldir  <- paste0(output_dir, "crop_calendars/plots/", gcm, "/", scen, "/")
ggdir  <- agmip_dir

if (!dir.exists(ncdir )) dir.create(ncdir,  recursive = T)
if (!dir.exists(pldir )) dir.create(pldir,  recursive = T)
if (!dir.exists(csvdir)) dir.create(csvdir, recursive = T)

# ------------------------------------------------------#
# Assemble crop calendar time series from individual DT files
# Filter discontinuities
# Save ncdf files

generateCropCalTSerie_isimip3(
    gcm        = gcm,
    scen       = scen,
    cro        = cro,
    irri       = irri,
    SYs        = SYs,
    EYs        = EYs,
    FYs        = FYs,
    LYs        = LYs,
    HYs        = HYs,
    years_nc   = YEARSnc,
    crop_ls    = crop_ls,
    irri_ls    = irri_ls,
    ggdir      = ggdir,
    output_dir = dtdir,
    ncdir      = ncdir,
    csvdir     = csvdir,
    pldir      = pldir,
    ncores     = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
    makeplot   = FALSE
)

# --------------------------#

print(summary(warnings()))


etime <- Sys.time()
cat("Run time:\n--------\n")
print(etime-stime)

cat("\n------\nDone !\n------\n")
