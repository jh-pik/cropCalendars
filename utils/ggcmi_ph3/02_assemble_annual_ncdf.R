# ---------------------------------------------------------------------------- #
# Step 02 (annual): assemble the publication-ready GGCMI/ISIMIP3b crop-calendar
# NetCDF directly from the annual sliding-window calendars (output of
# 01_compute_annual_calendars.R).
#
# Writes the DRS-compliant file in ONE pass (final variable names, time axis with
# reference year 1601, calendar "standard", ASCENDING latitude, _FillValue +
# missing_value = 1e20, per-timestep chunking, compression, DRS filename + publish
# path). This replaces the old stage 02 + the 04-07 NCO/CDO chain.
#
# Header conventions matched to the official ISIMIP3b reference
#   .../crop_calendar/v.../<soc>/ggcmi-crop-calendar_<gcm>_<soc>_<var>-<irr>_annual_<y0>_<y1>.nc
# Latitude is ASCENDING (-89.75..89.75) as in the official product (the old
# pipeline's cdo invertlat wrongly produced descending latitude).
#
# Args: GCM SCENARIO CROP(ggcmi) IRRI(ggcmi)
# ---------------------------------------------------------------------------- #

rm(list = ls(all.names = TRUE))
suppressMessages(library(ncdf4))
stime <- Sys.time(); print(stime)

work_dir <- getwd()
source(file.path(work_dir, "00_config.R"))

if (cluster_job == TRUE) { options(echo = FALSE); args <- commandArgs(trailingOnly = TRUE) } else {
  args <- c("GFDL-ESM4", "historical", "mai", "ir")
}
print(args)
gcm <- args[1]; scen <- args[2]; cro <- args[3]; irri <- args[4]

lons <- seq(-179.75, 179.75, by = 0.5); nlon <- length(lons)
lats <- seq( -89.75,  89.75, by = 0.5); nlat <- length(lats)   # ASCENDING (S->N), as official

cr <- which(crop_ls[["ggcmi"]] == cro)
ir <- which(irri_ls[["ggcmi"]] == irri)
rb <- crop_ls[["rb_cal"]][cr]

# Load the annual calendar first; the product year range comes from the file (set
# by stage 01), NOT hardcoded -- so any scenario works (ssp 2015-2100, GSWP3
# 1901-2019, ...). Rule-based crop: its own file; default crop (rb = NA): Maize's.
load_annual <- function(crop_name) {
  f <- Sys.glob(paste0(output_dir, "crop_calendars/annual/", scen, "/", gcm, "/",
                       "annual_calendar_", crop_name, "_", gcm, "_", scen, "_*.Rdata"))[1]
  if (is.na(f)) stop("annual calendar not found: ", crop_name)
  e <- new.env(); load(f, envir = e); e
}
e <- load_annual(if (!is.na(rb)) rb else "Maize")
grid_clm <- e$grid_clm; ncell <- nrow(grid_clm)
years_nc <- as.integer(e$emit_years); nyears <- length(years_nc)
y0 <- min(years_nc); y1 <- max(years_nc)

# DRS naming / publish path (soc specifier: file token vs directory).
# TODO: ISIMIP3a observational forcings (GSWP3-W5E5 obsclim/counterclim) use a
# different DRS soc convention -- confirm before publishing those.
if (scen == "historical") { soc_file <- "histsoc"; soc_dir <- "historical" } else {
  soc_file <- scen; soc_dir <- paste0(scen, "soc-adapt") }
irr_tok  <- if (irri == "ir") "firr" else "noirr"
gcm_lc   <- tolower(gcm)
publish_dir <- paste0(output_dir, "ISIMIP3b/InputData/socioeconomic/crop_calendar")
outdir   <- paste0(publish_dir, "/", gcm, "/", soc_dir, "/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
ncfname  <- paste0(outdir, "ggcmi-crop-calendar_", gcm_lc, "_", soc_file, "_",
                   cro, "-", irr_tok, "_annual_", y0, "_", y1, ".nc")
cat(sprintf("\n%s %s | crop %s (rb=%s) irri %s -> %s\n", gcm, scen, cro, rb, irri, basename(ncfname)))

# ------------------------------------ #
# Observed GGCMI/AgMIP dates, oriented to our ASCENDING lat / lon grid.
nf <- nc_open(paste0(agmip_dir, cro, "_", irri, "_ggcmi_crop_calendar_phase3_v1.01.nc4"))
sdggcmi <- ncvar_get(nf, "planting_day"); hdggcmi <- ncvar_get(nf, "maturity_day")
gg_lat <- as.numeric(nf$dim[["lat"]]$vals); gg_lon <- as.numeric(nf$dim[["lon"]]$vals)
nc_close(nf)
if (gg_lat[1] > gg_lat[length(gg_lat)]) { sdggcmi <- sdggcmi[, ncol(sdggcmi):1]; hdggcmi <- hdggcmi[, ncol(hdggcmi):1] }
if (gg_lon[1] > 0) { stop("Unexpected GGCMI lon origin; expected -179.75..179.75") }

# ------------------------------------ #
# Per-cell annual fields (rule-based crop) or all-default (rb = NA).
ilon  <- match(round(grid_clm$lon, 2), round(lons, 2))
ilat  <- match(round(grid_clm$lat, 2), round(lats, 2))
lin   <- ilon + (ilat - 1L) * nlon
yones <- rep(1, nyears)

if (!is.na(rb)) {
  cal  <- e$cal
  pd <- cal$sow; md <- cal[[paste0("maty_", irri)]]; gp <- cal[[paste0("gp_", irri)]]
  seas <- cal$seas; hr <- cal[[paste0("hr_", irri)]]; ss <- cal$ss
  isdef <- cal$dflag == 0
} else {
  pd <- md <- gp <- seas <- hr <- ss <- matrix(NA_real_, ncell, nyears)
  isdef <- matrix(TRUE, ncell, nyears)
}

# GGCMI default-date replacement (broadcast per-cell observed date over years).
sdg <- outer(sdggcmi[cbind(ilon, ilat)], yones); hdg <- outer(hdggcmi[cbind(ilon, ilat)], yones)
repl <- which(isdef)
pd[repl] <- ifelse(!is.na(sdg[repl]), sdg[repl], 0)
md[repl] <- ifelse(!is.na(hdg[repl]), hdg[repl], 0)
gp[repl] <- ifelse(pd[repl] <= md[repl], md[repl] - pd[repl], md[repl] + 365 - pd[repl])

to_grid <- function(field) { A <- matrix(NA_real_, nlon * nlat, nyears); A[lin, ] <- field
  dim(A) <- c(nlon, nlat, nyears); A }
AR <- list(planting_day = to_grid(pd), maturity_day = to_grid(md), growing_period = to_grid(gp),
           seasonality = to_grid(seas), harvest_reason = to_grid(hr), planting_season = to_grid(ss))

# ------------------------------------ #
# Write the DRS-compliant NetCDF (one pass).
FV <- 1.0e20
londim <- ncdim_def("lon", "degrees_east",  lons, longname = "longitude")
latdim <- ncdim_def("lat", "degrees_north", lats, longname = "latitude")
timdim <- ncdim_def("time", "years since 1601-1-1 00:00:00", as.double(years_nc - 1601L),
                    unlim = TRUE, longname = "time", calendar = "standard")
diml <- list(londim, latdim, timdim); chk <- c(nlon, nlat, 1L)
def <- function(name, units, long) ncvar_def(name, units, diml, missval = FV, longname = long,
                                             prec = "single", chunksizes = chk, compression = 6)
vdef <- list(
  planting_day    = def("planting_day",    "day of year", "Rule-based sowing date (annual, 30-yr sliding window)"),
  maturity_day    = def("maturity_day",    "day of year", "Rule-based harvest date (annual, 30-yr sliding window)"),
  growing_period  = def("growing_period",  "days",        "Rule-based growing period duration (annual)"),
  seasonality     = def("seasonality",     "-", "Climate seasonality type, (1=No Seas; 2=Prec; 3=PrecTemp; 4=Temp; 5=TempPrec)"),
  harvest_reason  = def("harvest_reason",  "-", "Rule triggering harvest (1=GPmin; 2=GPmed; 3=GPmax; 4=Wstress; 5=Topt; 6=Thigh)"),
  planting_season = def("planting_season", "days", "Sowing season (1=Winter; 2=Spring)"))

ncout <- nc_create(ncfname, vdef, force_v4 = TRUE, verbose = FALSE)
for (v in names(vdef)) {
  ncvar_put(ncout, vdef[[v]], AR[[v]])
  ncatt_put(ncout, v, "missing_value", FV, prec = "float")
}
ncatt_put(ncout, "lon", "standard_name", "longitude"); ncatt_put(ncout, "lon", "axis", "X")
ncatt_put(ncout, "lat", "standard_name", "latitude");  ncatt_put(ncout, "lat", "axis", "Y")
ncatt_put(ncout, "time", "standard_name", "time");     ncatt_put(ncout, "time", "axis", "T")
ncatt_put(ncout, 0, "Conventions", "CF-1.6")
ncatt_put(ncout, 0, "Note", paste("The unit of the time dimension is calendar year and it refers to the",
          "plant-day variable. Mind that this does not always coincide with the maty-day year.",
          "If in a year (t), plant-day > maty-day, maty-day occurs the following year (t+1)."))
ncatt_put(ncout, 0, "Crop", paste0(cro, "_", irri))
ncatt_put(ncout, 0, "Institution", "Potsdam Institute for Climate Impact Research (PIK), Germany")
ncatt_put(ncout, 0, "History", paste0("Created ", format(Sys.time(), "%Y-%m-%d"),
          " by cropCalendars annual sliding-window pipeline (30-yr window slid by 1 year)."))
nc_close(ncout)
cat("saved", ncfname, "\n"); print(Sys.time() - stime)
