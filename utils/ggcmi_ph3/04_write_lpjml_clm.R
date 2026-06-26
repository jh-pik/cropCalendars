# ---------------------------------------------------------------------------- #
# Step 04 (annual): write LPJmL CLM inputs (sdate, hdate, phu) from the pipeline's
# netCDF products. One CLM (version-2) binary per variable, 30 bands = 15 crops x 2
# irrigations (rainfed crops first, then irrigated), cell-major (all bands of a cell
# contiguous), int16, scalar 1.
#
# Sources (decadal -- the PHU and the dates LPJmL is driven with all change only per
# decade and are mutually consistent):
#   sdate <- planting_day-median   (stage-02 product)
#   hdate <- maturity_day-median   (stage-02 product)
#   phu   <- phu                   (stage-03 product, output/crop_calendars/ncdf)
#
# Args: GCM SCENARIO
# ---------------------------------------------------------------------------- #

rm(list = ls(all.names = TRUE))
suppressMessages(library(ncdf4))
suppressMessages(library(lpjmlkit))
stime <- Sys.time(); print(stime)

work_dir <- getwd()
source(file.path(work_dir, "00_config.R"))

if (cluster_job == TRUE) { options(echo = FALSE); args <- commandArgs(trailingOnly = TRUE) } else {
  args <- c("GFDL-ESM4", "ssp126")
}
print(args)
gcm <- args[1]; scen <- args[2]

# ------------------------------------ #
# Stage-02 product directory + soc token (mirror 02/03), and stage-03 PHU directory.
if (scen %in% isimip3a_scenarios) {
  soc_file <- if (scen == "counterclim") "countersoc" else "histsoc"
  cal_dir  <- paste0(output_dir, "ISIMIP3a/InputData/socioeconomic/crop_calendar/", soc_file, "/")
} else {
  soc_dir  <- if (scen == "historical") "historical" else paste0(scen, "soc-adapt")
  soc_file <- if (scen == "historical") "histsoc" else scen
  cal_dir  <- paste0(output_dir, "ISIMIP3b/InputData/socioeconomic/crop_calendar/", gcm, "/", soc_dir, "/")
}
phu_dir <- paste0(output_dir, "crop_calendars/ncdf/", gcm, "/", scen, "/")
clm_dir <- paste0(output_dir, "crop_calendars/clm/", gcm, "/", scen, "/")
if (!dir.exists(clm_dir)) dir.create(clm_dir, recursive = TRUE)

# Product year range from any stage-02 file (all crops share it).
gcm_lc  <- tolower(gcm)
irr_tok <- function(ir) if (ir == "ir") "firr" else "noirr"
hits <- Sys.glob(paste0(cal_dir, "ggcmi-crop-calendar_", gcm_lc, "_", soc_file, "_*_annual_*.nc"))
if (length(hits) == 0) stop("Stage-02 crop-calendar files not found (run stage 02 first): ", cal_dir)
m  <- regmatches(basename(hits[1]), regexec("_annual_([0-9]{4})_([0-9]{4})\\.nc$", basename(hits[1])))[[1]]
SY <- as.integer(m[2]); EY <- as.integer(m[3])
years <- SY:EY; NYEARS <- length(years)

# ------------------------------------ #
# 30 bands: 15 crops rainfed, then 15 crops irrigated (LPJmL band order).
crops  <- crop_ls[["ggcmi"]]
bands  <- rep(crops, 2)
irris  <- rep(c("rf", "ir"), each = length(crops))
NBANDS <- length(bands)

date_nc <- function(cro, ir) paste0(cal_dir, "ggcmi-crop-calendar_", gcm_lc, "_", soc_file,
                                    "_", cro, "-", irr_tok(ir), "_annual_", SY, "_", EY, ".nc")
phu_nc  <- function(cro, ir) paste0(phu_dir, cro, "_", ir, "_", gcm, "_", scen, "_",
                                    SY, "-", EY, "_ggcmi_ph3_rule_based_phu.nc4")
nc1_paths <- mapply(date_nc, bands, irris)   # dates (planting/maturity-median)
nc2_paths <- mapply(phu_nc,  bands, irris)   # phu
miss <- c(nc1_paths, nc2_paths)[!file.exists(c(nc1_paths, nc2_paths))]
if (length(miss) > 0) stop("Missing input(s) (run stages 02 & 03 first):\n  ", paste(miss, collapse = "\n  "))

# ------------------------------------ #
# LPJmL grid (the CLM cell order) and the 720x360 -> cell index mapping.
NCELLS  <- 67420L
grid_io <- suppressWarnings(read_io(grid_file, silent = TRUE))
grid_df <- data.frame(lon = round(as.numeric(grid_io$data[, 1, 1]), 2),
                      lat = round(as.numeric(grid_io$data[, 1, 2]), 2))
if (nrow(grid_df) != NCELLS) stop("Unexpected grid cell count: ", nrow(grid_df))

nc_tmp <- nc_open(nc1_paths[1]); lons <- ncvar_get(nc_tmp, "lon"); lats <- ncvar_get(nc_tmp, "lat")
pmask  <- !is.na(ncvar_get(nc_tmp, "planting_day", start = c(1, 1, 1), count = c(720, 360, 1)))
nc_close(nc_tmp)
lin_idx <- (match(grid_df$lat, lats) - 1L) * 720L + match(grid_df$lon, lons)
# Reconcile the LPJmL vs GGCMI grid mismatch (a 1-cell sub-antarctic-island artefact, as in
# stage 03): any LPJmL cell off the product grid inherits the nearest product land cell.
off <- which(!pmask[lin_idx])
if (length(off) > 0) {
  land     <- which(as.vector(pmask))
  land_lon <- lons[((land - 1L) %% 720L) + 1L]
  land_lat <- lats[((land - 1L) %/% 720L) + 1L]
  for (i in off) {
    dd <- (land_lon - grid_df$lon[i])^2 + (land_lat - grid_df$lat[i])^2
    lin_idx[i] <- land[which.min(dd)]
  }
  cat("Remapped", length(off), "LPJmL cell(s) off the product grid to nearest land cell\n")
}

# ------------------------------------ #
# CLM version-2 header writer (LPJmL binary: headername + 7 int32 + 2 float32).
fwriteheader2 <- function(con, headername, bands, firstyear, nyears,
                          ncells = NCELLS, resolution = 0.5, scalar = 1) {
  writeChar(headername, con, eos = NULL)
  writeBin(as.integer(c(2L, 1L, firstyear, nyears, 0L, ncells, bands)), con, size = 4, endian = .Platform$endian)
  writeBin(c(resolution, scalar), con, size = 4, endian = .Platform$endian)   # float32 CELLSIZE, SCALAR
}

FN <- function(v) paste0(clm_dir, v, "_", gcm, "_", scen, "_", SY, "_", EY,
                         "_ggcmi_ph3_rule_based_crop_calendar_30bands.clm")
sdfile <- file(FN("sdate"), "wb"); fwriteheader2(sdfile, "LPJSOWD", NBANDS, SY, NYEARS)
hdfile <- file(FN("hdate"), "wb"); fwriteheader2(hdfile, "LPJSOWD", NBANDS, SY, NYEARS)
hufile <- file(FN("phu"),   "wb"); fwriteheader2(hufile, "LPJmLHU", NBANDS, SY, NYEARS)

# Open all ncdfs once.
nc1_list <- lapply(nc1_paths, nc_open)
nc2_list <- lapply(nc2_paths, nc_open)

cat(sprintf("\n[%s] %s %s | %d bands | years %d-%d -> CLM (cell-major, int16)\n",
            format(Sys.time(), "%H:%M:%S"), gcm, scen, NBANDS, SY, EY))
t_all <- Sys.time()

# ------------------------------------ #
# Year by year: per band map 720x360 -> cells, write all bands of a cell contiguously.
for (yy in seq_len(NYEARS)) {
  if (yy %% 20L == 1L || yy == NYEARS)
    cat(sprintf("[%s]   year %d/%d (%d)  elapsed %.1f min\n", format(Sys.time(), "%H:%M:%S"),
                yy, NYEARS, years[yy], as.numeric(Sys.time() - t_all, units = "mins")))
  xsd <- matrix(0L, NCELLS, NBANDS); xhd <- matrix(0L, NCELLS, NBANDS); xph <- matrix(0L, NCELLS, NBANDS)
  for (bb in seq_len(NBANDS)) {
    sdate <- ncvar_get(nc1_list[[bb]], "planting_day-median", start = c(1, 1, yy), count = c(720, 360, 1))
    hdate <- ncvar_get(nc1_list[[bb]], "maturity_day-median", start = c(1, 1, yy), count = c(720, 360, 1))
    phu   <- ncvar_get(nc2_list[[bb]], "phu",                 start = c(1, 1, yy), count = c(720, 360, 1))
    xsd[, bb] <- sdate[lin_idx]; xhd[, bb] <- hdate[lin_idx]; xph[, bb] <- phu[lin_idx]
  }
  xsd[is.na(xsd)] <- 0L; xhd[is.na(xhd)] <- 0L; xph[is.na(xph)] <- 0L
  writeBin(as.integer(c(t(xsd))), sdfile, size = 2, endian = .Platform$endian)
  writeBin(as.integer(c(t(xhd))), hdfile, size = 2, endian = .Platform$endian)
  writeBin(as.integer(c(t(xph))), hufile, size = 2, endian = .Platform$endian)
}

for (bb in seq_len(NBANDS)) { nc_close(nc1_list[[bb]]); nc_close(nc2_list[[bb]]) }
close(sdfile); close(hdfile); close(hufile)
cat("\nwrote CLM files to:", clm_dir, "\n"); print(Sys.time() - stime)
