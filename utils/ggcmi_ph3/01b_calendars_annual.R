# ---------------------------------------------------------------------------- #
# Step 01b (annual): rule-based crop calendars from the cached per-year
# climatology written by 01a. For each year, loads that year's climatology
# (~0.6 GB), runs calcCropCalendars for every crop x cell via mclapply, and
# writes the per-crop annual calendars (same format as the old unified stage 01).
#
# Because the heavy ring buffer is NOT live here (only one year's climatology is
# in the parent at fork time), the 64-worker fork is cheap (~tens of GB, no OOM),
# and this stage is re-runnable in minutes on every sowing/harvest RULE change
# WITHOUT re-reading the raw climate.
#
# Args: GCM SCENARIO [NCORES]
# Env : YEARS subset for dev iteration, e.g. YEARS="2000:2014" (default: every
#       year for which 01a wrote a climatology file). NB: a sparse subset breaks
#       the per-cell wettest-window hysteresis continuity (only relevant when
#       wet_window_eps > 0).
# ---------------------------------------------------------------------------- #

rm(list = ls(all.names = TRUE))
suppressMessages(library(parallel))
starttime <- Sys.time(); print(starttime)

work_dir <- getwd()
source(file.path(work_dir, "00_config.R"))

if (cluster_job == TRUE) { options(echo = FALSE); args <- commandArgs(trailingOnly = TRUE) } else {
  args <- c("GFDL-ESM4", "historical", "8")
}
print(args)
gcm    <- args[1]
scen   <- args[2]
ncores <- if (length(args) >= 3) as.integer(args[3]) else
          as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

if (!exists("wet_window_eps"))        wet_window_eps        <- 0
if (!exists("wet_window_drift_gate")) wet_window_drift_gate <- 7L
if (!exists("cross_min_duration"))    cross_min_duration    <- 1L
if (!exists("seas_eps"))              seas_eps              <- 0

parse_years <- function(s) {
  if (is.null(s) || s == "") return(NULL)
  sort(unique(unlist(lapply(strsplit(s, ",")[[1]], function(p) {
    p <- trimws(p)
    if (grepl(":", p)) { ab <- as.integer(strsplit(p, ":")[[1]]); ab[1]:ab[2] } else as.integer(p)
  }))))
}
years_env <- parse_years(Sys.getenv("YEARS", ""))

clim_dir <- paste0(output_dir, "/crop_calendars/annual_climatology/", scen, "/", gcm, "/")
out_dir  <- paste0(output_dir, "/crop_calendars/annual/", scen, "/", gcm, "/")
if (!dir.exists(clim_dir)) stop("No climatology dir (run 01a first): ", clim_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Discover the cached climatology years; intersect with the YEARS subset.
cfiles <- Sys.glob(file.path(clim_dir, sprintf("climatology_%s_%s_*.Rdata", gcm, scen)))
cyears <- as.integer(sub(".*_(\\d{4})\\.Rdata$", "\\1", cfiles))
o <- order(cyears); cfiles <- cfiles[o]; cyears <- cyears[o]
if (!is.null(years_env)) { keep <- cyears %in% years_env; cfiles <- cfiles[keep]; cyears <- cyears[keep] }
if (length(cyears) == 0) stop("No climatology files match (run 01a; check YEARS).")
emit_years <- cyears
nE <- length(emit_years)
cat(sprintf("\n%s %s | %d climatology years (%d..%d) cores=%d | wet_eps=%g cross_min_dur=%d seas_eps=%g\n",
            gcm, scen, nE, min(emit_years), max(emit_years), ncores, wet_window_eps, cross_min_duration, seas_eps))

# Crops + pre-extracted parameters. CROPS env (rb_cal names, comma-separated, e.g.
# "Maize") restricts the crop set for fast dev/validation runs (default: all).
crops <- as.character(stats::na.omit(crop_ls[["rb_cal"]]))
crops_env <- trimws(strsplit(Sys.getenv("CROPS", ""), ",")[[1]])
crops_env <- crops_env[nzchar(crops_env)]
if (length(crops_env) > 0) {
  miss <- setdiff(crops_env, crops); if (length(miss)) stop("Unknown CROPS: ", paste(miss, collapse=", "))
  crops <- crops_env
}
cat("crops:", paste(crops, collapse=", "), "\n")
cpfile  <- system.file("extdata", "crop_parameters.csv", package = "cropCalendars", mustWork = TRUE)
cparams <- lapply(crops, function(cr) getCropParam(crops = cr, cropparam_file = cpfile))
names(cparams) <- crops

SEAS_LEV <- c("NO_SEASONALITY", "PREC", "PRECTEMP", "TEMP", "TEMPPREC")
HARV_LEV <- c("hd_first", "hd_maxrp", "hd_last", "hd_wetseas", "hd_temp_base", "hd_temp_opt")
FLDS     <- c("sow", "ss", "seas", "dflag", "maty_rf", "maty_ir",
              "gp_rf", "gp_ir", "hr_rf", "hr_ir")

# One year of calendars for all cells x crops (identical to the unified driver).
computeYear <- function(clim, prev_wet, prev_seas) {
  res <- mclapply(seq_len(NCELLS), function(j) {
    mcl <- list(mtemp = clim$mtemp[j, ], mprec = clim$mprec[j, ], mpet = clim$mpet[j, ],
                mppet = clim$mppet[j, ], mppet_diff = clim$mppet_diff[j, ],
                dtemp = clim$dtemp[j, ], dprec = clim$dprec[j, ], dpet = clim$dpet[j, ])
    M <- matrix(NA_real_, length(crops), length(FLDS), dimnames = list(NULL, FLDS))
    wd <- NA_integer_; st <- NA_character_
    for (ci in seq_along(crops)) {
      r <- calcCropCalendars(lon = land_lon[j], lat = land_lat[j], mclimate = mcl,
                             crop_parameters = cparams[[ci]],
                             prev_wet_doy = prev_wet[j], wet_window_eps = wet_window_eps,
                             wet_window_drift_gate = wet_window_drift_gate,
                             cross_min_duration = cross_min_duration,
                             prev_seas = prev_seas[j], seas_eps = seas_eps)
      if (ci == 1L) { wd <- attr(r, "wet_doy"); st <- attr(r, "seas_type") }
      M[ci, ] <- c(r$sowing_doy[1],
                   ifelse(r$sowing_season[1] == "winter", 1, 2),
                   match(r$seasonality_type[1], SEAS_LEV),
                   ifelse(r$sowing_month[1] == 0, 0, 1),
                   r$maturity_doy[1], r$maturity_doy[2],
                   r$growing_period[1], r$growing_period[2],
                   match(r$harvest_reason[1], HARV_LEV),
                   match(r$harvest_reason[2], HARV_LEV))
    }
    list(M = M, wd = wd, st = st)
  }, mc.cores = ncores)
  if (any(vapply(res, function(x) !is.list(x) || !is.matrix(x$M), logical(1))))
    stop("computeYear: a worker failed (retry with ncores = 1 to see the error).")
  list(arr = simplify2array(lapply(res, `[[`, "M")),
       wet = vapply(res, function(x) x$wd, integer(1)),
       seas = vapply(res, function(x) x$st, character(1)))
}

# ------------------------------------ #
# Load grid from the first climatology file to define the cell set.
e1 <- new.env(); load(cfiles[1], envir = e1)
grid_clm <- e1$grid_clm; land_lon <- grid_clm$lon; land_lat <- grid_clm$lat
NCELLS   <- nrow(grid_clm); rm(e1)
cat("GGCMI land cells:", NCELLS, "\n")

OUT <- lapply(crops, function(.) { z <- lapply(FLDS, function(.) matrix(NA_real_, NCELLS, nE)); names(z) <- FLDS; z })
names(OUT) <- crops
store <- function(e, arr) for (ci in seq_along(crops)) for (fi in seq_along(FLDS))
  OUT[[ci]][[FLDS[fi]]][, e] <<- arr[ci, fi, ]

prev_wet  <- rep(NA_integer_, NCELLS)
prev_seas <- rep(NA_character_, NCELLS)
for (e in seq_len(nE)) {
  t0 <- Sys.time()
  ce <- new.env(); load(cfiles[e], envir = ce); clim <- ce$clim; rm(ce)
  cy <- computeYear(clim, prev_wet, prev_seas); prev_wet <- cy$wet; prev_seas <- cy$seas
  store(e, cy$arr)
  rm(clim)
  cat(sprintf("  year %d: %.1fs  (rss %.1f GB)\n", emit_years[e],
              as.numeric(Sys.time() - t0, units = "secs"), as.numeric(gc()[2, 2]) / 1024))
}

# ------------------------------------ #
for (ci in seq_along(crops)) {
  cal  <- OUT[[ci]]
  crop <- crops[ci]
  fn <- paste0(out_dir, "annual_calendar_", crop, "_", gcm, "_", scen, "_",
               min(emit_years), "_", max(emit_years), ".Rdata")
  save(grid_clm, emit_years, cal, gcm, scen, crop, file = fn)
  cat("saved", fn, "\n")
}
endtime <- Sys.time(); print(endtime); print(endtime - starttime)
