# ---------------------------------------------------------------------------- #
# Deployment settings — single source of truth for this pipeline copy.
#
# Sourced by the *.sh job scripts (native bash) and parsed by 00_config.R for
# the R side. Use plain KEY=VALUE with no spaces around "=".
#
# Keys read by R (OUTPUT_DIR, CLIMATE_DIR, ISIMIP3B_PATH, AGMIP_DIR) must be
# LITERAL paths — the R parser does not expand shell variables. Bash-only keys
# further down (NCDF_DIR, PUBLISH_DIR) may reference earlier keys; they are only
# `source`d by bash, where the expansion happens.
#
# When moving this pipeline to a new location, edit THIS file only.
# ---------------------------------------------------------------------------- #

# Working directory: where these pipeline scripts live (no trailing slash).
WD=/p/projects/macmit/users/heinke/crop_calendars/cropCalendars/utils/ggcmi_ph3

# SLURM account for sbatch (-A).
ACCOUNT=landuse

# Output root: where results are written (keep trailing slash).
OUTPUT_DIR=/p/projects/macmit/users/heinke/crop_calendars/output/

# Climate inputs, read-only (keep trailing slash).
CLIMATE_DIR=/p/projects/macmit/data/GGCMI/phase3/input_land_only_v2/

# ISIMIP3b .clm climate path used by generatePHUTserie_isimip3 (keep trailing slash).
ISIMIP3B_PATH=/p/projects/lpjml/input/scenarios/ISIMIP3bv2/

# LPJmL grid file (defines the 67420 land cells / their lon-lat), read by 00_config.R.
GRID_BIN=/p/projects/lpjml/input/historical/input_VERSION2/grid.bin

# AgMIP reference crop-calendar input, read-only (keep trailing slash).
# 02_generate_crop_cal_timeseries.R sets ggdir <- agmip_dir; generateCropCalTSerie_isimip3()
# reads ggdir as a global (it is not a function argument), so this must stay set.
AGMIP_DIR=/p/projects/macmit/data/GGCMI/AgMIP.input/phase3/crop_calendar/

# --- Post-processing trees (04–07), derived from OUTPUT_DIR; bash-only ---------
# ncdf working tree written by step 2/3 (source for 04).
NCDF_DIR=${OUTPUT_DIR}crop_calendars/ncdf
# Published ISIMIP3b output tree: 04 renames into it; 05/06/07 fix files in it.
PUBLISH_DIR=${OUTPUT_DIR}ISIMIP3b/InputData/socioeconomic/crop_calendar
