#!/bin/bash

module purge
# Toolchain modules (nco + cdo) — single source of truth, see env.sh
source "$(dirname "$(readlink -f "$0")")/env.sh"
load_nco_cdo_env

# Deployment settings (PUBLISH_DIR, ...) — see settings.sh
source "$(dirname "$(readlink -f "$0")")/settings.sh"
BASE_DIR_ROOT=$PUBLISH_DIR

#for PERIOD in $periods; do
#BASE_DIR=$BASE_DIR_ROOT #/$PERIOD
GCMS="GFDL-ESM4 IPSL-CM6A-LR MPI-ESM1-2-HR MRI-ESM2-0 UKESM1-0-LL"
SPECS="ssp585soc-adapt ssp370soc-adapt ssp126soc-adapt historical"

REF_DATE="1601-01-01,00:00:00,1year"
CALENDAR="standard"

for GCM in $GCMS; do

  GCM_LC=$(echo $GCM | tr '[:upper:]' '[:lower:]')
  echo $GCM_LC

for SPEC in $SPECS;do
  
  BASE_DIR=$BASE_DIR_ROOT/$GCM/$SPEC
  echo $BASE_DIR

for FILE in $(find $BASE_DIR -maxdepth 3 -type f | sort );do

  echo "  "
  echo $FILE
    
  if [ ! -f ${FILE}.tmp ]; then

    START_YEAR=""
    START_YEAR=$(echo $FILE | awk -F"/" '{print $NF}' | awk -F"_" '{print $6}')
    TIME=$(echo ${START_YEAR}-01-01,00:00:00,1year)
    echo $START_YEAR

#    continue

    cdo -s --history -setreftime,$REF_DATE -settaxis,$TIME -setcalendar,$CALENDAR $FILE $FILE.tmp
#    [ -f ${FILE}.tmp ] && ncatted -O -h -a missing_value,,o,f,1e+20 ${FILE}.tmp
    # invert latidues and correct fill value as per request by Matthias
    cdo -f nc4c -z zip invertlat -setctomiss,NaNf $FILE.tmp $FILE
#        exit
  fi
done
  find "$BASE_DIR" -type f -name '*.tmp' -delete
done
done
