#!/bin/bash

# Toolchain modules (nco + cdo) — single source of truth, see env.sh
source "$(dirname "$(readlink -f "$0")")/env.sh"
load_nco_cdo_env

# Deployment settings (PUBLISH_DIR, ...) — see settings.sh
source "$(dirname "$(readlink -f "$0")")/settings.sh"
BASE_DIR_ROOT=$PUBLISH_DIR

#for PERIOD in $periods; do
#BASE_DIR=$BASE_DIR_ROOT #/$PERIOD
GCMS="GFDL-ESM4"
SPECS="historical"

REF_DATE="1601-01-01,00:00:00,1year"
CALENDAR="standard"

for GCM in $GCMS; do

  GCM_LC=$(echo $GCM | tr '[:upper:]' '[:lower:]')
  echo $GCM_LC



for SPEC in $SPECS;do
  
  BASE_DIR=${BASE_DIR_ROOT}/$GCM/$SPEC
  echo $BASE_DIR

for FILE in $(find $BASE_DIR -maxdepth 3 -type f | sort );do

  echo "  "
  echo $FILE

  
#  if [ ! -f ${FILE}.tmp ]; then

#    continue

  if [ -z "$(ncdump -hs $FILE | grep missing_value)" ]; then
#    echo $FILE $(ncdump -hs $FILE | grep missing_value)
    ncatted -O -h -a missing_value,,o,f,1e+20 $FILE 2> /dev/null
#    exit
  fi

#  VARS=$(cdo -s showname $FILE)
#  for VAR in $VARS; do 
#    echo $VAR
#    ncatted -O -h -a _FillValue,$VAR,o,f,1e+20 $FILE 2> /dev/null
#  done
#    mv $FILE.tmp $FILE
#      exit
#  fi
done
done
done
