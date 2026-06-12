#!/bin/bash

# Toolchain modules (nco + cdo) — single source of truth, see env.sh
source "$(dirname "$(readlink -f "$0")")/env.sh"
load_nco_cdo_env


# Deployment settings (PUBLISH_DIR, ...) — see settings.sh
source "$(dirname "$(readlink -f "$0")")/settings.sh"
BASE_DIR_ROOT=$PUBLISH_DIR

# --fix_rec_dmn $VAR $FILE $FILE.tmp
GCMS="GFDL-ESM4"
SPECS="historical"

REF_DATE="1601-01-01,00:00:00,1year"
CALENDAR="standard"

for GCM in $GCMS; do

  GCM_LC=$(echo $GCM | tr '[:upper:]' '[:lower:]')
  echo $GCM_LC

for SPEC in $SPECS;do
  
  BASE_DIR=$BASE_DIR_ROOT/$GCM/$SPEC
  echo $BASE_DIR

for FILE in $(find $BASE_DIR -maxdepth 3 -type f -name '*.nc' | sort );do

  echo "  "
  echo $FILE
    
  FILENAME=$(basename $FILE)
  EXTENSION=${FILENAME##*.}
  BASENAME=$(basename $FILENAME .${EXTENSION})
  echo $BASENAME

  OUTDIR=$(dirname $FILE)
  FILENEW=$OUTDIR/$BASENAME.tmp
  echo $FILENEW

#  continue

  if [ ! -f ${FILENEW} ]; then

    nccopy -k4 -c "time/1,lat/360,lon/720" $FILE $FILENEW
    rm $FILE
    mv $FILENEW $FILE
#    [ -f ${FILE}.tmp ] && ncatted -O -h -a missing_value,,o,f,1e+20 ${FILE}.tmp
#    mv $FILE $FILENEW
#        exit
  fi
done
done
done




