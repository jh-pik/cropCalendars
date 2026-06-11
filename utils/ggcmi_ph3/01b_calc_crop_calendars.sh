#!/bin/bash

# Step 1b: crop calendars from the cached monthly climate (run AFTER 01a completes).
# Light step (no climate I/O, just calcCropCalendars per pixel) — modest core count.

# Deployment settings (WD, ACCOUNT, ...) — single source of truth, see settings.sh
source "$(dirname "$(readlink -f "$0")")/settings.sh"
wd=$WD

# Toolchain modules (R + packages) — single source of truth, see env.sh
source "$(dirname "$(readlink -f "$0")")/env.sh"
load_r_env

# GCM, scenario, crops
# gcms=('GFDL-ESM4' 'IPSL-CM6A-LR' 'MPI-ESM1-2-HR' 'MRI-ESM2-0' 'UKESM1-0-LL')
gcms=('GFDL-ESM4')
scens=('historical')
crops=('Maize' 'Rice' 'Millet' 'Sorghum' 'Soybean' 'Spring_Wheat' 'Winter_Wheat')

# sbatch settings (light: 16 workers on one node, no I/O)
nnodes=1
ntasks=16

# MAIN
for gc in "${!gcms[@]}";do
  for sc in "${!scens[@]}";do

    # Select years for each scenario (must match 01a)
    if [ ${scens[sc]} = 'picontrol' ]
    then
      years=($(seq 1601 10 2091))
    elif [ ${scens[sc]} = 'historical' ]
    then
      years=($(seq 1851 10 2021))
    else
      years=($(seq 2011 10 2091))
    fi

    for cr in "${!crops[@]}";do
      for yy in "${!years[@]}";do

echo "GCM: ${gcms[gc]} --- SCENARIO: ${scens[sc]} --- CROP: ${crops[cr]} YEAR: ${years[yy]}"

sbatch --nodes=${nnodes} --ntasks=1 --cpus-per-task=${ntasks} \
-t 00:30:00 -J crop_cal -A ${ACCOUNT} --chdir=${wd} --qos=standby \
R -f 01b_calc_crop_calendars.R \
--args "${gcms[gc]}" "${scens[sc]}" "${crops[cr]}" "${years[yy]}" "${nnodes}" "${ntasks}"

      done # yy
    done # cr
  done # sc
done # gc
