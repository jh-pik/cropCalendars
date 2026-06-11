#!/bin/bash

# Step 1a: compute & cache the average monthly climate (once per GCM x scenario x year).
# Year-streaming + cell-vectorised FAO-56: single-threaded, bounded memory (~one
# variable-year). One job per GCM x scenario x year; parallelise via job array if needed.

# Deployment settings (WD, ACCOUNT, ...) — single source of truth, see settings.sh
source "$(dirname "$(readlink -f "$0")")/settings.sh"
wd=$WD

# Toolchain modules (R + packages) — single source of truth, see env.sh
source "$(dirname "$(readlink -f "$0")")/env.sh"
load_r_env

# GCM, scenario (no crop — the monthly climate is crop-independent)
# gcms=('GFDL-ESM4' 'IPSL-CM6A-LR' 'MPI-ESM1-2-HR' 'MRI-ESM2-0' 'UKESM1-0-LL')
gcms=('GFDL-ESM4')
scens=('historical')

# MAIN
for gc in "${!gcms[@]}";do
  for sc in "${!scens[@]}";do

    # Select years for each scenario
    if [ ${scens[sc]} = 'picontrol' ]
    then
      years=($(seq 1601 10 2091))
    elif [ ${scens[sc]} = 'historical' ]
    then
      years=($(seq 1851 10 2021))
    else
      years=($(seq 2011 10 2091))
    fi

    for yy in "${!years[@]}";do

echo "GCM: ${gcms[gc]} --- SCENARIO: ${scens[sc]} --- YEAR: ${years[yy]}"

sbatch --ntasks=1 --cpus-per-task=1 --mem=8G \
-t 02:00:00 -J monthly_clm -A ${ACCOUNT} --chdir=${wd} --qos=standby \
R -f 01a_calc_monthly_climate.R \
--args "${gcms[gc]}" "${scens[sc]}" "${years[yy]}"

    done # yy
  done # sc
done # gc
