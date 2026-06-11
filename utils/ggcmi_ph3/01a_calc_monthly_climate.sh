#!/bin/bash

# Step 1a: compute & cache the average monthly climate (once per GCM x scenario x year).
# Heavy step (reads 7 daily climate vars + FAO-56 PET) — runs on a full node.

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

# sbatch settings
nnodes=1
ntasks=120  # 120 (not 128): R caps simultaneous connections at 128, makeCluster needs headroom

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

sbatch --nodes=${nnodes} --ntasks-per-node=${ntasks} --exclusive \
-t 02:00:00 -J monthly_clm -A ${ACCOUNT} --chdir=${wd} --qos=standby \
R -f 01a_calc_monthly_climate.R \
--args "${gcms[gc]}" "${scens[sc]}" "${years[yy]}" "${nnodes}" "${ntasks}"

    done # yy
  done # sc
done # gc
