#!/bin/bash
# Stage 01b (annual): rule-based crop calendars from the 01a climatology cache for
# the GCM x scenario matrix. The heavy ring buffer is NOT live here (only one year's
# climatology at a time), so the 64-worker mclapply fork is cheap (no OOM) and the
# stage is re-runnable in minutes on every sowing/harvest RULE change.
#
# Env: CORES (default 64), QOS (default priority), YEARS subset (e.g. "2000:2014").

source "$(dirname "$(readlink -f "$0")")/settings.sh"
source "$(dirname "$(readlink -f "$0")")/env.sh"

gcms=('GFDL-ESM4' 'IPSL-CM6A-LR' 'MPI-ESM1-2-HR' 'MRI-ESM2-0' 'UKESM1-0-LL')
scens=('historical' 'ssp126' 'ssp370' 'ssp585')
CORES=${CORES:-64}
QOS=${QOS:-priority}
YEARS=${YEARS:-}
CROPS=${CROPS:-}                                    # rb_cal names, comma-sep (e.g. "Maize"); empty=all
WD="$(dirname "$(readlink -f "$0")")"
mkdir -p "${WD}/logs"

for g in "${gcms[@]}"; do
  for s in "${scens[@]}"; do
    sbatch --ntasks=1 --cpus-per-task=${CORES} -t 04:00:00 -J cc_cal_${g}_${s} \
      -A "${ACCOUNT}" --chdir="${WD}" --qos=${QOS} \
      -o logs/cc_cal_${g}_${s}-%j.out \
      --export=ALL,YEARS="${YEARS}",CROPS="${CROPS}" \
      --wrap="source ${WD}/env.sh; load_r_env; export R_GC_MEM_GROW=0; Rscript 01b_calendars_annual.R ${g} ${s} ${CORES}"
    echo "submitted cc_cal_${g}_${s} (${CORES} cores, ${QOS}, YEARS='${YEARS}')"
  done
done
