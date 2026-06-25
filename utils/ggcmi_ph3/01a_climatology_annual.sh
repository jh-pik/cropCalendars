#!/bin/bash
# Stage 01a (annual): build the per-year sliding-window climatology cache for the
# GCM x scenario matrix. I/O + RAM heavy (the ring buffer), but single-threaded and
# fork-free, so memory is bounded (~25 GB) and it runs once. 01b then computes the
# calendars from this cache cheaply and re-runnably.
#
# Env: CORES (default 16), QOS (default priority), YEARS subset (e.g. "1990:2014").
#
# Per GCM the historical leg runs with SAVE_SEED=1 (banks the end-of-historical ring
# `seed_ring_<gcm>_w30_end<lastHistYear>.rds`), and each SSP leg is submitted with a SLURM
# dependency on that historical job so it LOADS the seed instead of re-reading the 30
# historical seed years. (`historical` must be first in `scens` for this to wire up.)

source "$(dirname "$(readlink -f "$0")")/settings.sh"
source "$(dirname "$(readlink -f "$0")")/env.sh"

gcms=('GFDL-ESM4' 'IPSL-CM6A-LR' 'MPI-ESM1-2-HR' 'MRI-ESM2-0' 'UKESM1-0-LL')
scens=('historical' 'ssp126' 'ssp370' 'ssp585')
CORES=${CORES:-16}                                 # ring only; cpus mainly buy memory budget
QOS=${QOS:-priority}
YEARS=${YEARS:-}
WD="$(dirname "$(readlink -f "$0")")"
mkdir -p "${WD}/logs"

for g in "${gcms[@]}"; do
  hist_jid=""
  for s in "${scens[@]}"; do
    SAVE=0; DEP=""
    if [ "$s" = historical ]; then SAVE=1
    elif [[ "$s" == ssp* && -n "$hist_jid" ]]; then DEP="--dependency=afterok:${hist_jid}"; fi
    jid=$(sbatch --parsable ${DEP} --ntasks=1 --cpus-per-task=${CORES} -t 08:00:00 -J cc_clm_${g}_${s} \
      -A "${ACCOUNT}" --chdir="${WD}" --qos=${QOS} \
      -o logs/cc_clm_${g}_${s}-%j.out \
      --export=ALL,EMIT_STEP=1,YEARS="${YEARS}",SAVE_SEED="${SAVE}" \
      --wrap="source ${WD}/env.sh; load_r_env; Rscript 01a_climatology_annual.R ${g} ${s}")
    [ "$s" = historical ] && hist_jid=$jid
    echo "submitted cc_clm_${g}_${s} -> ${jid} (SAVE_SEED=${SAVE} ${DEP})"
  done
done
