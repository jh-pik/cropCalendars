#!/bin/bash
# Stage 02 (annual): assemble the annual crop-calendar NetCDFs (run AFTER
# 01_compute_annual_calendars.R). One lightweight job per crop x irrigation.

source "$(dirname "$(readlink -f "$0")")/settings.sh"
source "$(dirname "$(readlink -f "$0")")/env.sh"

gcms=('GFDL-ESM4' 'IPSL-CM6A-LR' 'MPI-ESM1-2-HR' 'MRI-ESM2-0' 'UKESM1-0-LL')
scens=('historical' 'ssp126' 'ssp370' 'ssp585')
crops=('wwh' 'swh' 'mai' 'ri1' 'ri2' 'soy' 'mil' 'sor' 'pea' 'sgb' 'cas' 'rap' 'sun' 'nut' 'sgc')
irrigs=('rf' 'ir')

for gc in "${gcms[@]}"; do
  for sc in "${scens[@]}"; do
    for cr in "${crops[@]}"; do
      for ir in "${irrigs[@]}"; do
        sbatch --ntasks=1 --cpus-per-task=1 --mem=8G -t 00:20:00 -J nc_annual \
          -A "${ACCOUNT}" --chdir="${WD}" --qos=priority \
          -o logs/nc_annual_${gc}_${sc}_${cr}_${ir}-%j.out \
          --wrap="source $(dirname "$(readlink -f "$0")")/env.sh; load_r_env; Rscript 02_assemble_annual_ncdf.R ${gc} ${sc} ${cr} ${ir}"
      done
    done
  done
done
