#!/bin/bash

# Run time per job: 15 min

# Deployment settings (WD, ACCOUNT, ...) — single source of truth, see settings.sh
source "$(dirname "$(readlink -f "$0")")/settings.sh"
wd=$WD

# Toolchain modules (R + packages) — single source of truth, see env.sh
source "$(dirname "$(readlink -f "$0")")/env.sh"
load_r_env

# GCM, scenario, crops, irrigations
gcms=('GFDL-ESM4')
scens=('historical')
crops=('wwh' 'swh' 'mai' 'ri1' 'ri2' 'soy' 'mil' 'sor' 'pea' 'sgb' 'cas' 'rap' 'sun' 'nut' 'sgc')
irrigs=('rf' 'ir')

#gcms=('GFDL-ESM4')
#scens=('historical')
#crops=('mai')
#irrigs=('rf')

for gc in "${!gcms[@]}";do
  for sc in "${!scens[@]}";do
    for cr in "${!crops[@]}";do
      for ir in "${!irrigs[@]}";do

        echo "GCM: ${gcms[gc]} --- SCENARIO: ${scens[sc]} --- CROP: ${crops[cr]}"
        echo "------------------------------------------------------------------"

sbatch --ntasks=1 --cpus-per-task=4 -J nc_${gc}_${sc}_${cr}_${ir} -A ${ACCOUNT} \
-t 01:00:00  --chdir=${wd} --qos=standby  \
R -f 02_generate_crop_cal_timeseries.R \
--args "${gcms[gc]}" "${scens[sc]}" "${crops[cr]}" "${irrigs[ir]}"

      done
    done
    #sleep 5m
  done
done

# -o out_err/nc_${gcms[gc]}_${scens[sc]}_${crops[cr]}_${irrigs[ir]}_3b.out -e out_err/nc_${gcms[gc]}_${scens[sc]}_${crops[cr]}_${irrigs[ir]}_3b.err