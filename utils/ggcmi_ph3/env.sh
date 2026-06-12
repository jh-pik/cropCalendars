# ---------------------------------------------------------------------------- #
# Centralized toolchain / module setup for this pipeline.
#
# The pipeline uses two toolchains, each loaded by the stage that needs it:
#   - R 4.3.2 (+ netcdf-c) from the PIK 'piam' module set   -> stages 01-03
#   - NCO + CDO from the standard module set                -> stages 04-07
# piam does not provide nco/cdo, so each stage sources this file and calls ONLY
# the function it needs.
#
# Off the PIK cluster: replace the function bodies with whatever puts R (and its
# packages) / nco + cdo on PATH in your environment (e.g. conda, spack, apt).
# ---------------------------------------------------------------------------- #

# R + package toolchain (stages 01-03): provides R 4.3.2 and netcdf-c.
load_r_env() {
  source /p/system/modulefiles/defaults/piam/module_load_piam_fast
  # HDF5 file locking fails ("Permission denied" in nc_create) on this shared
  # filesystem when writing NetCDF4 from a login node; disable it so stage 02/03
  # work both interactively and under sbatch.
  export HDF5_USE_FILE_LOCKING=FALSE
}

# NetCDF post-processing toolchain (stages 04-07): provides nco and cdo.
load_nco_cdo_env() {
  module load nco
  module load cdo
}
