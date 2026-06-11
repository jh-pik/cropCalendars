# cropCalendars <a href=''><img src='./inst/img/logo_0.1.0.png' align="right" height="139" /></a>

R-package for simulating crop calendars and their adaptation following the approaches from [Waha et al. (2012)](https://doi.org/10.1111/j.1466-8238.2011.00678.x) and [Minoli et al. (2019)](https://doi.org/10.1016/j.gloplacha.2018.12.013).

## Installation

```bash
git clone https://github.com/AgMIP-GGCMI/cropCalendars.git .
cd ..
R CMD build cropCalendars
# this generates e.g. cropCalendars_0.1.0.tar.gz
R CMD INSTALL [-l my/Rlib/path] cropCalendars_0.1.0.tar.gz
```

Alternative in R with `devtools`

```r
library(devtools)
devtools::build()
devtools::install()
```

## Usage

See examples in `./utils`

```r
library(cropCalendars)

# Main functions
climate <- calcMonthlyClimate()
crop_calendars <- calcCropCalendars()
phenological_heat_units <- calcPHU()
```

## What is portable, and what is not

This repository has **two layers**:

1. **The R package (`R/`, portable).** The core algorithm functions
   (`calcMonthlyClimate`, `calcCropCalendars`, `calcSowingDate`, `calcHarvestDate`,
   `calcPHU`, `calcPET` / `calcPET_FAO56`, …) are plain R that operate on in-memory
   vectors / data frames. Default crop parameters ship inside the package
   (`inst/extdata/crop_parameters.csv`, read via `system.file`). The only hard
   dependencies are CRAN packages (`ncdf4`, `data.table`). This layer runs on **any**
   machine with R — no cluster, scheduler, or fixed data paths required.

2. **The GGCMI / ISIMIP3b pipeline (`utils/ggcmi_ph3/`, PIK-specific).** A set of
   driver scripts that run the package over the global grid for the GGCMI phase-3
   protocol. This layer is tied to the PIK cluster (SLURM, environment modules) and to a
   fixed input-data layout. See below for how to run it, and what to change to run it
   elsewhere.

### Using just the package

```r
library(cropCalendars)

# Per-pixel daily/monthly climate -> crop calendar -> phenological heat units.
# Feed your own climate series (tas, pr, and for FAO-56 PET also rsds, rlds,
# huss, sfcwind, ps) as plain R objects; no files or paths are required.
clim <- calcMonthlyClimate(lat = lat, mtemp = tas, mprec = pr, ...)
ccal <- calcCropCalendars(lon = lon, lat = lat, ...)
phu  <- calcPHU(sdate = ccal$sdate, hdate = ccal$hdate, ...)
```

`?calcCropCalendars`, `?calcMonthlyClimate`, `?calcPHU` document the expected inputs.

## Running the GGCMI / ISIMIP3b pipeline

Scripts live in `utils/ggcmi_ph3/`. All deployment-specific values (paths, SLURM
account, derived output trees) are centralized in **`settings.sh`** — edit that one file
when moving the pipeline; the `.sh` jobs `source` it and `00_config.R` parses it.

Stages (run in order):

| Stage | Script | Does | Needs |
|---|---|---|---|
| config | `00_config.R` | sourced by 1a–3; reads `settings.sh`, crop lists | R + package |
| 1a | `01a_calc_monthly_climate.{R,sh}` | cache crop-independent monthly climate per GCM×scenario×window, for **all climate land cells** | R + package |
| 1b | `01b_calc_crop_calendars.{R,sh}` | crop calendars per crop from the 1a cache → `DT/*.Rdata` | R + package |
| 2 | `02_generate_crop_cal_timeseries.{R,sh}` | assemble time series → NetCDF | R + package, AgMIP ref. |
| 3 | `03_calc_phu_for_lpjml.{R,sh}` | PHUs for LPJmL → NetCDF (`.clm`) | R + package, `.clm` climate, **`lpjmlkit`** (LPJmL grid) |
| 4 | `04_move_and_rename.sh` | rename to ISIMIP3b DRS layout | **NCO** (`ncrename`) |
| 5–7 | `05_fix_missval.sh`, `06_fix_chunks.sh`, `07_fix_timeaxis.sh` | fix `missing_value`, chunking, time axis | **NCO / CDO** |

The R driver scripts (1a–3) are submitted with `sbatch` via the matching `.sh`, but the
`.R` files themselves run standalone — a single test case is two steps (1a caches the
monthly climate, 1b computes one crop from it):

```bash
cd utils/ggcmi_ph3
Rscript --vanilla 01a_calc_monthly_climate.R GFDL-ESM4 historical 1991        # -> cache
Rscript --vanilla 01b_calc_crop_calendars.R GFDL-ESM4 historical Maize 1991   # -> DT
# Run from this dir; work_dir = getwd().
```

Beyond the package, the pipeline R scripts also use: `abind`, `foreach`, `doParallel`,
`zoo`. NetCDF post-processing (stages 4–7) requires the **NCO** and **CDO** command-line
tools on `PATH`.

Toolchain setup is centralized in **`env.sh`**, which defines two functions:
`load_r_env` (R + packages; stages 01–03) and `load_nco_cdo_env` (NCO + CDO; stages
04–07). Each `.sh` sources `env.sh` and calls the one it needs. On the PIK cluster these
are two separate module sets — R comes from the `piam` module set, which does **not**
include `nco`/`cdo`, so those are loaded separately. To run elsewhere, edit the two
function bodies.

## Running outside the PIK cluster

The package layer needs nothing special. The **pipeline** assumes a PIK environment in a
few concrete places — to run it elsewhere, replace each:

- **Job scheduler.** `01/02/03_*.sh` submit with `sbatch` (`-A $ACCOUNT`, `--qos=standby`,
  `--chdir`). Without SLURM, run the `.R` files directly with `Rscript` (see above), or
  adapt the wrappers to your scheduler. Stages 4–7 are plain bash loops (no scheduler).
- **Environment modules.** Toolchain loading lives in `env.sh` (`load_r_env` /
  `load_nco_cdo_env`). Off-cluster, replace the function bodies with however you provide
  R and NCO/CDO (e.g. conda: `conda install -c conda-forge nco cdo`), or empty them if the
  tools are already on `PATH`.
- **Input datasets**, currently expected at the paths in `settings.sh` / `00_config.R`:
  - `CLIMATE_DIR` — ISIMIP3b daily climate NetCDF (`tas`, `pr`, `rsds`, `rlds`, `huss`,
    `sfcwind`, `ps`).
  - `ISIMIP3B_PATH` — ISIMIP3b `.clm` binary climate, read by `get.isimip.tas()` in the
    PHU stage. This reader assumes the LPJmL `.clm` format and file naming.
  - `AGMIP_DIR` — AgMIP reference crop calendars (NetCDF), used by stage 2.
  - `GRID_BIN` — LPJmL `grid.bin` defining the land cells, read **only by stage 03** (via
    `lpjmlkit::read_io`) to map the gridded product onto the LPJmL cell order for the `.clm`
    output. Stages 1a/1b/2 don't use it (they work on the climate land mask). Point it at your
    grid, or replace `grid_df` in 03 with your own `lon`/`lat` table.
- **Grid / resolution.** Stage 6 hardcodes `lat/360, lon/720` (0.5° global). Change for a
  different grid.

In short: the **science** is portable (the package); the **GGCMI driver** is an
ISIMIP3b/PIK harness that you re-point via `settings.sh` plus the data-layout assumptions
listed above.

## Contact

- Sara Minoli (sara.minoli@pik-potsdam.de)

```
Potsdam Institute for Climate Impact Research (PIK)
Member of the Leibniz Association
14412, Potsdam, Germany
```

## References

```
@article{Waha2012,
  title={Climate-driven simulation of global crop sowing dates},
  author={Waha, K and Van Bussel, LGJ and M{\"u}ller, C and Bondeau, Alberte},
  journal={Global Ecology and Biogeography},
  volume={21},
  number={2},
  pages={247--259},
  year={2012},
  doi={https://doi.org/10.1111/j.1466-8238.2011.00678.x},
  publisher={Wiley Online Library}
}

@article{Minoli2019,
  title={Modelling cropping periods of grain crops at the global scale},
  author={Minoli, Sara and Egli, Dennis B and Rolinski, Susanne and M{\"u}ller, Christoph},
  journal={Global and Planetary Change},
  volume={174},
  pages={35--46},
  year={2019},
  doi={https://doi.org/10.1016/j.gloplacha.2018.12.013},
  publisher={Elsevier}
}
```