# Crop Calendar Pipeline — Deployment Hand-off

> **CURRENT STATE (annual sliding-window pipeline, package v0.2.0).** The pipeline has
> been rebuilt from the 10-year-step scheme to an **annual 30-yr sliding window** (see
> `METHODOLOGY_AND_CHANGES.md` in the working dir for the full design + change log).
> Most of the "Context"/"Deployment Plan" sections below describe the **superseded**
> window-scheme pipeline and are kept only for history. The current pipeline is:
>
> - **Config** `00_config.R`: `gcms` + per-GCM `scenarios` matrix; tunables
>   `clm_avg_years=30`, `clm_emit_step=1`, `pet_method="fao56"`, `phu_smooth_window=1`,
>   and the oscillation-suppression knobs `clm_smooth_window=15`, `cross_min_duration=5`,
>   `wet_window_eps=0.5`, `wet_window_decay=0.3`, `seas_eps=0.25`. Crops derived from
>   `crop_ls`. `enms`/`syears`/`eyears`/`ccal_years` removed.
> - **Climate** `settings.sh` `CLIMATE_DIR`: colon-separated **search list** of official
>   ISIMIP roots (3b primary+secondary, 3a obsclim/spinclim). Files discovered by
>   globbing; ensemble member/scenario/year-range read from file names. Cells = the
>   67420 GGCMI mask in `ggcmi_landcells.csv`.
> - **Stage 01 — use the SPLIT** (`01a_climatology_annual.R` → `01b_calendars_annual.R`).
>   `01a`: stream climate through the 30-yr ring, write per-year **smoothed** climatology
>   to disk (fork-free, ~25 GB; ~0.6 GB/yr cache). `01b`: load one year's climatology,
>   `calcCropCalendars` via `mclapply` → annual per-crop calendars (no ring at the fork →
>   no OOM, full 64 cores; ~20 min Maize / ~2.5 h 7 crops; `YEARS=`/`CROPS=` dev subsets).
>   The single-pass `01_compute_annual_calendars.R` still works but **must** run with
>   `R_GC_MEM_GROW=0` and ≤40 workers (else the ring + output held live across the 64-way
>   fork OOMs at ~340 GB — see METHODOLOGY §6). Future scenarios seed from historical
>   (`SAVE_SEED`). Submit: `01a_*.sh` then `01b_*.sh`.
> - **Oscillation suppression (v0.2.0)** — three independent, default-off fixes calibrated
>   on GFDL-ESM4 (per-cell mean year-to-year sowing change 1.94→0.55 d): daily-climatology
>   smoothing + sustained-crossing (`clm_smooth_window`/`cross_min_duration`, temperature
>   branch); distance-weighted wettest-window selection (`wet_window_eps`, PREC near-tie);
>   seasonality-threshold deadband (`seas_eps`, class flips). See METHODOLOGY §6.
> - **Wet-window refinement (later in v0.2.0)** — (a) fixed the wettest-window hysteresis
>   **anchor bug**: `prev_wet` was taken from crop #1's sowing day, wrong for a winter crop
>   (`Winter_Wheat`) in the combined run (caused a spurious 1850→1851 jump in 21k cells);
>   `wet_doy` is now computed crop-independently. (b) the `calcDoyWetMonth` distance weight
>   is now **Gaussian** with a `1−eps` floor and a `wet_window_decay` knob (faster
>   mid-distance decline → stickier for bimodal "chronic flippers"; same near/far behaviour
>   as the old linear kernel). **NB: any combined-run calendars/NetCDFs produced before this
>   fix are contaminated by the anchor bug and must be regenerated.** See NEWS / METHODOLOGY §6.
> - **Stage 02** `02_assemble_annual_ncdf.R` (replaces old 02 + **04–07**): writes the
>   DRS-compliant NetCDF in one pass (final names, 1601 time axis, ascending lat, fill
>   values, chunking, publish path). Submit: `02_assemble_annual_ncdf.sh`.
> - **Stage 03** `03_calc_phu_for_lpjml.R`: PHU per year (reads the DRS file).
> - Open items: regenerate stage-02 NetCDF from the v0.2.1 GFDL calendars (the on-disk
>   NetCDFs are pre-fix); re-run ssp245 `01b` + its stage-02 with v0.2.1; re-confirm the
>   ssp245 trend-tracking calibration with the Gaussian kernel; roll calibrated settings
>   across the full GCM × scenario matrix; header-parity check vs the official ISIMIP
>   reference; ISIMIP3a DRS soc/naming for GSWP3.

## Context (historical — superseded window scheme)

The `cropCalendars` R package was substantially revised on branch `fix-alg-vectorize-phu`
during development at `/p/projects/landuse/LPJmL_for_MAgPIE/cropCalendars/`. The whole
package + pipeline has since been **moved to its deployment home**:

```
/p/projects/macmit/users/heinke/crop_calendars/cropCalendars/
```

The remaining tasks are to **fix the hardcoded paths in the pipeline scripts**, install the
revised package, and run a test case (GFDL-ESM4 historical) end-to-end before deploying.

---

## What Was Changed in the Package (branch `fix-alg-vectorize-phu`)

All changes are committed; `git log --oneline` shows:

```
e1ef691  Vectorize generatePHUTserie_isimip3: replace cell loop with matrix ops
8cbb7bb  Fix namespace errors in generatePHUTserie_isimip3
a84770f  Optimize memory and runtime in PET calculation and pipeline script
174574a  01_calc_crop_calendars: read rsds, rlds, huss, sfcwind, ps; use FAO-56 PET
fa81e49  Use actual daily data for wet-season and threshold-crossing DOYs
5518588  calcPET_FAO56: remove tmax/tmin, use tas directly
1e1381f  calcPET_FAO56: use 24 h fluxes for LW and aerodynamic term
26ad35f  Improve PET: observed Rn in calcPET, add FAO-56 option
```

### Algorithm fixes

| Function | What changed |
|---|---|
| `calcPET_FAO56` | Removed `tmax`/`tmin`; uses daily mean `tas` + specific humidity `huss`; FAO-56 Penman-Monteith |
| `calcPET` | Added observed-Rn branch (when `swdown`/`lwdown` supplied); `pmax` for vectorization |
| `calcDoyWetMonth` | Now takes 365-value daily P/PET climatology; 120-day circular rolling window (eliminates up to 351-day interpolation errors) |
| `calcDoyCrossThreshold` | Now takes 365-value daily temperature; circular Dec→Jan wrap via `c(x[365], x[1:364])` |
| `calcMonthlyClimate` | FAO-56 PET path vectorized; returns `dtemp` and `dppet` (daily DOY climatologies) |
| `calcSowingDate` | Signature updated: `daily_ppet`, `daily_temp` replace monthly equivalents |
| `calcCropCalendars` | Passes `dtemp`/`dppet` from `calcMonthlyClimate` to `calcSowingDate` |
| `generatePHUTserie_isimip3` | Fully vectorized: `lin_idx` mapping, year-by-year climate accumulation, matrix PHU ops |

### New pipeline inputs (01_calc_crop_calendars.R)

Script now reads 7 climate variables per pixel: `tas`, `pr`, `rsds`, `rlds`, `huss`,
`sfcwind`, `ps` and passes them to `calcMonthlyClimate(pet_method="fao56", ...)`.

---

## Two Pipelines — Do Not Confuse Them

| | Package-based pipeline | Standalone pipeline |
|---|---|---|
| Location | `cropCalendars/utils/ggcmi_ph3/` | `crop_calendars/compute_sdate_hdate/` + `compute_phu/` |
| Uses package | Yes — `library(cropCalendars)` | No — sources own functions |
| Has algorithm fixes | YES — this is what was revised | No — separate older code |
| Deploy for | Production with revised algorithms | Separate project (ISIMIP3bv2) |

**We want the package-based pipeline.**

---

## Deployment Plan

### 1. Deployment directory — already in place

The package and its pipeline scripts have already been moved here. No copy step is needed:
```
Package repo:  /p/projects/macmit/users/heinke/crop_calendars/cropCalendars/
Pipeline dir:  /p/projects/macmit/users/heinke/crop_calendars/cropCalendars/utils/ggcmi_ph3/
```
Below, `<deployment_dir>` means the pipeline dir above.

### 2. Paths — now centralized in `settings.sh` (DONE)

Deployment paths and the SLURM account are no longer hardcoded in each script. They
live in a single file, **`utils/ggcmi_ph3/settings.sh`** (plain `KEY=VALUE`), which the
`.sh` job scripts `source` and `00_config.R` parses. **To move this pipeline again, edit
only `settings.sh`.**

Current values:

| Key | Value | Notes |
|---|---|---|
| `WD` | `…/heinke/crop_calendars/cropCalendars/utils/ggcmi_ph3` | pipeline dir (no trailing slash) |
| `ACCOUNT` | `landuse` | SLURM `-A`; `heinke` is **not** in `macmit`, so the old `-A macmit` was rejected. Available: `isimiplp` (default), `landuse`, `lpjml`, `magpie`, … (`sacctmgr -nP show assoc user=$USER format=Account`) |
| `OUTPUT_DIR` | `/p/projects/macmit/data/GGCMI/phase3/GGCMI_ph3_adaptation_cropping_calendars/` | write target; group-writable, heinke can write |
| `CLIMATE_DIR` | `/p/projects/macmit/data/GGCMI/phase3/input_land_only_v2/` | read-only input |
| `ISIMIP3B_PATH` | `/p/projects/lpjml/input/scenarios/ISIMIP3bv2/` | `.clm` climate for PHU step |
| `GRID_BIN` | `/p/projects/lpjml/input/historical/input_VERSION2/grid.bin` | LPJmL grid (land cells), read by `00_config.R` |
| `AGMIP_DIR` | `/p/projects/macmit/data/GGCMI/AgMIP.input/phase3/crop_calendar/` | AgMIP reference cal., read-only; `02` sets `ggdir`, read as a global by `generateCropCalTSerie_isimip3()` |
| `NCDF_DIR` | `${OUTPUT_DIR}crop_calendars/ncdf` | derived; ncdf tree from step 2/3 (source for `04`) |
| `PUBLISH_DIR` | `${OUTPUT_DIR}ISIMIP3b/InputData/socioeconomic/crop_calendar` | derived; published ISIMIP3b tree (`04` writes it; `05/06/07` fix it) |

How it wires together:
- `00_config.R` parses `settings.sh` → `output_dir`, `climate_dir`, `isimip3b.path`, `agmip_dir`.
- `01/02/03_*.R` use `work_dir <- getwd()` (sbatch passes `--chdir=$WD`; run interactively
  from the dir) — no path strings in the R files anymore. `02` uses `ggdir <- agmip_dir`.
- `01/02/03_*.sh` `source settings.sh` → `$WD` / `-A ${ACCOUNT}`.
- `04_move_and_rename.sh` `source settings.sh` → `BASE_DIR_ROOT=$NCDF_DIR`, `BASE_DIR_OUT=$PUBLISH_DIR`.
- `05/06/07_*.sh` `source settings.sh` → `BASE_DIR_ROOT=$PUBLISH_DIR`.

Toolchain modules are centralized in **`env.sh`** (separate from `settings.sh`):
`load_r_env` (piam → R 4.3.2 + packages) for `01/02/03_*.sh`, and `load_nco_cdo_env`
(`module load nco`/`cdo`) for `04–07`. piam does not provide nco/cdo, so they stay
separate. Each `.sh` sources `env.sh` and calls the relevant function.

`NCDF_DIR`/`PUBLISH_DIR` use bash expansion of `OUTPUT_DIR`; they are only read by the
`.sh` scripts. The R parser ignores them (it reads only the four literal keys above), so
keep the R-consumed keys as literal paths.

Note: `02`'s `ggdir <- agmip_dir` is **required** — `generateCropCalTSerie_isimip3()`
(R/generateCropCalTSerie_isimip3.R:56) reads `ggdir` as a free/global variable, not as a
function argument. Same global-injection pattern as `generatePHUTserie_isimip3`. Do not delete.

### 3. Install the package

The pipeline scripts call `library(cropCalendars)`. The repo is already on branch
`fix-alg-vectorize-phu`, so install the working tree directly:

```r
devtools::install("/p/projects/macmit/users/heinke/crop_calendars/cropCalendars")
```

Or build and install as a tarball if devtools is not available on the cluster nodes:
```bash
R CMD build /p/projects/macmit/users/heinke/crop_calendars/cropCalendars
R CMD INSTALL cropCalendars_*.tar.gz
```

Verify after install:
```r
library(cropCalendars)
packageVersion("cropCalendars")
# Check key functions exist
exists("calcPET_FAO56")
exists("calcDoyWetMonth")
```

### 4. Create output directory structure

`01_calc_crop_calendars.R` creates `dfout_dir` and `plot_dir` automatically, but the
root `output_dir` must exist:
```bash
mkdir -p <output_dir>/crop_calendars/DT
mkdir -p <output_dir>/crop_calendars/ncdf
```

### 5. Test run — GFDL-ESM4 historical, single crop, single year

Before submitting the full grid, run one job interactively to verify:

```bash
cd <deployment_dir>
Rscript --vanilla 01_calc_crop_calendars.R \
  GFDL-ESM4 historical Maize 1991 1 1
```

Args: `GCM SC CROP YEAR NNODES NTASKS`  
Year 1991 uses the 1961-1990 climate average (30-year window ending 1990).

Check:
- No R errors
- Output file created: `<output_dir>/crop_calendars/DT/historical/GFDL-ESM4/DT_output_crop_calendars_Maize_GFDL-ESM4_historical_1961_1990.Rdata`
- Sowing/harvest dates look spatially reasonable (not all NA, not all 1)

### 6. Full test — GFDL-ESM4 historical

Once single-crop test passes, submit the full historical run using the `.sh` scripts
(update `gcms` and `scens` arrays to only `GFDL-ESM4` and `historical`):

**Step 1** (`01_calc_crop_calendars.sh`):
- One job per crop × year combination
- Historical needs years: 1851, 1861, ..., 2011 (as in original config: `seq(1851, 2011, by=10)`)
- ~7 crops × 17 years = ~119 jobs, each using 120 CPUs, runtime ~2 h

**Step 2** (`02_generate_crop_cal_timeseries.sh`):
- After step 1 completes
- 15 crops × 2 irri = 30 jobs, ~15 min each

**Step 3** (`03_calc_phu_for_lpjml.sh`):
- After step 2 completes
- 15 crops × 2 irri = 30 jobs, ~1 h each
- Uses `generatePHUTserie_isimip3()` from the package (now vectorized)

---

## Key Technical Notes for the Next Session

### SLURM `launch_failed_requeued_held`

Jobs sometimes go to this state immediately after starting (node-side prolog failure).
They will NOT restart on their own. Release with:
```bash
scontrol release $(squeue -u heinke --state=PD | awk 'NR>1 {print $1}')
```
Watch for jobs that go RUNNING → COMPLETING within seconds of starting — that is the
failure, not normal completion.

### Stage 02 seasonality / harv-reason encoding (FIXED)

`generateCropCalTSerie_isimip3` previously encoded the `seasonality` and
`harv-reason` ncdf variables with `as.numeric(as.factor(per-pixel vector))`. Because
`as.factor` was applied **per pixel**, the integer code was alphabetical among only
the levels present at that pixel — so every temporally-constant pixel mapped to code 1
regardless of its actual type (≈94 % spurious `NO_SEASONALITY`; `harv-reason` collapsed
to `GPmin`). The standalone avoided this by storing those DT columns as factors with
**globally-fixed levels**. Fix: convert with fixed global level vectors
(`season_levels`, `harvreason_levels`) so codes are consistent across the grid:
`seasonality` 1=NoSeas 2=Prec 3=PrecTemp 4=Temp 5=TempPrec; `harv-reason`
1=GPmin(hd_first) 2=GPmed(hd_maxrp) 3=GPmax(hd_last) 4=Wstress(hd_wetseas)
5=Topt(hd_temp_base) 6=Thigh(hd_temp_opt). This changes ONLY those two variables;
sowing/harvest dates are unaffected (jump detection depends on change-points, not the
absolute code). After the fix the GFDL-ESM4/historical maize distribution matches the
standalone (validated). Any ncdf produced before this fix must be regenerated (stage 02
→ re-publish 04-07; stage 03 PHU is unaffected as it reads only plant-day/maty-day).

### Stage 07 CDO segfault — must use `cdo -L`

`07_fix_timeaxis.sh` calls CDO (2.4.4) `setreftime`/`settaxis`/`invertlat` on the
NetCDF4 files. The cluster's NetCDF4/HDF5 library is **not thread-safe**, so CDO's
default multi-threaded I/O **segfaults** (`cdi error (cdf_enddef): NetCDF: HDF error`,
then a 2 KB stub file). The script masked this because a trailing `find -delete`
returned exit 0. Both `cdo` calls now pass `-L` to serialise HDF5 access — this is
mandatory here. Symptom of the bug: published files keep the stage-06 time axis
(no `since 1601` reftime) and are ~2 KB instead of ~1.3 GB.

### Climate units (ISIMIP3b NetCDF)

| Variable | Unit in file | Conversion needed |
|---|---|---|
| `tas` | K | subtract 273.15 → °C |
| `pr` | kg/m²/s | multiply by 86400 → mm/day |
| `rsds`, `rlds` | W/m² | none |
| `huss` | kg/kg | none |
| `sfcwind` | m/s at 10 m | none |
| `ps` | Pa | none |

These conversions are already in `01_calc_crop_calendars.R` (`k2deg()` and `×86400`).

### FAO-56 PET requirement

`calcMonthlyClimate(pet_method="fao56")` requires `rsds`, `rlds`, `huss`, `sfcwind`,
`ps` in addition to `tas` and `pr`. All are read in `01_calc_crop_calendars.R`.
The old `pet_method="pt"` (Priestley-Taylor) only needs `tas`/`pr` — do not revert.

### `generatePHUTserie_isimip3` globals

The function reads several variables from the calling script's environment:
`grid_df`, `NCELLS`, `years`, `nyears`, `crop_ls`, `irri_ls`, `work_dir`, `LYs`,
`isimip3b.path`. These are all set in `00_config.R` and `03_calc_phu_for_lpjml.R`.
The vectorized implementation also uses `get.isimip.tas()` which reads from CLM binary
files at `isimip3b.path` — not from the NetCDF climate input.

### ncdf variable names

The crop calendar ncdf files use `"plant-day"` and `"maty-day"` as variable names
(with hyphens). `generatePHUTserie_isimip3` reads these correctly. Do not rename to
`"planting_day"` / `"harvest_day"`.

---

## Files Changed in the Package (for reference)

```
R/calcPET_FAO56.R
R/calcPET.R
R/calcDoyWetMonth.R
R/calcDoyCrossThreshold.R
R/calcMonthlyClimate.R
R/calcSowingDate.R
R/calcCropCalendars.R
R/generatePHUTserie_isimip3.R
utils/ggcmi_ph3/01_calc_crop_calendars.R
```

---

## What Has NOT Been Done Yet

- ~~Path updates in the deployment copies of the pipeline scripts~~ DONE — centralized in `settings.sh` (step 2)
- Package installation in the deployment environment
- Validation of outputs against the old pipeline (expected: different sowing dates in wet/dry tropics due to 120-day rolling window fix; different PHU values for winter wheat/rapeseed due to corrected FAO-56 PET)
- Running the full GCM × scenario matrix (5 GCMs × 5 scenarios) — do GFDL-ESM4 historical first

---

## Future cleanups (non-blocking)

- **Cell set / LPJmL grid usage.** `01a` now processes **all climate land cells** (the non-NA
  cells of the first `tas` file), not the LPJmL grid — so the crop-calendar product (01a→01b→02)
  covers every cell with climate data. The LPJmL grid (`grid_df`, read via `lpjmlkit::read_io`
  in `00_config.R`) is only *used* by stage 03 to write the `.clm` files. For ISIMIP3b the two
  masks differ by exactly one sub-Antarctic island cell (climate has 178.75,−49.25; LPJmL has
  178.75,−49.75), so 01a yields 67420 cells vs the old 67419. **(DONE)** the grid is now read
  via `lpjmlkit::read_io` only in stage 03, and the superseded legacy `01_calc_crop_calendars.{R,sh}`
  has been removed (01a/01b are validated bit-identical). Remaining: stage 03 should gap-fill the
  one LPJmL cell the climate lacks (currently NA in the `.clm`).

- **(FIXED) `plotMap_ggplot` / `plotMapCropCalendars` namespacing + bugs.** Now qualify all
  external calls (`ggplot2::*`, `scales::squish`, `RColorBrewer::brewer.pal`,
  `ggplot2::facet_grid`), fixed the `fil = landFill` typo (→ `fill =`), and modernised
  `aes_string`→`aes(.data[[…]])` and `size`→`linewidth`. `scales`/`RColorBrewer` added to
  DESCRIPTION Suggests. Verified by rendering a full DT to PDF. `01b` still wraps the plot call
  in `tryCatch` (best-effort; the DT is saved first), but no longer needs to attach packages.


- **Remove hidden global dependencies from the package functions.**
  - **(DONE) `generateCropCalTSerie_isimip3()`** now takes `FYs`/`LYs`/`ggdir`/`ncdir`/`csvdir`/
    `pldir` as explicit args (and honours `years_nc`, previously dead). `findGlobals` shows no
    caller-env free variables remain.
  - **(TODO) `generatePHUTserie_isimip3()` cluster** still reads `grid_df`, `crop_ls`, `irri_ls`,
    `work_dir`, `LYs`, `years`/`nyears` (derivable from `FYnc:LYnc`) from the caller env, plus the
    nested helpers `get.isimip.tas()` (`isimip3b.path`) and `read.climate.input()` (8 defaults:
    `NCELLS`, `RYEAR`, `FYEAR`, `LYEAR`, `HEADER`, `NBANDS`, `DTYPE`, `SCALAR` — most undefined
    anywhere, so every call must pass them). This is a deep 3-function cluster and **stage 03
    cannot be validated without a full run** (needs the stage-02 ncdf + `.clm` climate), so it was
    deliberately deferred: do it as its own focused task *with* a stage-03 run to confirm. Same
    applies to the stage-03 gap-fill of the one LPJmL cell the climate lacks.

- **(FIXED) `date_to_doy(skip_feb29=TRUE)` leap-year fold point.**
  `R/zz_dates.R` used `doy1 > 28`, which mis-folded at **Jan 29** (DOY bin 28 got 2 days and
  late-Jan/Feb DOYs were shifted one day in leap years). Changed to `doy1 > 59` ("after Feb 28",
  day-of-year 59): Feb 29 now folds onto DOY 59 (shared with Feb 28), Jan 29 → DOY 29, and
  Mar 1–Dec 31 align between leap and non-leap years. This **changes outputs** — the daily
  climatologies `dtemp`/`dppet` in leap years and hence some sowing / threshold-crossing dates
  — so any pre-fix `master` outputs differ here. `01a` calls the same fixed `date_to_doy`, so
  the split stays faithful to `calcMonthlyClimate`. `01a` also guards the input `time:calendar`
  (must be leap-aware). (`skip_feb29=FALSE` was never the issue: it would yield 366 bins with a
  leap-only-sparse bin 366 and conflate Feb 29 with Mar 1; `calcMonthlyClimate` hardcodes `TRUE`.)

- **(DONE) Cache the monthly climate; split stage 01 into preprocess + per-crop steps.**
  `01a_calc_monthly_climate.{R,sh}` computes & caches the crop-independent monthly climate
  once per GCM×scenario×window → `crop_calendars/monthly_climate/<scen>/<gcm>/
  monthly_climate_<gcm>_<scen>_<sy>_<ey>.Rdata`; `01b_calc_crop_calendars.{R,sh}` loads the
  cache, runs `calcCropCalendars` per pixel, writes the same `DT_output_*.Rdata` stage 02
  consumes. Removes the ~N_crops× monthly-climate recompute (the reason the fused `01` was much
  slower than the standalone pipeline at `/p/projects/landuse/LPJmL_for_MAgPIE/crop_calendars/
  scripts/`), keeping FAO-56.
  Design (matches the standalone): `01a` **streams over years** (one year of the full grid in
  memory at a time) and **vectorises FAO-56 PET over all cells** — so it needs no spatial
  chunking and runs **single-threaded** (one job per GCM×scenario×year; parallelise via a SLURM
  job array if I/O-bound). `01b` is a **single-threaded** pixel loop (no chunking; parallelism
  is one job per crop×year). Verified identical to the per-pixel `calcMonthlyClimate` path:
  accumulation max abs diff 2e-14 (monthly fields exact after `round(,5)`); cell extraction
  exact; crop-calendar round-trip `identical==TRUE`. Original `01_calc_crop_calendars.{R,sh}`
  kept as reference; retire once 01a/01b are confirmed at full grid scale. (Standalone uses
  Priestley-Taylor / 2 vars vs FAO-56 / 7 vars — a separate intentional difference; keep FAO-56.)
