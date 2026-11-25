# Bulk_forcing_WRF_to_CROCO. [**STATUS: UNDER TESTING**]
Python interface that converts WRF Lambert-grid 3-hourly output into the ERA5 intermediate format read by croco_tools (make_ERA5.m). It rotates winds to east/north, takes exact 3-hour means of radiation and precipitation (handling WRF’s running-total overflows and resets), and regrids to a regular lon/lat grid. 


# WRF → croco_tools interface for the experiments

`wrf2crocotools.py` converts the SHAMAL WRF 3-hourly output into the intermediate
"ERA5 croco format" read by croco_tools' `Aforc_ERA5/make_ERA5.m`. The usual MATLAB
step then builds the monthly bulk forcing files (`croco_blk_*_YyyyyMmm.nc`) on any
CROCO grid (GULF_16R, GULF_48R, …).

```
WRF wrfout_3hr_d01_*  --(wrf2crocotools.py)-->  <out-dir>/<VAR>_Y<yyyy>M<m>.nc
<out-dir>             --(make_ERA5.m, ERA5_dir = <out-dir>)-->  croco_blk_*_Y<yyyy>M<mm>.nc
```

---

## 1. Input

| Item | Value |
|---|---|
| WRF files | `<wrf_root>/<YYYY>/wrfout_3hr/wrfout_3hr_d01_<YYYY-MM-DD>_00:00:00` (8 records/day, 00–21 UTC) |
| Historical root | `/scratch/access/production/datasets/atmos/shamal/history` |
| Scenario (PGW) root | **to be set** (same `<year>/wrfout_3hr/` layout assumed) |
| WRF grid | Lambert conformal, 1023 × 895, DX = DY = 4.8 km, TRUELAT1/2 = 5/35°N, STAND_LON = 55°E, centre 20.5°N 55°E |
| WRF fields used | `T2`, `Q2`, `U10`, `V10`, `ACSWDNB`, `ACSWUPB`, `ACLWDNB`, `ACLWUPB` (+ `I_AC*` buckets), `AFWA_TOTPRECIP`, `SWDNB`, `SWUPB`, `LWDNB`, `LWUPB`, `ACSHFLSM`, `ACGHFLSM`, `ACETLSM`, `SH2O` |

The 3-hourly files contain **no XLAT/XLONG/LANDMASK**; see §4 for how this is handled.

## 2. Output

Directories:

| Experiment | Output directory |
|---|---|
| historical (CTRL) | `/scratch/rup9745/ROMS_gulf_16R/CROCO_FILES/SHAMAL/historical` |
| scenario (PGW) | `/scratch/rup9745/ROMS_gulf_16R/CROCO_FILES/scenario` |

One file per variable per month, named `<VAR>_Y<yyyy>M<m>.nc` (e.g. `T2M_Y2008M1.nc`),
in the same layout as croco_tools' `ERA5_convert.py`:

- dimensions `lon`, `lat`, `time` (unlimited); 1-D increasing `lon`/`lat`
- data `<VAR>(time, lat, lon)`, float32, `missing_value = 9999`
- `time` in days since `Yorig`-1-1 (Yorig = 1950, as in `crocotools_param.m`)

| Variable | Meaning | Units | Derived from WRF |
|---|---|---|---|
| `LSM` | land–sea mask (0 = sea, >0 = land) | 0–1 | land mask, see §4 (single record) |
| `T2M` | 2 m air temperature | K | `T2` |
| `Q` | 2 m specific humidity | kg kg⁻¹ | `Q2/(1+Q2)` (`Q2` is a mixing ratio) |
| `U10M`, `V10M` | 10 m wind, eastward / northward | m s⁻¹ | `U10`, `V10` rotated from grid axes to east/north |
| `SSR` | net surface shortwave (downward +) | W m⁻² | Δ(`ACSWDNB` − `ACSWUPB`) / Δt |
| `STRD` | downward surface longwave | W m⁻² | Δ`ACLWDNB` / Δt |
| `STR` | net surface longwave (downward +) | W m⁻² | Δ(`ACLWDNB` − `ACLWUPB`) / Δt |
| `TP` | total precipitation | kg m⁻² s⁻¹ | Δ`AFWA_TOTPRECIP` / Δt |

These are exactly the names, units and conversions `interp_ERA5.m` expects
(T2M in K → °C, Q → relative humidity via `qsat`, TP → cm day⁻¹, winds rotated
to the CROCO grid angle). `STR` is not used by the current `interp_ERA5.m`
(and `radlw` is not needed with `BULK_LW`), but is written for older versions.

**Target grid:** regular 0.04° lon/lat over 47.0–59.2°E, 22.8–31.9°N
(the CROCO Gulf domain 47.55–58.65°E, 23.35–31.35°N plus ~0.5° margin), 306 × 228 points.
Size ≈ 69 MB per variable per month uncompressed (≈ 0.55 GB/month for all variables;
less with the default zlib compression).

## 3. Time convention

- Each output record is the **mean over one 3-h WRF interval**, stamped at the
  **interval midpoint** (01:30, 04:30, …, 22:30 UTC): 8 records/day, 248 for a 31-day month.
- Radiative fluxes and precipitation are exact interval means from WRF's accumulated
  fields. T2M, Q and winds are the mean of the two end-point snapshots.
- Stamping at the midpoint avoids the 1.5 h phase lag of the diurnal shortwave cycle
  that end-of-interval labelling would give at 3-hourly resolution. Since every record
  is an interval mean, CROCO's linear time interpolation conserves daily totals.
- A month's last record (21:00–24:00 of the last day) uses the 00:00 record of the
  next month's first file (in the next year's directory for December). If that file
  does not exist (end of the dataset), the last interval is persisted and a `note`
  attribute is added to the files.
- Month overlaps (`itolap_era5`) are added by `make_ERA5.m` from the neighbouring
  months' files, as for ERA5.

### Accumulated fields

- **Bucket overflow** (`BUCKET_J = 1e9`): undone using the `I_AC*` counters when present.
  `I_ACLWUPB` is absent in the SHAMAL files, so for `ACLWUPB` a drop of the value just
  below the bucket size is recognised as an overflow.
- **Resets** (e.g. a restart at 1 January): a negative increment that is not an
  overflow is treated as a reset. That interval is filled with the mean of the
  instantaneous fluxes (`SWDNB`, `SWUPB`, `LWDNB`, `LWUPB`) at its two ends; for
  precipitation the previous interval is persisted.
- Counts of overflows/resets/fallbacks are printed for each month. If any accumulator
  decreases at more than 1 % of points/intervals, a **WARNING** is printed: that field
  may not be a running total and the month should be checked.

## 4. Horizontal grid and land mask

- **Coordinates:** WRF mass-point lat/lon are rebuilt from the Lambert projection
  attributes in the file (`TRUELAT1/2`, `STAND_LON`, `CEN_LAT/LON`, `DX/DY`;
  sphere R = 6370 km, as WRF/WPS). Domain centre = mass index ((nx−1)/2, (ny−1)/2).
  For the SHAMAL domain this gives corners (−0.51°N, 33.52°E) … (38.25°N, 82.73°E).
- **Interpolation:** bilinear in the native projection plane (exact WRF grid geometry)
  onto the regular target grid; only the WRF sub-window covering the target box
  (j 496:711, i 343:601) is read.
- **Wind rotation:** u_east = u·cosα + v·sinα, v_north = v·cosα − u·sinα,
  α = cone × (lon − STAND_LON) (wrf-python `uvmet` convention).
- **Land mask (default):** land where any of the Noah-MP accumulators `ACSHFLSM`,
  `ACGHFLSM`, `ACETLSM` is non-zero (they stay 0 over water), checked against
  `SH2O == 1` over water; the agreement is printed. `make_ERA5.m` treats every
  `LSM ≠ 0` point as land and fills it from the nearest sea values.
- **`--geo` (to do later):** when `geo_em.d01.nc` or `wrfinput_d01` is available, pass
  it with `--geo`. The script then compares the rebuilt coordinates with its
  XLAT/XLONG (stops if they differ by more than 0.02°) and uses its `LANDMASK`.

## 5. Requirements

Python 3 with `numpy` and `netCDF4` (the `roms_env` environment has them).
NCO and xESMF are not needed.

## 6. How to run

### Step 1 — dry run (login node, seconds)

Sets up the grid and land mask and reports them; writes nothing.

```bash
python wrf2crocotools.py --exp historical --start 2008-01 --end 2008-01 --dry-run
```

Check the `[grid]` lines: land fraction, agreement with the SH2O test (should be ~100 %),
WRF corners and the sea-point fraction of the target grid.

### Step 2 — convert (compute node)

```bash
# historical / CTRL
python wrf2crocotools.py --exp historical --start 2008-01 --end 2008-12 --nproc 4

# scenario / PGW (set the root in EXPERIMENTS or pass it here)
python wrf2crocotools.py --exp scenario --wrf-root /path/to/pgw/history \
                         --start 2008-01 --end 2008-12 --nproc 4
```

- Months run in parallel with `--nproc`; each worker needs well under 1 GB of memory.
- Existing complete months are skipped (use `--overwrite` to redo them).
- Files are written as `*.part` and renamed only when the month is complete, so a
  failed month leaves no partial files.
- A month fails (and is reported at the end, exit code 1) if any of its daily WRF
  files is missing or the time step is irregular.

Each month prints a summary of sea-point means, e.g.

```
[2008-01] 248 records; sea-point means: T2M 17.61 C, Q 15.73 g/kg, WSPD 5.88 m/s,
          SSR 269.3, STRD 380.7, STR -69.3 W/m2, TP 2.000 mm/day
```

Compare these with ERA5 for the same month as a first sanity check.

### Step 3 — bulk files with croco_tools (MATLAB)

Copy `Aforc_ERA5/make_ERA5.m` to `make_SHAMAL.m` and, right after the
`crocotools_param` line, replace the two `'_ERA5_'` prefix lines with:

```matlab
ERA5_dir   = '/scratch/rup9745/ROMS_gulf_16R/CROCO_FILES/SHAMAL/historical/';  % trailing slash needed
frc_prefix = [frc_prefix,'_SHAMAL_hist_'];
blk_prefix = [blk_prefix,'_SHAMAL_hist_'];
makefrc = 0;  makeblk = 1;                 % BULK_FLUX only needs the blk files
Ymin = 2008; Ymax = 2008; Mmin = 1; Mmax = 12;
```

For the scenario use the scenario directory and e.g. `'_SHAMAL_pgw_'`. The resulting
`croco_blk_SHAMAL_hist_Y2008M01.nc` replaces `croco_blk_ERA5_Y2008M01.nc` in the
`bulk_forcing:` entry of `croco.in`. No change to `cppdefs.h` is needed
(`BULK_FLUX`, `BULK_LW`, `ONLINE` undefined — same as the ERA5 hindcast).

## 7. Command-line options

| Option | Default | Meaning |
|---|---|---|
| `--exp` | `historical` | experiment preset (`historical`, `scenario`) from `EXPERIMENTS` |
| `--wrf-root` | from preset | WRF history root (contains `<YYYY>/wrfout_3hr/`) |
| `--out-dir` | from preset | output directory |
| `--start`, `--end` | required | first and last month, `YYYY-MM` (inclusive) |
| `--geo` | none | `geo_em.d01.nc` / `wrfinput_d01` for coordinate check and LANDMASK |
| `--box` | `47.0 59.2 22.8 31.9` | target lon min/max, lat min/max |
| `--res` | `0.04` | target grid step (deg) |
| `--yorig` | `1950` | time origin, must match `Yorig` in `crocotools_param.m` |
| `--month-fmt` | `d` | month in file names: `d` → `M1`, `02d` → `M01` |
| `--vars` | all 9 | subset of `LSM T2M Q TP SSR STRD STR U10M V10M` |
| `--nproc` | `1` | months processed in parallel |
| `--overwrite` | off | redo months that already exist |
| `--no-compress` | off | write without zlib compression |
| `--dry-run` | off | report grid/land mask only |

Defaults (paths, box, resolution, Yorig, month format) are set in the
"User settings" block at the top of the script.

## 8. Checks before production

1. **Month naming:** `make_ERA5.m` in croco_tools (Feb 2023 version) opens
   `T2M_Y2008M1.nc` (`num2str(M)`, no zero padding). Confirm against the existing ERA5
   intermediate files:
   `ls /scratch/rup9745/ROMS_gulf_16R/DATA/ERA5_ACCESS_GULF16 | head`
   If they are named `…M01.nc`, run with `--month-fmt 02d`.
2. **Scenario WRF path:** fill `EXPERIMENTS["scenario"]["wrf_root"]` or pass `--wrf-root`.
3. **Scenario output path:** used as given (`CROCO_FILES/scenario`, without `SHAMAL/`).
4. **Precipitation:** `AFWA_TOTPRECIP` is the only precipitation field in the 3-hourly
   files; watch for the reset WARNING in the first converted month.
5. **Land mask:** check the SH2O agreement in the dry run; replace with `--geo` when the
   WRF static file is available.

## 9. Validation

Tested on synthetic WRF files with the SHAMAL projection and analytic fields
(geometry from pyproj and earth-relative winds projected onto grid axes independently
of the script), covering a normal month boundary, bucket overflows with and without
counters, a 1 January accumulator reset and a missing next-month file:

| Check | Result |
|---|---|
| Rebuilt lat/lon vs pyproj (full 1023 × 895 grid) | max error 9e-14° |
| Rotation angle vs pyproj meridian convergence | max difference 2e-10° |
| U10M / V10M after rotation | max error 5e-5 m s⁻¹ |
| T2M / Q | max error 1.5e-4 K / 5e-9 kg kg⁻¹ |
| SSR / STRD / STR (incl. bucket overflows) | max error 0.05 / 4e-4 / 5e-4 W m⁻² |
| TP | max error 2e-9 kg m⁻² s⁻¹ |
| Reset interval (1 Jan) | STRD within 0.2 W m⁻² (instantaneous fallback); TP persisted |
| Time axis | first record 2008-01-01 01:30, step 3 h, 248 records in January |

Not yet tested on the real SHAMAL files.

## 10. Reference

Output format checked against croco_tools `Aforc_ERA5/make_ERA5.m`, `interp_ERA5.m`
and `ERA5_convert.py` (Feb 2023 version, GitHub mirror `SAEON/croco_tools`).
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
