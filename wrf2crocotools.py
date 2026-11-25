#!/usr/bin/env python3

import argparse
import calendar
import datetime as dt
import functools
import multiprocessing as mp
import os
import sys
import traceback

import numpy as np
from netCDF4 import Dataset, chartostring

# =============================================================================
# User settings (edit here or override on the command line)
# =============================================================================
EXPERIMENTS = {
    "historical": {
        "wrf_root": "/scratch/access/production/datasets/atmos/shamal/history",
        "out_dir": "/scratch/rup9745/ROMS_gulf_16R/CROCO_FILES/SHAMAL/historical",
    },
    "scenario": {
        "wrf_root": None,  # TODO: PGW wrfout root (same <year>/wrfout_3hr/ layout)
        "out_dir": "/scratch/rup9745/ROMS_gulf_16R/CROCO_FILES/scenario",
    },
}
# <wrf_root>/<YYYY>/wrfout_3hr/wrfout_3hr_d01_<YYYY-MM-DD>_00:00:00
WRF_SUBDIR = "wrfout_3hr"
WRF_PATTERN = "wrfout_3hr_d01_{date:%Y-%m-%d}_00:00:00"
WRF_DT_HOURS = 3

# Target regular grid: CROCO Gulf domain (47.55-58.65E, 23.35-31.35N) + ~0.5 deg
BOX = (47.0, 59.2, 22.8, 31.9)  # lonmin, lonmax, latmin, latmax
RES = 0.04                      # deg (WRF dx = 4.8 km)
YORIG = 1950                    # crocotools_param.m: Yorig
MONTH_FMT = "d"                 # make_ERA5.m opens T2M_Y2008M1.nc (num2str(M))
VARIABLES = ["LSM", "T2M", "Q", "TP", "SSR", "STRD", "STR", "U10M", "V10M"]

# =============================================================================
# Constants
# =============================================================================
EARTH_R = 6370000.0  # WRF/WPS sphere radius [m]
FILL = 9999.0
VAR_META = {  # units, long_name (as ERA5_convert.py)
    "LSM": ("(0-1)", "land_sea_mask"),
    "T2M": ("K", "2m_temperature"),
    "Q": ("kg kg-1", "specific_humidity_2m"),
    "TP": ("kg m-2 s-1", "total_precipitation"),
    "SSR": ("W m-2", "surface_net_solar_radiation"),
    "STRD": ("W m-2", "surface_thermal_radiation_downwards"),
    "STR": ("W m-2", "surface_net_thermal_radiation"),
    "U10M": ("m s-1", "10m_u_component_of_wind"),
    "V10M": ("m s-1", "10m_v_component_of_wind"),
}
# key: (accumulated variable, bucket counter, bucket-size attribute, kind)
ACC = {
    "SWDN": ("ACSWDNB", "I_ACSWDNB", "BUCKET_J", "J"),
    "SWUP": ("ACSWUPB", "I_ACSWUPB", "BUCKET_J", "J"),
    "LWDN": ("ACLWDNB", "I_ACLWDNB", "BUCKET_J", "J"),
    "LWUP": ("ACLWUPB", "I_ACLWUPB", "BUCKET_J", "J"),  # no I_ACLWUPB in SHAMAL
    "PREC": ("AFWA_TOTPRECIP", None, None, "mm"),
}
INST = {"SWDN": "SWDNB", "SWUP": "SWUPB", "LWDN": "LWDNB", "LWUP": "LWUPB"}
STATE = ["T2", "Q2", "U10", "V10"]
MAX_INC_J = 3.0e7                 # > any 3-h radiative increment [J m-2]
TOL = {"J": 1.0e4, "mm": 1.0e-2}  # float32 round-off allowed in differences
LAND_VARS = ("ACSHFLSM", "ACGHFLSM", "ACETLSM")


def log(tag, msg):
    print(f"[{tag}] {msg}", flush=True)


# =============================================================================
# WRF Lambert-conformal grid
# =============================================================================
class LambertGrid:
    """WRF Lambert mass grid rebuilt from the global attributes of a wrfout file."""

    def __init__(self, attrs, nx, ny):
        if int(attrs["MAP_PROJ"]) != 1:
            raise NotImplementedError("only MAP_PROJ = 1 (Lambert conformal)")
        t1, t2 = float(attrs["TRUELAT1"]), float(attrs["TRUELAT2"])
        if t1 <= 0 or t2 <= 0:
            raise NotImplementedError("southern-hemisphere Lambert not handled")
        self.stand_lon = float(attrs["STAND_LON"])
        self.dx, self.dy = float(attrs["DX"]), float(attrs["DY"])
        self.nx, self.ny = nx, ny
        p1, p2 = np.radians(t1), np.radians(t2)
        if abs(t1 - t2) > 0.1:  # cone factor exactly as WRF/WPS
            self.n = np.log(np.cos(p1) / np.cos(p2)) / np.log(
                np.tan(np.pi / 4 + p2 / 2) / np.tan(np.pi / 4 + p1 / 2))
        else:
            self.n = np.sin(p1)
        self.F = np.cos(p1) * np.tan(np.pi / 4 + p1 / 2) ** self.n / self.n
        # domain centre (CEN_LAT, CEN_LON) sits at mass index ((nx-1)/2, (ny-1)/2)
        self.ic, self.jc = (nx - 1) / 2.0, (ny - 1) / 2.0
        self.xc, self.yc = self._xy(float(attrs["CEN_LAT"]), float(attrs["CEN_LON"]))

    def _dlon(self, lon):
        return (np.asarray(lon, float) - self.stand_lon + 180.0) % 360.0 - 180.0

    def _xy(self, lat, lon):
        r = EARTH_R * self.F / np.tan(np.pi / 4 + np.radians(lat) / 2) ** self.n
        th = self.n * np.radians(self._dlon(lon))
        return r * np.sin(th), -r * np.cos(th)

    def ll_to_ij(self, lat, lon):
        """lat/lon -> fractional 0-based mass-grid indices (i: west_east, j: south_north)."""
        x, y = self._xy(lat, lon)
        return (x - self.xc) / self.dx + self.ic, (y - self.yc) / self.dy + self.jc

    def ij_to_ll(self, i, j):
        x = (np.asarray(i, float) - self.ic) * self.dx + self.xc
        y = (np.asarray(j, float) - self.jc) * self.dy + self.yc
        r = np.hypot(x, y)
        lon = self.stand_lon + np.degrees(np.arctan2(x, -y) / self.n)
        lat = np.degrees(2 * np.arctan((EARTH_R * self.F / r) ** (1 / self.n)) - np.pi / 2)
        return lat, (lon + 180.0) % 360.0 - 180.0

    def rotation(self, lon):
        """cos/sin of the grid-to-earth rotation angle.
        u_east = u*cos + v*sin ; v_north = v*cos - u*sin   (wrf-python uvmet)"""
        a = self.n * np.radians(self._dlon(lon))
        return np.cos(a), np.sin(a)


class Bilinear:
    """Bilinear interpolation from a native slab onto fixed target points."""

    def __init__(self, fi, fj):
        self.i0 = np.floor(fi).astype(np.intp)
        self.j0 = np.floor(fj).astype(np.intp)
        wx, wy = fi - self.i0, fj - self.j0
        self.w = ((1 - wx) * (1 - wy), wx * (1 - wy), (1 - wx) * wy, wx * wy)

    def __call__(self, a):
        i0, j0 = self.i0, self.j0
        w00, w10, w01, w11 = self.w
        return (a[..., j0, i0] * w00 + a[..., j0, i0 + 1] * w10
                + a[..., j0 + 1, i0] * w01 + a[..., j0 + 1, i0 + 1] * w11)


# =============================================================================
# File helpers
# =============================================================================
def wrf_path(root, date):
    return os.path.join(root, f"{date:%Y}", WRF_SUBDIR, WRF_PATTERN.format(date=date))


def out_name(out_dir, var, year, month, month_fmt):
    return os.path.join(out_dir, f"{var}_Y{year}M{format(month, month_fmt)}.nc")


def read_times(nc):
    return [dt.datetime.strptime(str(s), "%Y-%m-%d_%H:%M:%S")
            for s in chartostring(nc.variables["Times"][:])]


def read_wrf(path, slab, only_first=False):
    """Read the variables needed for the forcing on the native slab."""
    names = STATE + [v for k in ACC for v in ACC[k][:2] if v] + list(INST.values())
    with Dataset(path) as nc:
        nc.set_auto_mask(False)
        times = read_times(nc)
        tsl = slice(0, 1) if only_first else slice(None)
        data = {v: np.asarray(nc.variables[v][tsl, slab[0], slab[1]])
                for v in names if v in nc.variables}
        buckets = {a: float(nc.getncattr(a)) for a in ("BUCKET_J", "BUCKET_MM")
                   if a in nc.ncattrs()}
    missing = [v for v in STATE + [ACC[k][0] for k in ACC] if v not in data]
    if missing:
        raise KeyError(f"{path}: missing variables {missing}")
    return (times[:1] if only_first else times), data, buckets


# =============================================================================
# Target grid and land mask
# =============================================================================
def build_target(grid, box, res):
    lon = np.round(np.arange(box[0], box[1] + 0.5 * res, res), 6)
    lat = np.round(np.arange(box[2], box[3] + 0.5 * res, res), 6)
    LON, LAT = np.meshgrid(lon, lat)
    fi, fj = grid.ll_to_ij(LAT, LON)
    if fi.min() < 1 or fj.min() < 1 or fi.max() > grid.nx - 2 or fj.max() > grid.ny - 2:
        raise ValueError("target box extends beyond the WRF domain")
    i0, i1 = max(int(np.floor(fi.min())) - 2, 0), min(int(np.ceil(fi.max())) + 3, grid.nx)
    j0, j1 = max(int(np.floor(fj.min())) - 2, 0), min(int(np.ceil(fj.max())) + 3, grid.ny)
    jj, ii = np.mgrid[j0:j1, i0:i1]
    slat, slon = grid.ij_to_ll(ii, jj)
    cosa, sina = grid.rotation(slon)
    return {"lon": lon, "lat": lat, "slab": (slice(j0, j1), slice(i0, i1)),
            "interp": Bilinear(fi - i0, fj - j0), "cosa": cosa, "sina": sina}


def land_mask(first_file, grid, slab, geo=None):
    """Native-grid land mask (1 land, 0 water) on the slab."""
    if geo:
        with Dataset(geo) as g:
            g.set_auto_mask(False)
            vlat = "XLAT_M" if "XLAT_M" in g.variables else "XLAT"
            vlon = "XLONG_M" if "XLONG_M" in g.variables else "XLONG"
            glat = np.asarray(g.variables[vlat][0])
            glon = np.asarray(g.variables[vlon][0])
            if "LANDMASK" in g.variables:
                lm = np.asarray(g.variables["LANDMASK"][0])
            else:
                lm = (np.asarray(g.variables["XLAND"][0]) < 1.5).astype(float)
        if glat.shape != (grid.ny, grid.nx):
            raise ValueError(f"--geo grid {glat.shape} != wrfout grid {(grid.ny, grid.nx)}")
        jj, ii = np.mgrid[0:grid.ny, 0:grid.nx]
        alat, alon = grid.ij_to_ll(ii, jj)
        err = max(np.abs(alat - glat).max(),
                  np.abs((alon - glon + 180) % 360 - 180).max())
        log("grid", f"rebuilt lat/lon vs {os.path.basename(geo)}: max error {err:.2e} deg")
        if err > 0.02:
            raise ValueError("projection rebuild disagrees with --geo coordinates")
        return (lm[slab] > 0.5).astype(np.float64)

    with Dataset(first_file) as nc:
        nc.set_auto_mask(False)
        k = len(nc.dimensions["Time"]) - 1  # last record: most accumulation
        land, used = None, []
        for v in LAND_VARS:
            if v in nc.variables:
                a = np.asarray(nc.variables[v][k, slab[0], slab[1]]) != 0
                land = a if land is None else (land | a)
                used.append(v)
        water_sh2o = None
        if "SH2O" in nc.variables:
            water_sh2o = np.asarray(nc.variables["SH2O"][0, 0, slab[0], slab[1]]) >= 0.999
    if land is None:
        if water_sh2o is None:
            raise KeyError("no variable to derive a land mask: give --geo")
        land, used = ~water_sh2o, ["SH2O"]
    msg = f"land mask from {'|'.join(used)} != 0: land fraction {land.mean():.3f}"
    if water_sh2o is not None and used != ["SH2O"]:
        msg += f"; agreement with SH2O==1 water test {np.mean(land == ~water_sh2o) * 100:.2f}%"
    log("grid", msg)
    return land.astype(np.float64)


# =============================================================================
# Interval means
# =============================================================================
def acc_increment(ac, cnt, bucket, kind):
    """Increments of a WRF accumulator between consecutive records.
    Bucket overflow is undone (counter if present, else wrap test); a negative
    increment that is not a wrap is a reset -> NaN (filled later)."""
    tot = ac.astype(np.float64)
    wrap_ok = bucket is not None and bucket > 0
    if cnt is not None and wrap_ok:
        tot = tot + cnt.astype(np.float64) * bucket
        wrap_ok = False
    d = np.diff(tot, axis=0)
    d[(d < 0) & (d > -TOL[kind])] = 0.0
    neg = d < 0
    n_wrap = 0
    if wrap_ok:
        wrap = neg & (tot[:-1] > bucket - MAX_INC_J)
        d[wrap] += bucket
        n_wrap = int(wrap.sum())
        neg &= ~wrap
    d[neg] = np.nan
    return d, n_wrap, int(neg.sum())


def mid(a):
    return 0.5 * (a[:-1] + a[1:])


def interval_fields(times, data, buckets, tgt, carry, stats):
    """Native-grid interval means for consecutive records in `data`."""
    dts = np.array([(b - a).total_seconds() for a, b in zip(times[:-1], times[1:])])
    if np.any(np.abs(dts - WRF_DT_HOURS * 3600) > 1):
        bad = [f"{a}->{b}" for a, b, d in zip(times[:-1], times[1:], dts)
               if abs(d - WRF_DT_HOURS * 3600) > 1]
        raise ValueError(f"irregular WRF time step: {bad[:5]}")
    dsec = dts[:, None, None]

    inc = {}
    for key, (vac, vcnt, battr, kind) in ACC.items():
        bucket = buckets.get(battr) if battr else None
        cnt = data.get(vcnt) if vcnt else None
        inc[key], nw, nr = acc_increment(data[vac], cnt, bucket, kind)
        stats["wrap"] += nw
        stats["reset_" + key] += nr
    stats["points"] += inc["PREC"].size

    def inst_mean(expr):
        try:
            return mid(expr())
        except KeyError:
            return None

    out = {}
    # shortwave / longwave: fall back to instantaneous end-point means on resets
    pairs = {
        "SSR": (lambda: (inc["SWDN"] - inc["SWUP"]),
                lambda: data[INST["SWDN"]].astype(np.float64) - data[INST["SWUP"]]),
        "STRD": (lambda: inc["LWDN"], lambda: data[INST["LWDN"]].astype(np.float64)),
        "STR": (lambda: (inc["LWDN"] - inc["LWUP"]),
                lambda: data[INST["LWDN"]].astype(np.float64) - data[INST["LWUP"]]),
    }
    for name, (acc_fn, inst_fn) in pairs.items():
        f = acc_fn() / dsec
        bad = ~np.isfinite(f)
        if bad.any():
            alt = inst_mean(inst_fn)
            if alt is None:
                alt = np.zeros_like(f)
            f[bad] = alt[bad]
            stats["fallback_" + name] += int(bad.sum())
        out[name] = f

    # precipitation: persistence of the previous interval on resets
    tp = inc["PREC"] / dsec
    for k in range(tp.shape[0]):
        bad = ~np.isfinite(tp[k])
        if bad.any():
            prev = carry.get("TP")
            tp[k][bad] = prev[bad] if prev is not None else 0.0
            stats["fallback_TP"] += int(bad.sum())
        carry["TP"] = tp[k].copy()
    out["TP"] = tp

    out["SSR"] = np.maximum(out["SSR"], 0.0)
    out["STRD"] = np.maximum(out["STRD"], 0.0)
    out["TP"] = np.maximum(out["TP"], 0.0)

    q2 = data["Q2"].astype(np.float64)
    out["Q"] = mid(q2 / (1.0 + q2))
    out["T2M"] = mid(data["T2"].astype(np.float64))
    u, v = mid(data["U10"].astype(np.float64)), mid(data["V10"].astype(np.float64))
    out["U10M"] = u * tgt["cosa"] + v * tgt["sina"]
    out["V10M"] = v * tgt["cosa"] - u * tgt["sina"]

    tmid = [a + (b - a) / 2 for a, b in zip(times[:-1], times[1:])]
    return tmid, out


# =============================================================================
# Output
# =============================================================================
class MonthWriter:
    def __init__(self, ctx, year, month, history):
        self.ctx, self.n, self.paths, self.files = ctx, 0, {}, {}
        lon, lat = ctx["tgt"]["lon"], ctx["tgt"]["lat"]
        for var in ctx["variables"]:
            final = out_name(ctx["out_dir"], var, year, month, ctx["month_fmt"])
            part = final + ".part"
            nc = Dataset(part, "w", format="NETCDF4")
            nc.createDimension("lon", len(lon))
            nc.createDimension("lat", len(lat))
            nc.createDimension("time", None)
            vlon = nc.createVariable("lon", "f4", ("lon",))
            vlat = nc.createVariable("lat", "f4", ("lat",))
            vt = nc.createVariable("time", "f8", ("time",))
            vlon.long_name, vlon.units = "longitude of RHO-points", "degree_east"
            vlat.long_name, vlat.units = "latitude of RHO-points", "degree_north"
            vt.long_name, vt.units = "Time", f"days since {ctx['yorig']}-1-1"
            vd = nc.createVariable(var, "f4", ("time", "lat", "lon"),
                                   zlib=ctx["compress"], complevel=1,
                                   chunksizes=(1, len(lat), len(lon)))
            vd.missing_value = np.float32(FILL)
            vd.units, vd.long_name = VAR_META[var]
            vlon[:], vlat[:] = lon, lat
            nc.source = "SHAMAL WRF " + ctx["wrf_root"]
            nc.history = history
            nc.time_convention = ("3-h interval means stamped at interval midpoints; "
                                  "fluxes from WRF accumulated fields")
            self.paths[var], self.files[var] = (part, final), nc

    def append(self, tdays, fields):
        n = len(tdays)
        for var, nc in self.files.items():
            if var == "LSM":
                if self.n == 0:  # single record
                    nc.variables["time"][0] = tdays[0]
                    nc.variables[var][0] = fields[var]
                continue
            nc.variables["time"][self.n:self.n + n] = tdays
            nc.variables[var][self.n:self.n + n] = fields[var]
        self.n += n

    def close(self, ok, note=None):
        for var, nc in self.files.items():
            if note:
                nc.note = note
            nc.close()
            part, final = self.paths[var]
            if ok:
                os.replace(part, final)
            elif os.path.exists(part):
                os.remove(part)


# =============================================================================
# One month
# =============================================================================
def days_since(t, yorig):
    return (t - dt.datetime(yorig, 1, 1)).total_seconds() / 86400.0


def process_month(ym, ctx):
    year, month = ym
    tag = f"{year}-{month:02d}"
    finals = [out_name(ctx["out_dir"], v, year, month, ctx["month_fmt"]) for v in ctx["variables"]]
    if not ctx["overwrite"] and all(os.path.exists(f) for f in finals):
        log(tag, "outputs exist, skipped (use --overwrite)")
        return tag, True

    ndays = calendar.monthrange(year, month)[1]
    t_start = dt.datetime(year, month, 1)
    t_end = t_start + dt.timedelta(days=ndays)
    day_files = [wrf_path(ctx["wrf_root"], t_start + dt.timedelta(days=d)) for d in range(ndays)]
    missing = [f for f in day_files if not os.path.exists(f)]
    if missing:
        log(tag, f"ERROR: {len(missing)} WRF file(s) missing, e.g. {missing[:3]}")
        return tag, False
    next_file = wrf_path(ctx["wrf_root"], t_end)

    tgt, interp = ctx["tgt"], ctx["tgt"]["interp"]
    lsm = ctx["lsm"]
    sea = lsm == 0
    stats = {k: 0 for k in ["points", "wrap", "fallback_SSR", "fallback_STRD", "fallback_STR",
                            "fallback_TP"] + ["reset_" + k for k in ACC]}
    nsea = max(int(sea.sum()), 1)
    sums = {"WSPD": 0.0}  # sea-point sums for the end-of-month summary
    nrec = 0
    note = None
    history = f"{dt.datetime.now():%Y-%m-%d %H:%M} wrf2crocotools.py"
    writer = MonthWriter(ctx, year, month, history)
    ok = False
    try:
        prev_t, prev_d, carry, last = None, None, {}, None
        files = [(f, False) for f in day_files] + [(next_file, True)]
        for path, is_next in files:
            if is_next and prev_t[-1] >= t_end:
                break
            if is_next and not os.path.exists(path):
                note = f"next-month file {os.path.basename(path)} missing: last interval persisted"
                log(tag, "WARNING: " + note)
                break
            times, data, buckets = read_wrf(path, tgt["slab"], only_first=is_next)
            keep = [k for k, t in enumerate(times)
                    if t_start <= t <= t_end and (prev_t is None or t > prev_t[-1])]
            if not keep:
                continue
            times = [times[k] for k in keep]
            data = {v: a[keep] for v, a in data.items()}
            if prev_t is None:
                if times[0] != t_start:
                    raise ValueError(f"first record is {times[0]}, expected {t_start}")
            else:
                times = prev_t + times
                data = {v: np.concatenate([prev_d[v], a]) for v, a in data.items()}
            prev_t, prev_d = times[-1:], {v: a[-1:] for v, a in data.items()}
            if len(times) < 2:
                continue
            tmid, nat = interval_fields(times, data, buckets, tgt, carry, stats)
            fields = {v: interp(nat[v]).astype(np.float32) for v in nat}
            fields["LSM"] = lsm.astype(np.float32)
            for v, a in fields.items():
                if v != "LSM" and not np.all(np.isfinite(a)):
                    raise ValueError(f"non-finite values in {v}")
            writer.append([days_since(t, ctx["yorig"]) for t in tmid], fields)
            last = (tmid[-1], {v: a[-1:] for v, a in fields.items() if v != "LSM"})
            for v, arr in fields.items():
                if v != "LSM":
                    sums[v] = sums.get(v, 0.0) + float(np.sum(arr[:, sea])) / nsea
            sums["WSPD"] += float(np.sum(np.hypot(fields["U10M"], fields["V10M"])[:, sea])) / nsea
            nrec += len(tmid)

        if prev_t[-1] < t_end:  # persist the final interval up to month end
            t_last, f_last = last
            f_last["LSM"] = lsm.astype(np.float32)
            while t_last + dt.timedelta(hours=WRF_DT_HOURS) < t_end:
                t_last += dt.timedelta(hours=WRF_DT_HOURS)
                writer.append([days_since(t_last, ctx["yorig"])], f_last)

        expected = ndays * 24 // WRF_DT_HOURS
        if writer.n != expected:
            raise ValueError(f"wrote {writer.n} records, expected {expected}")
        ok = True
    except Exception as exc:  # noqa: BLE001
        log(tag, f"ERROR: {exc}")
        traceback.print_exc()
    finally:
        writer.close(ok, note)
    if not ok:
        return tag, False

    m = {v: s / max(nrec, 1) for v, s in sums.items()}
    log(tag, f"{writer.n} records; sea-point means: T2M {m['T2M'] - 273.15:.2f} C, "
        f"Q {m['Q'] * 1e3:.2f} g/kg, WSPD {m['WSPD']:.2f} m/s, SSR {m['SSR']:.1f}, "
        f"STRD {m['STRD']:.1f}, STR {m['STR']:.1f} W/m2, TP {m['TP'] * 86400:.3f} mm/day")
    flags = {k: v for k, v in stats.items() if v and k != "points"}
    if flags:
        log(tag, f"accumulator events (grid points): {flags}")
    for key in ACC:
        frac = stats["reset_" + key] / max(stats["points"], 1)
        if frac > 0.01:
            log(tag, f"WARNING: {ACC[key][0]} decreased at {frac * 100:.1f}% of points/intervals; "
                "it may not be a running total - check before using this month")
    return tag, True


# =============================================================================
# Main
# =============================================================================
def month_range(start, end):
    y, m = start
    while (y, m) <= end:
        yield y, m
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)


def parse_ym(s):
    y, m = s.split("-")
    return int(y), int(m)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exp", choices=sorted(EXPERIMENTS), default="historical")
    p.add_argument("--wrf-root", help="override the experiment's WRF history root")
    p.add_argument("--out-dir", help="override the experiment's output directory")
    p.add_argument("--start", required=True, type=parse_ym, help="YYYY-MM")
    p.add_argument("--end", required=True, type=parse_ym, help="YYYY-MM (inclusive)")
    p.add_argument("--geo", help="geo_em.d01.nc or wrfinput_d01 (optional, for later)")
    p.add_argument("--box", nargs=4, type=float, default=BOX,
                   metavar=("LONMIN", "LONMAX", "LATMIN", "LATMAX"))
    p.add_argument("--res", type=float, default=RES, help="target grid step [deg]")
    p.add_argument("--yorig", type=int, default=YORIG)
    p.add_argument("--month-fmt", default=MONTH_FMT,
                   help="'d' -> M1 (make_ERA5.m default), '02d' -> M01")
    p.add_argument("--vars", nargs="+", default=VARIABLES, choices=VARIABLES)
    p.add_argument("--nproc", type=int, default=1, help="months processed in parallel")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--no-compress", action="store_true")
    p.add_argument("--dry-run", action="store_true",
                   help="set up grid and land mask, report, write nothing")
    a = p.parse_args()

    exp = EXPERIMENTS[a.exp]
    wrf_root = a.wrf_root or exp["wrf_root"]
    out_dir = a.out_dir or exp["out_dir"]
    if not wrf_root:
        sys.exit(f"no WRF root for '{a.exp}': set it in EXPERIMENTS or pass --wrf-root")

    months = list(month_range(a.start, a.end))
    first = wrf_path(wrf_root, dt.datetime(*a.start, 1))
    if not os.path.exists(first):
        sys.exit(f"first WRF file not found: {first}")
    with Dataset(first) as nc:
        attrs = {k: nc.getncattr(k) for k in nc.ncattrs()}
        nx, ny = len(nc.dimensions["west_east"]), len(nc.dimensions["south_north"])
    grid = LambertGrid(attrs, nx, ny)
    tgt = build_target(grid, a.box, a.res)
    lsm_native = land_mask(first, grid, tgt["slab"], a.geo)
    lsm = tgt["interp"](lsm_native)
    lsm[lsm < 1e-6] = 0.0

    corners = grid.ij_to_ll(np.array([0, nx - 1, 0, nx - 1]), np.array([0, 0, ny - 1, ny - 1]))
    log("grid", f"WRF {nx}x{ny}, dx {grid.dx / 1e3:.1f} km, cone {grid.n:.5f}; corners "
        + ", ".join(f"({la:.2f}N,{lo:.2f}E)" for la, lo in zip(*corners)))
    sj, si = tgt["slab"]
    log("grid", f"target {len(tgt['lon'])}x{len(tgt['lat'])} pts at {a.res} deg, "
        f"lon {tgt['lon'][0]}-{tgt['lon'][-1]}, lat {tgt['lat'][0]}-{tgt['lat'][-1]}; "
        f"WRF slab j {sj.start}:{sj.stop}, i {si.start}:{si.stop}; "
        f"sea points {(lsm == 0).mean() * 100:.1f}%")
    log("run", f"{a.exp}: {wrf_root} -> {out_dir}; {len(months)} month(s); vars {a.vars}")
    if a.dry_run:
        return

    os.makedirs(out_dir, exist_ok=True)
    ctx = {"wrf_root": wrf_root, "out_dir": out_dir, "tgt": tgt, "lsm": lsm,
           "variables": a.vars, "yorig": a.yorig, "month_fmt": a.month_fmt,
           "compress": not a.no_compress, "overwrite": a.overwrite}
    work = functools.partial(process_month, ctx=ctx)
    if a.nproc > 1:
        with mp.get_context("fork").Pool(a.nproc) as pool:
            results = pool.map(work, months, chunksize=1)
    else:
        results = [work(m) for m in months]
    failed = [t for t, ok in results if not ok]
    log("run", f"done: {len(results) - len(failed)} ok, {len(failed)} failed {failed}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
