#!/usr/bin/env python3
"""
Download SILO daily GeoTIFFs from s3://silo-open-data (public, unsigned) for
the Adelaide River catchment and extract area-weighted subcatchment means.

    python scripts/download_silo.py --start 1900 --end 2024
    python scripts/download_silo.py --start 1900 --end 2024 --skip-download
"""

import argparse
import datetime
import os
import sys
from pathlib import Path

import boto3
import geopandas as gpd
import numpy as np
from botocore import UNSIGNED
from botocore.config import Config


BUCKET = "silo-open-data"
S3_PREFIX = "Official/daily"
VARIABLES = {
    "daily_rain": {"label": "Rainfall", "unit": "millimeter", "header": "Rainfall[millimeter]:Step_Accumulated"},
    "evap_morton_lake": {"label": "Evap", "unit": "millimeter", "header": "Evaporation[millimeter]:Step_Accumulated"},
}

# Adelaide River catchment bounding box (WGS84), with 0.1° buffer
BBOX = {
    "lon_min": 130.78,
    "lon_max": 131.77,
    "lat_min": -13.74,
    "lat_max": -12.07,
}

# SILO grid resolution (degrees)
SILO_RES = 0.05


def parse_args():
    p = argparse.ArgumentParser(description="Download SILO data from AWS for Adelaide River catchment")
    p.add_argument("--start", type=int, default=1900, help="Start year (inclusive)")
    p.add_argument("--end", type=int, default=2024, help="End year (inclusive)")
    p.add_argument("--variables", nargs="+", default=list(VARIABLES.keys()),
                   help="SILO variables to download")
    p.add_argument("--shapefile", type=str, default="data/spatial/NAMCatchments.shp",
                   help="Subcatchment polygon shapefile for area-weighting")
    p.add_argument("--id-col", type=str, default="MUID",
                   help="Shapefile attribute with subcatchment IDs")
    p.add_argument("--outdir", type=str, default="data/forcing",
                   help="Output root directory")
    p.add_argument("--cache-dir", type=str, default="data/cache/silo_tif",
                   help="Local cache for downloaded GeoTIFFs")
    p.add_argument("--skip-download", action="store_true",
                   help="Skip S3 download, only run extraction from cached TIFs")
    p.add_argument("--skip-extract", action="store_true",
                   help="Only download TIFs, skip subcatchment extraction")
    p.add_argument("--workers", type=int, default=8,
                   help="Parallel download threads")
    return p.parse_args()



def download_tifs(variable, years, cache_dir, workers=8):
    """Download daily GeoTIFFs for one variable from S3 (unsigned)."""
    import concurrent.futures

    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED),
                      region_name="ap-southeast-2")

    var_cache = Path(cache_dir) / variable
    var_cache.mkdir(parents=True, exist_ok=True)

    tasks = []
    for year in years:
        d = datetime.date(year, 1, 1)
        end = datetime.date(year, 12, 31)
        while d <= end:
            fname = f"{d.strftime('%Y%m%d')}.{variable}.tif"
            s3_key = f"{S3_PREFIX}/{variable}/{year}/{fname}"
            local = var_cache / str(year) / fname
            if not local.exists():
                tasks.append((s3_key, local))
            d += datetime.timedelta(days=1)

    if not tasks:
        print(f"  [{variable}] All {sum(366 if y % 4 == 0 else 365 for y in years)} files already cached.")
        return

    print(f"  [{variable}] Downloading {len(tasks)} GeoTIFFs ...")

    def _download(item):
        key, path = item
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            s3.download_file(BUCKET, key, str(path))
            return True
        except Exception as e:
            print(f"    WARN: {key} -> {e}", file=sys.stderr)
            return False

    ok = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(_download, tasks):
            if result:
                ok += 1
                if ok % 1000 == 0:
                    print(f"    ... {ok}/{len(tasks)}")
    print(f"  [{variable}] Downloaded {ok}/{len(tasks)} files.")



def build_subcatchment_weights(shapefile, id_col, variable, cache_dir):
    """Fractional-coverage pixel weights per subcatchment via exactextract.

    Returns {subcatchment_id: [(row, col, weight), ...]} with weights
    normalised to sum to 1 within each polygon.
    """
    import rasterio
    from exactextract import exact_extract

    gdf = gpd.read_file(shapefile).to_crs(epsg=4326)

    # Reference TIF for grid definition
    var_cache = Path(cache_dir) / variable
    ref_tif = None
    for root, dirs, files in os.walk(var_cache):
        for f in sorted(files):
            if f.endswith(".tif"):
                ref_tif = os.path.join(root, f)
                break
        if ref_tif:
            break

    if ref_tif is None:
        raise FileNotFoundError(f"No cached TIF found for {variable} in {var_cache}")

    with rasterio.open(ref_tif) as src:
        transform = src.transform
        height, width = src.height, src.width

    coverage_results = exact_extract(
        ref_tif, gdf, ["cell_id", "coverage", "values"],
        include_cols=[id_col],
        output="pandas",
    )

    weights = {}
    for _, row in gdf.iterrows():
        sid = str(row[id_col])

        sub = coverage_results[coverage_results[id_col] == row[id_col]]
        if sub.empty:
            print(f"    WARN: subcatchment {sid} — no pixel coverage from exactextract", file=sys.stderr)
            # Centroid fallback
            cx, cy = row.geometry.centroid.x, row.geometry.centroid.y
            c = int((cx - transform.c) / transform.a)
            r = int((cy - transform.f) / transform.e)
            weights[sid] = [(r, c, 1.0)]
            continue

        cell_ids = sub["cell_id"].values
        coverages = sub["coverage"].values

        # Convert flat cell_id to (row, col) in the raster grid
        rows = cell_ids // width
        cols = cell_ids % width

        # Normalise so weights sum to 1 per subcatchment
        total_cov = coverages.sum()
        if total_cov > 0:
            normed = coverages / total_cov
        else:
            normed = np.ones_like(coverages) / len(coverages)

        weights[sid] = [(int(r), int(c), float(w))
                        for r, c, w in zip(rows, cols, normed)]

    total_pairs = sum(len(v) for v in weights.values())
    print(f"    Built fractional-coverage weights for {len(weights)} subcatchments "
          f"({total_pairs} pixel-subcatchment pairs)")
    return weights, transform, (height, width)


def extract_subcatchment_means(variable, years, cache_dir, weights, outdir, var_meta):
    """Read cached TIFs, compute area-weighted means, write per-subcatchment text files."""
    import rasterio

    var_cache = Path(cache_dir) / variable
    out_var = Path(outdir) / var_meta["label"].lower()
    out_var.mkdir(parents=True, exist_ok=True)

    # Collect all dates
    dates = []
    for year in years:
        d = datetime.date(year, 1, 1)
        end = datetime.date(year, 12, 31)
        while d <= end:
            dates.append(d)
            d += datetime.timedelta(days=1)

    # Pre-allocate: {subcatchment_id: [values]}
    sids = sorted(weights.keys())
    data = {sid: [] for sid in sids}

    print(f"  [{variable}] Extracting {len(dates)} days for {len(sids)} subcatchments ...")

    for i, d in enumerate(dates):
        fname = f"{d.strftime('%Y%m%d')}.{variable}.tif"
        tif_path = var_cache / str(d.year) / fname

        if not tif_path.exists():
            # Missing day — fill with NA
            for sid in sids:
                data[sid].append(float("nan"))
            continue

        with rasterio.open(str(tif_path)) as src:
            band = src.read(1)  # 2D array
            nodata = src.nodata

        for sid in sids:
            vals = []
            ws = []
            for r, c, w in weights[sid]:
                if 0 <= r < band.shape[0] and 0 <= c < band.shape[1]:
                    v = band[r, c]
                    if nodata is not None and v == nodata:
                        continue
                    vals.append(float(v))
                    ws.append(w)
            if vals:
                total_w = sum(ws)
                data[sid].append(sum(v * w for v, w in zip(vals, ws)) / total_w)
            else:
                data[sid].append(float("nan"))

        if (i + 1) % 5000 == 0:
            print(f"    ... {i + 1}/{len(dates)} days")

    # Write output files
    start_year = years[0]
    end_year = years[-1]
    label = var_meta["label"]
    header = var_meta["header"]
    time_suffix = "00:00:00" if variable == "daily_rain" else "09:00:00"

    for sid in sids:
        out_file = out_var / f"{label}_{start_year}-{end_year}_{sid}.txt"
        with open(out_file, "w") as f:
            f.write(f"{header}\n")
            f.write("Time\tCurrent\n")
            for j, d in enumerate(dates):
                val = data[sid][j]
                ts = f"{d.strftime('%Y-%m-%d')} {time_suffix}"
                f.write(f"{ts}\t{val:.6g}\n")

    print(f"  [{variable}] Wrote {len(sids)} files to {out_var}/")



def main():
    args = parse_args()
    years = list(range(args.start, args.end + 1))
    print(f"SILO downloader: {args.start}–{args.end}, variables: {args.variables}")
    print(f"Bounding box: lon [{BBOX['lon_min']}, {BBOX['lon_max']}], lat [{BBOX['lat_min']}, {BBOX['lat_max']}]")

    if not args.skip_download:
        print("\nDownloading GeoTIFFs from S3 ...")
        for var in args.variables:
            if var not in VARIABLES:
                print(f"  SKIP unknown variable: {var}", file=sys.stderr)
                continue
            download_tifs(var, years, args.cache_dir, workers=args.workers)
    else:
        print("  (skipped, --skip-download)")

    if not args.skip_extract:
        print("\nExtracting subcatchment area-weighted means ...")
        if not Path(args.shapefile).exists():
            print(f"  ERROR: Shapefile not found: {args.shapefile}", file=sys.stderr)
            print("  Place your subcatchment polygons shapefile at the path above,")
            print("  or use --shapefile to point to a different location.")
            sys.exit(1)

        for var in args.variables:
            if var not in VARIABLES:
                continue
            print(f"\n  Building pixel weights for {var} ...")
            weights, transform, shape = build_subcatchment_weights(
                args.shapefile, args.id_col, var, args.cache_dir
            )
            extract_subcatchment_means(
                var, years, args.cache_dir, weights, args.outdir, VARIABLES[var]
            )
    else:
        print("  (skipped, --skip-extract)")

    print("\nDone.")


if __name__ == "__main__":
    main()
