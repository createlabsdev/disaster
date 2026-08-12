#!/usr/bin/env python3
"""
Generate proxy hazard labels from existing physics outputs.

Use this when Google Earth Engine is unavailable to create real
inventory-based labels. These proxies let you test the full ML
pipeline end-to-end; swap in real labels later for production.

Chellanam (flood):
  - Labels low-elevation (<3m), low-slope (<5°), non-water pixels as flood-prone
  - Rationale: 2018 Kerala floods inundated low-lying coastal areas around Chellanam

Meppadi (landslide):
  - Labels pixels with Factor of Safety < 1.0 as landslide-prone
  - Rationale: FS<1 indicates slope failure under saturated conditions

Usage:
  python generate_proxy_labels.py chellanam
  python generate_proxy_labels.py meppadi
  python generate_proxy_labels.py both
"""

import os
import sys
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def generate_chellanam_flood_label():
    """Generate proxy flood label for Chellanam from elevation + slope."""
    site = 'chellanam'
    out_dir = os.path.join(BASE_DIR, 'output', site)
    dem_path = os.path.join(out_dir, 'dem_ref.tif')
    lc_path = os.path.join(BASE_DIR, f'{site}_landcover.tif')

    print(f"[{site}] Generating proxy flood label...")

    with rasterio.open(dem_path) as ds:
        dem = ds.read(1).astype(np.float32)
        profile = ds.profile.copy()
        transform = ds.transform

    # Compute slope
    x_res = abs(transform.a) * 111320 * np.cos(np.radians(9.8))  # approx m
    y_res = abs(transform.e) * 110540
    dzdx = np.gradient(dem, x_res, axis=1)
    dzdy = np.gradient(dem, y_res, axis=0)
    slope_deg = np.degrees(np.arctan(np.sqrt(dzdx**2 + dzdy**2)))

    # Load and align landcover
    if os.path.exists(lc_path):
        with rasterio.open(lc_path) as lc_ds:
            lc_data = lc_ds.read(1).astype(np.float32)
            lc_aligned = np.zeros_like(dem)
            reproject(
                source=lc_data,
                destination=lc_aligned,
                src_transform=lc_ds.transform,
                src_crs=lc_ds.crs,
                dst_transform=transform,
                dst_crs=profile['crs'],
                resampling=Resampling.nearest
            )
    else:
        lc_aligned = np.full_like(dem, 40)  # assume cropland

    # Proxy flood label:
    #   elevation < 3m AND slope < 5° AND not water (class 80)
    flood = ((dem < 3.0) & (slope_deg < 5.0) & (lc_aligned != 80)).astype(np.uint8)

    # Also include areas with elevation < 1m regardless of slope
    flood |= ((dem < 1.0) & (lc_aligned != 80)).astype(np.uint8)

    n_flood = np.sum(flood == 1)
    n_total = flood.size
    print(f"  Elevation range: {dem.min():.1f} to {dem.max():.1f} m")
    print(f"  Flood proxy pixels: {n_flood} / {n_total} ({100*n_flood/n_total:.1f}%)")

    if n_flood == 0:
        print("  WARNING: No flood pixels! Loosening threshold to elevation < 5m")
        flood = ((dem < 5.0) & (slope_deg < 10.0) & (lc_aligned != 80)).astype(np.uint8)
        n_flood = np.sum(flood == 1)
        print(f"  Flood proxy pixels (relaxed): {n_flood} / {n_total} ({100*n_flood/n_total:.1f}%)")

    out_path = os.path.join(out_dir, 'flood_extent_2018.tif')
    profile.update(dtype='uint8', nodata=255)
    with rasterio.open(out_path, 'w', **profile) as dst:
        dst.write(flood, 1)

    print(f"  Saved: {out_path}")
    return out_path


def generate_meppadi_landslide_label():
    """Generate proxy landslide label for Meppadi from Factor of Safety."""
    site = 'meppadi'
    out_dir = os.path.join(BASE_DIR, 'output', site)
    fs_path = os.path.join(out_dir, 'factor_of_safety.tif')
    lc_path = os.path.join(BASE_DIR, f'{site}_landcover.tif')

    print(f"\n[{site}] Generating proxy landslide label...")

    with rasterio.open(fs_path) as ds:
        fs = ds.read(1).astype(np.float32)
        profile = ds.profile.copy()
        transform = ds.transform

    # Load DEM for slope
    dem_path = os.path.join(out_dir, 'dem_ref.tif')
    with rasterio.open(dem_path) as ds:
        dem = ds.read(1).astype(np.float32)

    x_res = abs(transform.a) * 111320 * np.cos(np.radians(11.6))
    y_res = abs(transform.e) * 110540
    dzdx = np.gradient(dem, x_res, axis=1)
    dzdy = np.gradient(dem, y_res, axis=0)
    slope_deg = np.degrees(np.arctan(np.sqrt(dzdx**2 + dzdy**2)))

    # Load and align landcover
    if os.path.exists(lc_path):
        with rasterio.open(lc_path) as lc_ds:
            lc_data = lc_ds.read(1).astype(np.float32)
            lc_aligned = np.zeros_like(dem)
            reproject(
                source=lc_data,
                destination=lc_aligned,
                src_transform=lc_ds.transform,
                src_crs=lc_ds.crs,
                dst_transform=transform,
                dst_crs=profile['crs'],
                resampling=Resampling.nearest
            )
    else:
        lc_aligned = np.full_like(dem, 10)  # assume tree cover

    # Proxy landslide label:
    #   FS < 1.0 AND slope > 10° AND not water/urban
    landslide = (
        (fs < 1.0) &
        (fs > 0) &  # exclude nodata/zero
        (slope_deg > 10.0) &
        (lc_aligned != 80) &  # not water
        (lc_aligned != 50)    # not urban
    ).astype(np.uint8)

    n_ls = np.sum(landslide == 1)
    n_total = landslide.size
    print(f"  FS range: {fs[fs > 0].min():.2f} to {fs[fs < 999].max():.2f}")
    print(f"  Landslide proxy pixels: {n_ls} / {n_total} ({100*n_ls/n_total:.1f}%)")

    if n_ls == 0:
        print("  WARNING: No landslide pixels with FS<1! Using FS < 1.5 instead")
        landslide = (
            (fs < 1.5) & (fs > 0) & (slope_deg > 5.0) &
            (lc_aligned != 80) & (lc_aligned != 50)
        ).astype(np.uint8)
        n_ls = np.sum(landslide == 1)
        print(f"  Landslide proxy pixels (relaxed): {n_ls} / {n_total} ({100*n_ls/n_total:.1f}%)")

    out_path = os.path.join(out_dir, 'landslide_inventory_2024.tif')
    profile.update(dtype='uint8', nodata=255)
    with rasterio.open(out_path, 'w', **profile) as dst:
        dst.write(landslide, 1)

    print(f"  Saved: {out_path}")
    return out_path


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ('chellanam', 'meppadi', 'both'):
        print("Usage: python generate_proxy_labels.py <chellanam|meppadi|both>")
        sys.exit(1)

    site = sys.argv[1]

    print("=" * 60)
    print("Generating Proxy Hazard Labels")
    print("(Substitute for real GEE-derived labels)")
    print("=" * 60)

    if site in ('chellanam', 'both'):
        generate_chellanam_flood_label()
    if site in ('meppadi', 'both'):
        generate_meppadi_landslide_label()

    print("\n" + "=" * 60)
    print("Done. You can now run the ML pipeline:")
    if site == 'both':
        print("  python run_ml_pipeline.py chellanam --skip-gee")
        print("  python run_ml_pipeline.py meppadi --skip-gee")
    else:
        print(f"  python run_ml_pipeline.py {site} --skip-gee")
    print("\nNote: These are PROXY labels derived from physics outputs.")
    print("Replace with real inventory data (GEE, COOLR, EMS) for production.")
    print("=" * 60)


if __name__ == '__main__':
    main()
