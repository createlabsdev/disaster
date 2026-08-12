#!/usr/bin/env python3
"""
Sentinel-1 SAR change detection for 2018 Kerala flood extent at Chellanam.

Methodology:
  - Pre-flood composite:    Sentinel-1 GRD, VH polarization, DESCENDING orbit,
                            Jul 1-25, 2018 (before Kerala flood peak)
  - During-flood composite: Aug 15-25, 2018 (peak of 2018 Kerala floods)
  - Compute backscatter difference: pre_vh_dB - during_vh_dB
  - Threshold: difference > 3 dB => flooded (SAR backscatter drops over water)
  - Also implements Otsu automatic thresholding as alternative
  - Mask out permanent water bodies using JRC Global Surface Water

Fallback (if GEE unavailable):
  Download Copernicus Emergency Management Service activation EMSR294:
  https://emergency.copernicus.eu/mapping/list-of-components/EMSR294
  The delineation products provide flood extent polygons for Kerala 2018.
  Rasterize the shapefile to a binary GeoTIFF matching the DEM grid.

Requirements:
  pip install earthengine-api rasterio numpy
  Optional: pip install geemap  (for local export without Google Drive)

Usage:
  python gee_flood_extent.py                # export to Google Drive
  python gee_flood_extent.py --local        # save locally (needs geemap)
  python gee_flood_extent.py --threshold 4  # custom dB threshold
"""

import argparse
import os
import sys
import numpy as np

# Chellanam study area
BBOX = [76.28, 9.78, 76.35, 9.83]  # [west, south, east, north]
SITE = 'chellanam'
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output', SITE)
OUTPUT_FILE = 'flood_extent_2018.tif'

# Temporal windows
PRE_FLOOD_START  = '2018-07-01'
PRE_FLOOD_END    = '2018-07-25'
FLOOD_START      = '2018-08-15'
FLOOD_END        = '2018-08-25'


def authenticate_ee():
    """Authenticate and initialize Google Earth Engine."""
    try:
        import ee
    except ImportError:
        print("ERROR: earthengine-api not installed. Run: pip install earthengine-api")
        sys.exit(1)

    try:
        ee.Initialize(project='flood-502410')
        print("[OK] Earth Engine already authenticated and initialized.")
    except Exception:
        print("[INFO] Authenticating Earth Engine (browser will open)...")
        try:
            ee.Authenticate()
            ee.Initialize(project='flood-502410')
            print("[OK] Earth Engine authenticated and initialized.")
        except Exception as e:
            print(f"ERROR: Could not authenticate Earth Engine: {e}")
            print("\nFallback: Download Copernicus EMS EMSR294 flood polygons from:")
            print("  https://emergency.copernicus.eu/mapping/list-of-components/EMSR294")
            print("Rasterize them to a binary GeoTIFF and pass via --label-source to run_ml_pipeline.py")
            sys.exit(1)
    return ee


def get_s1_composite(ee, aoi, start_date, end_date):
    """Get median Sentinel-1 GRD VH backscatter composite."""
    collection = (ee.ImageCollection('COPERNICUS/S1_GRD')
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .select('VH'))

    count = collection.size().getInfo()
    print(f"  Found {count} Sentinel-1 scenes for {start_date} to {end_date}")

    if count == 0:
        print("  WARNING: No scenes found! Trying ASCENDING orbit...")
        collection = (ee.ImageCollection('COPERNICUS/S1_GRD')
            .filterBounds(aoi)
            .filterDate(start_date, end_date)
            .filter(ee.Filter.eq('instrumentMode', 'IW'))
            .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
            .select('VH'))
        count = collection.size().getInfo()
        print(f"  Found {count} ASCENDING scenes")

    if count == 0:
        raise RuntimeError(f"No Sentinel-1 scenes found for {start_date} to {end_date}")

    # Apply speckle filter (focal median, 50m radius)
    def speckle_filter(image):
        return image.focal_median(50, 'circle', 'meters').copyProperties(image, image.propertyNames())

    return collection.map(speckle_filter).median()


def compute_flood_mask(ee, pre_vh, during_vh, aoi, threshold_db=3.0, use_otsu=False):
    """Compute binary flood mask from SAR backscatter change."""
    # Difference: pre - during (positive = backscatter dropped = flooded)
    diff = pre_vh.subtract(during_vh)

    if use_otsu:
        # Otsu thresholding on the difference image
        print("  Computing Otsu threshold...")
        histogram = diff.reduceRegion(
            reducer=ee.Reducer.histogram(255, 0.1),
            geometry=aoi,
            scale=10,
            maxPixels=1e9
        )
        hist_info = histogram.getInfo()
        vh_hist = hist_info.get('VH', {})
        counts = vh_hist.get('histogram', [])
        bucket_means = vh_hist.get('bucketMeans', [])

        if counts and bucket_means:
            counts = np.array(counts)
            bucket_means = np.array(bucket_means)
            # Otsu's method: find threshold that maximizes between-class variance
            total = counts.sum()
            best_thresh = bucket_means[0]
            best_var = 0
            sum_total = (counts * bucket_means).sum()
            sum_bg = 0
            weight_bg = 0
            for i in range(len(counts)):
                weight_bg += counts[i]
                if weight_bg == 0:
                    continue
                weight_fg = total - weight_bg
                if weight_fg == 0:
                    break
                sum_bg += counts[i] * bucket_means[i]
                mean_bg = sum_bg / weight_bg
                mean_fg = (sum_total - sum_bg) / weight_fg
                var_between = weight_bg * weight_fg * (mean_bg - mean_fg) ** 2
                if var_between > best_var:
                    best_var = var_between
                    best_thresh = bucket_means[i]
            threshold_db = float(best_thresh)
            print(f"  Otsu threshold: {threshold_db:.2f} dB")
        else:
            print(f"  WARNING: Otsu histogram empty, falling back to {threshold_db} dB")
    else:
        print(f"  Using fixed threshold: {threshold_db} dB")

    # Apply threshold
    flood_raw = diff.gt(threshold_db)

    # Mask permanent water bodies using JRC Global Surface Water
    jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
    permanent_water = jrc.select('occurrence').gt(80)
    flood_mask = flood_raw.And(permanent_water.Not())

    # Clean up: remove isolated pixels (morphological opening)
    flood_mask = flood_mask.focal_min(1).focal_max(1)

    return flood_mask.rename('flood').toByte()


def export_to_drive(ee, flood_mask, aoi, scale=10):
    """Export flood mask to Google Drive as GeoTIFF."""
    task = ee.batch.Export.image.toDrive(
        image=flood_mask,
        description='chellanam_flood_extent_2018',
        folder='digital_twin',
        fileNamePrefix='flood_extent_2018',
        region=aoi,
        scale=scale,
        crs='EPSG:4326',
        maxPixels=1e9,
        fileFormat='GeoTIFF'
    )
    task.start()
    print(f"\n[EXPORT] Task started: {task.id}")
    print("  Check progress at: https://code.earthengine.google.com/tasks")
    print(f"  Output will appear in Google Drive folder 'digital_twin'")
    print(f"\n  After download, copy to: {os.path.join(OUTPUT_DIR, OUTPUT_FILE)}")
    print(f"  Then resample to 90m: gdal_translate -outsize <cols> <rows> -r average <src> <dst>")
    return task


def export_local(ee, flood_mask, aoi, output_path, scale=10):
    """Export flood mask directly to local GeoTIFF using geemap."""
    try:
        import geemap
    except ImportError:
        print("ERROR: geemap not installed for local export. Run: pip install geemap")
        print("  Falling back to Drive export...")
        return export_to_drive(ee, flood_mask, aoi, scale)

    import rasterio
    from rasterio.transform import from_bounds
    from rasterio.warp import reproject, Resampling

    print("[INFO] Downloading flood mask via geemap...")
    west, south, east, north = BBOX

    # Export at native 10m, then resample to 90m
    temp_path = output_path.replace('.tif', '_10m.tif')
    geemap.ee_export_image(
        flood_mask,
        filename=temp_path,
        scale=scale,
        region=aoi,
        crs='EPSG:4326',
        file_per_band=False
    )
    print(f"  Saved 10m flood mask: {temp_path}")

    # Resample to match DEM (90m)
    dem_ref = os.path.join(OUTPUT_DIR, 'dem_ref.tif')
    if os.path.exists(dem_ref):
        print("  Resampling to match DEM grid (90m)...")
        with rasterio.open(dem_ref) as ref:
            ref_transform = ref.transform
            ref_crs = ref.crs
            ref_width = ref.width
            ref_height = ref.height
            ref_profile = ref.profile.copy()

        with rasterio.open(temp_path) as src:
            src_data = src.read(1)
            src_transform = src.transform
            src_crs = src.crs

        dst_data_float = np.empty((ref_height, ref_width), dtype=np.float32)
        reproject(
            source=src_data.astype(np.float32),
            destination=dst_data_float,
            src_transform=src_transform,
            src_crs=src_crs,
            dst_transform=ref_transform,
            dst_crs=ref_crs,
            resampling=Resampling.average
        )
        dst_data = np.where(dst_data_float > 0.05, 1, 0).astype(np.uint8)

        ref_profile.update(dtype='uint8', nodata=255)
        with rasterio.open(output_path, 'w', **ref_profile) as dst:
            dst.write(dst_data, 1)
        print(f"  Saved 90m flood mask: {output_path}")

        # Clean up temp
        os.remove(temp_path)
    else:
        print(f"  WARNING: DEM reference not found at {dem_ref}, keeping 10m resolution")
        os.rename(temp_path, output_path)

    # Print summary stats
    with rasterio.open(output_path) as ds:
        data = ds.read(1)
        n_flood = np.sum(data == 1)
        n_total = np.sum(data != 255)
        pct = 100.0 * n_flood / max(n_total, 1)
        print(f"\n  Flood extent summary:")
        print(f"    Flooded pixels:  {n_flood}")
        print(f"    Total pixels:    {n_total}")
        print(f"    Flood fraction:  {pct:.1f}%")

    return None


def main():
    parser = argparse.ArgumentParser(
        description='Extract 2018 Kerala flood extent for Chellanam using Sentinel-1 SAR')
    parser.add_argument('--local', action='store_true',
                        help='Save output locally instead of Google Drive (requires geemap)')
    parser.add_argument('--threshold', type=float, default=3.0,
                        help='SAR backscatter drop threshold in dB (default: 3.0)')
    parser.add_argument('--otsu', action='store_true',
                        help='Use Otsu automatic thresholding instead of fixed threshold')
    parser.add_argument('--scale', type=int, default=10,
                        help='Export resolution in meters (default: 10)')
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print("Chellanam Flood Extent — Sentinel-1 SAR Change Detection")
    print("=" * 60)
    print(f"  Study area: {BBOX}")
    print(f"  Pre-flood:  {PRE_FLOOD_START} to {PRE_FLOOD_END}")
    print(f"  Flood:      {FLOOD_START} to {FLOOD_END}")

    # Authenticate
    ee = authenticate_ee()

    # Define AOI
    aoi = ee.Geometry.Rectangle(BBOX)

    # Get composites
    print("\n[1/3] Building pre-flood SAR composite...")
    pre_vh = get_s1_composite(ee, aoi, PRE_FLOOD_START, PRE_FLOOD_END)

    print("\n[2/3] Building during-flood SAR composite...")
    during_vh = get_s1_composite(ee, aoi, FLOOD_START, FLOOD_END)

    # Compute flood mask
    print("\n[3/3] Computing flood mask...")
    flood_mask = compute_flood_mask(
        ee, pre_vh, during_vh, aoi,
        threshold_db=args.threshold,
        use_otsu=args.otsu
    )

    # Export
    output_path = os.path.join(OUTPUT_DIR, OUTPUT_FILE)
    if args.local:
        export_local(ee, flood_mask, aoi, output_path, scale=args.scale)
    else:
        export_to_drive(ee, flood_mask, aoi, scale=args.scale)

    print("\n" + "=" * 60)
    print("Done.")
    print("=" * 60)


if __name__ == '__main__':
    main()
