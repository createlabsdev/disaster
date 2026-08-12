#!/usr/bin/env python3
"""
Sentinel-2 NDVI change detection for 2024 Wayanad landslide scars at Meppadi.

Methodology:
  - Pre-event composite:  Sentinel-2 L2A, Jun 1-30, 2024 (before Jul 30 event)
  - Post-event composite: Aug 1-31, 2024 (after Wayanad landslide disaster)
  - Cloud masking using QA60 band (< 20% cloud cover filter)
  - Compute NDVI = (B8 - B4) / (B8 + B4) for both composites
  - Delta NDVI = pre_ndvi - post_ndvi
  - Threshold: dNDVI > 0.3 => landslide scar (major vegetation loss)
  - Additional filter: slope > 10 degrees (exclude flat agricultural changes)
  - Mask out water bodies (WorldCover class 80) and urban (class 50)

Alternative data sources (if GEE unavailable):
  - NASA COOLR (Cooperative Open Online Landslide Repository):
    https://maps.nccs.nasa.gov/arcgis/apps/MapAndAppGallery/index.html
  - GSI Bhukosh (Geological Survey of India):
    https://bhukosh.gsi.gov.in/
  - Wayanad 2024 event is well-documented; landslide inventories may be
    available from KSDMA (Kerala State Disaster Management Authority) or
    published research papers.

Requirements:
  pip install earthengine-api rasterio numpy
  Optional: pip install geemap  (for local export without Google Drive)

Usage:
  python gee_landslide_scars.py                # export to Google Drive
  python gee_landslide_scars.py --local        # save locally (needs geemap)
  python gee_landslide_scars.py --threshold 0.25  # custom dNDVI threshold
"""

import argparse
import os
import sys
import numpy as np

# Meppadi study area
BBOX = [76.10, 11.50, 76.22, 11.60]  # [west, south, east, north]
SITE = 'meppadi'
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output', SITE)
OUTPUT_FILE = 'landslide_inventory_2024.tif'

# Temporal windows
PRE_EVENT_START  = '2024-01-01'
PRE_EVENT_END    = '2024-05-31'
POST_EVENT_START = '2024-10-01'
POST_EVENT_END   = '2024-12-31'


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
            print("\nFallback options for Meppadi landslide inventory:")
            print("  1. NASA COOLR: https://maps.nccs.nasa.gov/arcgis/apps/MapAndAppGallery/index.html")
            print("  2. GSI Bhukosh: https://bhukosh.gsi.gov.in/")
            print("  3. KSDMA reports for Wayanad 2024")
            print("  Rasterize polygons and pass via --label-source to run_ml_pipeline.py")
            sys.exit(1)
    return ee


def cloud_mask_s2(image):
    """Mask clouds in Sentinel-2 using QA60 band."""
    import ee as _ee
    qa = image.select('QA60')
    # Bits 10 and 11 are clouds and cirrus
    cloud_bit_mask = 1 << 10
    cirrus_bit_mask = 1 << 11
    mask = qa.bitwiseAnd(cloud_bit_mask).eq(0).And(
           qa.bitwiseAnd(cirrus_bit_mask).eq(0))
    return image.updateMask(mask)


def get_ndvi_composite(ee, aoi, start_date, end_date):
    """Get cloud-masked Sentinel-2 NDVI composite."""
    collection = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
        .map(cloud_mask_s2))

    count = collection.size().getInfo()
    print(f"  Found {count} Sentinel-2 scenes for {start_date} to {end_date}")

    if count == 0:
        # Relax cloud filter
        print("  WARNING: No scenes with <20% cloud. Trying <50%...")
        collection = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
            .filterBounds(aoi)
            .filterDate(start_date, end_date)
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 50))
            .map(cloud_mask_s2))
        count = collection.size().getInfo()
        print(f"  Found {count} scenes with relaxed cloud filter")

    if count == 0:
        raise RuntimeError(f"No Sentinel-2 scenes found for {start_date} to {end_date}")

    # Compute NDVI for each image, then take median composite
    def add_ndvi(image):
        ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
        return image.addBands(ndvi)

    return collection.map(add_ndvi).select('NDVI').median()


def compute_landslide_mask(ee, pre_ndvi, post_ndvi, aoi, threshold=0.3):
    """Compute binary landslide scar mask from NDVI change."""
    # Delta NDVI: positive means vegetation loss (landslide scar)
    dndvi = pre_ndvi.subtract(post_ndvi)

    print(f"  Using dNDVI threshold: {threshold}")

    # Threshold: dNDVI > threshold => landslide
    scar_raw = dndvi.gt(threshold)

    # Slope filter: only keep scars on slopes > 10 degrees
    # Use SRTM DEM available in GEE
    srtm = ee.Image('USGS/SRTMGL1_003')
    slope = ee.Terrain.slope(srtm)
    slope_mask = slope.gt(10)
    scar_sloped = scar_raw.And(slope_mask)

    # Mask out water bodies and urban areas using ESA WorldCover
    worldcover = ee.ImageCollection('ESA/WorldCover/v200').first()
    water = worldcover.eq(80)      # water bodies
    urban = worldcover.eq(50)      # built-up
    land_mask = water.Not().And(urban.Not())
    scar_masked = scar_sloped.And(land_mask)

    # Morphological cleaning: remove isolated pixels
    scar_clean = scar_masked.focal_min(1).focal_max(1)

    return scar_clean.rename('landslide').toByte()


def export_to_drive(ee, landslide_mask, aoi, scale=10):
    """Export landslide mask to Google Drive as GeoTIFF."""
    task = ee.batch.Export.image.toDrive(
        image=landslide_mask,
        description='meppadi_landslide_inventory_2024',
        folder='digital_twin',
        fileNamePrefix='landslide_inventory_2024',
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
    return task


def export_local(ee, landslide_mask, aoi, output_path, scale=10):
    """Export landslide mask directly to local GeoTIFF using geemap."""
    try:
        import geemap
    except ImportError:
        print("ERROR: geemap not installed for local export. Run: pip install geemap")
        print("  Falling back to Drive export...")
        return export_to_drive(ee, landslide_mask, aoi, scale)

    import rasterio
    from rasterio.warp import reproject, Resampling

    print("[INFO] Downloading landslide mask via geemap...")

    temp_path = output_path.replace('.tif', '_10m.tif')
    geemap.ee_export_image(
        landslide_mask,
        filename=temp_path,
        scale=scale,
        region=aoi,
        crs='EPSG:4326',
        file_per_band=False
    )
    print(f"  Saved 10m landslide mask: {temp_path}")

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
        print(f"  Saved 90m landslide mask: {output_path}")
        os.remove(temp_path)
    else:
        print(f"  WARNING: DEM reference not found at {dem_ref}, keeping 10m resolution")
        os.rename(temp_path, output_path)

    # Print summary stats
    with rasterio.open(output_path) as ds:
        data = ds.read(1)
        n_scar = np.sum(data == 1)
        n_total = np.sum(data != 255)
        pct = 100.0 * n_scar / max(n_total, 1)
        print(f"\n  Landslide inventory summary:")
        print(f"    Scar pixels:    {n_scar}")
        print(f"    Total pixels:   {n_total}")
        print(f"    Scar fraction:  {pct:.1f}%")

    return None


def main():
    parser = argparse.ArgumentParser(
        description='Extract 2024 Wayanad landslide scars for Meppadi using Sentinel-2 NDVI')
    parser.add_argument('--local', action='store_true',
                        help='Save output locally instead of Google Drive (requires geemap)')
    parser.add_argument('--threshold', type=float, default=0.3,
                        help='NDVI drop threshold (default: 0.3)')
    parser.add_argument('--scale', type=int, default=10,
                        help='Export resolution in meters (default: 10)')
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("============================================================")
    print("Meppadi Landslide Scars — Sentinel-2 NDVI Change Detection")
    print("============================================================")
    

    print(f"  Pre-event:  {PRE_EVENT_START} to {PRE_EVENT_END}")
    print(f"  Post-event: {POST_EVENT_START} to {POST_EVENT_END}")

    # Authenticate
    ee = authenticate_ee()

    # Define AOI
    aoi = ee.Geometry.Rectangle(BBOX)

    # Get NDVI composites
    print("\n[1/3] Building pre-event NDVI composite...")
    pre_ndvi = get_ndvi_composite(ee, aoi, PRE_EVENT_START, PRE_EVENT_END)

    print("\n[2/3] Building post-event NDVI composite...")
    post_ndvi = get_ndvi_composite(ee, aoi, POST_EVENT_START, POST_EVENT_END)

    # Compute landslide mask
    print("\n[3/3] Computing landslide scar mask...")
    landslide_mask = compute_landslide_mask(
        ee, pre_ndvi, post_ndvi, aoi,
        threshold=args.threshold
    )

    # Export
    output_path = os.path.join(OUTPUT_DIR, OUTPUT_FILE)
    if args.local:
        export_local(ee, landslide_mask, aoi, output_path, scale=args.scale)
    else:
        export_to_drive(ee, landslide_mask, aoi, scale=args.scale)

    print("\n" + "=" * 60)
    print("Done.")
    print("=" * 60)


if __name__ == '__main__':
    main()
