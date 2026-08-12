#!/usr/bin/env python3
"""
run_batch_training.py — Mass-Scale Multi-Site Training Pipeline

Processes 40 sites (10 landslides + 10 floods + 20 safe zones) to build
a globally-trained AI model for disaster risk prediction.

Usage:
    python run_batch_training.py                       # Process all 40 sites
    python run_batch_training.py --category landslide   # Only landslide sites
    python run_batch_training.py --skip-physics         # Skip physics (use existing)
    python run_batch_training.py --skip-download        # Skip GEE downloads
    python run_batch_training.py --train-only           # Only merge CSVs & train
"""

import argparse
import os
import sys
import subprocess
import time
import numpy as np
import pandas as pd
import rasterio

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)

from sites_registry import TRAINING_SITES, get_sites_by_category


def create_blank_label(site_name, out_dir):
    """Create an all-zero label raster matching the DEM grid for safe sites."""
    dem_path = os.path.join(out_dir, 'dem_ref.tif')
    label_path = os.path.join(out_dir, 'safe_label.tif')
    with rasterio.open(dem_path) as src:
        profile = src.profile.copy()
        profile.update(dtype='uint8', count=1, nodata=255)
        data = np.zeros((src.height, src.width), dtype=np.uint8)
    with rasterio.open(label_path, 'w', **profile) as dst:
        dst.write(data, 1)
    print(f"  Created blank label: {label_path}")
    return label_path


def fetch_disaster_labels(site_name, site_config, out_dir):
    """Download real disaster labels from GEE for disaster sites."""
    label_type = site_config['label_type']
    bbox = site_config['bbox']

    if label_type == 'sentinel1_sar':
        label_path = os.path.join(out_dir, 'flood_extent.tif')
        print(f"  Fetching Sentinel-1 SAR flood labels for {site_name}...")
        try:
            import ee
            import geemap
            try:
                ee.Initialize(project='flood-502410')
            except Exception:
                ee.Authenticate()
                ee.Initialize(project='flood-502410')

            aoi = ee.Geometry.Rectangle(list(bbox))
            pre_start, pre_end = site_config['pre_event']
            post_start, post_end = site_config['post_event']

            # Pre-flood Sentinel-1
            pre = (ee.ImageCollection('COPERNICUS/S1_GRD')
                   .filterBounds(aoi)
                   .filterDate(pre_start, pre_end)
                   .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
                   .filter(ee.Filter.eq('instrumentMode', 'IW'))
                   .select('VH')
                   .median())

            # During-flood Sentinel-1
            post = (ee.ImageCollection('COPERNICUS/S1_GRD')
                    .filterBounds(aoi)
                    .filterDate(post_start, post_end)
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
                    .filter(ee.Filter.eq('instrumentMode', 'IW'))
                    .select('VH')
                    .median())

            # Difference and threshold
            diff = pre.subtract(post)
            flooded = diff.gt(3.0).selfMask()

            # Mask permanent water
            jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence')
            perm_water = jrc.gt(80)
            flooded = flooded.updateMask(perm_water.Not())

            # Export at 90m to match DEM
            flooded_binary = flooded.unmask(0).clip(aoi)

            geemap.ee_export_image(
                flooded_binary,
                filename=label_path,
                scale=90,
                region=aoi,
                crs='EPSG:4326',
                file_per_band=False
            )
            print(f"  Saved flood label: {label_path}")
            return label_path
        except Exception as e:
            print(f"  WARNING: Failed to fetch flood labels: {e}")
            print(f"  Creating blank label as fallback.")
            return create_blank_label(site_name, out_dir)

    elif label_type == 'sentinel2_ndvi':
        label_path = os.path.join(out_dir, 'landslide_scars.tif')
        print(f"  Fetching Sentinel-2 NDVI landslide labels for {site_name}...")
        try:
            import ee
            import geemap
            try:
                ee.Initialize(project='flood-502410')
            except Exception:
                ee.Authenticate()
                ee.Initialize(project='flood-502410')

            aoi = ee.Geometry.Rectangle(list(bbox))
            pre_start, pre_end = site_config['pre_event']
            post_start, post_end = site_config['post_event']

            def mask_clouds(image):
                qa = image.select('QA60')
                cloud_mask = qa.bitwiseAnd(1 << 10).eq(0).And(qa.bitwiseAnd(1 << 11).eq(0))
                return image.updateMask(cloud_mask)

            # Pre-event NDVI
            pre = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                   .filterBounds(aoi)
                   .filterDate(pre_start, pre_end)
                   .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
                   .map(mask_clouds)
                   .median())
            pre_ndvi = pre.normalizedDifference(['B8', 'B4']).rename('ndvi')

            # Post-event NDVI
            post = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                    .filterBounds(aoi)
                    .filterDate(post_start, post_end)
                    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
                    .map(mask_clouds)
                    .median())
            post_ndvi = post.normalizedDifference(['B8', 'B4']).rename('ndvi')

            # dNDVI and threshold
            dndvi = pre_ndvi.subtract(post_ndvi)
            scars = dndvi.gt(0.15)  # 0.15 threshold (more sensitive)

            # Slope filter (>10 degrees)
            dem = ee.Image('CGIAR/SRTM90_V4')
            slope = ee.Terrain.slope(dem)
            scars = scars.And(slope.gt(10))

            # Mask water and urban
            worldcover = ee.ImageCollection('ESA/WorldCover/v200').first()
            water = worldcover.eq(80)
            urban = worldcover.eq(50)
            scars = scars.updateMask(water.Not()).updateMask(urban.Not())

            scars_binary = scars.unmask(0).clip(aoi)

            geemap.ee_export_image(
                scars_binary,
                filename=label_path,
                scale=90,
                region=aoi,
                crs='EPSG:4326',
                file_per_band=False
            )
            print(f"  Saved landslide label: {label_path}")
            return label_path
        except Exception as e:
            print(f"  WARNING: Failed to fetch landslide labels: {e}")
            print(f"  Creating blank label as fallback.")
            return create_blank_label(site_name, out_dir)

    else:
        return create_blank_label(site_name, out_dir)


def process_site(site_name, site_config, skip_download=False, skip_physics=False):
    """Process a single site: download data, run physics, build features."""
    bbox = site_config['bbox']
    category = site_config['category']
    out_dir = os.path.join(BASE_DIR, 'output', site_name)
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  Processing: {site_name} ({category})")
    print(f"  {site_config.get('description', '')}")
    print(f"  BBox: {bbox}")
    print(f"{'='*60}")

    start = time.time()

    # Step 1: Download DEM + Landcover
    if not skip_download:
        print(f"\n  [1/5] Downloading terrain data...")
        cmd = [sys.executable, os.path.join(BASE_DIR, 'init_site.py'),
               site_name, str(bbox[0]), str(bbox[1]), str(bbox[2]), str(bbox[3])]
        try:
            subprocess.run(cmd, check=True, cwd=BASE_DIR)
        except subprocess.CalledProcessError as e:
            print(f"  ERROR downloading data for {site_name}: {e}")
            return False
    else:
        print(f"\n  [1/5] Skipping download (--skip-download)")

    # Step 2: Compute soil params
    if not skip_physics:
        print(f"\n  [2/5] Computing soil parameters...")
        cmd = [sys.executable, os.path.join(BASE_DIR, 'compute_soil_params.py'), site_name]
        try:
            subprocess.run(cmd, check=True, cwd=BASE_DIR)
        except subprocess.CalledProcessError as e:
            print(f"  ERROR computing soil params for {site_name}: {e}")
            return False
    else:
        print(f"\n  [2/5] Skipping soil params (--skip-physics)")

    # Step 3: Run physics engines
    if not skip_physics:
        print(f"\n  [3/5] Running physics engines (RUSLE + FS + ANUGA)...")
        cmd = [sys.executable, os.path.join(BASE_DIR, 'run_pipeline_advanced.py'),
               '--site', site_name]
        try:
            subprocess.run(cmd, check=True, cwd=BASE_DIR)
        except subprocess.CalledProcessError as e:
            print(f"  ERROR running physics for {site_name}: {e}")
            return False
    else:
        print(f"\n  [3/5] Skipping physics (--skip-physics)")

    # Step 4: Get labels
    print(f"\n  [4/5] Getting labels...")
    if category == 'safe':
        label_path = create_blank_label(site_name, out_dir)
    else:
        if not skip_download:
            label_path = fetch_disaster_labels(site_name, site_config, out_dir)
        else:
            # Check if labels already exist
            label_path = None
            for candidate in ['flood_extent.tif', 'landslide_scars.tif',
                              'flood_extent_2018.tif', 'landslide_inventory_2024.tif',
                              'safe_label.tif']:
                candidate_path = os.path.join(out_dir, candidate)
                if os.path.exists(candidate_path):
                    label_path = candidate_path
                    break
            if label_path is None:
                print(f"  No existing labels found, creating blank.")
                label_path = create_blank_label(site_name, out_dir)

    # Step 5: Build training data
    print(f"\n  [5/5] Building feature stack and training data...")
    cmd = [sys.executable, os.path.join(BASE_DIR, 'build_training_data.py'),
           site_name, '--label-source', label_path, '--n-samples', '2000']
    try:
        subprocess.run(cmd, check=True, cwd=BASE_DIR)
    except subprocess.CalledProcessError as e:
        print(f"  ERROR building training data for {site_name}: {e}")
        return False

    elapsed = time.time() - start
    print(f"\n  ✓ {site_name} completed in {elapsed:.1f}s")
    return True


def merge_training_data(sites_to_process):
    """Merge all individual training_data.csv files into one global dataset."""
    print(f"\n{'='*60}")
    print(f"  MERGING TRAINING DATA FROM {len(sites_to_process)} SITES")
    print(f"{'='*60}")

    all_dfs = []
    for site_name in sites_to_process:
        csv_path = os.path.join(BASE_DIR, 'output', site_name, 'training_data.csv')
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            df['site'] = site_name
            df['site_id'] = list(sites_to_process).index(site_name)
            all_dfs.append(df)
            print(f"  {site_name}: {len(df)} samples "
                  f"(pos={int(df['label'].sum())}, neg={int((df['label']==0).sum())})")
        else:
            print(f"  WARNING: {csv_path} not found, skipping.")

    if not all_dfs:
        print("ERROR: No training data found!")
        sys.exit(1)

    merged = pd.concat(all_dfs, ignore_index=True)

    # Save merged dataset
    global_dir = os.path.join(BASE_DIR, 'output', 'global_model')
    os.makedirs(global_dir, exist_ok=True)
    merged_path = os.path.join(global_dir, 'training_data.csv')
    merged.to_csv(merged_path, index=False)

    n_pos = int(merged['label'].sum())
    n_neg = int((merged['label'] == 0).sum())
    n_sites = merged['site'].nunique()

    print(f"\n  Merged dataset: {len(merged)} total samples")
    print(f"  Positive (disaster): {n_pos}")
    print(f"  Negative (safe):     {n_neg}")
    print(f"  Sites represented:   {n_sites}")
    print(f"  Saved to: {merged_path}")

    return merged_path


def train_global_model():
    """Train the global model on merged multi-site data."""
    print(f"\n{'='*60}")
    print(f"  TRAINING GLOBAL AI MODEL")
    print(f"{'='*60}")

    global_dir = os.path.join(BASE_DIR, 'output', 'global_model')
    cmd = [sys.executable, os.path.join(BASE_DIR, 'train_susceptibility.py'),
           'global_model', '--model', 'both', '--no-map',
           '--output-dir', global_dir]
    subprocess.run(cmd, check=True, cwd=BASE_DIR)


def main():
    parser = argparse.ArgumentParser(
        description='Mass-scale multi-site training pipeline')
    parser.add_argument('--category', choices=['landslide', 'flood', 'safe', 'all'],
                        default='all', help='Which category of sites to process')
    parser.add_argument('--skip-download', action='store_true',
                        help='Skip GEE downloads (use existing data)')
    parser.add_argument('--skip-physics', action='store_true',
                        help='Skip physics simulations (use existing)')
    parser.add_argument('--train-only', action='store_true',
                        help='Only merge CSVs and train (skip all processing)')
    parser.add_argument('--sites', nargs='+', default=None,
                        help='Specific site names to process')
    args = parser.parse_args()

    total_start = time.time()

    # Determine which sites to process
    if args.sites:
        sites_to_process = {k: TRAINING_SITES[k] for k in args.sites if k in TRAINING_SITES}
    elif args.category == 'all':
        sites_to_process = TRAINING_SITES
    else:
        sites_to_process = get_sites_by_category(args.category)

    print(f"\n{'#'*60}")
    print(f"  MASS-SCALE MULTI-SITE TRAINING PIPELINE")
    print(f"  Sites to process: {len(sites_to_process)}")
    print(f"{'#'*60}")

    # Process each site
    if not args.train_only:
        successful = []
        failed = []
        for i, (site_name, site_config) in enumerate(sites_to_process.items()):
            print(f"\n\n>>> SITE {i+1}/{len(sites_to_process)}")
            ok = process_site(site_name, site_config,
                              skip_download=args.skip_download,
                              skip_physics=args.skip_physics)
            if ok:
                successful.append(site_name)
            else:
                failed.append(site_name)

        print(f"\n\n{'='*60}")
        print(f"  PROCESSING SUMMARY")
        print(f"{'='*60}")
        print(f"  Successful: {len(successful)}/{len(sites_to_process)}")
        if failed:
            print(f"  Failed: {', '.join(failed)}")

    # Merge and train
    merge_training_data(sites_to_process)
    train_global_model()

    total_elapsed = time.time() - total_start
    print(f"\n{'#'*60}")
    print(f"  GLOBAL TRAINING COMPLETE!")
    print(f"  Total time: {total_elapsed/60:.1f} minutes")
    print(f"  Global model saved to: output/global_model/")
    print(f"  ")
    print(f"  To predict risk for any new location, run:")
    print(f"  python generate_digital_twin.py <site> <W> <S> <E> <N> --model-site global_model")
    print(f"{'#'*60}")


if __name__ == '__main__':
    main()
