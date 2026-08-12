#!/usr/bin/env python3
"""
Build labeled training data for ML susceptibility modeling.

Extracts features from existing rasters, computes terrain derivatives,
and samples labeled points (event=1, non-event=0) with slope-stratified
negative sampling to avoid terrain distribution bias.

Outputs:
  output/{site}/training_data.csv     - labeled samples with all features
  output/{site}/feature_stack.tif     - multi-band GeoTIFF for prediction

Usage:
  python build_training_data.py chellanam
  python build_training_data.py meppadi --label-source path/to/mask.tif
  python build_training_data.py chellanam --n-samples 10000 --seed 123
"""

import argparse
import os
import sys
import warnings
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling
from scipy.ndimage import distance_transform_edt, binary_dilation, generate_binary_structure

warnings.filterwarnings('ignore', category=rasterio.errors.NotGeoreferencedWarning)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

SITES = {
    'chellanam': {
        'bbox': (76.28, 9.78, 76.35, 9.83),
        'default_label': 'flood_extent_2018.tif',
    },
    'meppadi': {
        'bbox': (76.05, 11.55, 76.15, 11.65),
        'default_label': 'landslide_inventory_2024.tif',
    },
}

# Features to extract from existing rasters (relative to output/{site}/)
RASTER_FEATURES = {
    'elevation':        'dem_ref.tif',
    'soil_depth':       'soil_depth.tif',
    'cohesion':         'cohesion.tif',
    'phi':              'phi.tif',
    'k_factor':         'k_factor.tif',
    'c_factor':         'c_factor.tif',
    'r_factor':         'r_factor.tif',
    'manning_n':        'manning_n.tif',
    'rusle_soil_loss':  'rusle_soil_loss.tif',
    'factor_of_safety': 'factor_of_safety.tif',
}

# Features that may be in UTM (need reprojection)
UTM_RASTER_FEATURES = {
    'avaflow_depth':    'anuga_depth_max.tif',
    'avaflow_velocity': 'anuga_velocity_max.tif',
}

# Feature band names (order matters for feature_stack.tif)
TERRAIN_FEATURES = ['elevation', 'slope', 'aspect', 'plan_curvature',
                    'profile_curvature', 'twi', 'dist_to_stream']
ALL_FEATURE_NAMES = (TERRAIN_FEATURES +
                     list(RASTER_FEATURES.keys())[1:] +  # skip elevation (already in terrain)
                     list(UTM_RASTER_FEATURES.keys()) +
                     ['landcover'])


# ────────────────────────── Terrain Derivatives ──────────────────────────

def compute_slope_aspect(dem, x_res, y_res):
    """Compute slope (degrees) and aspect (degrees 0-360) from DEM."""
    dzdx = np.gradient(dem, x_res, axis=1)
    dzdy = np.gradient(dem, y_res, axis=0)
    slope_rad = np.arctan(np.sqrt(dzdx**2 + dzdy**2))
    slope_deg = np.degrees(slope_rad)
    aspect = np.degrees(np.arctan2(-dzdy, dzdx))
    aspect = np.mod(aspect, 360.0)
    return slope_deg, aspect


def compute_curvatures(dem, x_res, y_res):
    """Compute plan and profile curvature from DEM."""
    dzdx = np.gradient(dem, x_res, axis=1)
    dzdy = np.gradient(dem, y_res, axis=0)
    d2zdx2 = np.gradient(dzdx, x_res, axis=1)
    d2zdy2 = np.gradient(dzdy, y_res, axis=0)
    d2zdxy = np.gradient(dzdx, y_res, axis=0)

    p = dzdx**2 + dzdy**2
    q = p + 1.0

    # Profile curvature (in direction of steepest slope)
    profile = np.where(p > 1e-10,
        -(dzdx**2 * d2zdx2 + 2*dzdx*dzdy*d2zdxy + dzdy**2 * d2zdy2) / (p * np.sqrt(q**3)),
        0.0)

    # Plan curvature (perpendicular to steepest slope)
    plan = np.where(p > 1e-10,
        -(dzdy**2 * d2zdx2 - 2*dzdx*dzdy*d2zdxy + dzdx**2 * d2zdy2) / (p**1.5),
        0.0)

    # Clip extreme values from numerical artifacts
    profile = np.clip(profile, -10, 10)
    plan = np.clip(plan, -10, 10)

    return plan.astype(np.float32), profile.astype(np.float32)


def d8_flow_accumulation(dem):
    """Simple D8 flow accumulation."""
    nrows, ncols = dem.shape
    acc = np.ones((nrows, ncols), dtype=np.float64)

    # Sort cells by elevation (highest first)
    flat_indices = np.argsort(-dem.ravel())

    # D8 neighbor offsets (row, col)
    neighbors = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]

    for idx in flat_indices:
        r, c = divmod(idx, ncols)
        elev = dem[r, c]
        if np.isnan(elev):
            continue

        # Find steepest downslope neighbor
        min_elev = elev
        min_r, min_c = -1, -1
        for dr, dc in neighbors:
            nr, nc_n = r + dr, c + dc
            if 0 <= nr < nrows and 0 <= nc_n < ncols:
                ne = dem[nr, nc_n]
                if not np.isnan(ne) and ne < min_elev:
                    min_elev = ne
                    min_r, min_c = nr, nc_n

        # Route flow to steepest neighbor
        if min_r >= 0:
            acc[min_r, min_c] += acc[r, c]

    return acc


def compute_twi(dem, slope_deg, x_res, y_res):
    """Compute Topographic Wetness Index = ln(a / tan(slope))."""
    print("    Computing D8 flow accumulation (may take a moment)...")
    acc = d8_flow_accumulation(dem)

    # Contributing area: acc * cell_area / cell_width
    cell_area = x_res * y_res
    cell_width = (x_res + y_res) / 2.0
    contributing_area = acc * cell_area / cell_width

    # TWI = ln(a / tan(slope))
    slope_rad = np.radians(slope_deg)
    tan_slope = np.maximum(np.tan(slope_rad), 0.001)  # avoid log(inf)
    twi = np.log(contributing_area / tan_slope)

    return np.clip(twi, -5, 30).astype(np.float32)


def compute_dist_to_stream(acc, x_res, y_res, threshold_pct=5.0):
    """Compute Euclidean distance to nearest stream channel.

    Streams defined as top threshold_pct of flow accumulation.
    """
    threshold = np.percentile(acc[acc > 0], 100 - threshold_pct)
    streams = acc >= threshold
    # Distance transform: distance from each non-stream cell to nearest stream
    # Invert: distance_transform works on 0-valued cells
    dist = distance_transform_edt(~streams, sampling=[y_res, x_res])
    return dist.astype(np.float32)


# ────────────────────────── Raster I/O ──────────────────────────

def load_reference(site):
    """Load DEM reference grid (defines output shape/CRS/transform)."""
    dem_path = os.path.join(BASE_DIR, 'output', site, 'dem_ref.tif')
    with rasterio.open(dem_path) as ds:
        dem = ds.read(1).astype(np.float32)
        profile = ds.profile.copy()
        transform = ds.transform
        crs = ds.crs
    return dem, profile, transform, crs


def load_and_align(raster_path, ref_profile, ref_transform, ref_crs):
    """Load a raster and reproject/resample to match reference grid."""
    ref_height = ref_profile['height']
    ref_width = ref_profile['width']

    if not os.path.exists(raster_path):
        print(f"    WARNING: {raster_path} not found, filling with 0")
        return np.zeros((ref_height, ref_width), dtype=np.float32)

    with rasterio.open(raster_path) as src:
        src_data = src.read(1).astype(np.float32)
        src_transform = src.transform
        src_crs = src.crs

        # Check if reprojection is needed
        if src_crs == ref_crs and src_data.shape == (ref_height, ref_width):
            return src_data

    dst_data = np.zeros((ref_height, ref_width), dtype=np.float32)
    reproject(
        source=src_data,
        destination=dst_data,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=ref_transform,
        dst_crs=ref_crs,
        resampling=Resampling.bilinear
    )
    return dst_data


def load_landcover(site, ref_profile, ref_transform, ref_crs):
    """Load and align ESA WorldCover landcover."""
    lc_path = os.path.join(BASE_DIR, f'{site}_landcover.tif')
    if not os.path.exists(lc_path):
        print(f"    WARNING: {lc_path} not found")
        return np.zeros((ref_profile['height'], ref_profile['width']), dtype=np.float32)

    ref_height = ref_profile['height']
    ref_width = ref_profile['width']

    with rasterio.open(lc_path) as src:
        src_data = src.read(1).astype(np.float32)
        src_transform = src.transform
        src_crs = src.crs

    dst_data = np.zeros((ref_height, ref_width), dtype=np.float32)
    reproject(
        source=src_data,
        destination=dst_data,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=ref_transform,
        dst_crs=ref_crs,
        resampling=Resampling.nearest  # nearest for categorical
    )
    return dst_data


# ────────────────────────── Sampling ──────────────────────────

def slope_stratified_sample(pos_rows, pos_cols, neg_candidates_mask,
                            slope, n_samples, rng):
    """Sample negatives matching the slope distribution of positives."""
    # Get slope values for positives
    pos_slopes = slope[pos_rows, pos_cols]

    # Create slope bins (10 quantile bins)
    bin_edges = np.percentile(pos_slopes[~np.isnan(pos_slopes)],
                              np.linspace(0, 100, 11))
    bin_edges[0] = -999
    bin_edges[-1] = 999

    # Count positives per bin → target counts for negatives
    pos_bins = np.digitize(pos_slopes, bin_edges) - 1
    bin_counts = np.bincount(pos_bins, minlength=10)
    target_per_bin = (bin_counts / max(bin_counts.sum(), 1) * n_samples).astype(int)

    # Get all negative candidate pixels
    neg_r, neg_c = np.where(neg_candidates_mask)
    neg_slopes = slope[neg_r, neg_c]
    neg_bins = np.digitize(neg_slopes, bin_edges) - 1

    sampled_r, sampled_c = [], []
    for b in range(10):
        mask_b = neg_bins == b
        indices_b = np.where(mask_b)[0]
        n_target = min(target_per_bin[b], len(indices_b))
        if n_target > 0:
            chosen = rng.choice(indices_b, size=n_target, replace=False)
            sampled_r.extend(neg_r[chosen])
            sampled_c.extend(neg_c[chosen])

    # If we didn't get enough, fill from remaining
    n_got = len(sampled_r)
    if n_got < n_samples:
        remaining_mask = neg_candidates_mask.copy()
        for r, c in zip(sampled_r, sampled_c):
            remaining_mask[r, c] = False
        rem_r, rem_c = np.where(remaining_mask)
        n_fill = min(n_samples - n_got, len(rem_r))
        if n_fill > 0:
            chosen = rng.choice(len(rem_r), size=n_fill, replace=False)
            sampled_r.extend(rem_r[chosen])
            sampled_c.extend(rem_c[chosen])

    return np.array(sampled_r[:n_samples]), np.array(sampled_c[:n_samples])


# ────────────────────────── Main ──────────────────────────

def build_features(site, out_dir):
    """Build complete feature stack from all rasters + terrain derivatives."""
    print(f"\n[1/4] Loading reference DEM...")
    dem, profile, transform, crs = load_reference(site)
    nrows, ncols = dem.shape
    print(f"  Grid: {nrows} x {ncols}, CRS: {crs}")

    # Replace nodata with NaN
    dem = np.where(dem == 0, np.nan, dem)

    # Cell size (approximate, in meters for geographic CRS)
    x_res_deg = abs(transform.a)
    y_res_deg = abs(transform.e)
    # Convert degrees to approximate meters at this latitude
    lat_center = transform.f + transform.e * nrows / 2
    x_res_m = x_res_deg * 111320 * np.cos(np.radians(lat_center))
    y_res_m = y_res_deg * 110540
    print(f"  Approximate cell size: {x_res_m:.0f} x {y_res_m:.0f} m")

    # ── Terrain derivatives ──
    print("\n[2/4] Computing terrain derivatives...")
    dem_filled = np.nan_to_num(dem, nan=np.nanmedian(dem))

    print("  Slope & aspect...")
    slope, aspect = compute_slope_aspect(dem_filled, x_res_m, y_res_m)

    print("  Curvatures...")
    plan_curv, prof_curv = compute_curvatures(dem_filled, x_res_m, y_res_m)

    print("  TWI (topographic wetness index)...")
    twi = compute_twi(dem_filled, slope, x_res_m, y_res_m)

    print("  Flow accumulation & distance to streams...")
    acc = d8_flow_accumulation(dem_filled)
    dist_stream = compute_dist_to_stream(acc, x_res_m, y_res_m)

    features = {
        'elevation': dem_filled.astype(np.float32),
        'slope': slope.astype(np.float32),
        'aspect': aspect.astype(np.float32),
        'plan_curvature': plan_curv,
        'profile_curvature': prof_curv,
        'twi': twi,
        'dist_to_stream': dist_stream,
    }

    # ── Existing rasters ──
    print("\n[3/4] Loading existing raster features...")
    for name, filename in RASTER_FEATURES.items():
        if name == 'elevation':
            continue  # already loaded
        raster_path = os.path.join(out_dir, filename)
        print(f"  {name}: {filename}")
        features[name] = load_and_align(raster_path, profile, transform, crs)

    # UTM rasters (need reprojection)
    for name, filename in UTM_RASTER_FEATURES.items():
        raster_path = os.path.join(out_dir, filename)
        print(f"  {name}: {filename} (UTM → reproject)")
        features[name] = load_and_align(raster_path, profile, transform, crs)

    # Landcover
    print(f"  landcover: {site}_landcover.tif")
    features['landcover'] = load_landcover(site, profile, transform, crs)

    # ── Save feature stack ──
    print("\n[4/4] Saving feature stack...")
    n_bands = len(ALL_FEATURE_NAMES)
    stack_profile = profile.copy()
    stack_profile.update(dtype='float32', count=n_bands, nodata=np.nan)

    stack_path = os.path.join(out_dir, 'feature_stack.tif')
    with rasterio.open(stack_path, 'w', **stack_profile) as dst:
        for i, name in enumerate(ALL_FEATURE_NAMES):
            band_data = features.get(name, np.zeros((nrows, ncols), dtype=np.float32))
            band_data = np.nan_to_num(band_data, nan=0.0)
            dst.write(band_data.astype(np.float32), i + 1)
            dst.set_band_description(i + 1, name)

    print(f"  Saved: {stack_path} ({n_bands} bands)")

    return features, slope, profile, transform, crs


def extract_training_data(site, features, slope, label_path, n_samples, seed, out_dir):
    """Extract labeled training samples with stratified negative sampling."""
    rng = np.random.default_rng(seed)
    _, profile, transform, crs = load_reference(site)
    nrows, ncols = profile['height'], profile['width']

    # Load label raster
    print(f"\n  Loading label raster: {label_path}")
    label = load_and_align(label_path, profile, transform, crs)
    label = np.where(label > 0.5, 1, 0).astype(np.uint8)

    n_event = np.sum(label == 1)
    n_total = label.size
    print(f"  Event pixels: {n_event} / {n_total} ({100*n_event/n_total:.1f}%)")

    # Create valid mask (exclude nodata, water bodies)
    landcover = features.get('landcover', np.zeros((nrows, ncols)))
    valid_mask = (landcover != 80)  # exclude water
    elev = features.get('elevation', np.zeros((nrows, ncols)))
    valid_mask &= ~np.isnan(elev) & (elev != 0)

    if n_event == 0:
        print("  INFO: No event pixels found — this is a SAFE site.")
        print("  Sampling only negative (non-event) pixels...")
        # For safe sites, sample random negatives across the entire valid area
        neg_candidate_mask = (label == 0) & valid_mask
        neg_r, neg_c = np.where(neg_candidate_mask)
        n_neg = min(n_samples, len(neg_r))
        if n_neg > 0:
            chosen = rng.choice(len(neg_r), size=n_neg, replace=False)
            neg_r, neg_c = neg_r[chosen], neg_c[chosen]
        print(f"  Sampled {len(neg_r)} negative (safe) pixels")
        
        all_r = neg_r
        all_c = neg_c
        labels = np.zeros(len(neg_r))
        
        data = {'row': all_r, 'col': all_c, 'label': labels.astype(int)}
        for name in ALL_FEATURE_NAMES:
            feat = features.get(name, np.zeros((nrows, ncols), dtype=np.float32))
            values = feat[all_r, all_c]
            data[name] = values.astype(np.float32)
        
        df = pd.DataFrame(data)
        for col in ALL_FEATURE_NAMES:
            mask = ~np.isfinite(df[col])
            if mask.any():
                median_val = df.loc[~mask, col].median() if (~mask).any() else 0.0
                df.loc[mask, col] = median_val
        
        csv_path = os.path.join(out_dir, 'training_data.csv')
        df.to_csv(csv_path, index=False)
        print(f"\n  Saved: {csv_path}")
        print(f"  Shape: {df.shape}")
        return df



    # Buffer around event boundaries (3 pixels) to avoid mixed pixels
    struct = generate_binary_structure(2, 2)  # 8-connectivity
    event_dilated = binary_dilation(label == 1, structure=struct, iterations=3)
    boundary_buffer = event_dilated & (label == 0)

    # ── Sample positives ──
    pos_mask = (label == 1) & valid_mask
    pos_r, pos_c = np.where(pos_mask)
    n_pos = min(n_samples, len(pos_r))
    if n_pos < len(pos_r):
        chosen = rng.choice(len(pos_r), size=n_pos, replace=False)
        pos_r, pos_c = pos_r[chosen], pos_c[chosen]
    print(f"  Sampled {len(pos_r)} positive (event) pixels")

    # ── Sample negatives (slope-stratified) ──
    neg_candidate_mask = (label == 0) & valid_mask & ~boundary_buffer
    n_neg = min(n_samples, np.sum(neg_candidate_mask))
    neg_r, neg_c = slope_stratified_sample(
        pos_r, pos_c, neg_candidate_mask, slope, n_neg, rng)
    print(f"  Sampled {len(neg_r)} negative (non-event) pixels (slope-stratified)")

    # ── Build DataFrame ──
    all_r = np.concatenate([pos_r, neg_r])
    all_c = np.concatenate([pos_c, neg_c])
    labels = np.concatenate([np.ones(len(pos_r)), np.zeros(len(neg_r))])

    data = {'row': all_r, 'col': all_c, 'label': labels.astype(int)}
    for name in ALL_FEATURE_NAMES:
        feat = features.get(name, np.zeros((nrows, ncols), dtype=np.float32))
        values = feat[all_r, all_c]
        data[name] = values.astype(np.float32)

    df = pd.DataFrame(data)

    # Replace NaN/inf with column median
    for col in ALL_FEATURE_NAMES:
        mask = ~np.isfinite(df[col])
        if mask.any():
            median_val = df.loc[np.isfinite(df[col]), col].median()
            df.loc[mask, col] = median_val if np.isfinite(median_val) else 0.0

    # Save
    csv_path = os.path.join(out_dir, 'training_data.csv')
    df.to_csv(csv_path, index=False)
    print(f"\n  Saved: {csv_path}")
    print(f"  Shape: {df.shape}")

    # Print feature statistics
    print("\n  Feature statistics:")
    print(f"  {'Feature':<22s} {'Min':>10s} {'Max':>10s} {'Mean':>10s} {'Std':>10s}")
    print("  " + "-" * 62)
    for name in ALL_FEATURE_NAMES:
        col = df[name]
        print(f"  {name:<22s} {col.min():>10.3f} {col.max():>10.3f} "
              f"{col.mean():>10.3f} {col.std():>10.3f}")

    return df


def main():
    parser = argparse.ArgumentParser(
        description='Build labeled training data for ML susceptibility model')
    parser.add_argument('site', help='Study site name')
    parser.add_argument('--label-source', type=str, default=None,
                        help='Path to custom label raster (binary GeoTIFF, 1=event)')
    parser.add_argument('--n-samples', type=int, default=5000,
                        help='Number of samples per class (default: 5000)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Override output directory')
    args = parser.parse_args()

    site = args.site
    out_dir = args.output_dir or os.path.join(BASE_DIR, 'output', site)
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 60)
    print(f"Building Training Data — {site}")
    print("=" * 60)

    # Determine label source
    if args.label_source:
        label_path = args.label_source
    else:
        label_path = SITES.get(site, {}).get('default_label')
        if label_path:
            label_path = os.path.join(out_dir, label_path)

    skip_labels = False
    if label_path == 'dummy' or not label_path or not os.path.exists(label_path):
        print(f"\nWarning: Label raster not found or dummy provided.")
        print(f"Skipping training data CSV generation. Will only build feature stack.")
        skip_labels = True

    # Build feature stack
    features, slope, profile, transform, crs = build_features(site, out_dir)

    if skip_labels:
        print("\n" + "=" * 60)
        print(f"Done. Feature stack ready at output/{site}/feature_stack.tif")
        print("Skipped training CSV generation (no labels).")
        print("=" * 60)
        return

    # Extract training data
    df = extract_training_data(
        site, features, slope, label_path,
        n_samples=args.n_samples, seed=args.seed, out_dir=out_dir)

    print("\n" + "=" * 60)
    print(f"Done. Training data ready for train_susceptibility.py")
    print("=" * 60)


if __name__ == '__main__':
    main()
