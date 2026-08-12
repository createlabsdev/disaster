#!/usr/bin/env python3
"""
FINAL fully automated pipeline – no GRASS interaction:
  RUSLE → Infinite‑slope landslide initiation → ANUGA flood (transmissive boundary)
r.avaflow must be run manually later (see instructions below).

ANUGA boundary fix (2024-07-16): replaced Dirichlet_boundary([0,0,0]) with
Transmissive_stage_zero_momentum_boundary at the outlet. The old Dirichlet
boundary forced stage=0, which for Chellanam's sub-sea-level terrain (-14 to
+18m) drained rainfall instantly → all-zero depth/velocity output.
"""

import os, sys, subprocess, warnings
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling
from scipy.interpolate import griddata, RegularGridInterpolator

warnings.filterwarnings('ignore', category=DeprecationWarning)

import argparse

# ---------- Configuration ----------
# Default configurations for backward compatibility
DEFAULT_SITES = {
    'chellanam': {'bbox': (76.28, 9.78, 76.35, 9.83)},
    'meppadi':   {'bbox': (76.10, 11.50, 76.22, 11.60)},
}
DEFAULT_DESIGN_RAIN = {'duration_h': 3, 'intensity_mmh': 100}  # 100 mm/h for 3 h
MESH_RES_DEM = 90   # resolution in metres

# ---------- Infinite slope factor of safety ----------
def infinite_slope_fs(slope_deg, depth_m, cohesion_kpa, phi_deg, saturated=True):
    gamma = 18.0      # kN/m³
    gamma_w = 9.81    # kN/m³
    beta = np.radians(slope_deg)
    c = cohesion_kpa
    phi = np.radians(phi_deg)
    z = depth_m
    cos_beta = np.cos(beta)
    sin_beta = np.sin(beta)
    sigma_n = gamma * z * cos_beta**2
    u = gamma_w * z * cos_beta**2 if saturated else 0.0
    tau = gamma * z * sin_beta * cos_beta
    FS = np.where(tau > 1e-6, (c + (sigma_n - u) * np.tan(phi)) / tau, 999.0)
    return FS

# ---------- Utility: align a raster to the DEM grid ----------
def align_to_ref(src_tif, ref_tif, dst_tif):
    with rasterio.open(ref_tif) as ref:
        ref_transform = ref.transform
        ref_crs = ref.crs
        ref_width = ref.width
        ref_height = ref.height
    with rasterio.open(src_tif) as src:
        src_data = src.read(1).astype(np.float32)
        src_transform = src.transform
        src_crs = src.crs
    dst_data = np.empty((ref_height, ref_width), dtype=np.float32)
    reproject(
        source=src_data,
        destination=dst_data,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=ref_transform,
        dst_crs=ref_crs,
        resampling=Resampling.bilinear
    )
    with rasterio.open(ref_tif) as ref:
        profile = ref.profile.copy()
    profile.update(dtype='float32')
    with rasterio.open(dst_tif, 'w', **profile) as dst:
        dst.write(dst_data, 1)

# ---------- Auto-run compute_soil_params if needed ----------
def ensure_soil_params(site):
    """Run compute_soil_params.py if soil_depth.tif doesn't exist."""
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output', site)
    soil_depth = os.path.join(out_dir, 'soil_depth.tif')
    if not os.path.exists(soil_depth):
        print(f"[{site}] Running compute_soil_params.py...")
        subprocess.run([sys.executable, 'compute_soil_params.py', site], check=True)
    else:
        print(f"[{site}] Soil params exist, skipping compute_soil_params.py")

def main():
    parser = argparse.ArgumentParser(description="Run advanced physics pipeline")
    parser.add_argument("--site", type=str, help="Specific site to run (if omitted, runs default sites)")
    parser.add_argument("--rain-intensity", type=float, default=100.0, help="Rain intensity in mm/h")
    parser.add_argument("--rain-duration", type=float, default=3.0, help="Rain duration in hours")
    args = parser.parse_args()

    design_rain = {'duration_h': args.rain_duration, 'intensity_mmh': args.rain_intensity}

    if args.site:
        sites_to_run = [args.site]
    else:
        sites_to_run = list(DEFAULT_SITES.keys())

    for site in sites_to_run:
        print(f"\n===== Processing site: {site} =====")
        bbox = DEFAULT_SITES.get(site, {}).get('bbox', (0,0,0,0))
        out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output', site)
        os.makedirs(out_dir, exist_ok=True)

        # --- 0a. Ensure soil params are computed ---
        ensure_soil_params(site)

        # --- 0b. Prepare master DEM and align all inputs ---
        dem_raw = f'{site}_dem_90m.tif'
        dem_ref = os.path.join(out_dir, 'dem_ref.tif')
        if not os.path.exists(dem_ref) and os.path.exists(dem_raw):
            os.system(f'cp {dem_raw} {dem_ref}')



        # ====== 1. RUSLE ======
        print(f"[{site}] RUSLE...")
        with rasterio.open(dem_ref) as dem_src:
            dem = dem_src.read(1).astype(np.float32)
            transform = dem_src.transform
            x_res = transform[0]
            y_res = -transform[4]
            dzdx = np.gradient(dem, x_res, axis=1)
            dzdy = np.gradient(dem, y_res, axis=0)
            slope_rad = np.arctan(np.sqrt(dzdx**2 + dzdy**2))
            slope_deg = np.degrees(slope_rad)
            profile = dem_src.profile.copy()

        LS = np.clip((slope_deg / 10.0) ** 1.3, 0, 100)

        with rasterio.open(os.path.join(out_dir, 'r_factor.tif')) as r_src:
            R = r_src.read(1).astype(np.float32)
        # Scale R-factor by live rain intensity (baseline = 100 mm/h)
        rain_scale = design_rain['intensity_mmh'] / 100.0
        R = R * rain_scale
        print(f"  R-factor scaled by {rain_scale:.2f}x (live rain: {design_rain['intensity_mmh']} mm/h)")
        with rasterio.open(os.path.join(out_dir, 'k_factor.tif')) as k_src:
            K = k_src.read(1).astype(np.float32)
        with rasterio.open(os.path.join(out_dir, 'c_factor.tif')) as c_src:
            C = c_src.read(1).astype(np.float32)
        P = 1.0
        A = R * K * LS * C * P

        profile.update(dtype='float32')
        rusle_out = os.path.join(out_dir, 'rusle_soil_loss.tif')
        with rasterio.open(rusle_out, 'w', **profile) as dst:
            dst.write(A.astype(np.float32), 1)
        print(f"  {rusle_out}")

        # ====== 2. Physics‑based landslide release (infinite‑slope FS) ======
        print(f"[{site}] Computing Factor of Safety...")
        cohesion_tif = os.path.join(out_dir, 'cohesion.tif')
        phi_tif = os.path.join(out_dir, 'phi.tif')
        with rasterio.open(cohesion_tif) as c_src:
            cohesion = c_src.read(1).astype(np.float32)
        with rasterio.open(phi_tif) as phi_src:
            phi = phi_src.read(1).astype(np.float32)
        with rasterio.open(os.path.join(out_dir, 'soil_depth.tif')) as d_src:
            depth = d_src.read(1).astype(np.float32)
            depth_profile = d_src.profile.copy()

        FS = infinite_slope_fs(slope_deg, depth, cohesion, phi, saturated=True)
        fs_out = os.path.join(out_dir, 'factor_of_safety.tif')
        with rasterio.open(fs_out, 'w', **depth_profile) as dst:
            dst.write(FS.astype(np.float32), 1)

        release = np.where(FS < 1.0, depth, 0.0).astype(np.float32)
        release_tif = os.path.join(out_dir, 'landslide_release.tif')
        with rasterio.open(release_tif, 'w', **depth_profile) as dst:
            dst.write(release, 1)

        # ====== 3. ANUGA flood routing ======
        print(f"[{site}] ANUGA...")
        import anuga

        dem_utm = os.path.join(out_dir, 'dem_utm.tif')
        subprocess.run(['gdalwarp', '-t_srs', 'EPSG:32643', dem_ref, dem_utm], check=True)

        with rasterio.open(dem_utm) as dem_utm_src:
            bounds_utm = dem_utm_src.bounds
            elev_data = dem_utm_src.read(1).astype(np.float32)
            utm_transform = dem_utm_src.transform
            nodata_val = dem_utm_src.nodata
        west_utm, south_utm, east_utm, north_utm = bounds_utm

        if nodata_val is not None:
            elev_data[elev_data == nodata_val] = np.nan
        elev_data[elev_data < -500] = np.nan

        bounding_polygon = [
            (west_utm, south_utm),
            (east_utm, south_utm),
            (east_utm, north_utm),
            (west_utm, north_utm)
        ]

        south_mean = np.nanmean(elev_data[0, :])
        north_mean = np.nanmean(elev_data[-1, :])
        west_mean = np.nanmean(elev_data[:, 0])
        east_mean = np.nanmean(elev_data[:, -1])
        edges = {'bottom': south_mean, 'top': north_mean, 'left': west_mean, 'right': east_mean}
        outlet_edge = min(edges, key=edges.get)
        seg_map = {'bottom': 0, 'right': 1, 'top': 2, 'left': 3}
        outlet_seg = seg_map[outlet_edge]
        wall_segs = [i for i in range(4) if i != outlet_seg]

        domain = anuga.create_domain_from_regions(
            bounding_polygon=bounding_polygon,
            boundary_tags={'outlet': [outlet_seg], 'wall': wall_segs},
            maximum_triangle_area=MESH_RES_DEM**2
        )

        domain.set_zone(43)
        domain.geo_reference.south = False

        ny, nx = elev_data.shape
        x = np.linspace(west_utm + utm_transform.a/2, east_utm - utm_transform.a/2, nx)
        y = np.linspace(south_utm + utm_transform.a/2, north_utm - utm_transform.a/2, ny)
        elev_median = np.nanmedian(elev_data)
        elev_filled = np.where(np.isnan(elev_data), elev_median, elev_data)

        interpolator = RegularGridInterpolator(
            (np.sort(y), x), elev_filled[::-1, :],
            bounds_error=False, fill_value=elev_median
        )

        def elev_func(x_coord, y_coord):
            return interpolator(np.column_stack((y_coord, x_coord)))

        domain.set_quantity('elevation', elev_func, location='centroids')
        domain.quantities['stage'].set_values(
            domain.quantities['elevation'].centroid_values.copy(),
            location='centroids'
        )
        domain.set_quantity('friction', 0.03)

        domain.set_boundary({
            'outlet': anuga.Transmissive_stage_zero_momentum_boundary(domain),
            'wall': anuga.Reflective_boundary(domain)
        })

        rain_rate_m_s = (design_rain['intensity_mmh'] / 1000.0) / 3600.0
        print(f"  Rain rate: {design_rain['intensity_mmh']/3600.0:.4f} mm/s = {design_rain['intensity_mmh']} mm/h")
        rain = anuga.Rate_operator(domain, rate=rain_rate_m_s)
        
        centroid_coords = domain.get_centroid_coordinates()
        max_depth = np.zeros(len(centroid_coords), dtype=np.float32)
        max_velocity = np.zeros_like(max_depth)

        final_time = design_rain['duration_h'] * 3600
        for t in domain.evolve(yieldstep=300, finaltime=final_time):
            stage = domain.get_quantity('stage').centroid_values
            elev = domain.get_quantity('elevation').centroid_values
            depth_curr = np.maximum(0.0, stage - elev)
            xmom = domain.get_quantity('xmomentum').centroid_values
            ymom = domain.get_quantity('ymomentum').centroid_values
            depth_safe = np.maximum(depth_curr, 1e-6)
            vel_curr = np.sqrt((xmom / depth_safe)**2 + (ymom / depth_safe)**2)
            vel_curr = np.where(depth_curr > 0.001, vel_curr, 0.0)
            max_depth = np.maximum(max_depth, depth_curr)
            max_velocity = np.maximum(max_velocity, vel_curr)
            print(f"  Time: {domain.get_time():.0f}s  wet_cells={np.sum(depth_curr>0.001)}")

        grid_x, grid_y = np.meshgrid(x, y)
        depth_grid = griddata(centroid_coords, max_depth, (grid_x, grid_y), method='linear', fill_value=0.0)
        velocity_grid = griddata(centroid_coords, max_velocity, (grid_x, grid_y), method='linear', fill_value=0.0)

        depth_out = os.path.join(out_dir, 'anuga_depth_max.tif')
        vel_out = os.path.join(out_dir, 'anuga_velocity_max.tif')
        with rasterio.open(dem_utm) as src:
            out_profile = src.profile.copy()
        out_profile.update(dtype='float32')

        with rasterio.open(depth_out, 'w', **out_profile) as dst:
            dst.write(depth_grid.astype(np.float32), 1)
        with rasterio.open(vel_out, 'w', **out_profile) as dst:
            dst.write(velocity_grid.astype(np.float32), 1)

        print(f"\n===== {site} completed =====\n")

# ---------- Main ----------
if __name__ == '__main__':
    main()
