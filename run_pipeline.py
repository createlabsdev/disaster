#!/usr/bin/env python3
"""
Final working pipeline: RUSLE → Landslide Release Map → ANUGA Flood
Outputs saved in output/<site>/
"""

import os, sys, subprocess, warnings
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling
from scipy.interpolate import griddata

warnings.filterwarnings('ignore', category=DeprecationWarning)

# ---------- Configuration ----------
SITES = {
    'chellanam': {'bbox': (76.28, 9.78, 76.35, 9.83)},
    'meppadi':   {'bbox': (76.05, 11.55, 76.15, 11.65)},
}
DESIGN_RAIN = {'duration_h': 24, 'intensity_mmh': 50}   # 50 mm/h
SLOPE_THRESH = 30.0          # degrees
DEPTH_THRESH = 0.5           # metres

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

# ---------- Main pipeline per site ----------
def run_site(site, bbox):
    west, south, east, north = bbox
    base = os.getcwd()
    out_dir = os.path.join(base, 'output', site)
    os.makedirs(out_dir, exist_ok=True)

    # --- 0. Prepare master DEM and align all inputs ---
    dem_raw = f'{site}_dem_90m.tif'
    dem_ref = os.path.join(out_dir, 'dem_ref.tif')
    os.system(f'cp {dem_raw} {dem_ref}')

    layers_to_align = {
        f'{site}_c_factor.tif':       os.path.join(out_dir, 'c_factor.tif'),
        f'{site}_manning_n.tif':      os.path.join(out_dir, 'manning_n.tif'),
        f'soil_{site}/k_factor.tif':  os.path.join(out_dir, 'k_factor.tif'),
        f'soil_{site}/soil_depth.tif': os.path.join(out_dir, 'soil_depth.tif'),
        f'soil_{site}/cohesion.tif':  os.path.join(out_dir, 'cohesion.tif'),
        f'soil_{site}/phi.tif':       os.path.join(out_dir, 'phi.tif'),
        f'{site}_rfactor.tif':        os.path.join(out_dir, 'r_factor.tif'),
    }
    print(f"[{site}] Aligning inputs...")
    for src, dst in layers_to_align.items():
        if os.path.exists(src):
            align_to_ref(src, dem_ref, dst)
        else:
            print(f"  WARNING: missing {src}")

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

    # ====== 2. Landslide release map (slope + depth threshold) ======
    print(f"[{site}] Landslide release map...")
    with rasterio.open(os.path.join(out_dir, 'soil_depth.tif')) as depth_src:
        depth_arr = depth_src.read(1).astype(np.float32)
        depth_profile = depth_src.profile.copy()
    release = np.where((slope_deg > SLOPE_THRESH) & (depth_arr > DEPTH_THRESH), depth_arr, 0.0).astype(np.float32)
    release_tif = os.path.join(out_dir, 'landslide_release.tif')
    with rasterio.open(release_tif, 'w', **depth_profile) as dst:
        dst.write(release, 1)
    print(f"  {release_tif}")

    # ====== 3. ANUGA flood routing ======
    print(f"[{site}] ANUGA...")
    import anuga

    # Project DEM to UTM 43N
    dem_utm = os.path.join(out_dir, 'dem_utm.tif')
    subprocess.run(['gdalwarp', '-t_srs', 'EPSG:32643', dem_ref, dem_utm], check=True)

    # Read UTM DEM: bounds, array, transform
    with rasterio.open(dem_utm) as dem_utm_src:
        bounds_utm = dem_utm_src.bounds
        elev_data = dem_utm_src.read(1).astype(np.float32)
        utm_transform = dem_utm_src.transform
        dem_crs = dem_utm_src.crs
    west_utm, south_utm, east_utm, north_utm = bounds_utm

    # Bounding polygon (list of points)
    bounding_polygon = [
        (west_utm, south_utm),
        (east_utm, south_utm),
        (east_utm, north_utm),
        (west_utm, north_utm)
    ]

    # Create domain
    domain = anuga.create_domain_from_regions(
        bounding_polygon=bounding_polygon,
        boundary_tags={'exterior': [0, 1, 2, 3]},
        maximum_triangle_area=8100
    )

    # Set UTM zone and hemisphere
    domain.set_zone(43)
    domain.geo_reference.south = False

    # Elevation interpolation function
    ny, nx = elev_data.shape
    x = np.linspace(west_utm + utm_transform.a/2, east_utm - utm_transform.a/2, nx)
    y = np.linspace(south_utm + utm_transform.a/2, north_utm - utm_transform.a/2, ny)
    X, Y = np.meshgrid(x, y)
    points = np.column_stack((X.ravel(), Y.ravel()))
    values = elev_data.ravel()

    def elev_func(x_coord, y_coord):
        interp_vals = griddata(points, values, (x_coord, y_coord), method='linear', fill_value=0.0)
        return interp_vals

    domain.set_quantity('elevation', elev_func, location='centroids')
    domain.set_quantity('friction', 0.03)
    domain.set_boundary({'exterior': anuga.Reflective_boundary(domain)})

    # Rain-on-grid
    rain_rate = DESIGN_RAIN['intensity_mmh'] / 3600000.0   # m/s
    rain = anuga.Rate_operator(domain, rate=rain_rate)

    # --- Run simulation & track maximum depth/velocity manually ---
    # Get centroid coordinates for output mapping
    centroid_coords = domain.get_centroid_coordinates()
    num_elements = len(centroid_coords)
    max_depth = np.zeros(num_elements, dtype=np.float32)
    max_velocity = np.zeros(num_elements, dtype=np.float32)

    for t in domain.evolve(yieldstep=300, finaltime=3600):
        print(f"  Time: {domain.get_time():.0f} s")
        # Current depth = stage - elevation
        stage = domain.get_quantity('stage').centroid_values
        elev = domain.get_quantity('elevation').centroid_values
        depth = np.maximum(0.0, stage - elev)
        # Velocity components
        xvel = domain.get_quantity('xmomentum').centroid_values
        yvel = domain.get_quantity('ymomentum').centroid_values
        vel = np.sqrt(xvel**2 + yvel**2)
        # Update maxima
        max_depth = np.maximum(max_depth, depth)
        max_velocity = np.maximum(max_velocity, vel)

    # Interpolate from unstructured centroids to regular grid
    grid_x, grid_y = np.meshgrid(x, y)
    depth_grid = griddata(centroid_coords, max_depth, (grid_x, grid_y), method='linear', fill_value=0.0)
    velocity_grid = griddata(centroid_coords, max_velocity, (grid_x, grid_y), method='linear', fill_value=0.0)

    # Write GeoTIFFs
    depth_out = os.path.join(out_dir, 'anuga_depth_max.tif')
    vel_out = os.path.join(out_dir, 'anuga_velocity_max.tif')
    with rasterio.open(dem_utm) as dem_utm_src:
        out_profile = dem_utm_src.profile.copy()
    out_profile.update(dtype='float32')

    with rasterio.open(depth_out, 'w', **out_profile) as dst:
        dst.write(depth_grid.astype(np.float32), 1)
    with rasterio.open(vel_out, 'w', **out_profile) as dst:
        dst.write(velocity_grid.astype(np.float32), 1)

    print(f"  ANUGA outputs saved.")
    print(f"\n===== {site} completed =====\n")

# ---------- Main ----------
if __name__ == '__main__':
    for site, info in SITES.items():
        run_site(site, info['bbox'])
    print("All simulations finished.")
