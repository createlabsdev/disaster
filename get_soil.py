#!/usr/bin/env python3
"""Download SoilGrids 2.0 layers and compute K factor, soil depth, cohesion, phi.
No login required. Works with SoilGrids 2020 WCS 2.0.1"""

import sys, os
import rasterio
import numpy as np
from owslib.wcs import WebCoverageService

if len(sys.argv) != 6:
    print("Usage: python get_soil.py site_name west south east north")
    sys.exit(1)

site = sys.argv[1]
west, south, east, north = map(float, sys.argv[2:])

out_dir = f"./soil_{site}"
os.makedirs(out_dir, exist_ok=True)

# Correct SoilGrids 2.0 WCS endpoint (no login)
url = "https://maps.isric.org/mapserv?map=/map/soilgrids2020.map"
wcs = WebCoverageService(url, version='2.0.1')

# Coverage IDs for 0–5 cm depth (topsoil)
layers = {
    'sand': 'sand_0-5cm_mean',
    'silt': 'silt_0-5cm_mean',
    'clay': 'clay_0-5cm_mean',
    'soc':  'soc_0-5cm_mean',      # soil organic carbon (g/kg)
    'depth':'bdticm',              # absolute depth to bedrock (cm)
}

# Download each layer as GeoTIFF (subset by bbox)
for var, layer in layers.items():
    print(f"Fetching {var}...")
    response = wcs.getCoverage(
        identifier=[layer],
        bbox=(west, south, east, north),
        crs='urn:ogc:def:crs:EPSG::4326',
        format='image/tiff',          # GeoTIFF
        width=40, height=40,           # rough pixel count (will be resampled)
        timeout=120
    )
    fname = f"{out_dir}/{site}_{var}.tif"
    with open(fname, 'wb') as f:
        f.write(response.read())

# Verify we got valid GeoTIFFs
sand_tif = f"{out_dir}/{site}_sand.tif"
try:
    with rasterio.open(sand_tif) as src:
        sand = src.read(1)
except:
    print("Error: downloaded sand.tif is not a valid GeoTIFF. "
          "The WCS request may have failed. Check the file content with "
          f"'head {sand_tif}'.")
    sys.exit(1)

# Read all raw rasters (units: % × 10 for sand/silt/clay, g/kg × 10 for soc, cm for depth)
with rasterio.open(sand_tif) as s, rasterio.open(f"{out_dir}/{site}_silt.tif") as si, \
     rasterio.open(f"{out_dir}/{site}_clay.tif") as c, rasterio.open(f"{out_dir}/{site}_soc.tif") as sc, \
     rasterio.open(f"{out_dir}/{site}_depth.tif") as d:

    sand_pct = s.read(1) / 10.0   # convert to percent (0–100)
    silt_pct = si.read(1) / 10.0
    clay_pct = c.read(1) / 10.0
    soc_pct  = sc.read(1) / 10.0   # g/kg -> %
    depth_m  = d.read(1) / 100.0   # cm -> m

    profile = s.profile
    profile.update(dtype='float32')

# ---------- K factor (Renard et al. 1997) ----------
M = (100 - clay_pct) * (silt_pct + sand_pct)   # M factor
OM = soc_pct * 1.72                             # organic matter (%)
# structure code (1 = fine granular, typical for silt/clay) and permeability code (3 = moderate)
s_code = 1
p_code = 3
K = 0.1317 * (2.1e-4 * (M**1.14) * (12 - OM) + 3.25 * (s_code - 2) + 2.5 * (p_code - 3)) / 100.0
K = np.clip(K, 0, 1)

with rasterio.open(f"{out_dir}/{site}_k_factor.tif", 'w', **profile) as dst:
    dst.write(K.astype(np.float32), 1)

# ---------- Cohesion and friction angle (simple texture lookup) ----------
cohesion = np.zeros_like(sand_pct, dtype=np.float32)
phi      = np.zeros_like(sand_pct, dtype=np.float32)
for i in range(sand_pct.shape[0]):
    for j in range(sand_pct.shape[1]):
        if clay_pct[i,j] > 35:
            cohesion[i,j] = 5.0   # kPa
            phi[i,j]      = 25.0  # degrees
        elif silt_pct[i,j] > 50:
            cohesion[i,j] = 3.0
            phi[i,j]      = 30.0
        else:
            cohesion[i,j] = 0.0
            phi[i,j]      = 35.0

with rasterio.open(f"{out_dir}/{site}_cohesion.tif", 'w', **profile) as dst:
    dst.write(cohesion, 1)
with rasterio.open(f"{out_dir}/{site}_phi.tif", 'w', **profile) as dst:
    dst.write(phi, 1)

# Save depth (already in metres)
with rasterio.open(f"{out_dir}/{site}_soil_depth.tif", 'w', **profile) as dst:
    dst.write(depth_m.astype(np.float32), 1)

print(f"Soil layers saved in {out_dir}/")
