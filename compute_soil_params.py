#!/usr/bin/env python3
"""Compute K factor, cohesion, phi, and soil depth (landcover+slope lookup) from SoilGrids + DEM + WorldCover."""
import sys, os
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling

site = sys.argv[1]
outdir = f"./output/{site}"
os.makedirs(outdir, exist_ok=True)

# First get profile and shape from DEM
dem_path = os.path.join(outdir, "dem_ref.tif")
if not os.path.exists(dem_path):
    # fallback to old location for backwards compatibility
    dem_path = f"{site}_dem_90m.tif"

with rasterio.open(dem_path) as src:
    profile = src.profile.copy()
    profile.update(dtype='float32', count=1)
    ref_transform = src.transform
    ref_crs = src.crs
    ref_shape = (src.height, src.width)
    dem = src.read(1)

# Check if soil grids exist
indir = f"./soil_{site}"
if os.path.exists(f"{indir}/sand.tif"):
    with rasterio.open(f"{indir}/sand.tif") as s, rasterio.open(f"{indir}/silt.tif") as si, \
         rasterio.open(f"{indir}/clay.tif") as c, rasterio.open(f"{indir}/soc.tif") as sc:
        sand_pct = s.read(1) / 10.0
        silt_pct = si.read(1) / 10.0
        clay_pct = c.read(1) / 10.0
        soc_pct  = sc.read(1) / 10.0
else:
    # Use global averages for tropical hilly regions (Kerala)
    print(f"Warning: SoilGrids not found in {indir}, using regional averages.")
    sand_pct = np.full(ref_shape, 40.0, dtype=np.float32)
    silt_pct = np.full(ref_shape, 30.0, dtype=np.float32)
    clay_pct = np.full(ref_shape, 30.0, dtype=np.float32)
    soc_pct  = np.full(ref_shape, 2.0, dtype=np.float32)

M = (100 - clay_pct) * (silt_pct + sand_pct)
OM = soc_pct * 1.72
K = 0.1317 * (2.1e-4 * (M**1.14) * (12 - OM) + 3.25 * (1 - 2) + 2.5 * (3 - 3)) / 100.0
K = np.clip(K, 0, 1)

# Ensure output goes to output dir
with rasterio.open(f"{outdir}/k_factor.tif", 'w', **profile) as dst:
    dst.write(K.astype(np.float32), 1)

cohesion = np.zeros_like(sand_pct, dtype=np.float32)
phi = np.zeros_like(sand_pct, dtype=np.float32)
for i in range(sand_pct.shape[0]):
    for j in range(sand_pct.shape[1]):
        if clay_pct[i,j] > 35:
            cohesion[i,j] = 5.0; phi[i,j] = 25.0
        elif silt_pct[i,j] > 50:
            cohesion[i,j] = 3.0; phi[i,j] = 30.0
        else:
            cohesion[i,j] = 0.0; phi[i,j] = 35.0

with rasterio.open(f"{outdir}/cohesion.tif", 'w', **profile) as dst:
    dst.write(cohesion, 1)
with rasterio.open(f"{outdir}/phi.tif", 'w', **profile) as dst:
    dst.write(phi, 1)

def resample_to_ref(src_path, resampling=Resampling.bilinear):
    with rasterio.open(src_path) as src:
        dst = np.empty(ref_shape, dtype=np.float32)
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=ref_transform,
            dst_crs=ref_crs,
            resampling=resampling
        )
    return dst

landcover = resample_to_ref(f"{site}_landcover.tif", resampling=Resampling.nearest)

px = abs(ref_transform.a)
py = abs(ref_transform.e)
dzdx = np.gradient(dem, px, axis=1)
dzdy = np.gradient(dem, py, axis=0)
slope_deg = np.degrees(np.arctan(np.sqrt(dzdx**2 + dzdy**2)))

steep = slope_deg >= 15.0

soil_depth = np.select(
    condlist=[
        landcover == 80,
        landcover == 90,
        landcover == 50,
        landcover == 60,
        (landcover == 10) & steep,
        (landcover == 10) & ~steep,
        landcover == 40,
        landcover == 20,
        landcover == 30,
        landcover == 95,
    ],
    choicelist=[0.0, 0.1, 0.3, 0.2, 1.1, 1.7, 1.4, 0.8, 0.7, 0.3],
    default=1.0
).astype(np.float32)

c_factor = np.select(
    condlist=[
        landcover == 10, landcover == 20, landcover == 30, landcover == 40,
        landcover == 50, landcover == 60, landcover == 80, landcover == 90,
        landcover == 95
    ],
    choicelist=[0.01, 0.05, 0.10, 0.30, 0.05, 1.00, 0.00, 0.00, 0.01],
    default=0.10
).astype(np.float32)

manning_n = np.select(
    condlist=[
        landcover == 10, landcover == 20, landcover == 30, landcover == 40,
        landcover == 50, landcover == 60, landcover == 80, landcover == 90,
        landcover == 95
    ],
    choicelist=[0.10, 0.07, 0.035, 0.04, 0.015, 0.02, 0.03, 0.05, 0.10],
    default=0.035
).astype(np.float32)

r_factor = np.full(ref_shape, 3000.0, dtype=np.float32)

with rasterio.open(f"{outdir}/soil_depth.tif", 'w', **profile) as dst:
    dst.write(soil_depth, 1)
with rasterio.open(f"{outdir}/c_factor.tif", 'w', **profile) as dst:
    dst.write(c_factor, 1)
with rasterio.open(f"{outdir}/manning_n.tif", 'w', **profile) as dst:
    dst.write(manning_n, 1)
with rasterio.open(f"{outdir}/r_factor.tif", 'w', **profile) as dst:
    dst.write(r_factor, 1)

print(f"All physical parameters generated and written to {outdir}/")
