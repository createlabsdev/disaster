#!/usr/bin/env python3
"""Calculate R factor from WorldClim monthly precipitation (no login).
Usage: python get_rfactor.py <site_name> <west> <south> <east> <north>"""

import sys, os
import numpy as np
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling

if len(sys.argv) != 6:
    print("Usage: python get_rfactor.py site_name west south east north")
    sys.exit(1)

site = sys.argv[1]
west, south, east, north = map(float, sys.argv[2:])

# WorldClim monthly precipitation base URL (30s = ~1km)
base_url = "/vsicurl/https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_prec_"

# We'll fetch each month, clip to bbox, and accumulate R factor components.
# Because the global files are large (~1GB each), we use gdal_translate with projwin to grab only our area.
# This is done via a temporary VRT and then read into memory.

months = range(1, 13)
annual_R = None
profile = None
annual_sum = None

for m in months:
    url = f"{base_url}{m:02d}.tif"
    print(f"Processing month {m}...")
    # Use gdal_translate to extract our bbox
    temp_file = f"/tmp/prec_{site}_{m}.tif"
    cmd = (f"gdal_translate -q -projwin {west} {north} {east} {south} "
           f"-projwin_srs EPSG:4326 {url} {temp_file}")
    os.system(cmd)
    if not os.path.exists(temp_file):
        print(f"  Failed to download month {m}")
        continue
    with rasterio.open(temp_file) as src:
        prec = src.read(1)
        if profile is None:
            profile = src.profile
        if annual_sum is None:
            annual_sum = np.zeros_like(prec, dtype=np.float32)
        annual_sum += prec

    # R factor contribution: 0.264 * p_i^2 / P (but P is unknown until end)
    # We'll accumulate p_i^2 and compute final R after the loop.
    if annual_R is None:
        annual_R = np.zeros_like(prec, dtype=np.float32)
    annual_R += (prec ** 2)

    os.remove(temp_file)

# After all months, annual_sum holds total precip P (mm)
# Final R = Σ (0.264 * p_i^2) / P
R = 0.264 * annual_R / annual_sum   # units MJ·mm/(ha·h·yr)
R[np.isnan(R) | np.isinf(R)] = 0

# Write R factor raster
profile.update(dtype='float32', driver='GTiff')
with rasterio.open(f"{site}_rfactor.tif", 'w', **profile) as dst:
    dst.write(R.astype(np.float32), 1)
print(f"R factor saved to {site}_rfactor.tif")
