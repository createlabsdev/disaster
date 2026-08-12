#!/usr/bin/env python3
"""
Reclassify ESA WorldCover 10m (2021) to RUSLE C factor and Manning's n.
Usage: python reclassify_worldcover.py <input_tif> <output_c_factor> <output_manning>
"""
import sys
import rasterio
import numpy as np

# Lookup table: class -> (C factor, Manning's n)
lookup = {
    10: (0.001, 0.100),   # Tree cover
    20: (0.01,  0.040),   # Shrubland
    30: (0.01,  0.040),   # Grassland
    40: (0.2,   0.035),   # Cropland
    50: (0.9,   0.015),   # Built-up
    60: (1.0,   0.025),   # Bare/sparse vegetation
    70: (0,     0.010),   # Snow and ice
    80: (0,     0.015),   # Permanent water bodies
    90: (0.05,  0.050),   # Herbaceous wetland
    95: (0.001, 0.080),   # Mangroves
    100:(0,     0.050)    # Moss and lichen
}

if len(sys.argv) != 4:
    print("Usage: python reclassify_worldcover.py input.tif c_factor.tif manning_n.tif")
    sys.exit(1)

input_tif = sys.argv[1]
c_out = sys.argv[2]
n_out = sys.argv[3]

with rasterio.open(input_tif) as src:
    lulc = src.read(1).astype(int)
    profile = src.profile
    profile.update(dtype='float32', nodata=-9999)

    c_factor = np.zeros_like(lulc, dtype=np.float32)
    manning = np.zeros_like(lulc, dtype=np.float32)

    for cls, (cf, mn) in lookup.items():
        mask = lulc == cls
        c_factor[mask] = cf
        manning[mask] = mn

    with rasterio.open(c_out, 'w', **profile) as dst:
        dst.write(c_factor, 1)

    with rasterio.open(n_out, 'w', **profile) as dst:
        dst.write(manning, 1)

    print(f"Done. Created {c_out} and {n_out}")
