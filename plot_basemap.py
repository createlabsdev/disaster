import matplotlib.pyplot as plt
import rasterio
from rasterio.plot import show
import contextily as cx
import geopandas as gpd
from shapely.geometry import box
import numpy as np
import sys
import os

site = sys.argv[1] if len(sys.argv) > 1 else 'iiit_kottayam'
tif_path = f'output/{site}/susceptibility_classes.tif'

if not os.path.exists(tif_path):
    print(f"Error: {tif_path} not found.")
    sys.exit(1)

with rasterio.open(tif_path) as src:
    data = src.read(1)
    
    # Replace nodata (or nan) with np.nan for plotting
    data = data.astype(float)
    data[data == src.nodata] = np.nan
    data[data < 1] = np.nan
    
    bounds = src.bounds
    crs = src.crs

# Plot the raster overlaid on a basemap
fig, ax = plt.subplots(figsize=(12, 10))

cmap = plt.cm.RdYlGn_r

show(data, transform=src.transform, ax=ax, cmap=cmap, alpha=0.6, vmin=1, vmax=5)

# Add basemap directly matching the raster's CRS
cx.add_basemap(ax, crs=crs.to_string(), source=cx.providers.OpenStreetMap.Mapnik, attribution_size=8)

# Add title and remove axis
ax.set_title(f'Landslide Susceptibility Map: {site.upper()}\n(Overlaid on OpenStreetMap)', fontsize=14)
ax.set_axis_off()

out_img = f'output/{site}/basemap_plot.png'
plt.savefig(out_img, dpi=200, bbox_inches='tight')
print(f"Saved basemap overlay to {out_img}")
