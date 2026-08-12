#!/usr/bin/env python3
"""
init_site.py — Dynamic Data Fetcher for the Digital Twin

Downloads SRTM DEM (90m) and ESA WorldCover (10m) for a given bounding box.
"""

import argparse
import os
import sys

try:
    import ee
    import geemap
except ImportError:
    print("ERROR: geemap or earthengine-api not installed.")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Initialize a new site by downloading DEM and Landcover.")
    parser.add_argument("site", help="Name of the site (e.g. munnar)")
    parser.add_argument("west", type=float, help="Min Longitude")
    parser.add_argument("south", type=float, help="Min Latitude")
    parser.add_argument("east", type=float, help="Max Longitude")
    parser.add_argument("north", type=float, help="Max Latitude")
    args = parser.parse_args()

    site = args.site
    bbox = [args.west, args.south, args.east, args.north]
    
    # Init GEE
    try:
        ee.Initialize(project='flood-502410')
    except Exception:
        print("[INFO] Authenticating Earth Engine...")
        ee.Authenticate()
        ee.Initialize(project='flood-502410')

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output', site)
    os.makedirs(out_dir, exist_ok=True)
    
    dem_path = os.path.join(out_dir, 'dem_ref.tif')
    lc_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'{site}_landcover.tif')

    aoi = ee.Geometry.Rectangle(bbox)

    print(f"\n[1/2] Downloading SRTM DEM (90m) for {site}...")
    dem = ee.Image('CGIAR/SRTM90_V4').clip(aoi)
    geemap.ee_export_image(
        dem,
        filename=dem_path,
        scale=90,
        region=aoi,
        crs='EPSG:4326',
        file_per_band=False
    )
    
    print(f"\n[2/2] Downloading ESA WorldCover (10m) for {site}...")
    lc = ee.ImageCollection("ESA/WorldCover/v200").first().clip(aoi)
    geemap.ee_export_image(
        lc,
        filename=lc_path,
        scale=10,
        region=aoi,
        crs='EPSG:4326',
        file_per_band=False
    )
    
    print(f"\n[OK] Site initialized. Data saved to {out_dir}/ and {lc_path}")

if __name__ == "__main__":
    main()
