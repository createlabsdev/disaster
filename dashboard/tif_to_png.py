#!/usr/bin/env python3
"""
tif_to_png.py — Convert single-band probability GeoTIFF to colorized RGBA PNG

Converts a susceptibility map (0.0–1.0 probability) into a QGIS-style
classified overlay with solid green/yellow/orange/red colors.
"""

import numpy as np
import rasterio
from PIL import Image
import os


def probability_to_rgba(prob_array, nodata=-9999.0, absolute=False, multiplier=1.0):
    """
    Convert a 2D probability array (0.0–1.0) to a QGIS-style classified
    RGBA image with solid, clearly visible colors.

    If absolute=False: Uses quantiles to show relative terrain vulnerability.
    If absolute=True: Multiplies prob_array by multiplier and uses fixed thresholds
                      for Active Risk mapping.
    """
    h, w = prob_array.shape
    rgba = np.zeros((h, w, 4), dtype=np.uint8)

    # Nodata mask
    nodata_mask = (prob_array == nodata) | np.isnan(prob_array) | np.isinf(prob_array)

    # Clamp raw probabilities (UNSCALED)
    p_raw = np.clip(prob_array, 0.0, 1.0)
    
    # Calculate quantiles on the UN-SCALED AI data to serve as the absolute baseline for this specific region
    valid_vals_unscaled = p_raw[~nodata_mask]
    if len(valid_vals_unscaled) == 0:
        return rgba

    q20 = float(np.percentile(valid_vals_unscaled, 20))
    q40 = float(np.percentile(valid_vals_unscaled, 40))
    q60 = float(np.percentile(valid_vals_unscaled, 60))
    q80 = float(np.percentile(valid_vals_unscaled, 80))

    # If generating the Active Risk Map, scale the probabilities down by the weather severity.
    # We will compare these scaled probabilities against the unscaled baseline thresholds to accurately
    # reflect real-world storm severity.
    if absolute:
        p = np.clip(p_raw * multiplier, 0.0, 1.0)
    else:
        p = p_raw

    # Class 1: Very Low → Green #2ecc71
    mask1 = (~nodata_mask) & (p <= q20)
    rgba[mask1] = [46, 204, 113, 140]

    # Class 2: Low → Light Green #a8e06c
    mask2 = (~nodata_mask) & (p > q20) & (p <= q40)
    rgba[mask2] = [168, 224, 108, 150]

    # Class 3: Moderate → Yellow #f1c40f
    mask3 = (~nodata_mask) & (p > q40) & (p <= q60)
    rgba[mask3] = [241, 196, 15, 170]

    # Class 4: High → Orange #e67e22
    mask4 = (~nodata_mask) & (p > q60) & (p <= q80)
    rgba[mask4] = [230, 126, 34, 185]

    # Class 5: Very High → Red #e74c3c
    mask5 = (~nodata_mask) & (p > q80)
    rgba[mask5] = [231, 76, 60, 200]

    # Nodata: fully transparent
    rgba[nodata_mask] = [0, 0, 0, 0]

    return rgba


def convert_tif_to_png(tif_path, png_path, nodata=-9999.0, absolute=False, multiplier=1.0):
    """
    Read a single-band GeoTIFF and write a colorized RGBA PNG.

    Parameters
    ----------
    tif_path : str
        Path to the input GeoTIFF (single-band, values 0.0–1.0)
    png_path : str
        Path to the output PNG
    nodata : float
        Nodata value in the GeoTIFF
    absolute : bool
        If true, use fixed color thresholds (0.2, 0.4, 0.6, 0.8)
    multiplier : float
        Multiply raw values by this factor before coloring (used with absolute=True)

    Returns
    -------
    bbox : tuple (west, south, east, north)
        Geographic bounding box of the raster
    """
    with rasterio.open(tif_path) as src:
        prob = src.read(1).astype(np.float32)
        bounds = src.bounds  # (left, bottom, right, top) = (west, south, east, north)
        if src.nodata is not None:
            nodata = src.nodata

    rgba = probability_to_rgba(prob, nodata=nodata, absolute=absolute, multiplier=multiplier)

    img = Image.fromarray(rgba, mode='RGBA')
    os.makedirs(os.path.dirname(png_path) or '.', exist_ok=True)
    img.save(png_path, 'PNG')

    bbox = (bounds.left, bounds.bottom, bounds.right, bounds.top)
    print(f"  Converted {tif_path} → {png_path}")
    print(f"  BBox: {bbox}")
    print(f"  Size: {img.width}x{img.height}")

    return bbox


def compute_overall_risk_score(tif_path, nodata=-9999.0):
    """
    Compute an overall risk score (0–100) from a susceptibility GeoTIFF.

    The score is the 90th percentile of valid pixel probabilities,
    weighted to emphasize high-risk areas.
    """
    with rasterio.open(tif_path) as src:
        prob = src.read(1).astype(np.float32)
        if src.nodata is not None:
            nodata = src.nodata

    valid = prob[(prob != nodata) & ~np.isnan(prob) & (prob >= 0) & (prob <= 1)]

    if len(valid) == 0:
        return 0.0

    # Weighted score: emphasize high-risk pixels
    p90 = float(np.percentile(valid, 90))
    p_max = float(np.max(valid))
    mean_high = float(np.mean(valid[valid > 0.5])) if np.any(valid > 0.5) else 0.0
    fraction_high = float(np.mean(valid > 0.5))

    # Composite score
    score = (p90 * 0.3 + p_max * 0.2 + mean_high * 0.3 + fraction_high * 0.2) * 100
    
    # Confidence score: how certain is the model? 
    # Certainty is highest near 0.0 or 1.0, lowest near 0.5.
    # We map p -> 2 * abs(p - 0.5), yielding [0, 1]. Mean over all pixels * 100
    certainty = 2 * np.abs(valid - 0.5)
    confidence = float(np.mean(certainty)) * 100

    return min(score, 100.0), min(confidence, 100.0)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print("Usage: python tif_to_png.py <input.tif> <output.png>")
        sys.exit(1)
        
    score, conf = compute_overall_risk_score(sys.argv[1])
    print(f"Overall Risk Score: {score:.1f}%")
    print(f"Model Confidence: {conf:.1f}%")
    convert_tif_to_png(sys.argv[1], sys.argv[2])
