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
from scipy.ndimage import label


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


def probability_to_siltation_rgba(prob_array, nodata=-9999.0, manning_path=None):
    """
    Convert a 2D probability array to a River Siltation & Sediment overlay.
    STRICTLY masks out land, displaying sediment colors ONLY across river channels, streams, lakes, and reservoirs.
    """
    h, w = prob_array.shape
    rgba = np.zeros((h, w, 4), dtype=np.uint8)

    nodata_mask = (prob_array == nodata) | np.isnan(prob_array) | np.isinf(prob_array)
    p_raw = np.clip(prob_array, 0.0, 1.0)

    # Extract exact water channel & riverbed mask using Manning roughness or top flow channels
    water_mask = None
    if manning_path and os.path.exists(manning_path):
        try:
            with rasterio.open(manning_path) as src:
                mn = src.read(1).astype(np.float32)
                # Manning n <= 0.020 represents open water channels & rivers
                raw_water = (mn <= 0.020) & (mn > 0)
                
                # Morphological component filter: eliminate isolated town noise dots (< 8 connected pixels)
                labeled, num_features = label(raw_water)
                sizes = np.bincount(labeled.ravel())
                water_mask = np.zeros_like(raw_water, dtype=bool)
                for i in range(1, num_features + 1):
                    if sizes[i] >= 8:
                        water_mask |= (labeled == i)
        except Exception:
            water_mask = None

    if water_mask is None or np.sum(water_mask) == 0:
        valid_vals = p_raw[~nodata_mask]
        if len(valid_vals) == 0:
            return rgba
        # Fallback water channel approximation (top 5% valley channels)
        q95 = float(np.percentile(valid_vals, 95))
        raw_water = (p_raw >= q95)
        labeled, num_features = label(raw_water)
        sizes = np.bincount(labeled.ravel())
        water_mask = np.zeros_like(raw_water, dtype=bool)
        for i in range(1, num_features + 1):
            if sizes[i] >= 8:
                water_mask |= (labeled == i)

    # Strictly mask out all LAND pixels (100% Transparent)
    combined_mask = (~nodata_mask) & water_mask
    valid_water_vals = p_raw[combined_mask]

    if len(valid_water_vals) == 0:
        return rgba

    q33 = float(np.percentile(valid_water_vals, 33))
    q66 = float(np.percentile(valid_water_vals, 66))

    # Class 1: Low Siltation Deposit → Bright Gold #f1c40f
    mask1 = combined_mask & (p_raw <= q33)
    rgba[mask1] = [241, 196, 15, 180]

    # Class 2: Moderate Siltation Deposit → Burnt Orange #e67e22
    mask2 = combined_mask & (p_raw > q33) & (p_raw <= q66)
    rgba[mask2] = [230, 126, 34, 210]

    # Class 3: Heavy Silt & Mud Accumulation → Deep Crimson #c0392b
    mask3 = combined_mask & (p_raw > q66)
    rgba[mask3] = [192, 57, 43, 235]

    # Land / Nodata: fully transparent
    rgba[~combined_mask] = [0, 0, 0, 0]

    return rgba


def convert_tif_to_png(tif_path, png_path, nodata=-9999.0, absolute=False, multiplier=1.0, siltation=False):
    """
    Read a single-band GeoTIFF and write a colorized RGBA PNG.
    """
    with rasterio.open(tif_path) as src:
        prob = src.read(1).astype(np.float32)
        bounds = src.bounds  # (left, bottom, right, top) = (west, south, east, north)
        if src.nodata is not None:
            nodata = src.nodata

    if siltation:
        manning_path = os.path.join(os.path.dirname(tif_path), 'manning_n.tif')
        rgba = probability_to_siltation_rgba(prob, nodata=nodata, manning_path=manning_path)
    else:
        rgba = probability_to_rgba(prob, nodata=nodata, absolute=absolute, multiplier=multiplier)

    img = Image.fromarray(rgba, mode='RGBA')
    os.makedirs(os.path.dirname(png_path) or '.', exist_ok=True)
    img.save(png_path, 'PNG')

    bbox = (bounds.left, bounds.bottom, bounds.right, bounds.top)
    print(f"  Converted {tif_path} → {png_path}")
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
