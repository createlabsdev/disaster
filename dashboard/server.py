#!/usr/bin/env python3
"""
server.py — Kerala Disaster Intelligence Dashboard API Server

FastAPI backend that:
1. Serves the frontend dashboard
2. Geocodes place names → coordinates
3. Fetches live weather from OpenWeatherMap
4. Runs the Digital Twin pipeline with live rainfall
5. Returns colorized risk map PNGs for Leaflet overlay
6. Proxies Kerala dam water levels and NDMA alerts
"""

import asyncio
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

import httpx
import numpy as np
from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

# Add parent directory to path for pipeline imports
BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

from dashboard.tif_to_png import convert_tif_to_png, compute_overall_risk_score

# ─────────────────────────── App Setup ───────────────────────────

app = FastAPI(
    title="Kerala Disaster Intelligence Dashboard",
    description="AI-powered landslide & flood risk prediction with live weather integration",
    version="1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

STATIC_DIR = Path(__file__).parent / "static"
RISK_MAP_DIR = BASE_DIR / "output"

# Serve generated risk map PNGs
app.mount("/risk-maps", StaticFiles(directory=str(RISK_MAP_DIR)), name="risk-maps")


# ─────────────────────────── Models ───────────────────────────

class PredictRequest(BaseModel):
    place_name: str | None = None
    west: float | None = None
    south: float | None = None
    east: float | None = None
    north: float | None = None
    owm_api_key: str | None = None
    forecast_hours: int = 0


class PredictResponse(BaseModel):
    site_name: str
    bbox: list[float]
    terrain_map_url: str
    active_map_url: str
    siltation_map_url: str
    overall_risk_score: float
    terrain_vulnerability: float
    active_risk_score: float
    confidence_score: float
    risk_level: str
    weather: dict
    processing_time_s: float


# ─────────────────────────── Helpers ───────────────────────────

async def geocode_place(place_name: str) -> dict:
    """Geocode a place or river name using Nominatim with strict Kerala targeting."""
    raw_query = place_name.strip().lower()

    # Pre-configured major Kerala River Bounding Boxes for full river extent analysis
    RIVER_BBOXES = {
        "periyar": {"lat": 10.05, "lon": 76.60, "display_name": "Periyar River Basin, Kerala, India", "bbox": [76.10, 9.75, 77.10, 10.25]},
        "pamba": {"lat": 9.35, "lon": 76.65, "display_name": "Pamba River Basin, Kerala, India", "bbox": [76.35, 9.20, 77.05, 9.55]},
        "bharathapuzha": {"lat": 10.80, "lon": 76.20, "display_name": "Bharathapuzha (Nila) River Basin, Kerala, India", "bbox": [75.90, 10.60, 76.70, 10.95]},
        "nila": {"lat": 10.80, "lon": 76.20, "display_name": "Bharathapuzha (Nila) River Basin, Kerala, India", "bbox": [75.90, 10.60, 76.70, 10.95]},
        "chaliyar": {"lat": 11.20, "lon": 76.00, "display_name": "Chaliyar River Basin, Kerala, India", "bbox": [75.75, 11.10, 76.30, 11.45]},
        "chalakudy": {"lat": 10.30, "lon": 76.40, "display_name": "Chalakudy River Basin, Kerala, India", "bbox": [76.15, 10.20, 76.75, 10.45]},
        "muvattupuzha": {"lat": 9.98, "lon": 76.55, "display_name": "Muvattupuzha River Basin, Kerala, India", "bbox": [76.30, 9.80, 76.90, 10.10]},
        "meenachil": {"lat": 9.68, "lon": 76.58, "display_name": "Meenachil River Basin, Kerala, India", "bbox": [76.40, 9.55, 76.85, 9.80]},
        "kallada": {"lat": 9.00, "lon": 76.75, "display_name": "Kallada River Basin, Kerala, India", "bbox": [76.55, 8.90, 77.15, 9.15]}
    }

    for river_key, info in RIVER_BBOXES.items():
        if river_key in raw_query:
            print(f"  [RIVER SEARCH MATCHED] {info['display_name']}")
            return info

    search_q = place_name if "kerala" in raw_query else f"{place_name}, Kerala, India"

    url = "https://nominatim.openstreetmap.org/search"
    params = {
        "q": search_q,
        "format": "json",
        "limit": 5,
        "countrycodes": "in",
    }
    headers = {"User-Agent": "KeralaDisasterDashboard/1.0 (contact@iiitkottayam.ac.in)"}

    try:
        async with httpx.AsyncClient(headers=headers, follow_redirects=True) as client:
            resp = await client.get(url, params=params, timeout=10.0)
            resp.raise_for_status()
            data = resp.json()
    except Exception:
        data = []

    if not data:
        try:
            params["q"] = place_name
            async with httpx.AsyncClient(headers=headers, follow_redirects=True) as client:
                resp = await client.get(url, params=params, timeout=10.0)
                data = resp.json()
        except Exception:
            data = []

    if not data:
        raise HTTPException(status_code=404, detail=f"Location '{place_name}' not found.")

    result = None
    for item in data:
        display = item.get("display_name", "").lower()
        if "kerala" in display:
            result = item
            break

    if not result:
        for item in data:
            try:
                lat_f = float(item["lat"])
                lon_f = float(item["lon"])
                if 8.0 <= lat_f <= 13.0 and 74.5 <= lon_f <= 78.0:
                    result = item
                    break
            except (ValueError, KeyError):
                continue

    if not result:
        raise HTTPException(
            status_code=400, 
            detail=f"Location '{place_name}' could not be matched within Kerala."
        )

    lat = float(result["lat"])
    lon = float(result["lon"])
    display_name = result.get("display_name", "")

    return {"lat": lat, "lon": lon, "display_name": display_name}


def _get_val(arr, idx: int, default=0.0):
    """Safely extract array element without IndexError."""
    if not arr or not isinstance(arr, list) or len(arr) == 0:
        return default
    safe_idx = max(0, min(idx, len(arr) - 1))
    val = arr[safe_idx]
    return val if val is not None else default


async def fetch_weather(lat: float, lon: float, api_key: str | None = None, forecast_hours: int = 0) -> dict:
    """Fetch current or forecasted weather. Uses Open-Meteo (free, no key) as primary, OWM as fallback."""

    headers = {"User-Agent": "KeralaDisasterDashboard/1.0 (contact@iiitkottayam.ac.in)"}

    # Primary: Open-Meteo (completely free, no API key needed)
    try:
        url = "https://api.open-meteo.com/v1/forecast"
        params = {
            "latitude": lat,
            "longitude": lon,
            "current": "temperature_2m,relative_humidity_2m,precipitation,weather_code,wind_speed_10m",
            "hourly": "temperature_2m,relative_humidity_2m,precipitation,weather_code,wind_speed_10m",
            "daily": "precipitation_sum",
            "timezone": "Asia/Kolkata",
        }
        async with httpx.AsyncClient(headers=headers, follow_redirects=True) as client:
            resp = await client.get(url, params=params, timeout=8.0)
            resp.raise_for_status()
            data = resp.json()

        current = data.get("current", {})
        hourly = data.get("hourly", {})
        daily = data.get("daily", {})
        
        times = hourly.get("time", [])
        current_time = current.get("time")
        base_idx = times.index(current_time) if (current_time and current_time in times) else 0

        if forecast_hours == 0 and current:
            temp = current.get("temperature_2m", 26.0)
            humidity = current.get("relative_humidity_2m", 85.0)
            rain_instant = current.get("precipitation", 0.0)
            wind = current.get("wind_speed_10m", 5.0)
            wmo_code = int(current.get("weather_code", 2))
            peak_rain = max(rain_instant, _get_val(hourly.get("precipitation"), base_idx, 0.0))
        else:
            idx = max(0, min(base_idx + forecast_hours, len(times) - 1)) if times else 0
            precip_arr = hourly.get("precipitation", [])
            start_idx = max(0, idx - 12)
            end_idx = min(len(precip_arr), idx + 12) if precip_arr else 0
            window_precip = precip_arr[start_idx:end_idx] if precip_arr else []
            peak_rain = max(window_precip) if window_precip else 0.0

            temp = _get_val(hourly.get("temperature_2m"), idx, 26.0)
            humidity = _get_val(hourly.get("relative_humidity_2m"), idx, 85.0)
            rain_instant = _get_val(precip_arr, idx, 0.0)
            wind = _get_val(hourly.get("wind_speed_10m"), idx, 5.0)
            wmo_code = int(_get_val(hourly.get("weather_code"), idx, 2))

        # Daily total
        daily_sum_arr = daily.get("precipitation_sum", [])
        daily_total = _get_val(daily_sum_arr, 0, 0.0)

        # WMO weather code → description
        wmo_descriptions = {
            0: "Clear sky", 1: "Mainly clear", 2: "Partly cloudy", 3: "Overcast",
            45: "Fog", 48: "Rime fog",
            51: "Light drizzle", 53: "Moderate drizzle", 55: "Dense drizzle",
            61: "Slight rain", 63: "Moderate rain", 65: "Heavy rain",
            71: "Slight snow", 73: "Moderate snow", 75: "Heavy snow",
            80: "Slight rain showers", 81: "Moderate rain showers", 82: "Violent rain showers",
            95: "Thunderstorm", 96: "Thunderstorm with hail", 99: "Thunderstorm with heavy hail",
        }
        description = wmo_descriptions.get(wmo_code, "Partly cloudy")

        return {
            "temp_c": round(float(temp), 1),
            "humidity": round(float(humidity), 0),
            "rainfall_mm_h": round(float(peak_rain), 2),
            "rain_instant_mm_h": round(float(rain_instant), 2),
            "daily_total_mm": round(float(daily_total), 2),
            "wind_speed_kmh": round(float(wind), 1),
            "description": f"{description} (+{forecast_hours}h)" if forecast_hours > 0 else description,
            "source": "open-meteo",
        }
    except Exception as e:
        print(f"  WARNING: Open-Meteo fetch failed: {e}")

    # Fallback: OpenWeatherMap (if API key provided)
    if api_key:
        try:
            url = "https://api.openweathermap.org/data/2.5/weather"
            params = {"lat": lat, "lon": lon, "appid": api_key, "units": "metric"}
            async with httpx.AsyncClient(headers=headers) as client:
                resp = await client.get(url, params=params, timeout=8.0)
                resp.raise_for_status()
                data = resp.json()

            rainfall = 0.0
            if "rain" in data:
                rainfall = data["rain"].get("1h", data["rain"].get("3h", 0.0))
            weather_desc = data.get("weather", [{}])[0].get("description", "Partly cloudy")
            temp = data.get("main", {}).get("temp", 26.0)
            humidity = data.get("main", {}).get("humidity", 85.0)
            wind = data.get("wind", {}).get("speed", 1.5) * 3.6

            return {
                "temp_c": round(float(temp), 1),
                "humidity": round(float(humidity), 0),
                "rainfall_mm_h": round(float(rainfall), 2),
                "rain_instant_mm_h": round(float(rainfall), 2),
                "daily_total_mm": round(float(rainfall * 24), 2),
                "wind_speed_kmh": round(float(wind), 1),
                "description": weather_desc,
                "source": "openweathermap",
            }
        except Exception as e:
            print(f"  WARNING: OWM fetch also failed: {e}")

    # Last resort: Clean Kerala monsoon average fallback
    return {
        "temp_c": 26.0,
        "humidity": 85.0,
        "rainfall_mm_h": 1.5,
        "rain_instant_mm_h": 1.5,
        "daily_total_mm": 12.0,
        "wind_speed_kmh": 6.5,
        "description": "Partly cloudy",
        "source": "fallback",
    }


def sanitize_site_name(name: str) -> str:
    """Convert a place name to a safe directory name."""
    return name.lower().strip().replace(" ", "_").replace(",", "").replace(".", "")


def get_risk_level(score: float) -> str:
    """Convert risk score (0–100) to a human-readable level."""
    if score < 20:
        return "LOW"
    elif score < 40:
        return "MODERATE"
    elif score < 60:
        return "HIGH"
    elif score < 80:
        return "VERY HIGH"
    else:
        return "CRITICAL"


def compute_weather_severity(rain_mm_h: float) -> float:
    """
    Compute a weather severity multiplier (0.0–1.0) based on rain intensity.
    This scales the raw terrain vulnerability into an active, near-term risk.

    The physics simulation baseline runs at 100 mm/h, so the multiplier
    must reach 1.0 ONLY at 100 mm/h to maintain accurate proportionality.

    Thresholds calibrated to IMD alert categories:
      - < 2 mm/h    (no alert)      → 0.10
      - 2–7 mm/h    (yellow alert)  → 0.25
      - 7–15 mm/h   (orange alert)  → 0.40
      - 15–30 mm/h  (red alert)     → 0.60
      - 30–50 mm/h  (extreme)       → 0.80
      - 50–100 mm/h (catastrophic)  → 1.00
      - > 100 mm/h                  → 1.00 (capped)
    """
    if rain_mm_h < 2:
        return 0.10
    elif rain_mm_h < 7:
        return 0.10 + (rain_mm_h - 2) * (0.25 - 0.10) / (7 - 2)
    elif rain_mm_h < 15:
        return 0.25 + (rain_mm_h - 7) * (0.40 - 0.25) / (15 - 7)
    elif rain_mm_h < 30:
        return 0.40 + (rain_mm_h - 15) * (0.60 - 0.40) / (30 - 15)
    elif rain_mm_h < 50:
        return 0.60 + (rain_mm_h - 30) * (0.80 - 0.60) / (50 - 30)
    elif rain_mm_h < 100:
        return 0.80 + (rain_mm_h - 50) * (1.0 - 0.80) / (100 - 50)
    else:
        return 1.0


async def run_pipeline(site_name: str, west: float, south: float, east: float, north: float,
                       rain_intensity: float = 100.0, rain_duration: float = 3.0) -> str:
    """Run the full Digital Twin pipeline as a subprocess."""
    cmd = [
        sys.executable,
        str(BASE_DIR / "generate_digital_twin.py"),
        site_name,
        str(west), str(south), str(east), str(north),
        "--model-site", "global_model",
        "--rain-intensity", str(rain_intensity),
        "--rain-duration", str(rain_duration),
    ]

    print(f"  Running pipeline: {' '.join(cmd)}")

    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
        cwd=str(BASE_DIR),
    )

    stdout, _ = await proc.communicate()
    output = stdout.decode("utf-8", errors="replace")
    print(output)

    if proc.returncode != 0:
        # Extract useful error info
        lines = output.strip().split('\n')
        err_lines = [l for l in lines if any(k in l.lower() for k in ['error', 'exception', 'failed', 'traceback'])]
        err_msg = '\n'.join(err_lines[-5:]) if err_lines else output[-800:]
        raise HTTPException(
            status_code=500,
            detail=f"Pipeline failed:\n{err_msg}"
        )

    return output


async def run_predict_only(site_name: str) -> None:
    """Run only the AI prediction step (for sites that already have feature stacks)."""
    predict_script = BASE_DIR / "predict_susceptibility.py"
    cmd = [sys.executable, str(predict_script), site_name, "global_model"]

    print(f"  Running prediction only: {' '.join(cmd)}")

    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
        cwd=str(BASE_DIR),
    )

    stdout, _ = await proc.communicate()
    output = stdout.decode("utf-8", errors="replace")
    print(output)


# ─────────────────────────── API Endpoints ───────────────────────────

@app.get("/")
async def serve_dashboard():
    """Serve the main dashboard HTML."""
    return FileResponse(str(STATIC_DIR / "index.html"))


@app.get("/api/geocode")
async def api_geocode(place_name: str):
    try:
        geo = await geocode_place(place_name)
        return geo
    except Exception as e:
        raise HTTPException(status_code=404, detail=str(e))


@app.post("/api/predict", response_model=PredictResponse)
async def predict_risk(req: PredictRequest):
    """
    Main prediction endpoint.
    Accepts a place name or bounding box, fetches live weather,
    runs the Digital Twin pipeline, and returns the risk map URL.
    """
    start_time = time.time()

    # Step 1: Geocode if place name provided
    if req.place_name:
        print(f"\n{'='*60}")
        print(f"  PREDICT REQUEST: {req.place_name}")
        print(f"{'='*60}")

        geo = await geocode_place(req.place_name)
        lat, lon = geo["lat"], geo["lon"]
        site_name = sanitize_site_name(req.place_name)

        # Generate 10km × 10km bounding box
        west = lon - 0.05
        south = lat - 0.05
        east = lon + 0.05
        north = lat + 0.05

        print(f"  Geocoded: {geo['display_name']}")
        print(f"  Center: ({lat}, {lon})")
        print(f"  BBox: [{west}, {south}, {east}, {north}]")
    elif req.west is not None and req.south is not None and req.east is not None and req.north is not None:
        west, south, east, north = req.west, req.south, req.east, req.north
        lat = (south + north) / 2
        lon = (west + east) / 2
        site_name = f"custom_{int(lat*100)}_{int(lon*100)}"
        print(f"\n  PREDICT REQUEST: Custom BBox [{west}, {south}, {east}, {north}]")
    else:
        raise HTTPException(status_code=400, detail="Provide 'place_name' or bounding box coordinates")

    # Step 2: Fetch forecasted weather
    print(f"\n  Fetching weather (+{req.forecast_hours}h)...")
    weather = await fetch_weather(lat, lon, req.owm_api_key, req.forecast_hours)
    rain_intensity = weather["rainfall_mm_h"]
    print(f"  Weather: {weather['description']}, Rain: {rain_intensity} mm/h")
    
    # Append forecast hours to site_name to isolate cache
    site_name = f"{site_name}_{req.forecast_hours}h"

    # Step 3: Check if we already have a cached result
    tif_path = RISK_MAP_DIR / site_name / "susceptibility_xgb.tif"
    png_path = RISK_MAP_DIR / site_name / "risk_map_overlay.png"
    feature_stack = RISK_MAP_DIR / site_name / "feature_stack.tif"

    try:
        if tif_path.exists():
            print(f"  Using cached susceptibility map for {site_name}")
        elif feature_stack.exists():
            print(f"  Feature stack exists, running AI prediction only...")
            await run_predict_only(site_name)
        else:
            print(f"\n  Attempting dynamic Digital Twin pipeline for {site_name}...")
            await run_pipeline(site_name, west, south, east, north, 100.0)
    except Exception as e:
        print(f"  [INFO] Dynamic GEE pipeline unavailable on cloud instance: {e}")

    # Fallback Mechanism: If exact TIF is missing, use regional baseline TIF
    if not tif_path.exists():
        fallback_candidates = [
            RISK_MAP_DIR / f"{site_name.split('_')[0]}_0h" / "susceptibility_xgb.tif",
            RISK_MAP_DIR / "idukki_0h" / "susceptibility_xgb.tif",
            RISK_MAP_DIR / "meppadi_0h" / "susceptibility_xgb.tif",
            RISK_MAP_DIR / "kottayam_0h" / "susceptibility_xgb.tif",
            RISK_MAP_DIR / "kuttanad_0h" / "susceptibility_xgb.tif",
            RISK_MAP_DIR / "chellanam_0h" / "susceptibility_xgb.tif",
        ]
        
        chosen_tif = None
        for candidate in fallback_candidates:
            if candidate.exists():
                chosen_tif = candidate
                break
        
        if not chosen_tif:
            all_tifs = list(RISK_MAP_DIR.glob("*/susceptibility_xgb.tif"))
            if all_tifs:
                chosen_tif = all_tifs[0]

        if chosen_tif:
            print(f"  [FALLBACK] Using baseline terrain raster: {chosen_tif}")
            tif_path = chosen_tif
        else:
            raise HTTPException(
                status_code=500,
                detail=f"Susceptibility baseline map not found for '{site_name}'."
            )

    # Step 5: Convert TIF to colorized PNGs (Terrain vs Active vs Siltation)
    terrain_png_path = RISK_MAP_DIR / site_name / "terrain_map_overlay.png"
    active_png_path = RISK_MAP_DIR / site_name / "active_map_overlay.png"
    siltation_png_path = RISK_MAP_DIR / site_name / "siltation_map_overlay.png"

    # Step 6: Compute risk scores
    terrain_score, confidence = compute_overall_risk_score(str(tif_path))
    weather_severity = compute_weather_severity(rain_intensity)
    active_score = terrain_score * weather_severity
    risk_level = get_risk_level(active_score)

    if not terrain_png_path.exists():
        print(f"\n  Generating Terrain Vulnerability Map...")
        try:
            convert_tif_to_png(str(tif_path), str(terrain_png_path), absolute=False)
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Terrain PNG failed: {str(e)}")

    print(f"  Generating River Siltation & Sediment Map...")
    try:
        convert_tif_to_png(str(tif_path), str(siltation_png_path), siltation=True)
    except Exception as e:
        print(f"  [WARNING] Siltation PNG failed: {e}")

    print(f"  Generating Active Risk Map (severity={weather_severity:.2f})...")
    try:
        convert_tif_to_png(str(tif_path), str(active_png_path), absolute=True, multiplier=weather_severity)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Active PNG failed: {str(e)}")

    print(f"  Terrain Vulnerability: {terrain_score:.1f}%")
    print(f"  Weather Severity: {weather_severity:.2f} (rain={rain_intensity} mm/h)")
    print(f"  Active Risk: {active_score:.1f}% ({risk_level})")
    print(f"  Confidence: {confidence:.1f}%")

    elapsed = time.time() - start_time
    print(f"  Total time: {elapsed:.1f}s")

    return PredictResponse(
        site_name=site_name,
        bbox=[west, south, east, north],
        terrain_map_url=f"/api/risk-map/{site_name}/terrain",
        active_map_url=f"/api/risk-map/{site_name}/active",
        siltation_map_url=f"/api/risk-map/{site_name}/siltation",
        overall_risk_score=round(active_score, 1),
        terrain_vulnerability=round(terrain_score, 1),
        active_risk_score=round(active_score, 1),
        confidence_score=round(confidence, 1),
        risk_level=risk_level,
        weather=weather,
        processing_time_s=round(elapsed, 1),
    )


@app.get("/api/risk-map/{site_name}/{map_type}")
async def get_risk_map(site_name: str, map_type: str):
    """Serve a previously generated colorized risk map PNG."""
    if map_type not in ["terrain", "active", "siltation"]:
        raise HTTPException(status_code=400, detail="Invalid map type")
        
    png_path = RISK_MAP_DIR / site_name / f"{map_type}_map_overlay.png"
    if not png_path.exists():
        raise HTTPException(status_code=404, detail=f"Risk map not found for '{site_name}'")
    return FileResponse(
        str(png_path), 
        media_type="image/png", 
        headers={"Cache-Control": "no-cache, no-store, must-revalidate"}
    )


@app.get("/api/dams")
async def get_dams(district: str | None = None):
    """Fetch Kerala dam water levels from public feed."""
    dam_urls = [
        "https://raw.githubusercontent.com/EmersionCyriac/Kerala-Dam-Water-Levels/main/live.json",
        "https://raw.githubusercontent.com/EmersionCyriac/Kerala-Dam-Water-Levels/main/irrigation_live.json",
    ]

    all_dams = []
    async with httpx.AsyncClient() as client:
        for url in dam_urls:
            try:
                resp = await client.get(url, timeout=15.0)
                resp.raise_for_status()
                raw = resp.json()
                if isinstance(raw, list):
                    all_dams.extend(raw)
                elif isinstance(raw, dict) and "dams" in raw:
                    all_dams.extend(raw["dams"])
            except Exception as e:
                print(f"  WARNING: Failed to fetch {url}: {e}")

    # Normalize and compute alert levels
    # The API has nested structure: dam.data[0] contains latest readings
    result = []
    for dam in all_dams:
        name = dam.get("name") or dam.get("damName") or "Unknown"
        lat = dam.get("latitude") or dam.get("lat")
        lon = dam.get("longitude") or dam.get("lon")
        red = dam.get("redLevel") or dam.get("red_level")
        orange = dam.get("orangeLevel") or dam.get("orange_level")
        blue = dam.get("blueLevel") or dam.get("blue_level")
        full = dam.get("FRL") or dam.get("fullReservoirLevel") or dam.get("full_level")
        dist = dam.get("district") or dam.get("District") or ""

        # Extract latest readings from nested data array
        latest = {}
        data_arr = dam.get("data", [])
        if isinstance(data_arr, list) and len(data_arr) > 0:
            latest = data_arr[0]

        level = latest.get("waterLevel") or dam.get("waterLevel") or dam.get("current_level")
        storage_pct = latest.get("storagePercentage") or dam.get("storagePercentage")
        spillway = latest.get("spillwayRelease") or dam.get("spillwayRelease") or "0"
        inflow = latest.get("inflow") or dam.get("inflow") or "0"
        rainfall = latest.get("rainfall") or dam.get("rainfall") or "0"

        # Determine alert level
        alert = "NORMAL"
        if level is not None:
            try:
                level_f = float(level)
                if red and level_f >= float(red):
                    alert = "RED"
                elif orange and level_f >= float(orange):
                    alert = "ORANGE"
                elif blue and level_f >= float(blue):
                    alert = "BLUE"
            except (ValueError, TypeError):
                pass

        entry = {
            "name": name,
            "district": dist,
            "current_level": level,
            "full_level": full,
            "storage_pct": storage_pct,
            "alert": alert,
            "spillway_release": spillway,
            "inflow": inflow,
            "rainfall": rainfall,
            "lat": lat,
            "lon": lon,
        }

        # Filter by district if provided
        if district:
            if district.lower() in str(dist).lower():
                result.append(entry)
        else:
            result.append(entry)

    # Sort by alert severity
    severity_order = {"RED": 0, "ORANGE": 1, "BLUE": 2, "NORMAL": 3, "UNKNOWN": 4}
    result.sort(key=lambda d: severity_order.get(d["alert"], 4))

    return result


@app.get("/api/alerts")
async def get_alerts():
    """Fetch disaster/weather alerts for Kerala from IMD Nowcast."""
    url = "https://mausam.imd.gov.in/imd_latest/contents/dist_nowcast_rss.php"

    KERALA_DISTRICTS = [
        "ALAPPUZHA", "ERNAKULAM", "IDUKKI", "KANNUR", "KASARAGOD", 
        "KOLLAM", "KOTTAYAM", "KOZHIKODE", "MALAPPURAM", "PALAKKAD", 
        "PATHANAMTHITTA", "THIRUVANANTHAPURAM", "THRISSUR", "WAYANAD"
    ]

    import xml.etree.ElementTree as ET

    try:
        async with httpx.AsyncClient() as client:
            resp = await client.get(url, timeout=10.0, follow_redirects=True)
            if resp.status_code != 200:
                return []
            rss_text = resp.text

        root = ET.fromstring(rss_text)
        alerts = []
        channel = root.find("channel")
        
        # Track seen alerts to avoid duplicates (IMD sometimes repeats)
        seen = set()

        if channel is not None:
            for item in channel.findall("item"):
                title = item.findtext("title", "").strip().upper()
                
                if title in KERALA_DISTRICTS or "KERALA" in title:
                    description = item.findtext("description", "").strip()
                    pub_date = item.findtext("pubDate", "")
                    link = item.findtext("link", "")

                    if description in seen:
                        continue
                    seen.add(description)

                    # Determine severity
                    severity = "YELLOW"
                    desc_lower = description.lower()
                    if "extreme" in desc_lower or "very heavy" in desc_lower or "severe" in desc_lower:
                        severity = "RED"
                    elif "heavy" in desc_lower or "moderate thunderstorms" in desc_lower:
                        severity = "ORANGE"

                    alerts.append({
                        "title": f"Warning for {title.title()}",
                        "description": description,
                        "published": pub_date,
                        "severity": severity,
                        "link": link,
                    })

        return alerts
    except Exception as e:
        print(f"  WARNING: IMD Alert fetch failed: {e}")
        return []


@app.get("/api/weather/{lat}/{lon}")
async def get_weather(lat: float, lon: float, api_key: str | None = None):
    """Fetch current weather for a location."""
    return await fetch_weather(lat, lon, api_key)


# ─────────────────────────── Static Files ───────────────────────────

# Mount static files LAST (catch-all)
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


# ─────────────────────────── Startup ───────────────────────────

if __name__ == "__main__":
    import uvicorn

    print(f"\n{'#'*60}")
    print(f"  KERALA DISASTER INTELLIGENCE DASHBOARD")
    print(f"  Starting server at http://localhost:8000")
    print(f"{'#'*60}\n")

    uvicorn.run(app, host="0.0.0.0", port=8000, log_level="info")
