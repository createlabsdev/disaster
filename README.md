# 🌍 Kerala Digital Twin & AI Disaster Predictor
**Complete Project Overview**

This document outlines the entire journey of building a state-of-the-art Regional Risk Predictor capable of forecasting landslides and floods in Kerala, combining physics-based simulations with Artificial Intelligence.

---

## Phase 1: Data Acquisition & The "Digital Twin"
To predict disasters, we first needed to recreate Kerala's physical environment inside a computer—a "Digital Twin."

1. **Topography (Elevation):** We used the Google Earth Engine (GEE) API to download the **SRTM 90m Digital Elevation Model (DEM)**. This gave us a 3D map of the hills, valleys, and slopes for any given coordinate in Kerala.
2. **Land Cover:** We downloaded the **ESA WorldCover 10m** dataset to understand what was on the ground (forests, bare soil, urban concrete, water bodies).
3. **Soil Mechanics:** We estimated critical soil parameters (soil depth, cohesion, friction angle) using geological averages for Kerala.

## Phase 2: Physics Engines (Simulating the Elements)
Before applying AI, we ran massive mathematical physics engines on the 3D terrain to see how water and soil would behave under heavy rain. In total, the pipeline computes dozens of spatial parameters and physical states for every single 10-meter block of land:

1. **Topographic Derivatives (The Shape of the Land):** We run algorithms over the DEM to calculate:
    * **Slope Steepness:** Exactly how steep a hillside is (in degrees).
    * **Aspect:** Which direction the slope faces (affects soil moisture from the sun).
    * **Plan Curvature & Profile Curvature:** Determines if the terrain is convex (sheds water) or concave (collects water).
    * **Topographic Wetness Index (TWI):** A mathematical calculation of how much water will naturally pool in a specific pixel based on upstream drainage.
    * **Distance to Stream:** How close a pixel is to a natural riverbed.
2. **Soil Mechanics & Friction:**
    * **Manning's Roughness Coefficient (n):** Based on the land cover (e.g., dense forest vs. smooth concrete), we calculate how much friction the ground applies to flowing water.
    * **Soil Depth, Cohesion, & Friction Angle:** Critical geotechnical parameters that define how well the dirt holds itself together.
3. **Soil Erosion (RUSLE - Revised Universal Soil Loss Equation):** By combining rainfall intensity (R-factor) with land cover (C-factor) and soil type (K-factor), we calculated exactly how many tons of soil would wash away.
4. **Slope Stability (Factor of Safety):** We calculate the exact mathematical threshold at which a steep, rain-soaked hillside will collapse under its own weight (Factor of Safety < 1.0 means imminent landslide).
5. **Hydrodynamics (ANUGA):** We integrated the **ANUGA Hydrodynamic Engine**. We pour virtual rain over the 3D terrain and simulate the *Shallow Water Equations* across a computational mesh to calculate:
    * **Maximum Flood Depth (meters)**
    * **Maximum Flood Velocity (meters/second)**

## Phase 3: Artificial Intelligence (Machine Learning)
Physics alone isn't enough because real-world disasters are complex. We needed an AI to learn the patterns of *when* those physical parameters result in a disaster.

1. **Building the Dataset:** We ran the physics pipeline on **40 locations** across Kerala. 
    * **10 Landslide Epicenters** (e.g., Chooralmala, Kavalappara).
    * **10 Flood Zones** (e.g., Chellanam, Kuttanad).
    * **20 Safe Zones** (to teach the AI what *not* to flag).
2. **Feature Stacking (The 27 Parameters):** We stacked physical parameters into a massive dataset containing **80,784 pixels** of data. The 27 model input features include:
    * **18 Continuous Physical Calculations:** Elevation, Slope, Aspect, Plan Curvature, Profile Curvature, TWI (Topographic Wetness Index), Distance to Stream, Soil Depth, Cohesion, Friction Angle (Phi), K-Factor, C-Factor, R-Factor, Manning's n, RUSLE Soil Loss, Factor of Safety, Flood Depth (ANUGA), Flood Velocity (ANUGA).
    * **9 Categorical Land Cover Parameters:** The ESA WorldCover raster is one-hot encoded into 9 distinct binary categories (Tree cover, Shrubland, Grassland, Cropland, Built-up, Bare / sparse vegetation, Snow and ice, Permanent water bodies, Herbaceous wetland) representing exactly what is on the ground.
3. **Training XGBoost & Confidence Levels:** We trained an **XGBoost Classifier** to study these 27 features. Because XGBoost outputs probabilistic values (from 0.0 to 1.0) rather than simple binary classifications, the AI is able to assign a precise **Confidence Level** to every pixel on the map, allowing us to quantify exactly how mathematically certain the AI is of a disaster occurring.
4. **The Result:** The model achieved an outstanding **94.7% Accuracy (AUC)**, proving it successfully learned to distinguish safe areas from high-risk zones.

## Phase 4: Time-Series Hydrological Predictions (24, 48, 72 Hours)
Real-world disasters like landslides depend heavily on cumulative rainfall. We built a **Hydrological Memory Engine** to predict risk for future timeframes (+24h, +48h, +72h).
* When predicting 72 hours into the future, the system calculates the **48-hour Antecedent Rainfall** leading up to that exact future hour.
* This past rainfall acts as a saturation multiplier, creating an **Effective Rain Intensity** that artificially boosts the water volume sent into the physics engine.
* This mathematically forces the ANUGA simulation to model saturated soil and pre-flooded terrain, generating highly accurate time-series forecasts rather than relying solely on isolated 1-hour cloudbursts.

## Phase 5: The Interactive Dashboard
With a trained AI model, we built a modern web application to make the tool accessible.

1. **The Backend (FastAPI):** We built a Python server that acts as the brain. It waits for a user to search for a location, orchestrates the entire pipeline (downloading terrain, running physics, applying AI), and generates a colorized GeoTIFF risk map.
2. **Live APIs:** We integrated external live data feeds:
    * **Open-Meteo:** To fetch real-time rainfall, temperature, and humidity for the searched location. This live rain drives the physics engines.
    * **Kerala Dams:** We linked to a live JSON feed to monitor all Kerala dam water levels and alert statuses.
    * **IMD Nowcast:** We integrated the official India Meteorological Department RSS feed to pull live weather warnings for Kerala districts.
3. **The Frontend (UI):** We designed a stunning, glassmorphic dark-mode interface using pure HTML, CSS, and vanilla JavaScript. 
    * **Interactive Map (Leaflet):** Allows users to view the predicted Risk Map overlaid on top of Dark, Satellite, or Street maps.
    * **Real-time Metrics:** Displays the live weather, dam levels, official IMD alerts, and a dynamic 0-100% Flood Risk gauge.

## How It Works Today (The User Flow)
1. You type a place (e.g., "Pala") into the dashboard.
2. The dashboard geocodes it into coordinates and fetches the live rainfall.
3. The server downloads the 3D terrain for Pala and runs the physics engines using that live rainfall.
4. The AI analyzes the physics results and draws a Risk Map.
5. The dashboard instantly displays the transparent Red/Orange/Yellow danger zones right over your map.
