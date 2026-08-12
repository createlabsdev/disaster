#!/usr/bin/env python3
"""
generate_digital_twin.py — One-Code-Run Generalized Pipeline

Takes any bounding box, downloads the terrain data, runs the physics simulations,
and applies a pre-trained machine learning model to generate a risk map (Transfer Learning).
"""

import argparse
import os
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(description="Generate a Digital Twin for any location")
    parser.add_argument("site", help="Name of the new site (e.g., munnar)")
    parser.add_argument("west", type=float, help="Min Longitude")
    parser.add_argument("south", type=float, help="Min Latitude")
    parser.add_argument("east", type=float, help="Max Longitude")
    parser.add_argument("north", type=float, help="Max Latitude")
    parser.add_argument("--model-site", default="meppadi", help="Pre-trained model to use (chellanam for flood, meppadi for landslide)")
    parser.add_argument("--rain-intensity", type=float, default=100.0, help="Live rain intensity in mm/h")
    parser.add_argument("--rain-duration", type=float, default=3.0, help="Rain duration in hours")
    args = parser.parse_args()

    print(f"============================================================")
    print(f"  GENERATING DIGITAL TWIN FOR: {args.site.upper()}")
    print(f"  Bounding Box: [{args.west}, {args.south}, {args.east}, {args.north}]")
    print(f"  Using pre-trained AI model from: {args.model_site}")
    print(f"  Rain intensity: {args.rain_intensity} mm/h for {args.rain_duration}h")
    print(f"============================================================\n")

    base_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. Initialize Site (Download DEM and Landcover via GEE)
    print(">>> STEP 1: Fetching Terrain Data...")
    cmd_init = [sys.executable, os.path.join(base_dir, 'init_site.py'), 
                args.site, str(args.west), str(args.south), str(args.east), str(args.north)]
    subprocess.run(cmd_init, check=True)

    # 2. Run Physics Engines (ANUGA, Infinite Slope, RUSLE)
    print("\n>>> STEP 2: Running Physics Engines...")
    cmd_physics = [sys.executable, os.path.join(base_dir, 'run_pipeline_advanced.py'), '--site', args.site,
                   '--rain-intensity', str(args.rain_intensity), '--rain-duration', str(args.rain_duration)]
    subprocess.run(cmd_physics, check=True)

    # 3. Build Feature Stack (Using build_training_data.py but ignoring labels)
    print("\n>>> STEP 3: Building AI Feature Stack...")
    cmd_features = [sys.executable, os.path.join(base_dir, 'build_training_data.py'), args.site, '--n-samples', '1']
    # We pass n-samples 1 because we only care about generating the feature_stack.tif, not the CSV
    subprocess.run(cmd_features, check=False) # May fail on label check, but feature_stack is built first

    # 4. Predict Risk Map (Transfer Learning)
    print("\n>>> STEP 4: Applying Artificial Intelligence...")
    
    # We need a small script to just apply the joblib model to the feature stack
    predict_script = os.path.join(base_dir, 'predict_susceptibility.py')
    if not os.path.exists(predict_script):
        create_predict_script(predict_script)

    cmd_predict = [sys.executable, predict_script, args.site, args.model_site]
    subprocess.run(cmd_predict, check=True)

    print(f"\n============================================================")
    print(f"  DIGITAL TWIN COMPLETE!")
    print(f"  Check output/{args.site}/susceptibility_xgb_transfer.tif")
    print(f"============================================================")

def create_predict_script(path):
    """Creates a lightweight script just for inference/transfer learning"""
    code = '''#!/usr/bin/env python3
import sys, os
import joblib
import warnings

# Suppress pandas performance warnings from one-hot encoding
warnings.filterwarnings('ignore')

base_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(base_dir)
from train_susceptibility import generate_susceptibility_maps

site = sys.argv[1]
model_site = sys.argv[2]

# Load pre-trained model and scaler
model_path = os.path.join(base_dir, 'output', model_site, 'model_xgb.joblib')
scaler_path = os.path.join(base_dir, 'output', model_site, 'scaler.joblib')

if not os.path.exists(model_path):
    print(f"ERROR: Pre-trained model not found at {model_path}")
    sys.exit(1)

clf = joblib.load(model_path)
scaler = joblib.load(scaler_path)

# Reconstruct the feature names list (27 features total)
raw_features = ['elevation', 'slope', 'aspect', 'plan_curvature', 'profile_curvature', 
                'twi', 'dist_to_stream', 'soil_depth', 'cohesion', 'phi', 'k_factor', 
                'c_factor', 'r_factor', 'manning_n', 'rusle_soil_loss', 'factor_of_safety', 
                'avaflow_depth', 'avaflow_velocity']
LANDCOVER_CLASSES = [10, 20, 30, 40, 50, 60, 70, 80, 90]
feature_names = raw_features + [f"landcover_{c}" for c in LANDCOVER_CLASSES]

# Generate maps using the robust function from train_susceptibility
stack_path = os.path.join(base_dir, 'output', site, 'feature_stack.tif')
out_dir = os.path.join(base_dir, 'output', site)

generate_susceptibility_maps(
    model_rf=None,
    model_xgb=clf,
    scaler=scaler,
    feature_names=feature_names,
    feature_stack_path=stack_path,
    output_dir=out_dir,
    models_to_run='xgb'
)
'''
    with open(path, 'w') as f:
        f.write(code)

if __name__ == "__main__":
    main()
