#!/usr/bin/env python3
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
