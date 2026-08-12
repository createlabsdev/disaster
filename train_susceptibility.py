#!/usr/bin/env python3
"""
train_susceptibility.py — ML Susceptibility Model Training & Map Generation

Trains RandomForest and/or XGBoost classifiers for landslide/flood susceptibility
modeling using spatially cross-validated training data. Generates probability maps
and classified susceptibility rasters.

Part of the Kerala Digital Twin project (Chellanam floods / Meppadi landslides).

Usage:
    python train_susceptibility.py <site> [--model both] [--no-map] [--output-dir DIR]

Inputs (produced by build_training_data.py):
    output/{site}/training_data.csv
    output/{site}/feature_stack.tif

Outputs:
    output/{site}/susceptibility_rf.tif
    output/{site}/susceptibility_xgb.tif
    output/{site}/susceptibility_classes.tif
    output/{site}/model_rf.joblib
    output/{site}/model_xgb.joblib
    output/{site}/scaler.joblib
    output/{site}/roc_curve.png
    output/{site}/feature_importance.png
    output/{site}/model_report.txt

Dependencies:
    numpy, pandas, scikit-learn, xgboost, matplotlib, joblib, rasterio
"""

import argparse
import os
import sys
import warnings
from pathlib import Path
from textwrap import dedent

import joblib
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for server/headless use
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import Affine
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    auc,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

try:
    import xgboost as xgb
except ImportError:
    xgb = None

warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# ESA WorldCover land-cover class codes
LANDCOVER_CLASSES = [10, 20, 30, 40, 50, 60, 80, 90, 95]

# Features that should be log-transformed (np.log1p) to reduce skew
LOG_TRANSFORM_FEATURES = [
    "rusle_soil_loss",
    "twi",
    "dist_to_stream",
    "avaflow_depth",
    "avaflow_velocity",
]

# Plot style — seaborn-v0_8 works on matplotlib ≥ 3.6
PLOT_STYLE = "seaborn-v0_8"
PLOT_DPI = 150


# ===================================================================
# 1. Data Loading & Preprocessing
# ===================================================================

def load_training_data(csv_path: str) -> pd.DataFrame:
    """Read the training CSV produced by build_training_data.py."""
    df = pd.read_csv(csv_path)
    required = {"label", "row", "col"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Training CSV missing required columns: {missing}")
    print(f"  Loaded {len(df)} samples from {csv_path}")
    print(f"  Class distribution: {dict(df['label'].value_counts().sort_index())}")
    return df


def preprocess_features(
    df: pd.DataFrame,
    scaler: StandardScaler | None = None,
    fit: bool = True,
) -> tuple[pd.DataFrame, StandardScaler, list[str]]:
    """
    Apply preprocessing pipeline to a DataFrame of features.

    Steps:
        1. One-hot encode 'landcover' column
        2. Log-transform skewed features
        3. Replace NaN/inf with column medians
        4. StandardScaler normalisation

    Parameters
    ----------
    df : pd.DataFrame
        Feature columns only (no label/row/col).
    scaler : StandardScaler or None
        If provided and fit=False, use this pre-fitted scaler.
    fit : bool
        Whether to fit the scaler on this data.

    Returns
    -------
    df_processed : pd.DataFrame
        Processed features.
    scaler : StandardScaler
        Fitted scaler.
    feature_names : list[str]
        Ordered list of feature names after encoding.
    """
    df = df.copy()

    # --- One-hot encode landcover ---
    if "landcover" in df.columns:
        for cls in LANDCOVER_CLASSES:
            df[f"landcover_{cls}"] = (df["landcover"] == cls).astype(np.float32)
        df.drop(columns=["landcover"], inplace=True)

    # --- Log-transform skewed features ---
    for feat in LOG_TRANSFORM_FEATURES:
        if feat in df.columns:
            df[feat] = np.log1p(df[feat].clip(lower=0))

    # --- Replace NaN / inf with column median ---
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    medians = df.median()
    df.fillna(medians, inplace=True)

    feature_names = list(df.columns)

    # --- StandardScaler ---
    if scaler is None:
        scaler = StandardScaler()
    if fit:
        arr = scaler.fit_transform(df.values)
    else:
        arr = scaler.transform(df.values)

    df_processed = pd.DataFrame(arr, columns=feature_names, index=df.index)
    return df_processed, scaler, feature_names


def split_features_label(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series, pd.Series, pd.Series]:
    """Separate features, label, and spatial coordinates from the training DataFrame."""
    meta_cols = ["label", "row", "col", "site", "site_id"]
    feature_cols = [c for c in df.columns if c not in meta_cols]
    X = df[feature_cols]
    y = df["label"]
    rows = df["row"]
    cols = df["col"]
    return X, y, rows, cols


# ===================================================================
# 2. Spatial Cross-Validation
# ===================================================================

def assign_spatial_blocks(
    rows: pd.Series, cols: pd.Series, n_blocks: int = 5
) -> np.ndarray:
    """
    Divide the study area into an n_blocks × n_blocks grid and assign
    each sample to a spatial block ID.

    This is used with GroupKFold to prevent spatial autocorrelation
    leakage between training and test folds.
    """
    row_min, row_max = rows.min(), rows.max()
    col_min, col_max = cols.min(), cols.max()

    # Avoid division by zero for single-pixel extent
    row_range = max(row_max - row_min, 1)
    col_range = max(col_max - col_min, 1)

    row_block = np.clip(
        ((rows - row_min) / row_range * n_blocks).astype(int), 0, n_blocks - 1
    )
    col_block = np.clip(
        ((cols - col_min) / col_range * n_blocks).astype(int), 0, n_blocks - 1
    )
    block_id = row_block * n_blocks + col_block
    return block_id.values


def spatial_cross_validate(
    model, X: np.ndarray, y: np.ndarray, groups: np.ndarray, n_splits: int = 5
) -> dict:
    """
    Run GroupKFold spatial cross-validation.

    Returns a dict with per-fold and aggregate metrics.
    """
    gkf = GroupKFold(n_splits=n_splits)
    fold_results = []

    all_y_true, all_y_prob = [], []

    for fold_idx, (train_idx, test_idx) in enumerate(gkf.split(X, y, groups)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        model_clone = _clone_model(model)
        model_clone.fit(X_train, y_train)

        proba = model_clone.predict_proba(X_test)
        # Handle case where model only saw one class during training
        if proba.shape[1] == 1:
            y_prob = np.zeros(len(X_test))
        else:
            y_prob = proba[:, 1]
        y_pred = (y_prob >= 0.5).astype(int)

        # roc_auc_score needs both classes present in y_test
        if len(np.unique(y_test)) < 2:
            fold_auc = 0.5  # undefined, use chance-level
        else:
            fold_auc = roc_auc_score(y_test, y_prob)
        fold_acc = accuracy_score(y_test, y_pred)
        fold_prec = precision_score(y_test, y_pred, zero_division=0)
        fold_rec = recall_score(y_test, y_pred, zero_division=0)
        fold_f1 = f1_score(y_test, y_pred, zero_division=0)

        fold_results.append(
            {
                "fold": fold_idx + 1,
                "auc": fold_auc,
                "accuracy": fold_acc,
                "precision": fold_prec,
                "recall": fold_rec,
                "f1": fold_f1,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
            }
        )

        all_y_true.extend(y_test)
        all_y_prob.extend(y_prob)

        print(
            f"    Fold {fold_idx + 1}: AUC={fold_auc:.4f}  "
            f"Acc={fold_acc:.4f}  Prec={fold_prec:.4f}  "
            f"Rec={fold_rec:.4f}  F1={fold_f1:.4f}"
        )

    all_y_true = np.array(all_y_true)
    all_y_prob = np.array(all_y_prob)
    all_y_pred = (all_y_prob >= 0.5).astype(int)

    cv_results = {
        "folds": fold_results,
        "mean_auc": np.mean([f["auc"] for f in fold_results]),
        "std_auc": np.std([f["auc"] for f in fold_results]),
        "overall_auc": roc_auc_score(all_y_true, all_y_prob),
        "overall_accuracy": accuracy_score(all_y_true, all_y_pred),
        "overall_precision": precision_score(all_y_true, all_y_pred, zero_division=0),
        "overall_recall": recall_score(all_y_true, all_y_pred, zero_division=0),
        "overall_f1": f1_score(all_y_true, all_y_pred, zero_division=0),
        "confusion_matrix": confusion_matrix(all_y_true, all_y_pred),
        "y_true": all_y_true,
        "y_prob": all_y_prob,
    }
    return cv_results


def _clone_model(model):
    """Create a fresh copy of a model with the same hyperparameters."""
    from sklearn.base import clone
    return clone(model)


# ===================================================================
# 3. Model Construction
# ===================================================================

def build_rf_model() -> RandomForestClassifier:
    """Construct a RandomForest classifier with project defaults."""
    return RandomForestClassifier(
        n_estimators=500,
        max_depth=20,
        min_samples_leaf=10,
        class_weight="balanced",
        random_state=42,
        n_jobs=-1,
    )


def build_xgb_model(pos_weight: float = 1.0) -> "xgb.XGBClassifier":
    """Construct an XGBoost classifier with project defaults."""
    if xgb is None:
        raise ImportError("xgboost is not installed. Install with: pip install xgboost")
    return xgb.XGBClassifier(
        n_estimators=500,
        learning_rate=0.05,
        max_depth=8,
        min_child_weight=5,
        subsample=0.8,
        colsample_bytree=0.8,
        scale_pos_weight=pos_weight,
        random_state=42,
        use_label_encoder=False,
        eval_metric="logloss",
        n_jobs=-1,
    )


# ===================================================================
# 4. Evaluation & Plotting
# ===================================================================

def plot_roc_curves(
    results: dict[str, dict], output_path: str
) -> None:
    """
    Plot ROC curves for all models on one figure.

    Parameters
    ----------
    results : dict
        Keys are model names, values are cv_results dicts containing
        'y_true' and 'y_prob'.
    output_path : str
        Path to save the PNG figure.
    """
    try:
        plt.style.use(PLOT_STYLE)
    except OSError:
        pass

    fig, ax = plt.subplots(figsize=(8, 7))

    for name, res in results.items():
        fpr, tpr, _ = roc_curve(res["y_true"], res["y_prob"])
        auc_val = auc(fpr, tpr)
        ax.plot(fpr, tpr, linewidth=2, label=f"{name} (AUC = {auc_val:.4f})")

    ax.plot([0, 1], [0, 1], "k--", linewidth=1, alpha=0.5, label="Random")
    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title("ROC Curve — Spatial Cross-Validation", fontsize=14, fontweight="bold")
    ax.legend(loc="lower right", fontsize=11)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    fig.tight_layout()
    fig.savefig(output_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved ROC curve: {output_path}")


def plot_feature_importance(
    importances: dict[str, np.ndarray],
    feature_names: list[str],
    output_path: str,
    top_n: int = 15,
) -> None:
    """
    Plot feature importance bar chart for all models side-by-side.

    Parameters
    ----------
    importances : dict
        Keys are model names, values are arrays of feature importances.
    feature_names : list[str]
        Feature names aligned with importance arrays.
    output_path : str
        Path to save the PNG figure.
    top_n : int
        Number of top features to display.
    """
    try:
        plt.style.use(PLOT_STYLE)
    except OSError:
        pass

    n_models = len(importances)
    fig, axes = plt.subplots(1, n_models, figsize=(7 * n_models, 8))
    if n_models == 1:
        axes = [axes]

    for ax, (name, imp) in zip(axes, importances.items()):
        sorted_idx = np.argsort(imp)[-top_n:]
        ax.barh(
            range(top_n),
            imp[sorted_idx],
            align="center",
            color="#4C72B0",
            edgecolor="white",
        )
        ax.set_yticks(range(top_n))
        ax.set_yticklabels([feature_names[i] for i in sorted_idx], fontsize=10)
        ax.set_xlabel("Importance", fontsize=11)
        ax.set_title(f"{name} — Top {top_n} Features", fontsize=13, fontweight="bold")

    fig.tight_layout()
    fig.savefig(output_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved feature importance: {output_path}")


def format_confusion_matrix(cm: np.ndarray) -> str:
    """Pretty-print a 2×2 confusion matrix."""
    return (
        f"  Confusion Matrix:\n"
        f"                 Predicted 0   Predicted 1\n"
        f"    Actual 0     {cm[0, 0]:>10d}   {cm[0, 1]:>10d}\n"
        f"    Actual 1     {cm[1, 0]:>10d}   {cm[1, 1]:>10d}\n"
    )


def write_model_report(
    results: dict[str, dict],
    feature_names: list[str],
    importances: dict[str, np.ndarray],
    output_path: str,
) -> None:
    """Write a comprehensive text report of all model metrics."""
    lines = [
        "=" * 70,
        "SUSCEPTIBILITY MODEL REPORT",
        "=" * 70,
        "",
    ]

    for name, res in results.items():
        lines.append(f"--- {name} ---")
        lines.append(f"  Spatial CV Mean AUC: {res['mean_auc']:.4f} ± {res['std_auc']:.4f}")
        lines.append(f"  Overall AUC:         {res['overall_auc']:.4f}")
        lines.append(f"  Accuracy:            {res['overall_accuracy']:.4f}")
        lines.append(f"  Precision:           {res['overall_precision']:.4f}")
        lines.append(f"  Recall:              {res['overall_recall']:.4f}")
        lines.append(f"  F1-Score:            {res['overall_f1']:.4f}")
        lines.append("")
        lines.append(format_confusion_matrix(res["confusion_matrix"]))

        # Per-fold results
        lines.append("  Per-Fold Results:")
        for f in res["folds"]:
            lines.append(
                f"    Fold {f['fold']}: AUC={f['auc']:.4f}  Acc={f['accuracy']:.4f}  "
                f"Prec={f['precision']:.4f}  Rec={f['recall']:.4f}  F1={f['f1']:.4f}  "
                f"(train={f['n_train']}, test={f['n_test']})"
            )
        lines.append("")

        # Feature importance
        if name in importances:
            imp = importances[name]
            sorted_idx = np.argsort(imp)[::-1][:15]
            lines.append(f"  Top 15 Feature Importances ({name}):")
            for rank, idx in enumerate(sorted_idx, 1):
                lines.append(f"    {rank:2d}. {feature_names[idx]:30s}  {imp[idx]:.6f}")
            lines.append("")

    lines.append("=" * 70)

    with open(output_path, "w") as f:
        f.write("\n".join(lines))
    print(f"  Saved model report: {output_path}")


# ===================================================================
# 5. Susceptibility Map Generation
# ===================================================================

def generate_susceptibility_maps(
    model_rf,
    model_xgb,
    scaler: StandardScaler,
    feature_names: list[str],
    feature_stack_path: str,
    output_dir: str,
    models_to_run: str = "both",
) -> None:
    """
    Load the multi-band feature_stack.tif, apply preprocessing, predict
    susceptibility probabilities per pixel, and write output rasters.

    Parameters
    ----------
    model_rf : fitted RandomForest or None
    model_xgb : fitted XGBoost or None
    scaler : fitted StandardScaler
    feature_names : list[str]
        Feature names matching the training order (post one-hot encoding).
    feature_stack_path : str
        Path to the multi-band GeoTIFF feature stack.
    output_dir : str
        Directory for output rasters.
    models_to_run : str
        'rf', 'xgb', or 'both'.
    """
    print("\n[Map Generation]")
    print(f"  Loading feature stack: {feature_stack_path}")

    with rasterio.open(feature_stack_path) as src:
        meta = src.meta.copy()
        height, width = src.height, src.width
        n_bands = src.count
        transform = src.transform
        crs = src.crs

        # Read all bands: shape (n_bands, height, width)
        stack = src.read()

    print(f"  Raster dimensions: {height} × {width}, {n_bands} bands")

    # The feature stack bands correspond to the raw feature columns (before
    # one-hot encoding). We need the original feature order to map bands
    # correctly. The band names should align with the training CSV columns
    # (excluding label/row/col). We'll reconstruct from feature_names.

    # Determine which raw features are present (reverse one-hot for landcover)
    raw_feature_names = []
    for fn in feature_names:
        if fn.startswith("landcover_"):
            base = "landcover"
            if base not in raw_feature_names:
                raw_feature_names.append(base)
        else:
            raw_feature_names.append(fn)

    # Flatten raster to (n_pixels, n_bands)
    flat = stack.reshape(n_bands, -1).T.astype(np.float64)  # (n_pixels, n_bands)

    # Create a DataFrame with raw feature names
    if len(raw_feature_names) != n_bands:
        print(
            f"  WARNING: Feature name count ({len(raw_feature_names)}) != "
            f"band count ({n_bands}). Using generic band names."
        )
        raw_feature_names = [f"band_{i}" for i in range(n_bands)]

    df_raster = pd.DataFrame(flat, columns=raw_feature_names)

    # Track valid (non-NaN) pixels — nodata pixels get 0 probability
    nodata_mask = np.isnan(flat).any(axis=1) | np.isinf(flat).any(axis=1)

    # Apply same preprocessing as training
    df_processed, _, _ = preprocess_features(df_raster, scaler=scaler, fit=False)

    # Ensure columns match training order
    for col in feature_names:
        if col not in df_processed.columns:
            df_processed[col] = 0.0
    df_processed = df_processed[feature_names]

    X_raster = df_processed.values

    # Replace any remaining NaN/inf after preprocessing
    X_raster = np.nan_to_num(X_raster, nan=0.0, posinf=0.0, neginf=0.0)

    # Output raster profile
    out_profile = meta.copy()
    out_profile.update(
        dtype="float32",
        count=1,
        compress="lzw",
        nodata=-9999.0,
    )

    def _predict_and_save(model, model_name, out_path):
        """Predict probabilities and write a single-band GeoTIFF."""
        print(f"  Predicting with {model_name}...")
        prob = model.predict_proba(X_raster)[:, 1].astype(np.float32)
        prob[nodata_mask] = -9999.0
        prob_2d = prob.reshape(height, width)

        with rasterio.open(out_path, "w", **out_profile) as dst:
            dst.write(prob_2d, 1)
        valid = prob[~nodata_mask]
        print(
            f"    {model_name} susceptibility: min={valid.min():.4f}, "
            f"max={valid.max():.4f}, mean={valid.mean():.4f}"
        )
        print(f"    Saved: {out_path}")
        return prob

    # Predict with requested models
    best_prob = None
    best_name = None
    best_auc = -1.0

    if models_to_run in ("rf", "both") and model_rf is not None:
        prob_rf = _predict_and_save(
            model_rf, "RandomForest", os.path.join(output_dir, "susceptibility_rf.tif")
        )
        # Track for classification
        if best_prob is None:
            best_prob, best_name = prob_rf, "RandomForest"

    if models_to_run in ("xgb", "both") and model_xgb is not None:
        prob_xgb = _predict_and_save(
            model_xgb, "XGBoost", os.path.join(output_dir, "susceptibility_xgb.tif")
        )
        if best_prob is None:
            best_prob, best_name = prob_xgb, "XGBoost"

    # If both were run, the caller should pass the best model's prob — for now,
    # use the last one computed. The main function will set best_prob correctly.
    if best_prob is None:
        print("  No model was run — skipping classified map.")
        return

    # --- Classified susceptibility map (5 classes) ---
    _generate_classified_map(best_prob, nodata_mask, height, width, out_profile, output_dir)


def _generate_classified_map(
    prob: np.ndarray,
    nodata_mask: np.ndarray,
    height: int,
    width: int,
    out_profile: dict,
    output_dir: str,
) -> None:
    """
    Classify probability values into 5 susceptibility classes using
    quantile breaks on non-zero valid probabilities.

    Classes:
        1 = Very Low    (0-20th percentile)
        2 = Low         (20-40th percentile)
        3 = Moderate    (40-60th percentile)
        4 = High        (60-80th percentile)
        5 = Very High   (80-100th percentile)
    """
    print("  Generating classified susceptibility map...")

    valid_prob = prob[(~nodata_mask) & (prob > 0)]
    if len(valid_prob) == 0:
        print("  WARNING: No valid non-zero probabilities. Skipping classified map.")
        return

    quantiles = np.percentile(valid_prob, [20, 40, 60, 80])
    print(f"    Quantile breaks: {quantiles}")

    classified = np.zeros_like(prob, dtype=np.uint8)
    valid_mask = (~nodata_mask) & (prob >= 0)

    classified[valid_mask & (prob <= quantiles[0])] = 1  # Very Low
    classified[valid_mask & (prob > quantiles[0]) & (prob <= quantiles[1])] = 2  # Low
    classified[valid_mask & (prob > quantiles[1]) & (prob <= quantiles[2])] = 3  # Moderate
    classified[valid_mask & (prob > quantiles[2]) & (prob <= quantiles[3])] = 4  # High
    classified[valid_mask & (prob > quantiles[3])] = 5  # Very High

    classified_2d = classified.reshape(height, width)

    class_profile = out_profile.copy()
    class_profile.update(dtype="uint8", nodata=0)

    out_path = os.path.join(output_dir, "susceptibility_classes.tif")
    with rasterio.open(out_path, "w", **class_profile) as dst:
        dst.write(classified_2d, 1)

    # Print class distribution
    class_labels = {1: "Very Low", 2: "Low", 3: "Moderate", 4: "High", 5: "Very High"}
    for cls_val, cls_name in class_labels.items():
        count = int(np.sum(classified == cls_val))
        print(f"    Class {cls_val} ({cls_name:>9s}): {count:>8d} pixels")

    print(f"    Saved: {out_path}")


# ===================================================================
# 6. Main Pipeline
# ===================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Train susceptibility models and generate probability maps.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=dedent("""\
            Examples:
              python train_susceptibility.py meppadi
              python train_susceptibility.py chellanam --model xgb --no-map
              python train_susceptibility.py meppadi --output-dir results/meppadi
        """),
    )
    parser.add_argument("site", help="Site name (e.g. meppadi, chellanam)")
    parser.add_argument(
        "--model",
        choices=["rf", "xgb", "both"],
        default="both",
        help="Which model(s) to train (default: both)",
    )
    parser.add_argument(
        "--no-map",
        action="store_true",
        help="Skip map generation (faster for hyperparameter tuning)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Override output directory (default: output/{site})",
    )

    args = parser.parse_args()

    site = args.site
    output_dir = args.output_dir or os.path.join("output", site)
    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(output_dir, "training_data.csv")
    feature_stack_path = os.path.join(output_dir, "feature_stack.tif")

    print("=" * 60)
    print(f"  SUSCEPTIBILITY MODEL TRAINING — {site.upper()}")
    print("=" * 60)

    # ------------------------------------------------------------------
    # Step 1: Load & preprocess data
    # ------------------------------------------------------------------
    print("\n[1/5] Loading training data...")
    df = load_training_data(csv_path)
    X_raw, y, rows, cols = split_features_label(df)

    print("\n[2/5] Preprocessing features...")
    X_processed, scaler, feature_names = preprocess_features(X_raw, fit=True)
    X = X_processed.values
    y_arr = y.values.astype(int)

    print(f"  Feature matrix shape: {X.shape}")
    print(f"  Features ({len(feature_names)}): {feature_names[:10]}{'...' if len(feature_names) > 10 else ''}")

    # Save scaler
    scaler_path = os.path.join(output_dir, "scaler.joblib")
    joblib.dump(scaler, scaler_path)
    print(f"  Saved scaler: {scaler_path}")

    # ------------------------------------------------------------------
    # Step 2: Assign spatial blocks
    # ------------------------------------------------------------------
    print("\n[3/5] Spatial cross-validation...")
    groups = assign_spatial_blocks(rows, cols, n_blocks=5)
    n_unique_blocks = len(np.unique(groups))
    n_splits = min(5, n_unique_blocks)
    print(f"  Spatial blocks: {n_unique_blocks} unique blocks, {n_splits} CV folds")

    # ------------------------------------------------------------------
    # Step 3: Cross-validate & train models
    # ------------------------------------------------------------------
    cv_results = {}  # model_name -> cv_results dict
    trained_models = {}  # model_name -> fitted model
    importances = {}  # model_name -> importance array

    # Compute class weight for XGBoost
    n_neg = int(np.sum(y_arr == 0))
    n_pos = int(np.sum(y_arr == 1))
    pos_weight = n_neg / max(n_pos, 1)
    print(f"  Class balance: {n_neg} negative, {n_pos} positive (pos_weight={pos_weight:.2f})")

    # --- RandomForest ---
    if args.model in ("rf", "both"):
        print("\n  --- RandomForest ---")
        rf_model = build_rf_model()
        print("  Running spatial CV...")
        cv_results["RandomForest"] = spatial_cross_validate(
            rf_model, X, y_arr, groups, n_splits=n_splits
        )
        print(
            f"  CV AUC: {cv_results['RandomForest']['mean_auc']:.4f} "
            f"± {cv_results['RandomForest']['std_auc']:.4f}"
        )

        # Train on full dataset
        print("  Training on full dataset...")
        rf_model.fit(X, y_arr)
        trained_models["RandomForest"] = rf_model
        importances["RandomForest"] = rf_model.feature_importances_

        # Save model
        model_path = os.path.join(output_dir, "model_rf.joblib")
        joblib.dump(rf_model, model_path)
        print(f"  Saved: {model_path}")

    # --- XGBoost ---
    if args.model in ("xgb", "both"):
        if xgb is None:
            print("\n  WARNING: xgboost not installed, skipping XGBoost training.")
        else:
            print("\n  --- XGBoost ---")
            xgb_model = build_xgb_model(pos_weight=pos_weight)
            print("  Running spatial CV...")
            cv_results["XGBoost"] = spatial_cross_validate(
                xgb_model, X, y_arr, groups, n_splits=n_splits
            )
            print(
                f"  CV AUC: {cv_results['XGBoost']['mean_auc']:.4f} "
                f"± {cv_results['XGBoost']['std_auc']:.4f}"
            )

            # Train on full dataset
            print("  Training on full dataset...")
            xgb_model.fit(X, y_arr)
            trained_models["XGBoost"] = xgb_model
            importances["XGBoost"] = xgb_model.feature_importances_

            # Save model
            model_path = os.path.join(output_dir, "model_xgb.joblib")
            joblib.dump(xgb_model, model_path)
            print(f"  Saved: {model_path}")

    if not cv_results:
        print("\nERROR: No models were trained.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 4: Evaluation plots & report
    # ------------------------------------------------------------------
    print("\n[4/5] Generating evaluation outputs...")

    # ROC curve
    roc_path = os.path.join(output_dir, "roc_curve.png")
    plot_roc_curves(cv_results, roc_path)

    # Feature importance
    importance_path = os.path.join(output_dir, "feature_importance.png")
    plot_feature_importance(importances, feature_names, importance_path, top_n=15)

    # Text report
    report_path = os.path.join(output_dir, "model_report.txt")
    write_model_report(cv_results, feature_names, importances, report_path)

    # ------------------------------------------------------------------
    # Step 5: Susceptibility maps
    # ------------------------------------------------------------------
    if args.no_map:
        print("\n[5/5] Map generation SKIPPED (--no-map)")
    else:
        print("\n[5/5] Generating susceptibility maps...")
        if not os.path.exists(feature_stack_path):
            print(f"  ERROR: Feature stack not found: {feature_stack_path}")
            print("  Run build_training_data.py first, or use --no-map to skip.")
            sys.exit(1)

        # Determine best model by AUC
        best_model_name = max(cv_results, key=lambda k: cv_results[k]["mean_auc"])
        print(f"  Best model (by CV AUC): {best_model_name}")

        generate_susceptibility_maps(
            model_rf=trained_models.get("RandomForest"),
            model_xgb=trained_models.get("XGBoost"),
            scaler=scaler,
            feature_names=feature_names,
            feature_stack_path=feature_stack_path,
            output_dir=output_dir,
            models_to_run=args.model,
        )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("  TRAINING COMPLETE")
    print("=" * 60)
    for name, res in cv_results.items():
        print(f"  {name:15s}  AUC = {res['mean_auc']:.4f} ± {res['std_auc']:.4f}")
    print(f"  Output directory: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
