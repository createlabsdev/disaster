#!/usr/bin/env python3
"""
run_ml_pipeline.py — Full ML Pipeline Orchestrator

Orchestrates the complete susceptibility modeling pipeline for the Kerala
Digital Twin project, running each step as a subprocess with timing and
error handling.

Pipeline steps:
    1. compute_soil_params.py   — derive soil parameters
    2. GEE label export         — download flood/landslide labels (optional)
    3. build_training_data.py   — assemble training CSV & feature stack
    4. train_susceptibility.py  — train models & generate maps
    5. Summary report           — list outputs with raster statistics

Usage:
    python run_ml_pipeline.py <site> [--skip-gee] [--label-source FILE]
                                      [--n-samples N] [--model both] [--seed 42]

Sites:
    chellanam — Coastal flood susceptibility (2018 Kerala floods)
    meppadi   — Landslide susceptibility (2024 Wayanad landslides)

Dependencies:
    numpy, rasterio, subprocess, argparse (plus all pipeline script deps)
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
from textwrap import dedent

# ---------------------------------------------------------------------------
# Configuration: maps site name → GEE script
# ---------------------------------------------------------------------------
GEE_SCRIPTS = {
    "chellanam": "gee_flood_extent.py",
    "meppadi": "gee_landslide_scars.py",
}

VALID_SITES = list(GEE_SCRIPTS.keys())


# ===================================================================
# Utility helpers
# ===================================================================

def _timestamp() -> str:
    """Return a human-readable timestamp for log output."""
    return time.strftime("%H:%M:%S")


def _format_elapsed(seconds: float) -> str:
    """Format elapsed time as a human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes = int(seconds // 60)
    secs = seconds % 60
    if minutes < 60:
        return f"{minutes}m {secs:.1f}s"
    hours = int(minutes // 60)
    mins = minutes % 60
    return f"{hours}h {mins}m {secs:.0f}s"


def run_step(
    step_num: int,
    step_name: str,
    cmd: list[str],
    cwd: str | None = None,
) -> float:
    """
    Run a pipeline step as a subprocess.

    Parameters
    ----------
    step_num : int
        Step number for display.
    step_name : str
        Human-readable step description.
    cmd : list[str]
        Command and arguments.
    cwd : str or None
        Working directory for the subprocess.

    Returns
    -------
    elapsed : float
        Wall-clock time in seconds.

    Raises
    ------
    SystemExit
        If the subprocess exits with a non-zero return code.
    """
    print(f"\n{'=' * 60}")
    print(f"  STEP {step_num}: {step_name}")
    print(f"  Command: {' '.join(cmd)}")
    print(f"  Started: {_timestamp()}")
    print("=" * 60)

    t0 = time.time()
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    except FileNotFoundError as e:
        elapsed = time.time() - t0
        print(f"\n  ERROR: Command not found: {e}")
        print(f"  Step {step_num} FAILED after {_format_elapsed(elapsed)}")
        sys.exit(1)

    elapsed = time.time() - t0

    # Stream the captured output
    if result.stdout:
        for line in result.stdout.rstrip("\n").split("\n"):
            print(f"  | {line}")

    if result.returncode != 0:
        print(f"\n  ERROR: Step {step_num} exited with code {result.returncode}")
        print(f"  Step {step_num} FAILED after {_format_elapsed(elapsed)}")
        sys.exit(result.returncode)

    print(f"\n  Step {step_num} completed in {_format_elapsed(elapsed)}")
    return elapsed


def print_raster_summary(path: str) -> None:
    """Read a GeoTIFF and print basic statistics."""
    try:
        import numpy as np
        import rasterio
    except ImportError:
        print(f"    {path} (rasterio not available for stats)")
        return

    if not os.path.exists(path):
        print(f"    {path} — NOT FOUND")
        return

    try:
        with rasterio.open(path) as src:
            data = src.read(1)
            nodata = src.nodata
            if nodata is not None:
                valid = data[data != nodata]
            else:
                valid = data[~np.isnan(data)]

            if len(valid) == 0:
                print(f"    {path} — {src.width}×{src.height}, NO VALID PIXELS")
            else:
                print(
                    f"    {path}\n"
                    f"      Size: {src.width}×{src.height}, CRS: {src.crs}\n"
                    f"      Min: {valid.min():.4f}, Max: {valid.max():.4f}, "
                    f"Mean: {valid.mean():.4f}, Valid pixels: {len(valid)}"
                )
    except Exception as e:
        print(f"    {path} — ERROR reading: {e}")


def print_file_summary(path: str) -> None:
    """Print file existence and size."""
    if os.path.exists(path):
        size_kb = os.path.getsize(path) / 1024
        if size_kb > 1024:
            size_str = f"{size_kb / 1024:.1f} MB"
        else:
            size_str = f"{size_kb:.1f} KB"
        print(f"    {path} ({size_str})")
    else:
        print(f"    {path} — NOT FOUND")


# ===================================================================
# Main pipeline
# ===================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Run the full ML susceptibility modeling pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=dedent("""\
            Examples:
              python run_ml_pipeline.py meppadi
              python run_ml_pipeline.py chellanam --skip-gee --model xgb
              python run_ml_pipeline.py meppadi --label-source labels.tif --n-samples 10000
        """),
    )
    parser.add_argument(
        "site",
        choices=VALID_SITES,
        help=f"Site name ({', '.join(VALID_SITES)})",
    )
    parser.add_argument(
        "--skip-gee",
        action="store_true",
        help="Skip GEE label export (use existing label rasters)",
    )
    parser.add_argument(
        "--label-source",
        default=None,
        metavar="FILE",
        help="Path to a pre-downloaded label raster (overrides GEE export)",
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=5000,
        help="Number of samples per class (default: 5000)",
    )
    parser.add_argument(
        "--model",
        choices=["rf", "xgb", "both"],
        default="both",
        help="Which model(s) to train (default: both)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )

    args = parser.parse_args()

    site = args.site
    output_dir = os.path.join("output", site)
    os.makedirs(output_dir, exist_ok=True)

    # Resolve script directory (assumes all scripts live alongside this one)
    script_dir = os.path.dirname(os.path.abspath(__file__))

    print("=" * 60)
    print(f"  DIGITAL TWIN ML PIPELINE — {site.upper()}")
    print(f"  Started: {_timestamp()}")
    print(f"  Output:  {os.path.abspath(output_dir)}")
    print("=" * 60)

    pipeline_start = time.time()
    step_times = {}

    # ==================================================================
    # Step 1: Compute soil parameters
    # ==================================================================
    soil_script = os.path.join(script_dir, "compute_soil_params.py")
    cmd = [sys.executable, soil_script, site]
    step_times["soil_params"] = run_step(
        1, "Compute Soil Parameters", cmd, cwd=script_dir
    )

    # ==================================================================
    # Step 2: GEE label export (optional)
    # ==================================================================
    if args.skip_gee or args.label_source:
        print(f"\n  STEP 2: GEE label export — SKIPPED")
        if args.label_source:
            print(f"    Using pre-downloaded labels: {args.label_source}")
        step_times["gee_export"] = 0.0
    else:
        gee_script_name = GEE_SCRIPTS.get(site)
        if gee_script_name is None:
            print(f"\n  WARNING: No GEE script configured for site '{site}'. Skipping.")
            step_times["gee_export"] = 0.0
        else:
            gee_script = os.path.join(script_dir, gee_script_name)
            cmd = [sys.executable, gee_script, site]
            step_times["gee_export"] = run_step(
                2, f"GEE Label Export ({gee_script_name})", cmd, cwd=script_dir
            )

    # ==================================================================
    # Step 3: Build training data
    # ==================================================================
    build_script = os.path.join(script_dir, "build_training_data.py")
    cmd = [
        sys.executable,
        build_script,
        site,
        "--n-samples",
        str(args.n_samples),
        "--seed",
        str(args.seed),
    ]
    if args.label_source:
        cmd.extend(["--label-source", args.label_source])

    step_times["build_data"] = run_step(
        3, "Build Training Data", cmd, cwd=script_dir
    )

    # ==================================================================
    # Step 4: Train susceptibility models
    # ==================================================================
    train_script = os.path.join(script_dir, "train_susceptibility.py")
    cmd = [
        sys.executable,
        train_script,
        site,
        "--model",
        args.model,
    ]

    step_times["train_models"] = run_step(
        4, "Train Susceptibility Models", cmd, cwd=script_dir
    )

    # ==================================================================
    # Step 5: Summary report
    # ==================================================================
    total_elapsed = time.time() - pipeline_start

    print("\n" + "=" * 60)
    print("  PIPELINE SUMMARY")
    print("=" * 60)

    # Timing breakdown
    print("\n  Step Timing:")
    for step_name, elapsed in step_times.items():
        print(f"    {step_name:20s}  {_format_elapsed(elapsed)}")
    print(f"    {'TOTAL':20s}  {_format_elapsed(total_elapsed)}")

    # Output file inventory
    print("\n  Output Files:")
    print("  " + "-" * 50)

    # GeoTIFF outputs — print with raster stats
    tif_outputs = [
        "susceptibility_rf.tif",
        "susceptibility_xgb.tif",
        "susceptibility_classes.tif",
        "feature_stack.tif",
    ]
    print("\n  Raster Outputs:")
    for fname in tif_outputs:
        fpath = os.path.join(output_dir, fname)
        if os.path.exists(fpath):
            print_raster_summary(fpath)

    # Model & data files
    other_outputs = [
        "training_data.csv",
        "model_rf.joblib",
        "model_xgb.joblib",
        "scaler.joblib",
        "model_report.txt",
        "roc_curve.png",
        "feature_importance.png",
    ]
    print("\n  Model & Data Files:")
    for fname in other_outputs:
        fpath = os.path.join(output_dir, fname)
        print_file_summary(fpath)

    print("\n" + "=" * 60)
    print(f"  Pipeline completed in {_format_elapsed(total_elapsed)}")
    print(f"  All outputs saved to: {os.path.abspath(output_dir)}")
    print("=" * 60)


if __name__ == "__main__":
    main()
