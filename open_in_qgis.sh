#!/bin/bash
# Opens QGIS with all final hazard pipeline outputs loaded for both sites

cd ~/digital_twin/output

qgis \
  chellanam/hillshade.tif \
  chellanam/rusle_soil_loss.tif \
  chellanam/factor_of_safety.tif \
  chellanam/landslide_release.tif \
  chellanam/avaflow_depth_max.tif \
  chellanam/avaflow_velocity_max.tif \
  chellanam/anuga_depth_max.tif \
  chellanam/anuga_velocity_max.tif \
  meppadi/hillshade.tif \
  meppadi/rusle_soil_loss.tif \
  meppadi/factor_of_safety.tif \
  meppadi/landslide_release.tif \
  meppadi/avaflow_depth_max.tif \
  meppadi/avaflow_velocity_max.tif \
  meppadi/anuga_depth_max.tif \
  meppadi/anuga_velocity_max.tif
