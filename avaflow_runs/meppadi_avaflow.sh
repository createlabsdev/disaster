#!/bin/bash
set -e
g.region raster=dem_meppadi_sq_filled -p

r.avaflow.40G -e -a \
  prefix=meppadi_dbf \
  phases=1 \
  elevation=dem_meppadi_sq_filled \
  hrelease=release_meppadi_sq_filled \
  density=1800 \
  friction=30,20,0.0 \
  cohesion=0.0 \
  deformation=1.0

r.out.gdal input=meppadi_dbf_results/meppadi_dbf_ascii/meppadi_dbf_hflow_max output=/home/emersion/digital_twin/output/meppadi/avaflow_depth_max.tif format=GTiff --overwrite
r.out.gdal input=meppadi_dbf_results/meppadi_dbf_ascii/meppadi_dbf_vflow_max output=/home/emersion/digital_twin/output/meppadi/avaflow_velocity_max.tif format=GTiff --overwrite
r.out.gdal input=meppadi_dbf_results/meppadi_dbf_ascii/meppadi_dbf_pflow_max output=/home/emersion/digital_twin/output/meppadi/avaflow_pressure_max.tif format=GTiff --overwrite
