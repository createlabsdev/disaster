#!/bin/bash
set -e
g.region raster=dem_chellanam_sq_filled -p

r.avaflow.40G -e -a \
  prefix=chellanam_dbf \
  phases=1 \
  elevation=dem_chellanam_sq_filled \
  hrelease=release_chellanam_sq_filled \
  density=1800 \
  friction=30,20,0.0 \
  cohesion=0.0 \
  deformation=1.0
