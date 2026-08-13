#!/bin/bash
grass /home/emersion/digital_twin/grassdata/chellanam_grass --exec bash <<'EOF_GRASS'
r.external input="/home/emersion/digital_twin/output/chellanam/dem_ref.tif" output=dem_chellanam
r.external input="/home/emersion/digital_twin/output/chellanam/landslide_release.tif" output=release_chellanam
r.avaflow.40G -e -k elevation=dem_chellanam hrelease=release_chellanam prefix=chellanam_ava time=1200 friction=0.2,200 cellsize=90
r.out.gdal input=chellanam_ava_H_max output="/home/emersion/digital_twin/output/chellanam/avaflow_h_max.tif" format=GTiff
r.out.gdal input=chellanam_ava_V_max output="/home/emersion/digital_twin/output/chellanam/avaflow_v_max.tif" format=GTiff
r.out.gdal input=chellanam_ava_P_max output="/home/emersion/digital_twin/output/chellanam/avaflow_p_max.tif" format=GTiff
EOF_GRASS
