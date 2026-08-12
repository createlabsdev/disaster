r.external input=output/chellanam/dem_ref.tif output=dem_chellanam
r.external input=output/chellanam/landslide_release.tif output=release_chellanam
r.avaflow.40G -e -k elevation=dem_chellanam hrelease=release_chellanam prefix=chellanam_ava time=1200 friction=0.2,200 cellsize=90
r.out.gdal input=chellanam_ava_H_max output=output/chellanam/avaflow_h_max.tif format=GTiff
r.out.gdal input=chellanam_ava_V_max output=output/chellanam/avaflow_v_max.tif format=GTiff
r.out.gdal input=chellanam_ava_P_max output=output/chellanam/avaflow_p_max.tif format=GTiff
r.avaflow.40G --help | grep friction
exit
