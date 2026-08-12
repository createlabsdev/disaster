r.external input=output/meppadi/dem_ref.tif output=dem_meppadi
r.external input=output/meppadi/landslide_release.tif output=release_meppadi
r.avaflow.40G -e -k elevation=dem_meppadi hrelease=release_meppadi prefix=meppadi_ava time=1200 phases=1 friction=0.2,200 cellsize=90
r.out.gdal input=meppadi_ava_H_max output=output/meppadi/avaflow_h_max.tif format=GTiff
r.out.gdal input=meppadi_ava_V_max output=output/meppadi/avaflow_v_max.tif format=GTiff
r.out.gdal input=meppadi_ava_P_max output=output/meppadi/avaflow_p_max.tif format=GTiff
r.avaflow.40G -e -k elevation=dem_meppadi hrelease=release_meppadi prefix=meppadi_ava time=1200 phases=1 friction=0.2,200,0.2,200,0.2,200 cellsize=90
r.avaflow.40G -e -k elevation=dem_meppadi hrelease=release_meppadi prefix=meppadi_ava time=1200 phases=1 density=2000 friction=0.2,200 cohesion=0 viscosity=0 deformation=0 cellsize=90
EXIT
exirt
exit
exit
