&inputfiles
 vtk_input_file='GEOM_geoqm.vtk',
/
&hdsgenparameters
 geometrical_type=2,
 limit_geobj_in_bin=NBIN,
 margin_type=2,
/
&btreeparameters
 tree_nxyz=1,2,2,
 tree_ttalg=2, ! use nxyz
 tree_type=3, ! multi-octree
/
&positionparameters
position_transform=1,
/
&plotselections
      plot_geobjq = .true.,
      plot_geoptq = .false.,
      plot_hdsm = .true.,
      plot_hdsq = .true.,
/
