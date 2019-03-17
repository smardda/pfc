&inputfiles
 vtk_input_file='../V/GEOM.vtk',
/
&hdsgenparameters
 geometrical_type=2,
 limit_geobj_in_bin=MAXBIN,
 margin_type=2,
 max_tolerance=0.01
/
&btreeparameters
 tree_type=3, ! multi-octree
 tree_ttalg=2, ! use nxyz
 tree_nxyz=2,2,1,
/
&positionparameters
position_transform=1,
/
&plotselections
      plot_hdsm = .true.,
      plot_hdsq = .true.,
      plot_geobjq = .true.,
      plot_geoptq = .true.,
/
