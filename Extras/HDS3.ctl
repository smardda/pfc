&inputfiles
 vtk_input_file='../S/SHAD_geoqm.vtk',
/
&hdsgenparameters
 geometrical_type=2,
 limit_geobj_in_bin=80,
 margin_type=2,
/
&btreeparameters
 tree_type=3, ! multi-octree
 tree_ttalg=2, ! use nxyz
 tree_nxyz=1,4,4,
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
