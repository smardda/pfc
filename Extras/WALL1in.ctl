&inputfiles
 vtk_input_file='WALL.vtk',
 eqdsk_input_file='EQDSK.eqdsk'
/
&miscparameters
/
&plotselections
      plot_geoqx = .true.,
      plot_geoqvolx = .false.,
      plot_geoqvolm = .false.,
      plot_gnum = .false.
      plot_geoqm = .true.,
      plot_geofldx = .false.,
      plot_gnu = .true.,
/
&beqparameters
      beq_rmove=-6.,
      beq_cenopt=4,
      beq_psiopt=1,
      beq_bdryopt=3,
      beq_thetaopt=1,
/
