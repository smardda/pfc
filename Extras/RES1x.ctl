&inputfiles
 vtk_input_file='RES.vtk',
 eqdsk_input_file='EQDSK.eqdsk'
/
&miscparameters
/
&plotselections
      plot_geoqx = .true.,
      plot_geoqvolx = .false.,
      plot_geoqm = .true.,
      plot_geofldx = .true.,
      plot_gnu = .true.,
      plot_gnusil = .true.,
/
&beqparameters
      beq_cenopt=2,
      beq_bdryopt=8,
      beq_psiopt=2,
      beq_thetaopt=2,
      beq_deltheta=0.,
      beq_rmove=-11.,
/
