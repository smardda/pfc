&inputfiles
 vtk_input_file='SHAD.vtk',
 eqdsk_input_file='EQDSK.eqdsk'
/
&miscparameters
/
&plotselections
      plot_geoqx = .true.,
      plot_geoqvolx = .false.,
      plot_geoqvolm = .false.,
      plot_gnum = .true.,
      plot_geoqm = .true.,
      plot_geofldx = .false.,
      plot_gnu = .true.,
      plot_gnusil = .false.,
      plot_gnusilm = .false.,
/
&beqparameters
      beq_cenopt=4,
      beq_bdryopt=2,
      beq_psiopt=2,
      beq_thetaopt=2,
      beq_deltheta=0.,
/
