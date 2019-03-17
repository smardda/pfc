&inputfiles
 vtk_input_file='GEOM.vtk',
 eqdsk_input_file='EQDSK.eqdsk'
/
&miscparameters
/
&plotselections
      plot_geoqx = .true.,
      plot_geoqvolx = .false.,
      plot_geoqm = .true.,
      plot_geofldx = .true.,
      plot_gnu = .false.,
      plot_gnusil = .false.,
/
&beqparameters
      beq_fldspec=FLDSPEC,
      beq_cenopt=4,
      beq_bdryopt=10,
      beq_thetaopt=2,
      beq_psiopt=2,
      beq_deltheta=0.01,
/
