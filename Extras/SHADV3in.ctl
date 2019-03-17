&inputfiles
 vtk_input_file='SHAD.vtk',
 eqdsk_input_file='EQDSK.eqdsk'
/
&miscparameters
/
&plotselections
      plot_geoqx = .true.,
      plot_geoqvolx = .false.,
      plot_geoqm = .true.,
      plot_geofldx = .false.,
      plot_gnu = .true.,
      plot_gnusil = .false.,
/
&beqparameters
      beq_fldspec=3,
      beq_cenopt=4,
      beq_bdryopt=7,
      beq_xiopt=2,
      beq_nzetap=1,
/
