&inputfiles
 vtk_input_file='RES.vtk',
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
      plot_gnu = .false.,
      plot_gnusil = .true.,
      plot_gnusilm = .true.,
/
&beqparameters
      beq_cenopt=4,
      beq_bdryopt=1,
      beq_psiref=PSIREF,
      beq_psiopt=2,
      beq_thetaopt=2, 
      beq_deltheta=0.,
/
