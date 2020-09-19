&inputfiles
 eqdsk_input_file='EQDSK.eqdsk'
 vtk_input_file='GEOM.vtk',
/
&miscparameters
/
&plotselections
      plot_geofldx = .true.,
      plot_geoqm = .true.,
      plot_geoqvolm = .false.,
      plot_geoqvolx = .false.,
      plot_geoqx = .true.,
      plot_gnu = .true.,
      plot_gnum = .true.,
      plot_gnusil = .true.,
      plot_gnusilm = .true.,
/
&beqparameters
      beq_bdryopt=BOPT,
      beq_cenopt=4,
      beq_fldspec=FLDS,
      beq_nzetap=NZETP,
      beq_psiopt=2,
      beq_thetaopt=2,
      beq_vacuum_field_file='VFLD'
      beq_arip=ARIP,
      beq_mrip=MRIP,
      beq_deltheta=0.,
      beq_psiref=0.,
      beq_rmove=0.,
/
