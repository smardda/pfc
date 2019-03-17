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
      beq_fldspec=2,
      beq_cenopt=2,
      beq_bdryopt=9,
      beq_psiref=0.,
      beq_xiopt=4,
      beq_nzetap=12,
      beq_vacuum_field_file='MLF.sp3'
/
