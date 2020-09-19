&inputfiles
 vtk_input_file='GEOM_geoqm.vtk', ! shadowing geometry
 vtkres_input_file='RES_geoqm.vtk', ! results geometry
 geoq_input_file='GEOM_geoq.out',
 hds_input_file='GEOMH_hds.hds',
/
&miscparameters
/
&plotselections
      plot_flinends = .false.,
      plot_flinm = .false.,
      plot_flinx = .false.,
      plot_powx = .true.,
/
&powcalparameters
      calculation_type='CAL'
      more_profiles = .true.
      shadow_control=1,
      termination_planes = .true.
/
