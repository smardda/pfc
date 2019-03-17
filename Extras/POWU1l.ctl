&inputfiles
 vtk_input_file='../S/SHAD_geoqm.vtk', ! shadowing geometry
 vtkres_input_file='../G/RES_geoqm.vtk', ! results geometry
 geoq_input_file='../S/SHAD_geoq.out',
 hds_input_file='../H/SHAD_hds.hds',
/
&miscparameters
/
&plotselections
      plot_powx = .true.,
      plot_flinx = .false.,
      plot_flinm = .false.,
      plot_flinends = .false.,
/
&powcalparameters
      shadow_control=1,
      calculation_type='local'
      more_profiles = .true.
/
&edgprofparameters
      decay_length=DECL,
      power_loss=POWe+06,
/
&odesparameters
      initial_dt=0.0001,
      abs_error=1.e-4,
      rel_error=1.e-4,
      max_zeta=6.,
/
