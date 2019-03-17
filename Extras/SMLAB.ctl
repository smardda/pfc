&analysisfiles
 vtk_data_file='RES.vtk',
 vtk_input_file='RES.vtk'
/
&miscparameters
/
&plotselections
      plot_smanalout = .true.,
      plot_anx = .true.,
      plot_ansmallx = .true.,
      plot_angnu = .true.,
/
&smanalparameters
     analysis_mode='regular'
     sort_key='Body'
     analyse_scalar='Body'
     required_statistics='max   '
     rekey=.false.
/
