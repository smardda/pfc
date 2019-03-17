&analysisfiles
 vtk_data_file='../G/RES.vtk',
 vtk_input_file='RES_powx.vtk'
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
     analyse_scalar='Q'
     required_statistics='max   ','intg  '
     new_key='angle'
     rekey=.true.
     number_of_clusters=36
/
