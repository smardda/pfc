&odesparameters
      abs_error=1.e-4,
      initial_dt=1.e-4,
      max_numsteps=100, 
      max_zeta=6., ! radians
      rel_error=1.e-4,
      termination_parameters(1)=1., ! avoid self-intersection, quantised units
      termination_parameters(2)=0.1, ! max step-size, fraction of domain units
/
