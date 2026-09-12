module atcontrol_h

  use const_kind_m

!! public types

!> run control parameters
  type, public :: atparams_t
     character(len=80) :: control !< option control parameter
     real(kr8) :: realpar !< real control parameter
     integer(ki4) :: intpar !< integer control parameter
     logical :: logicpar !< logical control parameter
  end type atparams_t

!> file names
  type, public :: atfiles_t
     character(len=80)  :: out       !< output data
     character(len=80)  :: log           !< log file
     character(len=80)  :: antestdata        !< antest input data file
     character(len=80)  :: antestout         !< antest output data file
     character(len=80)  :: vtk
     character(len=80)  :: powx
     character(len=80)  :: gnu !< gnuplot file
     character(len=80)  :: eqdsk         !< equilibrium type eqtype file name
     character(len=80)  :: equil         !< generic equilibrium type file name
     character(len=80)  :: antdata         !< antenna data
     character(len=80)  :: radata         !< radiation data
  end type atfiles_t

!> plot output selectors
  type, public :: atplots_t
     logical  :: antestout !< antest output data selector
     logical  :: vtk   !< vtk plot selector
     logical  :: gnu !< gnuplot plot selector
     logical  :: powx   !< powx plot selector
  end type atplots_t

end module atcontrol_h
