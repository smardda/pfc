module scontrol_h

  use const_kind_m

!! public types

!> run control parameters
  type, public :: snumerics_t
     character(len=80) :: control !< option control parameter
     character(len=80) :: mode !< control mode of operation
     character(len=80) :: namekey !< sort key identifier of scalar field
     character(len=80) :: namescal !< identifier of scalar field to analyse
     logical :: lurcen !< user sets nominal central \f$ R \f$ for poloidal angle analysis
     logical :: luzcen !< user sets nominal central \f$ Z \f$ for poloidal angle analysis
     real(kr8) :: urcen !< user defined value of central \f$ R \f$
     real(kr8) :: uzcen !< user defined value of central \f$ Z \f$
     integer(ki4) :: totstat !< number of active entries in namstat
     character(len=80) :: newkey !< name of new key (only angle allowed)
     logical :: rekey !< logical control parameter
     integer(ki4) :: nbin !< number of bins in position angle (sets cluster size)
     real(kr8) :: realpar !< real control parameter
     integer(ki4) :: intpar !< integer control parameter
     logical :: logicpar !< logical control parameter
     character(len=80), dimension(:), allocatable :: namstat !< statistics control parameter
  end type snumerics_t

!> file names
  type, public :: sfiles_t
     integer(ki4) :: nvtkdata !< number of vtkdta files
     character(len=80)  :: smanalout       !< output data
     character(len=80)  :: out       !< output data
     character(len=80)  :: log           !< log file
     character(len=80), dimension(:), allocatable  :: vtkdata !< smanal input data files
     character(len=80)  :: vtk   !< special input vtk file
     character(len=80)  :: vtkfull         !< smanal output data file
     character(len=80)  :: vtksmall         !< smanal output data file
     character(len=80)  :: gnu !< gnuplot file
     character(len=80)  :: gnure !< gnuplot file for rekey
  end type sfiles_t

!> plot output selectors
  type, public :: splots_t
     logical  :: smanalout !< smanal output data selector
     logical  :: vtk !< vtk plot selector
     logical  :: vtksmall   !< vtk plot selector
     logical  :: gnu !< gnuplot plot selector
  end type splots_t

end module scontrol_h
