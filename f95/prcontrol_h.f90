module prcontrol_h

  use const_kind_m
  use pcontrol_h

!! public types

!> run control parameters
  type, public, extends(pnumerics_t)  :: pxnumerics_t
     integer(ki4)  :: pclecalc = 1    !< default 1, if nonzero, perform particle power calculation
     integer(ki4)  :: radncalc = 0   !< default 0, if nonzero, perform radiation power calculation
     integer(ki4)  :: shadowr = 0 !< default 0, if nonzero, shadow radiation power
     integer(ki4)  :: mant = 0 !< maximum number of antennas to follow (<= hardwired parameter)
     integer(ki4)  :: nant = 0  !< if non-zero, number of antennas to follow
     integer(ki4), dimension(:), allocatable  :: antno !< array containing antenna/element numbers
     logical  :: plotrayg = .false.   !< rays plot selector (copy from prplots_t below)
  end type pxnumerics_t

!> file names
  type, public, extends(pfiles_t)  :: prfiles_t
     character(len=80)  :: antdata = ''         !< antenna data
     character(len=80)  :: radata = ''        !< radiation data
     character(len=80)  :: powcalrout  = ''     !< output data
     character(len=80)  :: rayg = ''       !< rays output
  end type prfiles_t

!> plot output selectors
  type, public, extends(pplots_t)  :: prplots_t
     logical  :: ant = .false. !< antenna plot selector
     logical  :: radn = .false.  !< radn plot selector
     logical  :: rayg = .false.   !< rays plot selector
  end type prplots_t

end module prcontrol_h
