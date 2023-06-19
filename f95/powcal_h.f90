module powcal_h

  use const_kind_m
  use geobjlist_h
  use geobj_m
  use position_h
  use fmesh_h
  use beq_h
  use odes_h
  use powres_h
  use edgprof_h
  use edgprof_m
  use spl2d_m
  use spl3d_m
  use pcontrol_h
  use termplane_h
  use smitermpi_h


! public types

!> type describing launch calculation data
  type, public :: powlau_t
     integer(ki4)  :: ndummy  !< type variable
  end type powlau_t

!> type describing all power calculation data
  type, public :: powcal_t
     type(powres_t) :: powres !< powres
     type(edgprof_t) :: edgprof !< edgprof
     type(powlau_t) :: powlau !< powlau
     type(odes_t) :: odes !< odes
     type(pnumerics_t) :: n !< control  parameters
     type(decomp_t) :: decomp !< MPI parameters
  end type powcal_t

end module powcal_h
