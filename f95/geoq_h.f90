module geoq_h

  use const_kind_m
  use position_h
  use geobjlist_h
  use geobj_m
  use fmesh_h
  use spl2d_m
  use spl3d_m
  use dcontrol_h
  use skyl_h
  use beq_h


!> data structure describing geometrical objects and equilibrium field
  type, public :: geoq_t
     type(geobjlist_t) :: objl !< geobj coord data and useful bits
     type(beq_t) :: beq !< equilibrium field
     type(skyl_t) :: skyl !< skylight definition(s)
  end type geoq_t

end module geoq_h
