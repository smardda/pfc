module apb_h


  use const_kind_m

! public types

!> data structure for 3D applied magnetic field
  type, public :: apb_t
     real(kr8), dimension(:,:,:,:), allocatable :: field !< 3D vector field
     real(kr8), dimension(:), allocatable :: pos1 !< list of positions, coordinate 1
     real(kr8), dimension(:), allocatable :: pos2 !< list of positions, coordinate 2
     real(kr8), dimension(:), allocatable :: pos3 !< list of positions, coordinate 3
     integer(ki4) :: n1   !< array dimension in coordinate 1
     integer(ki4) :: n2   !< array dimension in coordinate 2
     integer(ki4) :: n3   !< array dimension in coordinate 3
     integer(ki4) :: ncpt   !< number of components
  end type apb_t

end module apb_h
