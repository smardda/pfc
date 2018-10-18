module termplane_h

  use const_kind_m

! public types
!> coordinate directions \f$ (R,Z,\xi) \f$ correspond to resp. \f$ 1,2,3 \f$
!! Signed integer \f$ \pm 1-3 \f$ gives direction of travel in corresponding coordinate
  type, public :: termplane_t
   !> vector describing termination plane
   !! 1. position in coordinate given by termplanedir
   !! 2. Auxilliary test value for code 4
     real(kr4), dimension(:,:), allocatable :: termplane !< -
   !> Storage of vector for extremum test. 1.e+30 flags first time per track,
   !! thereafter becomes position of last but one step.
     real(kr4), dimension(:,:), allocatable :: termstore !< -
   !> termplane direction and operation, entries
   !! 1. coordinate direction of travel to be tested
   !! 2. sign of direction of travel
   !! 3. code determining whether
   !! - 1. going past a plane specified by termplane in a certain direction given by 1 and 2
   !! - 2. intersecting a plane specified by termplane in a certain coordinate given by 1
   !! - 3. reaching an extremum in a certain coordinate given by 1, max/min from 2
   !! - 4. intersecting a plane with auxilliary condition
   !! 4. direction in which auxilliary condition to be applied
   !! 5. sign of direction of auxilliary condition
     integer(ki4), dimension(:,:), allocatable :: termplanedir !< -
     integer(ki4)  :: ntermplane !< number of termination planes
     integer(ki4)  :: ntermactive !< number of active termination planes
  end type termplane_t

end module termplane_h
