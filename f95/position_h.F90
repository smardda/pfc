module position_h

  use const_kind_m

! public types

  logical, parameter :: POSITION_OVERRIDE_IFMIF=.FALSE. !< use values required by IFMIF project
  type, public :: posvecl_t
     real(kr4), dimension(3) :: posvec !< position vector
  end type posvecl_t

  type, public :: posveclis_t
     type(posvecl_t), dimension(:), allocatable :: pos !< type variable
     integer(ki4)  :: np !< number of position vectors
     integer(ki2par)  :: numnparpos=4 !< number of entries in nparpos
   !> Integer parameters describing transform of vectors rel to Cartesians in m
   !! Array index entries give
   !! 1. units (consistent with posang, -3 for mm, 0 for m)
   !! 2. coordinate system (after posang, q.v.)
   !! 3. quantised (if 1, else 0)
   !! 4. number of transforms applied - number of inverse transforms applied
     integer(ki2par), dimension(4)  :: nparpos =(/-3,0,0,0/) !< .
     integer(ki2par)  :: numrparpos=1 !< number of entries in rparpos
     real(kr4), dimension(1) :: rparpos=(/0._kr4/) !< real parameter(s) (NOT USED)
  end type posveclis_t

  type, public :: posvecb_t
     integer(ki2), dimension(3) :: posb !< binary position vector
  end type posvecb_t

  type, public :: quantfm_t
     real(kr4), dimension(3) :: hmin !< vector
     real(kr4), dimension(3) :: rhmin !< reciprocal of vector
     real(kr4), dimension(3) :: offvec !< offset vector
     integer(ki4) :: nqtfm !< type of quantum vector transform
  end type quantfm_t

  type, public :: tfmdata_t
     real(kr4), dimension(3) :: scale !< scaling vector of transform
     real(kr4), dimension(3) :: offset !< offset vector
     real(kr4), dimension(3,3) :: matrix  !< maxtrix of transform
     integer(ki4) :: ntfm !< type of position vector transform
   !    character(len=20) :: id !< identifies transformation
  end type tfmdata_t


!  type, public :: alltfmdata_t
!     type(quantfm) :: qtfm !! quantum scaling transform
!     type(tfmdata) :: tfm !! scaling transform
!  end type alltfmdata_t

!  type, public :: box_t
!     real(kr4), dimension(3) :: posvec !! position vector
!  end type box_t
!
end module position_h
