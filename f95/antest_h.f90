module antest_h

  use const_kind_m
  use antenna_h


!> parameters describing how to construct object
  type, public :: atnumerics_t
     type(posvecl_t), dimension(:,:), allocatable :: directions !< directions array
     integer(ki4) :: ntheta !< number_of_samples_in_theta
     integer(ki4) :: nphi !< number_of_samples_in_phi
     integer(ki4) :: nsample !< number of total sampled directions from triangle
     integer(ki4) :: nraysample !< total number of sampled points along ray
     real(kr8) :: termdistance !< terminating distance of ray
     integer(ki4) :: option !< case_select_for_sampling_pattern
     character(len=80) :: formula !< antenna variable
     character(len=80) :: pattern !< antenna sampling pattern
   !> type of equilibrium
   !! - 'eqdsk' for EQDSK input
   !! - 'equ' for FIESTA input
   !! - 'ana' for analytic
     character(len=80) :: eqtype !< -
     real(kr8) :: f !< power split (math variable name allowed)
     integer(ki4) :: nrpams !< number of real parameters
     integer(ki4) :: nipams !< number of integer parameters
     real(kr8), dimension(:), allocatable   :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable   :: npar !< general integer parameters
  end type atnumerics_t

  type, public :: antest_t
     real(kr8) :: totpow !< total power
     type(atnumerics_t) :: n !< control  parameters
     type(antenna_t) :: antenna !< antenna variables
     type(geobjlist_t) :: geobjl !< geometrical objects
     real(kr4), dimension(:), allocatable :: pow !< power
  end type antest_t

end module antest_h
