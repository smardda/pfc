module skyl_h

  use const_kind_m
  use dcontrol_h

  implicit none

!> parameters describing how to construct object
  type, public :: sknumerics_t
   !> extent in flux terms
   !! 1. lower, inner \n
   !! 2. lower, outer \n
   !! 3. upper, inner \n
   !! 4. upper, outer \n
     real(kr8), dimension(4) :: psiexts !< -
   !> choose how to calculate psiexts
   !! 1. lower, inner \n
   !! 2. lower, outer \n
   !! 3. upper, inner \n
   !! 4. upper, outer \n
     integer(ki4), dimension(4) :: extsopt !< -
   !> number of binning boxes in calculation of skylight(s)
   !! 1. lower, inner \n
   !! 2. lower, outer \n
   !! 3. upper, inner \n
   !! 4. upper, outer \n
     integer(ki4), dimension(4) :: nexts !< -
   !> size of binning boxes in calculation of skylight(s) (INERT as input)
   !! 1. lower, inner \n
   !! 2. lower, outer \n
   !! 3. upper, inner \n
   !! 4. upper, outer \n
     real(kr8), dimension(4) :: binsize !< -
   !> margin for flux extent calculation of skylight(s)
     real(kr8) :: extsdel !< -
     real(kr8) :: ngeoml !< number of typical lengthscales in sample extent
   !> encodes which if either skylight needed
   !! 0. no skylight (warning condition)
   !! 1. only lower divertor
   !! 2. lower and upper divertor
   !! 3. only upper divertor
     integer(ki4) :: skyltyp !< -
     logical :: lrext !< setting extent in \f$ R \f$
     logical :: lzext !< setting extent in \f$ Z \f$
     real(kr8) :: rextmin !< minimum extent in \f$ R \f$
     real(kr8) :: rextmax !< maximum extent in \f$ R \f$
     real(kr8) :: zextmin !< minimum extent in \f$ Z \f$
     real(kr8) :: zextmax !< maximum extent in \f$ Z \f$
     real(kr8), dimension(3) :: toli !< intersect tolerance vector
     integer(ki4) :: div !< number of divisions in angle
   !> special controls for skylight
   !! first entry=1, use launch points instead of geometry points
   !! second entry=1, force inner \f$ Z(\psi) \f$ monotone
   !! third entry=1, force outer \f$ Z(\psi) \f$ monotone
     integer(ki4), dimension(3) :: control !< -
     logical, dimension(2) :: lprovis !< provisional skylight output required, entry 1=>lower, 2=>upper
  end type sknumerics_t

! type which defines/instantiates the object
  type, public :: skyl_t
     real(kr8) :: geoml !< typical lengthscale of geometry
   !> limits on \f$ \psi \f$  for each of 4 = 2 x 2 skylights \n
   !! first index gives min and max \n
   !! second  index gives inner and outer \n
   !! third  index gives lower and upper
     real(kr8), dimension(2,2,2) :: psilts !< -
   !> \f$ \Delta \psi \f$  for each of 4 = 2 x 2 skylights \n
   !! first  index gives inner and outer \n
   !! second  index gives lower and upper
     real(kr8), dimension(2,2) :: psidelta !< -
     real(kr8), dimension(:,:), allocatable :: inboxr !< inner \f$ \psi \f$-box limits in \f$ R \f$
     real(kr8), dimension(:,:), allocatable :: inboxz !< inner \f$ \psi \f$-box limit in \f$ Z \f$
     real(kr8), dimension(:,:), allocatable :: ouboxr !< outer \f$ \psi \f$-box limits in \f$ R \f$
     real(kr8), dimension(:,:), allocatable :: ouboxz !< outer \f$ \psi \f$-box limit in \f$ Z \f$
   !> dimensions (1) of inboxr/z, (2) of ouboxr/z
     integer(ki4), dimension(2,2) :: dimbox !< -
     real(kr8) :: eps !< smallest distinguished lengthscale of geometry
     integer(ki4) :: nprovis1 !< lower provisional skylight dimension
     real(kr8), dimension(:,:), allocatable :: provis1 !< lower provisional skylight
     integer(ki4) :: nprovis2 !< upper provisional skylight dimension
     real(kr8), dimension(:,:), allocatable :: provis2 !< upper provisional skylight
     type(sknumerics_t) :: n !< control parameters
     type(dnumerics_t) :: dn !< control parameters from dnumerics
     integer(ki4), dimension(:), allocatable :: ntest !< test
     integer(ki4) :: debug !< output to debug file units if unity
     integer, dimension(10) :: ndskyl !< debug file units
     integer(ki4) :: ndskyln !< number of debug file units
  end type skyl_t

end module skyl_h
