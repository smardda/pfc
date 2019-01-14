module beqan_h

  use const_kind_m
  use const_numphys_h

! public types

! type storing analytic equilibrium description
  type, public :: eqparms_t
     character(len=80) :: formula !< equilibrium formula
     real(kr8) :: psi0 !< flux constant
     real(kr8) :: bphi !< toroidal field
     real(kr8) :: r1 !< toroidal inner radius
     real(kr8) :: r2 !< toroidal outer radius
     real(kr8) :: rm !< \f$ R \f$ of point of max height
     real(kr8) :: zm !< maximum \f$ Z \f$
     integer(ki4) :: nrpams !< number of real parameters
     integer(ki4) :: nipams !< number of integer parameters
     real(kr8), dimension(:), allocatable :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable :: npar !< general integer parameters
  end type eqparms_t

! type storing eqparms as well as other variables
  type, public :: beqan_t
     type(eqparms_t) :: eqparms !< type variable
     real(kr8) :: a !< minor radius
     real(kr8) :: rxsq !< separatrix \f$ R \f$ squared
     real(kr8) :: esq !< \f$ E^2 \f$
     real(kr8) :: ar !< cross section
     real(kr8) :: s !< torus surface
     real(kr8) :: v !< volume
     real(kr8) :: r0sq !< magnetic axis squared
     real(kr8) :: k !< elongation
     real(kr8) :: eps !< inverse aspect ratio
     real(kr8) :: del !< triangularity
     real(kr8) :: iphi !< toroidal current
     real(kr8) :: rc !< type variable
     real(kr8) :: psiaxis !< \f$ \psi \f$ on axis
     real(kr8) :: psibdry !< \f$ \psi \f$ at boundary
     real(kr8) :: rmin !< \f$ R_{\min} \f$
     real(kr8) :: rmax !< \f$ R_{\max} \f$
     integer(ki4) :: mr !< \f$ N_R \f$
     real(kr8) :: zmin !< \f$ Z_{\min} \f$
     real(kr8) :: zmax !< \f$ Z_{\max} \f$
     integer(ki4) :: mz !< \f$ N_Z \f$
  end type beqan_t

end module beqan_h

