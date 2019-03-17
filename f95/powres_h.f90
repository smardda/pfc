module powres_h

  use const_kind_m
  use geobjlist_h
  use geobj_m
  use position_h
  use fmesh_h
  use beq_h
  use odes_h
  use skyl_h


! public types

!> numerical control parameters for all power calculation
  type, public :: prnumerics_t
     real(kr8) :: f !< power split ion to electron direction
     real(kr8) :: lmid !< power decay length at outer midplane (m)
     real(kr8) :: ploss !< power crossing the separatrix (W)
     integer(ki4)  :: nlevel !< refinement level
     integer(ki4)  :: shadow !< shadowing calculation if positive
     real(kr8) :: qpara0 !< parallel reference flux (W/sqm)
  end type prnumerics_t

!> type storing power results data
  type, public :: powres_t
     type(posveclis_t) :: vecb !< list of B-vectors
     type(posveclis_t) :: vecx !< list of Cartesian position vectors
     type(geobjlist_t) :: geobjl !< geobj list
     type(beq_t) :: beq !< beq (only part defined)
     type(skyl_t) :: skyl !< skylight definition(s)
     integer(ki4)  :: npowe !< number of powelts (ie. level 1 triangles)
     integer(ki4)  :: npow !< number of entries in pow
     integer(ki4)  :: npows !< number of entries in pows (normally use npowe)
     real(kr8) :: rblfac !< Lomas' factor \f$ \frac{1}{2 \pi \lambda_{mid} R_m B_{pm}} \f$
     real(kr8) :: fpfac !< Lomas' factor \f$ \frac{F P_{loss}}}{2 \pi \lambda_{mid} R_m B_{pm}} \f$
     real(kr8) :: slfac !< factor \f$ \frac{\sigma}{2 \lambda_{mid}} \f$
     real(kr8) :: qfac !< quantising factor
     real(kr8), dimension(3) :: qfaca !< quantising factors
     real(kr4), dimension(:), allocatable :: pow !< power
     integer(ki4), dimension(:), allocatable :: pmask !< power mask
     real(kr4), dimension(:), allocatable :: pows !< power statistic
     real(kr4), dimension(:), allocatable :: powa !< power average statistic
     real(kr4), dimension(:), allocatable :: psista !< \f$ \psi \f$ at start of fieldline
     type(prnumerics_t) :: n !< control  parameters
     real(kr8) :: psimin !< minimum value of \f$ \psi \f$ for interpolation
     real(kr8) :: psimax !< maximum value of \f$ \psi \f$ for interpolation
     real(kr8) :: thetamin !< minimum value of \f$ \theta \f$ for interpolation
     real(kr8) :: thetamax !< maximum value of \f$ \theta \f$ for interpolation
     logical  :: flincart    !< plot switch for field line following
     logical  :: flinx    !< plot switch for field line following
     logical  :: flinends    !< gnuplot switch for field line ends
     logical  :: flinptz    !< plot switch for field line following in \f$ (\psi,\theta,\zeta) \f$
     logical  :: flinm    !< plot switch for field line following in \f$ (\psi,\theta,\zeta) \f$
     integer  :: nflends !< gnuplot file unit for field line ends
  end type powres_t

end module powres_h
