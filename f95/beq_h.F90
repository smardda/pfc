module beq_h

  use const_kind_m
  use fmesh_h
  use spl2d_m
  use spl3d_m

! public types

!> numerical control parameters for beq
  type, public :: bnumerics_t
   !> number of points for gnuplot
   !! display in direction of invariant coordinate \f$ \zeta \f$,
   !! if equal to unity, output plane at \f$ \zeta=0 \f$
     integer(ki4) :: mzetag !< -
   !> number of points for vtk
   !! display in direction of invariant coordinate \f$ \zeta \f$,
   !! if equal to unity, output plane at \f$ \zeta=0 \f$
     integer(ki4) :: mzetav !< -
   !> \f$ C_{opt} \f$. Set position of plasma centre \f$ R_{cen}, Z_{cen} \f$
   !! 1. User input of \f$ R_{cen}, Z_{cen} \f$
   !! 2. Find mesh-point with extremal psi value, start search at
   !!  \f$ (R_{\min}+R_{\max})/2, (Z_{\min}+Z_{\max})/2 \f$
   !! 3. Find mesh-point with extremal psi value, start search at
   !!  user-specified position rcen,zcen
   !! 4. Use values rqcen,zqcen from EQDSK file
     integer(ki4) :: cenopt !< -
   !> \f$ \psi{opt} \f$. Set range of flux to be used in flux mapping (if flux-mapped case)
   !! 1. User input of normalised \f$ \psi{\min}\f$ and \f$ \psi_{\max}\f$ for mapping.
   !! normalised by X-point or limiter value, no margin \f$ \delta \psi\f$.
   !! 2. \f$ \psi{\min}\f$ and \f$ \psi_{\max}\f$ are calculated from input geometry,
   !! with margin \f$ \delta \psi = \f$ beq_delpsi \f$ \times (\psi{\max}-\psi_{\min}) \f$
   !! 3. User input of \f$ \psi{\min}\f$ and \f$ \psi_{\max}\f$ for mapping.
     integer(ki4) :: psiopt !< -
   !> \f$ m_{opt} \f$. Set value of flux at plasma boundary psim.
   !! rsig is important rsig>0 -> plasma centre at a minimum of psi
   !! 1. User input of psim=psiref
   !! 2. Use value psiqbdry from EQDSK file unless psiltr is tighter
   !! 3. psim=psiltr calculated using silhouette of geometry (sampled with
   !! spatial rate deltal=10mm) unless psiqbdry is tighter.
   !! Tighter implies psim = max(psix, psiqbdry) if rsig<0, vv. if rsig>0
   !! 4. psim=psix is used psiqbdry is tighter.
   !! Tighter implies psim = max(psix, psiqbdry) if rsig<0, vv. if rsig>0
   !! 5. psim=psiref with inboard limiting assumed
   !! 6. Invalid
   !! 7. Invalid
   !! 8. as 4.
   !! 9. psim=psiref with outboard limiting assumed
     integer(ki4) :: bdryopt !< -
   !> \f$ N_{opt} \f$. Set numbers of points  \f$ N_{\psi}, N_{\theta} \f$ to be used in flux mapping
   !! 1. User input of \f$ N{\min}\f$ and \f$ N_{\max}\f$
   !! 2. Estimate sensible values of \f$ N_{\psi}, N_{\theta} \f$ proportional
   !! to \f$ N_R, N_Z \f$ from equilibrium file
     integer(ki4) :: nopt !< -
   !> \f$ \theta_{opt} \f$. Set range of \f$ \theta \f$ to be used in flux mapping.
   !! Origin of theta may be set using thetaref
   !! 1. User input of \f$ \theta{\min}\f$ and \f$ \theta_{\max}\f$ in coordinate
   !! system with zero at X-point, or relative to vertical, default unit radians,
   !! outboard is positive, no margin \f$ \delta \theta\f$.
   !! 2. \f$ \theta{\min}\f$ and \f$ \theta_{\max}\f$ are calculated from input geometry,
   !! with margin \f$ \delta \theta = \f$ beq_deltheta \f$ \times (\theta{\max}-\theta_{\min}) \f$
     integer(ki4) :: thetaopt !< -
   !> \f$ \zeta_{opt} \f$. Set range of \f$ \zeta \f$ to be used in mapping to cylindricals
   !! 1. User input of \f$ \zeta{\min}\f$ and \f$ \zeta_{\max}\f$, must be
   !! positive angles, default unit radians. Note \f$ \zeta = -\phi\f$.
   !! 2. \f$ \zeta{\min}\f$ and \f$ \zeta_{\max}\f$ are calculated as
   !! \f$ -\pi/m \f$ and \f$  +\pi/m \f$
   !! using m=nzetap (if nzetap>0), or m=mrip (if mrip>0) or m=1 (if nzetap=mrip=0).
   !! 3. \f$ \zeta{\min}\f$ and \f$ \zeta_{\max}\f$ are calculated as
   !! \f$ -\pi \f$ and \f$  +\pi \f$ regardless
   !! 4. \f$ \zeta{\min}\f$ and \f$ \zeta_{\max}\f$ are calculated as
   !! \f$ -3\pi/2 \f$ and \f$  +\pi/2 \f$ regardless
   !! 5. \f$ \zeta{\min}\f$ and \f$ \zeta_{\max}\f$ are calculated as
   !! \f$ -\pi/2 \f$ and \f$  +3\pi/2 \f$ regardless
     integer(ki4) :: zetaopt !< -
   !> \f$ \xi_{opt} \f$. Set range of \f$ \xi \f$ to be used in mapping to cylindricals
   !! 1. \f$ \xi{\min}\f$ and \f$ \xi_{\max}\f$ = \f$ \zeta{\min}\f$ and \f$ \zeta_{\max}\f$
   !! 2. \f$ \xi{\min}\f$ and \f$ \xi_{\max}\f$ are calculated as
   !! \f$ -\pi \f$ and \f$  +\pi \f$, also  nzets=nzetp (if nzetp>0) or nzets=mrip
   !! 3. Invalid
   !! 4. \f$ \xi{\min}\f$ and \f$ \xi_{\max}\f$ are calculated as
   !! \f$ -3\pi/2 \f$ and \f$  +\pi/2 \f$ regardless, also  nzets=nzetp (if nzetp>0) or nzets=mrip
   !! 5. \f$ \xi{\min}\f$ and \f$ \xi_{\max}\f$ are calculated as
   !! \f$ -\pi/2 \f$ and \f$  +3\pi/2 \f$ regardless, also  nzets=nzetp (if nzetp>0) or nzets=mrip
     integer(ki4) :: xiopt !< -
     real(kr8) :: rcen !< \f$ R_{cen} \f$
     real(kr8) :: zcen !< \f$ Z_{cen} \f$
     real(kr8) :: psimin !< \f$ \psi_{\min} \f$
     real(kr8) :: psimax !< \f$ \psi_{\max} \f$
     integer(ki4) :: npsi !< \f$ N_{\psi} \f$
     real(kr8) :: delpsi !< \f$ \delta\psi \f$
     real(kr8) :: thetamin !< \f$ \theta_{\min} \f$
     real(kr8) :: thetamax !< \f$ \theta_{\max} \f$
     real(kr8) :: deltheta !< \f$ \delta\theta \f$
     real(kr8) :: zetamin !< \f$ \zeta_{\min} \f$
     real(kr8) :: zetamax !< \f$ \zeta_{\max} \f$
     real(kr8) :: ximax !< maximum \f$ \xi \f$ value
     real(kr8) :: ximin !< minimum \f$ \xi \f$ value
     real(kr8), dimension(3,2) :: xbb !< limits on quantised vectors
     real(kr8) :: rmove !< \f$ R_{mov} \f$ displacement of equilibrium
     real(kr8) :: zmove !< \f$ Z_{mov} \f$ displacement of equilibrium
     real(kr8) :: fscale !< Scale toroidal field \f$ f \f$
     real(kr8) :: psiref !< \f$ \psi_X \f$
     real(kr8) :: thetaref !< \f$ \theta_X \f$
     integer(ki4) :: ntheta !< \f$ N_{\theta} \f$
     integer(ki4) :: nzetp !< \f$ N_{\zeta P} \f$
     logical :: skipb !< read B cpts from eqdsk file unless .TRUE.
     logical :: eqbdry !< produce files containing boundary and limiter points from eqdsk
     character(len=80) :: eqbdryfile !< file containing boundary points from eqdsk
     character(len=80) :: eqltrfile !< file containing limiter points from eqdsk
   !> specifies way equilibrium field described
   !! 1. R/J implies mapped case (axisymmetric) -> move
   !! 2. psi gradients
   !! 3. R/J x psi gradients (axisymmetric) -> move1
     integer(ki4) :: fldspec !< -
     integer(ki4):: fiesta !< flag equil file in modified FIESTA output format
     integer(ki4) :: psibig !< normally 0, 1 if \f$ \psi=\Psi \f$ has \f$ 2\pi \f$ factor
     integer(ki4) :: xsearch !< how to search for X-point (1 = within box)
     real(kr8) :: xrsta !< X-point box min R (m)
     real(kr8) :: xrend !< X-point box max R (m)
     real(kr8) :: xzsta !< X-point box min Z (m)
     real(kr8) :: xzend !< X-point box max Z (m)
     integer(ki4) :: search !< how to search for contact-point (1 = within box)
     real(kr8) :: lkrsta !< contact point box min R (m)
     real(kr8) :: lkrend !< contact point box max R (m)
     real(kr8) :: lkzsta !< contact point box min Z (m)
     real(kr8) :: lkzend !< contact point box max Z (m)

   !> type of equilibrium
   !! - 'eqdsk' for EQDSK input
   !! - 'equ' for FIESTA input
   !! - 'ana' for analytic
     character(len=80) :: eqtype !< -
     character(len=80) :: vacfile !< specifies file for spl3d format vacuum field
     integer(ki4) :: mrip !< \f$ N \f$ number of ripple coils
     real(kr8) :: irip !< unused parameter of ripple coils
     real(kr8) :: arip !< \f$ a \f$ for ripple coils
     logical :: leqok !< equilibrium helicity is ok, no need for override
     logical :: duct !< flag whether work in duct coordinates
     integer(ki4), dimension(0:10) :: objadd !< number of objects of given description to add to geobjlist
     logical :: skyl !< flag whether skylight present
     logical :: skylpsi !< flag whether flux limited skylight(s) present
     logical :: skylcen !< flag whether any skylight(s) defined by plasma centre line
     integer(ki4) :: skyldbg !< level of output to help understand skylight production
  end type bnumerics_t

!> data structure describing equilibrium field
  type, public :: beq_t
     real(kr8), dimension(:), allocatable :: f !< \f$ f(\psi) \f$ from EFIT
     real(kr8) :: rmin !< \f$ R_{\min} \f$ for equilibrium data
     real(kr8) :: rmax !< \f$ R_{\max} \f$ for equilibrium data
     integer(ki4) :: mr !< \f$ N_R \f$ for equilibrium data
     real(kr8) :: dr !< \f$ (R_{\max}-R_{\min})/N_R \f$
     real(kr8) :: zmin !< \f$ Z_{\min} \f$ for equilibrium data
     real(kr8) :: zmax !< \f$ Z_{\max} \f$ for equilibrium data
     integer(ki4) :: mz !< \f$ N_Z \f$ for equilibrium data
     real(kr8) :: dz !< \f$ (Z_{\max}-Z_{\min})/N_z \f$
     real(kr8), dimension(:), allocatable :: i !< \f$ I(\psi) \f$
     type(spl2d_t) :: psi !< \f$ \psi(R,Z) \f$ from input
     type(spl2d_t) :: dpsidr !< \f$ \frac{\partial\psi}{\partial R}(R,Z) \f$
     type(spl2d_t) :: dpsidz !< \f$ \frac{\partial\psi}{\partial Z}(R,Z) \f$
     type(spl2d_t) :: rispldr !< \f$ \frac{R}{I}\frac{\partial\psi}{\partial R}(R,Z) \f$
     type(spl2d_t) :: rispldz !< \f$ \frac{R}{I}\frac{\partial\psi}{\partial Z}(R,Z) \f$
     type(spl2d_t) :: r !< \f$ R(\psi,\theta) \f$
     type(spl2d_t) :: z !< \f$ Z(\psi,\theta) \f$
     type(spl2d_t) :: rjac !< \f$ R/J(\psi,\theta) \f$
     type(bnumerics_t) :: n !< control  parameters
     real(kr8) :: dpsi !< \f$ (\psi_{\max}-\psi_{\min})/N_\psi \f$
     real(kr8) :: psiaxis !< \f$ \psi \f$ on axis
     real(kr8) :: psibdry !< \f$ \psi \f$ at reference boundary (limiter, X-point or user)
     real(kr8) :: psiltr !< \f$ \psi \f$ at limiter
     real(kr8) :: rxpt !< \f$ R \f$ at X-point (\f$ R=0 \f$ implies no X-point defined)
     real(kr8) :: zxpt !< \f$ Z \f$ at X-point
     real(kr8) :: psixpt !< \f$ \psi \f$ at X-point
     real(kr8) :: thetaxpt !< \f$ \theta \f$ at X-point
     real(kr8) :: psiotr !< Extreme value of \f$ \psi \f$ for geometry which is not limiter extremum
     real(kr8) :: psinorm !< normalisation of \f$ \psi \f$, i.e. a representative magnitude
     real(kr8) :: psiqbdry !< \f$ \psi \f$ at boundary set by EQDSK
     real(kr8) :: thetagmax !< maximum \f$ \theta \f$ value of geometry
     real(kr8) :: thetagmin !< minimum \f$ \theta \f$ value of geometry
     real(kr8) :: rqcen !< \f$ R \f$ at central extremum of \f$ \psi \f$
     real(kr8) :: zqcen !< \f$ Z \f$ at central extremum of \f$ \psi \f$
     real(kr8) :: rbdry !< \f$ R \f$ at reference boundary
     real(kr8) :: bpbdry !< \f$ B_p \f$ at reference boundary, \f$ B_{pm} \f$
     real(kr8) :: btotbdry !< Total \f$ B \f$ at reference boundary \f$ B_m \f$
     real(kr8) :: ximaxm !< maximum \f$ \xi \f$ value \f$ -\epsilon \f$
     real(kr8) :: ximinp !< minimum \f$ \xi \f$ value \f$ +\epsilon \f$
     real(kr8), dimension(3) :: domlen !< dimensions of domain (quantised)
     integer(ki4) :: nzets !< \f$ N_{\zeta s} \f$ such that \f$ \xi = N_{\zeta s} \zeta \f$
     integer(ki4) :: ntmin !< minimum valid value of \f$ N_\theta \f$
     integer(ki4) :: ntmax !< maximum valid value of \f$ N_\theta \f$
     real(kr8) :: dtheta !< \f$ (\theta_{\max}-\theta_{\min})/N_\theta \f$
     real(kr8),dimension(:),allocatable :: srmin !< \f$ r_{\min}(\theta) \f$, set of \f$ r \f$  corresponding to \f$ \psi_{\min} \f$
     real(kr8),dimension(:),allocatable :: srmax !< \f$ r_{\max}(\theta) \f$, set of \f$ r \f$  corresponding to \f$ \psi_{\max} \f$
     real(kr8) :: ivac !< \f$ I \f$ for vacuum field
     real(kr8) :: psicen !<  Computed \f$ \psi \f$ on axis, compare PSIAXIS from eqdsk
     logical :: replasi !< If .TRUE. geoq is replacing psiaxis and psiqbdry from EQDSK with its estimates
     type(spl3d_t)  :: vacfld  !< vacuum field structure
     type(fmesh_t)  :: fmesh  !< field mesh for calculation in duct coordinates
     real(kr8), dimension(:,:), allocatable :: ctrackrz !< central track, of psi extremum through plasma centre
     integer(ki4) :: nctrack !<  bound for ctrackrz array
     real(kr8) :: rcenq !< quantised \f$ R_{cen} \f$
     real(kr8) :: zcenq !< quantised \f$ Z_{cen} \f$
  end type beq_t

! public variables
  integer(ki4), public, parameter :: beq_spline_order=4 !< set to four for cubic splines (use of tc01a assumes cubic splines)
  logical, public :: beq_nobinq !< flag whether \f$ B \f$ data in EQDSK file
  logical, parameter :: BEQ_OVERRIDE_ITER=.TRUE. !< use values required by ITER project (unless equil_helicity_ok set)
  logical, parameter :: BEQ_OVERRIDE_FTU=.FALSE. !< use values required by FTU project (unless equil_helicity_    ok set)

end module beq_h
