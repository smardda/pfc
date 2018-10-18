module pcontrol_h

  use const_kind_m
  use termplane_h


! public types

!> numerical control parameters for all power calculation
  type, public :: pnumerics_t
     real(kr8) :: f !< power split ion to electron direction
     real(kr8) :: lmid !< power decay length at outer midplane (m)
     real(kr8) :: ploss !< power crossing the separatrix (W)
     real(kr8) :: sigma !< diffusion length (m)
     integer(ki4)  :: nlevel !< refinement level
     integer(ki4)  :: shadow !< shadowing calculation if positive
     real(kr8) :: qpara0 !< parallel reference flux (W/sqm)
   ! additions for powcal
     character(len=80) :: caltype !< type of calculation ('afws','local', 'msus','global', 'msum','middle')
     integer(ki4)  :: mtrack !< maximum number of tracks to follow (<= hardwired parameter)
     integer(ki4)  :: ntrack !< if non-zero, number of tracks to follow
     integer(ki4), dimension(:), allocatable  :: trackno !< array containing track element numbers
     integer(ki4)  :: nanalau !< type of analytic definition of launch
     integer(ki4)  :: nlaupow !< calculate power deposited on launch geometry
     integer(ki4)  :: nshapow !< calculate power deposited on shadowing geometry
     logical  :: ltermplane !< termination planes present
     logical  :: ledgprof !< more edge profile functions available
     logical  :: lskyl !< use skylight(s) if true
     logical  :: lflagtol !< use changed skylight_tolerance in preference to tolerance from beqparameters
     real(kr8), dimension(3) :: tolin  !< input as skylight_tolerance
     logical  :: lpfpower !< if .TRUE., power in private flux region
     logical  :: lskylpfr !< if .TRUE., use skylights in private flux region
     type(termplane_t) :: termp !< termination planes
  end type pnumerics_t

  type, public :: pfiles_t
   !! file names
     character(len=80)  :: powcalout       !< output data
     character(len=80)  :: log            !< log file
     character(len=80)  :: geoq         !< geoq output format file name
     character(len=80)  :: vtkdata         !< vtk input data file
     character(len=80)  :: hdsdata         !< hds input data file
     character(len=80)  :: laudata         !< launch data
     character(len=80)  :: vtkres         !< vtk input res file
     character(len=80)  :: cartv   !< DUPLICATE  vtk plot file of power per big triangle in Cartesians
     character(len=80)  :: powstatx   !< vtk plot file of power per big triangle in Cartesians
     character(len=80)  :: allcartv   !< DUPLICATE  vtk plot file of power for all triangles in Cartesians
     character(len=80)  :: powx   !< vtk plot file of power for all triangles in Cartesians
     character(len=80)  :: wall   !< vtk plot file of power for shadowing wall triangles in Cartesians
     character(len=80)  :: ptzv    !< INERT vtk plot file of all powers in \f$ (\psi,\theta,\zeta) \f$
     character(len=80)  :: flincart  !< plot file for field line following
     character(len=80)  :: flinx  !< plot file for field line following
     character(len=80)  :: flinends  !< gnuplot file for field line ends
     character(len=80)  :: flinptz  !< plot file for field line following in \f$ (\psi,\theta,\zeta) \f$
     character(len=80)  :: flinm  !< plot file for field line following in \f$ (\psi,\theta,\zeta) \f$
     character(len=80)  :: allptzv    !< INERT vtk plot file of powers in \f$ (\psi,\theta,\zeta) \f$
     character(len=80)  :: geofld  !< INERT vtk plot file of geometry and field
     character(len=80)  :: geofldq  !< INERT vtk plot file of quantised geometry and field
     character(len=80)  :: gnu !< gnuplot file TBD
  end type pfiles_t

  type, public :: pplots_t
   !! vtk plot output selectors
     logical  :: cartv   !< DUPLICATE  power in Cartesians
     logical  :: powstatx   !< power in Cartesians
     logical  :: allcartv   !< DUPLICATE  all power in Cartesians
     logical  :: powx   !< all power in Cartesians
     logical  :: wall   !< wall power in Cartesians
     logical  :: ptzv    !< INERT vtk plot file of powers in \f$ (\psi,\theta,\zeta) \f$
     logical  :: flincart    !< plot file for field line following
     logical  :: flinx    !< plot file for field line following
     logical  :: flinends    !< gnuplot field line endpoints
     logical  :: flinptz    !< plot file for field line following in \f$ (\psi,\theta,\zeta) \f$
     logical  :: flinm    !< plot file for field line following in \f$ (\psi,\theta,\zeta) \f$
     logical  :: allptzv    !< INERT vtk plot file of all powers in \f$ (\psi,\theta,\zeta) \f$
     logical  :: geofld  !< INERT vtk plot file of geometry and field
     logical  :: geofldq  !< INERT vtk plot file of quantised geometry and field
     logical  :: gnu !< gnuplot file TBD
  end type pplots_t

end module pcontrol_h
