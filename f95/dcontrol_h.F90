module dcontrol_h

  use const_kind_m

!! public types

!> numerical run parameters
  type, public :: dnumerics_t
     character(len=80)  :: tfm       !< transformation
     real(kr8), dimension(:), allocatable :: r !< positions in 1 coordinate
     real(kr8), dimension(:), allocatable :: z !< positions in 2 coordinate
     integer(ki4) :: npos !< number of positions
     integer(ki4) :: div !< number of divisions in rotate/translate (should be even)
   !> number of divisions in straight line joining first
   !! and last \f$ (R,Z) \f$
     integer(ki4) :: ldiv !< -
     integer(ki2par) :: descode !< description of surface
     integer(ki2par) :: csys !< coordinate system positive as posang (q.v.), unity for polars
     integer(ki2par) :: cunits !< units of system (-3 for mm)
     real(kr8) :: stang  !< starting angle for surface generation
     real(kr8) :: finang  !< finishing angle for surface generation
     real(kr8), dimension(3) :: stpos !< starting position for surface generation
     real(kr8), dimension(3) :: finpos !< finishing position for surface generation
  end type dnumerics_t

!> file names
  type, public :: dfiles_t
     character(len=80)  :: out       !< output data
     character(len=80)  :: log           !< log file
     character(len=80)  :: rfile         !< file with position coordinate 1 data
     character(len=80)  :: zfile         !< file with position coordinate 2 data
     character(len=80)  :: rzfile         !< file with position coordinate 1 & 2 data
     character(len=80)  :: vtk   !< vtk file
     character(len=80)  :: gnu !< gnuplot file
  end type dfiles_t

!> plot output selectors
  type, public :: dplots_t
     logical  :: vtk   !< vtk plot selector
     logical  :: gnu !< gnuplot plot selector
  end type dplots_t

end module dcontrol_h
