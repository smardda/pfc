module posang_h

  use const_kind_m

! public types

!> vector (in polar coordinates) manipulation
  type, public :: posang_t
     real(kr8), dimension(3) :: vec !< vector (in polars by default)
     real(kr8), dimension(3) :: pos !< position vector (in polars by default)
   !> format currently used for vectors
   !! 0 - position in Cartesians \f$ (x,y,z) \f$
   !! 1 - position in cylindrical polars \f$ (R,Z,\zeta), |\zeta| \leq \pi \f$
   !! 2 - position in flux coordinates  \f$ (\psi,\theta,\zeta) \f$
   !! 16 - position and vector in Cartesians
   !! 17 - position and vector in cylindrical polars
   !! 18 - position and vector in flux coordinates \f$ (0,B_\theta,B_\zeta) \f$
     integer(ki4) :: opt !< local variable
     integer(ki4) :: units !<  (-3) for mm, 0 otherwise for metres
  end type posang_t

end module posang_h
