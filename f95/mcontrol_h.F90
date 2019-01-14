module mcontrol_h

  use const_kind_m
  use position_h

! public types

  type, public :: mnumerics_t
   !! numerical run parameters
     real(kr8) :: rax !<  \f$ R \f$  value at axis, point which determines vacuum field
     real(kr8) :: ireq !<  \f$ I = B R \f$ value required (determines scaling of input)
     integer(ki4) :: magtfm !< how vacuum field axis is determined
     integer(ki4) :: nord !< spline order to use (always 4)
     integer(ki4) :: n0 !< split for transform (always 2)
     integer(ki4), dimension(3) :: kmin !< minimum wavenumber for evaluation
     integer(ki4), dimension(3) :: kmax !< maximum wavenumber for evaluation
     integer(ki4), dimension(3) :: parity !< parity of component (1=even,2=odd)
     real(kr8) :: cutoff !<  cutoff for maximum mode number calculation
  end type mnumerics_t


! public types
  type, public :: mfiles_t
   !! file names
     character(len=80)  :: mout       !< output data
     character(len=80)  :: log            !< log file
     character(len=80)  :: magdata         !< mag input data file
     character(len=80)  :: magout         !< mag output data file
     character(len=80)  :: cartv   !< \f$ \bf{B} \f$ in Cartesians
     character(len=80)  :: allcartv   !< all \f$ \bf{B} \f$ in Cartesians
     character(len=80)  :: gnuv !< gnuplot file of fields as functions \f$ f(R,Z) \f$
     character(len=80)  :: gnu !< gnuplot file of \f$ fns(\zeta) \f$
     character(len=80)  :: modes !< gnuplot file of maximum mode numbers
  end type mfiles_t

  type, public :: mplots_t
   !! plot output selectors
     logical  :: cartv   !< \f$ \bf{B} \f$ in Cartesians
     logical  :: allcartv   !< all \f$ \bf{B} \f$ in Cartesians
     logical  :: gnuv !< gnuplot file of fields as functions \f$ f(R,Z) \f$
     logical  :: gnu !< gnuplot file of \f$ fns(\zeta) \f$
     logical  :: testtfm !< test Fourier transform in \f$ \zeta \f$ (NO plot)
     logical  :: modes !< gnuplot file of maximum mode numbers
  end type mplots_t


end module mcontrol_h
