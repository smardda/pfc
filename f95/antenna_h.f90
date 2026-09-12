module antenna_h

  use const_kind_m
  use position_h
  use geobjlist_h

!> parameters describing how to construct object

  type, public :: antnumerics_t
     integer(ki4) :: ntheta = 0 !< number_of_samples_in_theta
     integer(ki4) :: nphi = 0 !< number_of_samples_in_phi
     character(len=80) :: formula = '' !< antenna variable
     real(kr8) :: f = 0.0_kr8 !< power split (math variable name allowed)
     integer(ki4) :: nrpams = 0 !< number of real parameters
     integer(ki4) :: nipams = 0 !< number of integer parameters
     real(kr8), dimension(:), allocatable   :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable   :: npar !< general integer parameters
     logical :: plot_diagvtk  = .false.!< logical to state whether or not to produce paraview visualisation of antenna
  end type antnumerics_t

  type, public :: antenna_t
     real(kr8) :: totpow = 0.0_kr8  !< total power on antenna
     type(posvecl_t) :: tricen !< position vector of barycentre of triangle
     type(posvecl_t) :: trinorm !< normal vector to triangle
     type(geobjlist_t) :: hemisphere !< hemisphere data
     type(antnumerics_t) :: n !< control  parameters
     real(kr4), dimension(:), allocatable :: area !< area of each triangle in hemisphere
     type(posvecl_t), dimension(:), allocatable :: dir !< launch directions from triangulated hemisphere normals
     type(posvecl_t), dimension(:), allocatable :: posdir !< launch position of directions
  end type antenna_t

!! currently ray object data is just stored within antenna_t as current testing is all within antenna modules

end module antenna_h
