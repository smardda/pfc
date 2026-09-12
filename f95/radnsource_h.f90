module radnsource_h

  use const_kind_m
  use fmesh_h
  use fmesh_m
  use spl2d_m

! public types

!> parameters describing how to construct object
  type, public :: radnumerics_t
     character(len=80) :: formula = '' !< radiation formula
     integer(ki4) :: nrpams = 0 !< number of real parameters
     integer(ki4) :: nipams = 0 !< number of integer parameters
     real(kr8), dimension(:), allocatable   :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable   :: npar !< general integer parameters

     real(kr8) :: zcen = 0.0_kr8 !< \f$ Z_{cen} \f$
     real(kr8) :: rcen = 0.0_kr8  !< \f$ R_{cen} \f$
     real(kr8) :: psicen = 0.0_kr8!< value of psi at centre of plasma
     real(kr8) :: psibdry = 0.0_kr8 !< value of psi at plasma boundary

   !real(kr8) :: rcminor !< parameter for circular test formula
   !real(kr8) :: bcpar !< parameter for circular test formula
     real(kr8) :: apar = 0.0_kr8 !< parameter for analytic formula
     real(kr8) :: bpar = 0.0_kr8 !< parameter for analytic formula

   !! physics radiation mechanism parameters
     real(kr8) :: ne = 0.0_kr8  !< \f$ N_e \f$, number density of electrons in plasma (in 1/m^3)
     real(kr8) :: anpar = 0.0_kr8  !< \f$ a_n \f$, parameter for electron density analytic formula
     real(kr8) :: bnpar = 0.0_kr8  !< \f$ b_n \f$, parameter for electron density analytic formula
     real(kr8) :: te = 0.0_kr8  !< \f$ T_e \f$, temperature of electrons (in eV)
     real(kr8) :: atpar = 0.0_kr8  !< \f$ a_T \f$, parameter for electron temperature analytic formula
     real(kr8) :: btpar = 0.0_kr8  !< \f$ b_T \f$, parameter for electron temperature analytic formula
     real(kr8) :: bmag = 0.0_kr8  !< \f$ |B| \f$, magnitude of toroidal magnetic field in Tesla
     logical :: lte = .false. !< specify whether or not to use lte assumptions within code
     integer(ki4) :: nswbrems =0 !< switch for Bremsstrahlung (0=off, 1=on)
     integer(ki4) :: nswcyclo = 0!< switch for cyclotron (0=off, 1=on)
     integer(ki4) :: nswrecomb = 0!< switch for recombination (0=off, 1=on) (NOT USED)

   !! restrict source to box region of plasma
     logical :: lplasbox = .false. !< plasma box set
     real(kr8), dimension(4) :: plasbox = 0.0_kr8 !< plasma box corners \f$ (R_1,Z_1), (R_2,Z_2) \f$

     logical :: radndia = .false.!< switch for radiation diagnostics
     real(kr8) :: powtot = 0.0_kr8  !< total radiation power leaving plasma (positive is user specified)
     real(kr8) :: q0 = 0.0_kr8  !< total power normalisation factor (not used)
     integer(ki4) :: nraorder = 0 !< order of spline used in radn source
     type (fmesh_t) :: fmesh  !< fmesh

  end type radnumerics_t

! type which defines/instantiates the object
  type, public :: radnsource_t

   !! radnsource parameters
     type(spl2d_t) :: srcspl2d  !< spl2d representation of radn source
     real(kr8) :: normfac = 0.0_kr8  !< power normalisation factor (NOT USED)
     real(kr8) :: normfac2 = 0.0_kr8  !< power normalisation factor (NOT USED)
     real(kr8) :: prad = 0.0_kr8  !< total radiative power lost from plasma (NOT USED)
     integer(ki4) :: nspecies = 0 !< number of species in plasma to sum over (NOT USED)
     real(kr8),allocatable,dimension(:) :: noz !< \f$ N(Z) \f$, number density of particular ion species (NOT USED)
     real(kr8),allocatable,dimension(:) :: ezinf !< \f$ E^Z_inf \f$, ionization energies different species (NOT USED)
     real(kr8),allocatable,dimension(:) :: zspecies !< \f$ Z \f$, charge of each species in plasma (NOT USED)
     type(radnumerics_t) :: n !< control  parameters
  end type radnsource_t

! type which defines/instantiates the object
  type, public :: adas_data

   !! adas data parameters
     integer(ki4) :: n_impurities = 0 !< number of impurities in the plasma
     integer(ki4) :: n_temps = 0 !< number of temperatures for loss function data
     integer(ki4) :: dwell_time = 0 !< number of different times used in adas data
     real(kr8),allocatable,dimension(:) :: temps !< temperature data from adas files
     real(kr8),allocatable,dimension(:,:) :: loss_fun !< loss function data from adas files
     Character(len=2), dimension(:), allocatable :: imp_pre !< impurity prefix of file names
  end type adas_data

  type, public :: profile_data
   !! profile data parameters
     real (kind=8),dimension(10001) :: psi = 0.0_kr8 !< psi for defining density and temperature profiles
     real (kind=8),dimension(10001,14) :: density_profile = 0.0_kr8 !< density profile for each impurity
     real (kind=8),dimension(10001,14) :: temperature_profile = 0.0_kr8 !< temperature profile for each impurity
  end type profile_data

end module radnsource_h
