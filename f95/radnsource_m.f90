module radnsource_m

  use radnsource_h
  use misc_m
  use log_m
  use const_numphys_h
  use const_kind_m
  use geobjlist_h
  use geobjlist_m

  implicit none
  private

! public subroutines
  public :: &
  radnsource_init, & !< initialise source array from parametric or disc data
  radnsource_openfile,  & !< open file
  radnsource_readcon,  & !< read data from file
  radnsource_parameter,  & !< set radnsource emission function parameters
  radnsource_fn, &  !< general external function call
  radnsource_dia, &  !< object diagnostics to log file
  radnsource_initwrite, & !< open new file, making up name
  radnsource_write, &  !< write out object
  radnsource_writeg, &  !< write out object as gnuplot
  radnsource_writev, &  !< write out object as vtk
  radnsource_delete, & !< delete object
  radnsource_close, & !< close file
  radnsource_closewrite !< close write file
  real (kr8) :: frac_impurity(14) !< fraction of each impurity in plasma (use physics case)
  integer(ki4)  :: prof_custom=0 !< used to define the temperature and density profiles used to each impurity

! private subroutines
  private :: &
  radnsource_initfile, & !< initialise source array from disc data
  radnsource_initformul,& !< initialise source array from parametric data
  radnsource_userdefined,  & !< user-defined function
  radnsource_circular,  & !< user-defined function
  radnsource_physics,  & !< user-defined function
  radnsource_analytic,  & !< user-defined function
  radnsource_norm,  & !< normalisation of radiation source
  radnsource_cyclo, & !< calculate cyclotron radiation power density (watts/m^3)
  radnsource_recomb, & !< calculate recombination (free-bound) radiation power density (watts/m^3)
  radnsource_brems,& !< calculate bremsstrahlung radiation power density (watts/m^3)
  radnsource_adas,& !< subroutine for reading in adas data
  radnsource_profile,& !< subroutine for reading in profile_data
  radnsource_linebrem !< calculate impurity line and bremsstrahlung radiation power density (watts/m^3)

! private variables
  character(*), parameter :: m_name='radnsource_m' !< module name
  real(kr8), allocatable, dimension(:,:,:,:) :: rwork   !< input samples
  integer(ki4)  :: status   !< error status
  integer(ki4), save  :: ninra=5     !< control file unit number
  integer(ki4), save  :: noutra=6      !< output file unit number
  character(len=80) :: icbuff !< buffer for input/output
  character(len=80), save :: radnsourcefile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4)  :: ilog      !< for namelist dump after error
  real(kr8), parameter :: const_halfd=0.5_kr8 !< double precision half
  real(kr8), parameter :: const_twod=2.0_kr8 !< double precision two
  real(kr8), parameter :: const_zerod=0.0_kr8 !< double precision zero
  real(kr8), parameter :: const_oned=1.0_kr8 !< double precision unit
  real(kr8), parameter :: const_epsd=0.1_kr8 !< double precision tenth
  real(kr8), parameter :: const_eps2d=0.01_kr8 !< double precision hundredth
  real(kr8), parameter :: const_eps3d=0.001_kr8 !< double precision thousandth
  real(kr8), parameter :: const_eps4d=1.0e-4_kr8 !< double precision small
  real(kr8), parameter :: const_tend=10.0_kr8 !< double precision ten
  real(kr8), parameter :: const_ten8d=1.0e+8_kr8 !<  small infinity

  character(len=29), parameter :: cfmtbv0c='(g14.7,", ",g14.7,", ",g14.7)'
  character(len=39), parameter :: cfmtbv1c='(g16.7,"",g16.7,"",g16.7)'

  contains
!---------------------------------------------------------------------
!> initialise source array and normalise
subroutine radnsource_init(self,pivac,psi,psicen,psibdry,prcen,pzcen,file)

  !! arguments
  type(radnsource_t), intent(inout) :: self !< module object
  type(spl2d_t), intent(inout) :: psi !< \f$ \psi(R,Z) \f$ from input
  real(kr8), intent(in) :: pivac   !<  \f$ I \f$
  real(kr8), intent(in) :: psicen   !<  \f$ \psi_c \f$
  real(kr8), intent(in) :: psibdry   !<  \f$ \psi_b \f$
  real(kr8), intent(in) :: prcen   !<  \f$ R_c \f$
  real(kr8), intent(in) :: pzcen   !<  \f$ Z_c \f$
  character(*), intent(in) :: file !< file name
  !! local
  character(*), parameter :: s_name='radnsource_init' !< subroutine name

  self%n%psicen=psicen
  self%n%psibdry=psibdry
  self%n%rcen=prcen
  self%n%zcen=pzcen
  self%n%bmag=pivac/prcen

  formula_init: select case (self%n%formula)
  case('file_input')
     call radnsource_initfile(self,psi,file)

  case default
     call radnsource_initformul(self,psi,file)

  end select formula_init

  call radnsource_norm(self,psi)

end subroutine radnsource_init
!---------------------------------------------------------------------
!> radiation source from file
subroutine radnsource_initfile(self,psi,file)

  !! arguments
  type(radnsource_t), intent(inout) :: self !< module object
  type(spl2d_t), intent(inout) :: psi !< \f$ \psi(R,Z) \f$ from input
  character(*), intent(in) :: file !< file name
  !! local
  character(*), parameter :: s_name='radnsource_initfile' !< subroutine name
  integer(ki4) :: in1   !< first nominal dimension
  integer(ki4) :: in2   !< second nominal dimension
  integer(ki4) :: in3   !< third nominal dimension
  integer(ki4)  :: iorder   !< spline order`
  real(kr8) :: zorg1   !< minimum mesh value in coordinate 1
  real(kr8) :: zorg2   !< minimum mesh value in coordinate 2
  real(kr8) :: zh1   !< mesh spacing, coordinate 1
  real(kr8) :: zh2   !< mesh spacing, coordinate 2
  real(kr8) :: zpt1   !< position in coordinate 1
  real(kr8) :: zpt2   !< position in coordinate 2
  real(kr8) :: zpsi   !< flux value

  iorder=self%n%nraorder
  icbuff=' '

  call radnsource_openfile(file)
  if (ninra<0) then
     !! error opening file
     print '("Fatal error: Unable to open radn data file, ",a)',file
     call log_error(m_name,s_name,1,error_fatal,'Cannot open radn data file')
     stop
  end if
  call fmesh_readmlab(self%n%fmesh,ninra)
  in1=self%n%fmesh%nxf
  in2=self%n%fmesh%nyf
  in3=self%n%fmesh%nzf
  zorg1=self%n%fmesh%xf(1)
  zorg2=self%n%fmesh%yf(1)
  zh1=self%n%fmesh%xf(2)-self%n%fmesh%xf(1)
  zh2=self%n%fmesh%yf(2)-self%n%fmesh%yf(1)
  ! allocate 4D storage
  allocate(rwork(in3,1,in1,in2), stat=status)
  call log_alloc_check(m_name,s_name,5,status)
  call misc_readmlab(rwork,in1,in2,in3,0,0,1,1,icbuff,ninra)

  !! assign to field

  if (in3==1) then
     call spl2d_init(self%srcspl2d,rwork(1,1,:,:),in1-1,in2-1,zorg1,zorg2,zh1,zh2,iorder)
  else
     ! volume source : either use spl3d - need to apply F.T. as apb_convert
     !   call spl3d_init(self%vol3,rwork,1,in0,in1,in2,in3,self%xf,self%yf,self%zf,iorder,imin,imax)
     ! or code as look up in 3D volume after beqart
  end if

  deallocate(rwork)

end subroutine radnsource_initfile
!---------------------------------------------------------------------
!> open file
subroutine radnsource_openfile(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer(ki4), intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='radnsource_openfile' !< subroutine name
  logical :: unitused !< flag to test unit is available

  if (trim(file)=='null') then
     call log_error(m_name,s_name,1,log_info,'null filename ignored')
     ninra=-1
     if (present(channel)) channel=ninra
     return
  end if

  !! get file unit
  call misc_getfileunit(ninra)
  if (present(channel)) channel=ninra

  !! open file
  radnsourcefile=trim(file)
  call log_value("Radn source data file",trim(radnsourcefile))
  open(unit=ninra,file=radnsourcefile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open radn sourcefile, ",a)',radnsourcefile
     call log_error(m_name,s_name,2,error_fatal,'Cannot open radn source data file')
     stop
  end if

end subroutine radnsource_openfile
!---------------------------------------------------------------------
!> radiation source from formula
subroutine radnsource_initformul(self,psi,file)

  !! arguments
  type(radnsource_t), intent(inout) :: self !< module object
  type(spl2d_t), intent(inout) :: psi !< \f$ \psi(R,Z) \f$ from input
  character(*), intent(in) :: file !< file name
  !! local
  character(*), parameter :: s_name='radnsource_initformul' !< subroutine name
  type(spl2d_t) :: psin !< normalised \f$ \psi(R,Z) \f$
  integer(ki4) :: in1   !< first nominal dimension
  integer(ki4) :: in2   !< second nominal dimension
  integer(ki4) :: in3   !< third nominal dimension
  integer(ki4)  :: iorder   !< spline order`
  real(kr8) :: zorg1   !< minimum mesh value in coordinate 1
  real(kr8) :: zorg2   !< minimum mesh value in coordinate 2
  real(kr8) :: zh1   !< mesh spacing, coordinate 1
  real(kr8) :: zh2   !< mesh spacing, coordinate 2
  real(kr8) :: zpt1   !< position in coordinate 1
  real(kr8) :: zpt2   !< position in coordinate 2
  real(kr8) :: zpsi   !< flux value
  type(adas_data) :: imp !< contains adas data
  type(profile_data) :: profiles !< contain density and temperature profile data
  iorder=self%n%nraorder

  call spl2d_specoffscale(psi, &
 &self%n%psicen,self%n%psibdry,self%n%rcen,self%n%zcen,self%n%lplasbox,self%n%plasbox, &
 &psin)

  ! fmesh from psi spl2d structure
  self%n%fmesh%nxf=psi%n1p
  self%n%fmesh%nyf=psi%n2p
  self%n%fmesh%nzf=1
  self%n%fmesh%x0f=psi%org1
  self%n%fmesh%y0f=psi%org2
  self%n%fmesh%z0f=0
  self%n%fmesh%dxf=psi%h1
  self%n%fmesh%dyf=psi%h2
  self%n%fmesh%dzf=0
  self%n%fmesh%ndimf=2
  call fmesh_uniform(self%n%fmesh)
  call radnsource_parameter(self)
  ! now as locals
  in1=self%n%fmesh%nxf
  in2=self%n%fmesh%nyf
  in3=self%n%fmesh%nzf
  zorg1=self%n%fmesh%xf(1)
  zorg2=self%n%fmesh%yf(1)
  zh1=self%n%fmesh%xf(2)-self%n%fmesh%xf(1)
  zh2=self%n%fmesh%yf(2)-self%n%fmesh%yf(1)
  ! allocate 4D storage
  allocate(rwork(in3,1,in1,in2), stat=status)
  call log_alloc_check(m_name,s_name,6,status)
  ! loop over mesh
  formula_chosen: select case (self%n%formula) 
  ! reads in adas data if in physics mode
  case('physics')
     call radnsource_adas(imp)
     call radnsource_profile(profiles)
  end select formula_chosen

  do j=1,in2
     zpt2=self%n%fmesh%yf(j)
     do i=1,in1
        zpt1=self%n%fmesh%xf(i)
        call spl2d_eval(psin,zpt1,zpt2,zpsi)
        rwork(1,1,i,j)=radnsource_fn(self,zpsi,zpt1,zpt2,imp,profiles)
     end do
  end do

  !! assign to field

  if (in3==1) then
     call spl2d_init(self%srcspl2d,rwork(1,1,:,:),in1-1,in2-1,zorg1,zorg2,zh1,zh2,iorder)
  else
     ! volume source : either use spl3d - need to apply F.T. as apb_convert
     !   call spl3d_init(self%vol3,rwork,1,in0,in1,in2,in3,self%xf,self%yf,self%zf,iorder,imin,imax)
     ! or code as look up in 3D volumea after beqart
  end if

  deallocate(rwork)

end subroutine radnsource_initformul
!---------------------------------------------------------------------
!> read data from file
subroutine radnsource_readcon(selfn,channel)

  !! arguments
  type(radnumerics_t), intent(out) :: selfn !< type which data will be assigned to
  integer(ki4), intent(in),optional :: channel   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='radnsource_readcon' !< subroutine name
  character(len=80) :: radnsource_formula !< formula to be used
  integer(ki4), parameter :: max_number_of_parameters=10 !< maximum number of parameters allowed
  real(kr8), dimension(max_number_of_parameters) :: general_real_parameters  !< local variable
  integer(ki4), dimension(max_number_of_parameters) :: general_integer_parameters  !< local variable
  integer(ki4) :: number_of_real_parameters  !< local variable
  integer(ki4) :: number_of_integer_parameters  !< local variable
  integer(ki4) :: bremsstrahlung_switch !< switch for Bremsstrahlung
  integer(ki4) :: cyclotron_switch !< switch for cyclotron
  integer(ki4) :: recombination_switch !< switch for recombination
  integer(ki4) :: order_of_source_spline !< order of spline used in radn source
  real(kr8) :: a_exponent !<  \f$ a \f$ parameter for radnsource
  real(kr8) :: te_a_exponent !<  \f$ a_T \f$ parameter for radnsource
  real(kr8) :: ne_a_exponent !<  \f$ a_n \f$ parameter for radnsource
  real(kr8) :: b_exponent !<  \f$ b \f$ parameter for radnsource
  real(kr8) :: te_b_exponent !<  \f$ b_T \f$ parameter for radnsource
  real(kr8) :: ne_b_exponent !<  \f$ b_n \f$ parameter for radnsource
  real(kr8) :: te !< electron temperature
  real(kr8) :: ne !< electron number density
  logical :: plasma_box_set  !< plasma box limits set
  real(kr8) :: box_r_start  !< plasma box limit
  real(kr8) :: box_z_start  !< plasma box limit
  real(kr8) :: box_r_end  !< plasma box limit
  real(kr8) :: box_z_end  !< plasma box limit
  real(kr8) :: z1  !< local variable
  real(kr8) :: z2  !< local variable

  logical :: lte !< set whether lte assumptions
  real(kr8) :: total_power !< total radiative power leaving plasma
  logical :: radiation_diagnostics !< switch for plotting radiation heatmap

  !! radnsource parameters
  namelist /radnsourceparameters/ &
 &radnsource_formula, &
 &a_exponent, b_exponent, &
 &te_a_exponent, te_b_exponent, &
 &ne_a_exponent, ne_b_exponent, &
 &te, ne, lte, &
 &plasma_box_set, &
 &box_r_start, box_z_start, box_r_end, box_z_end, &
 &bremsstrahlung_switch, &
 &cyclotron_switch, &
 &recombination_switch, &
 &order_of_source_spline, &
 &total_power, &
 &radiation_diagnostics, &
 &prof_custom, &
 &frac_impurity, &
 &general_real_parameters, number_of_real_parameters, &
 &general_integer_parameters, number_of_integer_parameters

  !! set default radnsource parameters

  radnsource_formula='unset'
  general_real_parameters=0
  general_integer_parameters=0
  number_of_real_parameters=0
  number_of_integer_parameters=0

  a_exponent = 2
  te_a_exponent = 2
  ne_a_exponent = 2
  b_exponent = 1
  te_b_exponent = 1
  ne_b_exponent = 1
  te = 8800
  ne = 1.415e+20_kr8
  lte = .false.
  plasma_box_set=.false.
  box_r_start=0
  box_z_start=0
  box_r_end=0
  box_z_end=0

  bremsstrahlung_switch=1
  cyclotron_switch=0
  recombination_switch=0
  order_of_source_spline=4
  total_power=-1
  radiation_diagnostics = .false.
  frac_impurity=0

  if(present(channel).and.channel/=0) then
     !! assume unit already open and reading infile
     ninra=channel
  end if

  !!read radnsource parameters
  read(ninra,nml=radnsourceparameters,iostat=status)
  if(status/=0) then
     !!dump namelist contents to logfile to assist error location
     print '("Fatal error reading radnsource parameters")'
     call log_getunit(ilog)
     write(ilog,nml=radnsourceparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading radnsource parameters')
  end if

  call lowor(radnsource_formula,1,len_trim(radnsource_formula))
  !! check for valid data

  formula_chosen: select case (radnsource_formula)
  case('unset','file_input')

  case('circular')
     if (total_power<=0) &
 &   call log_error(m_name,s_name,10,error_fatal,'total_power must be greater than 0')
     if (a_exponent<=0) &
 &   call log_error(m_name,s_name,11,error_fatal,'a_exponent as minor radius must be positive')
     if (b_exponent<const_epsbdry) &
 &   call log_error(m_name,s_name,12,error_fatal,'b_exponent must be positive')

  case('physics')
     if (te_a_exponent<const_epsbdry) &
 &   call log_error(m_name,s_name,20,error_fatal,'te_a_exponent must be positive')
     if (ne_a_exponent<const_epsbdry) &
 &   call log_error(m_name,s_name,21,error_fatal,'ne_a_exponent must be positive')
     if (te_b_exponent<const_epsbdry) &
 &   call log_error(m_name,s_name,22,error_fatal,'te_b_exponent must be positive')
     if (ne_b_exponent<const_epsbdry) &
 &   call log_error(m_name,s_name,23,error_fatal,'ne_b_exponent must be positive')
     if (te<=0) &
 &   call log_error(m_name,s_name,22,error_fatal,'te must be positive')
     if (ne<=0) &
 &   call log_error(m_name,s_name,23,error_fatal,'ne must be positive')

  case('analytic')
     if (total_power<=0) &
 &   call log_error(m_name,s_name,30,error_fatal,'total_power must be greater than 0')
     if (a_exponent<const_epsbdry) &
 &   call log_error(m_name,s_name,31,error_fatal,'a_exponent must be positive')
     if (b_exponent<const_epsbdry) &
 &   call log_error(m_name,s_name,32,error_fatal,'b_exponent must be positive')

  case('userdefined')
     if(number_of_real_parameters<0) &
 &   call log_error(m_name,s_name,44,error_fatal,'number of real parameters must be >=0')
     if(number_of_real_parameters>max_number_of_parameters) then
        call log_value("max number of real parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,45,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters<0) &
 &   call log_error(m_name,s_name,46,error_fatal,'number of integer parameters must be >=0')
     if(number_of_integer_parameters>max_number_of_parameters) then
        call log_value("max number of integer parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,47,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters==0.and.number_of_real_parameters==0) &
 &   call log_error(m_name,s_name,48,error_fatal,'no parameters set')

  case default
     call log_error(m_name,s_name,50,error_fatal,'radnsource_formula not recognised')
  end select formula_chosen

  if(plasma_box_set) then
     ! could apply to any case
     z1=min(box_r_start,box_r_end)
     z2=max(box_r_start,box_r_end)
     if (abs(z1-z2)<const_eps3d) &
 &   call log_error(m_name,s_name,51,error_fatal,'box box too small in R')
     box_r_start=z1
     box_r_end=z2
     z1=min(box_z_start,box_z_end)
     z2=max(box_z_start,box_z_end)
     if (abs(z1-z2)<const_eps3d) &
 &   call log_error(m_name,s_name,51,error_fatal,'box box too small in Z')
     box_z_start=z1
     box_z_end=z2
  end if

  !! store values
  selfn%formula=radnsource_formula

  selfn%apar = a_exponent
  selfn%atpar = te_a_exponent
  selfn%anpar = ne_a_exponent
  selfn%bpar = b_exponent
  selfn%btpar = te_b_exponent
  selfn%bnpar = ne_b_exponent
  selfn%te = te
  selfn%ne = ne
  selfn%lte = lte
  selfn%lplasbox=plasma_box_set
  selfn%plasbox(1)=box_r_start
  selfn%plasbox(2)=box_z_start
  selfn%plasbox(3)=box_r_end
  selfn%plasbox(4)=box_z_end

  selfn%nswbrems=bremsstrahlung_switch
  selfn%nswcyclo=cyclotron_switch
  selfn%nswrecomb=recombination_switch
  selfn%nraorder=order_of_source_spline
  selfn%powtot=total_power
  selfn%radndia=radiation_diagnostics

  selfn%plasbox(1)=box_r_start
  selfn%plasbox(2)=box_z_start
  selfn%plasbox(3)=box_r_end
  selfn%plasbox(4)=box_z_end

  !! allocate arrays and assign

  selfn%nrpams=number_of_real_parameters
  selfn%nipams=number_of_integer_parameters

  formula_allocate: select case (radnsource_formula)

  case('userdefined')
     if (number_of_real_parameters>0) then
        allocate(selfn%rpar(number_of_real_parameters), stat=status)
        call log_alloc_check(m_name,s_name,65,status)
        selfn%rpar=general_real_parameters(:number_of_real_parameters)
     end if
     if (number_of_integer_parameters>0) then
        allocate(selfn%npar(number_of_integer_parameters), stat=status)
        call log_alloc_check(m_name,s_name,66,status)
        selfn%npar=general_integer_parameters(:number_of_integer_parameters)
     end if
  case default
  end select formula_allocate

end  subroutine radnsource_readcon
!---------------------------------------------------------------------
!> set radnsource emission function parameters
subroutine radnsource_parameter(self)

  !! arguments
  type(radnsource_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='radnsource_parameter' !< subroutine name
  real(kr8) :: zint !< value of integral

  zint = exp(const_oned) - 1
  self%prad = 100
  self%normfac2 = 1
  self%normfac = self%prad/(2*const_pid*zint)

  !! do not convert units of parameters for use in radiation formula

  !self%n%bmag = 10000*self%n%bmag ! set bmag to units of gauss from tesla
  !self%n%ne = 1.0e-6_kr8*self%n%ne ! set number Ne to cm^3 from m^3
  !self%n%te = 1000*self%n%te ! set Te to eV from KeV

  !  print*, "bmag=", self%n%bmag
  !  print*, "ne=", self%n%ne
  !  print*, "te=", self%n%te

  !! allocate arrays

  ! species include ions, alphas and any impurities

  allocate(self%noz(self%nspecies))
  allocate(self%ezinf(self%nspecies))
  allocate(self%zspecies(self%nspecies))


end subroutine radnsource_parameter
!---------------------------------------------------------------------
!> normalisation of radiation source
subroutine radnsource_norm(self,psi)

  !! arguments
  type(radnsource_t), intent(inout) :: self !< module object
  type(spl2d_t), intent(inout) :: psi !< \f$ \psi(R,Z) \f$ from input
  !! local
  character(*), parameter :: s_name='radnsource_norm' !< subroutine name
  real(kr8) :: zint !< integrated radiation
  real(kr8) :: scalefac !< scaling factor
  real(kr8) :: zintpow  !< local variable

  call spl2d_sumint(self%srcspl2d,zint,1)

  zintpow=2*const_pid*zint

  if (self%n%powtot>0) then
     scalefac=self%n%powtot/(4*const_pid)/zintpow
  else
     ! source is per steradian
     scalefac=1/(4*const_pid)
  end if
  call spl2d_scale(self%srcspl2d,scalefac)

  if (self%n%radndia) then
     call log_value("total power before normalisation", zintpow)
  end if

end subroutine radnsource_norm
!---------------------------------------------------------------------
!> output to log file
subroutine radnsource_dia(self)

  !! arguments
  type(radnsource_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='radnsource_dia' !< subroutine name


end subroutine radnsource_dia
!---------------------------------------------------------------------
!> general external function call
function radnsource_fn(self,psi,pos1,pos2,imp,prof_data)

  !! arguments
  type(radnsource_t), intent(in) :: self !< module object
  real(kr8) :: radnsource_fn !< function variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$
  real(kr8), intent(in) :: pos1 !< value of \f$ R \f$
  real(kr8), intent(in) :: pos2 !< value of \f$ Z \f$
  type(adas_data) :: imp !< contains adas data
  type(profile_data) :: prof_data !< contains density and temperature data for each impurity

  !! local variables
  character(*), parameter :: s_name='radnsource_fn' !< subroutine name
  real(kr8),dimension(14) :: prepow !< local variable
  real(kr8) :: pow !< local variable

  pow=0._kr8
  !! select radnsource
  formula_chosen: select case (self%n%formula)
  case('unset','userdefined')
     pow=radnsource_userdefined(self,psi,pos1,pos2)
  case('circular')
     pow=radnsource_circular(self,pos1,pos2)
  case('physics')
     pow=radnsource_physics(self,psi,imp,pos1,pos2,prof_data)
  case('analytic')
     prepow=radnsource_analytic(psi,self%n%apar,self%n%bpar,0,pos1,pos2,self)
     pow=prepow(1)
  end select formula_chosen

  !! return radnsource
  radnsource_fn=pow

end function radnsource_fn
!---------------------------------------------------------------------
!> userdefined function
function radnsource_userdefined(self,psi,pos1,pos2)

  !! arguments
  type(radnsource_t), intent(in) :: self !< module object
  real(kr8) :: radnsource_userdefined !< function variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$
  real(kr8), intent(in) :: pos1 !< value of \f$ R \f$
  real(kr8), intent(in) :: pos2 !< value of \f$ Z \f$

  !! local variables
  character(*), parameter :: s_name='radnsource_userdefined' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  integer(ki4) :: ilocal !< local integer variable

  zpos=psi
  pow=0._kr8
  !> user defines \tt{pow} here
  !! .....
  !! return radnsource
  radnsource_userdefined=max(pow,const_zerod)

end function radnsource_userdefined
!---------------------------------------------------------------------
!> circular function
function radnsource_circular(self,pos1,pos2)

  !! arguments
  type(radnsource_t), intent(in) :: self !< module object
  real(kr8) :: radnsource_circular !< function variable
  real(kr8), intent(in) :: pos1 !< value of \f$ R \f$
  real(kr8), intent(in) :: pos2 !< value of \f$ Z \f$

  !! local variables
  character(*), parameter :: s_name='radnsource_circular' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zposq !< position
  real(kr8) :: zrfac  !< local variable

  zposq=(pos1-self%n%rcen)**2+(pos2-self%n%zcen)**2
  zrfac=max( 1-zposq/self%n%apar**2 , const_zerod )
  pow=zrfac**self%n%bpar
  radnsource_circular=pow

end function radnsource_circular
!---------------------------------------------------------------------
!> physics function
function radnsource_physics(self,psi,imp,r,z,prof_data)

  !! arguments
  type(radnsource_t), intent(in) :: self !< module object
  real(kr8) :: radnsource_physics !< function variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$
  type(adas_data) :: imp !< contains adas data
  real(kr8) :: r !< value of \f$ R \f$
  real(kr8) :: z !< value of \f$ Z \f$
  type(profile_data):: prof_data !< contains density and temperature data for each impurity

  !! local variables
  character(*), parameter :: s_name='radnsource_physics' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  real(kr8),dimension(14) :: zne !<\f$ n_e \f$, local number density of electrons in plasma
  real(kr8),dimension(14) :: zte !<\f$ T_e \f$, local temperature of electrons (in eV)
  zpos=psi
  if(prof_custom /= 2) then
     zne=self%n%ne*radnsource_analytic(zpos,self%n%anpar,self%n%bnpar,0,r,z,self,prof_data)
     zte=self%n%te*radnsource_analytic(zpos,self%n%atpar,self%n%btpar,1,r,z,self,prof_data)
  else
     zne=radnsource_analytic(zpos,self%n%anpar,self%n%bnpar,0,r,z,self,prof_data)
     zte=radnsource_analytic(zpos,self%n%atpar,self%n%btpar,1,r,z,self,prof_data)
  end if
  pow=self%n%nswbrems*radnsource_linebrem(self,zne,zte,imp)!+&
  !  self%n%nswrecomb*radnsource_recomb(self,zne,zte)
  radnsource_physics=pow

end function radnsource_physics
!---------------------------------------------------------------------
!> adas data read in
subroutine radnsource_adas(impurities)
  !!arguments
  type(adas_data) :: impurities !< contains adas data
  !local
  integer (ki4) :: j !< loop counter
  integer :: ninz !< local unit number
  impurities%n_impurities=14
  allocate(impurities%imp_pre(impurities%n_impurities))
  !defines prefix for file names involving impurities
  impurities%imp_pre(1)='ar';   impurities%imp_pre(2)='c'; impurities%imp_pre(3)='h'
  impurities%imp_pre(4)='kr';   impurities%imp_pre(5)='ne';   impurities%imp_pre(6)='o'
  impurities%imp_pre(7)='w'; impurities%imp_pre(8)='be';   impurities%imp_pre(9)='fe'
  impurities%imp_pre(10)='he';impurities%imp_pre(11)='n';  impurities%imp_pre(12)='ni'
  impurities%imp_pre(13)='si';impurities%imp_pre(14)='xe'

  !reads adas data for all the impurities labelled in myfilename
  do j=1,impurities%n_impurities
     call misc_getfileunit(ninz)
     open (ninz, file = trim(impurities%imp_pre(j))//'_lz_tau.dat' , status = 'old')
     do i=1,12
        if(i==4) then
           read(ninz,*) impurities%n_temps , impurities%dwell_time
        else
           read(ninz,*)
        end if
     enddo
     if(j==1) then
        allocate(impurities%loss_fun(impurities%n_impurities,impurities%n_temps))
        allocate(impurities%temps(impurities%n_temps))
        impurities%loss_fun=0
        impurities%temps=0
     end if
     do i=0,impurities%n_temps/8 - 1
        Read(ninz,'(X,7(E9.3,X),E9.3)') impurities%temps(1 + 8*i : 8 + 8*i )
     end do
     do i=1,impurities%dwell_time + impurities%dwell_time*(impurities%n_temps/8) + 1
        read(ninz,*)
     enddo
     do i=0,impurities%n_temps/8 - 1
        Read(ninz,'(X,7(E9.3,X),E9.3)') impurities%loss_fun(j,1 + 8 *i : 8 + 8*i )
     end do
     close(ninz)
  enddo
end subroutine
!---------------------------------------------------------------------
!> temperature and density profile read in
subroutine radnsource_profile(profiles)
  !!arguments
  type(profile_data) :: profiles !< contains density and temperature data for each impurity
  integer :: ninz !< local unit number
  integer(ki4) :: j !< loop counter
  if(prof_custom==2) then
     call misc_getfileunit(ninz)
     open (ninz, file = 'profiles.txt' , status = 'old')
     do j=1,10001
        read(ninz,*) profiles%psi(j), profiles%density_profile(j,:),profiles%temperature_profile(j,:)
     end do
     close(ninz)
  end if
end subroutine
!---------------------------------------------------------------------
!> line + bremstrahlung radiation produced by plasma
function radnsource_linebrem(self,pne,pte,imp)

  !! arguments
  type(radnsource_t), intent(in) :: self !< module object
  real(kr8),dimension(14), intent(in) :: pne !< \f$ n_e \f$, local number density of electrons in plasma
  real(kr8),dimension(14), intent(in) :: pte !< \f$ t_e \f$, local temperature of electrons (in ev)
  real(kr8) :: radnsource_linebrem !< function variable
  type(adas_data) :: imp !< impurity data from adas data

  !! local
  character(*), parameter :: s_name='radnsource_linebrem' !< subroutine name
  integer (ki4) :: j !< local variable
  integer (ki4) :: id(1) !< local variable
  !  real (kr8), dimension(:), allocatable :: frac_imp
  real (kr8) :: plb !< local variable
  real (kr8) :: dldt !< local variable
  real (kr8) :: f !< local variable

  !  allocate(frac_imp(imp%n_impurities))
  plb=0
  if(prof_custom==2) frac_impurity=1
  do j=1,imp%n_impurities

     id=minloc(abs(pte(j)-imp%temps))
     !  if(id(1)==0) then
     !     dldt=(imp%loss_fun(j,1)-imp%loss_fun(j,0))/(imp%temps(1)-imp%temps(0))
     !     f=imp%loss_fun(j,0) + (dldt)*(pte(j)-imp%loss_fun(j,0))
     !     plb=plb+frac_impurity(j)*f*(pne(j)**2)
     !  else if(id(1)==imp%n_temps) then
     !     dldt=(imp%loss_fun(j,id(1))-imp%loss_fun(j,id-1))/(imp%temps(id)-imp%temps(id-1))
     !     f=imp%loss_fun(j,id(1)-1) + (dldt)*(pte(j)-imp%loss_fun(j,id-1))
     !     plb=plb+frac_impurity(j)*f*(pne(j)**2)
     !  else
     !     if(imp%temps(id(1)) .lt. pte(j) ) then
     !        dldt=(imp%loss_fun(j,id(1)+1)-imp%loss_fun(j,id))/(imp%temps(id+1)-imp%temps(id))
     !        f=imp%loss_fun(j,id(1)) + (dldt)*(pte(j)-imp%loss_fun(j,id))
     !     else
     !        dldt=(imp%loss_fun(j,id(1))-imp%loss_fun(j,id-1))/(imp%temps(id)-imp%temps(id-1))
     !        f=imp%loss_fun(j,id(1)-1) + (dldt)*(pte(j)-imp%loss_fun(j,id))
     !     end if
     !     plb=plb+frac_impurity(j)*f*(pne(j)**2)
     !  end if
     plb=plb+frac_impurity(j)*imp%loss_fun(j,id(1))*(pne(j)**2)
  enddo
  !  deallocate(frac_imp,imp%temps,imp%loss_fun)
  radnsource_linebrem=plb
end function radnsource_linebrem
!---------------------------------------------------------------------
!!> cyclotron radiation produced by plasma (version 2)
!function radnsource_cyclo2(self,pne,pte)
!
!  !! arguments
!  type(radnsource_t), intent(in) :: self !< module object
!  real(kr8), intent(in) :: pne !< \f$ n_e \f$, local number density of electrons in plasma
!  real(kr8), intent(in) :: pte !< \f$ t_e \f$, local temperature of electrons (in ev)
!  real(kr8) :: radnsource_cyclo2 !< function variable
!  !! local
!  character(*), parameter :: s_name='radnsource_cyclo2' !< subroutine name
!  real(kr8) :: pc !< cyclotron radiation power density (watt/m^3)
!  real(kr8) :: synch_const=9.645e-10,t0,n0,p_a0,s_m_t,k_x,x
!  t0=self%n%te/1000.0
!  n0=self%n%ne/(10.0d0**20.0d0)
!  p_a0=(6.04e+3)*(s_m_t*n0/self%n%bmag)
!! calculate cyclotron radiation
!  zne=self%n%ne*radnsource_analytic(zpos,self%n%anpar,self%n%bnpar)
!  zte=self%n%te*radnsource_analytic(zpos,self%n%atpar,self%n%btpar)
!     x=sqrt(self%n%zcen**2 + self%n%rcen**2)/s_m_t
!     pc = (10.0d0**6.0d0)*synch_const*sqrt(self%n%zcen**2 + self%n%rcen**2)&
!     *(self%n%bmag**2.62d0)*(s_minor_t**1.38d0)*(k_x**0.79d0)*&
!      (n0**0.38d0)*t0*((16.0d0 + t0)**2.61d0)*&
!      k(self%n%anpar,self%n%atpart,self%n%btpar)*&
!        g(x)/((1.0d0 + 0.12d0*t0/p_a0)**1.51d0)
!
! radnsource_cyclo2 = pc
!
!end function radnsource_cyclo2
!aspect ratio correction factor
!function g(x)
!  real (kr8) :: g,x
!  g=0.93d0*(1.0d0 + 0.85d0*exp(-0.82d0*x))
!
!end function
!profile correction factor
!function k(alpha_n_syn,alpha_t_syn,beta_t_syn)
!  real (kr8) :: k,alpha_n_syn,alpha_t_syn,beta_t_syn
!  k = ((alpha_n_syn + 3.87d0*alpha_t_syn + 1.46d0)**-0.79d0)*&
!     ((1.98d0 + alpha_t_syn)**1.36d0)*&
!     ((beta_t_syn**1.53d0 +1.87d0*alpha_t_syn -0.16d0)**-1.33d0)*&
!     (beta_t_syn**2.14d0)
!end function
!---------------------------------------------------------------------
!> cyclotron radiation produced by plasma
function radnsource_cyclo(self,pne,pte)

  !! arguments
  type(radnsource_t), intent(in) :: self !< module object
  real(kr8), intent(in) :: pne !< \f$ n_e \f$, local number density of electrons in plasma
  real(kr8), intent(in) :: pte !< \f$ T_e \f$, local temperature of electrons (in eV)
  real(kr8) :: radnsource_cyclo !< function variable
  !! local
  character(*), parameter :: s_name='radnsource_cyclo' !< subroutine name
  real(kr8) :: pc !< cyclotron radiation power density (watt/m^3)

  !! calculate cyclotron radiation

  if (self%n%lte) then
     ! if lte specified then calculate cyclotron rad from ne and te
     pc = 5.00e-44_kr8*(pne**2)*(pte**2)
  else
     ! otherwise calculate using magnitude of magnetic field B
     pc = 6.21e-20_kr8*(self%n%bmag**2)*pne*pte
     ! should the constant be 0.387554
     ! current constant comes out as (e^5)/(3 pi e_0 m_e^3 c^3)
     ! this is due to conversion of temperature to eV
     !print*, "no lte"
  end if

  !! return power density
  radnsource_cyclo = pc

end function radnsource_cyclo
!---------------------------------------------------------------------
!> recombination (free-bound) radiation
function radnsource_recomb(self,pne,pte,imp)

  !! arguments
  type(radnsource_t), intent(in) :: self !< module object
  real(kr8), intent(in) :: pne !< \f$ n_e \f$, local number density of electrons in plasma
  real(kr8), intent(in) :: pte !< \f$ T_e \f$, local temperature of electrons (in eV)
  real(kr8) :: radnsource_recomb !< function variable
  type(adas_data) :: imp !< contains adas data

  !! local
  character(*), parameter :: s_name='radnsource_recomb' !< subroutine name
  real(kr8) :: pr !< recombination radiation power density (watt/m^3)
  real(kr8) :: prefac !< factor before sum
  real(kr8) :: sum_z !< sum over ionization states, Z
  real(kr8) :: e_ion !< energy level within ion
  real (kr8), dimension(:), allocatable :: frac_imp !< fraction of each impurity
  real (kr8), dimension(:), allocatable :: frac_ion !< fraction of each impurity
  integer(ki4) :: n !< principal quantum number
  integer(ki4), dimension(:), allocatable :: z !< charge of impurity
  allocate(frac_imp(imp%n_impurities),frac_ion(imp%n_impurities),z(imp%n_impurities))
  frac_imp=const_eps2d
  radnsource_recomb=0
  do n=1,imp%n_impurities
     e_ion = const_rydberg*(z(n)**2)/n**2
     radnsource_recomb = radnsource_recomb + (4*sqrt(2*const_pid)/sqrt(3.0_kr8))&
    &*((const_charge**2)/(4*const_pid*const_epsilon0))*&
    &(4/(3*(const_clight**3)))*(1/(const_masse*sqrt(const_masse)))*&
    &(3*const_rydberg*(z(n)**4)/n**3)*&
    &(pne**2)*frac_imp(n)*frac_ion(n)*(1/pte)*exp(e_ion/pte)
  end do

end function radnsource_recomb
!---------------------------------------------------------------------
!> bremsstrahlung radiation
function radnsource_brems(self,pne,pte)

  !! arguments
  type(radnsource_t), intent(in) :: self !< module object
  real(kr8), intent(in) :: pne !< \f$ n_e \f$, local number density of electrons in plasma
  real(kr8), intent(in) :: pte !< \f$ T_e \f$, local temperature of electrons (in eV)
  real(kr8) :: radnsource_brems !< function variable
  !! local
  character(*), parameter :: s_name='radnsource_brems' !< subroutine name
  real(kr8) :: pbr !< bremsstrahlung radiation power density (watt/m^3)
  real(kr8) :: prefac !< factor before sum
  real(kr8) :: sum_z !< sum over ionization states, Z
  integer(ki4) :: z !< ionization state Z

  !! bremsstrahlung radiation
  ! As a first approximation we can assume that a DT plasma radiation = n^2 * t^1/2

  prefac = 1.69e-38_kr8*(sqrt(pte))*pne**2

  !do z = 1, self%self%nspecies ! sum over all species in plasma
  ! sum_z = (self%zpecies(z)**2)*self%noz(z)
  !end do
  !+ sum_z calculate total power density
  pbr = prefac

  !! return power density
  radnsource_brems = pbr

end function radnsource_brems
!---------------------------------------------------------------------
!> analytic function
function radnsource_analytic(psi,apar,bpar,choice,r2,z2,self,prof)

  !! arguments
  type(radnsource_t), intent(in) :: self !< module object
  real(kr8),dimension(14) :: radnsource_analytic !< function variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$
  real(kr8), intent(in) :: apar !< value of \f$ a \f$
  real(kr8), intent(in) :: bpar !< value of \f$ b \f$
  integer (ki4), intent(in) :: choice !< used to pick density or temperature allocation within subroutine
  real(kr8) :: r2 !< value of \f$ R \f$
  real(kr8) :: z2 !< value of \f$ Z \f$
  type(profile_data), optional :: prof !<contains density and temperature data for each impurity

  !! local variables
  character(*), parameter :: s_name='radnsource_analytic' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zfac  !< local variable
  integer (ki4) :: nlines !< local variable
  integer (ki4) :: io !< local variable
  integer (ki4) :: i !< local variable
  integer :: ninz !< local unit number
  real(kr8) :: dummy !< to be moved
  real(kr8) :: zte !< to be moved
  real(kr8) :: zti !< to be moved
  real(kr8) :: zne !< to be moved
  real(kr8) :: zni !< to be moved
  real(kr8) :: zr !< to be moved
  real(kr8) :: distance !< distance from centre of plasma
  integer (ki4) :: id(1) !< local variable
  distance=sqrt((r2-self%n%rcen)**2 + (z2-self%n%zcen)**2)

  if(prof_custom==0) then
     zfac=max(1-abs(psi)**apar,const_zerod)
     !default function for profiles
     pow=zfac**bpar
     radnsource_analytic=pow
  end if

  if(prof_custom==1) then
     !user input function for profiles

     if(choice==0) then
        !default function for profiles
        zfac=max(1-abs(distance)**apar,const_zerod)
        pow=zfac**bpar
     end if

     if(choice==1) then
        !default function for profiles
        zfac=max(1-abs(distance)**apar,const_zerod)
        pow=zfac**bpar
     end if
     radnsource_analytic=pow
  end if

  if(prof_custom==2) then
     !fortran generated profiles

     id=minloc(abs(psi-prof%psi))
     do i=1,14
        if(choice==0) then
           radnsource_analytic(i)=prof%density_profile(id(1),i)
        end if

        if(choice==1) then
           radnsource_analytic(i)=prof%temperature_profile(id(1),i)
        end if
     end do


  end if

  if(prof_custom==3) then
     !user read in function for profiles (transport code)
     call misc_getfileunit(ninz)
     open (unit = ninz, file = "custom_psi.txt")
     read(ninz,*)
     read(ninz,*)
     nlines=0
     do
        read(ninz,9002,iostat=io)
        if (io/=0) exit
        nlines = nlines + 1
     end do
     rewind(ninz)
     read(ninz,*)
     read(ninz,*)
     do i=1,nlines
        read(ninz,9002) dummy,zr,zte,dummy,zti,dummy,zne,zni,dummy,dummy,dummy,dummy
     end do
9002     format(4x,i4,3x,f7.3,1x,2(3x,1pe10.3,1x,a2),6(3x,1pe10.3))
     close(ninz)
  end if


end function radnsource_analytic
!---------------------------------------------------------------------
!> open new file, making up name
subroutine radnsource_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer(ki4), intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='radnsource_initwrite' !< subroutine name
  logical :: unitused !< flag to test unit is available
  character(len=80) :: outputfile !< output file name

  !! get file unit
  call misc_getfileunit(noutra)
  if (present(channel)) channel=noutra

  !! open file
  outputfile=trim(fileroot)//"_radnsource.out"
  call log_value("Control data file",trim(outputfile))
  open(unit=noutra,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutra,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if

end subroutine radnsource_initwrite
!---------------------------------------------------------------------
!> write radnsource data
subroutine radnsource_write(self,channel)

  !! arguments
  type(radnsource_t), intent(in) :: self   !< radnsource data structure
  integer(ki4), intent(in), optional :: channel   !< output channel for radnsource data structure

  !! local
  character(*), parameter :: s_name='radnsource_write' !< subroutine name
  integer(ki4) :: iout   !< output channel for radnsource data structure

  !! sort out unit
  if(present(channel)) then
     iout=channel
  else
     iout=noutra
  end if

  write(iout,*,iostat=status) 'radnsource_formula'
  call log_write_check(m_name,s_name,18,status)
  write(iout,*,iostat=status) self%n%formula
  call log_write_check(m_name,s_name,19,status)
  write(iout,*,iostat=status) 'nrpams'
  call log_write_check(m_name,s_name,46,status)
  write(iout,*,iostat=status) self%n%nrpams
  call log_write_check(m_name,s_name,47,status)
  if (self%n%nrpams>0) then
     write(iout,*,iostat=status) 'real_parameters'
     call log_write_check(m_name,s_name,48,status)
     write(iout,*,iostat=status) self%n%rpar
     call log_write_check(m_name,s_name,49,status)
  end if
  write(iout,*,iostat=status) 'nipams'
  call log_write_check(m_name,s_name,50,status)
  write(iout,*,iostat=status) self%n%nipams
  call log_write_check(m_name,s_name,51,status)
  if (self%n%nipams>0) then
     write(iout,*,iostat=status) 'integer_parameters'
     call log_write_check(m_name,s_name,52,status)
     write(iout,*,iostat=status) self%n%npar
     call log_write_check(m_name,s_name,53,status)
  end if

end subroutine radnsource_write
!---------------------------------------------------------------------
!> write object data as gnuplot
subroutine radnsource_writeg(self,select,channel)

  !! arguments
  type(radnsource_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for radnsource data structure

  !! local
  character(*), parameter :: s_name='radnsource_writeg' !< subroutine name
  integer(ki4) :: iout   !< output channel for radnsource data structure

  call log_error(m_name,s_name,1,log_info,'gnuplot file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine radnsource_writeg
!---------------------------------------------------------------------
!> write object data as vtk
subroutine radnsource_writev(self,select,channel)

  !! arguments
  type(radnsource_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for radnsource data structure

  !! local
  character(*), parameter :: s_name='radnsource_writev' !< subroutine name
  integer(ki4) :: iout   !< output channel for radnsource data structure

  call log_error(m_name,s_name,1,log_info,'vtk file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine radnsource_writev
!---------------------------------------------------------------------
!> close write file
subroutine radnsource_closewrite

  !! local
  character(*), parameter :: s_name='radnsource_closewrite' !< subroutine name

  !! close file
  close(unit=noutra,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine radnsource_closewrite
!---------------------------------------------------------------------
!> delete object
subroutine radnsource_delete(self)

  !! arguments
  type(radnsource_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='radnsource_delete' !< subroutine name

  formula_deallocate: select case (self%n%formula)
  case('userdefined')
     if (self%n%nrpams>0) deallocate(self%n%rpar)
     if (self%n%nipams>0) deallocate(self%n%npar)
  case default
  end select formula_deallocate


end subroutine radnsource_delete
!---------------------------------------------------------------------
!> close file
subroutine radnsource_close

  !! local
  character(*), parameter :: s_name='radnsource_close' !< subroutine name

  !! close file
  close(unit=ninra,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close radn source file, ",a)',radnsourcefile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close radn source data file')
     stop
  end if

end subroutine radnsource_close

end module radnsource_m
