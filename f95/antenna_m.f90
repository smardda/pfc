module antenna_m

  use misc_m
  use control_h
  use radnsource_h
  use radnsource_m
  use antenna_h
  use log_m
  use const_numphys_h
  use const_kind_m
  use posang_h
  use posang_m
  use spl2d_m
  use beq_h
  use beq_m
  use vfile_m
  use geobj_m
  use geobjlist_h
  use geobjlist_m
  use powelt_m
  use btree_m
  use ray_m

  implicit none
  private

! public subroutines
  public :: &
  antenna_init, &  !< initialise object
  antenna_initfile,  & !< open file
  antenna_readcon,  & !< read data from file
  antenna_normal, & !< calculate and store launch positions and directions of rays
  antenna_centre, & !< translate paraview generated antenna's centre to (0,0,0)
  antenna_align, & !< align antenna to CAD triangle
  antenna_walltest, & !< testing antenna functionality to produce antennas for each CAD triangle
  antenna_testelt, & !< testing antenna functionality to produce antennas for each CAD triangle
  antenna_powflux, & !< flux function that assigns power to ray segments
  antenna_userdefined,  & !< user-defined function
  antenna_fn, &  !< general external function call
  antenna_elt, & !< launch rays from antenna
  antenna_powelt, & !< launch rays from antenna/powelt
  antenna_dia, &  !< object diagnostics to log file
  antenna_initwrite, & !< open new file, making up name
  antenna_write, &  !< write out object
  antenna_writeg, &  !< write out object as gnuplot
  antenna_writev, &  !< write out object as vtk
  antenna_delete, & !< delete object
  antenna_close, & !< close file
  antenna_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='antenna_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4), save  :: ninant=5     !< control file unit number
  integer(ki4), save  :: noutant=6      !< output file unit number
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4)  :: ilog      !< for namelist dump after error
  integer(ki4) :: testunit !< general file unit for gnuplot quick testing puporses (value changes whenever misc_getfileunit called)

  character(len=29), parameter :: cfmtbv0c='(g14.7,", ",g14.7,", ",g14.7)'
  character(len=39), parameter :: cfmtbv1c='(g16.7,"",g16.7,"",g16.7)'
  real(kr8), parameter :: const_eps5=1.0e-5_kr4 !< single precision small
  real(kr8), parameter :: const_halfd=0.5_kr8 !< double precision half
  real(kr8), parameter :: const_twod=2.0_kr8 !< double precision two
  real(kr8), parameter :: const_zerod=0.0_kr8 !< double precision zero
  real(kr8), parameter :: const_oned=1.0_kr8 !< double precision unit
  real(kr8), parameter :: const_epsd=0.1_kr8 !< double precision tenth
  real(kr8), parameter :: const_eps2d=0.01_kr8 !< double precision hundredth
  real(kr8), parameter :: const_eps3d=0.001_kr8 !< double precision thousandth
  real(kr8), parameter :: const_tend=10.0_kr8 !< double precision ten
  real(kr8), parameter :: const_ten8d=1.0e+8_kr8 !<  small infinity

  contains
!---------------------------------------------------------------------
!> initialise object
subroutine antenna_init(self)

  !! arguments
  type(antenna_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='antenna_init' !< subroutine name

  call geobjlist_area(self%hemisphere,self%area)

  !! allocate arrays (needed for diagnostic)
  allocate(self%dir(self%hemisphere%ng), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  allocate(self%posdir(self%hemisphere%ng), stat=status)
  call log_alloc_check(m_name,s_name,3,status)

end subroutine antenna_init
!---------------------------------------------------------------------
!> open file
subroutine antenna_initfile(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer(ki4), intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='antenna_initfile' !< subroutine name
  logical :: unitused !< flag to test unit is available

  if (trim(file)=='null') then
     call log_error(m_name,s_name,1,log_info,'null filename ignored')
     ninant=-1
     if (present(channel)) channel=ninant
     return
  end if

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  call misc_getfileunit(ninant)
  open(unit=ninant,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,2,error_fatal,'Cannot open control data file')
     stop
  end if
  if (present(channel)) channel=ninant

end subroutine antenna_initfile
!---------------------------------------------------------------------
!> read data from file
subroutine antenna_readcon(selfn,channel)

  !! arguments
  type(antnumerics_t), intent(out) :: selfn !< type which data will be assigned to
  integer(ki4), intent(in),optional :: channel   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='antenna_readcon' !< subroutine name
  character(len=80) :: antenna_formula !< formula to be used
  character(len=80) :: pattern !< sampling pattern to be used
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  real(kr8) :: power_split !< variable with meaningful name

  real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters  !< local variable
  integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters  !< local variable
  integer(ki4) :: number_of_real_parameters  !< local variable
  integer(ki4) :: number_of_integer_parameters  !< local variable
  integer(ki4) :: ntheta !< number of samples in theta, normally overwritten by value in geometry file
  integer(ki4) :: nphi !< number of samples in phi, normally overwritten by value in geometry file
  logical :: plot_diagvtk !< logical to decide whether to visualise or not

  !! antenna parameters
  namelist /antennaparameters/ &
 &power_split, &
 &ntheta, nphi, &
 &general_real_parameters, number_of_real_parameters, &
 &general_integer_parameters, number_of_integer_parameters, &
 &plot_diagvtk

  !! set default antenna parameters
  power_split=0.5_kr8

  antenna_formula='unset'
  general_real_parameters=0
  general_integer_parameters=0
  number_of_real_parameters=0
  number_of_integer_parameters=0
  ntheta = 25
  nphi = 25
  plot_diagvtk = .false.

  if(present(channel).AND.channel/=0) then
     !! assume unit already open and reading infile
     ninant=channel
  end if

  !!read antenna parameters
  read(ninant,nml=antennaparameters,iostat=status)
  if(status/=0) then
     !!dump namelist contents to logfile to assist error location
     print '("Fatal error reading antenna parameters")'
     call log_getunit(ilog)
     write(ilog,nml=antennaparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading antenna parameters')
  end if

  call lowor(antenna_formula,1,len_trim(antenna_formula))
  !! check for valid data

  formula_chosen: select case (antenna_formula)
  case('unset','exp')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,11,error_fatal,'power_split must be >=0 and <=1')

  case('expdouble')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,21,error_fatal,'power_split must be >=0 and <=1')

  case('userdefined')
     if(number_of_real_parameters<0) &
 &   call log_error(m_name,s_name,44,error_fatal,'number of real parameters must be >=0')
     if(number_of_real_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of real parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,45,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters<0) &
 &   call log_error(m_name,s_name,46,error_fatal,'number of integer parameters must be >=0')
     if(number_of_integer_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of integer parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,47,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters==0.AND.number_of_real_parameters==0) &
 &   call log_error(m_name,s_name,48,error_fatal,'no parameters set')

  end select formula_chosen

  !! store values
  selfn%formula=antenna_formula

  selfn%f=power_split

  selfn%ntheta = ntheta
  selfn%nphi = nphi
  selfn%plot_diagvtk = plot_diagvtk

  !! allocate arrays and assign

  selfn%nrpams=number_of_real_parameters
  selfn%nipams=number_of_integer_parameters


  formula_allocate: select case (antenna_formula)

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

end  subroutine antenna_readcon
!---------------------------------------------------------------------
!> triagulate hemispherical surface
subroutine antenna_normal(self,beq,ray)

  !! arguments
  type(antenna_t), intent(inout) :: self !< module object
  type(beq_t), intent(inout) :: beq !< local variable
  type(ray_t), intent(inout) :: ray !< local variable
  !! local
  character(*), parameter :: s_name='antenna_normal' !< subroutine name
  real(kr4) :: area !< surface area of hemisphere
  type(geobj_t) :: igeobj !< local variable
  real(kr4), dimension(3) :: normal !< normal of triangle on hemisphere
  real(kr4) :: mag !< magnitude of normal
  integer(ki4) :: jtri !< loop counter
  integer :: kplot !< local variable
  integer :: unorm !< local variable
  integer :: ucentr !< local variable
  character(10) :: kchar !< local variable
  real(kr4), dimension(3) :: centre !< local variable


  !call misc_getfileunit(unorm)
  !open(unit = unorm, file = "trinorm.txt")

  !call misc_getfileunit(ucentr)
  !open(unit = ucentr, file = "tricentre.txt")


  ! loop over every triangle in antenna
  do jtri = 1, self%hemisphere%ng
     igeobj%geobj=self%hemisphere%obj2(jtri)%ptr
     igeobj%objtyp=self%hemisphere%obj2(jtri)%typ
     ! calculate tri centre
     call geobj_centre(igeobj,self%hemisphere%posl,self%hemisphere%nodl,centre) 
     ! store ray launch position as tri centre
     self%posdir(jtri)%posvec = centre 
     !write(ucentr,cfmtbv0c) self%posdir(jtri)%posvec
     ! calculate tri normal
     call geobj_normal(igeobj,self%hemisphere%posl,self%hemisphere%nodl,normal,mag) 
     self%dir(jtri)%posvec = normal 
     ! store ray launch direction as tri normal
     !write(unorm,cfmtbv0c)self%dir(jtri)%posvec
     !call antenna_launch(self,hemisphere,beq,self%posdir(jtri),self%dir(jtri))
     call ray_launch(ray,beq,self%posdir(jtri),self%dir(jtri))

  end do

  !! write out antenna as vtk for testing purposes (prompts error if ran with large mesh)
  !> provides visualisation capabilities to see if/how the hemispheres align to wall triangles

  if (self%n%plot_diagvtk) then
     kchar = 'geometry'
     call misc_getfileunit(kplot)
     call geobjlist_writev(self%hemisphere,kchar,kplot)
  end if

end subroutine antenna_normal
!---------------------------------------------------------------------
!> accumulate power along a ray
!subroutine antenna_poweray(ray,raypowe,beq)

!! arguments
!  type(poweray_t), intent(inout) :: ray !< module object
!  real(kr4), intent(out) :: raypowe
!  type(beq_t), intent(in) :: beq

!! local
!  character(*), parameter :: s_name='antenna_poweray' !< subroutine name
!  integer(ki4) :: j !< loop counter
!  real(kr4) :: intpow !< interim power
!  real(kr4) :: powtotal !< total power of ray
!  real(kr4) :: psix

!  psix = beq%n%psiref


!  powtotal = 0
!  do j = 1,ray%nsample

!    intpow = antenna_powflux(ray%fluxposl%pos(j)%posvec(1),psix)

!    powtotal = powtotal + intpow

!  end do

!  raypowe = powtotal

!end subroutine antenna_poweray
!---------------------------------------------------------------------

!> test whether hemisphere is centered on (0,0,0) stop code if even value of phi used
subroutine antenna_centre(hemisphere)

  !! arguments
  type(geobjlist_t), intent(inout) :: hemisphere !< trinagulated hemisphere representing antenna

  !! local
  character(*), parameter :: s_name='antenna_centre' !< subroutine name
  integer(ki4) :: jpos !< loop counter over postions in hemisphere geobjlist
  integer(ki4) :: transform  !< local variable
  type(posvecl_t) :: offset !< offset vector to transform positions
  real(kr4), dimension(3) :: scale !< local variable
  logical :: adjust = .false. !< decide whether or not to carry out adjustment
  integer :: kplot !< local variable
  character(10) :: kchar !< local variable


  !! search for if offset needed and how much by
  do j = 1, hemisphere%np
     if (hemisphere%posl%pos(j)%posvec(3)<0) then
        !offset%posvec(1) = 0
        !offset%posvec(2) = 0
        !offset%posvec(3) = - hemisphere%posl%pos(j)%posvec(3)
        adjust = .true.
        exit
     end if
  end do

  !  !! stop code and prompt user to rerun with odd number valued phi resolution for hemisphere

  if (adjust) then
     print '("Fatal error: Even valued phi resolution of hemisphere selected. Rerun with odd valued phi")'
     call log_error(m_name,s_name,1,error_fatal,'Even valued phi resolution selected')
     stop
  end if

  !! translate hemisphere base up to z=0 xy plane
  ! currently not needed as code is setup to stop if even value of phi for hemisphere selected

  !if (adjust) then
  !  transform = 1
  !  scale = (/1,1,1/)
  !  hemisphere%tfmdata%ntfm = 2 ! set transform 2 from position_tfm
  !  hemisphere%tfmdata%scale = scale
  !  hemisphere%tfmdata%offset = offset%posvec
  !  call geobjlist_tfm(hemisphere,transform)
  !  call log_error(m_name,s_name,70,log_info,'antenna hemisphere centered on origin')
  !end if

  !kchar = 'geometry'
  ! call misc_getfileunit(kplot)
  ! call geobjlist_writev(hemisphere,kchar,kplot)

end subroutine antenna_centre
!---------------------------------------------------------------------
!> align antenna to CAD triangle
subroutine antenna_align(self,tricentre,trinormal,beq,ray)

  !! arguments
  type(antenna_t), intent(inout) :: self !< module object
  type(posvecl_t), intent(inout) :: tricentre !< centre of triangle to position hemisphere on
  type(posvecl_t), intent(inout) :: trinormal !< normal of triangle to rotationally align to
  type(beq_t), intent(inout) :: beq !< local variable
  type(ray_t), intent(inout) :: ray !< local variable
  !! local
  character(*), parameter :: s_name='antenna_align' !< subroutine name
  real(kr4), dimension(3,3) :: tfm !< transformation matrix aligns polar axis to triangle normal
  real(kr4), dimension(3,3) :: trtfm !< matrix transpose
  real(kr4) :: n !< transformation matrix component
  real(kr4) :: n1 !< transformation matrix component
  real(kr4) :: n2 !< transformation matrix component
  real(kr4) :: n3 !< transformation matrix component
  real(kr4) :: nf !< normalisation factor if triangle normal not given as unit normal
  integer(ki4) :: jin !< loop counter
  integer(ki4) :: transform !< specify whether forwards or backwards transform
  real(kr4), dimension(3) :: scale !< local variable

  !! setup transformation matrix to align antenna to given triangle

  n1 = trinormal%posvec(1)
  n2 = trinormal%posvec(2)
  n3 = trinormal%posvec(3)
  n = sqrt(n1**2 + n2**2)

  if (n3<=0) then
     n = -n ! deal with n3 < 0 case specially
  end if

  if (abs(n)<const_eps5*abs(n3)) then
     call log_error(m_name,s_name,1,log_info,'Problematic alignment case')
     call log_value("parameter |n| ",abs(n))
     tfm(1,1) = n3
     tfm(1,2) = 0
     tfm(1,3) = 0
     tfm(2,1) = 0
     tfm(2,2) = n3
     tfm(2,3) = 0
     tfm(3,1) = 0
     tfm(3,2) = 0
     tfm(3,3) = n3
  else
     tfm(1,1) = n1*n3/n
     tfm(1,2) = -n2/n
     tfm(1,3) = n1
     tfm(2,1) = n2*n3/n
     tfm(2,2) = n1/n
     tfm(2,3) = n2
     tfm(3,1) = -n
     tfm(3,2) = 0
     tfm(3,3) = n3
  end if

  trtfm = transpose(tfm) ! calculate transpose of rotation matrix (to be used in inverse transform)

  scale = (/1,1,1/) ! set scale to
  self%hemisphere%tfmdata%ntfm = 2 ! set transform 2 from position_tfm
  self%hemisphere%tfmdata%scale = scale
  transform = 1 ! specify that this is forward transform
  self%hemisphere%tfmdata%offset = tricentre%posvec ! set offset vector for transform

  ! perform rotation of positions via rotation matrix tfm
  do jin = 1, self%hemisphere%np
     self%hemisphere%posl%pos(jin)%posvec = matmul(self%hemisphere%posl%pos(jin)%posvec,trtfm)
  end do
  call geobjlist_tfm(self%hemisphere,transform) !positional offset

  call antenna_normal(self,beq,ray) ! calculate antenna launch directions/positions
  transform = -1 ! specify inverse transform
  call geobjlist_tfm(self%hemisphere,transform) ! reset position via inverse transform

  ! perform inverse rotation of positions via inverse rotation matrix tfm
  do jin = 1, self%hemisphere%np
     self%hemisphere%posl%pos(jin)%posvec = matmul(self%hemisphere%posl%pos(jin)%posvec,tfm)
  end do


end subroutine antenna_align

!---------------------------------------------------------------------
!> loop antenna over all CAD triangles
subroutine antenna_walltest(self,wall,beq,ray)

  !! arguments
  type(antenna_t), intent(inout) :: self !< module object
  type(geobjlist_t), intent(inout) :: wall !< triangulated wall segment
  type(beq_t), intent(inout) :: beq !< equilibrium data
  type(ray_t), intent(inout) :: ray !< ray data

  !! local
  character(*), parameter :: s_name='antenna_walltest' !< subroutine name
  real(kr4), dimension(3) :: normal !< vector to hold output from geobj_normal
  real(kr4), dimension(3) :: centre !< vector to hold output from geobj_centre
  type(posvecl_t) :: tricentre !< centre of particular triangle to be passed to antenna_align
  type(posvecl_t) :: trinormal !< normal of particular triangle to be passed to antenna_align
  integer(ki4) :: jwall !< loop counter
  type(geobj_t) :: walligeobj !< geobj_t specifying ptr and geobj type for the CAD mesh
  real(kr4) :: mag !< magnitude of normal
  integer :: kplot !< local variable
  character(10) :: kchar !< local variable
  integer(ki4) :: transform !< = 1 for forward, anything else = inverse
  real(kr8) :: psi !< local variable
  integer(ki4) :: percent!< progress percentage indicators
  integer(ki4) :: marker!< progress percentage indicators
  integer(ki4) :: per !< progress percentage indicators

  per = 10
  marker = 1+int(wall%ng/100)
  percent = marker
  !print*, "wall%ng=", wall%ng
  !print*, "10% mark=", marker
  !print*, "self%hemisphere%ng=", self%hemisphere%ng

  !! an issue where geobjlists arent read in correctly causes code to not run so call error and stop
  if (wall%ng>500000) then
     ! figure out a way for code stop without arbitrary value exceed
     print '("Fatal error: Error reading geobjlist correctly, please rerun the code")'
     call log_error(m_name,s_name,1,error_fatal,'Error reading geobjlist correctly')
  end if

  !! allocate arrays (perhaps this should be moved elsewhere)

  allocate(self%dir(self%hemisphere%ng), stat=status)
  call log_alloc_check(m_name,s_name,15,status)

  allocate(self%posdir(self%hemisphere%ng), stat=status)
  call log_alloc_check(m_name,s_name,15,status)


  call antenna_centre(self%hemisphere) ! bring hemisphere up to z=0 so it sits flush on triangle surface
  if (self%n%plot_diagvtk) then
     if (wall%ng>=100) then
        ! mesh too large do not produce antenna vtks
        self%n%plot_diagvtk = .false.
        call log_error(m_name,s_name,1,log_info,'antenna vtks not produced as mesh limit exceeded')
        print '("Warning: Wall mesh too large for antenna visualisation check. Antenna vtks not produced")'
     end if
  end if


  ! loop over every triangle in test structure
  do jwall = 1, wall%ng
     walligeobj%geobj=wall%obj2(jwall)%ptr
     walligeobj%objtyp=wall%obj2(jwall)%typ

     !    kchar = 'geometry'
     !    call misc_getfileunit(kplot)
     !    call geobjlist_writev(self%hemisphere,kchar,kplot)

     ! calculate normal of current tri
     call geobj_normal(walligeobj,wall%posl,wall%nodl,normal,mag)
     trinormal%posvec = normal
     ! calculate centre of current tri
     call geobj_centre(walligeobj,wall%posl,wall%nodl,centre)
     tricentre%posvec = centre
     ! positional offset/rotational alignment and calculation/storage
     call antenna_align(self,tricentre,trinormal,beq,ray)

     !! print statements for diagnostics

     !print*, '----------------------------Triangle', jwall
     if (jwall==percent) then
        print*,per,'% of triangles processed'
        per = per + 10
        percent = int(marker*per)
     end if

  end do

  !    kchar = 'geometry'
  !    call misc_getfileunit(kplot)
  !    call geobjlist_writev(wall,kchar,kplot)

end subroutine antenna_walltest
!---------------------------------------------------------------------
!> loop antenna over all CAD triangles
subroutine antenna_testelt(self,wall,psi,radnsource,ray,power)

  !! arguments
  type(antenna_t), intent(inout) :: self !< module object
  type(geobjlist_t), intent(inout) :: wall !< triangulated wall segment
  type(spl2d_t), intent(inout) :: psi !< equilibrium flux data
  type(radnsource_t), intent(inout) :: radnsource !< radiation source array
  type(ray_t), intent(inout) :: ray !< ray data
  real(kr4), dimension(:), intent(inout) :: power  !< power on triangles

  !! local
  character(*), parameter :: s_name='antenna_testelt' !< subroutine name
  integer(ki4) :: jwall !< loop counter
  integer(ki4) :: transform !< = 1 for forward, anything else = inverse
  integer(ki4) :: percent!< progress percentage indicators
  integer(ki4) :: marker!< progress percentage indicators
  integer(ki4) :: per !< progress percentage indicators
  integer(ki4) :: ielt !< standin for element object

  if (wall%ng>500000) then
     ! figure out a way for code stop without arbitrary value exceed
     print '("Fatal error: Error reading geobjlist correctly, please rerun the code")'
     call log_error(m_name,s_name,1,error_fatal,'Error reading geobjlist correctly')
  end if

  if (self%n%plot_diagvtk) then
     if (wall%ng>=100) then
        ! mesh too large do not produce antenna vtks
        self%n%plot_diagvtk = .false.
        call log_error(m_name,s_name,4,log_info,'antenna vtks not produced as mesh limit exceeded')
        print '("Warning: Wall mesh too large for antenna visualisation check. Antenna vtks not produced")'
     end if
  end if

  per = 10
  marker = 1+int(wall%ng/100)
  ! loop over every triangle in results geometry
  do jwall = 1, wall%ng

     ielt=jwall
     write(*,*) '--------------elt--------------', ielt !WA2
     call antenna_elt(self,ielt,wall,psi,radnsource,ray,power(jwall))

     !! print statements for diagnostics
     if (jwall==percent) then
        print*,per,'% of triangles processed'
        per = per + 10
        percent = int(marker*per)
     end if

  end do
  print*,per,'% of triangles processed'

end subroutine antenna_testelt
!---------------------------------------------------------------------
!> launch rays from antenna
subroutine antenna_elt(self,ielt,wall,psi,radnsource,ray,power)

  !! arguments
  type(antenna_t), intent(inout) :: self !< module object
  integer(ki4), intent(in) :: ielt !< standin for element object
  type(geobjlist_t), intent(inout) :: wall !< triangulated wall segment
  type(spl2d_t), intent(inout) :: psi !< equilibrium flux data
  type(radnsource_t), intent(inout) :: radnsource !< radiation source array
  type(ray_t), intent(inout) :: ray !< ray data
  real(kr4),  intent(out) :: power  !< power on antenna elt

  !! local
  character(*), parameter :: s_name='antenna_elt' !< subroutine name
  integer(ki4) :: iwall !< standin for element object
  type(geobj_t) :: walligeobj !< geobj_t specifying ptr and geobj type for the CAD mesh
  real(kr4), dimension(3) :: normal !< vector to hold output from geobj_normal
  real(kr4), dimension(3) :: centre !< vector to hold output from geobj_centre
  type(posvecl_t) :: tricentre !< centre of wall geometry triangle
  type(posvecl_t) :: trinormal !< normal of wall geometry triangle
  type(posvecl_t) :: antricentre !< centre of antenna triangle
  type(posvecl_t) :: antrinormal !< normal of antenna triangle
  real(kr4), dimension(3,3) :: tfm !< transformation matrix aligns polar axis to triangle normal
  real(kr4), dimension(3,3) :: trtfm !< matrix transpose
  real(kr4), dimension(3) :: zscale !< local variable
  real(kr4) :: n !< transformation matrix component
  real(kr4) :: n1 !< transformation matrix component
  real(kr4) :: n2 !< transformation matrix component
  real(kr4) :: n3 !< transformation matrix component
  real(kr4) :: nf !< normalisation factor if triangle normal not given as unit normal
  integer(ki4) :: jin !< loop counter
  integer(ki4) :: transform !< specify whether forwards or backwards transform
  type(geobj_t)  :: igeobj   !< local geobj definition
  real(kr4) :: mag !< magnitude of normal
  integer(ki4) :: jtri !< loop counter
  integer :: kplot !< unit number
  integer :: unorm !< unit number
  integer :: ucentr !< unit number
  character(10) :: kchar !< local variable
  real(kr4) :: zcos1 !< cosine angle of incidence
  real(kr8) :: zrecdpow !< total power received by antenna
  real(kr8) :: powtot !< total power of ray
  type(posang_t) :: posang !< local variable
  type(geobjlist_t)  :: geobjlil   !< local geobjlist copy

  iwall=ielt
  walligeobj%geobj=wall%obj2(iwall)%ptr
  walligeobj%objtyp=wall%obj2(iwall)%typ

  ! calculate normal of current tri
  call geobj_normal(walligeobj,wall%posl,wall%nodl,normal,mag)
  trinormal%posvec = normal
  ! calculate centre of current tri
  call geobj_centre(walligeobj,wall%posl,wall%nodl,centre)
  tricentre%posvec = centre

  ! positional offset/rotational alignment and calculation/storage

  !! setup transformation matrix to align antenna to given triangle

  n1 = trinormal%posvec(1)
  n2 = trinormal%posvec(2)
  n3 = trinormal%posvec(3)
  n = sqrt(n1**2 + n2**2)

  if (n3<=0) then
     ! deal with n3 < 0 case specially
     n = -n
  end if

  if (abs(n)<const_eps5*abs(n3)) then
     call log_error(m_name,s_name,1,log_info,'Problematic alignment case')
     call log_value("parameter |n| ",abs(n))
     tfm(1,1) = n3
     tfm(1,2) = 0
     tfm(1,3) = 0
     tfm(2,1) = 0
     tfm(2,2) = n3
     tfm(2,3) = 0
     tfm(3,1) = 0
     tfm(3,2) = 0
     tfm(3,3) = n3
  else
     tfm(1,1) = n1*n3/n
     tfm(1,2) = -n2/n
     tfm(1,3) = n1
     tfm(2,1) = n2*n3/n
     tfm(2,2) = n1/n
     tfm(2,3) = n2
     tfm(3,1) = -n
     tfm(3,2) = 0
     tfm(3,3) = n3
  end if
  call geobjlist_copy(self%hemisphere,geobjlil,0)

  ! offset and rotate positions via rotation matrix tfm
  geobjlil%tfmdata%matrix=tfm
  ! set scale
  zscale = (/1,1,1/)
  geobjlil%tfmdata%scale = zscale
  ! set transform 3 for position_tfm
  geobjlil%tfmdata%ntfm = 3
  ! set offset vector for transform
  geobjlil%tfmdata%offset = tricentre%posvec

  call geobjlist_tfm(geobjlil,1)

  zrecdpow=0
  ! loop over every triangle in antenna
  do jtri = 1, geobjlil%ng
     !WA1     write(*,*) 'jtri===============',jtri !WA1
     igeobj%geobj=geobjlil%obj2(jtri)%ptr
     igeobj%objtyp=geobjlil%obj2(jtri)%typ
     ! calculate tri centre
     call geobj_centre(igeobj,geobjlil%posl,geobjlil%nodl,centre)
     !convert to metres
     posang%pos=centre
     posang%opt=0
     posang%units=-3
     call posang_units(posang,0)
     ! store ray launch position as tri centre
     antricentre%posvec = posang%pos
     self%posdir(jtri)%posvec = centre !DIAG
     !write(ucentr,cfmtbv0c) self%posdir(jtri)%posvec !DIAG
     ! calculate tri normal
     call geobj_normal(igeobj,geobjlil%posl,geobjlil%nodl,normal,mag)
     ! store ray launch direction as tri normal
     antrinormal%posvec = normal
     self%dir(jtri)%posvec = normal !DIAG
     !write(unorm,cfmtbv0c)self%dir(jtri)%posvec !DIAG
     ! track ray to end
     call ray_trackinlcfs(ray,psi,radnsource%srcspl2d,antricentre,antrinormal)
     ! calculate power contribution
     call ray_trackpower(ray,psi,radnsource%srcspl2d,powtot)
     ! accumulate over antenna solid angle
     zcos1=dot_product(antrinormal%posvec,trinormal%posvec)
     !DIA1 write(*,*) 'recdpow,pow,area,cos',zrecdpow,powtot,self%area(jtri),zcos1 !DIA1
     zrecdpow=zrecdpow + powtot*self%area(jtri)*zcos1

  end do

  !call geobjlist_tfm(geobjlil,-1)

  power=zrecdpow

end subroutine antenna_elt
!---------------------------------------------------------------------
!> launch rays from antenna/powelt
subroutine antenna_powelt(self,ielt,wall,psi,radnsource,ray,gshadl,powr,btree)

  !! arguments
  type(antenna_t), intent(inout) :: self !< module object
  type(powelt_t), intent(in) :: ielt !< element object
  type(geobjlist_t), intent(inout) :: wall !< triangulated wall segment
  type(spl2d_t), intent(inout) :: psi !< equilibrium flux data
  type(radnsource_t), intent(inout) :: radnsource !< radiation source array
  type(ray_t), intent(inout) :: ray !< ray data
  type(geobjlist_t), intent(inout) :: gshadl !< shadow geobj list data
  real(kr4), dimension(:), intent(out) :: powr  !< radiation power on elt
  type(btree_t), intent(in), optional :: btree   !< binary tree data

  !! local
  character(*), parameter :: s_name='antenna_powelt' !< subroutine name
  integer(ki4) :: iwall !< address of wall triangular element
  type(geobj_t) :: walligeobj !< geobj_t specifying ptr and geobj type for the CAD mesh
  real(kr4), dimension(3) :: normal !< vector to hold output from geobj_normal
  real(kr4), dimension(3) :: centre !< vector to hold output from geobj_centre
  type(posvecl_t) :: tricentre !< centre of wall geometry triangle
  type(posvecl_t) :: trinormal !< normal of wall geometry triangle
  type(posvecl_t) :: antricentre !< centre of antenna triangle
  type(posvecl_t) :: antrinormal !< normal of antenna triangle
  real(kr4), dimension(3,3) :: tfm !< transformation matrix aligns polar axis to triangle normal
  real(kr4), dimension(3,3) :: trtfm !< matrix transpose
  real(kr4), dimension(3) :: zscale !< local variable
  real(kr4) :: n !< transformation matrix component
  real(kr4) :: n1 !< transformation matrix component
  real(kr4) :: n2 !< transformation matrix component
  real(kr4) :: n3 !< transformation matrix component
  real(kr4) :: nf !< normalisation factor if triangle normal not given as unit normal
  integer(ki4) :: jin !< loop counter
  integer(ki4) :: transform !< specify whether forwards or backwards transform
  type(geobj_t)  :: igeobj   !< local geobj definition
  real(kr4) :: mag !< magnitude of normal
  integer(ki4) :: jtri !< loop counter
  integer :: kplot !< unit number
  integer :: unorm !< unit number
  integer :: ucentr !< unit number
  character(10) :: kchar !< local variable
  real(kr8) :: zrecdpow !< total power received by antenna
  real(kr8) :: powtot !< total power of ray
  type(posang_t) :: posang !< local variable
  type(geobjlist_t)  :: geobjlil   !< local geobjlist copy
  !integer(ki4) :: npowe   !< powelt address factor
  integer(ki4), parameter :: nshadowr=1 !< if nonzero, shadow radiation power

  iwall=ielt%ie
  walligeobj%geobj=wall%obj2(iwall)%ptr
  walligeobj%objtyp=wall%obj2(iwall)%typ

  ! calculate normal of current tri
  call geobj_normal(walligeobj,wall%posl,wall%nodl,normal,mag)
  trinormal%posvec = normal
  ! calculate centre of current tri
  call geobj_centre(walligeobj,wall%posl,wall%nodl,centre)
  tricentre%posvec = centre

  ! positional offset/rotational alignment and calculation/storage

  !! setup transformation matrix to align antenna to given triangle

  n1 = trinormal%posvec(1)
  n2 = trinormal%posvec(2)
  n3 = trinormal%posvec(3)
  n = sqrt(n1**2 + n2**2)

  if (n3<=0) then
     ! deal with n3 < 0 case specially
     n = -n
  end if

  if (abs(n)<const_eps5*abs(n3)) then
     call log_error(m_name,s_name,1,log_info,'Problematic alignment case')
     call log_value("parameter |n| ",abs(n))
     tfm(1,1) = n3
     tfm(1,2) = 0
     tfm(1,3) = 0
     tfm(2,1) = 0
     tfm(2,2) = n3
     tfm(2,3) = 0
     tfm(3,1) = 0
     tfm(3,2) = 0
     tfm(3,3) = n3
  else
     tfm(1,1) = n1*n3/n
     tfm(1,2) = -n2/n
     tfm(1,3) = n1
     tfm(2,1) = n2*n3/n
     tfm(2,2) = n1/n
     tfm(2,3) = n2
     tfm(3,1) = -n
     tfm(3,2) = 0
     tfm(3,3) = n3
  end if
  call geobjlist_copy(self%hemisphere,geobjlil,0)

  ! offset and rotate positions via rotation matrix tfm
  geobjlil%tfmdata%matrix=tfm
  ! set scale
  zscale = (/1,1,1/)
  geobjlil%tfmdata%scale = zscale
  ! set transform 3 for position_tfm
  geobjlil%tfmdata%ntfm = 3
  ! set offset vector for transform
  geobjlil%tfmdata%offset = tricentre%posvec

  call geobjlist_tfm(geobjlil,1)

  zrecdpow=0
  ! loop over every triangle in antenna
  do jtri = 1, geobjlil%ng
     !WA1     write(*,*) 'jtri===============',jtri !WA1
     igeobj%geobj=geobjlil%obj2(jtri)%ptr
     igeobj%objtyp=geobjlil%obj2(jtri)%typ
     ! calculate tri centre
     call geobj_centre(igeobj,geobjlil%posl,geobjlil%nodl,centre)
     !convert to metres
     posang%pos=centre
     posang%opt=0
     posang%units=-3
     call posang_units(posang,0)
     ! store ray launch position as tri centre
     antricentre%posvec = posang%pos
     self%posdir(jtri)%posvec = centre !DIAG
     !write(ucentr,cfmtbv0c) self%posdir(jtri)%posvec !DIAG
     ! calculate tri normal
     call geobj_normal(igeobj,geobjlil%posl,geobjlil%nodl,normal,mag)
     ! store ray launch direction as tri normal
     antrinormal%posvec = normal
     self%dir(jtri)%posvec = normal !DIAG
     !write(unorm,cfmtbv0c)self%dir(jtri)%posvec !DIAG
     ! track ray to end
     if (nshadowr/=0) then
        call ray_trackcoll(ray,psi,radnsource%srcspl2d,gshadl,btree,&
 &      antricentre,antrinormal,-1)
     else
        call ray_trackinlcfs(ray,psi,radnsource%srcspl2d,antricentre,antrinormal)
     end if
     ! calculate power contribution
     call ray_trackpower(ray,psi,radnsource%srcspl2d,powtot)
     ! accumulate over antenna solid angle
     zrecdpow=zrecdpow + powtot*self%area(jtri)

  end do

  !call geobjlist_tfm(geobjlil,-1)

  powr(iwall)=zrecdpow

end subroutine antenna_powelt
!---------------------------------------------------------------------
!> output to log file
subroutine antenna_dia(self)

  !! arguments
  type(antenna_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='antenna_dia' !< subroutine name

  call log_value("ntheta ",self%n%ntheta)
  call log_value("power ",self%totpow)

end subroutine antenna_dia
!---------------------------------------------------------------------
!> userdefined function
function antenna_powflux(psi,psix)

  !! arguments
  real(kr4) :: antenna_powflux !< local variable
  real(kr4), intent(in) :: psi !< position in \f$ \psi \f$
  real(kr4), intent(in) :: psix !< beq structure

  !! local variables
  character(*), parameter :: s_name='antenna_userdefined' !< subroutine name
  real(kr4) :: pow !< local variable
  real(kr4) :: psimin !< local variable
  real(kr4) :: psimax !< local variable
  real(kr4) :: pmax !< local variable

  pmax = 5

  pow = pmax - (psi - psix)**2
  antenna_powflux = pow

end function antenna_powflux

!---------------------------------------------------------------------
!> userdefined function
function antenna_userdefined(self,psi)

  !! arguments
  type(antenna_t), intent(in) :: self !< module object
  real(kr8) :: antenna_userdefined !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='antenna_userdefined' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  integer(ki4) :: ilocal !< local integer variable

  zpos=psi
  pow=0._kr8
  !> user defines \Tt{pow} here
  !! .....
  !! return antenna
  antenna_userdefined=pow

end function antenna_userdefined

!---------------------------------------------------------------------
!> general external function call
function antenna_fn(self,psi)

  !! arguments
  type(antenna_t), intent(in) :: self !< module object
  real(kr8) :: antenna_fn !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='antenna_fn' !< subroutine name
  real(kr8) :: pow !< local variable

  pow=0._kr8
  !! select antenna
  formula_chosen: select case (self%n%formula)
  case('userdefined')
     pow=antenna_userdefined(self,psi)
  end select formula_chosen

  !! return antenna
  antenna_fn=pow

end function antenna_fn
!---------------------------------------------------------------------
!> open new file, making up name
subroutine antenna_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer(ki4), intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='antenna_initwrite' !< subroutine name
  logical :: unitused !< flag to test unit is available
  character(len=80) :: outputfile !< output file name

  !! get file unit
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        if (present(channel)) channel=i
        exit
     end if
  end do
  noutant=i

  !! open file
  outputfile=trim(fileroot)//"_antenna.out"
  call log_value("Control data file",trim(outputfile))
  open(unit=noutant,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutant,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if

end subroutine antenna_initwrite
!---------------------------------------------------------------------
!> write antenna data
subroutine antenna_write(self,channel)

  !! arguments
  type(antenna_t), intent(in) :: self   !< antenna data structure
  integer(ki4), intent(in), optional :: channel   !< output channel for antenna data structure

  !! local
  character(*), parameter :: s_name='antenna_write' !< subroutine name
  integer(ki4) :: iout   !< output channel for antenna data structure

  !! sort out unit
  if(present(channel)) then
     iout=channel
  else
     iout=noutant
  end if

  write(iout,*,iostat=status) 'antenna_formula'
  call log_write_check(m_name,s_name,18,status)
  write(iout,*,iostat=status) self%n%formula
  call log_write_check(m_name,s_name,19,status)

  write(iout,*,iostat=status) 'ntheta'
  call log_write_check(m_name,s_name,20,status)                      !! ntheta
  write(iout,*,iostat=status)  self%n%ntheta
  call log_write_check(m_name,s_name,21,status)

  write(iout,*,iostat=status) 'nphi'
  call log_write_check(m_name,s_name,22,status)                      !! nphi
  write(iout,*,iostat=status)  self%n%nphi
  call log_write_check(m_name,s_name,23,status)


  !  write(iout,*,iostat=status) 'f'
  !  call log_write_check(m_name,s_name,26,status)
  !  write(iout,*,iostat=status) self%n%f
  !  call log_write_check(m_name,s_name,27,status)

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

end subroutine antenna_write
!---------------------------------------------------------------------
!> write object data as gnuplot
subroutine antenna_writeg(self,select,channel)

  !! arguments
  type(antenna_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for antenna data structure

  !! local
  character(*), parameter :: s_name='antenna_writeg' !< subroutine name
  integer(ki4) :: iout   !< output channel for antenna data structure

  call log_error(m_name,s_name,1,log_info,'gnuplot file produced')

end subroutine antenna_writeg
!---------------------------------------------------------------------
!> write object data as vtk
subroutine antenna_writev(self,select,channel)

  !! arguments
  type(antenna_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for antenna data structure

  !! local
  character(*), parameter :: s_name='antenna_writev' !< subroutine name
  integer(ki4) :: iout   !< output channel for antenna data structure

  call log_error(m_name,s_name,1,log_info,'vtk file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine antenna_writev
!---------------------------------------------------------------------
!> close write file
subroutine antenna_closewrite

  !! local
  character(*), parameter :: s_name='antenna_closewrite' !< subroutine name

  !! close file
  close(unit=noutant,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine antenna_closewrite
!---------------------------------------------------------------------
!> delete object
subroutine antenna_delete(self)

  !! arguments
  type(antenna_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='antenna_delete' !< subroutine name


  formula_deallocate: select case (self%n%formula)
  case('userdefined')
     if (self%n%nrpams>0) deallocate(self%n%rpar)
     if (self%n%nipams>0) deallocate(self%n%npar)
  case default
  end select formula_deallocate

  if (allocated(self%dir)) deallocate(self%dir)
  if (allocated(self%posdir)) deallocate(self%posdir)

end subroutine antenna_delete
!---------------------------------------------------------------------
!> close file
subroutine antenna_close

  !! local
  character(*), parameter :: s_name='antenna_close' !< subroutine name

  !! close file
  close(unit=ninant,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine antenna_close

end module antenna_m
