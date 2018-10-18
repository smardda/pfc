module fmesh_m

  use position_h
  use fmesh_h
  use log_m
  use misc_m
  use const_kind_m
  use position_m

  implicit none
  private

! public subroutines
  public :: &
  fmesh_init,  & !< open file
  fmesh_readcon,  & !< read data from file
  fmesh_uniform,  & !< calculate uniform field mesh
  fmesh_userdefined,  & !< userdefined analytic field mesh
  fmesh_generic,  & !< generic analytic call
  fmesh_initwrite, & !< open new output file
  fmesh_copy, &  !< copy object
  fmesh_write, &  !< write out object
  fmesh_read, &  !< read in object
  fmesh_delete, & !< delete object
  fmesh_close, & !< close file
  fmesh_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='fmesh_m' !< module name
  integer  :: status   !< error status
  integer,save  :: ninfm=5      !< control file unit number
  integer,save  :: noutfm=6      !< output file unit number
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  character(len=80) :: ibuff !< buffer for input/output

  contains
!---------------------------------------------------------------------
!> open file
subroutine fmesh_init(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer, intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='fmesh_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then channel=i exit end if end do ninfm=i

  !! open file
  controlfile=trim(file)
  call log_value("Mesh data file",trim(controlfile))
  call misc_getfileunit(ninfm)
  open(unit=ninfm,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open mesh data file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open mesh data data file')
     stop
  end if
  if (present(channel)) channel=ninfm

end  subroutine fmesh_init
!---------------------------------------------------------------------
!> read data from file
subroutine fmesh_readcon(self,kin)

  !! arguments
  type(fmesh_t), intent(out) :: self !< type which data will be assigned to
  integer, intent(in),optional :: kin   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='fmesh_readcon' !< subroutine name
  character(len=80) :: mesh_formula !< mesh formula
  logical :: mesh_transform !< mesh transform
  integer(ki4)  :: dimensionality !< dimension of mesh array (1, 2 or 3)
  integer(ki4)  :: nstaggered !< whether mesh staggered (1) or not (0)
  character(len=80)  :: length_units !< units, either me(tres) or mm
  character(len=80)  :: plot_units !< units, either me(tres) or mm
  real(kr8) :: length_x !< physical extent of mesh
  real(kr8) :: length_y !< physical extent of mesh
  real(kr8) :: length_z !< physical extent of mesh
  real(kr8) :: xstart !< origin of mesh in first direction
  real(kr8) :: ystart !< origin of mesh in second direction
  real(kr8) :: zstart !< origin of mesh in third direction
  integer(ki4) :: nx !< input field mesh size
  integer(ki4) :: ny !< input field mesh size
  integer(ki4) :: nz !< input field mesh size

  !     real(kr8), dimension (:), allocatable:: xvalues !< input field mesh
  !     real(kr8), dimension (:), allocatable:: yvalues !< input field mesh
  !     real(kr8), dimension (:), allocatable:: zvalues !< input field mesh
  real(kr8) :: dx !< input field mesh spacing size
  real(kr8) :: dy !< input field mesh spacing size
  real(kr8) :: dz !< input field mesh spacing size
  real(kr8) :: xend !< far point of mesh in first direction
  real(kr8) :: yend !< far point of mesh in second direction
  real(kr8) :: zend !< far point of mesh in third direction
  real(kr4) :: faclen=1. !< lengths default is metre

  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters !< general real parameters
  integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters !< general integer_parameters
  integer(ki4) :: number_of_real_parameters !< number of real parameters
  integer(ki4) :: number_of_integer_parameters !< number of integer parameters
  integer(ki4) :: flag  !< local variable

  !! mesh parameters
  namelist /fmeshparameters/ &
 &mesh_formula , mesh_transform , &
 &dimensionality , nstaggered, &
 &length_units, &
 &plot_units, &
 &length_x , length_y , length_z , &
 &xstart , ystart , zstart, xend , yend, zend, &
 &nx , ny , nz , &
 &dx , dy , dz
  ! &xvalues , yvalues , zvalues , &

  !! set default mesh parameters
  mesh_formula='regular'
  mesh_transform=.TRUE.
  dimensionality=3
  nstaggered=0
  length_units='metre'
  plot_units='millimetre'
  length_x=1.
  length_y=1.
  length_z=1.
  xstart=0.
  ystart=0.
  zstart=0.
  nx=26
  ny=26
  nz=26
  ! if mesh_formula='increment', use instead of length_s
  dx=.04
  dy=.04
  dz=.04
  ! if mesh_formula='limits', use instead of length_s
  xend=1.
  yend=1.
  zend=1.

  !! userdefined formula
  general_real_parameters=0
  general_integer_parameters=0
  number_of_real_parameters=0
  number_of_integer_parameters=0

  if(present(kin)) then
     !! assume unit already open and reading infile
     ninfm=kin
  end if

  !!read mesh parameters
  read(ninfm,nml=fmeshparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_warning,'Error reading mesh parameters')
     call log_getunit(ilog)
     write(ilog,nml=fmeshparameters)
     print '("Error reading mesh parameters")'
  end if

  !! check inputs
  call lowor(mesh_formula,1,len_trim(mesh_formula))
  if (mesh_formula=='regular') then

     if (dimensionality<=0.OR.dimensionality>3) then
        call log_error(m_name,s_name,2,error_fatal,'Mesh should be 1-D, 2-D or 3-D')
     end if
     if (nx<=0) then
        call log_error(m_name,s_name,4,error_fatal,'need positive number of mesh points in x')
     end if
     if (ny<=0.AND.dimensionality>1) then
        call log_error(m_name,s_name,5,error_fatal,'need positive number of mesh points in y')
     end if
     if (nz<=0.AND.dimensionality>2) then
        call log_error(m_name,s_name,6,error_fatal,'need positive number of mesh points in z')
     end if
     !Other tests !todo

  end if

  !! store mesh values
  self%formula=mesh_formula
  !! store field mesh values
  self%ndimf=dimensionality
  self%nstagf=nstaggered
  self%lunit=length_units(1:2)
  self%punit=plot_units(1:2)

  mesh_form: select case (mesh_formula)
  case ('regular')
     dx=length_x/(nx-1)
     dy=length_y/(ny-1)
     dz=length_z/(nz-1)
     xend=xstart+length_x
     yend=ystart+length_y
     zend=zstart+length_z
     mesh_formula='uniform'

  case ('increment')
     length_x=(nx-1)*dx
     length_y=(ny-1)*dy
     length_z=(nz-1)*dz
     xend=xstart+length_x
     yend=ystart+length_y
     zend=zstart+length_z
     mesh_formula='uniform'

  case ('limits')
     length_x= xend- xstart
     length_y= yend- ystart
     length_z= zend- zstart
     dx=length_x/(nx-1)
     dy=length_y/(ny-1)
     dz=length_z/(nz-1)
     mesh_formula='uniform'

  case ('userdefined')
     if(number_of_real_parameters<0) &
 &   call log_error(m_name,s_name,14,error_fatal,'number of real parameters must be >=0')
     if(number_of_real_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of real_p rameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,15,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters<0) &
 &   call log_error(m_name,s_name,16,error_fatal,'number of integer parameters must be >=0')
     if(number_of_integer_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of integer parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,17,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters==0.AND.number_of_real_parameters==0) &
 &   call log_error(m_name,s_name,18,error_fatal,'no parameters set')

  case default
     call log_error(m_name,s_name,20,error_fatal,'Mesh with this formula not implemented')
  end select mesh_form

  !! scale if required
  if (length_units(1:2)=='me') then
     faclen=1.
  else
     faclen=0.001
     xstart=faclen*xstart; ystart=faclen*ystart; zstart=faclen*zstart
     xend=faclen*xend; yend=faclen*yend; zend=faclen*zend
     dx=faclen*dx; dy=faclen*dy; dz=faclen*dz
     length_x=faclen*length_x; length_y=faclen*length_y; length_z=faclen*length_z
     length_units='millimetre'
  end if

  self%xlengthf=length_x
  self%ylengthf=length_y
  self%zlengthf=length_z
  self%x0f=xstart
  self%y0f=ystart
  self%z0f=zstart
  self%nxf=nx
  self%nyf=ny
  self%nzf=nz
  self%dxf=dx
  self%dyf=dy
  self%dzf=dz
  self%x1f=xend
  self%y1f=yend
  self%z1f=zend
  self%nrpams=number_of_real_parameters
  self%nipams=number_of_integer_parameters

  !! allocate input arrays and assign

  formula_allocate: select case (mesh_formula)
  case('userdefined')
     if (number_of_real_parameters>0) allocate(self%rpar(number_of_real_parameters), stat=status)
     call log_alloc_check(m_name,s_name,65,status)
     if (number_of_integer_parameters>0) allocate(self%npar(number_of_integer_parameters), stat=status)
     call log_alloc_check(m_name,s_name,66,status)
     self%rpar=general_real_parameters(:number_of_real_parameters)
     self%npar=general_integer_parameters(:number_of_integer_parameters)
     self%formula='userdefined'
  case default
  end select formula_allocate

  self%iltfm=mesh_transform
  if (mesh_transform) then
     !! read in transform
     call position_readcon(self%tfmdata,ninfm)
  else
     !! set to identity
     call position_unitfm(self%tfmdata)
  end if

  !! scale offset as required
  call position_scaleunits(self%tfmdata,faclen)

end  subroutine fmesh_readcon
!---------------------------------------------------------------------
!> calculate uniform mesh field
subroutine fmesh_uniform(self)

  !! arguments
  type(fmesh_t), intent(inout) :: self !< type containing field mesh definition parameters

  !! local variables
  character(*), parameter :: s_name='fmesh_uniform' !< subroutine name

  !! check mesh formula !todo

  allocate(self%xf(self%nxf), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  do i=1,self%nxf
     self%xf(i)=self%x0f+(i-1)*self%dxf
  end do

  if (self%ndimf>1) then
     allocate(self%yf(self%nyf), stat=status)
     call log_alloc_check(m_name,s_name,2,status)
     do j=1,self%nyf
        self%yf(j)=self%y0f+(j-1)*self%dyf
     end do
  end if

  if (self%ndimf>2) then
     allocate(self%zf(self%nzf), stat=status)
     call log_alloc_check(m_name,s_name,3,status)
     do k=1,self%nzf
        self%zf(k)=self%z0f+(k-1)*self%dzf
     end do
  end if

end subroutine fmesh_uniform
!---------------------------------------------------------------------
!> user defined analytic field mesh
subroutine fmesh_userdefined(self)

  !! arguments
  type(fmesh_t), intent(inout) :: self !< type containing field mesh definition parameters

  !! local variables
  character(*), parameter :: s_name='fmesh_userdefined' !< subroutine name

  !! user defines field mesh here

end subroutine fmesh_userdefined
!>---------------------------------------------------------------------
subroutine fmesh_generic(self)

  !! arguments
  type(fmesh_t), intent(inout) :: self !< type containing field mesh definition parameters

  !! local variables
  character(*), parameter :: s_name='fmesh_generic' !< subroutine name

  !! select formula
  formula_chosen: select case (self%formula)
  case('uniform')
     call fmesh_uniform(self)
  case('userdefined')
     call fmesh_userdefined(self)
  end select formula_chosen

end subroutine fmesh_generic
!---------------------------------------------------------------------
!> open new file
subroutine fmesh_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer, intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='fmesh_initwrite' !< subroutine name
  character(len=80) :: outputfile !< output file name
  !! logical :: unitused !< flag to test unit is available
  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then channel=i exit end if end do noutfm=i

  !! open file
  outputfile=trim(fileroot)//"_fmesh.out"
  call log_value("Mesh data file",trim(outputfile))
  call misc_getfileunit(noutfm)
  open(unit=noutfm,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutfm,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if
  if (present(channel)) channel=noutfm

end subroutine fmesh_initwrite
!---------------------------------------------------------------------
!> copy fmesh data
subroutine fmesh_copy(selfi,selfo)

  !! arguments
  type(fmesh_t), intent(in) :: selfi   !< input fmesh data structure
  type(fmesh_t), intent(out) :: selfo   !< output fmesh data structure

  !! local
  character(*), parameter :: s_name='fmesh_copy' !< subroutine name

  selfo=selfi

end subroutine fmesh_copy
!---------------------------------------------------------------------
!> write fmesh data
subroutine fmesh_write(self,kout)

  !! arguments
  type(fmesh_t), intent(in) :: self   !< fmesh data structure
  integer, intent(in), optional :: kout   !< output channel for fmesh data structure

  !! local
  character(*), parameter :: s_name='fmesh_write' !< subroutine name
  integer :: iout   !< output channel for fmesh data structure

  !! sort out unit
  if(present(kout)) then
     iout=kout
  else
     iout=noutfm
  end if

  write(iout,*,iostat=status) 'mesh_formula'
  call log_write_check(m_name,s_name,1,status)
  write(iout,*,iostat=status) self%formula
  call log_write_check(m_name,s_name,2,status)
  write(iout,*,iostat=status) 'ndimf'
  call log_write_check(m_name,s_name,3,status)
  write(iout,*,iostat=status) self%ndimf
  call log_write_check(m_name,s_name,4,status)
  write(iout,*,iostat=status) 'nstagf'
  call log_write_check(m_name,s_name,5,status)
  write(iout,*,iostat=status) self%nstagf
  call log_write_check(m_name,s_name,12,status)
  write(iout,*,iostat=status) 'nxf'
  call log_write_check(m_name,s_name,13,status)
  write(iout,*,iostat=status) self%nxf
  call log_write_check(m_name,s_name,14,status)
  write(iout,*,iostat=status) 'nyf'
  call log_write_check(m_name,s_name,15,status)
  write(iout,*,iostat=status) self%nyf
  call log_write_check(m_name,s_name,16,status)
  write(iout,*,iostat=status) 'nzf'
  call log_write_check(m_name,s_name,17,status)
  write(iout,*,iostat=status) self%nzf
  call log_write_check(m_name,s_name,18,status)
  write(iout,*,iostat=status) 'x0f'
  call log_write_check(m_name,s_name,23,status)
  write(iout,*,iostat=status) self%x0f
  call log_write_check(m_name,s_name,24,status)
  write(iout,*,iostat=status) 'y0f'
  call log_write_check(m_name,s_name,25,status)
  write(iout,*,iostat=status) self%y0f
  call log_write_check(m_name,s_name,26,status)
  write(iout,*,iostat=status) 'z0f'
  call log_write_check(m_name,s_name,27,status)
  write(iout,*,iostat=status) self%z0f
  call log_write_check(m_name,s_name,28,status)
  write(iout,*,iostat=status) 'xlengthf'
  call log_write_check(m_name,s_name,29,status)
  write(iout,*,iostat=status) self%xlengthf
  call log_write_check(m_name,s_name,30,status)
  write(iout,*,iostat=status) 'ylengthf'
  call log_write_check(m_name,s_name,31,status)
  write(iout,*,iostat=status) self%ylengthf
  call log_write_check(m_name,s_name,32,status)
  write(iout,*,iostat=status) 'zlengthf'
  call log_write_check(m_name,s_name,33,status)
  write(iout,*,iostat=status) self%zlengthf
  call log_write_check(m_name,s_name,34,status)
  write(iout,*,iostat=status) 'lunit'
  call log_write_check(m_name,s_name,35,status)
  write(iout,*,iostat=status) self%lunit
  call log_write_check(m_name,s_name,36,status)
  write(iout,*,iostat=status) 'dxf'
  call log_write_check(m_name,s_name,37,status)
  write(iout,*,iostat=status) self%dxf
  call log_write_check(m_name,s_name,38,status)
  write(iout,*,iostat=status) 'dyf'
  call log_write_check(m_name,s_name,39,status)
  write(iout,*,iostat=status) self%dyf
  call log_write_check(m_name,s_name,40,status)
  write(iout,*,iostat=status) 'dzf'
  call log_write_check(m_name,s_name,41,status)
  write(iout,*,iostat=status) self%dzf
  call log_write_check(m_name,s_name,42,status)
  write(iout,*,iostat=status) 'x1f'
  call log_write_check(m_name,s_name,43,status)
  write(iout,*,iostat=status) self%x1f
  call log_write_check(m_name,s_name,44,status)
  write(iout,*,iostat=status) 'y1f'
  call log_write_check(m_name,s_name,45,status)
  write(iout,*,iostat=status) self%y1f
  call log_write_check(m_name,s_name,46,status)
  write(iout,*,iostat=status) 'z1f'
  call log_write_check(m_name,s_name,47,status)
  write(iout,*,iostat=status) self%z1f
  call log_write_check(m_name,s_name,48,status)

  call position_writetfm(self%tfmdata,iout)

  if (self%formula=='uniform') return

  write(iout,*,iostat=status) 'nrpams'
  call log_write_check(m_name,s_name,56,status)
  write(iout,*,iostat=status) self%nrpams
  call log_write_check(m_name,s_name,57,status)
  if (self%nrpams>0) then
     write(iout,*,iostat=status) 'real_parameters'
     call log_write_check(m_name,s_name,58,status)
     write(iout,*,iostat=status) self%rpar
     call log_write_check(m_name,s_name,59,status)
  end if

  write(iout,*,iostat=status) 'nipams'
  call log_write_check(m_name,s_name,60,status)
  write(iout,*,iostat=status) self%nipams
  call log_write_check(m_name,s_name,61,status)
  if (self%nipams>0) then
     write(iout,*,iostat=status) 'integer_parameters'
     call log_write_check(m_name,s_name,62,status)
     write(iout,*,iostat=status) self%npar
     call log_write_check(m_name,s_name,63,status)
  end if

end subroutine fmesh_write
!---------------------------------------------------------------------
!> read fmesh data
subroutine fmesh_read(self,infile,channel)

  !! arguments
  type(fmesh_t), intent(out) :: self   !< fmesh data structure
  character(len=80),intent(in) :: infile !< name of input file
  integer, intent(in),optional :: channel   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='fmesh_read' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  if(present(channel).AND.channel/=0) then
     !! assume unit already open and reading infile
     ninfm=channel
  else if(present(channel).AND.channel==0) then
     !! assume unit already open and know unit
     continue
  else
     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then ninfm=i exit end if end do

     !! open file
     call misc_getfileunit(ninfm)
     open(unit=ninfm,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        call log_error(m_name,s_name,1,error_fatal,'Error opening data structure file')
     else
        call log_error(m_name,s_name,2,log_info,'data structure file opened')
     end if
  end if

  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,1,status)
  read(ninfm,*,iostat=status) self%formula
  call log_read_check(m_name,s_name,2,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,3,status)
  read(ninfm,*,iostat=status) self%ndimf
  call log_read_check(m_name,s_name,4,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,5,status)
  read(ninfm,*,iostat=status) self%nstagf
  call log_read_check(m_name,s_name,12,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,13,status)
  read(ninfm,*,iostat=status) self%nxf
  call log_read_check(m_name,s_name,14,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,15,status)
  read(ninfm,*,iostat=status) self%nyf
  call log_read_check(m_name,s_name,16,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,17,status)
  read(ninfm,*,iostat=status) self%nzf
  call log_read_check(m_name,s_name,18,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,23,status)
  read(ninfm,*,iostat=status) self%x0f
  call log_read_check(m_name,s_name,24,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,25,status)
  read(ninfm,*,iostat=status) self%y0f
  call log_read_check(m_name,s_name,26,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,27,status)
  read(ninfm,*,iostat=status) self%z0f
  call log_read_check(m_name,s_name,28,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,29,status)
  read(ninfm,*,iostat=status) self%xlengthf
  call log_read_check(m_name,s_name,30,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,31,status)
  read(ninfm,*,iostat=status) self%ylengthf
  call log_read_check(m_name,s_name,32,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,33,status)
  read(ninfm,*,iostat=status) self%zlengthf
  call log_read_check(m_name,s_name,34,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,35,status)
  read(ninfm,*,iostat=status) self%lunit
  call log_read_check(m_name,s_name,36,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,37,status)
  read(ninfm,*,iostat=status) self%dxf
  call log_read_check(m_name,s_name,38,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,39,status)
  read(ninfm,*,iostat=status) self%dyf
  call log_read_check(m_name,s_name,40,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,41,status)
  read(ninfm,*,iostat=status) self%dzf
  call log_read_check(m_name,s_name,42,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,43,status)
  read(ninfm,*,iostat=status) self%x1f
  call log_read_check(m_name,s_name,44,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,45,status)
  read(ninfm,*,iostat=status) self%y1f
  call log_read_check(m_name,s_name,46,status)
  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,47,status)
  read(ninfm,*,iostat=status) self%z1f
  call log_read_check(m_name,s_name,48,status)

  !    write(*,*) 'self%dxf',self%dxf
  call position_readtfm(self%tfmdata,ninfm)

  if (self%formula=='uniform') return

  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,56,status)
  read(ninfm,*,iostat=status) self%nrpams
  call log_read_check(m_name,s_name,57,status)
  if (self%nrpams>0) then
     read(ninfm,*,iostat=status) ibuff
     call log_read_check(m_name,s_name,58,status)
     read(ninfm,*,iostat=status) self%rpar
     call log_read_check(m_name,s_name,59,status)
  end if

  read(ninfm,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,60,status)
  read(ninfm,*,iostat=status) self%nipams
  call log_read_check(m_name,s_name,61,status)
  if (self%nipams>0) then
     read(ninfm,*,iostat=status) ibuff
     call log_read_check(m_name,s_name,62,status)
     read(ninfm,*,iostat=status) self%npar
     call log_read_check(m_name,s_name,63,status)
  end if

end subroutine fmesh_read
!---------------------------------------------------------------------
!> delete object
subroutine fmesh_delete(self)

  !! arguments
  type(fmesh_t), intent(inout) :: self !< type containing analytic mesh parameters
  !! local
  character(*), parameter :: s_name='fmesh_delete' !< subroutine name

  formula_deallocate: select case (self%formula)
  case('uniform')
     deallocate(self%xf)
     if (self%ndimf>1) deallocate(self%yf)
     if (self%ndimf>2) deallocate(self%zf)
  case('userdefined')
     if (self%nrpams>0) deallocate(self%rpar)
     if (self%nipams>0) deallocate(self%npar)
  case default
  end select formula_deallocate

end subroutine fmesh_delete
!---------------------------------------------------------------------
!> close file
subroutine fmesh_close

  !! local
  character(*), parameter :: s_name='fmesh_close' !< subroutine name

  !! close file
  close(unit=ninfm,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close mesh data file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close mesh data data file')
     stop
  end if

end subroutine fmesh_close
!---------------------------------------------------------------------
!> close write file
subroutine fmesh_closewrite

  !! local
  character(*), parameter :: s_name='fmesh_closewrite' !< subroutine name

  !! close file
  close(unit=noutfm,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine fmesh_closewrite

end module fmesh_m
