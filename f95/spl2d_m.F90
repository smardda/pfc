module spl2d_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use misc_m
  use position_h
  use position_m

  implicit none
  private

! public subroutines
  public :: spl2d_read, & !< read in data structure format
  spl2d_init,   & !< create spl2d data structure and knots
  spl2d_initpart,   & !< create create dummy part of spl2d data structure
  spl2d_ptlimits,   & !< get \f$ \psi,\theta \f$ limits of interpolation
  spl2d_ptscale,  &  !< scale points at which spline defined
  spl2d_scale,  &  !< scale spline values
  spl2d_delete,   & !< delete  spl2d data structure
  spl2d_writeg, & !< write in gnuplot format
  spl2d_writev, & !< write in paraview format
  spl2d_write, & !< write in data structure format
  spl2d_locpos, & !< find integer coordinates of point
  spl2d_deriv, & !< differentiate with respect to coordinate
  spl2d_coeff, & !< calculate coefficient array of splines
  spl2d_eval,  & !< evaluate splines
  spl2d_evaln   !< evaluate splines reusing coefficients

! public types
!> data structure for 2D splines
  type, public :: spl2d_t
     real(kr8), dimension(:,:), allocatable :: sampl !< values sampled on mesh
     real(kr8), dimension(:,:), allocatable :: coeff !< coefficients of 2-D splines
     real(kr8), dimension(:), allocatable :: pos1 !< padded list of positions, coordinate 1
     real(kr8), dimension(:), allocatable :: pos2 !< padded list of positions, coordinate 2
     real(kr8), dimension(:), allocatable :: knot1 !< padded list of knots, coordinate 1
     real(kr8), dimension(:), allocatable :: knot2 !< padded list of knots, coordinate 2
     real(kr8), dimension(:), allocatable :: wv1 !< 1-D spline work vector
     real(kr8), dimension(:), allocatable :: wv2 !< 1-D spline work vector
     real(kr8), dimension(:), allocatable :: wv3 !< 1-D spline work vector
     real(kr8), dimension(:,:), allocatable :: wa1k !< 1-D \f$ x \f$ (order) spline preserved work vector
     real(kr8), dimension(:,:), allocatable :: wa2k !< 1-D \f$ x \f$ (order) spline preserved work vector
     real(kr8), dimension(:), allocatable :: wv2k !< 1-D \f$ x \f$ (order) spline work vector
     real(kr8), dimension(:), allocatable :: val1 !< (order) spline work vector
     real(kr8), dimension(:), allocatable :: val2 !< (order) spline work vector
     integer(ki4), dimension(:), allocatable :: iwa1 !< 1-D spline preserved work vector
     integer(ki4), dimension(:), allocatable :: iwa2 !< 1-D spline preserved work vector
     integer(ki4) :: lunif !< if unity, knots are quasi-uniform
     real(kr8) :: h1 !< mesh spacing, coordinate 1
     real(kr8) :: h2 !< mesh spacing, coordinate 2
     real(kr8) :: rh1 !< reciprocal mesh spacing, coordinate 1
     real(kr8) :: rh2 !< reciprocal mesh spacing, coordinate 2
     integer(ki4) :: n1   !< unpadded array dimension in coordinate 1
     integer(ki4) :: n2   !< unpadded array dimension in coordinate 2
     integer(ki4) :: n1p   !< padded array dimension in coordinate 1
     integer(ki4) :: n2p   !< padded array dimension in coordinate 2
     integer(ki4) :: pad1=0   !< array padding in coordinate 1
     integer(ki4) :: pad2=0   !< array padding in coordinate 2
     integer(ki4) :: nord=4   !< order of splines
     integer(ki4) :: noff=2   !< fixup offset (might need to be four)
     real(kr8) :: org1 !< origin (minimum value) in coordinate 1
     real(kr8) :: org2 !< origin (minimum value) in coordinate 2
  end type spl2d_t

! public variables
! private types

! private variables
  character(*), parameter :: m_name='spl2d_m' !< module name
  real(kr8), dimension(:,:), allocatable :: work2 !< 2D work array
  integer   :: status   !< error status
  integer :: nin   !< input channel for spl2d data
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: i12m !< max 1D index
  integer(ki4) :: isw !< max 1D index
  integer(ki4) :: iswm !< max 1D index
  integer(ki4) :: idum !< dummy integer
  character(len=80) :: ibuff !< buffer for input/output
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> read in data format
subroutine spl2d_read(self,infile,kin)

  !! arguments
  type(spl2d_t), intent(out) :: self   !< object data structure
  character(len=80),intent(in) :: infile !< name of input file
  integer, intent(in),optional :: kin   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='spl2d_read' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  if(present(kin).AND.kin/=0) then
     !! assume unit already open and reading infile
     nin=kin
  else
     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

     !! open file
     call misc_getfileunit(nin)
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        call log_error(m_name,s_name,1,error_fatal,'Error opening data structure file')
     else
        call log_error(m_name,s_name,2,log_info,'data structure file opened')
     end if
  end if


  read(nin,*,iostat=status) ibuff
  !write(*,*) ibuff
  read(nin,*,iostat=status) self%lunif
  if(status/=0) then
     call log_error(m_name,s_name,7,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%h1
  if(status/=0) then
     call log_error(m_name,s_name,8,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%h2
  if(status/=0) then
     call log_error(m_name,s_name,9,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%rh1
  if(status/=0) then
     call log_error(m_name,s_name,10,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%rh2
  if(status/=0) then
     call log_error(m_name,s_name,11,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n1
  if(status/=0) then
     call log_error(m_name,s_name,12,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n2
  if(status/=0) then
     call log_error(m_name,s_name,13,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n1p
  if(status/=0) then
     call log_error(m_name,s_name,14,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n2p
  if(status/=0) then
     call log_error(m_name,s_name,15,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%pad1
  if(status/=0) then
     call log_error(m_name,s_name,16,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%pad2
  if(status/=0) then
     call log_error(m_name,s_name,17,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%nord
  if(status/=0) then
     call log_error(m_name,s_name,18,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%noff
  if(status/=0) then
     call log_error(m_name,s_name,19,error_fatal,'Error reading object data')
  end if

  !-----------------------------------------------------------------------
  !              Allocate 2D storage and read
  !! data array
  allocate(self%sampl(self%n1p,self%n2p), stat=status)
  call log_alloc_check(m_name,s_name,40,status)

  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%sampl
  if(status/=0) then
     call log_error(m_name,s_name,41,error_fatal,'Error reading object data')
  end if

  !! spline coefficients
  allocate(self%coeff(self%n1p,self%n2p), stat=status)
  call log_alloc_check(m_name,s_name,42,status)

  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%coeff
  if(status/=0) then
     call log_error(m_name,s_name,43,error_fatal,'Error reading object data')
  end if

  !-----------------------------------------------------------------------
  !              Allocate 1D storage and read
  ! position 1 array
  !! allocate position 1 storage
  allocate(self%pos1(self%n1p), stat=status)
  call log_alloc_check(m_name,s_name,50,status)

  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%pos1
  if(status/=0) then
     call log_error(m_name,s_name,51,error_fatal,'Error reading object data')
  end if
  ! position 2 array
  !! allocate position 2 storage
  allocate(self%pos2(self%n2p), stat=status)
  call log_alloc_check(m_name,s_name,52,status)

  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%pos2
  if(status/=0) then
     call log_error(m_name,s_name,53,error_fatal,'Error reading object data')
  end if
  ! knot 1 array
  !! allocate knot 1 storage
  allocate(self%knot1(self%n1p+self%nord), stat=status)
  call log_alloc_check(m_name,s_name,54,status)

  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%knot1
  if(status/=0) then
     call log_error(m_name,s_name,55,error_fatal,'Error reading object data')
  end if
  ! knot 2 array
  !! allocate knot 2 storage
  allocate(self%knot2(self%n2p+self%nord), stat=status)
  call log_alloc_check(m_name,s_name,56,status)

  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%knot2
  if(status/=0) then
     call log_error(m_name,s_name,57,error_fatal,'Error reading object data')
  end if

end subroutine spl2d_read
!---------------------------------------------------------------------
!> create spl2d data structure and knots
subroutine spl2d_init(self,pwork,kn1,kn2,porg1,porg2,ph1,ph2,korder)

  !! arguments
  type(spl2d_t), intent(inout) :: self   !< object data structure
  real(kr8), dimension(:,:), intent(in) :: pwork   !< input samples
  integer(ki4), intent(in) :: kn1   !< first nominal dimension
  integer(ki4), intent(in) :: kn2   !< second nominal dimension
  real(kr8), intent(in) :: porg1   !< minimum mesh value in coordinate 1
  real(kr8), intent(in) :: porg2   !< minimum mesh value in coordinate 2
  real(kr8), intent(in) :: ph1   !< mesh spacing, coordinate 1
  real(kr8), intent(in) :: ph2   !< mesh spacing, coordinate 2
  integer(ki4), intent(in) :: korder   !< spline order

  !! local
  character(*), parameter :: s_name='spl2d_init' !< subroutine name
  integer(ki4) :: i12k !< max 1D index times spline order
  !-----------------------------------------------------------------------
  !              Initialise scalar variables
  !self%pad1=0
  !self%pad2=0
  !self%nord=4

  self%lunif=2

  self%nord=korder
  self%noff=1+(korder-1)/2

  self%n1=kn1
  self%n2=kn2
  self%h1=ph1
  self%h2=ph2

  self%n1p=self%n1+2*self%pad1+1
  self%n2p=self%n2+2*self%pad2+1

  self%org1=porg1
  self%org2=porg2

  !-----------------------------------------------------------------------
  !            Input data array
  !! allocate 2D storage
  allocate(self%sampl(self%n1p,self%n2p), stat=status)
  call log_alloc_check(m_name,s_name,40,status)
  !-----------------------------------------------------------------------
  !              Initialise position arrays
  ! position 1 array
  !! allocate position 1 storage
  allocate(self%pos1(self%n1p), stat=status)
  call log_alloc_check(m_name,s_name,1,status)

  do i=1,self%n1p
     self%pos1(i)=(i-1-self%pad1)*self%h1+self%org1
  end do

  ! position 2 array
  !! allocate position 2 storage
  allocate(self%pos2(self%n2p), stat=status)
  call log_alloc_check(m_name,s_name,2,status)

  do j=1,self%n2p
     self%pos2(j)=(j-1-self%pad2)*self%h2+self%org2
  end do


  !-----------------------------------------------------------------------
  !              Initialise knots
  ! knot 1 array
  !! allocate knot 1 storage
  allocate(self%knot1(self%n1p+self%nord), stat=status)
  call log_alloc_check(m_name,s_name,3,status)

  do i=1,self%n1p-self%nord
     self%knot1(i+self%nord)=self%pos1(i+self%noff+self%pad1)
  end do

  do k=1,self%nord
     self%knot1(k)=self%pos1(1)
     self%knot1(self%n1p+k)=self%pos1(self%n1p)
  end do

  ! knot 2 array
  !! allocate knot 2 storage
  allocate(self%knot2(self%n2p+self%nord), stat=status)
  call log_alloc_check(m_name,s_name,4,status)

  do j=1,self%n2p-self%nord
     self%knot2(j+self%nord)=self%pos2(j+self%noff+self%pad2)
  end do

  do k=1,self%nord
     self%knot2(k)=self%pos2(1)
     self%knot2(self%n2p+k)=self%pos2(self%n2p)
  end do

  !-----------------------------------------------------------------------
  !              Other stuff
  self%rh1=1/self%h1
  self%rh2=1/self%h2


  !-----------------------------------------------------------------------
  !              Allocate 2D storage and work space
  !! allocate 2D storage
  allocate(self%coeff(self%n1p,self%n2p), stat=status)
  call log_alloc_check(m_name,s_name,10,status)

  i12m=max(self%n1p,self%n2p)
  !! allocate work storage
  allocate(self%wv1(i12m), stat=status)
  call log_alloc_check(m_name,s_name,11,status)
  self%wv1=0
  allocate(self%wv2(i12m), stat=status)
  call log_alloc_check(m_name,s_name,12,status)
  self%wv2=0
  allocate(self%iwa1(self%n1p), stat=status)
  call log_alloc_check(m_name,s_name,14,status)
  self%iwa1=0
  allocate(self%iwa2(self%n2p), stat=status)
  call log_alloc_check(m_name,s_name,15,status)
  self%iwa2=0
  allocate(self%wa1k(self%n1p,self%nord), stat=status)
  call log_alloc_check(m_name,s_name,16,status)
  allocate(self%wa2k(self%n2p,self%nord), stat=status)
  call log_alloc_check(m_name,s_name,17,status)
  iswm=2*i12m+3*self%nord+1
  allocate(self%wv3(iswm), stat=status)
  call log_alloc_check(m_name,s_name,18,status)
  i12k=i12m*self%nord
  allocate(self%wv2k(i12k), stat=status)
  call log_alloc_check(m_name,s_name,20,status)
  !-----------------------------------------------------------------------
  !              Initialise spline look up arrays
  ! value arrays
  !! allocate value 1 storage
  allocate(self%val1(self%nord), stat=status)
  call log_alloc_check(m_name,s_name,30,status)
  self%val1=0
  !! allocate value 2 storage
  allocate(self%val2(self%nord), stat=status)
  call log_alloc_check(m_name,s_name,31,status)
  !-----------------------------------------------------------------------
  !            Finally copy over input data and evaluate coefficients

  self%sampl=pwork

  call spl2d_coeff(self)

end subroutine spl2d_init
!---------------------------------------------------------------------
!> create dummy part of spl2d data structure
subroutine spl2d_initpart(self)

  !! arguments
  type(spl2d_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='spl2d_initpart' !< subroutine name
  integer(ki4) :: idum !< dummy allocation

  !-----------------------------------------------------------------------
  !              Initialise other scalar variables
  self%org1=self%pos1(1)
  self%org2=self%pos2(1)

  !-----------------------------------------------------------------------
  !              Allocate 2D storage and work space
  !! allocate work storage
  idum=1
  allocate(self%wv1(idum), stat=status)
  call log_alloc_check(m_name,s_name,11,status)
  self%wv1=0
  allocate(self%wv2(idum), stat=status)
  call log_alloc_check(m_name,s_name,12,status)
  self%wv2=0
  allocate(self%iwa1(idum), stat=status)
  call log_alloc_check(m_name,s_name,14,status)
  self%iwa1=0
  allocate(self%iwa2(idum), stat=status)
  call log_alloc_check(m_name,s_name,15,status)
  self%iwa2=0
  allocate(self%wa1k(idum,idum), stat=status)
  call log_alloc_check(m_name,s_name,16,status)
  allocate(self%wa2k(idum,idum), stat=status)
  call log_alloc_check(m_name,s_name,17,status)
  allocate(self%wv3(idum), stat=status)
  call log_alloc_check(m_name,s_name,18,status)
  allocate(self%wv2k(idum), stat=status)
  call log_alloc_check(m_name,s_name,20,status)
  !-----------------------------------------------------------------------
  !              Allocate value arrays
  !! allocate value 1 storage
  allocate(self%val1(self%nord), stat=status)
  call log_alloc_check(m_name,s_name,30,status)
  !! allocate value 2 storage
  allocate(self%val2(self%nord), stat=status)
  call log_alloc_check(m_name,s_name,31,status)

end subroutine spl2d_initpart
!---------------------------------------------------------------------
!> get \f$ \psi,\theta \f$ limits of interpolation
subroutine spl2d_ptlimits(self,p1min,p1max,p2min,p2max)

  !! arguments
  type(spl2d_t), intent(in) :: self   !< object data structure
  real(kr8), intent(out) :: p1min !< minimum value of position-1
  real(kr8), intent(out) :: p1max !< maximum value of position-1
  real(kr8), intent(out) :: p2min !< minimum value of position-2
  real(kr8), intent(out) :: p2max !< maximum value of position-2

  !! local
  character(*), parameter :: s_name='spl2d_ptlimits' !< subroutine name

  p1min=self%pos1(1)
  p1max=self%pos1(self%n1p)
  p2min=self%pos2(1)
  p2max=self%pos2(self%n2p)

end subroutine spl2d_ptlimits
!---------------------------------------------------------------------
!> scale points at which spline defined
subroutine spl2d_ptscale(self,qtfmdata)

  !! arguments
  type(spl2d_t), intent(inout) :: self   !< object data structure
  type(quantfm_t), intent(inout) :: qtfmdata   !< data defining transform

  !! local
  character(*), parameter :: s_name='spl2d_ptscale' !< subroutine name
  type(posveclis_t) :: zposl !< local variable
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zposq !< local variable
  integer(ki4) :: inqtfm   !< save transform type
  integer(ki4) :: inkn1   !< number of knots in 1-direction
  integer(ki4) :: inkn2   !< number of knots in 2-direction
  integer(ki4) :: inkn   !< maximum number of knots in either direction

  ! origins
  allocate(zposl%pos(1),stat=status)
  call log_alloc_check(m_name,s_name,10,status)
  zposl%pos(1)%posvec=(/self%org1, self%org2, 0._kr8/)
  zposl%np=1
  call position_qtfmlis(zposl,qtfmdata)
  self%org1=zposl%pos(1)%posvec(1)
  self%org2=zposl%pos(1)%posvec(2)
  deallocate(zposl%pos)

  ! spacings
  zpos%posvec=(/self%h1, self%h2 , 0._kr8/)
  inqtfm=qtfmdata%nqtfm
  ! case 1 just scales
  qtfmdata%nqtfm=1
  zposq=position_qtfm(zpos,qtfmdata)
  self%h1=zposq%posvec(1)
  self%rh1=1/zposq%posvec(1)
  self%rh2=1/zposq%posvec(2)
  self%h2=zposq%posvec(2)
  ! restore
  qtfmdata%nqtfm=inqtfm

  ! knots
  inkn1=self%n1p+self%nord
  inkn2=self%n2p+self%nord

  inkn=max(inkn1,inkn2)
  allocate(zposl%pos(inkn),stat=status)
  call log_alloc_check(m_name,s_name,20,status)
  do j=1,inkn
     zposl%pos(j)%posvec=0
  end do
  do j=1,inkn1
     zposl%pos(j)%posvec(1)=self%knot1(j)
  end do
  do j=1,inkn2
     zposl%pos(j)%posvec(2)=self%knot2(j)
  end do
  zposl%np=inkn
  call position_qtfmlis(zposl,qtfmdata)
  do j=1,inkn1
     self%knot1(j)=zposl%pos(j)%posvec(1)
  end do
  do j=1,inkn2
     self%knot2(j)=zposl%pos(j)%posvec(2)
  end do
  deallocate(zposl%pos)

end subroutine spl2d_ptscale
!---------------------------------------------------------------------
!> scale spline values
subroutine spl2d_scale(self,pfac)

  !! arguments
  type(spl2d_t), intent(inout) :: self   !< object data structure
  real(kr8), intent(in) :: pfac   !< scale factor

  !! local
  character(*), parameter :: s_name='spl2d_scale' !< subroutine name

  self%sampl=pfac*self%sampl
  self%coeff=pfac*self%coeff

end subroutine spl2d_scale
!---------------------------------------------------------------------
!> delete  spl2d data structure
subroutine spl2d_delete(self)

  !! arguments
  type(spl2d_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='spl2d_delete' !< subroutine name

  deallocate(self%sampl)
  deallocate(self%pos1)
  deallocate(self%pos2)
  deallocate(self%knot1)
  deallocate(self%knot2)
  deallocate(self%coeff)
  deallocate(self%wv1)
  deallocate(self%wv2)
  deallocate(self%iwa1)
  deallocate(self%iwa2)
  deallocate(self%wv3)
  deallocate(self%wa1k)
  deallocate(self%wa2k)
  deallocate(self%wv2k)
  deallocate(self%val1)
  deallocate(self%val2)

end subroutine spl2d_delete
!---------------------------------------------------------------------
!> write in gnuplot format
subroutine spl2d_writeg(self,kchar,kout)

  !! arguments
  type(spl2d_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: kchar  !< case
  integer, intent(in) :: kout   !< output channel for object data structure


  !! local
  character(*), parameter :: s_name='spl2d_writeg' !< subroutine name

  plot_type: select case (kchar)
  case('sampl')

     do j=1,self%n2p
        do i=1,self%n1p
           write(kout,'(2(1x,i4),'//cfmt2v,iostat=status) &
 &         i-1,j-1,self%pos1(i),self%pos2(j),self%sampl(i,j)
        end do
        write(kout,*,iostat=status) ' '
     end do
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Error writing sampl')
     end if

  case default

  end select plot_type

end subroutine spl2d_writeg
!---------------------------------------------------------------------
!> write in paraview format
subroutine spl2d_writev(self,kchar,kout)

  !! arguments
  type(spl2d_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: kchar  !< case
  integer, intent(in) :: kout   !< output channel for object data structure


  !! local
  character(*), parameter :: s_name='spl2d_writev' !< subroutine name

  plot_type: select case (kchar)
  case('sampl')

     do j=1,self%n2p
        do i=1,self%n1p
           write(kout,cfmtbs1,iostat=status) self%sampl(i,j)
        end do
     end do
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Error writing sampl')
     end if

  case default

  end select plot_type

end subroutine spl2d_writev
!---------------------------------------------------------------------
!> write in data format
subroutine spl2d_write(self,kout)

  !! arguments
  type(spl2d_t), intent(in) :: self   !< object data structure
  integer, intent(in) :: kout   !< output channel for object data structure


  !! local
  character(*), parameter :: s_name='spl2d_write' !< subroutine name

  write(kout,*,iostat=status) 'lunif'
  write(kout,*,iostat=status) self%lunif
  if(status/=0) then
     call log_error(m_name,s_name,7,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'h1'
  write(kout,*,iostat=status) self%h1
  if(status/=0) then
     call log_error(m_name,s_name,8,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'h2'
  write(kout,*,iostat=status) self%h2
  if(status/=0) then
     call log_error(m_name,s_name,9,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'rh1'
  write(kout,*,iostat=status) self%rh1
  if(status/=0) then
     call log_error(m_name,s_name,10,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'rh2'
  write(kout,*,iostat=status) self%rh2
  if(status/=0) then
     call log_error(m_name,s_name,11,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'n1'
  write(kout,*,iostat=status) self%n1
  if(status/=0) then
     call log_error(m_name,s_name,12,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'n2'
  write(kout,*,iostat=status) self%n2
  if(status/=0) then
     call log_error(m_name,s_name,13,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'n1p'
  write(kout,*,iostat=status) self%n1p
  if(status/=0) then
     call log_error(m_name,s_name,14,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'n2p'
  write(kout,*,iostat=status) self%n2p
  if(status/=0) then
     call log_error(m_name,s_name,15,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'pad1'
  write(kout,*,iostat=status) self%pad1
  if(status/=0) then
     call log_error(m_name,s_name,16,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'pad2'
  write(kout,*,iostat=status) self%pad2
  if(status/=0) then
     call log_error(m_name,s_name,17,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'nord'
  write(kout,*,iostat=status) self%nord
  if(status/=0) then
     call log_error(m_name,s_name,18,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'noff'
  write(kout,*,iostat=status) self%noff
  if(status/=0) then
     call log_error(m_name,s_name,19,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'sampl'
  write(kout,*,iostat=status) self%sampl
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'coeff'
  write(kout,*,iostat=status) self%coeff
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'pos1'
  write(kout,*,iostat=status) self%pos1
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'pos2'
  write(kout,*,iostat=status) self%pos2
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'knot1'
  write(kout,*,iostat=status) self%knot1
  if(status/=0) then
     call log_error(m_name,s_name,5,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'knot2'
  write(kout,*,iostat=status) self%knot2
  if(status/=0) then
     call log_error(m_name,s_name,6,error_fatal,'Error writing object data')
  end if

end subroutine spl2d_write
!---------------------------------------------------------------------
!> calculate coefficient array of splines
subroutine spl2d_coeff(self)

  !! arguments
  type(spl2d_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='spl2d_coeff' !< subroutine name

  !-----------------------------------------------------------------------
  !              spline  coefficients via direct product
  !       note that arrays iwa,wa aka IL,AN are saved from one call to the next,
  !       only recalculated when last calling argument changes
  !
  !! first work in 1 coordinate
  isw=2*self%n1p+3*self%nord+1
  do j=1,self%n2p
     self%wv1(:self%n1p)=self%sampl(:,j)
     call tb06a(self%n1p,self%pos1,self%wv1,self%nord,self%n1p, &
 &   self%knot1(self%nord+1),self%iwa1,self%wa1k,self%wv2k,isw,self%wv3, &
 &   self%wv2,1)
     self%coeff(:self%n1p,j)=self%wv2(:self%n1p)
  end do

  !!  then in 2 coordinate
  isw=2*self%n2p+3*self%nord+1
  do i=1,self%n1p
     self%wv1(:self%n2p)=self%coeff(i,:)
     call tb06a(self%n2p,self%pos2,self%wv1,self%nord,self%n2p, &
 &   self%knot2(self%nord+1),self%iwa2,self%wa2k,self%wv2k,isw,self%wv3, &
 &   self%wv2,2)
     self%coeff(i,:self%n2p)=self%wv2(:self%n2p)
  end do

end subroutine spl2d_coeff
!---------------------------------------------------------------------
!> evaluate spline
subroutine spl2d_eval(self,p1,p2,pe)

  !! arguments
  type(spl2d_t), intent(inout) :: self   !< object data structure
  real(kr8), intent(in) :: p1   !< first coordinate of point
  real(kr8), intent(in) :: p2   !< second coordinate of point
  real(kr8), intent(out) :: pe   !< value of spline at (p1,p2)

  !! local
  character(*), parameter :: s_name='spl2d_eval' !< subroutine name
  integer(ki4), save :: ipp=0   !<  first integer coordinate of point (poss. offset)
  integer(ki4), save :: iqq=0   !<  second integer coordinate of point (poss. offset)
  integer(ki4) :: ip   !<  first integer coordinate of point
  integer(ki4) :: iq   !<  second integer coordinate of point
  integer(ki4) :: ii   !<  integer loop work
  integer(ki4) :: ij   !<  integer loop work
  integer(ki4) :: iflag   !<  warning flag
  real(kr8) :: zz    !< local variable
  real(kr8) :: p1f   !< fractional part of p1
  real(kr8) :: p2f   !< fractional part of p2
  !dbgwval1  real(kr8), dimension(4) :: zwork   !< workspace !dbgwval1

  zz=0
  if (self%lunif==0) then
     ! assumes points non-uniformly distributed
     call interv(self%knot1,self%n1p+self%nord,p1,ipp,iflag)
     if (iflag>0) then
        call log_error(m_name,s_name,1,error_warning,'Point not in range of spline')
        iflag=0
     end if
     call interv(self%knot2,self%n2p+self%nord,p2,iqq,iflag)
     if (iflag>0) then
        call log_error(m_name,s_name,2,error_warning,'Point not in range of spline')
        iflag=0
     end if
     ip=ipp-self%nord
     iq=iqq-self%nord
     call bsplvn(self%knot1,self%nord,1,p1,ipp,self%val1)
     do i=1,self%nord
        ii=i+iq
        self%val2(i)=0
        do j=1,self%nord
           self%val2(i)=self%val2(i)+self%val1(j)*self%coeff(j+ip,ii)
        end do
     end do

     call bsplvn(self%knot2,self%nord,1,p2,iqq,self%val1)
     do i=1,self%nord
        zz=zz+self%val2(i)*self%val1(i)
     end do

  else if (self%lunif==1) then
     ! assume uniform distribution for point location, not for final evaluation
     call interu(self%rh1,self%n1p,p1-self%org1,ipp,self%pad1,self%noff,self%nord)
     call interu(self%rh2,self%n2p,p2-self%org2,iqq,self%pad2,self%noff,self%nord)
     ip=ipp-self%nord
     iq=iqq-self%nord
     call bsplvn(self%knot1,self%nord,1,p1,ipp,self%val1)
     !dbgwval1     write(121,*) 'val1 array 1',ip,(self%val1(i),i=1,4) !dbgwval1
     do i=1,self%nord
        ii=i+iq
        self%val2(i)=0
        do j=1,self%nord
           self%val2(i)=self%val2(i)+self%val1(j)*self%coeff(j+ip,ii)
        end do
     end do

     call bsplvn(self%knot2,self%nord,1,p2,iqq,self%val1)
     !dbgwval1     write(121,*) 'val1 array 2',iq,(self%val1(i),i=1,4) !dbgwval1
     do i=1,self%nord
        zz=zz+self%val2(i)*self%val1(i)
     end do

  else if (self%lunif==2) then
     ! assume uniform distribution for both point location and final evaluation
     call interw(self%rh1,self%n1p,p1-self%org1,self%pad1,self%noff,self%nord,ipp,p1f)
     call interw(self%rh2,self%n2p,p2-self%org2,self%pad2,self%noff,self%nord,iqq,p2f)
     ip=ipp-self%nord
     iq=iqq-self%nord
     call bsplwn(ipp,p1f,self%n1p,self%val1)
     !dbgwval1     zwork(1:4)=self%val1 !dbgwval1
     !dbgwval1     write(122,*) 'val1 array 1',ip,(zwork(i),i=1,4) !dbgwval1
     !optw     ii=iq
     !optw     do i=1,self%nord
     !optw        ii=ii+1
     !optw        self%val2(i)=0
     !optw        ij=ip
     !optw        do j=1,self%nord
     !optw           ij=ij+1
     !optw           self%val2(i)=self%val2(i)+self%val1(j)*self%coeff(ij,ii)
     !optw        end do
     !optw     end do
     self%val2=matmul(self%val1,self%coeff(ip+1:ipp,iq+1:iqq))

     call bsplwn(iqq,p2f,self%n2p,self%val1)
     !dbgwval1     zwork(1:4)=self%val1 !dbgwval1
     !dbgwval1     write(122,*) 'val1 array 2',iq,(zwork(i),i=1,4) !dbgwval1
     !optw     do i=1,self%nord
     !optw        zz=zz+self%val2(i)*self%val1(i)
     !optw     end do
     zz=dot_product(self%val2,self%val1)
  end if

  pe=zz

end subroutine spl2d_eval
!---------------------------------------------------------------------
!> evaluate spline reusing coefficients.
subroutine spl2d_evaln(self,p1,p2,kcall,pe)

  !! arguments
  type(spl2d_t), intent(inout) :: self   !< object data structure
  real(kr8), intent(in) :: p1   !< first coordinate of point
  real(kr8), intent(in) :: p2   !< second coordinate of point
  integer(ki4), intent(in) :: kcall   !< first or second call
  real(kr8), intent(out) :: pe   !< value of spline at (p1,p2)

  !! local
  character(*), parameter :: s_name='spl2d_evaln' !< subroutine name
  integer(ki4), save :: ipp=0   !<  first integer coordinate of point (poss. offset)
  integer(ki4), save :: iqq=0   !<  second integer coordinate of point (poss. offset)
  integer(ki4) :: ip   !<  first integer coordinate of point
  integer(ki4) :: iq   !<  second integer coordinate of point
  integer(ki4) :: ii   !<  integer loop work
  integer(ki4) :: ij   !<  integer loop work
  integer(ki4) :: iflag   !<  warning flag
  real(kr8) :: zz    !< local variable
  real(kr8) :: p1f   !< fractional part of p1
  real(kr8) :: p2f   !< fractional part of p2
  !dbgwval1  real(kr8), dimension(4) :: zwork   !< workspace !dbgwval1
  integer(ki4), save :: sip   !< reusable index
  integer(ki4), save :: sipp  !< reusable index
  integer(ki4), save :: siq   !< reusable index
  integer(ki4), save :: siqq  !< reusable index
  real(kr8), save, dimension(4) :: scoef1   !< reusable coefficients, NB only cubic case
  real(kr8), save, dimension(4) :: scoef2   !< reusable coefficients

  zz=0
  if (self%lunif==0) then
     ! assumes points non-uniformly distributed
     call interv(self%knot1,self%n1p+self%nord,p1,ipp,iflag)
     if (iflag>0) then
        call log_error(m_name,s_name,1,error_warning,'Point not in range of spline')
        iflag=0
     end if
     call interv(self%knot2,self%n2p+self%nord,p2,iqq,iflag)
     if (iflag>0) then
        call log_error(m_name,s_name,2,error_warning,'Point not in range of spline')
        iflag=0
     end if
     ip=ipp-self%nord
     iq=iqq-self%nord
     call bsplvn(self%knot1,self%nord,1,p1,ipp,self%val1)
     do i=1,self%nord
        ii=i+iq
        self%val2(i)=0
        do j=1,self%nord
           self%val2(i)=self%val2(i)+self%val1(j)*self%coeff(j+ip,ii)
        end do
     end do

     call bsplvn(self%knot2,self%nord,1,p2,iqq,self%val1)
     do i=1,self%nord
        zz=zz+self%val2(i)*self%val1(i)
     end do

  else if (self%lunif==1) then
     ! assume uniform distribution for point location, not for final evaluation
     call interu(self%rh1,self%n1p,p1-self%org1,ipp,self%pad1,self%noff,self%nord)
     call interu(self%rh2,self%n2p,p2-self%org2,iqq,self%pad2,self%noff,self%nord)
     ip=ipp-self%nord
     iq=iqq-self%nord
     call bsplvn(self%knot1,self%nord,1,p1,ipp,self%val1)
     !dbgwval1     write(121,*) 'val1 array 1',ip,(self%val1(i),i=1,4) !dbgwval1
     do i=1,self%nord
        ii=i+iq
        self%val2(i)=0
        do j=1,self%nord
           self%val2(i)=self%val2(i)+self%val1(j)*self%coeff(j+ip,ii)
        end do
     end do

     call bsplvn(self%knot2,self%nord,1,p2,iqq,self%val1)
     !dbgwval1     write(121,*) 'val1 array 2',iq,(self%val1(i),i=1,4) !dbgwval1
     do i=1,self%nord
        zz=zz+self%val2(i)*self%val1(i)
     end do

  else if (self%lunif==2) then
     ! assume uniform distribution for both point location and final evaluation
     if (kcall==1) then
        call interw(self%rh1,self%n1p,p1-self%org1,self%pad1,self%noff,self%nord,ipp,p1f)
        call interw(self%rh2,self%n2p,p2-self%org2,self%pad2,self%noff,self%nord,iqq,p2f)
        ip=ipp-self%nord
        iq=iqq-self%nord
        call bsplwn(ipp,p1f,self%n1p,scoef1)
        call bsplwn(iqq,p2f,self%n2p,scoef2)
        sip=ip;sipp=ipp;siq=iq;siqq=iqq
     end if
     self%val2=matmul(scoef1,self%coeff(sip+1:sipp,siq+1:siqq))
     zz=dot_product(self%val2,scoef2)
  end if

  pe=zz

end subroutine spl2d_evaln
!---------------------------------------------------------------------
!> find integer coordinates of point
subroutine spl2d_locpos(self,p1,p2,kp1,kp2)

  !! arguments
  type(spl2d_t), intent(inout) :: self   !< object data structure
  real(kr8), intent(in) :: p1   !< first coordinate of point
  real(kr8), intent(in) :: p2   !< second coordinate of point
  integer(ki4), intent(out) :: kp1   !<  first integer coordinate of point
  integer(ki4), intent(out) :: kp2   !<  second integer coordinate of point

  !! local
  character(*), parameter :: s_name='spl2d_locpos' !< subroutine name
  integer(ki4) :: ipp   !<  first integer coordinate of point (poss. offset)
  integer(ki4) :: iqq   !<  second integer coordinate of point (poss. offset)
  integer(ki4) :: iflag   !<  warning flag

  if (self%lunif==0) then
     call interv(self%knot1,self%n1p+self%nord,p1,ipp,iflag)
     if (iflag>0) then
        call log_error(m_name,s_name,1,error_warning,'Point not in range of spline')
        iflag=0
     end if
     call interv(self%knot2,self%n2p+self%nord,p2,iqq,iflag)
     if (iflag>0) then
        call log_error(m_name,s_name,2,error_warning,'Point not in range of spline')
        iflag=0
     end if
     kp1=ipp-self%nord
     kp2=iqq-self%nord
  else if (self%lunif>=1) then
     call interu(self%rh1,self%n1p,p1-self%org1,ipp,0,0,0)
     call interu(self%rh2,self%n2p,p2-self%org2,iqq,0,0,0)
     if (ipp*iqq==0) then
        call log_error(m_name,s_name,3,error_warning,'Point not in range of spline')
        ipp=max(1,ipp)
        iqq=max(1,iqq)
     end if
     kp1=ipp
     kp2=iqq
  end if

end subroutine spl2d_locpos
!---------------------------------------------------------------------
!> differentiate with respect to coordinate
subroutine spl2d_deriv(self,dself,kcoo)

  !! arguments
  type(spl2d_t), intent(inout) :: self   !< object data structure
  type(spl2d_t), intent(out) :: dself   !< derivative object data structure
  integer(ki4), intent(in) :: kcoo !< coordinate direction

  !! local
  character(*), parameter :: s_name='spl2d_deriv' !< subroutine name
  integer(ki4) :: iflag   !<  error flag
  real(kr8) :: z   !< point to evaluate derivative
  real(kr8) :: zd   !< derivative value

  allocate(work2(self%n1p,self%n2p), stat=status)
  call log_alloc_check(m_name,s_name,1,status)

  if (kcoo==1) then
     ! derivative in 1-direction
     do j=1,self%n2p
        self%wv1(:self%n1p)=self%coeff(:,j)
        do i=1,self%n1p
           z=(i-1)*self%h1+self%org1
           call bvalue(self%knot1,self%wv1,self%n1p,self%nord,z,1,self%lunif,iflag,zd)
           if (iflag>0) then
              call log_error(m_name,s_name,10,error_fatal,'Point not in range of spline')
              iflag=0
           end if
           self%wv2(i)=zd
        end do
        work2(:self%n1p,j)=self%wv2(:self%n1p)
     end do
  else

     ! derivative in 2-direction
     do i=1,self%n1p
        self%wv1(:self%n2p)=self%coeff(i,:)
        do j=1,self%n2p
           z=(j-1)*self%h2+self%org2
           call bvalue(self%knot2,self%wv1,self%n2p,self%nord,z,1,self%lunif,iflag,zd)
           if (iflag>0) then
              call log_error(m_name,s_name,20,error_fatal,'Point not in range of spline')
              iflag=0
           end if
           self%wv2(j)=zd
        end do
        work2(i,:self%n2p)=self%wv2(:self%n2p)
     end do
  end if

  call spl2d_init(dself,work2,self%n1,self%n2,&
 &self%org1,self%org2,self%h1,self%h2,self%nord)
  deallocate(work2)

end subroutine spl2d_deriv

end module spl2d_m
