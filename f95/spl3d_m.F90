module spl3d_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use misc_m
  use position_h
  use position_m
  use spl2d_m

  implicit none
  private

! public subroutines
  public :: spl3d_read, & !< read in data structure format
  spl3d_init,   & !< create spl3d data structure
  spl3d_mask,   & !< calculate spline masking
  spl3d_delete,   & !< delete  spl3d data structure
  spl3d_writeg,   &  !< write in gnuplot structure format
  spl3d_write,  &  !< write in data structure format
  spl3d_initpart,   & !< create create dummy part of spl3d data structure
  spl3d_ptlimits,   & !< get \f$ \psi,\theta \f$ limits of interpolation
  spl3d_ptscale,  &  !< scale points at which spline defined
  spl3d_scale,  &  !< scale spline values
  spl3d_deriv, & !< differentiate with respect to coordinate
  spl3d_coeff, & !< calculate coefficient array of splines
  spl3d_testtfm, & !< test transform by inverting
  spl3d_eval, &  !< evaluate splines
  spl3d_evalm   !< evaluate splines using mask

! public types
!> data structure for 3D splines
  type, public :: spl3d_t
     real(kr8), dimension(:), allocatable :: pos3 !< list of positions, coordinate 3
     real(kr8), dimension(:), allocatable :: posk !< k-space list of positions, coordinate 3
     real(kr8), dimension(:), allocatable :: poskof !< k times offset, coordinate 3
     real(kr8) :: h3 !< mesh spacing, coordinate 3 (not used)
     real(kr8) :: rh3 !< reciprocal mesh spacing, coordinate 3  (not used)
     integer(ki4) :: n3   !< array dimension in coordinate 3 (=n0*nk)
     integer(ki4) :: n0   !< first dimension in 3-D spline definition (usually 2 for cos/sin)
     integer(ki4) :: ncpt   !< number of components in spline definition
     integer(ki4) :: nk   !< number of wave numbers in spline definition
     integer(ki4), dimension(3) :: kmin   !< minimum wave number in spline evaluation
     integer(ki4), dimension(3) :: kmax   !< maximum wave number in spline evaluation
     integer(ki4), dimension(3) :: parity !< parity of component (1=even,2=odd,3=both)
     integer(ki4) :: nstate   !< transformed state of spline definition
     type(spl2d_t), dimension(:,:,:), allocatable :: farm !< 3-D vector spline field
     integer(ki4), dimension(:,:,:,:), allocatable :: mask !< sets mode number limits for parity and component
  end type spl3d_t

! public variables
  integer(ki4), public, parameter :: spl3d_rpower=1 !< power of \f$ R \f$ multiplying \f$ \bf{B} \f$

! private types

! private variables
  character(*), parameter :: m_name='spl3d_m' !< module name
  real(kr8), dimension(:,:), allocatable :: work2 !< 2D work array
  real(kr8), dimension(:,:,:), allocatable :: work3 !< 3D work array
  real(kr8), dimension(:,:,:,:), allocatable :: work4 !< 4D work array
  integer   :: status   !< error status
  integer :: nin   !< input channel for spl3d data
  integer(ki4) :: i !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: jp !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: idum !< dummy integer
  character(len=80) :: ibuff !< buffer for input/output
  logical :: iltest !< logical flag
  logical, parameter :: lplog2d=.false. !< logical flag for logging by 2d
  logical, parameter :: lpddtest=.false. !< logical flag for testing dodgy data
  logical, parameter :: lpmask=.true. !< logical flag for masking
  logical, parameter :: lpmask3=.false. !< logical flag for masking using 3-cpt values only
  logical, parameter :: lpmaskr=.true. !< logical flag for masking using reference field value
  integer(ki4), parameter :: ip1diag=75 !< diagnostic in 1 coordinate
  integer(ki4), parameter :: ip2diag=8 !< diagnostic in 2 coordinate
  integer(ki4), parameter :: ip3diag=32 !< diagnostic in 3 coordinate
  integer(ki4), parameter :: ipcpt=3 !< diagnostic in component
  integer(ki4), parameter :: ipcpt2=2 !< diagnostic 2 in component

  contains
!---------------------------------------------------------------------
!> read in data format
subroutine spl3d_read(self,infile,kin)

  !! arguments
  type(spl3d_t), intent(out) :: self   !< object data structure
  character(len=80),intent(in) :: infile !< name of input file
  integer, intent(in),optional :: kin   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='spl3d_read' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer(ki4) :: in1   !< first nominal dimension
  integer(ki4) :: in2   !< second nominal dimension

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
  call log_read_check(m_name,s_name,10,status)
  read(nin,*,iostat=status) self%n0
  call log_read_check(m_name,s_name,11,status)
  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,12,status)
  read(nin,*,iostat=status) self%ncpt
  call log_read_check(m_name,s_name,13,status)
  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,14,status)
  read(nin,*,iostat=status) self%nk
  call log_read_check(m_name,s_name,15,status)
  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,16,status)
  read(nin,*,iostat=status) self%nstate
  call log_read_check(m_name,s_name,17,status)
  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,18,status)
  read(nin,*,iostat=status) self%n3
  call log_read_check(m_name,s_name,19,status)
  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,20,status)
  read(nin,*,iostat=status) self%h3
  call log_read_check(m_name,s_name,21,status)
  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,22,status)
  read(nin,*,iostat=status) self%rh3
  call log_read_check(m_name,s_name,23,status)
  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,24,status)
  read(nin,*,iostat=status) self%kmin
  call log_read_check(m_name,s_name,25,status)
  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,26,status)
  read(nin,*,iostat=status) self%kmax
  call log_read_check(m_name,s_name,27,status)
  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,28,status)
  read(nin,*,iostat=status) self%parity
  call log_read_check(m_name,s_name,29,status)

  !! check for valid data
  if(self%n0<=0.OR.self%n0>=5) &
 &call log_error(m_name,s_name,30,error_fatal,'n0 must be small positive integer')
  if(self%ncpt<=0.OR.self%ncpt>=4) &
 &call log_error(m_name,s_name,31,error_fatal,'ncpt must be small positive integer')
  if(self%nk<=0) &
 &call log_error(m_name,s_name,32,error_fatal,'nk must be positive integer')
  if(self%nstate<=0.OR.self%nstate>3) &
 &call log_error(m_name,s_name,33,error_fatal,'nstate must be small positive integer')

  !-----------------------------------------------------------------------
  !              Allocate 1D storage and read
  ! position 3 array
  !! allocate position 3 storage
  allocate(self%pos3(self%n3), stat=status)
  call log_alloc_check(m_name,s_name,50,status)

  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,51,status)
  read(nin,*,iostat=status) self%pos3
  if(status/=0) then
     call log_error(m_name,s_name,52,error_fatal,'Error reading object data')
  end if
  ! position k array
  !! allocate position k storage
  allocate(self%posk(0:self%nk), stat=status)
  call log_alloc_check(m_name,s_name,53,status)

  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,54,status)
  read(nin,*,iostat=status) self%posk
  if(status/=0) then
     call log_error(m_name,s_name,55,error_fatal,'Error reading object data')
  end if
  ! k offset array
  !! allocate k offset storage
  allocate(self%poskof(0:self%nk), stat=status)
  call log_alloc_check(m_name,s_name,56,status)

  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,57,status)
  read(nin,*,iostat=status) self%poskof
  if(status/=0) then
     call log_error(m_name,s_name,58,error_fatal,'Error reading object data')
  end if

  !-----------------------------------------------------------------------
  !              Allocate 3D storage and read
  !! data array
  allocate(self%farm(self%ncpt,self%n0,0:self%nk), stat=status)
  call log_alloc_check(m_name,s_name,60,status)

  ij=0
  do l=0,self%nk
     do k=1,self%n0
        do j=1,self%ncpt
           call spl2d_read( self%farm(j,k,l),infile,nin )
           if (lplog2d) print '("2D spline read index = ",3(i9,1x))',j,k,l
           ij=ij+1
           if (lplog2d) call log_value("number of 2D splines read ",ij)
           call spl2d_initpart( self%farm(j,k,l) )
        end do
     end do
  end do

  !-----------------------------------------------------------------------
  !              Allocate 4D storage and read

  in1=self%farm(1,1,0)%n1p
  in2=self%farm(1,1,0)%n2p
  !! mask array
  allocate(self%mask(self%ncpt,self%n0,in1,in2), stat=status)
  call log_alloc_check(m_name,s_name,70,status)

  read(nin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,80,status)
  read(nin,*,iostat=status) self%mask
  if(status/=0) then
     call log_error(m_name,s_name,81,error_fatal,'Error reading mask data')
  end if

  print '("finished reading 3D spline")'
  call log_error(m_name,s_name,99,log_info,'finished reading 3D spline')

end subroutine spl3d_read
!---------------------------------------------------------------------
!> create spl3d data structure
subroutine spl3d_init(self,pwork,kcpt,kn0,kn1,kn2,kn3,pos1,pos2,pos3,korder,kmin,kmax)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure
  real(kr8), dimension(:,:,:,:), intent(in) :: pwork   !< input samples
  integer(ki4), intent(in) :: kcpt   !< number of vector components
  integer(ki4), intent(in) :: kn0   !< zero nominal dimension
  integer(ki4), intent(in) :: kn1   !< first nominal dimension
  integer(ki4), intent(in) :: kn2   !< second nominal dimension
  integer(ki4), intent(in) :: kn3   !< third nominal dimension
  real(kr8), dimension(:), intent(in) :: pos1   !< mesh array, coordinate 1
  real(kr8), dimension(:), intent(in) :: pos2   !< mesh array, coordinate 2
  real(kr8), dimension(:), intent(in) :: pos3   !< mesh array, coordinate 3
  integer(ki4), intent(in) :: korder   !< spline order
  integer(ki4), dimension(3), intent(in) :: kmin   !< min wavenumber to be used in evaluation
  integer(ki4), dimension(3), intent(in) :: kmax   !< max wavenumber to be used in evaluation

  !! local
  character(*), parameter :: s_name='spl3d_init' !< subroutine name
  integer(ki4) :: in1   !< first nominal dimension
  integer(ki4) :: in2   !< second nominal dimension
  integer(ki4) :: ik   !< third nominal dimension
  real(kr8) :: zorg1   !< minimum mesh value in coordinate 1
  real(kr8) :: zorg2   !< minimum mesh value in coordinate 2
  real(kr8) :: zh1   !< mesh spacing, coordinate 1
  real(kr8) :: zh2   !< mesh spacing, coordinate 2
  real(kr8) :: zh3   !< mesh spacing, coordinate 3
  real(kr8) :: znorm   !< normalisation in wavenumber
  real(kr8) :: zkoff   !< offset in wavenumber
  integer(ki4) :: il   !< counter
  real(kr8) :: znk   !< maximum in coordinate 3
  !-----------------------------------------------------------------------
  !              Initialise scalar variables
  !self%pad1=0
  !self%pad2=0
  !self%nord=4

  self%n0=kn0
  self%n3=kn3
  self%ncpt=kcpt

  in1=kn1
  in2=kn2
  ik=kn3/kn0
  znk=max(1,ik)
  self%nk=ik
  self%kmin=kmin
  self%kmax=kmax

  zh1=pos1(2)-pos1(1)
  zh2=pos2(2)-pos2(1)
  zh3=pos3(2)-pos3(1)
  self%h3=zh3

  zorg1=pos1(1)
  zorg2=pos2(1)

  !-----------------------------------------------------------------------
  !              Initialise position arrays
  ! position 3 array
  !! allocate position 3 storage
  allocate(self%pos3(self%n3), stat=status)
  call log_alloc_check(m_name,s_name,1,status)

  do i=1,self%n3
     self%pos3(i)=pos3(i)
  end do

  ! position k array (for k space values of in 3 coordinate)
  !! allocate position k storage
  allocate(self%posk(0:ik), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  !! allocate k offset storage
  allocate(self%poskof(0:ik), stat=status)
  call log_alloc_check(m_name,s_name,3,status)

  if (lpddtest) &
 &call log_error(m_name,s_name,4,log_info,'Assuming equil. field duplicated in data')
  if (lpmask) then
     call log_error(m_name,s_name,5,log_info,'Masking in testtfm')
     if (lpmask3) &
 &   call log_error(m_name,s_name,6,log_info,'Masking using field 3-component even')
     if (lpmaskr) &
 &   call log_error(m_name,s_name,7,log_info,'Masking using reference amplitude')
  end if
  znorm=2*const_pid/(self%pos3(self%n3)-self%pos3(1))
  zkoff=0
  !      else
  !! these are to be scaled later
  !         znorm=1
  !         zkoff=0
  !      end if

  do l=0,ik
     self%posk(l)=l*znorm
     self%poskof(l)=l*zkoff
  end do

  !-----------------------------------------------------------------------
  !            Input data array
  !! allocate 3D storage
  if (kn0*kcpt*ik>0) then
     allocate(self%farm(kcpt,kn0,0:ik), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
  else
     call log_error(m_name,s_name,11,error_fatal,'no data')
  end if
  !! allocate 2D storage
  if (in1*in2>0) then
     allocate(work2(in1,in2), stat=status)
     call log_alloc_check(m_name,s_name,12,status)
  else
     call log_error(m_name,s_name,13,error_fatal,'no data')
  end if

  !! define
  ij=0
  il=0
  do l=0,ik
     do k=1,kn0
        il=il+1
        do j=1,kcpt
           il=min(il,self%n3)
           work2=pwork(il,j,1:in1,1:in2)
           call  spl2d_init(self%farm(j,k,l),work2, &
 &         in1-1,in2-1,zorg1,zorg2,zh1,zh2,korder)
           if (lplog2d) print '("2D spline created index = ",3(i9,1x))',j,k,l
           ij=ij+1
           if (lplog2d) call log_value("created 2D spline number ",ij)
        end do
     end do
  end do

  !              Other stuff
  self%rh3=1/self%h3
  self%nstate=1
  deallocate(work2)

  call log_error(m_name,s_name,99,log_info,'finished initialising 3D spline')

end subroutine spl3d_init
!---------------------------------------------------------------------
!> calculate spline masking
subroutine spl3d_mask(self,cutoff,parity,refamp)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure
  real(kr8), intent(in) :: cutoff !<  for maximum mode number calculation
  integer(ki4), dimension(3), intent(in) :: parity !< parity of component (1=even,2=odd,3=both)
  real(kr8), intent(in) :: refamp !<  reference amplitude for cutoff

  !! local
  character(*), parameter :: s_name='spl3d_mask' !< subroutine name
  real(kr8) :: za   !< a part of complex mode amplitude
  real(kr8) :: zampt  !< sum of all amplitudes of complex mode
  integer(ki4) :: iorder   !< spline order
  integer(ki4) :: in1   !< first dimension
  integer(ki4) :: in2   !< second dimension
  integer(ki4) :: j1   !< loop counter for first dimension
  integer(ki4) :: j2   !< loop counter for second dimension
  integer(ki4) :: jj1   !< second loop counter for first dimension
  integer(ki4) :: jj2   !< second loop counter for second dimension
  integer(ki4) :: ij1min    !< local variable
  integer(ki4) :: ij2min    !< local variable
  integer(ki4), dimension(:,:), allocatable :: iwork2 !< 2D work array
  integer(ki4) :: ikmin   !< min wavenumber to be used in evaluation
  integer(ki4) :: lmax !< maximum mode number for point
  real(kr8) :: zampc  !< critical amplitude for mode

  ! record parity
  self%parity=parity
  ! record mask set
  self%nstate=3

  iorder=self%farm(1,1,0)%nord
  !! allocate storage
  in1=self%farm(1,1,0)%n1p
  in2=self%farm(1,1,0)%n2p
  if (in1*in2>0) then
     allocate(iwork2(in1,in2), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,3,error_fatal,'no data')
  end if
  !! mask
  allocate(self%mask(self%ncpt,self%n0,in1,in2), stat=status)
  call log_alloc_check(m_name,s_name,4,status)
  self%mask=0

  ! analyse 1-D k dependence for all R,Z
  do j=1,self%ncpt

     k=min(2,self%parity(j))
     looparity: do jp=1,max(1,self%parity(j)-1)
        ikmin=k-1

        loopRZ: do j2=1,in2
           do j1=1,in1

              ! first total amplitude
              zampt=0
              do l=ikmin,self%nk
                 za=self%farm(j,k,l)%sampl(j1,j2)
                 zampt=zampt+abs(za)
              end do
              ! now find largest k exceeding fraction given by cut-off
              iwork2(j1,j2)=0
              if (zampt>=cutoff) then
                 ! ignore very small total amps
                 if (lpmaskr) then
                    ! critical amplitude based on reference
                    zampc=cutoff*refamp
                 else
                    ! critical amplitude is relative
                    zampc=cutoff*zampt
                 end if
                 do l=self%nk,ikmin,-1
                    za=self%farm(j,k,l)%sampl(j1,j2)
                    if (abs(za)>zampc) then
                       iwork2(j1,j2)=l
                       exit
                    end if
                 end do
              end if

           end do
        end do loopRZ

        ! now maximum k over spline evaluation
        loopRZ2: do j2=1,in2
           ij2min=max(1,j2-iorder)
           do j1=1,in1
              ij1min=max(1,j1-iorder)

              lmax=0
              loopspline: do jj2=ij2min,j2
                 do jj1=ij1min,j1
                    lmax=max(lmax,iwork2(jj1,jj2))
                 end do
              end do loopspline
              ! assign to mask
              ! note lmax for (j1-1,j2-1) saved in (j1,j2)
              self%mask(j,k,j1,j2)=lmax

              !     write(*,*) 'mask',j,k,j1,j2, self%mask(j,k,j1,j2)

           end do
        end do loopRZ2

        k=k-1
     end do looparity

  end do

  deallocate(iwork2)

end subroutine spl3d_mask
!---------------------------------------------------------------------
!> delete  spl3d data structure
subroutine spl3d_delete(self)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='spl3d_delete' !< subroutine name

  ij=0
  do l=0,self%nk
     do k=1,self%n0
        do j=1,self%ncpt
           call  spl2d_delete(self%farm(j,k,l))
           ij=ij+1
           if (lplog2d) print '("deleted 2D spline number = ",i10)',ij
           if (lplog2d) call log_value("deleted 2D spline number ",ij)
        end do
     end do
  end do

  call log_error(m_name,s_name,99,log_info,'finished deleting 3D spline')

end subroutine spl3d_delete
!---------------------------------------------------------------------
!> write in gnuplot format
subroutine spl3d_writeg(self,kchar,kout)

  !! arguments
  type(spl3d_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: kchar  !< case
  integer, intent(in) :: kout   !< output channel for object data structure

  !! local
  character(*), parameter :: s_name='spl3d_writeg' !< subroutine name
  integer(ki4), parameter :: ipcpt=3 !< component of field
  integer(ki4), parameter :: ip1=80 !< 1-coordinate (R) index
  integer(ki4), parameter :: ip2=6 !< 2-coordinate (Z) index
  integer(ki4), parameter :: ip3=16 !< 3-coordinate (\f$ \zeta \f$) index
  integer(ki4), parameter :: ipk=1 !< \f$ k \f$-coordinate (\f$ \zeta \f$ wavenumber) index
  integer(ki4), parameter :: ip0=2 !< \f$ a \f$ or \f$ b \f$ part (\f$ \zeta \f$ wavenumber) index
  real(kr8), parameter :: cutoff=1.e-4 !< cutoff for relative mode amp
  real(kr8) :: za   !< a part of complex mode amplitude
  real(kr8) :: zb   !< b part of complex mode amplitude
  real(kr8) :: zamp   !< amplitude of complex mode
  real(kr8) :: zampt  !< sum of all amplitudes of complex mode
  integer(ki4) :: in1   !< first dimension
  integer(ki4) :: in2   !< second dimension
  integer(ki4) :: j1   !< loop counter for first dimension
  integer(ki4) :: j2   !< loop counter for second dimension
  integer(ki4) :: ipm   !< end of parity loop
  type(spl2d_t) :: spl2d   !< object data structure

  plot_type: select case (kchar)
  case('sampl')
     ! save 2-D dependence for preselected
     j=ipcpt
     !DALL     do j=1,ipcpt !DALL
     write(kout,'(1x,a,3(1x,i3))',iostat=status) &
     '# cpt,ip0,ipk=',j,ip0,ipk
     call spl2d_writeg(self%farm(j,ip0,ipk),'sampl',kout)
     !DALL     end do !DALL

  case('strip')
     ! save 1-D dependence for preselected
     write(kout,'(1x,a,'//cfmt2v,iostat=status) &
     '# R,Z=',(ip1-1)*self%farm(1,1,0)%h1+self%farm(1,1,0)%org1,  &
     (ip2-1)*self%farm(1,1,0)%h2+self%farm(1,1,0)%org2
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Error writing strip')
     end if
     do j=self%ncpt,1,-1
        do l=0,self%nk
           za=self%farm(j,1,l)%sampl(ip1,ip2)
           zb=self%farm(j,2,l)%sampl(ip1,ip2)
           zamp=sqrt(za**2+zb**2)
           write(kout,'(1x,i4,'//cfmt2v,iostat=status) l,zamp,za,zb
           if(status/=0) then
              call log_error(m_name,s_name,2,error_fatal,'Error writing strip')
           end if
        end do
        write(kout,*,iostat=status)
        call log_write_check(m_name,s_name,3,status)
     end do

  case('modecount')
     !! allocate 2D storage
     in1=self%farm(1,1,0)%n1p
     in2=self%farm(1,1,0)%n2p
     if (in1*in2>0) then
        allocate(work2(in1,in2), stat=status)
        call log_alloc_check(m_name,s_name,12,status)
     else
        call log_error(m_name,s_name,13,error_fatal,'no data')
     end if

     ! analyse 1-D k dependence for all R,Z
     do j=1,self%ncpt

        write(kout,'(a,i1)') 'Analysis for component = ',j

        do j2=1,in2
           do j1=1,in1

              ! first total amplitude
              zampt=0
              do l=0,self%nk
                 za=self%farm(j,1,l)%sampl(j1,j2)
                 zb=self%farm(j,2,l)%sampl(j1,j2)
                 zamp=sqrt(za**2+zb**2)
                 zampt=zampt+zamp
              end do
              ! now find largest k exceeding fraction given by cut-off
              work2(j1,j2)=0
              if (zampt>=cutoff) then
                 ! ignore very small total amps
                 do l=self%nk,0,-1
                    za=self%farm(j,1,l)%sampl(j1,j2)
                    zb=self%farm(j,2,l)%sampl(j1,j2)
                    zamp=sqrt(za**2+zb**2)
                    if (zamp>cutoff*zampt) then
                       work2(j1,j2)=l
                       exit
                    end if
                 end do
              end if

           end do
        end do

        ! copy to spline data structure
        call spl2d_init(spl2d,work2,&
 &      in1-1,in2-1,&
 &      self%farm(1,1,0)%pos1(1),self%farm(1,1,0)%pos2(1),&
 &      self%farm(1,1,0)%h1,self%farm(1,1,0)%h2,self%farm(1,1,0)%nord)
        !     call spl2d_initpart(spl2d)
        ! output
        call spl2d_writeg(spl2d,'sampl',kout)

        call spl2d_delete(spl2d)

     end do
     deallocate(work2)

  case('countout')

     !! allocate 2D storage
     in1=self%farm(1,1,0)%n1p
     in2=self%farm(1,1,0)%n2p
     if (in1*in2>0) then
        allocate(work2(in1,in2), stat=status)
        call log_alloc_check(m_name,s_name,12,status)
     else
        call log_error(m_name,s_name,13,error_fatal,'no data')
     end if

     ! output mask k limits for all R,Z
     do j=1,self%ncpt

        write(kout,'(a,i1)') 'Analysis for component = ',j

        k=min(2,self%parity(j))
        ipm=max(1,self%parity(j)-1)
        looparity: do jp=1,ipm

           do j2=1,in2
              do j1=1,in1
                 work2(j1,j2)=self%mask(j,k,j1,j2)
              end do
           end do

           ! copy to spline data structure
           call spl2d_init(spl2d,work2,&
 &         in1-1,in2-1,&
 &         self%farm(1,1,0)%pos1(1),self%farm(1,1,0)%pos2(1),&
 &         self%farm(1,1,0)%h1,self%farm(1,1,0)%h2,self%farm(1,1,0)%nord)
           ! output
           call spl2d_writeg(spl2d,'sampl',kout)

           call spl2d_delete(spl2d)

           k=k-1
           if (jp<ipm) write(kout,'(a)') ' '
        end do looparity

     end do

     deallocate(work2)

  case default

  end select plot_type

end subroutine spl3d_writeg
!---------------------------------------------------------------------
!> write in data format
subroutine spl3d_write(self,kout)

  !! arguments
  type(spl3d_t), intent(in) :: self   !< object data structure
  integer, intent(in) :: kout   !< output channel for object data structure

  !! local
  character(*), parameter :: s_name='spl3d_write' !< subroutine name

  write(kout,*,iostat=status) 'n0'
  call log_write_check(m_name,s_name,10,status)
  write(kout,*,iostat=status) self%n0
  call log_write_check(m_name,s_name,11,status)
  write(kout,*,iostat=status) 'ncpt'
  call log_write_check(m_name,s_name,12,status)
  write(kout,*,iostat=status) self%ncpt
  call log_write_check(m_name,s_name,13,status)
  write(kout,*,iostat=status) 'nk'
  call log_write_check(m_name,s_name,14,status)
  write(kout,*,iostat=status) self%nk
  call log_write_check(m_name,s_name,15,status)
  write(kout,*,iostat=status) 'nstate'
  call log_write_check(m_name,s_name,16,status)
  write(kout,*,iostat=status) self%nstate
  call log_write_check(m_name,s_name,17,status)
  write(kout,*,iostat=status) 'n3'
  call log_write_check(m_name,s_name,18,status)
  write(kout,*,iostat=status) self%n3
  call log_write_check(m_name,s_name,19,status)
  write(kout,*,iostat=status) 'h3'
  call log_write_check(m_name,s_name,20,status)
  write(kout,*,iostat=status) self%h3
  call log_write_check(m_name,s_name,21,status)
  write(kout,*,iostat=status) 'rh3'
  call log_write_check(m_name,s_name,22,status)
  write(kout,*,iostat=status) self%rh3
  call log_write_check(m_name,s_name,23,status)
  write(kout,*,iostat=status) 'kmin'
  write(kout,*,iostat=status) self%kmin
  call log_write_check(m_name,s_name,24,status)
  write(kout,*,iostat=status) 'kmax'
  write(kout,*,iostat=status) self%kmax
  call log_write_check(m_name,s_name,25,status)
  write(kout,*,iostat=status) 'parity'
  write(kout,*,iostat=status) self%parity
  call log_write_check(m_name,s_name,26,status)

  !-----------------------------------------------------------------------
  !              Write 1D storage
  ! position 3 array
  write(kout,*,iostat=status) 'pos3'
  write(kout,*,iostat=status) self%pos3
  call log_write_check(m_name,s_name,51,status)
  ! position k array
  write(kout,*,iostat=status) 'posk'
  write(kout,*,iostat=status) self%posk
  call log_write_check(m_name,s_name,53,status)
  ! k offset array
  write(kout,*,iostat=status) 'poskof'
  write(kout,*,iostat=status) self%poskof
  call log_write_check(m_name,s_name,54,status)

  !-----------------------------------------------------------------------
  !              Write 3D storage

  !! data array
  ij=0
  do l=0,self%nk
     do k=1,self%n0
        do j=1,self%ncpt
           call spl2d_write( self%farm(j,k,l),kout )
           if (lplog2d) print '("2D spline written index = ",3(i9,1x))',j,k,l
           ij=ij+1
           if (lplog2d) call log_value("number of 2D splines written ",ij)
        end do
     end do
  end do

  !-----------------------------------------------------------------------
  !              Write 4D storage

  write(kout,*,iostat=status) 'mask'
  write(kout,*,iostat=status) self%mask
  call log_write_check(m_name,s_name,60,status)

  print '("finished writing 3D spline")'
  call log_error(m_name,s_name,99,log_info,'finished writing 3D spline')

end subroutine spl3d_write
!---------------------------------------------------------------------
!> create dummy part of spl3d data structure
subroutine spl3d_initpart(self)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='spl3d_initpart' !< subroutine name

  !-----------------------------------------------------------------------

  ij=0
  do l=0,self%nk
     do k=1,self%n0
        do j=1,self%ncpt
           call  spl2d_initpart(self%farm(j,k,l))
           if (lplog2d) print '("created dummy 2D spline index = ",3(i9,1x))',j,k,l
           ij=ij+1
           if (lplog2d) call log_value("created dummy 2D spline number ",ij)
        end do
     end do
  end do

end subroutine spl3d_initpart
!---------------------------------------------------------------------
!> get \f$ \psi,\theta \f$ limits of interpolation
subroutine spl3d_ptlimits(self,p1min,p1max,p2min,p2max)

  !! arguments
  type(spl3d_t), intent(in) :: self   !< object data structure
  real(kr8), intent(out) :: p1min !< minimum value of position-1
  real(kr8), intent(out) :: p1max !< maximum value of position-1
  real(kr8), intent(out) :: p2min !< minimum value of position-2
  real(kr8), intent(out) :: p2max !< maximum value of position-2

  !! local
  character(*), parameter :: s_name='spl3d_ptlimits' !< subroutine name

  ij=0
  do l=0,0 !self%nk
     do k=1,1 !self%n0
        do j=1,1 !self%ncpt
           call  spl2d_ptlimits(self%farm(j,k,l),p1min,p1max,p2min,p2max)
           if (lplog2d) print '("got limits 2D spline index = ",3(i9,1x))',j,k,l
           ij=ij+1
           if (lplog2d) call log_value("got limits 2D spline number ",ij)
        end do
     end do
  end do

end subroutine spl3d_ptlimits
!---------------------------------------------------------------------
!> scale points at which spline defined
subroutine spl3d_ptscale(self,qtfmdata)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure
  type(quantfm_t), intent(inout) :: qtfmdata   !< data defining transform

  !! local
  character(*), parameter :: s_name='spl3d_ptscale' !< subroutine name
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zposq !< local variable
  integer(ki4) :: inqtfm   !< save transform type

  ij=0
  do l=0,self%nk
     do k=1,self%n0
        do j=1,self%ncpt
           call  spl2d_ptscale(self%farm(j,k,l),qtfmdata)
           if (lplog2d) print '("scaled point 2D spline index = ",3(i9,1x))',j,k,l
           ij=ij+1
           if (lplog2d) call log_value("scaled point 2D spline number ",ij)
        end do
     end do
  end do
  ! spacings
  inqtfm=qtfmdata%nqtfm
  ! case 1 is pure scaling
  qtfmdata%nqtfm=1
  ! don't bother to scale h3 or pos3, not used
  !     zpos%posvec=(/0._kr8,0._kr8,self%h3/)
  !     zposq=position_qtfm(zpos,qtfmdata)
  !     self%h3=zposq%posvec(3)
  !     self%rh3=1/zposq%posvec(3)
  ! k values are set rather than scaled
  zpos%posvec=1
  zposq=position_invqtfm(zpos,qtfmdata)
  do l=0,self%nk
     self%posk(l)=l*zposq%posvec(3)
  end do

  qtfmdata%nqtfm=2
  ! k offset values
  zpos%posvec=0
  zposq=position_invqtfm(zpos,qtfmdata)
  do l=0,self%nk
     self%poskof(l)=l*zposq%posvec(3)
  end do
  ! restore
  qtfmdata%nqtfm=inqtfm

end subroutine spl3d_ptscale
!---------------------------------------------------------------------
!> scale spline values
subroutine spl3d_scale(self,pfac)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure
  real(kr8), intent(in) :: pfac   !< scale factor

  !! local
  character(*), parameter :: s_name='spl3d_scale' !< subroutine name

  ij=0
  do l=0,self%nk
     do k=1,self%n0
        do j=1,self%ncpt
           call  spl2d_scale(self%farm(j,k,l),pfac)
           if (lplog2d) print '("scaled value 2D spline index = ",3(i9,1x))',j,k,l
           ij=ij+1
           if (lplog2d) call log_value("scaled value 2D spline number ",ij)
        end do
     end do
  end do
  print '("3D spline scaled by factor = ",1pg15.8)',pfac
  call log_value("3D spline scaled by factor ",pfac)

end subroutine spl3d_scale
!---------------------------------------------------------------------
!> calculate coefficient array of splines
subroutine spl3d_coeff(self)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='spl3d_coeff' !< subroutine name

  !-----------------------------------------------------------------------
  !              spline  coefficients via direct product
  ij=0
  do l=0,self%nk
     do k=1,self%n0
        do j=1,self%ncpt
           call  spl2d_coeff(self%farm(j,k,l))
           if (lplog2d) print '("created coeff 2D spline index = ",3(i9,1x))',j,k,l
           ij=ij+1
           if (lplog2d) call log_value("created coeff 2D spline number ",ij)
        end do
     end do
  end do

end subroutine spl3d_coeff
!---------------------------------------------------------------------
!> test zeta transform
subroutine spl3d_testtfm(self,field,pfac,phchange)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure
  real(kr8), dimension(:,:,:,:), intent(in) :: field   !< original field
  real(kr8), intent(in) :: pfac   !< scale factor
  integer(ki4), intent(in) :: phchange !< correct phase of \f$ \zeta \f$

  !! local
  character(*), parameter :: s_name='spl3d_testtfm' !< subroutine name
  integer(ki4) :: j1   !< loop counter for first dimension \f$ R \f$
  integer(ki4) :: j2   !< loop counter for second dimension \f$ Z \f$
  integer(ki4) :: jj   !< loop counter for component
  integer(ki4) :: jk   !< loop counter for third dimension
  integer(ki4) :: jl   !< loop counter for third dimension
  integer(ki4) :: in1   !< first nominal dimension
  integer(ki4) :: in2   !< second nominal dimension
  integer(ki4) :: il   !< third nominal dimension
  integer(ki4) :: ilw   !< third nominal dimension wrapped
  real(kr8), dimension(3) :: zpt !< vector point
  real(kr8), dimension(3) :: zinit !< original vector value
  real(kr8), dimension(3) :: zeval !< calculated vector value
  real(kr8), dimension(3) :: zdiff !< calculated vector difference
  real(kr8) :: zeta !< angle
  real(kr8) :: zmdiff !< calculated maximum difference
  real(kr8) :: zdiffsum !< calculated summed maximum difference
  real(kr8) :: zrfac  !< local variable
  real(kr8), dimension(3) :: zswap !< swap
  integer(ki4),dimension(4) :: impt   !< point of maximum difference
  real(kr8), dimension(:,:,:,:), allocatable :: work4 !< 4D work array
  real(kr8), dimension(:), allocatable :: pos3o !< offset pos3

  ! allocate 4D storage
  ! could refactor to avoid this, but would be tedious
  in1=self%farm(1,1,0)%n1p
  in2=self%farm(1,1,0)%n2p
  allocate(work4(self%n3,self%ncpt,in1,in2), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  allocate(pos3o(self%n3), stat=status)
  call log_alloc_check(m_name,s_name,2,status)

  ! set up pos3 corresponding to field arrangement
  pos3o=self%pos3-self%pos3(1)

  ! analyse 1-D wavenumber dependence for all R,Z,zeta
  zmdiff=0
  zdiffsum=0
  do j2=1,in2
     zpt(2)=self%farm(1,1,0)%pos2(j2)
     do j1=1,in1
        zpt(1)=self%farm(1,1,0)%pos1(j1)
        zrfac=1/(zpt(1)**spl3d_rpower)

        il=0
        do jl=0,self%nk
           do jk=1,self%n0
              il=il+1
              !f! added debugging clause
              !f      if (j1==ip1diag.AND.j2==ip2diag.AND.il==ip3diag) then
              il=min(il,self%n3)
              do jj=1,self%ncpt
                 zinit(jj)=pfac*field(il,jj,j1,j2)
              end do

              if (phchange==1) then
                 ! disentangle zeta for pi phase change
                 if (il<=self%n3/2) then
                    ilw=il+self%n3/2
                 else
                    ilw=il-self%n3/2
                 end if
                 ! pos3o got mapped
                 zeta=pos3o(il)
              else if (phchange==2.OR.phchange==4) then
                 ! disentangle zeta for phi -> zeta map
                 ilw=self%n3+1-il
                 zeta=pos3o(ilw)
              else
                 ilw=il
                 zeta=pos3o(il)
              end if

              zpt(3)=pos3o(ilw)
              if (lpmask) then
                 call spl3d_evalm(self,zpt,zeval)
              else
                 call spl3d_eval(self,zpt,zeval)
              end if
              ! convert back to Cartesians
              zswap(3)=zeval(2)
              if (phchange==1) then
                 zswap(1)=zeval(1)*cos(zeta)+zpt(1)*zeval(3)*sin(zeta)
                 zswap(2)=-zeval(1)*sin(zeta)+zpt(1)*zeval(3)*cos(zeta)
              else if (phchange==2) then
                 zswap(1)=zeval(1)*cos(zeta)-zpt(1)*zeval(3)*sin(zeta)
                 zswap(2)=zeval(1)*sin(zeta)+zpt(1)*zeval(3)*cos(zeta)
              else
                 zswap(1)=zeval(1)*cos(zeta)-zpt(1)*zeval(3)*sin(zeta)
                 zswap(2)=-zeval(1)*sin(zeta)-zpt(1)*zeval(3)*cos(zeta)
                 if (phchange==3) zswap=-zswap ! only works if no axisymm BR,BZ
              end if
              do jj=1,self%ncpt
                 work4(il,jj,j1,j2)=zswap(jj)*zrfac
              end do
              zdiff=abs(zswap*zrfac-zinit)
              zdiffsum=zdiffsum+sum(zdiff)
              !!    write(*,*) 'j1,j2,il-init',j1,j2,il,zinit
              !!    write(*,*) 'j1,j2,il-swap',j1,j2,il,zswap
              if (maxval(zdiff)>zmdiff) then
                 zmdiff=maxval(zdiff)
                 impt=(/0,j1,j2,il/)
              end if
              !f      end if
              !f! end debugging clause

           end do
        end do

     end do
  end do

  call log_value("impt(1) ",impt(1))
  call log_value("impt(2) ",impt(2))
  call log_value("impt(3) ",impt(3))
  call log_value("impt(4) ",impt(4))
  print '("maximum difference after transform = ",1pg15.8)',zmdiff
  call log_value("maximum difference after transform",zmdiff)
  zdiffsum=zdiffsum/(self%ncpt*in1*in2*self%n3)
  print '("mean difference after transform = ",1pg15.8)',zdiffsum
  call log_value("mean difference after transform",zdiffsum)

  !      call spl3d_init(selfnew,work4,self%ncpt,self%n0,&
  !     &in1,in2,self%n3,&
  !     &self%farm(1,1,0)%pos1,self%farm(1,1,0)%pos2,self%pos3,self%farm(1,1,0)%nord,&
  !     &self%kmin,self%kmax)

  deallocate(work4)

end subroutine spl3d_testtfm
!---------------------------------------------------------------------
!> evaluate spline
subroutine spl3d_eval(self,pt,pe)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure
  real(kr8), dimension(3), intent(in) :: pt   !< coordinates of point i(third coordinate of point on \f$ (0,2\pi) \f$)
  real(kr8), dimension(3), intent(out) :: pe   !< value of spline vector at pt

  !! local
  character(*), parameter :: s_name='spl3d_eval' !< subroutine name
  integer(ki4) :: jk   !<  variable for loop over wavenumber
  real(kr8),dimension(3) :: ze   !< point value
  real(kr8),dimension(3) :: zz   !< cumulated point value
  real(kr8) :: zcosk  !< intermediate variable
  real(kr8) :: zsink  !< intermediate variable
  real(kr8) :: zxi  !< \f$ \xi \f$
  integer(ki4) :: ikmin   !< min wavenumber to be used in evaluation
  integer(ki4) :: ikmax   !< max wavenumber to be used in evaluation
  integer(ki4) :: icall   !< if unity, evaluate coeffs, else reuse previous

  allocate(work3(self%ncpt,self%n0,0:self%nk), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  ikmin=minval(self%kmin)
  ikmax=maxval(self%kmax)
  if (ikmax<0) then
     !! fix up for negative kmax
     self%kmax=self%nk
     ikmax=maxval(self%kmax)
  end if
  !f     write(*,*) 'posk', self%posk
  !f      write(*,*) 'pt', pt
  !  get Fourier coeffts
  icall=0
  loopget: do l=ikmin,ikmax
     do k=1,self%n0
        do j=1,self%ncpt
           if (l>=self%kmin(j).AND.l<=self%kmax(j)) then
              icall=icall+1
              call  spl2d_evaln(self%farm(j,k,l),pt(1),pt(2),icall,ze(j))
           else
              ze(j)=0
           end if
        end do
        !f     write(*,*) 'j,k,l', j,k,l,ze
        work3(:,k,l)=ze
     end do
  end do loopget
  !  evaluate sum
  zz=work3(:,1,0)
  do jk=max(ikmin,1),min(ikmax,self%nk-1)
     zxi=self%posk(jk)*pt(3)+self%poskof(jk)
     zcosk=cos(zxi)
     zsink=sin(zxi)
     zz=zz+work3(:,1,jk)*zcosk+work3(:,2,jk)*zsink
  end do
  if (ikmax==self%nk) zz=zz+work3(:,1,ikmax)*&
 &cos(self%posk(ikmax)*pt(3)+self%poskof(ikmax))
  ! return
  pe=zz
  ! tidy up
  deallocate(work3)

end subroutine spl3d_eval
!---------------------------------------------------------------------
!> evaluate spline using mask
subroutine spl3d_evalm(self,pt,pe)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure
  real(kr8), dimension(3), intent(in) :: pt   !< coordinates of point (third coordinate of point on \f$ (0,2\pi) \f$)
  real(kr8), dimension(3), intent(out) :: pe   !< value of spline vector at pt

  !! local
  character(*), parameter :: s_name='spl3d_evalm' !< subroutine name
  integer(ki4) :: jk   !<  variable for loop over wavenumber
  real(kr8) :: ze   !< point value
  real(kr8) :: zz   !< cumulated point value
  real(kr8) :: zxi  !< \f$ \xi \f$
  integer(ki4) :: ikmin   !< min wavenumber to be used in evaluation
  integer(ki4) :: ikmax   !< max wavenumber to be used in evaluation
  integer(ki4) :: i1   !< largest integer index positioned less then point 1 coordinate
  integer(ki4) :: i2   !< largest integer index positioned less then point 2 coordinate
  integer(ki4) :: icall   !< if unity, evaluate coeffs, else reuse previous
  real(kr8), dimension(:,:), allocatable :: zwork2 !< 2D work array
  real(kr8), dimension(:,:), allocatable :: workc !< work array for coeffts

  allocate(zwork2(self%n0,0:self%nk), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  allocate(workc(self%n0,0:self%nk), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  !PE      write(*,*) 'R,Z,xi=', pt(1),pt(2),self%posk(1)*pt(3)+self%poskof(1) !PE
  !f     write(*,*) 'posk', self%posk
  !f      write(*,*) 'pt', pt
  !  get Fourier coeffts
  call  spl2d_locpos(self%farm(1,1,0),pt(1),pt(2),i1,i2)
  ! set up coefficient array
  !     ikmin=minval(self%parity)-1
  ikmin=0
  if (lpmask3) then
     ikmax=self%mask(3,1,i1,i2)
  else
     ikmax=maxval(self%mask(:,:,i1,i2))
  end if
  do jk=ikmin,ikmax
     zxi=self%posk(jk)*pt(3)+self%poskof(jk)
     workc(1,jk)=cos(zxi)
     workc(2,jk)=sin(zxi)
  end do

  !WAtest. test data has spurious even-parity content, hence lpddtest
  ! outer loop over field components
  icall=0
  do j=1,self%ncpt
     zwork2(1:2,0:ikmax)=0

     ! first evaluate mode coefficients
     k=min(2,self%parity(j))
     looparity: do jp=1,max(1,self%parity(j)-1)
        ikmin=k-1
        if (lpmask3) then
           ikmax=self%mask(3,1,i1,i2)
        else
           ikmax=self%mask(j,k,i1,i2)
        end if
        !KMAX      write(*,*) 'mask, j,k,i1,i2,ikmax=',j,k,i1,i2,ikmax !KMAX
        do l=ikmin,ikmax
           icall=icall+1
           call  spl2d_evaln(self%farm(j,k,l),pt(1),pt(2),icall,ze)
           !f     write(*,*) 'j,k,l', j,k,l,icall,ze
           zwork2(k,l)=ze
        end do
        k=k-1
     end do looparity

     if (lpddtest) then
        if(self%parity(j)==2)then
           icall=icall+1
           !WAtest include zero mode amplitude even if only odd parity set
           call  spl2d_evaln(self%farm(j,1,0),pt(1),pt(2),icall,ze)
           zwork2(1,0)=ze
        end if
     end if

     ! then evaluate field component using mode coefficients
     zz=zwork2(1,0)
     k=min(2,self%parity(j))
     looparity2: do jp=1,max(1,self%parity(j)-1)
        ikmin=k-1
        if (lpmask3) then
           ikmax=self%mask(3,1,i1,i2)
        else
           ikmax=self%mask(j,k,i1,i2)
        end if
        !  evaluate sum
        do jk=max(ikmin,1),min(ikmax,self%nk-1)
           zz=zz+zwork2(k,jk)*workc(k,jk)
        end do
        ! special for even parity
        if (k==1.AND.ikmax==self%nk) zz=zz+zwork2(1,ikmax)*workc(1,ikmax)
        k=k-1
     end do looparity2

     ! return
     pe(j)=zz
  end do

  !PE      write(*,*) 'pe=',pe  !PE
  ! tidy up
  deallocate(zwork2)
  deallocate(workc)

end subroutine spl3d_evalm
!---------------------------------------------------------------------
!> differentiate with respect to coordinate
subroutine spl3d_deriv(self,dself,kcoo)

  !! arguments
  type(spl3d_t), intent(inout) :: self   !< object data structure
  type(spl3d_t), intent(out) :: dself   !< derivative object data structure
  integer(ki4), intent(in) :: kcoo !< coordinate direction

  !! local
  character(*), parameter :: s_name='spl3d_deriv' !< subroutine name
  integer(ki4) :: il   !<  integer work
  integer(ki4) :: in1   !< first dimension
  integer(ki4) :: in2   !< second dimension
  type(spl2d_t) :: spl2d   !< derivative spline object

  ! spline inner 2D dimensions
  in1=self%farm(1,1,0)%n1p
  in2=self%farm(1,1,0)%n2p
  ! could refactor to avoid this, but would be tedious
  allocate(work4(self%n3,self%ncpt,in1,in2), stat=status)
  call log_alloc_check(m_name,s_name,2,status)

  ij=0
  il=0
  do l=0,self%nk
     do k=1,self%n0
        il=il+1

        if (il<=self%n3) then
           do j=1,self%ncpt
              call  spl2d_deriv(self%farm(j,k,l),spl2d,kcoo)
              work4(il,j,1:in1,1:in2)=spl2d%sampl
              if (lplog2d) print '("derivative 2D spline index = ",3(i9,1x))',j,k,l
              ij=ij+1
              if (lplog2d) call log_value("derivative 2D spline number ",ij)
           end do
        end if

     end do
  end do

  call spl3d_init(dself,work4,self%ncpt,self%n0,&
 &in1,in2,self%n3,&
 &self%farm(1,1,0)%pos1,self%farm(1,1,0)%pos2,self%pos3,self%farm(1,1,0)%nord,&
 &self%kmin,self%kmax)

  deallocate(work4)

end subroutine spl3d_deriv

end module spl3d_m
