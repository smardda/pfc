module position_m

  use const_kind_m
  use log_m
  use misc_m
  use position_h

  implicit none
  private

! public subroutines
  public :: &
 &position_inbox, & !< test whether within real box
 &position_innbox, & !< test whether within normalised box
 &position_asvecb, & !< encode as binary vector
 &position_tfm, & !< transform position vector
 &position_invtfm, & !< inverse transform position vector
 &position_qtfm, & !< quantise position vector
 &position_invqtfm, & !< inverse quantise position vector
 &position_scaleunits, & !< change units of transform in position space
 &position_readv, & !< read in visualisation format
 &position_writev, & !< write in visualisation format
 &position_readcon, & !< read position control data
 &position_unitfm, & !< set identity position transform
 &position_unitfmq, & !< set identity quantising position transform
 &position_readtfm, & !< read position control data
 &position_tfmquery, & !< interrogate position control data
 &position_writetfm, & !< write position control data
 &position_readlis, & !< read position list data
 &position_copylis, & !< copy position list data
 &position_readonlylis, & !< read (general vtk format) only list of position coordinates
 &position_readveclis, & !< read position list data as vector
 &position_deleteveclis, & !< delete list of vectors
 &position_writelis, & !< write position list data
 &position_lenlis, & !< physical length of position list data
 &position_tfmlis, & !< transform position list
 &position_invtfmlis, & !< inverse transform position list
 &position_qtfmlis, & !< quantise position list
 &position_invqtfmlis ! inverse quantise position list


!public variables

! private types

! private variables
  character(*), parameter :: m_name='position_m' !< module name
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  integer   :: status   !< error status
  integer :: nin   !< input channel for position list data
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  logical :: iltest !< logical flag


  contains
!---------------------------------------------------------------------
!> is position in real box
  pure logical function position_inbox(self,box)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  real(kr4), dimension(3,2), intent(in) :: box !< bounding box corners


!! local
  character(*), parameter :: s_name='position_inbox' !< subroutine name
  logical :: linint   !< position in interval
  logical :: linbox   !< position in interval
  integer(ki4) :: jj   !< loop counter

  linbox=.true.
  do jj=1,3
     linint=( box(jj,1)<=self%posvec(jj) )
     linint=linint.and.( self%posvec(jj)<box(jj,2) )
     linbox=linbox.and.linint
  end do

  position_inbox=linbox

end function position_inbox
!---------------------------------------------------------------------
!> is position in logical box
  pure logical function position_innbox(self,kcorn,kext)

!! arguments
  type(posvecl_t), intent(in) :: self   !< quantised position data
  integer(ki2), dimension(3), intent(in) :: kcorn !< bounding box corner
  integer(ki2), dimension(3), intent(in) :: kext !< bounding box etents


!! local
  character(*), parameter :: s_name='position_innbox' !< subroutine name
  logical :: linbox   !< position in interval
  integer(ki4) :: jj   !< loop counter
  integer(ki2), dimension(3) :: ivec !< integer position data
  integer(ki2) :: ie !< integer element
  integer(ki2) :: iesum !< integer sum

  if (any(ivec<0)) then
     linbox=.false.
  else
     ivec=int(self%posvec)
     iesum=0
     do jj=1,3
        ie=ishft( ieor(kcorn(jj),ivec(jj)), -kext(jj) )
        iesum=iesum+ie
     end do

     linbox=(iesum==0)
  end if
  position_innbox=linbox

end function position_innbox
!---------------------------------------------------------------------
!> convert transformed position to binary representation (3 X integer*2)
  pure type(posvecb_t) function position_asvecb(self)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data


!! local
  character(*), parameter :: s_name='position_asvecb' !< subroutine name
  integer(ki2), dimension(3) :: ivecb   !< binary rep of position data
  ivecb=int(self%posvec)

  position_asvecb%posb=ivecb

end function position_asvecb
!---------------------------------------------------------------------
!> transform in position space
  pure type(posvecl_t) function position_tfm(self,tfmdata)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  type(tfmdata_t), intent(in) :: tfmdata   !< data defining transform


!! local
  character(*), parameter :: s_name='position_tfm' !< subroutine name
  real(kr4), dimension(3) :: zvec   !< transformed position data
  integer(ki4) :: intfm   !< number of transform
  integer(ki4) :: jj   !< loop

  intfm=tfmdata%ntfm
  transform_number :  select case (intfm)
  case(1)
     zvec=self%posvec*tfmdata%scale
  case(2)
     zvec=self%posvec*tfmdata%scale+tfmdata%offset
  case(3)
     do jj=1,3
        zvec(jj)=dot_product(self%posvec,tfmdata%matrix(jj,:))+tfmdata%offset(jj)
     end do
  case(4)
     do jj=1,3
        zvec(jj)=tfmdata%offset(jj)+dot_product(self%posvec-tfmdata%offset,tfmdata%matrix(jj,:))
     end do
  case(5)
     zvec=tfmdata%scale*(self%posvec+tfmdata%offset)
  case(6)
     zvec=self%posvec+tfmdata%scale*tfmdata%offset
  end select transform_number

  position_tfm%posvec=zvec

end function position_tfm
!---------------------------------------------------------------------
!> inverse transform in position space
  pure type(posvecl_t) function position_invtfm(self,tfmdata)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  type(tfmdata_t), intent(in) :: tfmdata   !< data defining transform


!! local
  character(*), parameter :: s_name='position_invtfm' !< subroutine name
  real(kr4), dimension(3) :: zvec   !< transformed position data
  integer(ki4) :: intfm   !< number of transform
  integer(ki4) :: jj   !< loop

  intfm=tfmdata%ntfm
  transform_number :  select case (intfm)
  case(1)
     zvec=self%posvec/tfmdata%scale
  case(2)
     zvec=(self%posvec-tfmdata%offset)/tfmdata%scale
  case(3)
! assume transform is orthogonal matrix, i.e. inverse=transpose
     do jj=1,3
        zvec(jj)=dot_product(self%posvec-tfmdata%offset,tfmdata%matrix(:,jj))
     end do
  case(4)
! assume transform is orthogonal matrix, i.e. inverse=transpose
     do jj=1,3
        zvec(jj)=tfmdata%offset(jj)+dot_product(self%posvec-tfmdata%offset,tfmdata%matrix(:,jj))
     end do
  case(5)
     zvec=self%posvec/tfmdata%scale-tfmdata%offset
  case(6)
     zvec=self%posvec-tfmdata%scale*tfmdata%offset
  end select transform_number

  position_invtfm%posvec=zvec

end function position_invtfm
!---------------------------------------------------------------------
!> quantise in position space
  pure type(posvecl_t) function position_qtfm(self,qtfmdata)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  type(quantfm_t), intent(in) :: qtfmdata   !< data defining transform


!! local
  character(*), parameter :: s_name='position_qtfm' !< subroutine name
  real(kr4), dimension(3) :: zvec   !< quantised position data
  integer(ki4) :: inqtfm   !< number of transform

  inqtfm=qtfmdata%nqtfm
  transform_number :  select case (inqtfm)
  case(1)
     zvec=self%posvec*qtfmdata%rhmin
  case(2)
     zvec=self%posvec*qtfmdata%rhmin+qtfmdata%offvec
  end select transform_number

  position_qtfm%posvec=zvec

end function position_qtfm
!---------------------------------------------------------------------
!> inverse quantise in position space
  pure type(posvecl_t) function position_invqtfm(self,qtfmdata)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  type(quantfm_t), intent(in) :: qtfmdata   !< data defining transform


!! local
  character(*), parameter :: s_name='position_invqtfm' !< subroutine name
  real(kr4), dimension(3) :: zvec   !< quantised position data
  integer(ki4) :: inqtfm   !< number of transform

  inqtfm=qtfmdata%nqtfm
  transform_number :  select case (inqtfm)
  case(1)
     zvec=self%posvec*qtfmdata%hmin
  case(2)
     zvec=(self%posvec-qtfmdata%offvec)*qtfmdata%hmin
  end select transform_number

  position_invqtfm%posvec=zvec

end function position_invqtfm
!---------------------------------------------------------------------
!> change units of transform in position space
subroutine position_scaleunits(tfmdata,pfac)

  !! arguments
  type(tfmdata_t), intent(inout) :: tfmdata   !< data defining transform
  real(kr4), intent(in) :: pfac  !< scales only positional data

  !! local
  character(*), parameter :: s_name='position_scaleunits' !< subroutine name
  real(kr4), dimension(3) :: zoffset   !< change offset position data

  zoffset=tfmdata%offset*pfac

  tfmdata%offset=zoffset

end subroutine position_scaleunits
!---------------------------------------------------------------------
!> read in position data
subroutine position_readv(self,kin)

  !! arguments
  type(posvecl_t), intent(out) :: self   !< position data
  integer, intent(in) :: kin   !< input channel for position data


  !! local
  character(*), parameter :: s_name='position_readv' !< subroutine name

  read(kin,*,iostat=status) self%posvec
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading Particle position')
  end if

end subroutine position_readv
!---------------------------------------------------------------------
!> output position vector
subroutine position_writev(self,kplot)

  !! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  integer, intent(in) :: kplot   !< output channel for position data


  !! local
  character(*), parameter :: s_name='position_writev' !< subroutine name

  write(kplot,cfmtbv1,iostat=status) self%posvec
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error writing Particle position')
  end if

end subroutine position_writev
!---------------------------------------------------------------------
!> read position control data
subroutine position_readcon(ztfmdata,kin,flag)

  !! arguments
  type(tfmdata_t), intent(out) :: ztfmdata   !< position transform numeric controls
  integer :: kin  !< unit for input
  integer(ki4), optional :: flag  !< local variable

  !! local
  character(*), parameter :: s_name='position_readcon' !< subroutine name
  real(kr4), dimension(3,3):: position_matrix  !< transformation matrix 3x3
  real(kr4), dimension(9):: position_vector0  !< work transformation matrix 3x3 as vector
  real(kr4), dimension(9):: position_vector  !< transformation matrix 3x3 as vector
  real(kr4), dimension(3):: position_scale  !< transformation scale 3-vector
  real(kr4), dimension(3):: position_offset  !< transformation offset 3-vector
  integer(ki4):: position_transform  !< transformation type
  character(len=20) :: transform_desc !< describes transformation type
  character(len=20) :: transform_id !< identifies transformation
  integer(ki4), parameter :: npdict = 10 !< size of dictionary
  !> dictionary of long name transformations
  character(len=24), dimension(npdict), parameter :: dictlong = &  !< .
 &(/'cartesian_scale         ','scale+offset            ',&
 &'full_matrix_mult+offset ','offset+full_matrix_mult ',&
 &'five                    ','poloidal_tilt           ',&
 &'toroidal_tilt           ','toroidal_rotate         ',&
 &'poloidal_rotate         ','minor_radial            '/)
  !> dictionary of short name transformations
  character(len=6), dimension(npdict), parameter :: dictshort = &  !< .
 &(/'scale ','offsca','matoff','offmat','five  ',&
 &'poltil','tortil','torrot','polrot',&
 &'radial' /)
  !> integer dictionary of transformations
  !! 2-digit representations, first denotes coordinate system
  !! 0 - Cartesian, 1 - (R,Z,zeta), 2 - (psi,theta,zeta), 3&4 - (r,theta,zeta)
  !! second denotes
  !! 1 - scale, 2 - scale+offset, 3 - full matrix mult+offset, 4 - offset+full matrix mult
  !! Only solid body transformations (and not five) are actually implemented
  integer(ki4), dimension(npdict), parameter :: dictn = &  !< .
 &(/1,2,3,4,5, &
 &6,7,12,22, &
 &42/)
  logical :: ilmatch !< local variable

  !! position parameters
  namelist /positionparameters/ &
 &position_matrix , &
 &position_vector , &
 &position_scale , &
 &position_offset , &
 &position_transform, &
 &transform_desc, &
 &transform_id


  !! set default position parameters
  ! default identity
  position_transform=0
  transform_desc='cartesian_scale'
  position_vector0=(/1.,0.,0.,0.,1.,0.,0.,0.,1./)
  position_vector=position_vector0
  position_matrix(1,:)=(/1.,  0.,  0. /)
  position_matrix(2,:)=(/0.,  1.,  0. /)
  position_matrix(3,:)=(/0.,  0.,  1. /)
  position_scale=(/1.,  1.,  1. /)
  position_offset=(/0.,  0.,  0. /)
  transform_id='unit'
  ilmatch=.FALSE.

  !!read position parameters
  read(kin,nml=positionparameters,iostat=status)
  if (present(flag)) flag=status
  if(status/=0) then
     if (present(flag)) return
     print '("Fatal error reading position parameters")'
     call log_getunit(ilog)
     write(ilog, nml=positionparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading position parameters')
  end if


  !! check for valid data
  if(all(position_scale==0)) &
 &call log_error(m_name,s_name,2,error_fatal,'position_scale must be > 0')
  if(position_transform<1) then
     !! try for a dictionary match
     do j=1,npdict
        ilmatch=(transform_desc==trim(dictshort(j))).OR.&
 &      (transform_desc==trim(dictlong(j)))
        if (ilmatch) then
           position_transform=dictn(j)
           exit
        end if
     end do
     if (.NOT.ilmatch) &
 &   call log_error(m_name,s_name,3,error_fatal,'Position transform not recognised')
  else if(position_transform>99) then
     call log_error(m_name,s_name,4,error_fatal,'Position transform must be between 1 and 99')
  end if

  if (POSITION_OVERRIDE_IFMIF.AND.position_transform==3) then
     ! FZK's transform
     position_matrix(1,:)=(/0.,  0.,  1. /)
     position_matrix(2,:)=(/-1.,  0.,  0. /)
     position_matrix(3,:)=(/0.,  -1.,  0. /)
     position_scale=(/1.,  1.,  1. /)
     ! working in velocity space, position offset can be ignored
     !     position_offset=(/50.,  0.,  -390. /)
     position_offset=(/0.,  0.,  0. /)
  end if
  position_vector0=position_vector-position_vector0
  if (abs(maxval(position_vector0))+abs(minval(position_vector0))>1.e-6_kr8) then
     ! position vector has changed, use to set matrix
     do j=1,3
        do i=1,3
           position_matrix(i,j)=position_vector(i+(j-1)*3)
        end do
     end do
  end if
  !! store values
  ztfmdata%matrix=position_matrix
  ztfmdata%scale=position_scale
  ztfmdata%offset=position_offset
  ztfmdata%ntfm=position_transform
  !ID  ztfmdata%id=transform_id

end subroutine position_readcon
!---------------------------------------------------------------------
!> read position control data
subroutine position_readtfm(tfmdata,kin)

  !! arguments
  type(tfmdata_t), intent(out) :: tfmdata   !< position transform numeric controls
  integer, intent(in) :: kin   !< output channel for tfmdata

  !! local
  character(*), parameter :: s_name='position_readtfm' !< subroutine name

  read(kin,*,iostat=status) ibuf1
  call log_read_check(m_name,s_name,1,status)
  read(kin,*,iostat=status) tfmdata%scale
  call log_read_check(m_name,s_name,2,status)
  read(kin,*,iostat=status)  tfmdata%offset
  call log_read_check(m_name,s_name,3,status)
  read(kin,*,iostat=status)  tfmdata%matrix
  call log_read_check(m_name,s_name,4,status)
  read(kin,*,iostat=status)  tfmdata%ntfm
  call log_read_check(m_name,s_name,5,status)

end subroutine position_readtfm
!---------------------------------------------------------------------
!> set identity position transform
subroutine position_unitfm(tfmdata)

  !! arguments
  type(tfmdata_t), intent(out) :: tfmdata   !< position transform numeric controls

  !! local
  character(*), parameter :: s_name='position_unitfm' !< subroutine name

  tfmdata%scale=1
  tfmdata%offset=0
  tfmdata%matrix(1,:)=(/1,0,0/)
  tfmdata%matrix(2,:)=(/0,1,0/)
  tfmdata%matrix(3,:)=(/0,0,1/)
  tfmdata%ntfm=1

end subroutine position_unitfm
!---------------------------------------------------------------------
!> set identity quantising position transform
subroutine position_unitfmq(tfmdata)

  !! arguments
  type(quantfm_t), intent(out) :: tfmdata   !< position transform numeric controls

  !! local
  character(*), parameter :: s_name='position_unitfmq' !< subroutine name

  tfmdata%hmin=1
  tfmdata%rhmin=1
  tfmdata%offvec=0
  tfmdata%nqtfm=1

end subroutine position_unitfmq
!---------------------------------------------------------------------
!> get position control data
subroutine position_tfmquery(tfmdata,kchar,object)

  !! arguments
  type(tfmdata_t), intent(in) :: tfmdata   !< position transform numeric controls
  character(*), intent(in) :: kchar  !< case, specifies object
  real(kr4), intent(out), dimension(*) :: object  !< output of object

  !! local
  character(*), parameter :: s_name='position_tfmquery' !< subroutine name

  select case (kchar)
  case ('scale')
     object(1:3)=tfmdata%scale
  case ('offset')
     object(1:3)=tfmdata%offset
  end select

end subroutine position_tfmquery
!---------------------------------------------------------------------
!> write position control data
subroutine position_writetfm(tfmdata,kout)

  !! arguments
  type(tfmdata_t), intent(in) :: tfmdata   !< position transform numeric controls
  integer, intent(in) :: kout   !< output channel for tfmdata

  !! local
  character(*), parameter :: s_name='position_writetfm' !< subroutine name

  write(kout,'(a)',iostat=status) 'TFMDA'
  call log_write_check(m_name,s_name,1,status)
  write(kout,*,iostat=status) tfmdata%scale
  call log_write_check(m_name,s_name,2,status)
  write(kout,*,iostat=status)  tfmdata%offset
  call log_write_check(m_name,s_name,3,status)
  write(kout,*,iostat=status)  tfmdata%matrix
  call log_write_check(m_name,s_name,4,status)
  write(kout,*,iostat=status)  tfmdata%ntfm
  call log_write_check(m_name,s_name,5,status)

end subroutine position_writetfm
!---------------------------------------------------------------------
!> read (vtk) list of position coordinates
subroutine position_readlis(self,infile,channel,kopt)
  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  character(*),intent(in) :: infile !< name of input file
  integer, intent(out) :: channel   !< input channel for position list data
  integer(ki4), intent(in), optional :: kopt   !< options

  !! local
  character(*), parameter :: s_name='position_readlis' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer(ki4) :: inpos    !< local variable

  logical :: isnumb !< local variable
  external isnumb

  if(present(kopt)) then
     !! assume unit already open and reading infile
     nin=channel
  else

     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then channel=i exit end if end do nin=channel

     !! open file
     call misc_getfileunit(nin)
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        call log_error(m_name,s_name,1,error_fatal,'Error opening position list data file')
     else
        call log_error(m_name,s_name,1,log_info,'Position list data file opened')
     end if
     channel=nin

  end if
  !! File unit sorted out
  !! read header information
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,1,error_fatal,'Error reading header data')
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:6)=='POINTS') then
           iltest=isnumb(ibuf2,inpos,7)
           self%np=inpos
           exit
        end if
     end if
  end do

  !! allocate position storage
  if(self%np>0) then
     allocate(self%pos(self%np), &
 &   stat=status)
     !! check successful allocation
     if(status/=0) then
        call log_error(m_name,s_name,2,error_fatal,'Unable to allocate memory')
     end if
  else
     call log_error(m_name,s_name,3,error_fatal,'No position data')
  end if

  !! read coordinates
  do j=1,self%np
     call position_readv(self%pos(j),nin)
  end do
  print '("number of position coordinates read = ",i10)',self%np
  call log_value("number of position coordinates read ",self%np)

end subroutine position_readlis
!---------------------------------------------------------------------
!> copy list of position coordinates
subroutine position_copylis(selfin,selfout,kopt)
  !! arguments
  type(posveclis_t), intent(in) :: selfin !< position list data
  type(posveclis_t), intent(inout) :: selfout !< position list data
  integer(ki4), intent(in), optional :: kopt   !< options

  !! local
  character(*), parameter :: s_name='position_copylis' !< subroutine name

  ! need to check allocation of selfout (TODO)
  !! copy coordinates
  selfout%np=selfin%np
  do j=1,selfout%np
     selfout%pos(j)%posvec=selfin%pos(j)%posvec
  end do
  selfout%nparpos=selfin%nparpos
  selfout%rparpos=selfin%rparpos

end subroutine position_copylis
!---------------------------------------------------------------------
!> read (general vtk format) only list of position coordinates
subroutine position_readonlylis(self,kin,kfmt)
  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  integer, intent(in) :: kin   !< input channel for position list data
  integer(ki4), intent(in) :: kfmt   !< options

  !! local
  character(*), parameter :: s_name='position_readonlylis' !< subroutine name
  real(kr4), dimension(12):: positions  !< transformation offset 3-vector
  integer(ki4) :: inp   !< number of 3-vectors in line

  logical :: isnumb !< local variable
  external isnumb

  !! assume unit already open, reading infile at start of position vectors

  if (kfmt<=1) then
     !! read coordinates one point/line
     do j=1,self%np
        call position_readv(self%pos(j),kin)
     end do
  else
     i=0
     do
        if (i>=self%np) go to 10
        inp=min(kfmt,self%np-i)
        read(kin,*,iostat=status) (positions(k),k=1,3*inp)
        call log_read_check(m_name,s_name,1,status)
        l=1
        do j=1,inp
           self%pos(i+j)%posvec=positions(l:l+2)
           l=l+3
        end do
        i=i+kfmt
     end do
10       continue
  end if
  ! write(*,'(9G12.5)') (self%pos(j)%posvec, j=1,self%np)

end subroutine position_readonlylis
!---------------------------------------------------------------------
!> read (vtk) list of vectors
subroutine position_readveclis(self,infile,kcname,kin,kfmt,kopt)
  !! arguments
  type(posveclis_t), intent(inout) :: self !< vector list data
  character(*),intent(in) :: infile !< name of input file
  character(*),intent(in) :: kcname !< name of field required
  integer, intent(inout) :: kin   !< input channel for vector list data
  integer(ki4), intent(in) :: kfmt   !< numeric data format in file
  integer(ki4), intent(in), optional :: kopt   !< options

  !! local
  character(*), parameter :: s_name='position_readveclis' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer(ki4) :: invec   !< number of vectors
  character(len=80) :: vname !< name of vector
  integer(ki4) :: islen   !< length of vector field name
  integer(ki4) :: islen2   !< length of required vector field name

  logical :: isnumb !< local variable
  external isnumb

  ibuf1=adjustl(kcname)
  islen2=max(2,scan(ibuf1,' '))-1
  if(present(kopt)) then
     !! assume unit already open and reading infile
     if (kin==0) then
        inquire(file=infile,number=kin,iostat=status)
        if(status/=0.OR.kin==-1)then
           !! error opening file
           call log_error(m_name,s_name,1,error_fatal,'Error opening vector list data file')
        else
           call log_error(m_name,s_name,1,log_info,'Vector list data file opened')
           nin=kin
        end if
     else
        nin=kin
     end if
  else

     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kin=i exit end if end do nin=kin

     !! open file
     call misc_getfileunit(nin)
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        call log_error(m_name,s_name,2,error_fatal,'Error opening vector list data file')
     else
        call log_error(m_name,s_name,2,log_info,'Vector list data file opened')
     end if
     kin=nin

  end if

  !! File unit now sorted, get to where point data begin
  if (self%np==0) then
     !! read local header information to get number of vectors
     do
        read(nin,fmt='(a)',iostat=status) ibuf1
        !!eof
        if(status<0) then
           exit
           !! error
        else if (status>0) then
           call log_error(m_name,s_name,3,error_fatal,'Error reading header data')
        else
           ibuf2=adjustl(ibuf1)
           if(ibuf2(1:10)=='POINT_DATA') then
              iltest=isnumb(ibuf2,invec,11)
              self%np=invec
              exit
           end if
        end if
     end do
  end if

  !! find vector header
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,4,error_fatal,'Error reading header data')
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:7)=='VECTORS') then
           ibuf1=adjustl(ibuf2(8:))
           islen=max(2,scan(ibuf1,' '))-1
           vname=ibuf1(:islen)
           if (vname(:islen)==kcname(:islen2)) then
              call log_value("Vector field found ",vname)
              exit
           else
              call log_value("Skipped vector field ",vname)
           end if
           !         else if(ibuf2(1:10)=='POINT_DATA') then
           !            iltest=isnumb(ibuf2,invec,11)
           !            self%np=invec
        end if
     end if
  end do

  !! allocate vector storage
  if(self%np>0) then
     allocate(self%pos(self%np),  stat=status)
     !! check successful allocation
     if(status/=0) then
        call log_error(m_name,s_name,5,error_fatal,'Unable to allocate memory')
     end if
  else
     call log_error(m_name,s_name,6,error_fatal,'No vector data')
  end if

  !! read coordinates
  call position_readonlylis(self,nin,kfmt)
  print '("number of vectors read = ",i10)',self%np
  call log_value("number of vectors read ",self%np)

end subroutine position_readveclis
!---------------------------------------------------------------------
!> delete list of vectors
subroutine position_deleteveclis(self)
  !! arguments
  type(posveclis_t), intent(inout) :: self !< vector list data

  !! local
  character(*), parameter :: s_name='position_deleteveclis' !< subroutine name

  if (allocated(self%pos)) deallocate(self%pos)

end subroutine position_deleteveclis
!---------------------------------------------------------------------
!> write (vtk)  list of position coordinates
subroutine position_writelis(self,kchar,kplot)

  !! arguments
  type(posveclis_t), intent(in) :: self !< position list data
  character(*), intent(in) :: kchar  !< case, specifies objects
  integer, intent(in) :: kplot   !< output channel for position list data
  !     type(beq_t), intent(inout), optional :: pbeq   !> beq data structure

  !! local
  character(*), parameter :: s_name='position_writelis' !< subroutine name
  integer(ki4) :: iseg   !< number of segments in track

  !      if(present(pbeq)) then
  !! convert positions to cartesians
  !      end if

  !! plot list of all positions
  write(kplot,'(''DATASET UNSTRUCTURED_GRID'')')
  write(kplot,'(''POINTS '',I8, '' float'')') self%np
  do j=1,self%np
     call position_writev(self%pos(j),kplot)
  end do

  write(kplot, '('' '')')
  !      if(present(pbeq)) then
  !! convert positions back to psi-theta-zeta
  !      end if

  plot_type: select case (kchar)
  case('track')
     ! output positions as a track
     iseg=self%np-1
     write(kplot,'(''CELLS '',I8,1X,I8)') iseg, 3*iseg
     do i=1,iseg
        write(kplot,'(''2 '',I8,1X,I8)')  i-1,i
     end do
     write(kplot, '('' '')')
     write(kplot,'(''CELL_TYPES '',I8)') iseg
     do i=1,iseg
        write(kplot,'(''3'')')
     end do
     write(kplot, '('' '')')

  case default
     ! output positions as points
     write(kplot,'(''CELLS 1 '',I8)') self%np+1
     write(kplot,'(I8)') self%np,  (i-1,i=1,self%np)
     write(kplot,'(''CELL_TYPES 1'')')
     write(kplot,'(''2'')')
     write(kplot, '('' '')')

  end select plot_type

end subroutine position_writelis
!---------------------------------------------------------------------
!> physical length of list of position coordinates
subroutine position_lenlis(self,length)

  !! arguments
  type(posveclis_t), intent(in) :: self !< position list data
  real(kr4), intent(out) :: length   !< physical length

  !! local
  character(*), parameter :: s_name='position_lenlis' !< subroutine name
  real(kr8) :: zlensum   !< sum of physical length (working)
  real(kr8) :: zlen   !< physical length (working)
  real(kr8), dimension(3)  :: zlenvec !< physical length vector (working)


  !! sum
  zlensum=0
  do j=1,self%np-1
     zlenvec=self%pos(j+1)%posvec-self%pos(j)%posvec
     zlen=sqrt(max(0._kr8,zlenvec(1)**2+zlenvec(2)**2+zlenvec(3)**2))
     zlensum=zlensum+zlen
  end do

  length=zlensum

end subroutine position_lenlis
!---------------------------------------------------------------------
!> transform list of position coordinates
subroutine position_tfmlis(self,tfmdata)

  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  type(tfmdata_t), intent(in) :: tfmdata   !< data defining transform

  !! local
  character(*), parameter :: s_name='position_tfmlis' !< subroutine name
  type(posvecl_t) :: zpos !< local variable

  !! transform list of all positions
  do j=1,self%np
     zpos=position_tfm(self%pos(j),tfmdata)
     self%pos(j)=zpos
  end do
  self%nparpos(4)=self%nparpos(4)+1

end subroutine position_tfmlis
!---------------------------------------------------------------------
!> inverse transform list of position coordinates
subroutine position_invtfmlis(self,tfmdata)

  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  type(tfmdata_t), intent(in) :: tfmdata   !< data defining transform

  !! local
  character(*), parameter :: s_name='position_invtfmlis' !< subroutine name
  type(posvecl_t) :: zpos !< local variable


  !! transform list of all positions
  do j=1,self%np
     zpos=position_invtfm(self%pos(j),tfmdata)
     self%pos(j)=zpos
  end do
  self%nparpos(4)=self%nparpos(4)-1

end subroutine position_invtfmlis
!---------------------------------------------------------------------
!> transform quantised list of position coordinates
subroutine position_qtfmlis(self,qtfmdata)

  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  type(quantfm_t), intent(in) :: qtfmdata   !< data defining transform

  !! local
  character(*), parameter :: s_name='position_qtfmlis' !< subroutine name
  !      type(posvecl_t) position_qtfm
  type(posvecl_t) :: zpos !< local variable


  !! transform list of all positions
  do j=1,self%np
     zpos=position_qtfm(self%pos(j),qtfmdata)
     self%pos(j)=zpos
  end do
  self%nparpos(3)=self%nparpos(3)+1

end subroutine position_qtfmlis
!---------------------------------------------------------------------
!> inverse transform quantised list of position coordinates
subroutine position_invqtfmlis(self,qtfmdata)

  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  type(quantfm_t), intent(in) :: qtfmdata   !< data defining transform

  !! local
  character(*), parameter :: s_name='position_invqtfmlis' !< subroutine name
  !      type(posvecl_t) position_invqtfm
  type(posvecl_t) :: zpos !< local variable


  !! transform list of all positions
  do j=1,self%np
     zpos=position_invqtfm(self%pos(j),qtfmdata)
     self%pos(j)=zpos
  end do
  self%nparpos(3)=self%nparpos(3)-1

end subroutine position_invqtfmlis

end module position_m

