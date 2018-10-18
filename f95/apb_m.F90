module apb_m

  use const_kind_m
!      use const_numphys_h
  use log_m
  use misc_m
  use mcontrol_h
  use mcontrol_m
  use date_time_m
  use apb_h

  implicit none
  private

! public subroutines
  public :: &
 &apb_readmlab, & !< read in special matlab format
 &apb_checkperiod, & !< check \f$ \phi \f$ periodicity of  data structure
 &apb_copy, & !< copy data structure
 &apb_scale, & !< scale data structure
 &apb_convert !< Fourier transform in \f$ \zeta \f$

! public types
  integer :: nin !< input file unit

! private subroutines
  private :: &
 &apb_fileafter !< move to point in file after line beginning with string

! private variables
  character(*), parameter :: m_name='apb_m' !< module name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: idum !< dummy integer
  integer :: status   !< status variable
  character(len=80) :: ibuff !< buffer for input/output
  integer(ki4), parameter :: ip1diag=64 !< diagnostic in 1 coordinate
  integer(ki4), parameter :: ip2diag=64 !< diagnostic in 2 coordinate
  integer(ki4), parameter :: ip3diag=32 !< diagnostic in 3 coordinate
  integer(ki4), parameter :: ipcpt=3 !< diagnostic in component
  integer(ki4), parameter :: ipcpt2=2 !< diagnostic 2 in component

  contains
!---------------------------------------------------------------------
!> read in special matlab format
subroutine apb_readmlab(self,infile)

  !! arguments
  type(apb_t), intent(inout) :: self   !< object data structure
  character(*),intent(in) :: infile !< name of input file

  !! local
  character(*), parameter :: s_name='apb_readmlab' !< subroutine name
  character(13), parameter :: cfmt0='(1x,f14.7)'
  character(13), parameter :: cfmt1='(5(f14.7,1x))'
  !! logical :: unitused !< flag to test unit is available
  integer :: iin   !< input channel for object data structure

  integer(ki4), parameter :: geoffdata=0 !< flag for Geoff's data format
  real(kr8), dimension(:,:), allocatable :: workv2 !< vector work array
  real(kr8) :: zalph !< rotation angle (clockwise)
  real(kr8) :: zcosa !< cosine of rotation angle
  real(kr8) :: zsina !< sine of rotation angle
  real(kr8) :: zbx !< field components
  real(kr8) :: zby !< field components
  real(kr8) :: zbz !< field components

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then iin=i exit end if end do

  call misc_getfileunit(iin)
  open(unit=iin, file=infile,status='OLD',form='FORMATTED',iostat=status)
  !!eof or error
  if (status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error opening file')
  end if
  ! skip header data
  call apb_fileafter('Zeta grid',iin)

  read(iin,*,iostat=status) self%n3
  if(status/=0) then
     call log_error(m_name,s_name,11,error_fatal,'Error reading object data')
  end if
  ! position 3 array
  !! allocate position 3 storage
  if(self%n3>0) then
     allocate(self%pos3(self%n3), stat=status)
     call log_alloc_check(m_name,s_name,12,status)
  else
     call log_error(m_name,s_name,13,error_fatal,'No 1D data')
  end if

  read(iin,*,iostat=status)(self%pos3(i),i=1,self%n3)
  call log_alloc_check(m_name,s_name,14,status)
  print '("number of zeta values read = ",i10)',self%n3
  call log_value("number of zeta values read ",self%n3)

  call apb_fileafter('R grid',iin)

  read(iin,*,iostat=status) self%n1
  call log_read_check(m_name,s_name,21,status)
  ! position 1 array
  !! allocate position 1 storage
  if(self%n1>0) then
     allocate(self%pos1(self%n1), stat=status)
     call log_alloc_check(m_name,s_name,22,status)
  else
     call log_error(m_name,s_name,23,error_fatal,'No 1D data')
  end if

  read(iin,*,iostat=status)(self%pos1(i),i=1,self%n1)
  call log_alloc_check(m_name,s_name,24,status)
  print '("number of R values read = ",i10)',self%n1
  call log_value("number of R values read ",self%n1)

  call apb_fileafter('Z grid',iin)

  read(iin,*,iostat=status) self%n2
  call log_read_check(m_name,s_name,33,status)
  ! position 2 array
  !! allocate position 2 storage
  if(self%n2>0) then
     allocate(self%pos2(self%n2), stat=status)
     call log_alloc_check(m_name,s_name,32,status)
  else
     call log_error(m_name,s_name,33,error_fatal,'No 1D data')
  end if

  read(iin,*,iostat=status)(self%pos2(i),i=1,self%n2)
  call log_alloc_check(m_name,s_name,34,status)
  print '("number of Z values read = ",i10)',self%n2
  call log_value("number of Z values read ",self%n2)

  !!     write(*,*)(self%pos1(i),i=1,self%n1)
  !!     write(*,*)(self%pos2(i),i=1,self%n2)
  !!     write(*,*)(self%pos3(i),i=1,self%n3)

  ! skip blank and non-blank
  read(iin,fmt='(a)',iostat=status) ibuff
  read(iin,fmt='(a)',iostat=status) ibuff

  !! number of components
  self%ncpt=3
  !-----------------------------------------------------------------------
  ! fixup
  if (geoffdata/=0) then
     call log_error(m_name,s_name,1,log_info,'Assuming Geoffs data format')
     zalph=self%pos3(self%n3)
     zcosa=cos(zalph)
     zsina=sin(zalph)
  end if

  !              Allocate 2D storage for read
  allocate(workv2(self%n3,self%ncpt), stat=status)
  call log_alloc_check(m_name,s_name,50,status)
  !              Allocate 3D storage and read
  allocate(self%field(self%n3,self%ncpt,self%n1,self%n2), stat=status)
  call log_alloc_check(m_name,s_name,51,status)

  ij=0
  do l=1,self%n2
     loopk: do k=1,self%n1

        do j=1,self%n3
           read(iin,*,iostat=status)(workv2(j,i),i=1,self%ncpt)
           if(status/=0) then
              call log_error(m_name,s_name,12,error_fatal,'Error reading field values')
           end if
           ij=ij+self%ncpt
        end do
        !     write(*,*) 'per. test',self%field(1,3,k,l),self%field(self%n3,3,k,l)
        self%field(:,:,k,l)=workv2

        if (geoffdata==1) then
           ! transformations it would be nice not to have to do
           do j=1,self%n3
              zbx=workv2(j,1)
              zby=workv2(j,2)
              zbz=workv2(j,3)
              self%field(self%n3+1-j,1,k,l)=zcosa*zbx+zsina*zby
              self%field(self%n3+1-j,2,k,l)=-zsina*zbx+zcosa*zby
              self%field(self%n3+1-j,3,k,l)=zbz
           end do
        else if (geoffdata==1) then
           ! transformations it would be nice not to have to do
           do j=1,self%n3
              zbx=workv2(j,1)
              zby=workv2(j,2)
              zbz=workv2(j,3)
              self%field(j,1,k,l)=zcosa*zbx+zsina*zby
              self%field(j,2,k,l)=-zsina*zbx+zcosa*zby
              self%field(j,3,k,l)=zbz
           end do
        end if

     end do loopk
  end do
  print '("number of field values read = ",i10)',ij
  call log_value("number of field values read ",ij)

  !d      if (ip1diag>0.AND.ip2diag>0) then
  !d!  diagnostics
  !d      write(*,*) 'ipdiag',(self%field(j,ipcpt,ip1diag,ip2diag),j=1,ip3diag)
  !d      end if
  close(iin)

end subroutine apb_readmlab
!---------------------------------------------------------------------
!> check phi periodicity of  input
subroutine apb_checkperiod(self)

  !! arguments
  type(apb_t), intent(in) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='apb_checkperiod' !< subroutine name
  real(kr8) :: zrphi0    !<  \f$ R \f$ cpt of vector field at \f$ \phi=0 \f$
  real(kr8) :: zrphin    !<  \f$ R \f$ cpt of vector field at \f$ \phi=1 \f$
  real(kr8) :: zrref    !<  reference \f$ R \f$ cpt of vector field
  real(kr8) :: zrerr    !<  relative error in periodicity of \f$ R \f$ cpt of vector field
  real(kr8) :: zrerrmax    !<  maximum relative error in periodicity of \f$ R \f$ cpt of vector field
  real(kr8) :: zzphi0    !<  \f$ Z \f$ cpt of vector field at \f$ \phi=0 \f$
  real(kr8) :: zzphin    !<  \f$ Z \f$ cpt of vector field at \f$ \phi=1 \f$
  real(kr8) :: zzref    !<  reference \f$ Z \f$ cpt of vector field
  real(kr8) :: zzerr    !<  relative error in periodicity of \f$ Z \f$ cpt of vector field
  real(kr8) :: zerrmax    !<  maximum error in periodicity of \f$ Z \f$ cpt of vector field
  real(kr8) :: zzerrmax    !<  maximum relative error in periodicity of \f$ Z \f$ cpt of vector field

  zrerrmax=0
  zzerrmax=0
  zerrmax=0
  !-----------------------------------------------------------------------
  !              Compare
  do l=1,self%n2
     do k=1,self%n1
        zrphi0=sqrt(self%field(1,1,k,l)**2+self%field(1,2,k,l)**2)
        zrphin=sqrt(self%field(self%n3,1,k,l)**2+self%field(self%n3,2,k,l)**2)
        zrref=(abs(zrphi0)+abs(zrphin))/2
        zrerr=(abs(zrphi0-zrphin))/zrref
        zrerrmax=max(zrerrmax,zrerr)
        zrerrmax=max(zrerrmax,zrerr)
        zzphi0=self%field(1,3,k,l)
        zzphin=self%field(self%n3,3,k,l)
        zzref=(1.e-5+abs(zzphi0)+abs(zzphin))/2
        zzerr=(abs(zzphi0-zzphin))/zzref
        !      write(*,*) 'l,k,zrerr,zzerr ',l,k,zrerr,zzerr,abs(zrphi0-zrphin),abs(zzphi0-zzphin)
        zzerrmax=max(zzerrmax,zzerr)
        zerrmax=max(zerrmax,abs(zzphi0-zzphin))
     end do
     !      write(*,*) ' '
  end do

  call log_value("maximum rel error of R-field periodicity ",zrerrmax)
  call log_value("maximum rel error of Z-field periodicity ",zzerrmax)
  call log_value("maximum error of Z-field periodicity ",zerrmax)

end subroutine apb_checkperiod
!---------------------------------------------------------------------
!> copy input
subroutine apb_copy(self,selfcopy)

  !! arguments
  type(apb_t), intent(in) :: self   !< object data structure
  type(apb_t), intent(out) :: selfcopy   !< copyd object data structure

  !! local
  character(*), parameter :: s_name='apb_copy' !< subroutine name

  ! position 1 array
  !! allocate position 1 storage
  if(self%n1>0) then
     allocate(selfcopy%pos1(self%n1), stat=status)
     call log_alloc_check(m_name,s_name,22,status)
  else
     call log_error(m_name,s_name,23,error_fatal,'No 1D data')
  end if

  ! position 2 array
  !! allocate position 2 storage
  if(self%n2>0) then
     allocate(selfcopy%pos2(self%n2), stat=status)
     call log_alloc_check(m_name,s_name,32,status)
  else
     call log_error(m_name,s_name,33,error_fatal,'No 1D data')
  end if

  ! position 3 array
  !! allocate position 3 storage
  if(self%n3>0) then
     allocate(selfcopy%pos3(self%n3), stat=status)
     call log_alloc_check(m_name,s_name,12,status)
  else
     call log_error(m_name,s_name,13,error_fatal,'No 1D data')
  end if

  !! copy dimensions, number of components
  selfcopy%n1=self%n1
  selfcopy%n2=self%n2
  selfcopy%n3=self%n3
  selfcopy%ncpt=self%ncpt
  !-----------------------------------------------------------------------
  !              Allocate 3D storage and read
  if(self%ncpt>0) then
     allocate(selfcopy%field(selfcopy%n3,selfcopy%ncpt,selfcopy%n1,selfcopy%n2), stat=status)
     call log_alloc_check(m_name,s_name,50,status)
  else
     call log_error(m_name,s_name,51,error_fatal,'No components')
  end if
  ! to work on ifort, have to explicitly code copy
  !     selfcopy%field=self%field
  do l=1,self%n2
     do k=1,self%n1
        do i=1,self%ncpt
           do j=1,self%n3
              selfcopy%field(j,i,k,l)=self%field(j,i,k,l)
              !     write(*,*) 'j,i,k,l test',j,i,k,l
           end do
        end do
     end do
  end do

end subroutine apb_copy
!---------------------------------------------------------------------
!> scale field
subroutine apb_scale(self,numerics,pfac)

  !! arguments
  type(apb_t), intent(inout) :: self   !< object data structure
  type(mnumerics_t), intent(in) :: numerics   !< object data structure
  real(kr8) :: pfac    !<  scale factor

  !! local
  character(*), parameter :: s_name='apb_scale' !< subroutine name

  !     self%field=pfac*self%field
  do l=1,self%n2
     do k=1,self%n1
        do i=1,self%ncpt
           do j=1,self%n3
              self%field(j,i,k,l)=pfac*self%field(j,i,k,l)
           end do
        end do
     end do
  end do

end subroutine apb_scale
!---------------------------------------------------------------------
!> Fourier transform in zeta
subroutine apb_convert(self,numerics,kinit,rpower,phchange)

  !! arguments
  type(apb_t), intent(inout) :: self   !< object data structure
  type(mnumerics_t), intent(in) :: numerics   !< object data structure
  integer(ki4), intent(in) :: kinit   !< zero to initialise FFT
  integer(ki4), intent(in) :: rpower !< power of \f$ R \f$ multiplying \f$ \bf{B} \f$
  integer(ki4), intent(in) :: phchange !< change phase of \f$ \zeta \f$

  !! local
  character(*), parameter :: s_name='apb_convert' !< subroutine name
  real(kr8), dimension(:), save, allocatable :: ywork !< work array
  complex(kr8), dimension(:), allocatable :: ctfm !< complex work array
  real(kr8), dimension(:),  allocatable :: rwork !< work array
  real(kr8), dimension(:,:), allocatable :: workv !< vector work array
  real(kr8), dimension(:,:), allocatable :: workv2 !< vector work array
  integer(ki4), save :: itfm   !< length of tfm
  integer(ki4), save :: itfmp1  !< length of tfm + 1
  real(kr8), save :: ztfm  !< real itfm
  real(kr8), save :: znorm  !< 1./ztfm
  integer(ki4), save :: itfmb2  !< itfm/2
  integer(ki4), save :: itfmbp1  !< itfmb2+1
  real(kr8) :: zeta    !<  \f$ \zeta \f$
  real(kr8) :: zetaoff    !<  \f$ \zeta \f$ offset
  real(kr8) :: zcos    !<  \f$ \cos(\zeta) \f$
  real(kr8) :: zsin    !<  \f$ \sin(\zeta) \f$
  real(kr8) :: zr    !<  radius \f$ R \f$
  real(kr8) :: zrp    !<  \f$ R^p \f$
  real(kr8) :: zrr    !<  \f$ R^(p-1) \f$
  real(kr8) :: zfr    !<  radial cpt of vector field
  real(kr8) :: zfze    !<  \f$ \zeta \f$ cpt of vector field
  real(kr8) :: zfz    !<  \f$ z \f$ cpt of vector field

  if (kinit==0) then
     itfm=self%n3-1
     itfmp1=itfm+1
     itfmb2=itfm/2
     itfmbp1=itfmb2+1
     ztfm=itfm
     znorm=1/ztfm
     ! initialise FFT
     allocate(ywork(4*itfm+15), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     call zffti(itfm,ywork)
     ! initialise copy arrays
     allocate(rwork(0:itfm), stat=status)
     call log_alloc_check(m_name,s_name,2,status)
     allocate(workv(itfm,3), stat=status)
     call log_alloc_check(m_name,s_name,3,status)
     allocate(workv2(itfm,3), stat=status)
     call log_alloc_check(m_name,s_name,4,status)
     allocate(ctfm(itfm), stat=status)
     call log_alloc_check(m_name,s_name,5,status)
     ctfm=0

  else if (kinit==2) then
     ! destructor
     deallocate(ywork)
     deallocate(rwork)
     deallocate(ctfm)
     return

  end if


  do l=1,self%n2
     do k=1,self%n1
        zr=self%pos1(k)
        zrp=zr**rpower
        zrr=zrp/zr

        !      write(*,*) 'ends',self%field(1,1,k,l)**2+self%field(1,2,k,l)**2,&
        !     &                  self%field(self%n3,1,k,l)**2+self%field(self%n3,2,k,l)**2
        ! select vector
        iloop0: do i=1,self%ncpt
           do j=1,itfm
              workv(j,i)=self%field(j,i,k,l)
           end do
        end do iloop0
        !DD     write(*,*) workv2(:,1)
        ! convert vector to R,Z,Zeta
        do j=1,itfm
           zeta=self%pos3(j)
           zsin=sin(zeta)
           zcos=cos(zeta)
           zfz=zrp*workv(j,3)
           if (phchange==1) then
              ! phase change of pi in zeta
              zfr=zrp*(workv(j,1)*zcos-workv(j,2)*zsin)
              zfze=zrr*(workv(j,1)*zsin+workv(j,2)*zcos)
           else if (phchange==2) then
              ! do map phi -> zeta
              zfr=zrp*(workv(j,1)*zcos+workv(j,2)*zsin)
              zfze=zrr*(-workv(j,1)*zsin+workv(j,2)*zcos)
           else
              zfr=zrp*(workv(j,1)*zcos-workv(j,2)*zsin)
              zfze=-zrr*(workv(j,1)*zsin+workv(j,2)*zcos)
           end if
           workv(j,1)=zfr
           workv(j,2)=zfz
           workv(j,3)=zfze
        end do

        if (phchange==1) then
           ! phase change of pi in zeta
           iloop1: do i=1,self%ncpt
              do j=1,itfmb2
                 workv2(j,i)=workv(itfmb2+j,i)
                 workv2(itfmb2+j,i)=workv(j,i)
              end do
           end do iloop1
        else if (phchange==2.OR.phchange==4) then
           ! do map phi -> zeta
           iloop12: do i=1,self%ncpt
              do j=1,itfm
                 workv2(j,i)=workv(itfmp1-j,i)
              end do
           end do iloop12
        else
           workv2=workv
        end if
        !DD     write(*,*) workv(:,1)

        !      write(*,*) 'approx per. test',workv(1,1),workv(itfm,1)
        !      write(*,*) 'approx per. test',workv(1,2),workv(itfm,2)
        !      write(*,*) 'approx per. test',workv(1,3),workv(itfm,3)
        !d      if (k==ip1diag.AND.l==ip2diag) then
        !d!  diagnostics
        !d      write(*,*) 'workv',(workv(j,ipcpt2),j=1,ip3diag)
        !d      write(*,*) 'workv2',(workv2(j,ipcpt2),j=1,ip3diag)
        !d      end if

        ! transform
        iloop2: do i=1,self%ncpt
           do j=1,itfm
              ctfm(j)=workv2(j,i)
           end do
           !     write(*,*) 'before',ctfm
           call zfftf(itfm,ctfm,ywork)
           !     write(*,*) 'after',ctfm
           !     put into work into a:b  order and normalise
           !   a  and b entries
           rwork(0)=real(ctfm(1))
           rwork(1)=aimag(ctfm(1))
           ij=1
           do j=2,itfmb2
              ij=ij+1
              rwork(ij)=2*real(ctfm(j))
              ij=ij+1
              rwork(ij)=-2*aimag(ctfm(j))
           end do
           ij=ij+1
           rwork(ij)=real(ctfm(itfmbp1))
           !     write(*,*) 'rwork',rwork
           ! normalise
           rwork=znorm*rwork
           ! rwork back to field array
           self%field(1:,i,k,l)=rwork(0:)
        end do iloop2

        !d      if (k==ip1diag.AND.l==ip2diag) then
        !d!  diagnostics
        !d      write(*,*) 'self%field',(self%field(j,ipcpt2,k,l),j=1,ip3diag)
        !d      end if

     end do
  end do

  if (phchange==1) then
     ! complete phase change of pi in zeta
     zetaoff=-self%pos3(itfmbp1)
     self%pos3=self%pos3+zetaoff
  else if (phchange==2) then
     ! complete map phi -> zeta by doing nothing  to pos3 array
  else if (phchange==3) then
     ! reflect field in zeta direction
     do l=1,self%n2
        do k=1,self%n1
           iloop3: do i=1,self%ncpt
              rwork(0:)=self%field(1:,i,k,l)
              !   negate b entries
              ij=1
              do j=1,itfmb2
                 rwork(ij)=-rwork(ij)
                 ij=ij+2
              end do
              !   negate zeta-component
              if (i==self%ncpt) rwork=-rwork
              ! rwork back to field array
              self%field(1:,i,k,l)=rwork(0:)
           end do iloop3
        end do
     end do
  end if

end subroutine apb_convert
!>---------------------------------------------------------------------
subroutine apb_fileafter(kfind,kin)
  !< move to point in file after line beginning with string
  !! arguments
  character(len=*), intent(in) :: kfind   !< character string to find
  integer, intent(in) :: kin   !< file unit for reading

  !! local
  character(*), parameter :: s_name='apb_fileafter' !< subroutine name
  character(len=80) :: ibuff   !< character string

  do
     read(kin,fmt='(a)',iostat=status) ibuff
     call log_read_check(m_name,s_name,1,status)
     if (adjustl(ibuff)==kfind) exit
  end do

  if (adjustl(ibuff)/=kfind) then
     ! jm string not found
     call log_error(m_name,s_name,10,error_fatal,'Error reading object data')
  end if

end subroutine apb_fileafter

end module apb_m
