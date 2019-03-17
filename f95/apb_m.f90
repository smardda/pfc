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
 &apb_qartwrite, & !< data structure to beqart compatible format and write
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
  real(kr8) :: zdum !< dummy real
  integer :: status   !< status variable
  character(len=80) :: ibuff !< buffer for input/output
  integer(ki4), parameter :: ip1diag=64 !< diagnostic in 1 coordinate
  integer(ki4), parameter :: ip2diag=64 !< diagnostic in 2 coordinate
  integer(ki4), parameter :: ip3diag=32 !< diagnostic in 3 coordinate
  integer(ki4), parameter :: ip1=6 ! =80 !< 1-coordinate (R) index  !dbg20
  integer(ki4), parameter :: ip2=77 ! =6 !< 2-coordinate (Z) index  !dbg20
  integer(ki4), parameter :: ipcpt=1 !< diagnostic in component
  integer(ki4), parameter :: ipcpt2=2 !< diagnostic 2 in component
  integer(ki4), parameter :: ipcpt3=3 !< diagnostic 3 in component

  contains
!---------------------------------------------------------------------
!> read in special matlab format
subroutine apb_readmlab(self,numerics,infile)

  !! arguments
  type(apb_t), intent(inout) :: self   !< object data structure
  type(mnumerics_t), intent(in) :: numerics   !< object data structure
  character(*),intent(in) :: infile !< name of input file

  !! local
  character(*), parameter :: s_name='apb_readmlab' !< subroutine name
  character(13), parameter :: cfmt0='(1x,f14.7)'
  character(13), parameter :: cfmt1='(5(f14.7,1x))'
  !! logical :: unitused !< flag to test unit is available
  integer :: iin   !< input channel for object data structure

  real(kr8), dimension(:,:), allocatable :: workv2 !< vector work array
  real(kr8), dimension(:,:,:), allocatable :: workm2 !< vector work array
  real(kr8) :: zalph !< rotation angle (clockwise), from n3 array entries
  real(kr8) :: zdalph !< mean angular separation between, from n3 array entries
  real(kr8) :: zcosa !< cosine of rotation angle
  real(kr8) :: zsina !< sine of rotation angle
  real(kr8) :: zbx !< field components
  real(kr8) :: zby !< field components
  real(kr8) :: zbz !< field components
  integer(ki4) :: in3 !< compressed number of points in 3 direction
  integer(ki4) :: inp !< number of points per period
  integer(ki4) :: inj !< angular period selector
  integer(ki4) :: ioff !< offset within cell
  integer(ki4) :: jj !< local variable
  integer(ki4) :: jseg !< segment loop variable
  integer(ki4) :: iseg !< segment offset

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
  !! default layout
  ! number of spaces in pseudo-period, ie each inter-coil toroidal interval
  inp=self%n3/numerics%nzetap-1
  ! number of spaces between points
  in3=self%n3-1
  ioff=0
  layout_type: select case (numerics%layout)
  case(1)
     call log_error(m_name,s_name,numerics%layout,log_info,'Assuming Geoffs periodic data format')
     zalph=self%pos3(self%n3)
  case(2,3,12)
     call log_error(m_name,s_name,numerics%layout,log_info,'Assuming ITER periodic data format')
     zalph=abs(self%pos3(self%n3)-self%pos3(1))
  case(23,33,34,42)
     call log_error(m_name,s_name,numerics%layout,log_info,'Assuming ITER special global data format')
     zalph=0
     ioff=inp/2
     in3=numerics%nzetap*inp
  case default
     zalph=0
  end select layout_type
  zcosa=cos(zalph)
  zsina=sin(zalph)

  !              Allocate 2D storage for read
  allocate(workv2(self%n3,self%ncpt), stat=status)
  call log_alloc_check(m_name,s_name,50,status)
  !              Allocate 3D storage and read
  allocate(self%field(in3+1,self%ncpt,self%n1,self%n2), stat=status)
  call log_alloc_check(m_name,s_name,51,status)

  ij=0
  if (numerics%layout==33) then
  ! stacked segments
     iseg=inp+2
     loopseg: do jseg=1,numerics%nzetap
        do l=1,self%n2
           loopk0: do k=1,self%n1

              do j=1,inp+1
                 read(iin,*,iostat=status)(workv2(j,i),i=1,self%ncpt)
                 if(status/=0) then
                    call log_error(m_name,s_name,52,error_fatal,'Error reading field values')
                 end if
                 ij=ij+self%ncpt
              end do
              !dbgper write(*,*) 'per. test',self%field(1,3,k,l),self%field(self%n3,3,k,l) ! dbgper
              ! loop ignores first entry in workv2, because repeated
              do j=1,inp
                 inj=iseg-j
                 jj=mod( inj+in3-ioff-1,in3 ) + 1
                 !dbg19 if (k==1.AND.l==1) write(*,*) 'j,jj',j,jj  !dbg19
                 zbx=workv2(j,1)
                 zby=workv2(j,2)
                 zbz=workv2(j,3)
                 self%field(jj,1,k,l)=zcosa*zbx+zsina*zby
                 self%field(jj,2,k,l)=-zsina*zbx+zcosa*zby
                 self%field(jj,3,k,l)=zbz
              end do

           end do loopk0
        end do
        iseg=iseg+inp
     end do loopseg
     ! wrap last points
     self%field(in3+1,:,:,:)=self%field(1,:,:,:)
     deallocate(workv2)

  else if (numerics%layout==34) then
  ! stacked segments
     iseg=0
     loopseg4: do jseg=1,numerics%nzetap
        do l=1,self%n2
           loopk4: do k=1,self%n1

              do j=1,inp+1
                 read(iin,*,iostat=status)(workv2(j,i),i=1,self%ncpt)
                 if(status/=0) then
                    call log_error(m_name,s_name,53,error_fatal,'Error reading field values')
                 end if
                 ij=ij+self%ncpt
              end do
              !dbgper write(*,*) 'per. test',self%field(1,3,k,l),self%field(self%n3,3,k,l) ! dbgper
              ! loop ignores last entry in workv2, because repeated
              do j=1,inp
                 inj=iseg+j
                 jj=mod( inj+in3-ioff-1,in3 ) + 1
                 !dbg19 if (k==1.AND.l==1) write(*,*) 'j,jj',j,jj  !dbg19
                 zbx=workv2(j,1)
                 zby=workv2(j,2)
                 zbz=workv2(j,3)
                 self%field(jj,1,k,l)=zcosa*zbx+zsina*zby
                 self%field(jj,2,k,l)=-zsina*zbx+zcosa*zby
                 self%field(jj,3,k,l)=zbz
              end do

           end do loopk4
        end do
        iseg=iseg-inp
     end do loopseg4
     ! wrap last points
     self%field(in3+1,:,:,:)=self%field(1,:,:,:)
     deallocate(workv2)

  else if (numerics%layout==42) then
     ! assume j,k index ordering reversed
     allocate(workm2(self%n1,self%n3,self%ncpt), stat=status)
     call log_alloc_check(m_name,s_name,53,status)

     do l=1,self%n2

        do j=1,self%n3
           do k=1,self%n1
              read(iin,*,iostat=status)(workm2(k,j,i),i=1,self%ncpt)
              if(status/=0) then
                 call log_error(m_name,s_name,11,error_fatal,'Error reading field values')
              end if
              ij=ij+self%ncpt
           end do
        end do
        ! unravelling stacked periodic cases
        do j=1,in3
           inj=(j-2+ioff)/inp
           jj=mod(j+self%n3+ioff+inj-1,self%n3)+1
           !dbg19 write(*,*) 'j,jj',j,jj  !dbg19
           do k=1,self%n1
              zbx=workm2(k,jj,1)
              zby=workm2(k,jj,2)
              zbz=workm2(k,jj,3)
              self%field(j,1,k,l)=zcosa*zbx+zsina*zby
              self%field(j,2,k,l)=-zsina*zbx+zcosa*zby
              self%field(j,3,k,l)=zbz
           end do
        end do
        ! wrapping point
        do k=1,self%n1
           self%field(in3+1,:,k,l)=self%field(1,:,k,l)
        end do

     end do
     deallocate(workm2)
  else
     do l=1,self%n2
        loopk: do k=1,self%n1

           do j=1,self%n3
              read(iin,*,iostat=status)(workv2(j,i),i=1,self%ncpt)
              if(status/=0) then
                 call log_error(m_name,s_name,12,error_fatal,'Error reading field values')
              end if
              ij=ij+self%ncpt
           end do
           !dbgper  write(*,*) 'per. test',self%field(1,3,k,l),self%field(self%n3,3,k,l)! dbgper

           layout_type2: select case (numerics%layout)
           case(1,2,12)
              ! assumes data input on phi periodic cell, not zeta cell (phi=-zeta)
              do j=1,self%n3
                 zbx=workv2(j,1)
                 zby=workv2(j,2)
                 zbz=workv2(j,3)
                 self%field(self%n3+1-j,1,k,l)=zcosa*zbx+zsina*zby
                 self%field(self%n3+1-j,2,k,l)=-zsina*zbx+zcosa*zby
                 self%field(self%n3+1-j,3,k,l)=zbz
              end do
           case(3)
              ! another option for unravelling periodic cases
              do j=1,self%n3
                 zbx=workv2(j,1)
                 zby=workv2(j,2)
                 zbz=workv2(j,3)
                 self%field(j,1,k,l)=zcosa*zbx+zsina*zby
                 self%field(j,2,k,l)=-zsina*zbx+zcosa*zby
                 self%field(j,3,k,l)=zbz
              end do
           case(23)
              ! option for unravelling stacked periodic cases
              do j=1,in3
                 inj=(j-2+ioff)/inp
                 jj=mod(j+self%n3+ioff+inj-1,self%n3)+1
                 !dbg19 write(*,*) 'j,jj',j,jj  !dbg19
                 zbx=workv2(jj,1)
                 zby=workv2(jj,2)
                 zbz=workv2(jj,3)
                 self%field(j,1,k,l)=zcosa*zbx+zsina*zby
                 self%field(j,2,k,l)=-zsina*zbx+zcosa*zby
                 self%field(j,3,k,l)=zbz
              end do
              ! wrapping point
              self%field(in3+1,:,k,l)=self%field(1,:,k,l)
           case default
              self%field(:,:,k,l)=workv2
           end select layout_type2

        end do loopk
     end do
     deallocate(workv2)
  end if
  print '("number of field values read = ",i10)',ij
  call log_value("number of field values read ",ij)

  !! rearrange entries in pos3
  layout_type3: select case (numerics%layout)
  case(2,3,12)
     zdalph=abs(self%pos3(self%n3)-self%pos3(1))/in3
     do j=1,in3+1
        self%pos3(j)=(j-1)*zdalph
     end do
  case(23,33,34,42)
     zdalph=abs(self%pos3(self%n3)-self%pos3(1))/in3
     do j=1,in3+1
        self%pos3(j)=(j-1)*zdalph
     end do
     self%n3=in3+1
  end select layout_type3
  !d      if (ip1diag>0.AND.ip2diag>0) then
  !d!  diagnostics
  !d      write(*,*) 'ipdiag',(self%field(j,ipcpt,ip1diag,ip2diag),j=1,ip3diag)
  !d      end if
  !dbg20 if (ipcpt3==3) then !dbg20
  !dbg20 write(203,*) '#', 'ipcpt3,ip1,ip2=',ipcpt3,ip1,ip2 !dbg20
  !dbg20 write(203,'(3(1x,1pg14.7))') ((self%field(j,l,ip1,ip2),l=1,self%ncpt),j=1,self%n3) !dbg20
  !dbg20 write(203,*) '#---------------#' !dbg20
  !dbg20 else !dbg20
  !dbg20 write(200,*) '#', 'ipcpt,ip1,ip2=',ipcpt,ip1,ip2 !dbg20
  !dbg20 write(200,'((1x,1pg14.7))') (self%field(j,ipcpt,ip1,ip2), j=1,self%n3) !dbg20
  !dbg20 write(200,*) '#---------------#' !dbg20
  !dbg20 end if !dbg20
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
        !dbg  write(*,*) 'l,k,zrerr,zzerr ',l,k,zrerr,zzerr,abs(zrphi0-zrphin),abs(zzphi0-zzphin)  !dbg
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
  !! contents
  selfcopy%pos1=self%pos1
  selfcopy%pos2=self%pos2
  selfcopy%pos3=self%pos3
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
!> data structure to beqart compatible format and write
subroutine apb_qartwrite(self,numerics,kout)

  !! arguments
  type(apb_t), intent(in) :: self   !< object data structure
  type(mnumerics_t), intent(in) :: numerics   !< object data structure
  integer, intent(in) :: kout   !< output channel

  !! local
  character(*), parameter :: s_name='apb_qartwrite' !< subroutine name
  integer(ki4) :: in3 !< compressed number of points in 3 direction
  integer(ki4) :: inp !< number of points per period
  integer(ki4) :: ioff !< offset within cell
  integer(ki4) :: jj !< local variable
  real(kr8) :: zalph !< toroidal angle extent from n3 array entries

  in3=self%n3-1
  ioff=in3/2
  zalph=abs(self%pos3(self%n3)-self%pos3(1))
  zdum=0
  ! pseudo namelist output
  write(kout,'(A)') '&fmeshparameters'
  write(kout,'(A)') '  mesh_transform=.FALSE., mesh_formula=''limits'''
  write( kout,'(a8,1pg14.7,a6,g14.7,a4,i6)' ) &
  ' xstart=',self%pos1(1),',xend=',self%pos1(self%n1),',nx=',self%n1
  write( kout,'(a8,1pg14.7,a6,g14.7,a4,i6)' ) &
  ' ystart=',self%pos2(1),',yend=',self%pos2(self%n2),',ny=',self%n2
  write( kout,'(a8,1pg14.7,a6,g14.7,a4,i6)' ) &
  ' zstart=',-0.5*zalph,',zend=',0.5*zalph,',nz=',self%n3
  write(kout,'(A)') '/'

  write(kout,'(i6)') self%n1
  write(kout,'(i6)') self%n2
  write(kout,'(i6)') self%n3
  do j=1,self%n3
     ! map from (0,2pi) to (-pi,pi)
     jj=mod(j+ioff-1,in3)+1
     !dbg19 write(*,*) 'j,jj=',j,jj !dbg19
     do l=1,self%n2
        do k=1,self%n1
           write(kout,'(3(1x,1pg14.7))') (self%field(jj,i,k,l), i=1,self%ncpt)
        end do
     end do
  end do

   !dbg20 write(201,*) '#', 'ipcpt2,ip1,ip2=',ipcpt2,ip1,ip2 !dbg20
   !dbg20 write(201,'((1x,1pg14.7))') (self%field(jj,ipcpt2,ip1,ip2), jj=1,self%n3) !dbg20

end subroutine apb_qartwrite
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
        !DD     write(*,*) workv2(:,1) !DD
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

        !dbgper  write(*,*) 'approx per. test',workv(1,1),workv(itfm,1)! dbgper
        !dbgper  write(*,*) 'approx per. test',workv(1,2),workv(itfm,2)! dbgper
        !dbgper  write(*,*) 'approx per. test',workv(1,3),workv(itfm,3)! dbgper
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
