module smanal_m

  use log_m
  use misc_m
  use const_numphys_h
  use const_kind_m
  use smanal_h
  use scontrol_h
  use geobj_m
  use geobjlist_h
  use geobjlist_m
  use vfile_m
  use indict_m

  implicit none
  private

! public subroutines
  public :: &
  smanal_read,  & !< read data from file
  smanal_rekey, & !< reassign key
  smanal_stats,  & !< coordinate statistics subroutine
  smanal_initwrite, & !< open new file, making up name
  smanal_write, &  !< write out object
  smanal_writeg, &  !< write out object as gnuplot
  smanal_writev, &  !< write out object as vtk
  smanal_delete, & !< delete object
  smanal_close, & !< close file
  smanal_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='smanal_m' !< module name
  integer  :: status   !< error status
  integer, save  :: ninsm=0     !< geometry file unit number
  integer, save  :: ninscal=0     !< scalx file unit number
  integer, save  :: noutsm=0      !< output file unit number
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: iopt !< option
  integer(ki4), dimension(:), allocatable :: iwork !< work array
  integer(ki4) :: imaxstat   !< number of outputs for smanal data structure
  real(kr8), dimension(:), allocatable :: work !< 1D work array
  type(geobjlist_t)  :: igeobjl      !< geometrical object

  contains
!---------------------------------------------------------------------
!> initialise object with geobjl data
subroutine smanal_read(self,file,numerics)

  !! arguments
  type(smanal_t), intent(out) :: self !< type which data will be assigned to
  type(sfiles_t), intent(in)      :: file      !< names of files
  type(snumerics_t), intent(in)   :: numerics  !< numerical parameters
  !! local
  character(*), parameter :: s_name='smanal_read' !< subroutine name
  character(len=80),save :: iched !< vtk field file descriptor

  self%geobjl%ngtype=2
  if (file%nvtkdata==1) then
     ! just one file with Body data
     call geobjlist_read(self%geobjl,file%vtkdata(1),iched)
     ninsm=0
     iopt=1
     call vfile_iscalarread(self%key,self%nkey,file%vtkdata(1),numerics%namekey,ninsm,iopt)
     call indict_fromlist(self%dict,self%key,self%ndict)
     call indict_sort(self%dict,iwork,self%ndict,1)
     call geobjlist_centroids(self%geobjl,self%key,self%dict,self%ndict,self%centr)
  end if

  !dbg write(*,*) numerics%namekey, numerics%namescal !dbg
  if (trim(adjustl(numerics%namekey))==trim(adjustl(numerics%namescal))) then
     !! scalar is in first file, ie. Body or Surf
     self%nscal=self%nkey
     allocate(self%scal(self%nscal),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     self%scal=self%key
  else
     !! scalar typically pow
     igeobjl%ngtype=2
     call geobjlist_read(igeobjl,file%vtk,iched)
     ninscal=0
     iopt=1
     call vfile_dscalarread(self%scal,self%nscal,file%vtk,numerics%namescal,ninscal,iopt)
     call geobjlist_delete(igeobjl)
  end if

  self%n=numerics

end subroutine smanal_read
!---------------------------------------------------------------------
!> rekey geobjl data based on position
subroutine smanal_rekey(self)

  !! arguments
  type(smanal_t), intent(inout) :: self !< type
  !! local
  character(*), parameter :: s_name='smanal_rekey' !< subroutine name
  real(kr8), dimension(3) :: zcentr !< average of centroids
  real(kr8) :: zrcen    !<   \f$ R_C \f$
  real(kr8) :: zzcen    !<   \f$ Z_C \f$
  real(kr8) :: zr    !<   \f$ R \f$
  real(kr8) :: zx    !<   \f$ X \f$
  real(kr8) :: zy    !<   \f$ Y \f$
  real(kr8) :: zz    !<   \f$ Z-Z_C \f$
  real(kr8) :: zangmin !< minimum angle
  real(kr8) :: zangmax !< maximum angle
  real(kr8) :: zbinsiz !< angle increment
  integer(ki4) :: ikey   !< local variable
  integer(ki4) :: indx !< local variable
  integer(ki4) :: iindict !< local variable

  iindict=self%ndict

  ! allow for user defined central values
  if (.NOT.self%n%lurcen.OR..NOT.self%n%luzcen) then
     ! find centroid unless user is setting both
     zcentr=0
     do i=1,iindict
        zcentr(:)=zcentr(:)+self%centr(:,i)
     end do
     zcentr(:)=zcentr(:)/iindict
  end if
  if (self%n%lurcen) then
     zrcen=self%n%urcen
  else
     zrcen=sqrt( zcentr(1)**2+zcentr(2)**2 )
  end if
  if (self%n%luzcen) then
     zzcen=self%n%uzcen
  else
     zzcen=zcentr(3)
  end if

  call log_value('Nominal R-centre',zrcen)
  call log_value('Nominal Z-centre',zzcen)

  self%inbin=iindict
  new_key: select case (self%n%newkey)

  case('null','angle','poloidal')
     ! new sort key based on poloidal angle
     allocate(self%bin(iindict), self%oldict(iindict), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     do i=1,iindict
        zr=sqrt( self%centr(1,i)**2+self%centr(2,i)**2 )
        zz=self%centr(3,i)-zzcen
        ! zang(i)=atan2(zz,zr-zrcen) might be expected, but to get branch-cut down below
        self%bin(i)=atan2(zr-zrcen,zz)
        !write(*,*) self%dict(i),self%bin(i),self%centr(:,i)
     end do
     self%oldict=self%dict
     ! bin size
     zangmin=minval(self%bin)
     zangmax=maxval(self%bin)
     zbinsiz=(zangmax-zangmin)/self%n%nbin

  case('toroidal')
     ! new sort key based on toroidal angle
     allocate(self%bin(iindict), self%oldict(iindict), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     do i=1,iindict
        zx=self%centr(1,i)
        zy=self%centr(2,i)
        self%bin(i)=atan2(zy,zx)
     end do
     self%oldict=self%dict
     ! bin size
     zangmin=minval(self%bin)
     zangmax=maxval(self%bin)
     zbinsiz=(zangmax-zangmin)/self%n%nbin

  end select new_key

  ! create dictionary
  call indict_bycluster(self%rdict,self%bin,iindict,zbinsiz,self%nrdict)

  ! new key
  allocate(self%newkey(iindict), stat=status)
  call log_alloc_check(m_name,s_name,3,status)
  do i=1,iindict
     self%newkey(i)=indict(self%rdict,self%bin(i),self%nrdict)
  end do

  !write(*,*) self%newkey
  if (self%n%rekey) then
     ! rekey
     do j=1,self%geobjl%ng
        ikey=self%key(j)
        indx=indict(self%dict,ikey,self%ndict)
        self%key(j)=self%newkey(indx)
     end do
     call indict_fromlist(self%dict,self%key,self%ndict)
     call indict_sort(self%dict,iwork,self%ndict,1)
     deallocate(self%centr,self%geobjl%obj)
     call geobjlist_centroids(self%geobjl,self%key,self%dict,self%ndict,self%centr)
     self%n%namekey='Angle '
     call log_error(m_name,s_name,10,log_info,'Data rekeyed on Angle')
  end if

end subroutine smanal_rekey
!---------------------------------------------------------------------
!> coordinate statistics subroutine
subroutine smanal_stats(self)

  !! arguments
  type(smanal_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='smanal_stats' !< subroutine name
  character(len=6) :: icstat !< statistics name
  integer(ki4) :: ikey !< local variable
  integer(ki4) :: indx !< local variable
  integer(ki4) :: indim !< local variable
  logical :: ilarea !< control area calculation

  ! initialise statistics
  allocate(self%nstat(self%n%totstat), &
 &self%stat(self%ndict,self%n%totstat),stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  allocate(self%areab(self%ndict), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  self%areab=0
  ilarea=.FALSE.
  do i=1,self%n%totstat
     self%nstat(i)=self%ndict
     icstat=self%n%namstat(i)
     stat_init_type: select case (icstat)
     case('intg  ', 'intgsq')
        self%stat(:,i)=0
        ilarea=.TRUE.
     case('maxabs')
        self%stat(:,i)=0
     case('max   ')
        self%stat(:,i)=-const_infty
     case('min   ', 'minabs')
        self%stat(:,i)=const_infty
     end select stat_init_type
  end do
  self%larea=ilarea

  if (ilarea) call geobjlist_area(self%geobjl,self%area)

  ! control calculation
  do i=1,self%n%totstat
     icstat=self%n%namstat(i)
     stat_type: select case (icstat)

     case ('intg  ')
        do k=1,self%geobjl%ng
           ikey=self%key(k)
           indx=indict(self%dict,ikey,self%ndict)
           self%stat(indx,i)=self%stat(indx,i)+self%area(k)*self%scal(k)
           self%areab(indx)=self%areab(indx)+self%area(k)
        end do

     case ('intgsq')
        do k=1,self%geobjl%ng
           ikey=self%key(k)
           indx=indict(self%dict,ikey,self%ndict)
           self%stat(indx,i)=self%stat(indx,i)+self%area(k)*self%scal(k)**2
        end do

     case ('max   ')
        do k=1,self%geobjl%ng
           ikey=self%key(k)
           indx=indict(self%dict,ikey,self%ndict)
           self%stat(indx,i)=max(self%stat(indx,i),self%scal(k))
        end do

     case ('min   ')
        do k=1,self%geobjl%ng
           ikey=self%key(k)
           indx=indict(self%dict,ikey,self%ndict)
           self%stat(indx,i)=min(self%stat(indx,i),abs(self%scal(k)))
        end do

     case ('maxabs')
        do k=1,self%geobjl%ng
           ikey=self%key(k)
           indx=indict(self%dict,ikey,self%ndict)
           self%stat(indx,i)=max(self%stat(indx,i),self%scal(k))
        end do

     case ('minabs')
        do k=1,self%geobjl%ng
           ikey=self%key(k)
           indx=indict(self%dict,ikey,self%ndict)
           self%stat(indx,i)=min(self%stat(indx,i),abs(self%scal(k)))
        end do

     end select stat_type
  end do

end subroutine smanal_stats
!---------------------------------------------------------------------
!> open new file, making up name
subroutine smanal_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer, intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='smanal_initwrite' !< subroutine name
  character(len=80) :: outputfile !< output file name
  !! logical :: unitused !< flag to test unit is available
  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then
  !! if (present(channel)) channel=i exit end if end do noutsm=i

  !! open file
  outputfile=trim(fileroot)//"_smanal.out"
  call log_value("Control data file",trim(outputfile))
  call misc_getfileunit(noutsm)
  open(unit=noutsm,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutsm,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if
  if (present(channel)) channel=noutsm

end subroutine smanal_initwrite
!---------------------------------------------------------------------
!> write smanal data
subroutine smanal_write(self,channel)

  !! arguments
  type(smanal_t), intent(in) :: self   !< smanal data structure
  integer, intent(in), optional :: channel   !< output channel for smanal data structure

  !! local
  character(*), parameter :: s_name='smanal_write' !< subroutine name
  integer(ki4) :: inewkey   !< value of new key
  integer :: iout   !< output channel for smanal data structure

  !! sort out unit
  if(present(channel)) then
     iout=channel
  else
     iout=noutsm
  end if

  if (allocated(self%rdict)) then
     call indict_write(self%rdict,self%nrdict,iout)
  end if
  if (allocated(self%dict)) then
     call indict_write(self%dict,self%ndict,iout)
  end if
  if (allocated(self%bin)) then
     write(iout,'(3(A18,1X))',iostat=status) 'new_key_values', self%n%namekey
     call log_write_check(m_name,s_name,10,status)
     do i=1,self%inbin
        inewkey=indict(self%rdict,self%bin(i),self%nrdict)
        write(iout,*,iostat=status) inewkey,self%bin(i),self%oldict(i)
        call log_write_check(m_name,s_name,11,status)
     end do
  end if
  ! skip rest
  return
  write(iout,*,iostat=status) 'smanal_formula'
  call log_write_check(m_name,s_name,18,status)
  write(iout,*,iostat=status) self%formula
  call log_write_check(m_name,s_name,19,status)
  write(iout,*,iostat=status) 'f'
  call log_write_check(m_name,s_name,20,status)
  write(iout,*,iostat=status) self%f
  call log_write_check(m_name,s_name,21,status)
  write(iout,*,iostat=status) 'nrpams'
  call log_write_check(m_name,s_name,46,status)
  write(iout,*,iostat=status) self%nrpams
  call log_write_check(m_name,s_name,47,status)
  if (self%nrpams>0) then
     write(iout,*,iostat=status) 'real_parameters'
     call log_write_check(m_name,s_name,48,status)
     write(iout,*,iostat=status) self%rpar
     call log_write_check(m_name,s_name,49,status)
  end if
  write(iout,*,iostat=status) 'nipams'
  call log_write_check(m_name,s_name,50,status)
  write(iout,*,iostat=status) self%nipams
  call log_write_check(m_name,s_name,51,status)
  if (self%nipams>0) then
     write(iout,*,iostat=status) 'integer_parameters'
     call log_write_check(m_name,s_name,52,status)
     write(iout,*,iostat=status) self%npar
     call log_write_check(m_name,s_name,53,status)
  end if

end subroutine smanal_write
!---------------------------------------------------------------------
!> write object data as gnuplot
subroutine smanal_writeg(self,select,channel)

  !! arguments
  type(smanal_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer, intent(in), optional :: channel   !< output channel for smanal data structure

  !! local
  character(*), parameter :: s_name='smanal_writeg' !< subroutine name
  integer(ki4) :: inacros   !< number of outputs for smanal data structure
  integer(ki4) :: idummy=0   !< dummy output
  character(len=80) :: icfmt !< output header format
  character(len=80) :: ircfmt !< output header format
  character(len=80) :: icrcfmt !< output data format
  real(kr8), dimension(:), allocatable :: zstat !< 1D work array
  real(kr8), dimension(:), allocatable :: zextr !< 1D work array
  real(kr8), dimension(3) :: zcentr !< average of centroids
  real(kr8) :: zareasum !< total area
  real(kr8) :: zareamax !< maximum area
  integer(ki4) :: inewkey   !< value of new key

  inacros=self%n%totstat+4
  write(icfmt,'("(A6,",i2,"(3X,A6,4X))")') inacros
  write(ircfmt,'("(I6,",i2,"(1PG12.5,1X))")') inacros
  write(icrcfmt,'("(A1,I5,",i2,"(1PG12.5,1X))")') inacros

  plot_type: select case(select)
  case('gnuanal')
     write(channel,icfmt,iostat=status) '#'//self%n%namekey, &
 &   (self%n%namstat(j), j=1,self%n%totstat),'x-cen','y-cen','z-cen','area'
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Error writing header')
     end if

     imaxstat=maxval(self%nstat)
     do i=1,imaxstat
        write(channel,ircfmt,iostat=status) self%dict(i), &
 &      (self%stat(i,j), j=1,self%n%totstat),(self%centr(j,i),j=1,3), &
 &      self%areab(i)
     end do
     if(status/=0) then
        call log_error(m_name,s_name,2,error_fatal,'Error writing data')
     end if

     allocate(zstat(self%n%totstat), zextr(self%n%totstat),stat=status)
     call log_alloc_check(m_name,s_name,3,status)
     zstat=0
     zextr=-const_infty
     zcentr=0
     zareasum=0
     zareamax=-const_infty
     do i=1,imaxstat
        zstat(:)=zstat(:)+self%stat(i,:)
        zextr(:)=max(zextr(:),self%stat(i,:))
        zcentr(:)=zcentr(:)+self%centr(:,i)
        zareasum=zareasum+self%areab(i)
        zareamax=max(zareamax,self%areab(i))
     end do
     zcentr(:)=zcentr(:)/imaxstat
     write(channel,'(A)',iostat=status) '# Totals / Centroid'
     write(channel,icrcfmt,iostat=status) '#',idummy, &
 &   (zstat(j), j=1,self%n%totstat),(zcentr(j),j=1,3), &
 &   zareasum
     if(status/=0) then
        call log_error(m_name,s_name,4,error_fatal,'Error writing data')
     end if
     write(channel,'(A)',iostat=status) '# Maxima / Centroid'
     write(channel,icrcfmt,iostat=status) '#',idummy, &
 &   (zextr(j), j=1,self%n%totstat),(zcentr(j),j=1,3), &
 &   zareamax
     if(status/=0) then
        call log_error(m_name,s_name,4,error_fatal,'Error writing data')
     end if

     deallocate(zstat,zextr)

  case('gnureanal')
     write(channel,'(A18,A18)',iostat=status) '#New_key:', self%n%newkey
     write(channel,'(A)',iostat=status) '#New_key Key_value     Old_key'
     if(status/=0) then
        call log_error(m_name,s_name,10,error_fatal,'Error writing header')
     end if
     do i=1,self%inbin
        inewkey=indict(self%rdict,self%bin(i),self%nrdict)
        write(channel,'(I8,1X,1PG12.5,1X,I8)',iostat=status) inewkey,self%bin(i),self%oldict(i)
        call log_write_check(m_name,s_name,11,status)
     end do

  case default

  end select plot_type

  call log_error(m_name,s_name,1,log_info,'gnuplot file produced')

end subroutine smanal_writeg
!---------------------------------------------------------------------
!> write object data as vtk
subroutine smanal_writev(self,select,channel)

  !! arguments
  type(smanal_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer, intent(in) :: channel   !< output channel for smanal data structure

  !! local
  character(*), parameter :: s_name='smanal_writev' !< subroutine name
  integer :: iout   !< output channel for smanal data structure
  integer(ki4) :: imaxstat   !< number of outputs for smanal data structure
  integer(ki4) :: ikey !< local variable
  integer(ki4) :: indx !< local variable

  iout=channel

  plot_type: select case(select)
  case('full')
     call geobjlist_writev(self%geobjl,'geometry',iout)
     iopt=1
     ! construct full array from statistic
     allocate(work(self%geobjl%ng), stat=status)
     call log_alloc_check(m_name,s_name,2,status)
     do k=1,self%n%totstat
        do l=1,self%geobjl%ng
           ikey=self%key(l)
           indx=indict(self%dict,ikey,self%ndict)
           work(l)=self%stat(indx,k)
        end do
        call vfile_dscalarwrite(work,self%geobjl%ng, &
 &      self%n%namstat(k),'CELL',iout,iopt)
        iopt=0
     end do
     call vfile_close

  case('small')
     ! construct small geobjlist
     imaxstat=maxval(self%nstat)
     allocate(igeobjl%posl%pos(imaxstat), igeobjl%obj2(imaxstat), &
 &   igeobjl%nodl(imaxstat), stat=status)
     call log_alloc_check(m_name,s_name,5,status)

     do j=1,imaxstat
        do k=1,3
           igeobjl%posl%pos(j)%posvec(k)=self%centr(k,j)
        end do
        igeobjl%obj2(j)%typ=VTK_VERTEX
        igeobjl%nodl(j)=j
     end do
     igeobjl%ng=imaxstat
     igeobjl%np=imaxstat
     ! and write geometry
     call geobjlist_writev(igeobjl,'geometry',iout)
     iopt=1
     do k=1,self%n%totstat
        call vfile_dscalarwrite(self%stat(:,k),self%nstat(k), &
 &      self%n%namstat(k),'CELL',iout,iopt)
        iopt=0
     end do
     call vfile_close

  case default

  end select plot_type

  call log_error(m_name,s_name,10,log_info,'vtk file produced')

end subroutine smanal_writev
!---------------------------------------------------------------------
!> close write file
subroutine smanal_closewrite

  !! local
  character(*), parameter :: s_name='smanal_closewrite' !< subroutine name

  !! close file
  close(unit=noutsm,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine smanal_closewrite
!---------------------------------------------------------------------
!> delete object
subroutine smanal_delete(self)

  !! arguments
  type(smanal_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='smanal_delete' !< subroutine name

  mode_deallocate: select case (self%n%mode)
  case('regular')
     call geobjlist_delete(self%geobjl)
     deallocate(self%key)
     deallocate(self%scal)
     deallocate(self%dict)
     deallocate(self%centr)
     deallocate(self%areab)
     deallocate(self%stat,self%nstat)
     if (self%larea) deallocate(self%area)
     if (allocated(self%rdict)) deallocate(self%rdict)
     if (allocated(self%bin)) deallocate(self%bin)
  case default
  end select mode_deallocate

end subroutine smanal_delete
!---------------------------------------------------------------------
!> close file
subroutine smanal_close

  !! local
  character(*), parameter :: s_name='smanal_close' !< subroutine name

  !! close file
  close(unit=ninsm,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine smanal_close

end module smanal_m
