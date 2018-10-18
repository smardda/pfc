module vfile_m

  use const_kind_m
  use log_m
  use misc_m

  implicit none
  private

! public subroutines
  public :: &
  vfile_init, & !< initialise vtk file
  vfile_rscalarread, & !< read vtk real scalars
  vfile_rvectorread, & !< read vtk real vectors
  vfile_dscalarread, & !< read vtk real kr8 scalars
  vfile_iscalarread, & !< read vtk integer scalars
  vfile_rfieldread, & !< read vtk real scalar from vtk file with field consisting only of scalar arrays
  vfile_rscalarwrite, & !< write vtk real scalars
  vfile_dscalarwrite, & !< write vtk real kr8 scalars
  vfile_dvectorwrite, & !< write vtk real kr8 vectors
  vfile_iscalarwrite, & !< write vtk integer scalars
  vfile_close, & !< close vtk file
  vfile_getfmt,&   !< find vtk file numerical format
  vfile_skiptolabel !< find line label in vtk file

! public types

  integer(ki4), parameter, public :: vfile_made_up_data = 1 !< local variable

! private variables
  character(*), parameter :: m_name='vfile_m' !< module name
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  character(len=80) :: icfile !< file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer, save :: nin  !< file unit for input
  integer, save :: nplot  !< file unit for output
  integer :: status  !< status flag
  character(len=7) :: cwrite !< write status of file
  logical :: filefound !< true if file exists
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> initialise vtk file
subroutine vfile_init(fplot,descriptor,kplot)

  !! arguments
  character(len=*), intent(in) :: fplot !< file name root
  character(len=*), intent(in) :: descriptor !< dataset descriptor
  integer, intent(inout) :: kplot   !< unit number

  !! local
  character(*), parameter :: s_name='vfile_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  !! open file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kplot=i exit end if end do

  ! stop under user control if write-protected file already exists
  icfile=trim(fplot)//'.vtk'
  inquire(file=trim(icfile),exist=filefound,write=cwrite)
  !DBG write(*,*) 'cwrite=',cwrite !DBG
  if (filefound) then
     call log_value("vtk data file ",icfile)
     call log_error(m_name,s_name,1,error_warning,'File already exists and may be write protected')
  end if
  call misc_getfileunit(kplot)
  open(unit=kplot,file=icfile,iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,2,error_fatal,'Error opening vtk data file')
  else
     call log_error(m_name,s_name,2,log_info,'vtk data file opened')
  end if

  !! write vtk header
  write(kplot,'(''# vtk DataFile Version 2.0'')')
  write(kplot,'(A)') descriptor
  write(kplot,'(''ASCII'')')
  write(kplot,'(''           '')')

  nplot=kplot
  !DBGwrite(*,*) 'nplot=',nplot !DBG

end subroutine vfile_init
!---------------------------------------------------------------------
!> read vtk real scalars
subroutine vfile_rscalarread(self,kp,infile,kcname,kin,kopt)
  !! arguments
  real(kr4), dimension(:), allocatable, intent(inout) :: self !< real scalar list data
  integer(ki4), intent(inout) :: kp   !< size of scalar list data
  character(*),intent(in) :: infile !< name of input file
  character(*),intent(in) :: kcname !< name of field required
  integer, intent(inout) :: kin   !< input channel for scalar list data
  !> if positive on input, assume unit open and do not terminate execution if
  !! data is missing, return positive code if trouble, else zero
  !! if zero on input, missing data is a fatal error
  integer(ki4), intent(inout) :: kopt   !< .

  !! local
  character(*), parameter :: s_name='vfile_rscalarread' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer :: ierror   !<  whether error is to be fatal or not
  integer(ki4) :: insca   !< number of scalars
  character(80) :: sname !< name of scalar
  integer(ki4) :: islen   !< length of scalar field name
  integer(ki4) :: islen2   !< length of required scalar field name

  logical :: isnumb !< local variable
  external isnumb

  ibuf1=adjustl(kcname)
  islen2=max(2,scan(ibuf1,' '))-1
  if(kopt>0) then
     ierror=error_warning
     !! assume unit already open and reading infile
     if (kin==0) then
        inquire(file=infile,number=kin,iostat=status)
        if(status/=0.OR.kin==-1)then
           !! error opening file
           call log_error(m_name,s_name,1,error_fatal,'Error opening scalar list data file')
        else
           call log_error(m_name,s_name,1,log_info,'Scalar list data file opened')
        end if
     end if
     nin=kin
  else
     ierror=error_fatal
     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kin=i exit end if end do nin=kin

     !! open file
     call misc_getfileunit(nin)
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        call log_error(m_name,s_name,2,error_fatal,'Error opening scalar list data file')
     else
        call log_error(m_name,s_name,2,log_info,'Scalar list data file opened')
     end if
     kin=nin

  end if

  !! File unit now sorted, get to where point data begin
  !! read local header information
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
           iltest=isnumb(ibuf2,insca,11)
           kp=insca
           exit
        else if(ibuf2(1:8)=='POLYGONS') then
           iltest=isnumb(ibuf2,insca,9)
           kp=insca
           exit
        else if(ibuf2(1:9)=='CELL_DATA') then
           iltest=isnumb(ibuf2,insca,10)
           kp=insca
           exit
        end if
     end if
  end do

  !! find scalar header
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,4,ierror,'Error reading header data')
        kopt=4
        return
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:7)=='SCALARS') then
           ibuf1=adjustl(ibuf2(8:))
           islen=max(2,scan(ibuf1,' '))-1
           sname=ibuf1(:islen)
           if (sname(:islen)==kcname(:islen2)) then
              call log_value("Scalar field found ",sname)
              exit
           else
              call log_value("Skipped scalar field ",sname)
           end if
        end if
     end if
  end do

  !! allocate scalar storage
  if(kp>0) then
     allocate(self(kp),stat=status)
     call log_alloc_check(m_name,s_name,5,status)
  else
     call log_error(m_name,s_name,6,ierror,'No scalar data')
     kopt=6
     return
  end if
  !! skip LOOKUP table default
  read(nin,fmt='(a)',iostat=status) ibuf1
  !! read coordinates
  read(nin,*,iostat=status) (self(j),j=1,kp)
  call log_read_check(m_name,s_name,7,status)
  print '("number of scalars read = ",i10)',kp
  call log_value("number of scalars read ",kp)
  kopt=0

end subroutine vfile_rscalarread
!---------------------------------------------------------------------
!> read vtk real vectors
subroutine vfile_rvectorread(self,kp,kadim,infile,kcname,kin,kopt)
  !! arguments
  real(kr4), dimension(:), allocatable, intent(inout) :: self !< real vector list data
  integer(ki4), intent(inout) :: kp   !< size of vector list data
  integer(ki4), dimension(3), intent(out) :: kadim   !< dimensions of vector list data
  character(*),intent(in) :: infile !< name of input file
  character(*),intent(in) :: kcname !< name of field required
  integer, intent(inout) :: kin   !< input channel for vector list data
  !> if positive on input, assume unit open and do not terminate execution if
  !! data is missing, return positive code if trouble, else zero
  !! if zero on input, missing data is a fatal error
  integer(ki4), intent(inout) :: kopt   !< .

  !! local
  character(*), parameter :: s_name='vfile_rvectorread' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer :: ierror   !<  whether error is to be fatal or not
  integer(ki4) :: insca   !< number of vectors
  character(80) :: sname !< name of vector
  integer(ki4) :: islen   !< length of vector field name
  integer(ki4) :: islen2   !< length of required vector field name

  logical :: isnumb !< local variable
  external isnumb

  kadim=0
  ibuf1=adjustl(kcname)
  islen2=max(2,scan(ibuf1,' '))-1
  if(kopt>0) then
     ierror=error_warning
     !! assume unit already open and reading infile
     if (kin==0) then
        inquire(file=infile,number=kin,iostat=status)
        if(status/=0.OR.kin==-1)then
           !! error opening file
           call log_error(m_name,s_name,1,error_fatal,'Error opening vector list data file')
        else
           call log_error(m_name,s_name,1,log_info,'Vector list data file opened')
        end if
     end if
     nin=kin
  else
     ierror=error_fatal
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
  !! read local header information
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
        if(ibuf2(1:10)=='DIMENSIONS') then
           read(ibuf2(11:),*,iostat=status) kadim
           call log_read_check(m_name,s_name,4,status)
        else if(ibuf2(1:10)=='POINT_DATA') then
           iltest=isnumb(ibuf2,insca,11)
           kp=insca
           exit
        else if(ibuf2(1:8)=='POLYGONS') then
           iltest=isnumb(ibuf2,insca,9)
           kp=insca
           exit
        else if(ibuf2(1:9)=='CELL_DATA') then
           iltest=isnumb(ibuf2,insca,10)
           kp=insca
           exit
        end if
     end if
  end do

  if (kadim(1)/=0) then
     if (kp-kadim(1)*kadim(2)*kadim(3)/=0) then
        call log_error(m_name,s_name,3,error_warning,'Mismatch in vector array dimensions')
     end if
  end if

  !! find vector header
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,4,ierror,'Error reading header data')
        kopt=4
        return
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:7)=='VECTORS') then
           ibuf1=adjustl(ibuf2(8:))
           islen=max(2,scan(ibuf1,' '))-1
           sname=ibuf1(:islen)
           if (sname(:islen)==kcname(:islen2)) then
              call log_value("Vector field found ",sname)
              exit
           else
              call log_value("Skipped vector field ",sname)
           end if
        end if
     end if
  end do

  !! allocate vector storage
  if(kp>0) then
     allocate(self(3*kp),stat=status)
     call log_alloc_check(m_name,s_name,5,status)
  else
     call log_error(m_name,s_name,6,ierror,'No vector data')
     kopt=6
     return
  end if
  !! read coordinates
  read(nin,*,iostat=status) (self(j),j=1,3*kp)
  call log_read_check(m_name,s_name,7,status)
  print '("number of vectors read = ",i10)',kp
  call log_value("number of vectors read ",kp)
  kopt=0

end subroutine vfile_rvectorread
!---------------------------------------------------------------------
!> read vtk real scalars
subroutine vfile_dscalarread(self,kp,infile,kcname,kin,kopt)
  !! arguments
  real(kr8), dimension(:), allocatable, intent(inout) :: self !< real scalar list data
  integer(ki4), intent(inout) :: kp   !< size of scalar list data
  character(*),intent(in) :: infile !< name of input file
  character(*),intent(in) :: kcname !< name of field required
  integer, intent(inout) :: kin   !< input channel for scalar list data
  !> if positive on input, assume unit open and do not terminate execution if
  !! data is missing, return positive code if trouble, else zero
  !! if zero on input, missing data is a fatal error
  integer(ki4), intent(inout) :: kopt   !< .

  !! local
  character(*), parameter :: s_name='vfile_dscalarread' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer :: ierror   !<  whether error is to be fatal or not
  integer(ki4) :: insca   !< number of scalars
  character(80) :: sname !< name of scalar
  integer(ki4) :: islen   !< length of scalar field name
  integer(ki4) :: islen2   !< length of required scalar field name

  logical :: isnumb !< local variable
  external isnumb

  ibuf1=adjustl(kcname)
  islen2=max(2,scan(ibuf1,' '))-1
  if(kopt>0) then
     ierror=error_warning
     !! assume unit already open and reading infile
     if (kin==0) then
        inquire(file=infile,number=kin,iostat=status)
        if(status/=0.OR.kin==-1)then
           !! error opening file
           call log_error(m_name,s_name,1,error_fatal,'Error opening scalar list data file')
        else
           call log_error(m_name,s_name,1,log_info,'Scalar list data file opened')
        end if
     end if
     nin=kin
  else
     ierror=error_fatal
     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kin=i exit end if end do nin=kin

     !! open file
     call misc_getfileunit(nin)
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        call log_error(m_name,s_name,2,error_fatal,'Error opening scalar list data file')
     else
        call log_error(m_name,s_name,2,log_info,'Scalar list data file opened')
     end if
     kin=nin

  end if

  !! File unit now sorted, get to where point data begin
  !! read local header information
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
           iltest=isnumb(ibuf2,insca,11)
           kp=insca
           exit
        else if(ibuf2(1:8)=='POLYGONS') then
           iltest=isnumb(ibuf2,insca,9)
           kp=insca
           exit
        else if(ibuf2(1:9)=='CELL_DATA') then
           iltest=isnumb(ibuf2,insca,10)
           kp=insca
           exit
        end if
     end if
  end do

  !! find scalar header
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,4,ierror,'Error reading header data')
        kopt=4
        return
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:7)=='SCALARS') then
           ibuf1=adjustl(ibuf2(8:))
           islen=max(2,scan(ibuf1,' '))-1
           sname=ibuf1(:islen)
           if (sname(:islen)==kcname(:islen2)) then
              call log_value("Scalar field found ",sname)
              exit
           else
              call log_value("Skipped scalar field ",sname)
           end if
        end if
     end if
  end do

  !! allocate scalar storage
  if(kp>0) then
     allocate(self(kp),stat=status)
     call log_alloc_check(m_name,s_name,5,status)
  else
     call log_error(m_name,s_name,6,ierror,'No scalar data')
     kopt=6
     return
  end if
  !! skip LOOKUP table default
  read(nin,fmt='(a)',iostat=status) ibuf1
  !! read coordinates
  read(nin,*,iostat=status) (self(j),j=1,kp)
  call log_read_check(m_name,s_name,7,status)
  print '("number of scalars read = ",i10)',kp
  call log_value("number of scalars read ",kp)
  kopt=0

end subroutine vfile_dscalarread
!---------------------------------------------------------------------
!> read vtk integer scalars
subroutine vfile_iscalarread(kself,kp,infile,kcname,kin,kopt,kcdata)
  !! arguments
  integer(ki4), dimension(:), allocatable, intent(inout) :: kself !< integer scalar list data
  integer(ki4), intent(inout) :: kp   !< size of scalar list data
  character(*),intent(in) :: infile !< name of input file
  character(*),intent(in) :: kcname !< name of field required
  integer, intent(inout) :: kin   !< input channel for scalar list data
  !> if positive on input, assume unit open and do not terminate execution if
  !! data is missing, return positive code if trouble, else zero
  !! if zero on input, missing data is a fatal error
  integer(ki4), intent(inout) :: kopt   !< .
  character(*),intent(in), optional :: kcdata !< nonblank, must match POINT_ or CELL_ DATA

  !! local
  character(*), parameter :: s_name='vfile_iscalarread' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer :: ierror   !<  whether error is to be fatal or not
  integer(ki4) :: insca   !< number of scalars
  character(80) :: sname !< name of scalar
  integer(ki4) :: islen   !< length of scalar field name
  integer(ki4) :: islen2   !< length of required scalar field name
  logical :: ilfound !< flag whether field name found
  logical :: ilpointd !< point data is acceptable
  logical :: ilcelld !< cell data is acceptable

  logical :: isnumb !< local variable
  external isnumb

  ibuf1=adjustl(kcname)
  islen2=max(2,scan(ibuf1,' '))-1
  if(kopt>0) then
     ierror=error_warning
     !! assume unit already open and reading infile
     if (kin==0) then
        inquire(file=infile,number=kin,iostat=status)
        if(status/=0.OR.kin==-1)then
           !! error opening file
           call log_error(m_name,s_name,1,error_fatal,'Error opening scalar list data file')
        else
           call log_error(m_name,s_name,1,log_info,'Scalar list data file opened')
        end if
     end if
     nin=kin
  else
     ierror=error_fatal
     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kin=i exit end if end do nin=kin

     !! open file
     call misc_getfileunit(nin)
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        call log_error(m_name,s_name,2,error_fatal,'Error opening scalar list data file')
     else
        call log_error(m_name,s_name,2,log_info,'Scalar list data file opened')
     end if
     kin=nin

  end if

  ilpointd=.TRUE.
  ilcelld=.TRUE.
  if (present(kcdata)) then
     ! to avoid optional arg, replace with if (kcdata(1)==' ')
     ilpointd=(kcdata(1:5)=='POINT')
     ilcelld=(kcdata(1:4)=='CELL')
     if ( .NOT.(ilpointd.OR.ilcelld) ) then
        call log_error(m_name,s_name,3,error_warning,'Neither point nor cell data requested')
     end if
  end if

  !! File unit now sorted, get to where point/cell data begin
  !! read local header information
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
           if (ilpointd) then
              iltest=isnumb(ibuf2,insca,11)
              kp=insca
              exit
           end if
        else if(ibuf2(1:8)=='POLYGONS') then
           iltest=isnumb(ibuf2,insca,9)
           kp=insca
           exit
        else if(ibuf2(1:9)=='CELL_DATA') then
           if (ilcelld) then
              iltest=isnumb(ibuf2,insca,10)
              kp=insca
              exit
           end if
        end if
     end if
  end do

  !! find scalar header
  ilfound=.FALSE.
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,4,ierror,'Error reading header data')
        kopt=4
        return
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:7)=='SCALARS') then
           ibuf1=adjustl(ibuf2(8:))
           islen=max(2,scan(ibuf1,' '))-1
           sname=ibuf1(:islen)
           if (sname(:islen)==kcname(:islen2)) then
              call log_value("Scalar field found ",sname)
              ilfound=.TRUE.
              exit
           else
              call log_value("Skipped scalar field ",sname)
           end if
        end if
     end if
  end do

  !! allocate scalar storage
  if(ilfound) then
     if(kp>0) then
        allocate(kself(kp),stat=status)
        call log_alloc_check(m_name,s_name,5,status)
     else
        call log_error(m_name,s_name,6,ierror,'No scalar data')
        kopt=6
        return
     end if
  else if(vfile_made_up_data/=0) then
     if(kp>0) then
        allocate(kself(kp),stat=status)
        call log_alloc_check(m_name,s_name,7,status)
     else
        call log_error(m_name,s_name,8,ierror,'No scalar data')
        kopt=8
        return
     end if
  else
     call log_error(m_name,s_name,9,ierror,'No scalar data')
     kopt=9
     return
  end if
  !! data
  if(ilfound) then
     !! skip LOOKUP table default
     read(nin,fmt='(a)',iostat=status) ibuf1
     !! read data
     read(nin,*,iostat=status) (kself(j),j=1,kp)
     call log_read_check(m_name,s_name,10,status)
     print '("number of scalars read = ",i10)',kp
     call log_value("number of scalars read ",kp)
  else if(vfile_made_up_data/=0) then
     kself=vfile_made_up_data
     call log_error(m_name,s_name,10,error_warning,'Made up data')
     call log_value("number of scalars made up ",kp)
  end if
  kopt=0

end subroutine vfile_iscalarread
!---------------------------------------------------------------------
!> read vtk real scalar from vtk file with field consisting only of scalar arrays
subroutine vfile_rfieldread(self,kp,infile,kcname,kin,kopt)
  !! arguments
  real(kr4), dimension(:), allocatable, intent(inout) :: self !< real scalar data
  integer(ki4), intent(inout) :: kp   !< size of scalar data
  character(*),intent(in) :: infile !< name of input file
  character(*),intent(in) :: kcname !< name of scalar array required
  integer, intent(inout) :: kin   !< input channel for field data
  !> if unity on input, assume unit open and do not terminate execution if
  !! data is missing, return positive code if trouble, else zero
  !! kopt=2, assume at start point of field
  !! if zero on input, missing data is a fatal error
  integer(ki4), intent(inout) :: kopt   !< .

  !! local
  character(*), parameter :: s_name='vfile_rfieldread' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer :: ierror   !<  whether error is to be fatal or not
  integer(ki4) :: insca   !< number of scalars in array
  integer(ki4) :: infld   !< number of scalar arrays
  character(80) :: sname !< name of scalar
  integer(ki4) :: inlastbl   !< position of last blank in buffer
  integer(ki4) :: islen   !< length of scalar array name
  integer(ki4) :: islen2   !< length of required scalar array name
  integer(ki4) :: ilen   !< length of line string

  logical :: ilok !< flag whether named array found
  logical :: isnumb !< local variable
  external isnumb

  ibuf1=adjustl(kcname)
  islen2=max(2,scan(ibuf1,' '))-1
  if(kopt>0) then
     ierror=error_warning
     !! assume unit already open and reading infile
     if (kin==0) then
        inquire(file=infile,number=kin,iostat=status)
        if(status/=0.OR.kin==-1)then
           !! error opening file
           call log_error(m_name,s_name,1,error_fatal,'Error opening field data file')
        else
           call log_error(m_name,s_name,1,log_info,'Field data file opened')
        end if
     end if
     nin=kin
  else
     ierror=error_fatal
     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kin=i exit end if end do nin=kin

     !! open file
     call misc_getfileunit(nin)
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        call log_error(m_name,s_name,2,error_fatal,'Error opening field data file')
     else
        call log_error(m_name,s_name,2,log_info,'Field data file opened')
     end if
     kin=nin

  end if

  if (kopt<=1) then
     !! File unit now sorted, get to where point data begin
     !! read local header information
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
              iltest=isnumb(ibuf2,insca,11)
              kp=insca
              exit
           else if(ibuf2(1:8)=='POLYGONS') then
              iltest=isnumb(ibuf2,insca,9)
              kp=insca
              exit
           else if(ibuf2(1:9)=='CELL_DATA') then
              iltest=isnumb(ibuf2,insca,10)
              kp=insca
              exit
           end if
        end if
     end do
  end if


  !! find field header
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     ! write(*,*) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,4,ierror,'Error reading header data')
        kopt=4
        return
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:5)=='FIELD') then
           ilen=len_trim(ibuf2)
           inlastbl=max( 1,scan(ibuf2(1:ilen),' ',.TRUE.) )
           iltest=isnumb(ibuf1,infld,inlastbl)
           if (infld>0) then
              call log_value("Number of scalar arrays found ",infld)
              exit
           else
              call log_value("No arrays found in file infld ",infld)
              call log_error(m_name,s_name,5,ierror,'Error looking for scalar array')
              kopt=5
           end if
        end if
     end if
  end do

  !! allocate array storage
  if(kp>0) then
     allocate(self(kp),stat=status)
     call log_alloc_check(m_name,s_name,5,status)
  else
     call log_error(m_name,s_name,6,ierror,'No array data')
     kopt=6
     return
  end if
  !! loop to find right array
  !! (no LOOKUP table default)
  ilok=.FALSE.
  do l=1,infld
     read(nin,fmt='(a)',iostat=status) ibuf1
     read(nin,*,iostat=status) (self(j),j=1,kp)
     call log_read_check(m_name,s_name,7,status)
     islen=max(2,scan(ibuf1,' '))-1
     sname=ibuf1(:islen)
     if (sname(:islen)==kcname(:islen2)) then
        call log_value("Scalar array found ",sname)
        ilok=.TRUE.
        exit
     else
        call log_value("Skipped scalar array ",sname)
     end if
  end do

  if (ilok) then
     !! log data read in
     print '("number of scalars read = ",i10)',kp
     call log_value("number of scalars read ",kp)
     kopt=0
  else
     call log_value("Scalar array not found in field file, name ",kcname)
     call log_error(m_name,s_name,8,ierror,'Error looking for scalar array')
     kopt=8
  end if

end subroutine vfile_rfieldread
!---------------------------------------------------------------------
!> write vtk real scalars
subroutine vfile_rscalarwrite(self,kp,kcname,kctyp,kplot,kheader)
  !! arguments
  real(kr4), dimension(kp), intent(in) :: self !< real scalar list data
  integer(ki4), intent(in) :: kp   !< size of scalar list data
  character(*),intent(in) :: kcname !< name of field
  character(*),intent(in) :: kctyp !< type of data
  integer, intent(in) :: kplot   !< output channel for scalar list data
  integer(ki4), intent(in) :: kheader   !< header options

  !! local
  character(*), parameter :: s_name='vfile_rscalarwrite' !< subroutine name
  integer(ki4) :: islen   !< length of scalar field name
  integer(ki4) :: islen2   !< length of required scalar field name

  ibuf1=adjustl(kctyp)
  islen=len_trim(ibuf1)
  ibuf1=adjustl(kcname)
  islen2=len_trim(ibuf1)

  !! output scalar header if kheader is unity
  if(kheader==1) then
     write(kplot,'(A,''_DATA'',1X,I9)',iostat=status) kctyp(1:islen),kp
  end if
  call log_write_check(m_name,s_name,1,status)
  write(kplot,'(''SCALARS '',A,'' float 1'')',iostat=status) kcname(1:islen2)
  call log_write_check(m_name,s_name,2,status)
  write(kplot,'(''LOOKUP_TABLE default'')',iostat=status)
  call log_write_check(m_name,s_name,3,status)

  !! write scalars
  do j=1,kp
     write(kplot,fmt=cfmtbs,iostat=status) self(j)
     call log_write_check(m_name,s_name,5,status)
  end do

  call log_value("number of scalars written ",kp)

end subroutine vfile_rscalarwrite
!---------------------------------------------------------------------
!> write vtk real scalars
subroutine vfile_dscalarwrite(self,kp,kcname,kctyp,kplot,kheader)
  !! arguments
  real(kr8), dimension(kp), intent(in) :: self !< real scalar list data
  integer(ki4), intent(in) :: kp   !< size of scalar list data
  character(*),intent(in) :: kcname !< name of field
  character(*),intent(in) :: kctyp !< type of data
  integer, intent(in) :: kplot   !< output channel for scalar list data
  integer(ki4), intent(in) :: kheader   !< header options

  !! local
  character(*), parameter :: s_name='vfile_dscalarwrite' !< subroutine name
  integer(ki4) :: islen   !< length of scalar field name
  integer(ki4) :: islen2   !< length of required scalar field name

  ibuf1=adjustl(kctyp)
  islen=len_trim(ibuf1)
  ibuf1=adjustl(kcname)
  islen2=len_trim(ibuf1)

  !! output scalar header if kheader is unity
  if(kheader==1) then
     write(kplot,'(A,''_DATA'',I8)',iostat=status) kctyp(1:islen),kp
  end if
  call log_write_check(m_name,s_name,1,status)
  write(kplot,'(''SCALARS '',A,'' float 1'')',iostat=status) kcname(1:islen2)
  call log_write_check(m_name,s_name,2,status)
  write(kplot,'(''LOOKUP_TABLE default'')',iostat=status)
  call log_write_check(m_name,s_name,3,status)

  !! write scalars
  do j=1,kp
     write(kplot,fmt=cfmtbs,iostat=status) self(j)
     call log_write_check(m_name,s_name,5,status)
  end do

  call log_value("number of scalars written ",kp)

end subroutine vfile_dscalarwrite
!---------------------------------------------------------------------
!> write vtk double kr8 vectors
subroutine vfile_dvectorwrite(self,kp,kcname,kctyp,kplot,kheader)
  !! arguments
  real(kr8), dimension(3,kp), intent(in) :: self !< real vector list data
  integer(ki4), intent(in) :: kp   !< size of vector list data
  character(*),intent(in) :: kcname !< name of field
  character(*),intent(in) :: kctyp !< type of data
  integer, intent(in) :: kplot   !< output channel for vector list data
  integer(ki4), intent(in) :: kheader   !< header options

  !! local
  character(*), parameter :: s_name='vfile_dvectorwrite' !< subroutine name
  integer(ki4) :: islen   !< length of vector field name
  integer(ki4) :: islen2   !< length of required vector field name

  ibuf1=adjustl(kctyp)
  islen=len_trim(ibuf1)
  ibuf1=adjustl(kcname)
  islen2=len_trim(ibuf1)

  !! output vector header if kheader is unity
  if(kheader==1) then
     write(kplot,'(A,''_DATA'',I8)',iostat=status) kctyp(1:islen),kp
  end if
  call log_write_check(m_name,s_name,1,status)
  write(kplot,'(''VECTORS '',A,'' float'')',iostat=status) kcname(1:islen2)
  call log_write_check(m_name,s_name,2,status)

  !! write vectors
  do j=1,kp
     write(kplot,cfmtbv1,iostat=status) (self(k,j),k=1,3)
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Error writing double vector')
     end if
  end do

  call log_value("number of vectors written ",kp)

end subroutine vfile_dvectorwrite
!---------------------------------------------------------------------
!> write vtk real scalars
subroutine vfile_iscalarwrite(self,kp,kcname,kctyp,kplot,kheader)
  !! arguments
  integer(ki4), dimension(kp), intent(in) :: self !< integer scalar list data
  integer(ki4), intent(in) :: kp   !< size of scalar list data
  character(*),intent(in) :: kcname !< name of field
  character(*),intent(in) :: kctyp !< type of data
  integer, intent(in) :: kplot   !< output channel for scalar list data
  integer(ki4), intent(in) :: kheader   !< header options

  !! local
  character(*), parameter :: s_name='vfile_iscalarwrite' !< subroutine name
  integer(ki4) :: islen   !< length of scalar field name
  integer(ki4) :: islen2   !< length of required scalar field name

  ibuf1=adjustl(kctyp)
  islen=len_trim(ibuf1)
  ibuf1=adjustl(kcname)
  islen2=len_trim(ibuf1)

  !! output scalar header
  if(kheader==1) then
     write(kplot,'(A,''_DATA'',I8)',iostat=status) kctyp(1:islen),kp
  end if
  call log_write_check(m_name,s_name,1,status)
  write(kplot,'(''SCALARS '',A,'' int '')',iostat=status) kcname(1:islen2)
  call log_write_check(m_name,s_name,2,status)
  write(kplot,'(''LOOKUP_TABLE default'')',iostat=status)
  call log_write_check(m_name,s_name,3,status)

  !! write scalars
  do j=1,kp
     write(kplot,fmt='(I9)',iostat=status) self(j)
     call log_write_check(m_name,s_name,5,status)
  end do

  call log_value("number of scalars written ",kp)

end subroutine vfile_iscalarwrite
!---------------------------------------------------------------------
!> close vis vtk file on unit nplot
subroutine vfile_close

  close(nplot,iostat=status)
  !DBG write(*,*) 'iostat=',status !DBG

end subroutine vfile_close
!---------------------------------------------------------------------
!> find vtk file numerical format
subroutine vfile_getfmt(infile,kfmta,kin)
  !! arguments
  character(*),intent(in) :: infile !< name of input file
  integer(ki4), dimension(2), intent(out) :: kfmta   !< format as number of entries per line
  integer, intent(inout), optional :: kin   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='vfile_getfmt' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer(ki4) :: ir   !< counter for file header read
  integer(ki4) :: ifmt   !< formatting of position vectors
  character(len=132) :: bigbuf !<big buffer for input/output
  integer(ki4) :: inpt   !< number of points

  logical :: isnumb !< local variable
  external isnumb

  if(present(kin).AND.kin>0) then
     !! assume unit already open and reading infile
     nin=kin
     rewind(nin)
  else
     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

     !! open file
     call misc_getfileunit(nin)
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        write(*,*) 'infile=',infile
        call log_error(m_name,s_name,1,error_fatal,'Error opening vtk data file')
     else
        call log_error(m_name,s_name,2,log_info,'vtk data file opened')
     end if
     if (present(kin)) kin=nin
  end if

  !! first set of reads to determine format of position vectors
  ir=0
  do
     ir=ir+1
     read(nin,fmt='(a)',iostat=status) ibuf1
     call log_read_check(m_name,s_name,10,status)
     !! look for key at start of line
     ibuf2=adjustl(ibuf1)
     if(ibuf2(1:6)=='POINTS') exit
  end do
  if (ir>20) then
     call log_error(m_name,s_name,11,error_fatal,'Error reading header data')
  end if
  iltest=isnumb(ibuf2,inpt,7)
  read(nin,fmt='(a)',iostat=status) bigbuf
  call log_read_check(m_name,s_name,12,status)
  call misc_countnos(bigbuf,ifmt)
  if (ifmt<0.OR.ifmt>4) then
     call log_error(m_name,s_name,13,error_fatal,'Error cannot count points on input line')
  end if
  kfmta(1)=ifmt
  kfmta(2)=inpt

  rewind(nin)

end subroutine vfile_getfmt
!---------------------------------------------------------------------
!> find line label in vtk file
subroutine vfile_skiptolabel(infile,kclabel,kin)
  !! arguments
  character(*),intent(in) :: infile !< name of input file, ignored if kin set positive
  character(*),intent(in) :: kclabel !< find this label in file
  integer, intent(inout), optional :: kin   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='vfile_skiptolabel' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer(ki4) :: ir   !< counter for file header read
  integer(ki4) :: ilen   !< length of string in file

  if(present(kin).AND.kin>0) then
     !! assume unit already open and reading infile
     nin=kin
     rewind(nin)
  else
     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

     !! open file
     call misc_getfileunit(nin)
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        write(*,*) 'infile=',infile
        call log_error(m_name,s_name,1,error_fatal,'Error opening vtk data file')
     else
        call log_error(m_name,s_name,2,log_info,'vtk data file opened')
     end if
     if (present(kin)) kin=nin
  end if

  !! reads to find label
  ir=0
  do
     ir=ir+1
     read(nin,fmt='(a)',iostat=status) ibuf1
     call log_read_check(m_name,s_name,10,status)
     !! look for key at start of line
     ibuf2=adjustl(ibuf1)
     ilen=scan(ibuf2,' ')
     if(ibuf2(1:ilen-1)==kclabel) exit
  end do
  if (ir>20) then
     call log_error(m_name,s_name,11,error_fatal,'Error reading header data')
  end if
  call log_value("Successfully found label ",kclabel)

end subroutine vfile_skiptolabel

subroutine misc_countnos(bigbuf,kfmt)
  character(len=*),intent(in) :: bigbuf !< buffer for input
  integer(ki4), intent(out) :: kfmt !< format of buffer - number of 3-vectors
  character(len=132) :: ibuf !< buffer for input/output
  integer(ki4) :: ilen !< length of string
  integer(ki4) :: iblan !< number of bank substrings
  integer(ki4) :: isw !< switch on if last character was not blank
  integer(ki4) :: ji !< loop variable
  iblan=0
  ibuf=adjustl(bigbuf)
  ilen=len_trim(ibuf)
  isw=1
  do ji=1,ilen
     if (ibuf(ji:ji)==' ') then
        if (isw/=0) then
           iblan=iblan+1
           isw=0
        end if
     else
        isw=1
     end if
  end do
  kfmt=(iblan+1)/3
end subroutine misc_countnos

end module vfile_m
