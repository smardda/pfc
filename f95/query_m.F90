module query_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use misc_m
  use control_h
  use position_h
  use position_m

  implicit none
  private

! public subroutines
  public :: query_readv, & ! read in visualisation format
  query_writev, & ! write in visualisation format
  query_read, & ! read all query data
  query_readcon, & ! read query position control data
  query_set, & ! set up query array
  query_delete ! deallocate query array


  type, public :: queryset_t
     type(posveclis_t) :: posl !< list of positions
     real(kr4), dimension(:), allocatable :: en !< energy at position
     real(kr4), dimension(:), allocatable :: res !< local variable
     real(kr4), dimension(:), allocatable :: rest !< local variable
     real(kr4), dimension(:), allocatable :: resa !< local variable
     real(kr4), dimension(:), allocatable :: resb !< local variable
     integer(ki4):: np !< number of entries in posl
     integer(ki4):: ni !< index of res
     integer(ki4):: nj !< index of res
     integer(ki4):: nk !< index of res
     integer(ki4):: nie !< same as nextra
     character(len=80) :: which !< which quantity is required
  end type queryset_t

  type, public :: querydata_t
     real(kr4):: scalen !< units for scale energy (1=MeV)
     integer(ki4):: nptype !< type of query vector transform
     integer(ki4):: ntype !< type of query vector transform
     integer(ki4):: nenergy !< local variable
     integer(ki4):: nextra !< extra number of energies (eg if \f$ \gamma \f$ in database)
     integer(ki4):: ninline !< local variable
     integer(ki4):: nlong !< local variable
     integer(ki4):: nlat !< local variable
     type(tfmdata_t) :: tfmdata !< query position transform numeric controls
  end type querydata_t

!public variables
  integer :: nqury   !< output channel for query position data
  integer :: ninr   !< input channel for query density data

! private types

! private variables
  character(*), parameter :: m_name='query_m' !< module name
  integer   :: status   !< error status
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  character(len=80) :: ibuf1 !< local variable
  character(len=80) :: ibuf2 !< local variable
!! logical :: unitused !< flag to test unit is available
  logical :: iltest !< logical flag


  contains
!---------------------------------------------------------------------
!> read in (vtk) query density data
subroutine query_readv(self,infile)

  !! arguments
  type(queryset_t), intent(inout) :: self   !< query position data
  character(*),intent(in) :: infile !< name of input file
  logical :: isnumb !< local variable
  external isnumb


  !! local
  character(*), parameter :: s_name='query_readv' !< subroutine name
  integer(ki4) :: inpos !< number of position vectors read
  integer(ki4) :: inres !< number of result vectors read
  integer(ki4) :: ip !< address in array
  type(queryset_t) :: zquery   !< local query position data

  call position_readlis(zquery%posl,infile,ninr)

  ! need to check zquery%posl matches self%posl
  ! only check have same number of data
  inpos=zquery%posl%np

  ! skip to POINT_DATA, ignore CELL and CELL_TYPES
  do
     read(ninr,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,10,error_fatal,'Error reading header data')
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:10)=='POINT_DATA') then
           iltest=isnumb(ibuf2,inres,11)
           if (inres/=inpos) then
              call log_error(m_name,s_name,12,error_fatal,'Numbers of data disagree')
           end if
           exit
        end if
     end if
  end do
  self%np=inres
  !! check res storage
  !      if(self%np>0) then
  !             deallocate(self%res)
  !         allocate(self%res(self%np), &
  !     &   stat=status)
  !       !! check successful allocation
  !         if(status/=0) then
  !            call log_error(m_name,s_name,22,error_fatal,'Unable to allocate memory')
  !         end if
  !      else
  !         call log_error(m_name,s_name,23,error_fatal,'No res data')
  !      end if
  ! skip 2 records
  do j=1,2
     read(ninr,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,30,error_fatal,'Error reading skip data')
     else
        cycle
     end if
  end do

  do k=1,self%nk
     do j=1,self%nj
        do i=1,self%ni
           ip=i+(j-1)*self%ni+(k-1)*self%ni*self%nj
           read(ninr,*,iostat=status) self%res(ip)
           if(status/=0) then
              if(status<0) then
                 call log_error(m_name,s_name,14,error_fatal,'EOF encountered')
              else
                 call log_error(m_name,s_name,15,error_fatal,'Error in res read')
              end if
           end if
        end do
     end do
  end do

  call log_value("number of res read ",self%np)

end subroutine query_readv
!---------------------------------------------------------------------
!> read in query position data
subroutine query_read(self,querydata,infile)

  !! arguments
  type(queryset_t), intent(out) :: self   !< query data set
  type(querydata_t), intent(out) :: querydata   !< query position transform numeric controls
  character(*),intent(in) :: infile !< name of input file

  !! local
  character(*), parameter :: s_name='query_read' !< subroutine name
  integer(ki4):: i1 !< local variable
  integer(ki4):: i2 !< local variable
  integer(ki4):: ji !< local variable
  integer(ki4):: inline !< local variable
  integer(ki4):: ien !< local variable
  integer(ki4):: ilines !< local variable
  real(kr4) :: zscalen !< local variable


  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nqury=i exit end if end do

  !! open file
  call misc_getfileunit(nqury)
  open(unit=nqury,file=infile,status='OLD',form='FORMATTED',iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,1,error_fatal,'Error opening query data file')
  else
     call log_error(m_name,s_name,2,log_info,'Query data file opened ')
  end if

  call query_readcon(querydata,nqury)

  ! energy levels
  !! allocate energy storage
  if(querydata%nenergy>0) then
     allocate(self%en(querydata%nenergy), &
 &   stat=status)
     !! check successful allocation
     if(status/=0) then
        call log_error(m_name,s_name,3,error_fatal,'Unable to allocate memory')
     end if
  else
     call log_error(m_name,s_name,4,error_fatal,'No energy data')
  end if

  querydata_pos_type: select case (querydata%nptype)
  case(1)
     ! regular lattice in sphericals
     inline=querydata%ninline
     ien=querydata%nenergy
     ilines=(ien+inline-1)/inline

     do j=1,ilines
        i1=(j-1)*inline+1
        i2=min(ien, j*inline)
        read(nqury,*,iostat=status) (self%en(ji),ji=i1,i2)
        if(status/=0) then
           call log_error(m_name,s_name,5,error_fatal,'Error reading query energies')
        end if
     end do

  case(2)
     ien=2
     ! regular lattice in cartesians use only en(1) and en(2)
     read(nqury,*,iostat=status) (self%en(ji),ji=1,ien)
     if(status/=0) then
        call log_error(m_name,s_name,6,error_fatal,'Error reading query energies')
     end if
  case(3)
     ien=2
     ! random distribution in cartesians use only en(1) and en(2)
     read(nqury,*,iostat=status) (self%en(ji),ji=1,ien)
     if(status/=0) then
        call log_error(m_name,s_name,7,error_fatal,'Error reading query energies')
     end if
  case default
     call log_error(m_name,s_name,8,error_fatal,'Invalid type of position query')
  end select querydata_pos_type

  ! energy gets scaled
  zscalen=querydata%scalen
  do j=1,ien
     self%en(j)=zscalen*self%en(j)
  end do

  querydata_type: select case (querydata%ntype)
  case(1)
     ! regular density
     self%which='density'
  case(2)
     ! flux density
     self%which='adaptive density'
  case(3)
     ! flux density
     self%which='flux density'
  case default
     call log_error(m_name,s_name,20,error_fatal,'Invalid type of query')
  end select querydata_type

  self%ni=querydata%nenergy
  self%nj=querydata%nlat
  self%nk=querydata%nlong
  self%nie=querydata%nextra

end subroutine query_read
!---------------------------------------------------------------------
!> set up in query query set data
subroutine query_set(self,querydata)

  !! arguments
  type(queryset_t), intent(inout) :: self   !< query position set
  type(querydata_t), intent(in) :: querydata   !< query position transform numeric controls


  !! local
  character(*), parameter :: s_name='query_set' !< subroutine name
  integer(ki4) :: ip !< local variable
  integer(ki4) :: isize    !< local variable
  real(kr4) :: phi !< local variable
  real(kr4) :: cosp !< local variable
  real(kr4) :: sinp !< local variable
  real(kr4) :: theta !< local variable
  real(kr4) :: sinth !< local variable
  real(kr4) :: costh !< local variable
  real(kr4) :: zdeli !< local variable
  real(kr4) :: zdelj !< local variable
  real(kr4) :: zdelk !< local variable
  real(kr4) :: zwidi !< local variable
  real(kr4) :: zwidj !< local variable
  real(kr4) :: zwidk !< local variable
  real :: randlx !< local variable

  isize=querydata%nenergy*querydata%nlat*querydata%nlong
  ! energy levels
  !! allocate energy storage
  if(isize>0) then
     !        allocate(self%posl(querydata%nenergy,querydata%nlat,querydata%nlong), &
     allocate(self%posl%pos(isize), &
 &   stat=status)
     !! check successful allocation
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Unable to allocate memory')
     end if
  else
     call log_error(m_name,s_name,2,error_fatal,'No energy data')
  end if

  ! allocate space for results
  !        allocate(self%res(querydata%nenergy,querydata%nlat,querydata%nlong), &
  allocate(self%res(isize), self%rest(isize),&
 &stat=status)
  !! check successful allocation
  if(status/=0) then
     call log_error(m_name,s_name,10,error_fatal,'Unable to allocate memory')
  end if

  !        allocate(self%resa(querydata%nenergy,querydata%nlat,querydata%nlong), &
  allocate(self%resa(isize), self%resb(isize),&
 &stat=status)
  !! check successful allocation
  if(status/=0) then
     call log_error(m_name,s_name,11,error_fatal,'Unable to allocate memory')
  end if

  querydata_type: select case (querydata%nptype)
  case(1)
     ! regular lattice in sphericals
     self%np=isize
     self%posl%np=isize
     ! set up sample point array
     do k=1,querydata%nlong
        phi=(k-1)*2*const_pi/querydata%nlong
        cosp=cos(phi)
        sinp=sin(phi)
        do j=1,querydata%nlat
           theta=(j-1)*const_pi/(querydata%nlat-1)
           sinth =  sin(theta)
           costh =  cos(theta)
           do i=1,querydata%nenergy
              ip=i+(j-1)*self%ni+(k-1)*self%ni*self%nj
              self%posl%pos(ip)%posvec(1) = self%en(i)*sinth*cosp
              self%posl%pos(ip)%posvec(2) = self%en(i)*sinth*sinp
              self%posl%pos(ip)%posvec(3) = self%en(i)*costh
           end do
        end do
     end do

  case(2)
     ! regular lattice in cartesians (use only en(1) and en(2))
     self%np=isize
     self%posl%np=isize
     ! set up sample point array
     zdelk=(self%en(2)-self%en(1))/querydata%nlong
     zdelj=(self%en(2)-self%en(1))/querydata%nlat
     zdeli=(self%en(2)-self%en(1))/querydata%nenergy
     do k=1,querydata%nlong
        do j=1,querydata%nlat
           do i=1,querydata%nenergy
              ip=i+(j-1)*self%ni+(k-1)*self%ni*self%nj
              self%posl%pos(ip)%posvec(1) = self%en(1)+i*zdeli
              self%posl%pos(ip)%posvec(2) = self%en(1)+j*zdelj
              self%posl%pos(ip)%posvec(3) = self%en(1)+k*zdelk
           end do
        end do
     end do

  case(3)
     ! random samples in cartesians (use only en(1) and en(2))
     self%np=isize
     self%posl%np=isize
     ! set up sample point array
     zwidk=(self%en(2)-self%en(1))
     zwidj=(self%en(2)-self%en(1))
     zwidi=(self%en(2)-self%en(1))
     do k=1,querydata%nlong
        do j=1,querydata%nlat
           do i=1,querydata%nenergy
              ip=i+(j-1)*self%ni+(k-1)*self%ni*self%nj
              self%posl%pos(ip)%posvec(1) = self%en(1)+randlx()*zwidi
              self%posl%pos(ip)%posvec(2) = self%en(1)+randlx()*zwidj
              self%posl%pos(ip)%posvec(3) = self%en(1)+randlx()*zwidk
           end do
        end do
     end do

  end select querydata_type

  !      deallocate(self%en)


end subroutine query_set

!!---------------------------------------------------------------------
!> output query position vector
subroutine query_writev(self,numerics,kchar,kout)

  !! arguments
  type(queryset_t), intent(in) :: self   !< query position data
  type(numerics_t), intent(in) :: numerics   !< numeric controls
  character(*),intent(in):: kchar !< control output (dummy)
  integer, intent(in) :: kout   !< output channel for query position data


  !! local
  character(*), parameter :: s_name='query_writev' !< subroutine name
  type(posvecl_t) :: zpos !< local variable

  !! plot list of all positions
  write(kout,'(''DATASET UNSTRUCTURED_GRID'')')
  write(kout,'(''POINTS '',I8, '' float'')') self%np
  do j=1,self%np
     pos_type: select case (kchar)
     case('real space')
        zpos=position_invqtfm(self%posl%pos(j),numerics%geobj_coord_tfm)
     case default
        !     case('quantised space')
        zpos=self%posl%pos(j)
     end select pos_type
     call position_writev(zpos,kout)
  end do
  write(kout, '('' '')')

  write(kout,'(''CELLS 1 '',I8)') self%np+1
  write(kout,'(I8)') self%np,  (i-1,i=1,self%np)
  write(kout,'(''CELL_TYPES 1'')')
  write(kout,'(''2'')')
  write(kout, '('' '')')

  !! density
  write(kout,'(''POINT_DATA'',I8)') self%np
  write(kout,'(''SCALARS density float 1'')')
  write(kout,'(''LOOKUP_TABLE default'')')
  plot_type: select case (kchar)
  case('restored')
     write(kout,'((g14.7,1x))') ( self%rest(i),  i=1,self%np)
  case default
     write(kout,'((g14.7,1x))') ( self%res(i),  i=1,self%np)
  end select plot_type


end subroutine query_writev

subroutine query_readcon(querydata,kin)

  !! arguments
  type(querydata_t), intent(out) :: querydata   !< query position transform numeric controls
  integer(ki4) :: kin  !< local variable

  !! local
  character(*), parameter :: s_name='query_readcon' !< subroutine name
  real(kr4), dimension(3,3):: query_matrix  !< local variable
  real(kr4), dimension(3):: query_scale  !< local variable
  real(kr4), dimension(3):: query_offset  !< local variable
  character(len=80) :: query_file !< local variable
  real(kr4):: query_scale_energy  !< local variable
  integer(ki4):: query_transform  !< local variable
  integer(ki4):: nquery_position_type !< local variable
  integer(ki4):: nquery_type !< local variable
  integer(ki4):: nquery_energies !< local variable
  integer(ki4):: nquery_extra_energies !< local variable
  integer(ki4):: nquery_energies_per_line !< local variable
  integer(ki4):: nquery_in_longitude !< local variable
  integer(ki4):: nquery_in_latitude !< local variable

  !! query position parameters
  namelist /querypositionparameters/ &
 &query_file , &
 &nquery_position_type , &
 &nquery_type , &
 &nquery_energies , &
 &nquery_extra_energies , &
 &nquery_energies_per_line , &
 &nquery_in_longitude , &
 &nquery_in_latitude , &
 &query_matrix , &
 &query_scale , &
 &query_scale_energy , &
 &query_offset , &
 &query_transform

  !! set default query position parameters
  nquery_position_type=1
  nquery_type=1
  nquery_energies=150
  nquery_extra_energies=0
  nquery_energies_per_line=1
  nquery_in_longitude=10
  nquery_in_latitude=10
  ! default identity
  query_transform=1
  query_matrix(1,:)=(/1.,  0.,  0. /)
  query_matrix(2,:)=(/0.,  1.,  0. /)
  query_matrix(3,:)=(/0.,  0.,  1. /)
  query_scale=(/1.,  1.,  1. /)
  query_scale_energy=1.
  query_offset=(/0.,  0.,  0. /)

  !!read query position parameters
  read(kin,nml=querypositionparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading query position parameters")'
     call log_getunit(ilog)
     write(ilog,nml=querypositionparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading query position parameters')
  end if


  !! check for valid data
  if(nquery_position_type<=0) &
 &call log_error(m_name,s_name,2,error_fatal,'nquery_position_type must be > 0')
  if(nquery_energies<=0) &
 &call log_error(m_name,s_name,3,error_fatal,'nquery_energies must be > 0')
  if(nquery_energies_per_line<=0) &
 &call log_error(m_name,s_name,4,error_fatal,'nquery_energies_per_line must be > 0')
  if(nquery_in_longitude<=0) &
 &call log_error(m_name,s_name,5,error_fatal,'nquery_in_longitude must be > 0')
  if(nquery_in_latitude<=0) &
 &call log_error(m_name,s_name,6,error_fatal,'nquery_in_latitude must be > 0')
  if(all(query_scale<=0)) &
 &call log_error(m_name,s_name,7,error_fatal,'query_scale must be > 0')
  if(query_transform<1)  &
 &call log_error(m_name,s_name,8,error_fatal,'query_transform must be > 0')
  if(nquery_type<=0) &
 &call log_error(m_name,s_name,9,error_fatal,'nquery_type must be > 0')
  if(nquery_extra_energies<=0) &
 &call log_error(m_name,s_name,10,error_fatal,'nquery_extra_energies must be > 0')

  !! store values
  querydata%nptype=nquery_position_type
  querydata%ntype=nquery_type
  querydata%nenergy=nquery_energies
  querydata%nextra=nquery_extra_energies
  querydata%ninline=nquery_energies_per_line
  querydata%nlong=nquery_in_longitude
  querydata%nlat=nquery_in_latitude+1
  querydata%scalen=query_scale_energy
  querydata%tfmdata%matrix=query_matrix
  querydata%tfmdata%scale=query_scale
  querydata%tfmdata%offset=query_offset
  querydata%tfmdata%ntfm=query_transform

end subroutine query_readcon
!---------------------------------------------------------------------
!! delete queryset_t
subroutine query_delete(self)

  !! arguments
  type(queryset_t), intent(inout) :: self !< queryset data

  deallocate(self%res)
  deallocate(self%rest)
  deallocate(self%en)
  deallocate(self%posl%pos)
  deallocate(self%resa)
  deallocate(self%resb)

end subroutine query_delete

end module query_m
