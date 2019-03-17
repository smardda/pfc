module termplane_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use termplane_h

  implicit none
  private

! public subroutines
  public :: &
  termplane_readcon, & !< read termination planes control data
  termplane_delete !< deallocate termination planes arrays

!public variables

! private types

! private variables
  character(*), parameter :: m_name='termplane_m' !< module name
  integer   :: status   !< error status
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  character(len=80) :: ibuf1 !< local variable
  character(len=80) :: ibuf2 !< local variable
  logical :: unitused !< flag to test unit is available
  logical :: iltest !< logical flag
  integer(ki4), parameter :: maxinp=100 !< max number of inputs


  contains
!---------------------------------------------------------------------
!> read termination planes control data
subroutine termplane_readcon(self,kin)

  !! arguments
  type(termplane_t), intent(out) :: self   !< termination planes numeric controls
  integer :: kin !< input unit number

  !! local
  character(*), parameter :: s_name='termplane_readcon' !< subroutine name
  integer(ki4), dimension(maxinp):: termplane_direction  !< namelist variable
  integer(ki4), dimension(maxinp):: termplane_intersection  !< namelist variable
  integer(ki4), dimension(maxinp):: termplane_extremum  !< namelist variable
  real(kr8), dimension(maxinp):: termplane_position  !< namelist variable
  real(kr8), dimension(maxinp):: termplane_condition  !< namelist variable
  integer(ki4), dimension(maxinp):: termplane_condition_dir  !< namelist variable
  integer(ki4):: termplane_max  !< namelist variable
  integer(ki4):: nactive  !< local variable

  !! termplane position parameters
  namelist /termplaneparameters/ &
 &termplane_direction, &
 &termplane_intersection, &
 &termplane_extremum, &
 &termplane_position, &
 &termplane_condition, &
 &termplane_condition_dir, &
 &termplane_max

  !! set default termplane position parameters
  termplane_position=0
  termplane_condition=0
  termplane_condition_dir=0
  termplane_direction=0
  termplane_intersection=0
  termplane_extremum=0
  termplane_max=0

  !!read termplane position parameters
  read(kin,nml=termplaneparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading termplane parameters")'
     call log_getunit(ilog)
     write(ilog,nml=termplaneparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading termplane parameters')
  end if

  nactive=0
  do j=1,maxinp
     !! check for valid data and count
     if(abs(termplane_direction(j))+abs(termplane_intersection(j))+abs(termplane_extremum(j))==0) exit
     if(abs(termplane_direction(j))>3) then
        call log_error(m_name,s_name,2,error_warning,'termplane_direction must be <= 3')
     else if (termplane_direction(j)/=0) then
        nactive=nactive+1
     end if
     if(abs(termplane_intersection(j))>3) then
        call log_error(m_name,s_name,3,error_warning,'termplane_intersection must be <= 3')
     else if (termplane_intersection(j)/=0) then
        nactive=nactive+1
     end if
     if(abs(termplane_extremum(j))>3) then
        call log_error(m_name,s_name,4,error_warning,'termplane_extremum must be <= 3')
     else if (termplane_extremum(j)/=0) then
        nactive=nactive+1
     end if
     if(abs(termplane_condition_dir(j))>3) then
        call log_error(m_name,s_name,4,error_warning,'termplane_condition_dir must be <= 3')
     end if
  end do

  ! allow one extra space in arrays for possible dynamical test
  allocate(self%termplanedir(nactive+1,5),self%termplane(nactive+1,2),&
 &self%termstore(nactive+1,3),stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  self%ntermplane=nactive
  self%ntermactive=nactive
  self%termplanedir=0
  self%termstore=0
  self%termplane=0

  nactive=0
  do j=1,maxinp
     if(abs(termplane_direction(j))+abs(termplane_intersection(j))+abs(termplane_extremum(j))==0) exit
     !! store values
     if(abs(termplane_direction(j))>3.OR.termplane_direction(j)==0) then
     else
        nactive=nactive+1
        self%termplanedir(nactive,1)=abs(termplane_direction(j))
        self%termplanedir(nactive,2)=sign(1,termplane_direction(j))
        self%termplanedir(nactive,3)=0
        self%termplane(nactive,1)=termplane_position(j)
     end if
     if(abs(termplane_intersection(j))>3.OR.termplane_intersection(j)==0) then
     else
        nactive=nactive+1
        self%termplanedir(nactive,1)=abs(termplane_intersection(j))
        self%termplanedir(nactive,2)=sign(1,termplane_intersection(j))
        self%termplanedir(nactive,3)=1
        self%termplane(nactive,1)=termplane_position(j)
     end if
     if(abs(termplane_extremum(j))>3.OR.termplane_extremum(j)==0) then
     else
        nactive=nactive+1
        self%termplanedir(nactive,1)=abs(termplane_extremum(j))
        self%termplanedir(nactive,2)=sign(1,termplane_extremum(j))
        self%termplanedir(nactive,3)=2
        self%termplane(nactive,1)=termplane_position(j)
     end if
     if(abs(termplane_condition_dir(j))>0) then
        self%termplanedir(nactive,3)=3
        self%termplanedir(nactive,4)=termplane_condition_dir(j)
        self%termplanedir(nactive,5)=sign(1,termplane_condition_dir(j))
        self%termplane(nactive,2)=termplane_condition(j)
     else
        ! ensure non-zero direction when scaling
        self%termplanedir(nactive,4)=self%termplanedir(nactive,1)
        self%termplanedir(nactive,5)=self%termplanedir(nactive,2)
     end if
  end do
  ! dynamical values (skylight for lower or upper divertor, here lower)
  self%termplanedir(nactive+1,1)=2
  self%termplanedir(nactive+1,2)=1

end subroutine termplane_readcon
!---------------------------------------------------------------------
!> delete termplane_t
subroutine termplane_delete(self)

  !! arguments
  type(termplane_t), intent(inout) :: self !< termplane data

  deallocate(self%termplane)
  deallocate(self%termplanedir)

end subroutine termplane_delete

end module termplane_m
