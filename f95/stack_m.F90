module stack_m

  use const_kind_m
  use const_numphys_h
  use log_m

  implicit none
  private

! public subroutines
  public :: stack_add, & ! add to stack
  stack_get, & ! get data from top of stack and remove
  stack_empty, & ! is stack empty
  stack_init, & ! set up stack arrays
  stack_delete ! deallocate stack arrays


  interface stack_add
   module procedure stack_add_ki4
   module procedure stack_add_kr8
  end interface

  interface stack_get
   module procedure stack_get_ki4
   module procedure stack_get_kr8
  end interface

  interface stack_init
   module procedure stack_init_ki4
   module procedure stack_init_kr8
  end interface

  interface stack_delete
   module procedure stack_delete_ki4
   module procedure stack_delete_kr8
  end interface

  type, public :: stackset_t
     integer(ki4), dimension(:,:), allocatable :: iarra !< integer array one
     integer(ki4), dimension(:,:), allocatable :: iarrb !< integer array two
     real(kr8), dimension(:,:), allocatable :: arra !< real array one
     real(kr8), dimension(:,:), allocatable :: arrb !< real array two
     integer(ki4) :: head=0 !< head of stack
     integer(ki4) :: ndim=1 !< first index size of stack
     integer(ki4) :: nsdim=30 !< starting dimension of stack
     integer(ki4) :: ndfac=3 !< factor to increase dimension
     character(len=80) :: which !< which type of quantity on stack
  end type stackset_t

!public variables

! private types

! private variables
  character(*), parameter :: m_name='stack_m' !< module name
  type(stackset_t) :: stack !< stackset data
  integer   :: status   !< error status
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter

  contains
!---------------------------------------------------------------------
!> add to stack
subroutine stack_add_ki4(vself,kdim,ksave)

  !! arguments
  integer(ki4), intent(in) :: vself   !< stack type
  integer(ki4), intent(in) :: kdim   !< stack dimension
  integer(ki4), intent(in), dimension(kdim) :: ksave !< stack entry

  !! local
  character(*), parameter :: s_name='stack_add_ki4' !< subroutine name
  integer(ki4) :: isdim !< old size

  if (stack%which/='ki4'.OR.kdim/=stack%ndim) then
     call log_error(m_name,s_name,1,error_fatal,'stack of wrong type')
  end if

  stack%head=stack%head+1
  if (stack%head>stack%nsdim) then
     !! increase stack size
     isdim=stack%nsdim
     stack%nsdim=stack%ndfac*isdim
     if (allocated(stack%iarra)) then
        allocate(stack%iarrb(kdim,stack%nsdim), stat=status)
        call log_alloc_check(m_name,s_name,1,status)
        stack%iarrb(:,1:isdim)=stack%iarra(:,1:isdim)
        deallocate(stack%iarra)
     else
        allocate(stack%iarra(kdim,stack%nsdim), stat=status)
        call log_alloc_check(m_name,s_name,2,status)
        stack%iarra(:,1:isdim)=stack%iarrb(:,1:isdim)
        deallocate(stack%iarrb)
     end if
  end if

  if (allocated(stack%iarra)) then
     stack%iarra(:,stack%head)=ksave(:)
  else
     stack%iarrb(:,stack%head)=ksave(:)
  end if

end subroutine stack_add_ki4
subroutine stack_add_kr8(vself,kdim,ksave)

  !! arguments
  real(kr8), intent(in) :: vself   !< stack type
  integer(ki4), intent(in) :: kdim   !< stack dimension
  integer(ki4), intent(in), dimension(kdim) :: ksave !< stack entry

  !! local
  character(*), parameter :: s_name='stack_add_kr8' !< subroutine name
  integer(ki4) :: isdim !< old size

  if (stack%which/='kr8'.OR.kdim/=stack%ndim) then
     call log_error(m_name,s_name,1,error_fatal,'stack of wrong type')
  end if

  stack%head=stack%head+1
  if (stack%head>stack%nsdim) then
     !! increase stack size
     isdim=stack%nsdim
     stack%nsdim=stack%ndfac*isdim
     if (allocated(stack%arra)) then
        allocate(stack%arrb(kdim,stack%nsdim), stat=status)
        call log_alloc_check(m_name,s_name,1,status)
        stack%arrb(:,1:isdim)=stack%arra(:,1:isdim)
        deallocate(stack%arra)
     else
        allocate(stack%arra(kdim,stack%nsdim), stat=status)
        call log_alloc_check(m_name,s_name,2,status)
        stack%arra(:,1:isdim)=stack%arrb(:,1:isdim)
        deallocate(stack%arrb)
     end if
  end if

  if (allocated(stack%arra)) then
     stack%arra(:,stack%head)=ksave(:)
  else
     stack%arrb(:,stack%head)=ksave(:)
  end if

end subroutine stack_add_kr8
!---------------------------------------------------------------------
!> get data from top of stack and remove
subroutine stack_get_ki4(vself,kdim,ksave)

  !! arguments
  integer(ki4), intent(inout) :: vself   !< stack type
  integer(ki4), intent(in) :: kdim   !< stack dimension
  integer(ki4), intent(out), dimension(kdim) :: ksave !< stack entry

  !! local
  character(*), parameter :: s_name='stack_get_ki4' !< subroutine name

  if (stack%head<=0) then
     vself=stack%head-1
     call log_error(m_name,s_name,1,error_warning,'ki4 stack found empty')
     return
  else
     vself=1
     if (stack%which/='ki4') then
        call log_error(m_name,s_name,2,error_fatal,'stack of wrong type')
     end if
  end if

  if (allocated(stack%iarra)) then
     ksave=stack%iarra(:,stack%head)
  else
     ksave=stack%iarrb(:,stack%head)
  end if

  stack%head=stack%head-1

end subroutine stack_get_ki4
subroutine stack_get_kr8(vself,kdim,ksave)

  !! arguments
  real(kr8), intent(inout) :: vself   !< stack type
  integer(ki4), intent(in) :: kdim   !< stack dimension
  integer(ki4), intent(out), dimension(kdim) :: ksave !< stack entry

  !! local
  character(*), parameter :: s_name='stack_get_kr8' !< subroutine name

  if (stack%head<=0) then
     vself=stack%head-1
     call log_error(m_name,s_name,1,error_warning,'kr8 stack found empty')
     return
  else
     vself=1
     if (stack%which/='kr8') then
        call log_error(m_name,s_name,2,error_fatal,'stack of wrong type')
     end if
  end if

  if (allocated(stack%arra)) then
     ksave=stack%arra(:,stack%head)
  else
     ksave=stack%arrb(:,stack%head)
  end if

  stack%head=stack%head-1

end subroutine stack_get_kr8
!---------------------------------------------------------------------
! is stack empty
function stack_empty(vself,kdim)

  !! arguments
  integer(ki4), intent(in) :: vself   !< stack type
  integer(ki4), intent(in) :: kdim   !< stack dimension

  !! local
  character(*), parameter :: s_name='stack_empty' !< subroutine name
  logical :: stack_empty !< function return

  stack_empty=(stack%head<=0)

end function stack_empty
!---------------------------------------------------------------------
!> set up stack array
subroutine stack_init_ki4(vself,kdim)

  !! arguments
  integer(ki4), intent(in) :: vself   !< stack position type
  integer(ki4), intent(in) :: kdim   !< stack dimension

  !! local
  character(*), parameter :: s_name='stack_init' !< subroutine name

  stack%ndim=kdim
  stack%which='ki4'
  allocate(stack%iarra(kdim,stack%nsdim), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  call log_error(m_name,s_name,2,log_info,'ki4 stack initialised')

end subroutine stack_init_ki4
subroutine stack_init_kr8(vself,kdim)

  !! arguments
  real(kr8), intent(in) :: vself   !< stack position type
  integer(ki4), intent(in) :: kdim   !< stack dimension

  !! local
  character(*), parameter :: s_name='stack_init' !< subroutine name

  stack%ndim=kdim
  stack%which='kr8'
  allocate(stack%arra(kdim,stack%nsdim), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  call log_error(m_name,s_name,2,log_info,'kr8 stack initialised')

end subroutine stack_init_kr8
!---------------------------------------------------------------------
!! delete stack storage
subroutine stack_delete_ki4(vself)

  !! arguments
  integer(ki4), intent(in) :: vself !< stack type

  stack%head=0
  stack%ndim=1
  stack%which=' '
  if (allocated(stack%iarra)) deallocate(stack%iarra)
  if (allocated(stack%iarrb)) deallocate(stack%iarrb)

end subroutine stack_delete_ki4
subroutine stack_delete_kr8(vself)

  !! arguments
  real(kr8), intent(in) :: vself !< stack type

  stack%head=0
  stack%ndim=1
  stack%which=' '
  if (allocated(stack%arra)) deallocate(stack%arra)
  if (allocated(stack%arrb)) deallocate(stack%arrb)

end subroutine stack_delete_kr8

end module stack_m
