module clock_m

  use const_kind_m
  use date_time_m
  use log_m

  implicit none
  private

! public subroutines
  public :: &
  clock_init ,   & ! initialise clock storage
  clock_start,   & ! start a clock
  clock_restart, & ! restart a clock
  clock_stop,    & ! stop a clock
  clock_log  ,   & ! print time on clock to log
  clock_summary, & ! print clock time summary
  clock_delete

! public variables


! private variables
  character(*), parameter :: m_name='clock_m' !< module name
  integer :: status    !< error status flag
  integer(ki4) :: noclock   !< number of active clocks
  integer(ki4) :: maxclocks !< max no of clocks
  integer(ki4) :: iclock    !< index of current clock
  integer(ki4)  :: i         !< loop index
  integer(ki4), dimension(:), allocatable :: mapclock !< cpu clock numbers
  character(len=80), dimension(:), allocatable :: clocknames !< clock names
  real(kr8), dimension(:), allocatable :: starttime !< start time of cpu clock
  real(kr8), dimension(:), allocatable :: elapsetime !<  time on cpu clock
  real(kr8) :: timeoff !< time clock stopped
  logical :: clockfound !< flag for clock search
  contains
!---------------------------------------------------------------------
!> initialise cpu clocks
subroutine clock_init(noclocks)

  !! arguments
  integer(ki4), intent(in) :: noclocks !< number of clocks
  !! local
  character(*), parameter :: s_name='clock_init' !< subroutine name

  if(noclocks<1) then
     call log_error(m_name,s_name,1,error_warning,'cpu clock not enabled')
  else
     allocate(mapclock(noclocks), &
 &   clocknames(noclocks), &
 &   starttime(noclocks), &
 &   elapsetime(noclocks), &
     stat=status)
     if(status/=0)then
        call log_error(m_name,s_name,2,error_fatal,'allocation failed')
     end if
  end if
  maxclocks=noclocks
  noclock=0
  elapsetime=0.0

end subroutine clock_init
!---------------------------------------------------------------------
!> deallocate clocks
subroutine clock_delete
  deallocate(mapclock,clocknames,starttime,elapsetime)
end subroutine clock_delete
!---------------------------------------------------------------------
!> start cpu clock
subroutine clock_start(clockno,clockname)

  !! arguments
  integer(ki4), intent(in) :: clockno !< clock number
  character(*), intent(in) :: clockname !< clock name
  !! local
  character(*), parameter :: s_name='clock_start' !< subroutine name

  noclock=noclock+1
  mapclock(noclock)=clockno
  clocknames(noclock)=clockname

  call cpu_time(starttime(noclock))

end subroutine clock_start
!---------------------------------------------------------------------
!> stop cpu clock
subroutine clock_stop(clockno)

  !! arguments
  integer(ki4), intent(in) :: clockno !< clock number
  !! local
  character(*), parameter :: s_name='clock_stop' !< subroutine name

  !! get time
  call cpu_time(timeoff)

  !! find clock
  clockfound=.false.
  do i=1,noclock
     if(clockno==mapclock(i)) then
        clockfound=.true.
        iclock=i
        exit
     end if
  end do

  !! add time to elapse counter
  if(clockfound) then
     elapsetime(iclock)=elapsetime(iclock)+(timeoff-starttime(iclock))
  else
     call log_error(m_name,s_name,1,error_warning,'cpu clock not found')
     call log_value('clockno =',clockno)
  end if
end subroutine clock_stop
!---------------------------------------------------------------------
!> print time on cpu clock to log
subroutine clock_log(clockno)
  !! arguments
  integer(ki4), intent(in) :: clockno !< clock number
  !! local
  character(*), parameter :: s_name='clock_print' !< subroutine name

  !! find clock
  clockfound=.false.
  do i=1,noclock
     if(clockno==mapclock(i)) then
        clockfound=.true.
        iclock=i
        exit
     end if
  end do

  !! add time to log file
  if(clockfound) then
     call log_value(trim(clocknames(iclock)),elapsetime(iclock),'secs')
     print '(a, "=", t35,g12.5," secs")',trim(clocknames(iclock)),elapsetime(iclock)
  else
     call log_error(m_name,s_name,1,error_warning,'cpu clock not found')
     call log_value('clockno',clockno)
  end if
end subroutine clock_log
!---------------------------------------------------------------------
!> print cpu clock summary
subroutine clock_summary
  do i=1,noclock
     call clock_log(mapclock(i))
  end do

end subroutine clock_summary
!---------------------------------------------------------------------
!> restart cpu clock
subroutine clock_restart(clockno)

  !! arguments
  integer(ki4), intent(in) :: clockno !< clock number
  !! local
  character(*), parameter :: s_name='clock_restart' !< subroutine name

  call cpu_time(starttime(clockno))

end subroutine clock_restart

end module clock_m
