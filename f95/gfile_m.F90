module gfile_m

  use const_kind_m
  use log_m
  use misc_m

  implicit none
  private

! public subroutines
  public :: &
 &gfile_init, & !< find new unit and insert gnuplot header
 &gfile_iwrite, & !< write of labelled integer array pair
 &gfile_rwrite, & !< write of labelled real array pair
 &gfile_close !< close gnuplot print file on unit nprint
! public types


! private variables
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  character(*), parameter :: m_name='gfile_m' !< module name
  integer :: nprint  !< file unit for output
  integer :: status  !< status flag

  contains

!---------------------------------------------------------------------
!> find new unit and insert gnuplot header
subroutine gfile_init(fprint,descriptor,kprint)

  !! arguments
  character(len=*), intent(in) :: fprint !< file name root
  character(len=*), intent(in) :: descriptor !< dataset descriptor
  integer, intent(inout) :: kprint   !< unit number

  !! local
  character(*), parameter :: s_name='gfile_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  !! open file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kprint=i exit end if end do

  call misc_getfileunit(kprint)
  open(unit=kprint,file=trim(fprint)//'.gnu')

  !! write gnuplot header
  write(kprint,'(''#gnuplot data file'')')
  write(kprint,'(''# '',a)') descriptor
  write(kprint,'(''           '')')

  nprint=kprint

end subroutine gfile_init
!---------------------------------------------------------------------
!> write of labelled integer array pair
subroutine gfile_iwrite(xout,yout,nout,descriptor,channel)

  !! arguments
  integer(ki4), dimension(nout), intent(in) :: xout !< horizontal axis values
  integer(ki4), dimension(nout), intent(in) :: yout !< vertical axis values
  integer(ki4), intent(in) :: nout   !< number of entries in xout and yout
  character(len=*), intent(in) :: descriptor !< dataset descriptor
  integer, intent(in), optional :: channel   !< output channel for object data structure

  !! local
  character(*), parameter :: s_name='gfile_iwrite' !< subroutine name
  integer :: iprint !< unit number used

  if (present(channel)) then
     iprint=channel
  else
     iprint=nprint
  end if

  write(iprint,'(''#'',a)',iostat=status) descriptor
  call log_write_check(m_name,s_name,1,status)
  write(iprint,'(1P, 2(1X, G13.5))',iostat=status) (xout(k),yout(k),k=1,nout)
  call log_write_check(m_name,s_name,2,status)

end subroutine gfile_iwrite
!---------------------------------------------------------------------
!> write of labelled real array pair
subroutine gfile_rwrite(xout,yout,nout,descriptor,channel)

  !! arguments
  real(kr8), dimension(nout), intent(in) :: xout !< horizontal axis values
  real(kr8), dimension(nout), intent(in) :: yout !< vertical axis values
  integer(ki4), intent(in) :: nout   !< number of entries in xout and yout
  character(len=*), intent(in) :: descriptor !< dataset descriptor
  integer, intent(in), optional :: channel   !< output channel for object data structure

  !! local
  character(*), parameter :: s_name='gfile_rwrite' !< subroutine name
  integer :: iprint !< unit number used

  if (present(channel)) then
     iprint=channel
  else
     iprint=nprint
  end if

  write(iprint,'(''#'',a)',iostat=status) descriptor
  call log_write_check(m_name,s_name,1,status)
  write(iprint,'(1P, 2(1X, G13.5))',iostat=status) (xout(k),yout(k),k=1,nout)
  call log_write_check(m_name,s_name,2,status)

end subroutine gfile_rwrite
!---------------------------------------------------------------------
!> close gnuplot print file on unit nprint
subroutine gfile_close

  close(nprint)

end subroutine gfile_close

end module gfile_m
