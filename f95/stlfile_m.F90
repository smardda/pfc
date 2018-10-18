module stlfile_m

  use const_kind_m
  use log_m
  use misc_m

  implicit none
  private

! public subroutines
  public :: &
  stlfile_init, & !< initialise stl file
  stlfile_close  !< close stl file

! public types

  integer(ki4), parameter, public :: stlfile_made_up_data = 1 !< local variable

! private variables
  character(*), parameter :: m_name='stlfile_m' !< module name
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  character(len=80), save :: keepdesc !< save descriptor
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer :: nin  !< file unit for input
  integer :: nplot  !< file unit for output
  integer :: status  !< status flag
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> initialise stl file
subroutine stlfile_init(fplot,descriptor,kplot)

  !! arguments
  character(len=*), intent(in) :: fplot !< file name root
  character(len=*), intent(in) :: descriptor !< dataset descriptor
  integer, intent(inout) :: kplot   !< unit number

  !! local
  character(*), parameter :: s_name='stlfile_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  !! open file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kplot=i exit end if end do

  call misc_getfileunit(kplot)
  open(unit=kplot,file=trim(fplot)//'.stl')

  keepdesc=descriptor

  !! optionally write stl header
  if (trim(keepdesc)=='triangles') then
     write(kplot,'(''solid stlfile'')')
  end if

  nplot=kplot

end subroutine stlfile_init
!---------------------------------------------------------------------
!> close vis stl file on unit nplot
subroutine stlfile_close

  !! optionally write stl last line
  if (trim(keepdesc)=='triangles') then
     write(nplot,'(''endsolid stlfile'')')
  end if

  close(nplot)

end subroutine stlfile_close

end module stlfile_m
