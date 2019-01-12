module moutfile_m

  use const_kind_m
  use log_m
  use misc_m
  use mcontrol_h
  use date_time_m

  implicit none
  private

! public subroutines
  public :: &
 &moutfile_init, &
 &moutfile_read, &
 &moutfile_close

! public variables

! private variables
  character(*), parameter :: m_name='moutfile_m' !< module name
  integer :: nout !< output file unit
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: idum !< dummy integer
  integer   :: status   !< error status
  character(len=80) :: ibuff !< buffer for input/output

  contains
!---------------------------------------------------------------------
!> open output file
subroutine moutfile_init(fout,numerics,descriptor,kout)

  !! arguments
  character(len=*), intent(in) :: fout !< file name root
  type(mnumerics_t), intent(in)  :: numerics  !< numerical control parameters
  character(len=*), intent(in) :: descriptor !< dataset descriptor
  integer, intent(inout) :: kout   !< unit number

  !! local
  character(*), parameter :: s_name='moutfile_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  !! get file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kout=i exit end if end do

  call misc_getfileunit(kout)
  open(unit=kout,file=trim(fout),status='new',form='FORMATTED',iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,1,error_fatal,'Error opening data structure file')
  else
     call log_error(m_name,s_name,2,log_info,'data structure file opened')
  end if
  nout=kout

  !! header information

  write(kout,*,iostat=status) 'ireq'
  call log_write_check(m_name,s_name,1,status)
  write(kout,*,iostat=status) numerics%ireq
  call log_write_check(m_name,s_name,2,status)

end  subroutine moutfile_init
!---------------------------------------------------------------------
!> read output file
subroutine moutfile_read(infile,pvar,kin)

  !! arguments
  character(len=*), intent(in) :: infile !< file name
  real(kr8), intent(out) :: pvar !< real variable to read
  integer, intent(inout) :: kin   !< unit number

  !! local
  character(*), parameter :: s_name='moutfile_read' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  !! get file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kin=i exit end if end do

  !! open file
  call misc_getfileunit(kin)
  open(unit=kin,file=trim(infile),status='OLD',form='FORMATTED',iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,1,error_fatal,'Error opening magnetic field file')
  else
     call log_error(m_name,s_name,2,log_info,'magnetic field file opened')
  end if

  !! header information

  read(kin,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,10,status)
  read(kin,*,iostat=status) pvar
  call log_read_check(m_name,s_name,11,status)

end  subroutine moutfile_read
!---------------------------------------------------------------------
!> close output files
subroutine moutfile_close

  close(unit=nout)

end  subroutine moutfile_close

end module moutfile_m
