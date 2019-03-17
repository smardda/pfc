module soutfile_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use misc_m
  use scontrol_h
  use date_time_m
  use geobj_m
  use geobjlist_h
  use geobjlist_m
  use smanal_h
  use smanal_m

  implicit none
  private

! public subroutines
  public :: &
 &soutfile_init, &
 &soutfile_write, &
 &soutfile_close
! public types
  integer :: nout !< output file unit

! private variables
  character(*), parameter :: m_name='soutfile_m' !< module name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer   :: status   !< error status
  integer(ki4) :: idum !< dummy integer

  contains
!---------------------------------------------------------------------
!> open output file
subroutine soutfile_init(file,timestamp)

  !! argument
  type(sfiles_t), intent(in) :: file !< file names
  type(date_time_t), intent(in) :: timestamp !< timestamp of run


  !! local
  character(*), parameter :: s_name='soutfile_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  !! get file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nout=i exit end if end do

  call misc_getfileunit(nout)
  open(unit=nout,file=trim(file%smanalout),status='new',form='FORMATTED',iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,1,error_fatal,'Error opening data structure file')
  else
     call log_error(m_name,s_name,2,log_info,'data structure file opened')
  end if

  !! header information
  write(nout,'(a)') trim(file%smanalout)
  write(nout,'(a)') trim(timestamp%long)
  write(nout,'(" vtk_input_file = ",a,/)') trim(file%vtkdata(1))
  write(nout,'(" vtk_powx_file = ",a,/)') trim(file%vtk)

end  subroutine soutfile_init
!---------------------------------------------------------------------
!> write output from smanal
subroutine soutfile_write(smanal)

  !! arguments
  type(smanal_t), intent(in) :: smanal !< smanal object data structure

  ! local variables

  !!write
  call smanal_write(smanal,nout)

end  subroutine soutfile_write
!---------------------------------------------------------------------
!> close output files
subroutine soutfile_close

  close(unit=nout)

end  subroutine soutfile_close

end module soutfile_m
