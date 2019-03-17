module poutfile_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use misc_m
  use pcontrol_h
  use date_time_m
  use position_h
  use position_m
  use spl2d_m
  use spl3d_m
  use powcal_h
  use powres_h
  use edgprof_h
  use edgprof_m
  use powcal_m

  implicit none
  private

! public subroutines
  public :: &
 &poutfile_init, &
 &poutfile_write, &
 &poutfile_close
! public types
  integer :: nout !< output file unit

! private variables
  character(*), parameter :: m_name='poutfile_m' !< module name
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
subroutine poutfile_init(file,timestamp)

  !! argument
  type(pfiles_t), intent(in) :: file !< file names
  type(date_time_t), intent(in) :: timestamp !< timestamp of run


  !! local
  character(*), parameter :: s_name='poutfile_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  !! get file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nout=i exit end if end do

  call misc_getfileunit(nout)
  open(unit=nout,file=trim(file%powcalout),status='new',form='FORMATTED',iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,1,error_fatal,'Error opening data structure file')
  else
     call log_error(m_name,s_name,2,log_info,'data structure file opened')
  end if

  !! header information
  write(nout,'(a)') trim(file%powcalout)
  write(nout,'(a)') trim(timestamp%long)
  write(nout,'(" vtk_input_file = ",a,/)') trim(file%vtkdata)
  write(nout,'(" vtkres_input_file = ",a,/)') trim(file%vtkres)
  write(nout,'(" hds_input_file = ",a,/)') trim(file%hdsdata)
  write(nout,'(" geoq_input_file = ",a,/)') trim(file%geoq)
  write(nout,'(" lau_input_file = ",a,/)') trim(file%laudata)

end  subroutine poutfile_init
!---------------------------------------------------------------------
!> write output from powcal
subroutine poutfile_write(powcal,timestamp)

  !! arguments
  type(powcal_t), intent(in) :: powcal !< powcal object data structure
  type(date_time_t), intent(in) :: timestamp !< timestamp of run

  ! local variables

  !!write
  call powcal_write(powcal,nout)

end  subroutine poutfile_write
!---------------------------------------------------------------------
!> close output files
subroutine poutfile_close

  close(unit=nout)

end  subroutine poutfile_close

end module poutfile_m
