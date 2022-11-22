module boutfile_m

  use const_kind_m
  use log_m
  use misc_m
  use bcontrol_m
  use date_time_m
  use position_h
  use fmesh_h
  use beq_h
  use beq_m
  use spl2d_m

  implicit none
  private

! public subroutines
  public :: &
 &boutfile_init, & !< find unit, open file and write header
 &boutfile_write, & !< write data structure
 &boutfile_getunit,  & !< get unit number
 &boutfile_close

! private variables
  character(*), parameter :: m_name='boutfile_m' !< module name
  integer, save :: nout=-1 !< output file unit
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: idum !< dummy integer
  integer   :: status   !< error status

  contains
!---------------------------------------------------------------------
!> open output file
subroutine boutfile_init(file,timestamp)

  !! argument
  type(bfiles_t), intent(in) :: file !< file names
  type(date_time_t), intent(in) :: timestamp !< timestamp of run


  !! local
  character(*), parameter :: s_name='boutfile_init' !< subroutine name
   !! logical :: unitused !< flag to test unit is available
!! get file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nout=i exit end if end do

  call misc_getfileunit(nout)
  open(unit=nout,file=trim(file%geoqout),status='new',form='FORMATTED',iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,1,error_fatal,'Error opening data structure file')
  else
     call log_error(m_name,s_name,2,log_info,'data structure file opened')
  end if

  !! header information
  write(nout,'(a)') trim(file%geoqout)
  write(nout,'(a)') trim(timestamp%long)
  write(nout,'(" vtk_input_file = ",a,/)') trim(file%vtkdata)
  write(nout,'(" eqdsk_input_file = ",a,/)') trim(file%eqdsk)
  write(nout,'(" equil_input_file = ",a,/)') trim(file%equil)

end  subroutine boutfile_init
!---------------------------------------------------------------------
!> write output from beq
subroutine boutfile_write(beq,timestamp)

  !! arguments
  type(beq_t), intent(in) :: beq !< beq object data structure
  type(date_time_t), intent(in) :: timestamp !< timestamp of run

  ! local variables

  !!write
  fld_specn: select case (beq%n%fldspec)
  case(1)
     call beq_writepart(beq,nout)
  case default
     call beq_writeplus(beq,nout)
  end select fld_specn
  ! next is for debugging
  !dbg     call beq_write(beq,nout)

end  subroutine boutfile_write
!---------------------------------------------------------------------
!> get unit number for output
subroutine boutfile_getunit(kunit)

  !! arguments
  integer, intent(out) :: kunit    !< log unit number

  kunit=nout

end subroutine boutfile_getunit
!---------------------------------------------------------------------
!> close output files
subroutine boutfile_close

  close(unit=nout)

end  subroutine boutfile_close

end module boutfile_m
