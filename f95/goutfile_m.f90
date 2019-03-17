module goutfile_m

  use const_kind_m
  use log_m
  use misc_m
  use const_numphys_h
  use dcontrol_h
  use bcontrol_m
  use date_time_m
  use position_h
  use fmesh_h
  use skyl_h
  use beq_h
  use geoq_h
  use geobj_m
  use spl2d_m
  use geobjlist_h
  use geobjlist_m
  use beq_m
  use vfile_m
  use dcontrol_m
  use skyl_m

  implicit none
  private

! public subroutines
  public :: &
 &goutfile_init, & !< find unit, open file and write header
 &goutfile_write, & !< write data structure
 &goutfile_getunit,  & !< get unit number
 &goutfile_close

! private variables
  character(*), parameter :: m_name='goutfile_m' !< module name
  integer, save :: nout=-1 !< output file unit
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: idum !< dummy integer
  integer   :: status   !< error status
  integer(ki4) :: ifldspec !< field spec
  integer(ki4) :: iextra !< extra info extracted

  contains
!---------------------------------------------------------------------
!> open output file
subroutine goutfile_init(file,timestamp)

  !! argument
  type(bfiles_t), intent(in) :: file !< file names
  type(date_time_t), intent(in) :: timestamp !< timestamp of run


  !! local
  character(*), parameter :: s_name='goutfile_init' !< subroutine name
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

end  subroutine goutfile_init
!---------------------------------------------------------------------
!> write output from geoq
subroutine goutfile_write(geoq,timestamp)

  !! arguments
  type(geoq_t), intent(inout) :: geoq !< geoq object data structure
  type(date_time_t), intent(in) :: timestamp !< timestamp of run

  ! local variables

  !ifldspec=geoq%beq%n%fldspec
  !iextra=ifldspec/10
  !ifldspec=ifldspec-10*iextra
  ifldspec=mod(geoq%beq%n%fldspec,10)
  !! mark output file (should not have both duct and skylight)
  if (geoq%beq%n%duct) geoq%beq%n%fldspec=ifldspec+10
  if (geoq%beq%n%skylpsi) geoq%beq%n%fldspec=ifldspec+20
  if (geoq%beq%n%objadd(GEOBJ_SKYLIT)>0) geoq%beq%n%fldspec=ifldspec+30
  if (geoq%beq%n%skylpsi.AND.geoq%beq%n%objadd(GEOBJ_SKYLIT)>0) geoq%beq%n%fldspec=ifldspec+40

  !!write
  fld_specn: select case (ifldspec)
  case(1)
     ! writing part, only output fldspec in range 1 to 3
     geoq%beq%n%fldspec=ifldspec
     call beq_writepart(geoq%beq,nout)
  case default
     call beq_writeplus(geoq%beq,nout)
  end select fld_specn

  if (geoq%beq%n%skyl) call skyl_write(geoq%skyl,nout)

end  subroutine goutfile_write
!---------------------------------------------------------------------
!> get unit number for output
subroutine goutfile_getunit(kunit)

  !! arguments
  integer, intent(out) :: kunit    !< log unit number

  kunit=nout

end subroutine goutfile_getunit
!---------------------------------------------------------------------
!> close output files
subroutine goutfile_close

  close(unit=nout)

end  subroutine goutfile_close

end module goutfile_m
