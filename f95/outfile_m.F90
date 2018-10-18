module outfile_m

  use const_kind_m
  use log_m
  use misc_m
  use control_h
  use date_time_m
  use geobj_m
  use ls_m
  use geobjlist_h
  use btree_m

  implicit none
  private

! public subroutines
  public :: &
 &outfile_init, &
 &outfile_dinit, &
 &outfile_linit, &
 &outfile_minit, &
 &outfile_write, &
 &outfile_dwrite, &
 &outfile_mwrite, &
 &outfile_close
! public types
  integer :: nout !< output file unit

! private variables
  character(*), parameter :: m_name='outfile_m' !< module name
  integer(ki4) :: i!< loop counters
  integer(ki4) :: j!< loop counters
  integer(ki4) :: k!< loop counters
  integer(ki4) :: nc !< loop counters
  integer(ki4) :: iemptycell  !< local variable
  integer(ki4) :: iin  !< local variable
  integer(ki4) :: imaxincell  !< local variable
  integer(ki4) :: imaxincellu  !< local variable
  integer(ki4) :: itotal !< number of output lines in block
  integer(ki4) :: itodu  !< local variable

  contains


!---------------------------------------------------------------------
!! open output files

subroutine outfile_init(file,timestamp)

  !! argument
  type(files_t), intent(in) :: file !< file names
  type(date_time_t), intent(in) :: timestamp !< timestamp of run


  !! local
  character(*), parameter :: s_name='outfile_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  !! get file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nout=i exit end if end do

  call misc_getfileunit(nout)
  open(unit=nout,file=trim(file%hdsgenout),status='new')

  !! header information
  write(nout,'(a)') trim(file%hdsgenout)
  write(nout,'(a)') trim(timestamp%long)
  write(nout,'(" vtk_input_file = ",a,/)') trim(file%vtkdata)

end  subroutine outfile_init

subroutine outfile_dinit(file,timestamp)

  !! argument
  type(files_t), intent(in) :: file !< file names
  type(date_time_t), intent(in) :: timestamp !< timestamp of run


  !! local
  character(*), parameter :: s_name='outfile_dinit' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  !! get file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nout=i exit end if end do

  call misc_getfileunit(nout)
  open(unit=nout,file=trim(file%dengenout),status='new')

  !! header information
  write(nout,'(a)') trim(file%dengenout)
  write(nout,'(a)') trim(timestamp%long)
  write(nout,'(" vtk_input_file = ",a,/)') trim(file%vtkdata)
  write(nout,'(" hds_input_file = ",a,/)') trim(file%hdsdata)
  write(nout,'(" query_input_file = ",a,/)') trim(file%qrydata)

end  subroutine outfile_dinit

subroutine outfile_linit(file,timestamp,kout)

  !! argument
  type(files_t), intent(in) :: file !< file names
  type(date_time_t), intent(in) :: timestamp !< timestamp of run
  integer, intent(out) :: kout   !< output channel


  !! local
  character(*), parameter :: s_name='outfile_linit' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  !! get file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nout=i exit end if end do
  kout=nout

  call misc_getfileunit(nout)
  open(unit=nout,file=trim(file%legtfmout),status='new')

  !! header information
  write(nout,'(a)') trim(file%legtfmout)
  write(nout,'(a)') trim(timestamp%long)
  write(nout,'(" vtk_input_file = ",a,/)') trim(file%vtkdata)
  write(nout,'(" query_input_file = ",a,/)') trim(file%qrydata)

end  subroutine outfile_linit


subroutine outfile_minit(file,timestamp)

  !! argument
  type(files_t), intent(in) :: file !< file names
  type(date_time_t), intent(in) :: timestamp !< timestamp of run


  !! local
  character(*), parameter :: s_name='outfile_minit' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  !! get file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nout=i exit end if end do

  call misc_getfileunit(nout)
  open(unit=nout,file=trim(file%moveout),status='new')

  !! header information
  write(nout,'(a)') trim(file%moveout)
  write(nout,'(a)') trim(timestamp%long)
  write(nout,'(" vtk_input_file = ",a,/)') trim(file%vtkdata)
  write(nout,'(" hds_input_file = ",a,/)') trim(file%hdsdata)
  write(nout,'(" query_input_file = ",a,/)') trim(file%qrydata)

end  subroutine outfile_minit

!---------------------------------------------------------------------
!! write sorted geobj data

subroutine outfile_write(geobjl,btree)

  !! arguments
  type(geobjlist_t), intent(in) :: geobjl !< point data
  type(btree_t), intent(inout) :: btree !< btree data

  ! local variables
  integer(ki4) :: itotal  !< local variable
  integer(ki4) :: iadr  !< local variable

  !!dimensions
  !! count empty cells
  iemptycell=0
  imaxincell=0
  imaxincellu=0
  itotal=0
  itodu=0
  do i=1,btree%nt
     j=btree%pter(3,i)
     if (j==-2) then
        ! empty
        iemptycell=iemptycell+1
     else if (j==-1) then
        ! fully processed
        iadr=btree%pter(2,i)
        iin=btree%objectls%list(iadr,2)
        imaxincell=max(iin,imaxincell)
        itotal=itotal+iin
     else if (j==0) then
        ! still to be processed
        iadr=btree%pter(2,i)
        iin=btree%objectls%list(iadr,2)
        imaxincellu=max(iin,imaxincellu)
        itotal=itotal+iin
        itodu=itodu+1
     end if
  end do

  !!write
  write(nout,'("statistics")')
  write(nout,'("  total number of tree leaves   ",i10)') btree%nt
  write(nout,'("  actual depth of tree   ",i10)') btree%ndepth
  write(nout,'("  number of empty cells     ",i10)') iemptycell
  write(nout,'("  number of cells unfinished     ",i10)') itodu
  write(nout,'("  max number in finished cell ",i10)') imaxincell
  write(nout,'("  max number in unfinished cell ",i10)') imaxincellu
  write(nout,'("  total number of objects in cells ",i10)')itotal
  write(nout,'("  number of unassigned objects ",i10)') geobjl%ngunassigned
  write(nout,'("end_statistics",/)')


  !!geobj listed by cell
  write(nout,'("bycell")')
  !!format is
  !!    cell <no> <index>
  !!    p1 p2 ...
  do i=1,btree%nt
     j=btree%pter(3,i)
     if (j>0.OR.j==-2) then
        cycle
     else
        iadr=btree%pter(2,i)
        write(nout,'("cell ",2i8)') iadr, i
        call ls_write1(btree%objectls,iadr,nout)
     end if
  end do
  write(nout,'("bycell_end")')

end  subroutine outfile_write

!---------------------------------------------------------------------
!! write output from dengen

subroutine outfile_dwrite(timestamp)

  !! arguments
  type(date_time_t), intent(in) :: timestamp !< timestamp of run

  ! local variables

  !!write
  write(nout,'("input HDS file date stamp")')
  write(nout,'(a)') trim(timestamp%long)

end  subroutine outfile_dwrite


!---------------------------------------------------------------------
!! write output from move

subroutine outfile_mwrite(timestamp)

  !! arguments
  type(date_time_t), intent(in) :: timestamp !< timestamp of run

  ! local variables

  !!write
  write(nout,'("input HDS file date stamp")')
  write(nout,'(a)') trim(timestamp%long)

end  subroutine outfile_mwrite


!---------------------------------------------------------------------
!! close output files

subroutine outfile_close

  close(unit=nout)

end  subroutine outfile_close

end module outfile_m
