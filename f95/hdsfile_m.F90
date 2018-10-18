module hdsfile_m

  use const_kind_m
  use date_time_m
  use log_m
  use misc_m
  use control_h
  use btree_m
  use geobjlist_h

  implicit none
  private

!! public subroutines
  public :: &
 &hdsfile_init, & !< initialise HDS file for writing
 &hdsfile_dinit, & !< initialise HDS file for reading
 &hdsfile_write, & !< write HDS file
 &hdsfile_read, & !< read HDS file
 &hdsfile_close !> close HDS file

!! public types

!! private variables
  character(*), parameter :: m_name='hdsfile_m' !< module name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer :: nouth  !< file unit for output
  integer :: status  !< status flag

  contains

!> write HDS file header
subroutine hdsfile_init(fouth,timestamp,kouth)

  !! arguments
  character(len=*), intent(in) :: fouth !< file name root
  type(date_time_t), intent(in) :: timestamp !< timestamp of run
  integer, intent(inout) :: kouth   !< unit number

  !! local
  character(*), parameter :: s_name='hdsfile_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  !! open file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kouth=i exit end if end do

  call misc_getfileunit(kouth)
  open(unit=kouth,file=trim(fouth)//'.hds')

  !! write HDS header
  write(kouth,'(A)') trim(timestamp%long)
  write(kouth,'(''           '')')

  !! save unit

  nouth=kouth

end subroutine hdsfile_init
!---------------------------------------------------------------------
!> read HDS file header
subroutine hdsfile_dinit(infile,timestamp)

  !! arguments
  character(len=*), intent(in) :: infile !< file name root
  type(date_time_t), intent(out) :: timestamp !< timestamp of run

  !! local
  character(*), parameter :: s_name='hdsfile_dinit' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  !! open file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nouth=i exit end if end do

  call misc_getfileunit(nouth)
  open(unit=nouth,file=infile)

  !! read HDS header and rewind
  read(nouth,'(A)') timestamp%long
  !     read(nouth,'(''           '')')

  rewind(nouth)

end subroutine hdsfile_dinit
!---------------------------------------------------------------------
!> write HDS file
subroutine hdsfile_write(geobjl,numerics,btree)

  !! arguments
  type(geobjlist_t), intent(in) :: geobjl !< geometry data
  type(numerics_t), intent(inout) :: numerics !< geobjlist controls
  type(btree_t), intent(in) :: btree !< btree data

  !! local
  character(*), parameter :: s_name='hdsfile_write' !< subroutine name

  call btree_write(btree,numerics,nouth)

  !! write transform arrays

  ! write tfmdata label
  write(nouth,'(a)',iostat=status) 'TFMDA'
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error writing  tfmdata label')
  end if

  !     write(nouth,*,iostat=status) geobjl%tfmdata
  write(nouth,*,iostat=status) geobjl%tfmdata%scale
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error writing  tfmdata ')
  end if
  write(nouth,*,iostat=status)  geobjl%tfmdata%offset
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error writing  tfmdata ')
  end if
  write(nouth,*,iostat=status)  geobjl%tfmdata%matrix
  if(status/=0) then
     call log_error(m_name,s_name,5,error_fatal,'Error writing  tfmdata ')
  end if
  write(nouth,*,iostat=status)  geobjl%tfmdata%ntfm
  if(status/=0) then
     call log_error(m_name,s_name,6,error_fatal,'Error writing  tfmdata ')
  end if


  ! write quantfm label
  write(nouth,'(a)',iostat=status) 'QTFMD'
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error writing  quantfm label')
  end if

  !     write(nouth,*,iostat=status) geobjl%quantfm
  write(nouth,*,iostat=status) geobjl%quantfm%hmin
  if(status/=0) then
     call log_error(m_name,s_name,10,error_fatal,'Error writing  quantfm ')
  end if
  write(nouth,*,iostat=status)  geobjl%quantfm%rhmin
  if(status/=0) then
     call log_error(m_name,s_name,11,error_fatal,'Error writing  quantfm ')
  end if
  write(nouth,*,iostat=status)  geobjl%quantfm%offvec
  if(status/=0) then
     call log_error(m_name,s_name,12,error_fatal,'Error writing  quantfm ')
  end if
  write(nouth,*,iostat=status)  geobjl%quantfm%nqtfm
  if(status/=0) then
     call log_error(m_name,s_name,13,error_fatal,'Error writing  quantfm ')
  end if

  write(nouth,'(a)',iostat=status) 'NUMPAR'
  call log_write_check(m_name,s_name,20,status)

  write(nouth,*,iostat=status)  numerics%geomtype
  call log_write_check(m_name,s_name,21,status)
  write(nouth,*,iostat=status)  numerics%maxtolerance
  call log_write_check(m_name,s_name,22,status)
  write(nouth,*,iostat=status)  numerics%mintolerance
  call log_write_check(m_name,s_name,23,status)
  write(nouth,*,iostat=status)  numerics%mingeobjinbin
  call log_write_check(m_name,s_name,24,status)
  write(nouth,*,iostat=status)  numerics%nquante
  call log_write_check(m_name,s_name,25,status)
  write(nouth,*,iostat=status)  numerics%position_coord_tfm
  call log_write_check(m_name,s_name,26,status)

  write(nouth,'(a)',iostat=status) 'ADDCU'
  call log_write_check(m_name,s_name,27,status)

  write(nouth,*,iostat=status)  numerics%nbdcub
  call log_write_check(m_name,s_name,30,status)
  write(nouth,*,iostat=status)  numerics%dilen
  call log_write_check(m_name,s_name,31,status)
  write(nouth,*,iostat=status)  numerics%dolen
  call log_write_check(m_name,s_name,32,status)
  write(nouth,*,iostat=status)  numerics%cornflag
  call log_write_check(m_name,s_name,33,status)
  write(nouth,*,iostat=status)  numerics%lowcorner
  call log_write_check(m_name,s_name,34,status)
  write(nouth,*,iostat=status)  numerics%upcorner
  call log_write_check(m_name,s_name,35,status)

  write(nouth,'(a)',iostat=status) 'BOUBOX'
  call log_write_check(m_name,s_name,37,status)

  write(nouth,*,iostat=status)  geobjl%coordbb
  call log_write_check(m_name,s_name,40,status)
  write(nouth,*,iostat=status)  geobjl%binbb
  call log_write_check(m_name,s_name,41,status)

end subroutine hdsfile_write
!---------------------------------------------------------------------
!> read HDS file
subroutine hdsfile_read(geobjl,numerics,btree)

  !! arguments
  type(geobjlist_t), intent(inout) :: geobjl !< point data
  type(numerics_t), intent(inout) :: numerics !< geobjlist controls
  type(btree_t), intent(inout) :: btree !< btree data

  !! local
  character(*), parameter :: s_name='hdsfile_read' !< subroutine name
  character(len=5) :: ichar !< local variable

  ! skip header
  do j=1,2
     read(nouth,'(a)',iostat=status) ichar
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Error reading HDS file header')
     end if
  end do

  do
     read(nouth,'(a)',iostat=status) ichar
     !write(*,*) 'ichar=',ichar
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,2,error_fatal,'Error reading HDS file')
     end if
     if (ichar=='BTREE'.OR.ichar=='EXTEN'.OR.ichar=='LSLAB') then

        call btree_read(btree,numerics,ichar,nouth)

     else if (ichar=='TFMDA') then
        ! read tfmdata label
        read(nouth,*,iostat=status) numerics%position_coord_tfm%scale
        if(status/=0) then
           call log_error(m_name,s_name,3,error_fatal,'Error reading  tfmdata ')
        end if
        read(nouth,*,iostat=status)  numerics%position_coord_tfm%offset
        if(status/=0) then
           call log_error(m_name,s_name,4,error_fatal,'Error reading  tfmdata ')
        end if
        read(nouth,*,iostat=status)  numerics%position_coord_tfm%matrix
        if(status/=0) then
           call log_error(m_name,s_name,5,error_fatal,'Error reading  tfmdata ')
        end if
        read(nouth,*,iostat=status)  numerics%position_coord_tfm%ntfm
        if(status/=0) then
           call log_error(m_name,s_name,6,error_fatal,'Error reading  tfmdata ')
        end if
        geobjl%tfmdata=numerics%position_coord_tfm

        ! read quantfm label
     else if (ichar=='QTFMD') then
        read(nouth,*,iostat=status) numerics%geobj_coord_tfm%hmin
        if(status/=0) then
           call log_error(m_name,s_name,10,error_fatal,'Error reading quantfm ')
        end if
        read(nouth,*,iostat=status)  numerics%geobj_coord_tfm%rhmin
        if(status/=0) then
           call log_error(m_name,s_name,11,error_fatal,'Error reading quantfm ')
        end if
        read(nouth,*,iostat=status)  numerics%geobj_coord_tfm%offvec
        if(status/=0) then
           call log_error(m_name,s_name,12,error_fatal,'Error reading quantfm ')
        end if
        read(nouth,*,iostat=status)  numerics%geobj_coord_tfm%nqtfm
        if(status/=0) then
           call log_error(m_name,s_name,13,error_fatal,'Error reading quantfm ')
        end if
        geobjl%quantfm=numerics%geobj_coord_tfm

        ! read numerical parameters label
     else if (ichar=='NUMPA') then

        read(nouth,*,iostat=status)  numerics%geomtype
        if(status/=0) then
           call log_error(m_name,s_name,21,error_fatal,'Error reading numpar ')
        end if
        read(nouth,*,iostat=status)  numerics%maxtolerance
        if(status/=0) then
           call log_error(m_name,s_name,22,error_fatal,'Error reading numpar ')
        end if
        read(nouth,*,iostat=status)  numerics%mintolerance
        if(status/=0) then
           call log_error(m_name,s_name,23,error_fatal,'Error reading numpar ')
        end if
        read(nouth,*,iostat=status)  numerics%mingeobjinbin
        if(status/=0) then
           call log_error(m_name,s_name,24,error_fatal,'Error reading numpar ')
        end if
        read(nouth,*,iostat=status)  numerics%nquante
        if(status/=0) then
           call log_error(m_name,s_name,25,error_fatal,'Error reading numpar ')
        end if
        read(nouth,*,iostat=status)  numerics%position_coord_tfm
        if(status/=0) then
           call log_error(m_name,s_name,26,error_fatal,'Error reading numpar ')
        end if

     else if (ichar=='ADDCU') then
        read(nouth,*,iostat=status)  numerics%nbdcub
        if(status/=0) then
           call log_error(m_name,s_name,30,error_fatal,'Error reading addcubes ')
        end if
        read(nouth,*,iostat=status)  numerics%dilen
        if(status/=0) then
           call log_error(m_name,s_name,31,error_fatal,'Error reading addcubes ')
        end if
        read(nouth,*,iostat=status)  numerics%dolen
        if(status/=0) then
           call log_error(m_name,s_name,32,error_fatal,'Error reading addcubes ')
        end if
        read(nouth,*,iostat=status)  numerics%cornflag
        if(status/=0) then
           call log_error(m_name,s_name,33,error_fatal,'Error reading addcubes ')
        end if
        read(nouth,*,iostat=status)  numerics%lowcorner
        if(status/=0) then
           call log_error(m_name,s_name,34,error_fatal,'Error reading addcubes ')
        end if
        read(nouth,*,iostat=status)  numerics%upcorner
        if(status/=0) then
           call log_error(m_name,s_name,35,error_fatal,'Error reading addcubes ')
        end if

     else if (ichar=='BOUBO') then
        read(nouth,*,iostat=status)  geobjl%coordbb
        !write(*,*) 'geobjl%coordbb',geobjl%coordbb
        if(status/=0) then
           call log_error(m_name,s_name,40,error_fatal,'Error reading coordbb ')
        end if
        read(nouth,*,iostat=status)  geobjl%binbb
        if(status/=0) then
           call log_error(m_name,s_name,41,error_fatal,'Error reading binbb ')
        end if
     end if

  end do

end subroutine hdsfile_read
!---------------------------------------------------------------------
!> close HDS file on unit nouth
subroutine hdsfile_close

  close(nouth)

end subroutine hdsfile_close

end module hdsfile_m
