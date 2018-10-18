module ccontrol_m

  use const_kind_m
  use log_m
  use misc_m

  implicit none

! public subroutines
  public :: &
 &ccontrol_init, & !< open input ctlin format data file and output ctl format file
 &ccontrol_readctlin, & !< read ctlin format data for this code
 &ccontrol_writectl !< write  ctl data for this code

! private variables
  character(*), parameter :: m_name='ccontrol_m' !< module name
  integer(ki4), parameter :: maximum_number_of_files=100 !< maximum number of files allowed
  integer(ki4), parameter :: maximum_number_of_panels=100 !< maximum number of panels allowed

  integer(ki4), parameter :: maximum_number_of_lines=100 !< maximum number of lines in files allowed
  integer(ki4), parameter :: maximum_number_of_namelists=20 !< maximum number of namelists in files allowed
  integer(ki4), parameter :: maximum_number_of_bodies_in_file=100 !< maximum number of bodies in files allowed
  character(len=3), dimension(maximum_number_of_lines), save :: filekey  !< labels data input line
  character(len=132), dimension(maximum_number_of_lines), save :: filecon !< contents of data input line
  integer(ki4) :: nincon !< number of entries in filecon
  integer(ki4), dimension(maximum_number_of_lines), save :: nmlno !< number of namelist
  character(len=40), dimension(maximum_number_of_namelists), save :: namel !< name of namelist
  integer(ki4), save :: inpfile      !< number of input files

  integer  :: status   !< error status
  integer, save  :: nin      !< input file unit number
  integer, save  :: nout      !< output file unit number
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: islen !< string length
  integer(ki4), dimension(:), allocatable :: iwork1 !< work array
  integer(ki4), dimension(:), allocatable :: iwork2 !< work array
  character(len=132) :: buff !< workspace
  character(len=80) :: root !< file root

  contains
!---------------------------------------------------------------------
!> open input ctlin format data file and output ctl format file
subroutine ccontrol_init(fileroot,code)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  character(len=80), intent(out) :: code   !< code input file intended for
  !! local
  character(*), parameter :: s_name='ccontrol_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  character(len=80) :: inputfile !< input file name
  character(len=80) :: outputfile !< output file name
  character(len=3) :: key !< key
  character(len=80) :: icode !< code input file intended for
  logical :: lok !< local variable

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

  !! open file
  inputfile=trim(fileroot)//".ctlin"
  root=fileroot
  call log_value("Control data input file",trim(inputfile))
  call misc_getfileunit(nin)
  open(unit=nin,file=inputfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open input file, ",a)',inputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open input data file')
     stop
  end if
  !---------------------------------------------------------------------
  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nout=i exit end if end do

  !! open file
  outputfile=trim(fileroot)//".ctl"
  call log_value("Control data output file",trim(outputfile))
  call misc_getfileunit(nout)
  open(unit=nout,file=outputfile,iostat=status,recl=1610,delim='apostrophe')
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,2,error_fatal,'Cannot open output data file')
     stop
  end if

  read(nin,'(A3,A)',iostat=status) key, icode
  call log_read_check(m_name,s_name,3,status)
  if (key/='STA') then
     call log_error(m_name,s_name,2,error_fatal,'Input data must start with STA <progname>')
  end if
  code=trim(adjustl(icode))
  rewind(nin)

end  subroutine ccontrol_init
!---------------------------------------------------------------------
!> read ctlin format data for this code
subroutine ccontrol_readctlin(code)

  !! arguments
  character(len=80), intent(inout) :: code   !< code name

  !!local
  character(*), parameter :: s_name='ccontrol_readctlin' !< subroutine name
  character(len=3) :: key !< key
  character(len=80) :: icode !< code input file intended for
  integer(ki4) :: isub !< subtract
  logical :: lok !< local variable

  !! namelist names (not used at present)
  namel(1)='miscparameters'
  namel(2)='vtkfiles'
  namel(3)='panelarrayparameters'
  namel(4)='positionparameters'

  read(nin,'(A3,A)',iostat=status) key, icode
  call log_read_check(m_name,s_name,1,status)
  lok=(key=='STA'.AND.trim(adjustl(icode))==trim(code))
  if (.NOT.lok) then
     call log_error(m_name,s_name,2,error_fatal,'Input data must start with STA <matching progname>')
  end if

  i=2
  do
     read(nin,'(A3,A)',iostat=status) filekey(1), filecon(1)

     if(status<0) then
        !!eof
        exit
     else if (status>0) then
        !! error
        call log_error(m_name,s_name,3,error_fatal,'Error reading input data')
     else if (filekey(1)=='STA') then
        icode=filecon(1)
        exit
     else
        isub=0
        filekey(i)=filekey(1)
        filecon(i)=filecon(1)
        look_for_match: do j=2,i-1
           if (filekey(j)==filekey(1)) then
              filecon(j)=filecon(1)
              isub=1
              exit
           end if
        end do  look_for_match
        i=i-isub+1
     end if

  end do

  nincon=i
  filekey(1)=key
  filecon(1)=code
  code=trim(adjustl(icode))

  close(unit=nin)

end subroutine ccontrol_readctlin
!---------------------------------------------------------------------
!> write  ctl data for this code
subroutine ccontrol_writectl

  !local
  character(*), parameter :: s_name='ccontrol_writectl' !< subroutine name
  integer(ki4) :: iscan !< position of blank
  character(len=80),dimension(:), allocatable:: vtk_input_file !< vtk format geometry data input file
  integer(ki4),dimension(:), allocatable:: number_of_copies !< number of copies to make of input file
  character(len=80),dimension(:), allocatable:: vtk_label  !< labels vtk format geometry data input file
  integer(ki4),dimension(:), allocatable:: vtk_label_number !< numeric label for input file (not used)
  character(len=20) :: angle_units  !< units, either radian(s) or degree(s)
  character(len=20) :: option !< for defining panels and their transforms
  character(len=80) :: vtk_output_file  !< vtk format geometry data output file
  integer(ki4) :: number_of_panels !< number of panels for which bodies defined
  integer(ki4) :: number_of_transforms !< number of panels for which transform defined
  integer(ki4) :: max_number_of_files !< number >= no. of files for which transform defined
  integer(ki4) :: max_number_of_panels !< number >= no. of panels for which transform defined
  integer(ki4) :: max_number_of_transforms !< number >= no. of panels for which transform defined
  integer(ki4), dimension(:), allocatable :: panel_transform !< number of transform to apply
  character(len=20), dimension(:), allocatable :: transform_id !< INACTIVE id of transform to apply
  integer(ki4), dimension(:), allocatable :: panel_bodies !< bodies defining the geometry
  real(kr4), dimension(3,3):: position_matrix  !< transformation matrix 3x3
  real(kr4), dimension(3):: position_scale  !< transformation scale 3-vector
  real(kr4), dimension(3):: position_offset  !< transformation offset 3-vector
  real(kr4) :: lenfac  !< convert to m if set
  integer(ki4):: position_transform  !< transformation type
  character(len=20) :: transform_desc !< describes transformation type
  integer(ki4):: ntfm  !< number of transformations

  !! misc parameters
  namelist /miscparameters/ &
 &option, &
 &angle_units, &
 &max_number_of_panels,max_number_of_transforms ,&
 &max_number_of_files !, &
  ! &number_of_panels,number_of_transforms

  !! file names
  namelist /vtkfiles/ &
 &vtk_input_file, number_of_copies!, &
  ! &vtk_label, vtk_label_number, &
  ! &vtk_output_file

  !! panelarray parameters
  namelist /panelarrayparameters/ &
 &panel_bodies,panel_transform !, &
  !&transform_id !omitted

  !! position parameters
  namelist /positionparameters/ &
 &position_scale , &
 &position_offset , &
 &position_transform !, &
  ! &position_matrix , &  ! omitted
  ! &transform_desc, & !omitted
  ! &transform_id !omitted

  ! determine useful globals and assign namelist numbers
  ntfm=0
  do j=2,nincon
     nml_key: select case (filekey(j))
     case('DEF')
        nmlno(j)=1
     case('FIL')
        nmlno(j)=2
        !count  file entries
        iscan=1
        inpfile=1
        buff=filecon(j)
        do k=1,maximum_number_of_files
           buff=adjustl(buff(iscan:))
           islen=len_trim(buff)
           iscan=scan(buff(:islen),' ')
           if (iscan>0.AND.iscan<islen) then
              inpfile=inpfile+1
           else
              exit
           end if
        end do
     case('TFM')
        ntfm=1
        nmlno(j)=4
     case('SCA','OFS')
        nmlno(j)=4
     end select nml_key
  end do

  !! allocate arrays
  allocate(vtk_input_file(inpfile),number_of_copies(inpfile), &
  panel_bodies(inpfile),panel_transform(inpfile), stat=status)
  call log_alloc_check(m_name,s_name,55,status)
  allocate(vtk_label(inpfile),vtk_label_number(inpfile), &
  transform_id(inpfile), stat=status)
  call log_alloc_check(m_name,s_name,56,status)


  lenfac=1.
  !! set default miscparameters
  max_number_of_files = inpfile !maximum_number_of_files
  max_number_of_panels = inpfile !maximum_number_of_panels
  max_number_of_transforms = 1
  angle_units = 'degree'
  option = 'panel'
  !! these not written out
  number_of_panels = 0
  number_of_transforms = 0
  do j=2,nincon
     if (nmlno(j)==1) then
        ! always DEF
        buff=filecon(j)
        if (scan(buff,'R')>0) then
           angle_units = 'radians'
        end if
        if (scan(buff,'M')>0) then
           lenfac=1000
        end if
        exit
     end if
  end do
  write(nout,nml=miscparameters, iostat=status)
  call log_write_check(m_name,s_name,10,status)

  !! file names
  vtk_input_file='null'
  number_of_copies=1
  !! these not written out
  vtk_output_file='null'
  vtk_label='null'
  vtk_label_number=0
  do j=2,nincon
     if (nmlno(j)==2) then
        ! always FIL
        buff=filecon(j)
        iscan=1
        do k=1, inpfile
           buff=adjustl(buff(iscan:))
           islen=len_trim(buff)
           iscan=scan(buff(:islen),' ')
           if (iscan>0) islen=iscan-1
           vtk_input_file(k)=trim(buff(:islen))//'.vtk'
        end do
        exit
     end if
  end do
  number_of_copies=1
  write(nout,nml=vtkfiles, iostat=status)
  call log_write_check(m_name,s_name,11,status)

  !! panelarray parameters
  do j=1,inpfile
     panel_bodies(j)=maximum_number_of_bodies_in_file*j+1
  end do
  if (inpfile==1) panel_bodies(1)=1
  panel_transform=ntfm
  write(nout,nml=panelarrayparameters, iostat=status)
  call log_write_check(m_name,s_name,12,status)

  !! set default position parameters
  ! default identity
  position_transform=0
  position_matrix(1,:)=(/1.,  0.,  0. /)
  position_matrix(2,:)=(/0.,  1.,  0. /)
  position_matrix(3,:)=(/0.,  0.,  1. /)
  position_scale=(/1.,  1.,  1. /)
  position_offset=(/0.,  0.,  0. /)

  do j=2,nincon
     if (nmlno(j)==4) then

        var_key: select case (filekey(j))
        case('TFM')
           read(filecon(j),*) ntfm
           position_transform=ntfm
        case('SCA')
           read(filecon(j),*) position_scale
        case('OFS')
           read(filecon(j),*) position_offset
        case('MAT')
           read(filecon(j),*) position_matrix
        end select var_key

     end if
  end do
  select  case (ntfm)
  case(2,3,22,42)
     position_offset=lenfac*position_offset
  end select
  write(nout,nml=positionparameters, iostat=status)
  call log_write_check(m_name,s_name,14,status)

  close(unit=nout)

end subroutine ccontrol_writectl

end module ccontrol_m
