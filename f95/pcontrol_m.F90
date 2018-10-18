module pcontrol_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use misc_m
  use position_h
  use position_m
  use spl2d_m
  use spl3d_m
  use powcal_h
  use powres_h
  use powcal_m
  use odes_h
  use odes_m
!  use termplane_h
  use pcontrol_h

  implicit none
  private

! public subroutines
  public :: &
 &pcontrol_init, &
 &pcontrol_read, &
 &pcontrol_fix  !< create consistent controls


! public types

! private variables
  character(*), parameter :: m_name='pcontrol_m' !< module name
  integer  :: status   !< error status
  integer  :: nin      !< control file unit number
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  character(len=80) :: root !< file root

  contains
!---------------------------------------------------------------------
!> open input control data file
subroutine pcontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='pcontrol_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  character(len=80) :: controlfile !< control file name


  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

  !! open file
  controlfile=trim(fileroot)//".ctl"
  root=fileroot
  call log_value("Control data file",trim(controlfile))
  call misc_getfileunit(nin)
  open(unit=nin,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open control data file')
     stop
  end if

end  subroutine pcontrol_init
!---------------------------------------------------------------------
!> read data for this run
subroutine pcontrol_read(file,numerics,onumerics,edgprof,plot)

  !! arguments
  type(pfiles_t), intent(out) :: file !< file names
  type(pnumerics_t), intent(out) :: numerics !< input numerical parameters
  type(onumerics_t), intent(out) :: onumerics !< input ODE numerical parameters
  type(edgprof_t), intent(out) :: edgprof !< input ODE numerical parameters
  type(pplots_t), intent(out) :: plot !< vtk plot selectors

  !!local
  character(*), parameter :: s_name='pcontrol_read' !< subroutine name
  character(len=80) :: vtk_input_file  !< vtk format geometry data input file
  character(len=80) :: vtkres_input_file  !< vtk format `results' geometry data input file
  character(len=80) :: hds_input_file  !< vtk format HDS data input file
  character(len=80) :: geoq_input_file  !< geoq output format data input file
  character(len=80) :: lau_input_file  !< launch data input file
  logical :: filefound !< true of file exists
  integer(ki4) :: max_number_of_tracks !< maximum number of tracks
  logical :: plot_cartv !< DUPLICATE  vtk plot selector
  logical :: plot_powstatx !< vtk plot selector
  logical :: plot_allcartv !< DUPLICATE  vtk plot selector
  logical :: plot_powx !< vtk plot selector
  logical :: plot_wall !< vtk plot selector
  logical :: plot_ptzv !< INERT vtk plot selector
  logical :: plot_flincart !< DUPLICATE  vtk plot selector
  logical :: plot_flinx !< vtk plot selector
  logical :: plot_flinends !< gnu plot selector
  logical :: plot_flinptz !< DUPLICATE  vtk plot selector
  logical :: plot_flinm !< vtk plot selector
  logical :: plot_allptzv !< INERT vtk plot selector
  logical :: plot_geofld !< INERT vtk plot selector
  logical :: plot_geofld_quantised !< INERT vtk plot selector
  logical :: plot_gnu !< gnu plot selector

  !! file names
  namelist /inputfiles/ &
 &vtk_input_file, &
 &vtkres_input_file, &
 &hds_input_file, &
 &geoq_input_file, &
 &lau_input_file

  !! misc parameters
  namelist /miscparameters/ &
 &max_number_of_tracks

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_cartv, &
 &plot_powstatx, &
 &plot_allcartv, &
 &plot_powx, &
 &plot_wall, &
 &plot_ptzv, &
 &plot_flincart, &
 &plot_flinx, &
 &plot_flinptz, &
 &plot_flinm, &
 &plot_allptzv, &
 &plot_geofld, &
 &plot_geofld_quantised, &
 &plot_flinends, &
 &plot_gnu

  !! read input file names
  vtk_input_file='null'
  vtkres_input_file='null'
  hds_input_file='null'
  geoq_input_file='null'
  lau_input_file='null'
  read(nin,nml=inputfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_getunit(ilog)
     write(ilog,nml=inputfiles)
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
  end if

  file%vtkdata   = vtk_input_file
  file%vtkres   = vtkres_input_file
  file%hdsdata   = hds_input_file
  file%geoq   = geoq_input_file
  file%laudata   = lau_input_file

  !!check files exist

  call log_value("geometry data file, vtk_input_file",trim(file%vtkdata))
  if(file%vtkdata/='null') then
     inquire(file=vtk_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find geometry data file, ",a)',vtk_input_file
        call log_error(m_name,s_name,2,error_fatal,'geoq output data file not found')
     end if
  end if

  call log_value("results geometry data file, vtkres_input_file",trim(file%vtkres))
  inquire(file=vtkres_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find results geometry data file, ",a)',vtkres_input_file
     call log_error(m_name,s_name,3,error_fatal,'geoq output data file not found')
  end if

  if(file%hdsdata/='null') then
     call log_value("HDS data file, hds_input_file",trim(file%hdsdata))
     inquire(file=hds_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find hds data file, ",a)',hds_input_file
        call log_error(m_name,s_name,4,error_fatal,'HDS data file not found')
     end if
  end if

  call log_value("output from geoq, now geoq_input_file",trim(file%geoq))
  inquire(file=geoq_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find geoq output data file, ",a)',geoq_input_file
     call log_error(m_name,s_name,5,error_fatal,'geoq output data file not found')
  end if

  if(file%laudata/='null') then
     call log_value("launch data file, lau_input_file",trim(file%laudata))
     inquire(file=lau_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find lau data file, ",a)',lau_input_file
        call log_error(m_name,s_name,6,error_fatal,'launch data file not found')
     end if
  end if

  !! create output file names from root

  !! output file
  file%powcalout = trim(root)//"_powcal.out"
  inquire(file=file%powcalout,exist=filefound)
  call log_value("powcal output file",trim(file%powcalout))
  if(filefound) then
     !! error opening file
     print '("Fatal error: Output file ",a, " already exists")',trim(file%powcalout)
     print '("Remove it and restart run")'
     call log_error(m_name,s_name,6,error_fatal,'Output file already exists')
  end if

  !!vtk file roots
  file%cartv     =trim(root)//"_cartv"
  file%powstatx     =trim(root)//"_powstatx"
  file%allcartv     =trim(root)//"_allcartv"
  file%powx     =trim(root)//"_powx"
  file%wall     =trim(root)//"_wall"
  file%ptzv    =trim(root)//"_ptzv"
  file%flincart    =trim(root)//"_flincart"
  file%flinx    =trim(root)//"_flinx"
  file%flinptz    =trim(root)//"_flinptz"
  file%flinm    =trim(root)//"_flinm"
  file%allptzv    =trim(root)//"_allptzv"
  file%geofld =trim(root)//"_geofld"
  file%geofldq    =trim(root)//"_geofldq"
  !!gnu file roots
  file%gnu     =trim(root)//"_gnu"
  file%flinends    =trim(root)//"_flinends"

  !! set default misc parameters
  max_number_of_tracks = 100

  !!read misc parameters
  read(nin,nml=miscparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading misc parameters")'
     call log_getunit(ilog)
     write(ilog,nml=miscparameters)
     call log_error(m_name,s_name,10,error_fatal,'Error reading misc parameters')
  end if

  !! check for valid data
  if(max_number_of_tracks<0) &
 &call log_error(m_name,s_name,11,error_fatal,'max_number_of_tracks must be >=0')


  numerics%mtrack = max_number_of_tracks

  !! set default plot selections
  plot_cartv = .false.
  plot_powstatx = .false.
  plot_allcartv = .false.
  plot_powx = .false.
  plot_wall = .false.
  plot_ptzv = .false.
  plot_flincart = .false.
  plot_flinx = .false.
  plot_flinptz = .false.
  plot_flinm = .false.
  plot_allptzv = .false.
  plot_geofld = .false.
  plot_geofld_quantised = .false.
  plot_gnu = .false.
  plot_flinends = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     print '("Fatal error reading plot selections")'
     call log_getunit(ilog)
     write(ilog,nml=plotselections)
     call log_error(m_name,s_name,20,error_fatal,'Error reading plot selections')
  end if

  !! store values
  plot%powstatx = plot_cartv
  plot%powstatx   = plot_powstatx
  plot%powx = plot_allcartv
  plot%powx   = plot_powx
  plot%wall   = plot_wall
  plot%ptzv    = plot_ptzv
  plot%flinx    = plot_flinx
  plot%flinm    = plot_flinm
  plot%allptzv     = plot_allptzv
  plot%geofld    = plot_geofld
  plot%geofldq     = plot_geofld_quantised
  plot%gnu     = plot_gnu
  plot%flinends    = plot_flinends
  if (plot_flincart) then
     plot%flinx = plot_flincart
     call log_error(m_name,s_name,30,error_warning,'Obsolete file handling feature activated')
  end if
  if (plot_flinptz) then
     plot%flinm = plot_flinptz
     call log_error(m_name,s_name,31,error_warning,'Obsolete file handling feature activated')
  end if

  call powcal_readcon(numerics,nin)

  call edgprof_readcon(edgprof,numerics,nin)

  call odes_readcon(onumerics,nin)

end  subroutine pcontrol_read
!---------------------------------------------------------------------
!> create consistent controls
subroutine pcontrol_fix(numerics,onumerics,ifldspec)

  !! arguments
  type(pnumerics_t), intent(inout) :: numerics !< input numerical parameters
  type(onumerics_t), intent(inout) :: onumerics !< input ODE numerical parameters
  integer(ki4), intent(in) :: ifldspec !< field specification

  !!local
  character(*), parameter :: s_name='pcontrol_fix' !< subroutine name

  calcn_type: select case (numerics%caltype)
  case('afws','local')
     numerics%nanalau=0
     numerics%nlaupow=1
     numerics%nshapow=0
     onumerics%nstartcon=0
     onumerics%ntermcon=0
  case('msus','global')
     numerics%nanalau=0
     numerics%nlaupow=1
     numerics%nshapow=0
     onumerics%nstartcon=1
     onumerics%ntermcon=sign(1,onumerics%ntermcon)
  case default ! 'msum','middle'
  end select calcn_type

  fldspec_type: select case (ifldspec)
  case(1)
     onumerics%npdim=1
  case(2)
     onumerics%npdim=3
  case(3)
     onumerics%npdim=2
  end select fldspec_type

end  subroutine pcontrol_fix

end module pcontrol_m
