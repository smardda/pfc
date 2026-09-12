module atcontrol_m

  use const_kind_m
  use log_m
  use misc_m
  use beq_h
  use atcontrol_h
  use antest_h
  use antest_m
  use antenna_h
  use antenna_m
  use radnsource_h
  use radnsource_m
  use ray_m
  use beq_m

  implicit none
  private

! public subroutines
  public :: &
 &atcontrol_init, & !< open input control data file
 &atcontrol_getunit, & !< get unit number for atcontrol reading
 &atcontrol_close, & !< close input control data file
 &atcontrol_read !< read data for this run

! private variables
  character(*), parameter :: m_name='atcontrol_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4), save  :: nin      !< control file unit number
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer :: ierror   !<  whether error is to be fatal or not
  character(len=80) :: root !< file root

  contains
!---------------------------------------------------------------------
!> open input control data file
subroutine atcontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='atcontrol_init' !< subroutine name
  logical :: unitused !< flag to test unit is available
  character(len=80) :: controlfile !< control file name


  !! get file unit
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        nin=i
        exit
     end if
  end do

  !! open file
  controlfile=trim(fileroot)//".ctl"
  root=fileroot
  call log_value("Control data file",trim(controlfile))
  open(unit=nin,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open control data file')
     stop
  end if


end  subroutine atcontrol_init
!---------------------------------------------------------------------
!> get unit number for atcontrol reading
subroutine atcontrol_getunit(kunit)

  !! arguments
  integer, intent(out) :: kunit    !< log unit number

  kunit=nin

end subroutine atcontrol_getunit
!---------------------------------------------------------------------
!> close input control data file
subroutine atcontrol_close

  !! local
  character(*), parameter :: s_name='atcontrol_close' !< subroutine name

  !! close file unit
  close(unit=nin)

end  subroutine atcontrol_close
!---------------------------------------------------------------------
!> read data for this run
subroutine atcontrol_read(file,param,atnumerics,antnumerics,radnumerics,raynumerics,numerics,plot)

  !! arguments
  type(atfiles_t), intent(out) :: file !< file names
  type(atparams_t), intent(out) :: param !< control parameters
  type(atnumerics_t), intent(out) :: atnumerics !< controls for antest
  type(antnumerics_t), intent(out) :: antnumerics !< controls for antenna variables
  type(radnumerics_t), intent(out) :: radnumerics !< radnsource parameters
  type(raynumerics_t), intent(out) :: raynumerics !< ray parameters
  type(bnumerics_t)  :: numerics  !< numerical control parameters for equilibrium B
  type(atplots_t), intent(out) :: plot !< plot selectors

  !!local
  character(*), parameter :: s_name='atcontrol_read' !< subroutine name
  logical :: filefound !< true if file exists


  character(len=80) :: prog_control  !< option control parameter
  character(len=80) :: vtk_input_file  !< wall geometry data input file
  character(len=80) :: eqdsk_input_file  !< EQDSK G format data input file
  character(len=80) :: equil_input_file  !< generic equilibrium input file
  character(len=80) :: antenna_input_file  !< launch data input file
  character(len=80) :: radiation_input_file  !< radiation data input file
  character(len=80) :: prog_input_file  !< antest data input file
  character(len=80) :: prog_output_file  !< antest data output file
  real(kr8) :: prog_realpar !< real control parameter
  integer(ki4) :: prog_intpar !< integer control parameter
  logical :: prog_logicpar !< logical control parameter
  character(len=80) :: icsuf  !< eqdsk file suffix

  logical :: plot_antestout !< antest output data selector
  logical :: plot_vtk !< vtk plot selector
  logical :: plot_gnu !< gnuplot plot selector
  logical :: plot_powx !< powx vtk plot selector

  !! file names
  namelist /progfiles/ &
 &vtk_input_file, &
 &equil_input_file, &
 &eqdsk_input_file, &
 &antenna_input_file, &
 &radiation_input_file, &
 &prog_input_file, &
 &prog_output_file

  !! misc parameters
  namelist /miscparameters/ &
 &prog_control,&
 &prog_realpar, prog_intpar, &
 &prog_logicpar

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_antestout, &
 &plot_vtk, &
 &plot_gnu, &
 &plot_powx

  !---------------------------------------------------------------------
  !! read input file names
  vtk_input_file='wall.vtk'
  equil_input_file='null'
  eqdsk_input_file='ndunr.eqdsk'
  antenna_input_file='ant0_vtk.vtk'
  radiation_input_file='null'
  !!antest file root used
  prog_input_file = 'null'
  prog_output_file = trim(root)//"_antestout"
  read(nin,nml=progfiles,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
     print '("Fatal error reading input filenames")'
  end if

  if ( equil_input_file=='null') then
     equil_input_file   = eqdsk_input_file
  end if

  !!check files exist

  if(prog_input_file/='null') then
     inquire(file=prog_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find antest data file, ",a)',prog_input_file
        call log_error(m_name,s_name,2,error_fatal,'antest data file not found')
     end if
  end if

  ! Check files exist
  !! radiation_input_file is optional, not tested here
  if(vtk_input_file/='null') then
     call log_value("geometry field data file, vtk_input_file",trim(vtk_input_file))
     inquire(file=vtk_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Error : Unable to find Geometry data file, ",a)',vtk_input_file
        call log_error(m_name,s_name,3,error_warning,'Geometry data file not found')
     end if
  end if

  if(equil_input_file/='null') then
     call log_value("Beq field data file, equil_input_file",trim(equil_input_file))
     inquire(file=equil_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Error : Unable to find Beq field data file, ",a)',equil_input_file
        call log_error(m_name,s_name,4,error_warning,'Beq field data file not found')
     end if
  end if
  ! set equilibrium type depending on suffix
  call misc_fsuffixget(equil_input_file,icsuf,ierror)
  atnumerics%eqtype=icsuf

  if(antenna_input_file/='null') then
     call log_value("Antenna definition data file, antenna_input_file",trim(antenna_input_file))
     inquire(file=antenna_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Error : Unable to find Antenna definition data file, ",a)',antenna_input_file
        call log_error(m_name,s_name,5,error_warning,'Antenna definition data file not found')
     end if
  end if

  file%vtk   = vtk_input_file
  file%eqdsk   = equil_input_file
  file%equil   = equil_input_file
  file%antdata   = antenna_input_file
  file%radata   = radiation_input_file
  file%antestdata = prog_input_file
  file%antestout = prog_output_file
  file%out = prog_output_file
  !call log_value("antest data file, prog_input_file",trim(file%antestdata))

  !! create output file names from root
  !!gnu file root
  file%gnu = trim(root)//"_gnu"
  file%powx = trim(root)//"_powx"

  !---------------------------------------------------------------------
  !! set default misc parameters
  prog_control = 'standard'
  prog_realpar = .0001_kr8
  prog_intpar = 1
  prog_logicpar = .false.

  !!read misc parameters
  read(nin,nml=miscparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,10,error_fatal,'Error reading misc parameters')
     print '("Fatal error reading misc parameters")'
  end if

  !! check for valid data
  if(prog_control/='standard') then
     call log_error(m_name,s_name,11,error_fatal,'only standard control allowed')
  end if
  ! positive real parameter
  if(prog_realpar<0._kr8) then
     call log_error(m_name,s_name,20,error_fatal,'realpar must be non-negative')
  end if
  ! positive integer parameter
  if(prog_intpar<=0) then
     call log_error(m_name,s_name,30,error_fatal,'intpar must be positive')
  end if
  if(prog_logicpar) &
 &call log_error(m_name,s_name,40,error_warning,'logicpar is true')

  param%control=prog_control
  param%realpar=prog_realpar
  param%intpar=prog_intpar
  param%logicpar=prog_logicpar

  !---------------------------------------------------------------------
  !! set default plot selections
  plot_antestout = .false.
  plot_vtk = .false.
  plot_gnu = .false.
  plot_powx = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,50,error_fatal,'Error reading plot selections')
     print '("Fatal error reading plot selections")'
  end if

  call antest_readcon(atnumerics,nin)

  call antenna_readcon(antnumerics,nin)

  call radnsource_readcon(radnumerics,nin)

  call ray_readcon(raynumerics,nin)

  call beq_readcon(numerics,nin)

  !! store values
  antnumerics%ntheta = atnumerics%ntheta
  antnumerics%nphi = atnumerics%nphi
  plot%antestout = plot_antestout
  plot%vtk = plot_vtk
  plot%gnu = plot_gnu
  plot%powx = plot_powx

end  subroutine atcontrol_read

end module atcontrol_m
