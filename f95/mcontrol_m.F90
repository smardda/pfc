module mcontrol_m

  use const_kind_m
  use log_m
  use misc_m
  use position_h
  use mcontrol_h
  use position_m

  implicit none
  private

! public subroutines
  public :: &
 &mcontrol_init, & !< open input control data file
 &mcontrol_close, & !< close input control data file
 &mcontrol_read !< read data for this run

! private variables
  character(*), parameter :: m_name='mcontrol_m' !< module name
  integer  :: status   !< error status
  integer,save  :: nin      !< control file unit number
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
subroutine mcontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='mcontrol_init' !< subroutine name
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


end  subroutine mcontrol_init
!---------------------------------------------------------------------
!> close input control data file
subroutine mcontrol_close

  !! arguments
  !! local
  character(*), parameter :: s_name='mcontrol_close' !< subroutine name

  !! close file unit
  close(unit=nin)

end  subroutine mcontrol_close
!---------------------------------------------------------------------
!> read data for this run
subroutine mcontrol_read(file,numerics,plot)

  !! arguments
  type(mfiles_t), intent(out) :: file !< file names
  type(mnumerics_t), intent(out) :: numerics !< input numerical parameters
  type(mplots_t), intent(out) :: plot !< plot selectors

  !!local
  character(*), parameter :: s_name='mcontrol_read' !< subroutine name
  character(len=80) :: mag_input_file  !< mag format magnetic field data input file
  character(len=80) :: mag_output_file  !< mag format magnetic field data output file
  logical :: filefound !< true of file exists
  real(kr8) :: radius_of_axis !< \f$ R \f$  axis value
  real(kr8) :: i_requested !< \f$ B R \f$  to give to output field
  integer(ki4) :: vacuum_field_axis_calculation    !< how vacuum field is determined
  integer(ki4) :: spline_order    !< order of spline
  integer(ki4) :: split_transform    !< order of transform splitting
  integer(ki4), dimension(3) :: min_wavenumber !< minimum wavenumber for evaluation
  integer(ki4), dimension(3) :: max_wavenumber !< maximum wavenumber for evaluation
  integer(ki4), dimension(3) :: mode_parity !< parity (1=even,2=odd,3=both)
  real(kr8) :: mode_cutoff !< relative moe cut-off amplitude
  logical :: plot_cartv !< vtk plot selector
  logical :: plot_allcartv !< vtk plot selector
  logical :: plot_testtfm !< do test of FFT analysis
  logical :: plot_gnuv !< gnu 2-D plot selector
  logical :: plot_modes !< gnu 2-D plot selector
  logical :: plot_gnu !< gnu 1-D plot selector

  !! file names
  namelist /magfiles/ &
 &mag_input_file, &
 &mag_output_file

  !! misc parameters
  namelist /miscparameters/ &
 &radius_of_axis,vacuum_field_axis_calculation,spline_order,split_transform,&
 &min_wavenumber,max_wavenumber,&
 &mode_parity, mode_cutoff, &
 &i_requested

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_cartv, &
 &plot_allcartv, &
 &plot_testtfm, &
 &plot_gnuv, &
 &plot_modes, &
 &plot_gnu

  !! read input file names
  mag_input_file='null'
  read(nin,nml=magfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_getunit(ilog)
     write(ilog,nml=magfiles)
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
  end if

  file%magdata   = mag_input_file
  file%magout   = mag_output_file

  !!check file exists

  call log_value("magnetic field data file, mag_input_file",trim(file%magdata))
  if(file%magdata/='null') then
     inquire(file=mag_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find field data file, ",a)',mag_input_file
        call log_error(m_name,s_name,3,error_fatal,'field data file not found')
     end if
  end if

  !! create output file names from root

  !! output file
  file%mout = trim(root)//"_v.out"
  filefound=.false.
  !     inquire(file=file%mout,exist=filefound)
  !     call log_value("magtfm output file",trim(file%mout))
  if(filefound) then
     !! error opening file
     print '("Fatal error: Output file ",a, " already exists")',trim(file%mout)
     print '("Remove it and restart run")'
     call log_error(m_name,s_name,4,error_fatal,'Output file already exists')
  end if

  !!vtk file roots
  file%cartv     =trim(root)//"_cartv"
  file%allcartv     =trim(root)//"_allcartv"
  !!gnu file roots
  file%modes     =trim(root)//"_modes"
  file%gnuv     =trim(root)//"_gnuv"
  file%gnu     =trim(root)//"_gnu"

  !! set default misc parameters
  radius_of_axis = 0.8
  vacuum_field_axis_calculation = 1
  spline_order = 4
  split_transform = 2
  min_wavenumber = (/1,1,0/)
  max_wavenumber = 1
  mode_parity = (/2,2,1/)
  mode_cutoff = .0001_kr8
  i_requested = 1._kr8

  !!read misc parameters
  read(nin,nml=miscparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading misc parameters")'
     call log_getunit(ilog)
     write(ilog,nml=miscparameters)
     call log_error(m_name,s_name,10,error_fatal,'Error reading misc parameters')
  end if

  !! check for valid data
  if(vacuum_field_axis_calculation<=0) &
 &call log_error(m_name,s_name,11,error_fatal,'vacuum field calculation type must be >0')
  if(vacuum_field_axis_calculation==1 .AND.  radius_of_axis<=0) &
 &call log_error(m_name,s_name,12,error_fatal,'radius of axis must be >0')

  if(spline_order/=4) &
 &call log_error(m_name,s_name,13,error_warning,'only cubic splines allowed')
  !    &call log_error(m_name,s_name,13,error_fatal,'only cubic splines allowed')
  if(split_transform/=2) &
 &call log_error(m_name,s_name,14,error_fatal,'only double split allowed')
  if(minval(max_wavenumber)>=0) then
     ! assume explicit wavenumber limit
     if(minval(min_wavenumber)<0) &
 &   call log_error(m_name,s_name,15,error_fatal,'min wave number must be non-negative')
     if(minval(max_wavenumber-min_wavenumber)<0) &
 &   call log_error(m_name,s_name,16,error_warning,'max wave number less than min')
  else
     ! controlled by cutoff parameter
     if(minval(mode_parity)<1) &
 &   call log_error(m_name,s_name,17,error_fatal,'parities must be positive')
     if(maxval(mode_parity)>3) &
 &   call log_error(m_name,s_name,18,error_fatal,'parities must be positive less than three')
     if(mode_cutoff<0) &
 &   call log_error(m_name,s_name,19,error_fatal,'mode cutoff must be  non-negative')
  end if
  if(i_requested==0.) &
 &call log_error(m_name,s_name,20,error_warning,'expected non-zero I = B R product')

  numerics%rax=radius_of_axis
  numerics%magtfm=vacuum_field_axis_calculation
  numerics%nord=spline_order
  numerics%n0=split_transform
  numerics%kmin=min_wavenumber
  numerics%kmax=max_wavenumber
  numerics%parity=mode_parity
  numerics%cutoff=mode_cutoff
  numerics%ireq=i_requested

  !! set default plot selections
  plot_cartv = .false.
  plot_allcartv = .false.
  plot_gnuv = .false.
  plot_modes = .false.
  plot_gnu = .false.
  plot_testtfm = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     print '("Fatal error reading plot selections")'
     call log_getunit(ilog)
     write(ilog,nml=plotselections)
     call log_error(m_name,s_name,22,error_fatal,'Error reading plot selections')
  end if

  !! store values
  plot%cartv   = plot_cartv
  plot%allcartv   = plot_allcartv
  plot%gnuv     = plot_gnuv
  plot%modes     = plot_modes
  plot%gnu     = plot_gnu
  plot%testtfm     = plot_testtfm

end  subroutine mcontrol_read

end module mcontrol_m
