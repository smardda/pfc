module scontrol_m

  use const_kind_m
  use log_m
  use misc_m
  use scontrol_h

  implicit none
  private

! public subroutines
  public :: &
 &scontrol_init, & !< open input control data file
 &scontrol_close, & !< close input control data file
 &scontrol_read, & !< read data for this run
 &scontrol_readcon !< read data for this object

! public variables
  integer(ki4), parameter, public :: RECOGNISED_KEYS=4 !< number of rekeys recognised by code
  character(len=6),  dimension(RECOGNISED_KEYS), parameter, public :: reckey = & !< input variable must match
 &(/'null  ','angle ','poloid','toroid'/)

! private variables
  character(*), parameter :: m_name='scontrol_m' !< module name
  integer  :: status   !< error status
  integer, save  :: nin      !< control file unit number
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  character(len=80) :: root !< file root
  integer(ki4), parameter :: RECOGNISED_STATISTICS=6 !< number of statistics recognised by code
  character(len=6), dimension(RECOGNISED_STATISTICS), parameter, public :: statistic = & !< input variable must match
 &(/'intg  ','intgsq','max   ','min   ','maxabs','minabs'/)

  contains
!---------------------------------------------------------------------
!> open input control data file
subroutine scontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='scontrol_init' !< subroutine name
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

end  subroutine scontrol_init
!---------------------------------------------------------------------
!> close input control data file
subroutine scontrol_close

  !! local
  character(*), parameter :: s_name='scontrol_close' !< subroutine name

  !! close file unit
  close(unit=nin)

end  subroutine scontrol_close
!---------------------------------------------------------------------
!> read data for this run
subroutine scontrol_read(param,file,plot)

  !! arguments
  type(snumerics_t), intent(out) :: param !< control parameters
  type(sfiles_t), intent(out) :: file !< file names
  type(splots_t), intent(out) :: plot !< plot selectors

  !!local
  character(*), parameter :: s_name='scontrol_read' !< subroutine name
  logical :: filefound !< true if file exists

  integer(ki4), parameter  :: maximum_number_of_files=100 !< maximum number of files allowed
  character(len=80), dimension(maximum_number_of_files) :: vtk_data_file !< vtk format geometry data input file
  character(len=80) :: vtk_input_file  !< smanal data input file _powx format
  character(len=80) :: vtk_output_file  !< smanal data output file
  character(len=80) :: gnu_output_file  !< smanal gnu/csv data output file
  character(len=80) :: gnure_output_file  !< smanal gnu/csv data output file
  character(len=80) :: vtk_small_output_file  !< smanal data output at centroids file

  character(len=80) :: prog_control  !< INERT option control parameter
  real(kr8) :: prog_realpar !< INERT real control parameter
  integer(ki4) :: infil  !< for checking files

  logical :: plot_smanalout !< smanal output data selector
  logical :: plot_anx !< smanal vtk plot selector
  logical :: plot_ansmallx !< smanal vtk plot at centroids only
  logical :: plot_angnu !< smanal  gnuplot plot selector

  !! file names
  namelist /analysisfiles/ &
 &vtk_data_file, &
 &vtk_input_file, &
 &gnu_output_file, &
 &gnure_output_file, &
 &vtk_small_output_file, &
 &vtk_output_file

  !! misc parameters
  namelist /miscparameters/ &
 &prog_control,&
 &prog_realpar

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_smanalout, &
 &plot_anx, &
 &plot_ansmallx, &
 &plot_angnu

  !---------------------------------------------------------------------
  !! read input file names
  vtk_data_file='null'
  vtk_input_file='null'
  vtk_output_file='null'
  vtk_small_output_file='null'
  gnu_output_file='null'
  gnure_output_file='null'
  read(nin,nml=analysisfiles,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
     print '("Fatal error reading input filenames")'
  end if
  ! check number of files against actuals
  do j=1,maximum_number_of_files
     if (trim(vtk_data_file(j))=='null') exit
     infil=j
  end do
  file%nvtkdata=infil
  call log_value("actual number of files",infil)

  file%vtkfull = vtk_output_file
  file%vtksmall = vtk_small_output_file

  !!check files exist
  do j=1,infil
     call log_value("geometry data file number",j)
     call log_value("geometry data file, vtk_data_file",trim(vtk_data_file(j)))
     if(vtk_data_file(j)/='null') then
        inquire(file=vtk_data_file(j),exist=filefound)
        if(.not.filefound) then
           !! error opening file
           print '("Fatal error: Unable to find  data file, ",a)',vtk_data_file(j)
           call log_error(m_name,s_name,22,error_fatal,' data file not found')
        end if
     end if
  end do

  ! allocate and assign files
  allocate(file%vtkdata(infil), stat=status)
  call log_alloc_check(m_name,s_name,23,status)
  do j=1,infil
     file%vtkdata(j)   = vtk_data_file(j)
  end do

  inquire(file=vtk_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find smanal data file, ",a)',vtk_input_file
     call log_error(m_name,s_name,2,error_fatal,'smanal data file not found')
  end if
  file%vtk = vtk_input_file

  !! create output file names from root where necessary

  if (vtk_output_file=='null') then
     file%vtkfull = trim(root)//"_anx"
  else
     file%vtkfull = vtk_output_file
  end if
  if (vtk_small_output_file=='null') then
     file%vtksmall = trim(root)//"_ansmx"
  else
     file%vtksmall = vtk_small_output_file
  end if
  !!gnu file roots
  if (gnu_output_file=='null') then
     file%gnu = trim(root)//"_angnu"
  else
     file%gnu = gnu_output_file
  end if
  if (gnure_output_file=='null') then
     file%gnure = trim(root)//"_angnure"
  else
     file%gnure = gnure_output_file
  end if

  !!smanal file root
  file%smanalout = trim(root)//"_smanal.out"
  filefound=.false.
  inquire(file=file%smanalout,exist=filefound)
  call log_value("smanal output file",trim(file%smanalout))
  if(filefound) then
     !! error opening file
     print '("Fatal error: Output file ",a, " already exists")',trim(file%smanalout)
     print '("Remove it and restart run")'
     call log_error(m_name,s_name,4,error_fatal,'Output file already exists')
  end if

  !---------------------------------------------------------------------
  !! set default misc parameters
  prog_control = 'standard'
  prog_realpar = .0001_kr8

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

  param%control=prog_control
  param%realpar=prog_realpar

  !---------------------------------------------------------------------
  !! set default plot selections
  plot_smanalout = .false.
  plot_anx = .false.
  plot_ansmallx = .false.
  plot_angnu = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,50,error_fatal,'Error reading plot selections')
     print '("Fatal error reading plot selections")'
  end if

  !! store values
  plot%smanalout = plot_smanalout
  plot%vtk = plot_anx
  plot%vtksmall = plot_ansmallx
  plot%gnu = plot_angnu

  call scontrol_readcon(param,nin)

end  subroutine scontrol_read
!---------------------------------------------------------------------
!> read data for this object
subroutine scontrol_readcon(self,channel)

  !! arguments
  type(snumerics_t), intent(inout) :: self !< object control data structure
  integer, intent(in), optional :: channel   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='scontrol_readcon' !< subroutine name
  character(len=80) :: analysis_mode !< only one mode available
  character(len=80) :: sort_key !< key to be used to sort scalar field
  character(len=80) :: analyse_scalar !< identifier of scalar field to analyse
  logical :: set_r_central !< user sets nominal central \f$ R \f$ for poloidal angle analysis
  logical :: set_z_central !< user sets nominal central \f$ Z \f$ for poloidal angle analysis
  real(kr8) :: nominal_r_central !< user defined value of central \f$ R \f$
  real(kr8) :: nominal_z_central !< user defined value of central \f$ Z \f$
  integer(ki4) :: number_of_clusters !< sets cluster size by dividing range
  character(len=80) :: new_key !<  new key for sorting (only angle allowed)
  logical :: rekey !<  only rekey, i.e. use new key, if .true.
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=100 !< maximum number of parameters allowed

  character(len=80), dimension(MAX_NUMBER_OF_PARAMETERS) :: required_statistics  !< statistics to be calculated
  integer(ki4) :: itotstat  !< actual number of statistics needed
  integer(ki4) :: ierr  !< for checking statistics
  character(len=80) :: icstat  !< local variable

  !! smanal parameters
  namelist /smanalparameters/ &
 &analysis_mode, sort_key, &
 &analyse_scalar, number_of_clusters, &
 &set_r_central, set_z_central, &
 &nominal_r_central,nominal_z_central, &
 &new_key, &
 &rekey, &
 &required_statistics

  !! set default smanal parameters
  analysis_mode='regular'
  sort_key='Body'
  analyse_scalar='Q'
  set_r_central=.false.
  set_z_central=.false.
  nominal_r_central=5000.
  nominal_z_central=0.
  number_of_clusters = 18
  new_key = 'null'
  rekey = .false.
  do i=1,MAX_NUMBER_OF_PARAMETERS
     required_statistics(i)='unset '
  end do
  itotstat=0

  if(present(channel).AND.channel/=0) then
     !! assume unit already open and reading infile
     nin=channel
  end if

  !!read smanal parameters
  read(nin,nml=smanalparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading smanal parameters")'
     call log_error(m_name,s_name,1,error_fatal,'Error reading smanal parameters')
  end if

  ! positive integer parameter
  if(number_of_clusters<=0) then
     call log_error(m_name,s_name,2,error_fatal,'(number_of_clusters  must be positive')
  end if
  ! Check new key
  call lowor(new_key,1,len_trim(icstat))
  ierr=1
  do j=1,RECOGNISED_KEYS
     if (new_key(1:6)==reckey(j)) then
        ierr=0
        exit
     end if
  end do
  if (ierr==1) then
     call log_value("new key ",new_key)
     call log_error(m_name,s_name,3,error_fatal,'Unrecognised new key')
  end if

  ! work out how many statistics required
  do i=1,MAX_NUMBER_OF_PARAMETERS
     if (required_statistics(i)=='unset ') exit
     itotstat=itotstat+1
  end do

  do i=1,itotstat
     icstat=required_statistics(i)
     call lowor(icstat,1,len_trim(icstat))
     ierr=1
     do j=1,RECOGNISED_STATISTICS
        if (icstat(1:6)==statistic(j)) then
           ierr=0
           exit
        end if
     end do
     if (ierr==1) then
        call log_value("Statistic ",required_statistics(i))
        call log_error(m_name,s_name,4,error_fatal,'Unrecognised statistic')
     else
        required_statistics(i)=icstat
     end if
  end do

  !! store values
  self%namekey=sort_key
  self%namescal=analyse_scalar
  self%lurcen=set_r_central
  self%luzcen=set_z_central
  if (set_r_central) self%urcen=nominal_r_central
  if (set_z_central) self%uzcen=nominal_z_central
  self%nbin=number_of_clusters
  self%newkey=new_key
  self%rekey=rekey
  self%mode=analysis_mode

  !! allocate arrays and assign

  self%totstat=itotstat

  mode_selection: select case (analysis_mode)

  case('regular')
     if (itotstat>0) then
        allocate(self%namstat(itotstat), stat=status)
        call log_alloc_check(m_name,s_name,10,status)
        do i=1,itotstat
           self%namstat(i)=required_statistics(i)
        end do
     end if
  case default
  end select mode_selection

end  subroutine scontrol_readcon

end module scontrol_m
