module control_m

  use const_kind_m
  use log_m
  use misc_m
  use position_h
  use control_h
  use position_m

  implicit none
  private

! public subroutines
  public :: &
 &control_init, &
 &control_read, &
 &control_dread, &
 &control_lread, &
 &control_mread, &
 &control_btree


! private variables
  character(*), parameter :: m_name='control_m' !< module name
  integer  :: status   !< error status
  integer  :: nin      !< control file unit number
  integer  :: ilog      !< for namelist dump after error
  integer(ki4)  :: i!< loop counter
  integer(ki4)  :: j        !< loop counter
  character(len=80) :: root !< file root
  contains


!---------------------------------------------------------------------
!! open input control data file

subroutine control_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='control_init' !< subroutine name
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


end  subroutine control_init

!---------------------------------------------------------------------
!! read data for this run

subroutine control_read(file,numerics,plot)

  !! arguments
  type(files_t), intent(out) :: file !< file names
  type(numerics_t), intent(out) :: numerics !< input numerical parameters
  type(plots_t), intent(out) :: plot !< vtk plot selectors

  !!local
  character(*), parameter :: s_name='control_read' !< subroutine name
  character(len=80) :: vtk_input_file  !< geobj data input file
  logical :: filefound !< true of file exists
  integer(ki4) :: geometrical_type !< type of geometry
  integer(ki4) :: min_geobj_in_bin !< numerical parameter
  integer(ki4) :: limit_geobj_in_bin !< numerical parameter
  integer(ki4) :: nbins_child  !< numerical parameter
  real(kr4) :: min_tolerance !< numerical parameter
  real(kr4) :: max_tolerance !< numerical parameter
  integer(ki4) :: boundary_cubes !< number of extra bounding cubes
  real(kr4), dimension(3) :: lower_corner !< vector defining corner with lesser component values
  real(kr4), dimension(3) :: upper_corner !< vector defining corner with greater component values
  real(kr4) :: delta_inner_length   !< inner cube separation from geometry
  real(kr4) :: delta_outer_length   !< outer cube separation from inner
  real(kr4) :: zrhmin !< reciprocal of estimated hmin
  integer(ki4) :: type_geobj_coord_scaling !< type of scaling of geobj x value
  integer(ki4) :: quantising_number !< local variable
  integer(ki4)  :: i2        !< power of two
  integer(ki4) :: no_geobj_records !< local variable
  integer(ki4) :: margin_type !< local variable
  integer(ki4) :: corner_flag !< local variable
  logical :: plot_hds !< DUPLICATE vtk plot selector
  logical :: plot_hdsm !< vtk plot selector
  logical :: plot_hdsbins !< DUPLICATE vtk plot selector
  logical :: plot_hdsq !< vtk plot selector
  logical :: plot_geobj !< vtk plot selector
  logical :: plot_geobjq !< vtk plot selector
  logical :: plot_lostgeobj !< vtk plot selector
  logical :: plot_allgeobj_quantised !< DUPLICATE vtk plot selector
  logical :: plot_geoptq !< vtk plot selector
  logical :: plot_densitygeobj !< vtk plot selector

  !! file names
  namelist /inputfiles/ &
 &vtk_input_file

  !! numerical parameters for hdsgen
  namelist /hdsgenparameters/ &
 &min_geobj_in_bin, &
 &limit_geobj_in_bin, &
 &nbins_child, &
 &geometrical_type, &
 &min_tolerance, &
 &max_tolerance, &
 &boundary_cubes, &
 &lower_corner, &
 &upper_corner, &
 &delta_inner_length, &
 &delta_outer_length, &
 &quantising_number, &
 &no_geobj_records, &
 &margin_type, &
 &type_geobj_coord_scaling

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_hds, &
 &plot_hdsm, &
 &plot_hdsbins, &
 &plot_hdsq, &
 &plot_geobj, &
 &plot_geobjq, &
 &plot_lostgeobj, &
 &plot_allgeobj_quantised, &
 &plot_geoptq, &
 &plot_densitygeobj

  !! read input file names
  read(nin,nml=inputfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_getunit(ilog)
     write(ilog,nml=inputfiles)
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
  end if

  file%vtkdata   = vtk_input_file

  !!check files exist

  call log_value("Objectg data file, vtk_input_file",trim(file%vtkdata))
  inquire(file=vtk_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find geobj data file, ",a)',vtk_input_file
     call log_error(m_name,s_name,2,error_fatal,'Objectg data file not found')
  end if

  !! create output file names from root

  !! output file
  file%hdsgenout = trim(root)//"_hdsgen.out"
  inquire(file=file%hdsgenout,exist=filefound)
  call log_value("geobjgen output file",trim(file%hdsgenout))
  if(filefound) then
     !! error opening file
     print '("Fatal error: Output file ",a, " already exists")',trim(file%hdsgenout)
     print '("Remove it and restart run")'
     call log_error(m_name,s_name,4,error_fatal,'Output file already exists')
  end if

  !!vtk file roots
  file%hdslist     =trim(root)//"_hds"
  file%hdsv  =trim(root)//"_hdsv"
  file%hdsm  =trim(root)//"_hdsm"
  file%hdsq  =trim(root)//"_hdsq"
  file%geobjq =trim(root)//"_geobjq"
  file%lostgeobj    =trim(root)//"_lostgeobj"
  file%allgeobjq    =trim(root)//"_allgeobjq"
  file%geoptq    =trim(root)//"_geoptq"
  file%densitygeobj    =trim(root)//"_densitygeobj"

  !! set default numerical parameters
  min_geobj_in_bin = 10
  limit_geobj_in_bin = 10
  nbins_child = 10
  geometrical_type = 1
  min_tolerance = 1.0e-3
  max_tolerance = 1.0e-1
  boundary_cubes = 0
  lower_corner = 0
  upper_corner = 0
  delta_inner_length = -1.
  delta_outer_length = 0.
  quantising_number =1024
  no_geobj_records = 0
  margin_type = 0
  type_geobj_coord_scaling=2  ! allows for offset

  !!read numerical parameters
  read(nin,nml=hdsgenparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading hdsgen parameters")'
     call log_getunit(ilog)
     write(ilog,nml=hdsgenparameters)
     call log_error(m_name,s_name,5,error_fatal,'Error reading hdsgen parameters')
  end if

  !! check for valid data
  if(geometrical_type<1) &
 &call log_error(m_name,s_name,16,error_fatal,'geometrical_type must be > 0')
  limit_geobj_in_bin=max(limit_geobj_in_bin,min_geobj_in_bin)
  if(limit_geobj_in_bin<0) &
 &call log_error(m_name,s_name,6,error_fatal,'limit_geobj_in_bin must be >=0')
  if(nbins_child<1) &
 &call log_error(m_name,s_name,7,error_fatal,'nbins_child must be > 0')
  if(min_tolerance<0.0 .or. min_tolerance>max_tolerance) then
     call log_error(m_name,s_name,8,error_warning,'invalid min_tolerance, value reset')
     call log_value("min_tolerance",min_tolerance)
  end if
  if(max_tolerance<=0.0)then
     call log_value("max_tolerance",max_tolerance)
     call log_error(m_name,s_name,9,error_fatal,'Invalid max_tolerance, need value > 0')
  end if
  corner_flag = 0
  if(boundary_cubes<0 .OR. boundary_cubes>2 ) &
 &call log_error(m_name,s_name,17,error_fatal,'boundary_cubes must be >= 0 and <=2')
  if(boundary_cubes/=0 ) then
     if(delta_inner_length<0.0)then
        boundary_cubes=1
        corner_flag=1
        call log_value("Assuming number of cuboids defined by corners is",boundary_cubes)
        !call log_value("delta_inner_length",delta_inner_length)
        !call log_error(m_name,s_name,18,error_warning,'Invalid delta_inner_length, need value >= 0')
     end if
     if(delta_outer_length<0.0)then
        call log_value("delta_outer_length",delta_outer_length)
        call log_error(m_name,s_name,19,error_warning,'Invalid delta_outer_length, need value >= 0')
     end if
  end if
  if(quantising_number<2) &
 &call log_error(m_name,s_name,10,error_fatal,'quantising_number must be > 1')
  if(type_geobj_coord_scaling<0) &
 &call log_error(m_name,s_name,11,error_fatal,'type_geobj_coord_scaling must be positive')

  !! store values
  call control_quantise(numerics,quantising_number,type_geobj_coord_scaling)


  numerics%geomtype = geometrical_type
  numerics%mingeobjinbin = limit_geobj_in_bin
  numerics%ngeobj  = no_geobj_records
  numerics%mtype  = margin_type
  numerics%mintolerance = min_tolerance
  numerics%maxtolerance = max_tolerance
  numerics%nbdcub = boundary_cubes
  numerics%cornflag = corner_flag
  numerics%lowcorner = lower_corner
  numerics%upcorner = upper_corner
  numerics%dilen = delta_inner_length
  numerics%dolen = delta_outer_length
  numerics%geobj_coord_tfm%nqtfm=type_geobj_coord_scaling

  call control_btree(numerics,nin)

  call position_readcon(numerics%position_coord_tfm,nin)

  !! set default plot selections
  plot_hds = .false.
  plot_hdsm = .false.
  plot_hdsbins = .false.
  plot_hdsq = .false.
  plot_geobj = .false.
  plot_geobjq = .false.
  plot_lostgeobj = .false.
  plot_allgeobj_quantised = .false.
  plot_geoptq = .false.
  plot_densitygeobj = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     print '("Fatal error reading plot selections")'
     call log_getunit(ilog)
     write(ilog,nml=plotselections)
     call log_error(m_name,s_name,12,error_fatal,'Error reading plot selections')
  end if

  !! store values
  plot%hdsm   = plot_hdsm
  plot%hdsbin    = plot_hdsbins
  plot%hdsq    = plot_hdsq
  plot%geobjq    = plot_geobjq
  plot%lostgeobj     = plot_lostgeobj
  plot%allgeobjq     = plot_allgeobj_quantised
  plot%geoptq     = plot_geoptq
  plot%densitygeobj     = plot_densitygeobj
  if (plot_hds) then
     plot%hdsm = plot_hds
     call log_error(m_name,s_name,20,error_warning,'Obsolete plot selection feature activated')
  end if
  if (plot_geobj) then
     plot%geobjq = plot_geobj
     call log_error(m_name,s_name,21,error_warning,'Obsolete plot selection feature activated')
  end if

end  subroutine control_read
!!---------------------------------------------------------------------
!> control dengen reading of controls
subroutine control_dread(file,numerics,plot)

  !! arguments
  type(files_t), intent(out) :: file !< file names
  type(numerics_t), intent(inout) :: numerics !< input numerical parameters
  type(plots_t), intent(out) :: plot !< vtk plot selectors

  !!local
  character(*), parameter :: s_name='control_dread' !< subroutine name
  character(len=80) :: vtk_input_file  !< geobj data input file
  character(len=80) :: hds_input_file  !< hds data input file
  character(len=80) :: query_input_file  !< qry data input file
  logical :: filefound !< true of file exists
  logical :: plot_density_quantised !< vtk plot selector
  logical :: plot_allgeobj !< vtk plot selector
  logical :: plot_hdsbins !< DUPLICATE vtk plot selector
  logical :: plot_hdsq !< vtk plot selector

  integer(ki4) :: geometrical_type !< type of geometry
  real(kr4) :: min_tolerance !< numerical parameter
  real(kr4) :: max_tolerance !< numerical parameter
  integer(ki4) :: type_geobj_coord_scaling !< type of scaling of geobj x value
  integer(ki4) :: quantising_number !< local variable
  integer(ki4) :: no_geobj_records !< local variable
  integer(ki4) :: margin_type !< local variable

  !! file names
  namelist /inputfiles/ &
 &vtk_input_file, &
 &hds_input_file, &
 &query_input_file


  !! numerical parameters (aka hdsgenparameters)
  namelist /numericalparameters/ &
 &geometrical_type, &
 &min_tolerance, &
 &max_tolerance, &
 &quantising_number, &
 &no_geobj_records, &
 &margin_type, &
 &type_geobj_coord_scaling

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_density_quantised, &
 &plot_allgeobj, &
 &plot_hdsbins, &
 &plot_hdsq

  !! read input file names
  read(nin,nml=inputfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_getunit(ilog)
     write(ilog,nml=inputfiles)
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
  end if

  file%vtkdata   = vtk_input_file
  file%hdsdata   = hds_input_file
  file%qrydata   = query_input_file

  !!check files exist

  call log_value("Objectg data file, vtk_input_file",trim(file%vtkdata))
  inquire(file=vtk_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find geobj data file, ",a)',vtk_input_file
     call log_error(m_name,s_name,2,error_fatal,'Objectg data file not found')
  end if

  call log_value("HDS data file, hds_input_file",trim(file%hdsdata))
  inquire(file=hds_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find hds data file, ",a)',hds_input_file
     call log_error(m_name,s_name,3,error_fatal,'HDS data file not found')
  end if

  call log_value("Query data file, query_input_file",trim(file%qrydata))
  inquire(file=query_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find query data file, ",a)',query_input_file
     call log_error(m_name,s_name,4,error_fatal,'Mesh data file not found')
  end if

  !! set default numerical parameters
  geometrical_type = 1
  min_tolerance = 1.0e-3
  max_tolerance = 1.0e-1
  quantising_number =1024
  no_geobj_records = 0
  margin_type = 0
  type_geobj_coord_scaling=2  ! allows for offset

  !!read numerical parameters
  read(nin,nml=numericalparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading numerical parameters")'
     call log_getunit(ilog)
     write(ilog,nml=numericalparameters)
     call log_error(m_name,s_name,5,error_fatal,'Error reading numerical parameters')
  end if

  !! check for valid data
  if(geometrical_type<1) &
 &call log_error(m_name,s_name,16,error_fatal,'geometrical_type must be > 0')
  if(min_tolerance<0.0 .or. min_tolerance>max_tolerance) then
     call log_error(m_name,s_name,8,error_warning,'invalid min_tolerance, value reset')
     call log_value("min_tolerance",min_tolerance)
  end if
  if(max_tolerance<=0.0)then
     call log_value("max_tolerance",max_tolerance)
     call log_error(m_name,s_name,9,error_fatal,'Invalid max_tolerance, need value > 0')
  end if
  if(quantising_number<32) &
 &call log_error(m_name,s_name,10,error_fatal,'quantising_number must be > 31')

  call control_quantise(numerics,quantising_number,type_geobj_coord_scaling)

  if(type_geobj_coord_scaling<0) &
 &call log_error(m_name,s_name,11,error_fatal,'type_geobj_coord_scaling must be positive')

  !! store values
  numerics%geomtype = geometrical_type
  numerics%ngeobj  = no_geobj_records
  numerics%mintolerance = min_tolerance
  numerics%mtype  = margin_type
  numerics%maxtolerance = max_tolerance
  numerics%geobj_coord_tfm%nqtfm=type_geobj_coord_scaling

  !! create output file names from root

  !! output file
  file%dengenout = trim(root)//"_dengen.out"
  inquire(file=file%dengenout,exist=filefound)
  call log_value("dengen output file",trim(file%dengenout))
  if(filefound) then
     !! error opening file
     print '("Fatal error: Output file ",a, " already exists")',trim(file%dengenout)
     print '("Remove it and restart run")'
     call log_error(m_name,s_name,5,error_fatal,'Output file already exists')
  end if

  !!vtk file roots
  file%den     =trim(root)//"_den"
  file%denq     =trim(root)//"_denq"
  file%allgeobj  =trim(root)//"_allgeobj"
  file%hdsv  =trim(root)//"_hdsv"
  file%hdsm  =trim(root)//"_hdsm"

  !! set default plot selections
  plot_density_quantised = .false.
  plot_allgeobj = .false.
  plot_hdsbins = .false.
  plot_hdsq = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     print '("Fatal error reading plot selections")'
     call log_getunit(ilog)
     write(ilog,nml=plotselections)
     call log_error(m_name,s_name,6,error_fatal,'Error reading plot selections')
  end if

  !! store values
  plot%denq   = plot_density_quantised
  plot%allgeobj    = plot_allgeobj
  plot%hdsbin    = plot_hdsbins
  plot%hdsq    = plot_hdsq

end  subroutine control_dread

!!---------------------------------------------------------------------
!! control legtfm reading of controls

subroutine control_lread(file,numbers,plot)

  !! arguments
  type(files_t), intent(out) :: file !< file names
  type(numbers_t), intent(inout) :: numbers !< input numerical parameters
  type(plots_t), intent(out) :: plot !< vtk plot selectors

  !!local
  character(*), parameter :: s_name='control_lread' !< subroutine name
  character(len=80) :: vtk_input_file  !< result data input file
  character(len=80) :: query_input_file  !< qry data input file
  logical :: filefound !< true if file exists
  logical :: plot_transformed !< vtk plot selector
  logical :: plot_restored !< vtk plot selector

  integer(ki4) :: start_energy_level  !< local variable
  integer(ki4) :: stop_energy_level  !< local variable
  real(kr4) :: max_energy !< numerical parameter
  real(kr4) :: min_energy !< numerical parameter
  integer(ki4) :: sph_degree !< numerical parameter
  real(kr4), dimension(3)  :: src_coord !< source coord
  real(kr4), dimension(3)  :: u_vector !< u coord vector
  real(kr4), dimension(3)  :: v_vector !< v coord vector
  real(kr4)    :: src_strength !< source strength
  character(len=80)  :: src_name  !< source name

  !! file names
  namelist /inputfiles/ &
 &vtk_input_file, &
 &query_input_file


  !! numerical parameters
  namelist /numericalparameters/ &
 &start_energy_level, &
 &stop_energy_level, &
 &max_energy, &
 &min_energy, &
 &sph_degree, &
 &src_coord, &
 &u_vector, &
 &v_vector, &
 &src_strength, &
 &src_name

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_transformed, &
 &plot_restored

  !! read input file names
  read(nin,nml=inputfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_getunit(ilog)
     write(ilog,nml=inputfiles)
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
  end if

  file%vtkdata   = vtk_input_file
  file%qrydata   = query_input_file

  !!check files exist

  call log_value("Result data file, vtk_input_file",trim(file%vtkdata))
  inquire(file=vtk_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find result data file, ",a)',vtk_input_file
     call log_error(m_name,s_name,2,error_fatal,'Result data file not found')
  end if

  call log_value("Query data file, query_input_file",trim(file%qrydata))
  inquire(file=query_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find query data file, ",a)',query_input_file
     call log_error(m_name,s_name,4,error_fatal,'Mesh data file not found')
  end if

  !! set default numerical parameters
  start_energy_level = 1
  stop_energy_level=255  !
  max_energy = 1.0e+8 ! inert
  min_energy = 1.0e-3 ! inert
  sph_degree = 3
  src_coord=(/24.82,0.,-390./) ! cm
  u_vector=(/0.,0.,1./)
  v_vector=(/1.,0.,0./)
  src_strength=1.
  src_name='FZK001'

  !!read numerical parameters
  read(nin,nml=numericalparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading numerical parameters")'
     call log_getunit(ilog)
     write(ilog,nml=numericalparameters)
     call log_error(m_name,s_name,5,error_fatal,'Error reading numerical parameters')
  end if

  !! check for valid data
  if(start_energy_level<1) &
 &call log_error(m_name,s_name,16,error_fatal,'start_energy_level must be > 0')
  if(stop_energy_level<1) &
 &call log_error(m_name,s_name,11,error_fatal,'stop_energy_level must be positive')
  if(max_energy<=0.0)then
     call log_value("max_energy",max_energy)
     call log_error(m_name,s_name,9,error_fatal,'Invalid max_energy, need value > 0')
  end if
  if(min_energy<0.0 .or. min_energy>max_energy) then
     call log_error(m_name,s_name,8,error_warning,'invalid min_energy, value reset')
     call log_value("min_energy",min_energy)
  end if
  if(sph_degree<1) &
 &call log_error(m_name,s_name,16,error_fatal,'sph_degree must be > 0')

  !! store values
  numbers%nestart = start_energy_level
  numbers%nestop = stop_energy_level
  numbers%maxenergy = max_energy
  numbers%minenergy = min_energy
  numbers%degree = sph_degree
  numbers%coord=src_coord
  numbers%uvec=u_vector
  numbers%vvec=v_vector
  numbers%strength=src_strength
  numbers%name=src_name

  !! create output file names from root

  !! output file
  file%legtfmout = trim(root)//"_legtfm.out"
  inquire(file=file%legtfmout,exist=filefound)
  call log_value("legtfm output file",trim(file%legtfmout))
  if(filefound) then
     !! error opening file
     print '("Fatal error: Output file ",a, " already exists")',trim(file%legtfmout)
     print '("Remove it and restart run")'
     call log_error(m_name,s_name,5,error_fatal,'Output file already exists')
  end if

  !! rtt output file
  file%rtt = trim(root)//"_legtfm.rtt"
  inquire(file=file%rtt,exist=filefound)
  call log_value("legtfm output file",trim(file%rtt))

  !!vtk output file root
  file%tfm     =trim(root)//"_tfm"
  file%rest     =trim(root)//"_rest"

  !!gnuplot output file root
  file%gnu     =trim(root)//"_tfm"

  !! set default plot selections
  plot_transformed = .false.
  plot_restored = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     print '("Fatal error reading plot selections")'
     call log_getunit(ilog)
     write(ilog,nml=plotselections)
     call log_error(m_name,s_name,6,error_fatal,'Error reading plot selections')
  end if

  !! store values
  plot%tfm   = plot_transformed
  plot%rest   = plot_restored

end  subroutine control_lread


!!---------------------------------------------------------------------
!! control move reading of controls

subroutine control_mread(file,numerics,plot)

  !! arguments
  type(files_t), intent(out) :: file !< file names
  type(numerics_t), intent(inout) :: numerics !< input numerical parameters
  type(plots_t), intent(out) :: plot !< vtk plot selectors

  !!local
  character(*), parameter :: s_name='control_mread' !< subroutine name
  character(len=80) :: vtk_input_file  !< geobj data input file
  character(len=80) :: hds_input_file  !< hds data input file
  character(len=80) :: query_input_file  !< qry data input file
  logical :: filefound !< true of file exists
  logical :: plot_move_positions !< vtk plot selector
  logical :: plot_points !< vtk plot selector
  logical :: plot_allgeobj !< vtk plot selector
  logical :: plot_hdsbins !< DUPLICATE vtk plot selector
  logical :: plot_hdsq !< vtk plot selector

  integer(ki4) :: geometrical_type !< type of geometry
  real(kr4) :: min_tolerance !< numerical parameter
  real(kr4) :: max_tolerance !< numerical parameter
  integer(ki4) :: type_geobj_coord_scaling !< type of scaling of geobj x value
  integer(ki4) :: quantising_number !< local variable
  integer(ki4) :: no_geobj_records !< local variable
  integer(ki4) :: margin_type !< local variable

  !! file names
  namelist /inputfiles/ &
 &vtk_input_file, &
 &hds_input_file, &
 &query_input_file


  !! numerical parameters aka hdsgenparameters
  namelist /numericalparameters/ &
 &geometrical_type, &
 &min_tolerance, &
 &max_tolerance, &
 &quantising_number, &
 &no_geobj_records, &
 &margin_type, &
 &type_geobj_coord_scaling

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_move_positions, &
 &plot_points, &
 &plot_allgeobj, &
 &plot_hdsbins, &
 &plot_hdsq

  !! read input file names
  read(nin,nml=inputfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_getunit(ilog)
     write(ilog,nml=inputfiles)
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
  end if

  file%vtkdata   = vtk_input_file
  file%hdsdata   = hds_input_file
  file%qrydata   = query_input_file

  !!check files exist

  call log_value("Objectg data file, vtk_input_file",trim(file%vtkdata))
  inquire(file=vtk_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find geobj data file, ",a)',vtk_input_file
     call log_error(m_name,s_name,2,error_fatal,'Objectg data file not found')
  end if

  call log_value("HDS data file, hds_input_file",trim(file%hdsdata))
  inquire(file=hds_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find hds data file, ",a)',hds_input_file
     call log_error(m_name,s_name,3,error_fatal,'HDS data file not found')
  end if

  call log_value("Query data file, query_input_file",trim(file%qrydata))
  inquire(file=query_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find query data file, ",a)',query_input_file
     call log_error(m_name,s_name,4,error_fatal,'Mesh data file not found')
  end if

  !! set default numerical parameters
  geometrical_type = 1
  min_tolerance = 1.0e-3
  max_tolerance = 1.0e-1
  quantising_number =1024
  no_geobj_records = 0
  margin_type = 0
  type_geobj_coord_scaling=2  ! allows for offset

  !!read numerical parameters
  read(nin,nml=numericalparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading numerical parameters")'
     call log_getunit(ilog)
     write(ilog,nml=numericalparameters)
     call log_error(m_name,s_name,5,error_fatal,'Error reading numerical parameters')
  end if

  !! check for valid data
  if(geometrical_type<1) &
 &call log_error(m_name,s_name,16,error_fatal,'geometrical_type must be > 0')
  if(min_tolerance<0.0 .or. min_tolerance>max_tolerance) then
     call log_error(m_name,s_name,8,error_warning,'invalid min_tolerance, value reset')
     call log_value("min_tolerance",min_tolerance)
  end if
  if(max_tolerance<=0.0)then
     call log_value("max_tolerance",max_tolerance)
     call log_error(m_name,s_name,9,error_fatal,'Invalid max_tolerance, need value > 0')
  end if
  if(quantising_number<32) &
 &call log_error(m_name,s_name,10,error_fatal,'quantising_number must be > 31')

  call control_quantise(numerics,quantising_number,type_geobj_coord_scaling)

  if(type_geobj_coord_scaling<0) &
 &call log_error(m_name,s_name,11,error_fatal,'type_geobj_coord_scaling must be positive')

  !! store values
  numerics%geomtype = geometrical_type
  numerics%ngeobj  = no_geobj_records
  numerics%mintolerance = min_tolerance
  numerics%mtype  = margin_type
  numerics%maxtolerance = max_tolerance
  numerics%geobj_coord_tfm%nqtfm=type_geobj_coord_scaling

  !! create output file names from root

  !! output file
  file%moveout = trim(root)//"_move.out"
  inquire(file=file%moveout,exist=filefound)
  call log_value("move output file",trim(file%moveout))
  if(filefound) then
     !! error opening file
     print '("Fatal error: Output file ",a, " already exists")',trim(file%moveout)
     print '("Remove it and restart run")'
     call log_error(m_name,s_name,5,error_fatal,'Output file already exists')
  end if

  !!vtk file roots
  file%mov     =trim(root)//"_mov"
  file%movq     =trim(root)//"_movq"
  file%allgeobj  =trim(root)//"_allgeobj"
  file%hdsv  =trim(root)//"_hdsv"
  file%hdsm  =trim(root)//"_hdsm"
  file%pts     =trim(root)//"_pts"
  file%ptsq     =trim(root)//"_ptsq"

  !! set default plot selections
  plot_move_positions = .false.
  plot_points = .false.
  plot_allgeobj = .false.
  plot_hdsbins = .false.
  plot_hdsq = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     print '("Fatal error reading plot selections")'
     call log_getunit(ilog)
     write(ilog,nml=plotselections)
     call log_error(m_name,s_name,6,error_fatal,'Error reading plot selections')
  end if

  !! store values
  plot%movq   = plot_move_positions
  plot%allgeobj    = plot_allgeobj
  plot%hdsbin    = plot_hdsbins
  plot%hdsq    = plot_hdsq
  plot%ptsq    = plot_points

end  subroutine control_mread

subroutine control_btree(numerics,kin)

  !! arguments
  type(numerics_t), intent(out) :: numerics   !< binary tree numeric controls
  integer(ki4) :: kin  !< local variable


  !! local
  character(*), parameter :: s_name='control_btree' !< subroutine name
  integer(ki4):: btree_size !< size of binary tree
  integer(ki4):: btree_sizep !< size of binary tree pter
  integer(ki4):: btree_sizee !< size of exten array
  integer(ki4):: btree_sizeh !< size of list array hoc
  integer(ki4):: btree_sizel !< size of list array
  integer(ki4):: btree_depth !< max depth of tree
  integer(ki4):: tree_ttalg !< type of tree algorithm
  integer(ki4), dimension(3) :: tree_nxyz !< top of tree children
  real(kr4), dimension(3) :: tree_hxyz !< top of tree spacings
  integer(ki4):: tree_type !< type of tree

  !! btree parameters
  namelist /btreeparameters/ &
 &btree_size , &
 &btree_sizep , &
 &btree_sizee , &
 &btree_sizeh , &
 &btree_sizel , &
 &btree_depth , &
 &tree_ttalg , &
 &tree_nxyz , &
 &tree_hxyz , &
 &tree_type

  !! set default btree parameters
  btree_size=2000000 !max_size_tree_array
  btree_sizep=3 !! BSP
  btree_sizee=1000 !max_size_exten_array)
  btree_sizeh=2
  btree_sizel=2000000 ! max_size_list_array
  btree_depth=30
  tree_type=1 !! BSP 1 -BSP 2-octree 3-octree with special top
  tree_ttalg=1 !! for special top: 0 -default 1- use hxyz 2- use nxyz
  tree_nxyz=(/2,2,2/) !! for special top octree default
  tree_hxyz=(/0.1,0.1,0.1/) !! for special top octree default

  !!read btree parameters
  read(kin,nml=btreeparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading btree parameters")'
     call log_getunit(ilog)
     write(ilog,nml=btreeparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading btree parameters')
  end if


  !! check for valid data
  if(btree_size<1) &
 &call log_error(m_name,s_name,2,error_fatal,'btree_size must be > 0')
  if(btree_sizep<1) &
 &call log_error(m_name,s_name,3,error_fatal,'btree_sizep must be > 0')
  if(btree_sizee<1) &
 &call log_error(m_name,s_name,4,error_fatal,'btree_sizee must be > 0')
  if(btree_sizeh<1) &
 &call log_error(m_name,s_name,5,error_fatal,'btree_sizeh must be > 0')
  if(btree_sizel<1) &
 &call log_error(m_name,s_name,6,error_fatal,'btree_sizel must be > 0')
  if(btree_depth<1) &
 &call log_error(m_name,s_name,6,error_fatal,'btree_depth must be > 0')
  if(tree_type<1.OR.tree_type>3) &
 &call log_error(m_name,s_name,6,error_fatal,'tree_type must be > 0 and <=3')
  if(tree_ttalg<0.OR.tree_ttalg>2) &
 &call log_error(m_name,s_name,7,error_fatal,'tree_ttalg must be >= 0 and <=2')
  if(any(tree_nxyz<1)) &
 &call log_error(m_name,s_name,9,error_fatal,'all tree_nxyz must be > 0')
  if(any(tree_hxyz<0.)) &
 &call log_error(m_name,s_name,10,error_fatal,'all tree_hxyz must be > 0')

  if (tree_type==1) then
     ! BSP
     if (btree_depth>3*numerics%nquante) then
        call log_error(m_name,s_name,11,error_warning,'depth too great for quantising_number,reset')
        btree_depth=3*numerics%nquante
     end if
  else if (tree_type==2) then
     ! standard octree
     ! test depth not too great
     if (btree_depth>numerics%nquante) then
        call log_error(m_name,s_name,12,error_warning,'depth too great for quantising_number, reset')
        btree_depth=numerics%nquante
     end if
  end if

  !! store values
  numerics%nsize=btree_size
  numerics%nsizep=btree_sizep
  numerics%nsizee=btree_sizee
  numerics%nsizeh=btree_sizeh
  numerics%nsizel=btree_sizel
  numerics%ndepth=btree_depth
  numerics%nttype=tree_type
  numerics%nttalg=tree_ttalg
  numerics%nxyz=tree_nxyz
  numerics%hxyz=tree_hxyz

end subroutine control_btree

subroutine control_quantise(numerics,kquante,kqtfm)

  !! arguments
  type(numerics_t), intent(out) :: numerics   !< numeric controls
  integer(ki4), intent(in):: kquante !< quantising number
  integer(ki4), intent(in):: kqtfm  !< local variable

  !! local
  character(*), parameter :: s_name='control_quantise' !< subroutine name
  real(kr4):: zrhmin  !< local variable
  integer(ki4):: i2  !< local variable
  integer(ki2):: iquante !< log base two of quantising_number

  ! quantisation
  i2=1
  do j=1,31
     i2=2*i2
     if (i2>=kquante) then
        ! log2 of quantising_number
        iquante=j
        zrhmin=float(i2)
        exit
     end if
  end do

  if (iquante.gt.ki2bits-1) then
     ! too large for type ki2 integers
     call log_error(m_name,s_name,1,error_warning,'quantising number too large, reset')
     iquante=ki2bits-2
  end if

  !! store values
  numerics%nquante  = iquante

end subroutine control_quantise


end module control_m
