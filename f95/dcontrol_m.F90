module dcontrol_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use misc_m
  use geobj_m
  use dcontrol_h

  implicit none
  private

! public subroutines
  public :: &
 &dcontrol_init, & !< open input control data file
 &dcontrol_closex, & !< close input control data file and  unit
 &dcontrol_close, & !< close input control data file, keep unit number
 &dcontrol_read, & !< read data for this run
 &dcontrol_readnum, & !< read dnumerics namelist data
 &dcontrol_readprogfiles, & !< read progfiles namelist data
 &dcontrol_readatfile, & ! < read dnumerics datafile
 &dcontrol_lines2d, & !< wrapper for misc_lines2d
 &dcontrol_getunit,  & !< get unit number
 &dcontrol_delete !< delete object

! private subroutines

! private variables
  character(*), parameter :: m_name='dcontrol_m' !< module name
  integer  :: status   !< error status
  integer, save  :: nin      !< control file unit number
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  real(kr8) :: zdum !< dummy real
  character(len=80) :: root !< file root

  contains
!---------------------------------------------------------------------
!> open input control data file
subroutine dcontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='dcontrol_init' !< subroutine name
  character(len=80) :: controlfile !< control file name

  call  misc_getfileunit(nin)

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

end  subroutine dcontrol_init
!---------------------------------------------------------------------
!> close input control data file and flag by setting unit number negative
subroutine dcontrol_closex

  !! local
  character(*), parameter :: s_name='dcontrol_closex' !< subroutine name

  !! close file unit
  close(unit=nin,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close file unit, ",i5)',nin
     call log_error(m_name,s_name,1,error_fatal,'Cannot close data file')
     stop
  end if

  nin=-2

end  subroutine dcontrol_closex
!---------------------------------------------------------------------
!> close input control data file, keep unit number
subroutine dcontrol_close

  !! local
  character(*), parameter :: s_name='dcontrol_close' !< subroutine name

  !! close file unit
  close(unit=nin,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close file unit, ",i5)',nin
     call log_error(m_name,s_name,1,error_fatal,'Cannot close data file')
     stop
  end if

end  subroutine dcontrol_close
!---------------------------------------------------------------------
!> read data for this run
subroutine dcontrol_read(file,numerics,plot)

  !! arguments
  type(dfiles_t), intent(out) :: file !< file names
  type(dnumerics_t), intent(out) :: numerics !< control numerics
  type(dplots_t), intent(out), optional :: plot !< plot selectors

  !!local
  character(*), parameter :: s_name='dcontrol_read' !< subroutine name
  logical, parameter :: noplotnml=.TRUE. !< shut off plot namelist read
  logical :: filefound !< true if file exists

  logical :: filedata !< data file specifies geometry

  logical :: plot_vtk !< vtk plot selector
  logical :: plot_gnu !< gnuplot ploteselector

  real(kr8), dimension(3) :: zr    !<   \f$ R \f$
  real(kr8), dimension(3) :: zz    !<   \f$ Z \f$
  real(kr4) :: lenfac=1. !< length scale factor default is m

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_vtk, &
 &plot_gnu

  !---------------------------------------------------------------------
  !! prologue
  !! create output file names from root
  !!vtk file
  file%vtk = trim(root)//"_vtk"
  !!gnu file roots
  file%gnu = trim(root)//"_gnu"

  !---------------------------------------------------------------------
  !!  read namelist containing input file names
  call dcontrol_readprogfiles(file,nin)

  !---------------------------------------------------------------------
  !! read transformation parameters in numerics namelist
  !! may change flag whether there is additional data to be read from file
  filedata=.TRUE.
  call dcontrol_readnum(numerics,nin,filedata)

  !---------------------------------------------------------------------
  !! set default plot selections
  plot_vtk = .false.
  plot_gnu = .false.

  if (noplotnml) then
     !!read plot selections
     read(nin,nml=plotselections,iostat=status)
     if(status/=0) then
        print '("Fatal error reading plot selections")'
        call log_getunit(ilog)
        write(ilog,nml=plotselections)
        call log_error(m_name,s_name,50,error_fatal,'Error reading plot selections')
     end if
  end if

  if (present(plot)) then
     !! store values
     plot%vtk = plot_vtk
     plot%gnu = plot_gnu
  end if

  call dcontrol_close

  !---------------------------------------------------------------------
  !! read silhouette files and process
  if (filedata) then

     ! read
     call dcontrol_readatfile(file,numerics)

     !dbg write(*,*) numerics%r,numerics%z,numerics%npos,numerics%ldiv
     ! process
     if (numerics%ldiv>numerics%npos+1) then
        !! replace with uniformly divided line between (r(1),z(1)) and (r(npos),z(npos))
        zr(1)=numerics%r(1)
        zr(2)=numerics%r(numerics%npos)
        zz(1)=numerics%z(1)
        zz(2)=numerics%z(numerics%npos)
        deallocate(numerics%r,numerics%z)
        call misc_line2d(numerics%r,numerics%z,zr,zz,numerics%ldiv)
        numerics%npos=numerics%ldiv+1
        call log_error(m_name,s_name,1,error_warning,'r/z/rz data replaced with straight line')
     else
        !! divide up each segment of line
        call misc_lines2d(numerics%r,numerics%z,numerics%npos,numerics%ldiv)
     end if
     !dbg write(*,*) numerics%r,numerics%z,numerics%npos,numerics%ldiv

     !! scale if required
     if (numerics%cunits==0) then
        !! input was in m, but matching to CAD in mm
        lenfac=1000.
        do j=1,numerics%npos
           numerics%r(j)=lenfac*numerics%r(j)
           numerics%z(j)=lenfac*numerics%z(j)
        end do
        numerics%cunits=-3
     end if

  end if

  call dcontrol_close

end  subroutine dcontrol_read
!---------------------------------------------------------------------
!> read only dnumerics data
subroutine dcontrol_readnum(numerics,kunit,klfile)

  !! arguments
  type(dnumerics_t), intent(out) :: numerics !< control numerics
  integer, intent(in)  :: kunit     !< file unit number
  logical, intent(inout) :: klfile !< data file specifies geometry

  !!local
  character(*), parameter :: s_name='dcontrol_readnum' !< subroutine name

  character(len=80) :: description  !< surface interaction with particles
  logical :: filedata  !< data file specifies geometry
  character(len=80) :: transform_type  !< transform type
  real(kr8) :: start_angle !< start angle
  real(kr8) :: finish_angle !< finish angle
  real(kr8), dimension(3) :: start_position !< start position
  real(kr8), dimension(3) :: finish_position !< finish position
  integer(ki4) :: number_of_divisions !< number of divisions
  integer(ki4) :: coordinate_system !< coordinate system positive as posang (q.v.), unity for polars
  integer(ki4) :: line_divisions !< number of divisions
  character(len=20) :: angle_units  !< units, either radian(s) or degree(s)
  character(len=20) :: length_units  !< units, either me(tres) or mm
  real(kr4) :: angfac=1. !< angles default is radians
  real(kr4) :: lenfac=1. !< use to scale lengths to mm if input im metres
  integer(ki2par) :: igcode   !< code for geometry

  !! transformation parameters
  namelist /datvtkparameters/ &
 &angle_units,&
 &description,&
 &filedata,&
 &length_units,&
 &transform_type,&
 &start_angle, finish_angle, &
 &start_position, finish_position, &
 &number_of_divisions,&
 &coordinate_system,&
 &line_divisions

  !---------------------------------------------------------------------
  !! set default datvtk parameters
  !polars
  angle_units='degree'
  description='absorb'
  filedata=klfile
  length_units='mm'
  transform_type = 'rotate'
  start_angle = 0
  finish_angle = 90
  start_position = (/1,0,0/)
  finish_position = (/1,1,0/)
  line_divisions = 0
  number_of_divisions = 10
  coordinate_system = 1

  !!read datvtk parameters
  read(kunit,nml=datvtkparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading datvtk parameters")'
     call log_getunit(ilog)
     write(ilog,nml=datvtkparameters)
     call log_error(m_name,s_name,9,error_fatal,'Error reading datvtk parameters')
  end if

  igcode=-99
  surf_desc: select case (description(1:6))
  case('absorb')
     igcode=GEOBJ_ABSORB !  first wall (absorbing boundary)
  case('invisi')
     igcode=GEOBJ_INVISI !  transparent geometry, totally ignored
  case('skylit')
     igcode=GEOBJ_SKYLIT !  skylight
  case('escape')
     igcode=GEOBJ_ESCAPE !  escape boundary (loss)
  case('errlos')
     igcode=GEOBJ_ERRLOS !  beancan (error loss)
  case('cutout')
     igcode=GEOBJ_CUTOUT !  cutout (expected loss)
  end select surf_desc
  if(igcode==-99) then
     call log_error(m_name,s_name,10,error_warning,'description not recognised')

     igcode=GEOBJ_ABSORB
  end if

  !! check for valid data
  if(transform_type/='rotate'.AND.transform_type/='translate') then
     call log_error(m_name,s_name,11,error_fatal,'transform type not recognised')
  end if
  ! Check angles
  !! set up factor for angles
  if (angle_units(1:6)=='degree') angfac=const_degrad
  if(start_angle*angfac<-2.000001_kr8*const_pid.OR.start_angle*angfac>2.000001_kr8*const_pid) then
     call log_error(m_name,s_name,20,error_fatal,'start angle not in range -360 to +360')
  end if
  if(finish_angle*angfac<-2.000001_kr8*const_pid.OR.finish_angle*angfac>2.000001_kr8*const_pid) then
     call log_error(m_name,s_name,21,error_fatal,'finish angle not in range -360 to +360')
  end if
  !! set up factor for lengths
  if (length_units(1:2)=='me') then
     lenfac=1000.
     numerics%cunits=0
  else
     numerics%cunits=-3
  end if
  ! positive integer parameter
  if(number_of_divisions<=0) then
     call log_error(m_name,s_name,30,error_fatal,'number of divisions must be positive')
  end if
  if(coordinate_system/=1) then
     call log_error(m_name,s_name,31,error_fatal,'coordinate system must be positive, unity for polars')
  end if

  numerics%descode=igcode
  numerics%tfm=transform_type
  numerics%stang=start_angle*angfac
  numerics%finang=finish_angle*angfac
  numerics%stpos=start_position*lenfac
  numerics%finpos=finish_position*lenfac
  numerics%div=2*((number_of_divisions+1)/2)
  numerics%csys=coordinate_system

  if (.NOT.filedata) then
     ! line defined by namelist, check here and set up
     if(line_divisions<=0) then
        call log_error(m_name,s_name,32,error_fatal,'number of line divisions must be positive')
     end if
     call misc_line2d(numerics%r,numerics%z,numerics%stpos,numerics%finpos,line_divisions)
  end if

  numerics%ldiv=line_divisions
  !! next will be overwritten when there is dat file input
  numerics%npos=line_divisions+1
  klfile=filedata

end  subroutine dcontrol_readnum
!---------------------------------------------------------------------
!> read progfiles namelist data
subroutine dcontrol_readprogfiles(file,kunit)

  !! arguments
  type(dfiles_t), intent(out) :: file !< file names
  integer, intent(in)  :: kunit     !< file unit number

  !!local
  character(*), parameter :: s_name='dcontrol_readprogfiles' !< subroutine name
  logical :: filefound !< true if file exists

  character(len=80) :: r_input_file  !< position coordinate 1 data input file
  character(len=80) :: z_input_file !< position coordinate 2 data input file
  character(len=80) :: rz_input_file  !< position coordinate 1 & 2 data input file

  !! file names
  namelist /progfiles/ &
 &r_input_file, &
 &z_input_file, &
 &rz_input_file

  !---------------------------------------------------------------------
  !! read input file names
  r_input_file='null'
  z_input_file='null'
  rz_input_file='null'
  read(kunit,nml=progfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_getunit(ilog)
     write(ilog,nml=progfiles)
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
  end if

  file%rfile = r_input_file
  file%zfile = z_input_file
  file%rzfile = rz_input_file

  !!check files exist

  call log_value("datvtk data file, r_input_file",trim(file%rfile))
  if(file%rfile/='null') then
     inquire(file=r_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find datvtk data file, ",a)',r_input_file
        call log_error(m_name,s_name,2,error_fatal,'datvtk data file not found')
     end if
  end if
  call log_value("datvtk data file, z_input_file",trim(file%zfile))
  if(file%zfile/='null') then
     inquire(file=z_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find datvtk data file, ",a)',z_input_file
        call log_error(m_name,s_name,3,error_fatal,'datvtk data file not found')
     end if
  end if
  if(file%rzfile/='null') then
     inquire(file=rz_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find datvtk data file, ",a)',rz_input_file
        call log_error(m_name,s_name,4,error_fatal,'datvtk data file not found')
     end if
  end if

end  subroutine dcontrol_readprogfiles
!---------------------------------------------------------------------
!> read dnumerics datafile
subroutine dcontrol_readatfile(file,numerics)

  !! arguments
  type(dfiles_t), intent(inout) :: file !< file names
  type(dnumerics_t), intent(inout) :: numerics !< control numerics

  !!local
  character(*), parameter :: s_name='dcontrol_readatfile' !< subroutine name
  integer  :: ninrz     !< file unit number

  if (file%rzfile=='null') then
     call log_value("r input file",trim(file%rfile))
     call misc_getfileunit(ninrz)
     open(unit=ninrz,file=file%rfile,status='OLD',iostat=status)
     call log_open_check(m_name,s_name,10,status)

     i=0
     do
        read(ninrz,*,iostat=status) zdum
        if(status/=0) exit
        i=i+1
     end do
     if (i<=1) call log_error(m_name,s_name,1,error_fatal,'Insufficient data in file')
     allocate(numerics%r(i),numerics%z(i),stat=status)
     call log_alloc_check(m_name,s_name,11,status)
     numerics%npos=i

     !! read into store values
     rewind(ninrz,iostat=status)
     call log_read_check(m_name,s_name,12,status)
     do j=1,i
        read(ninrz,*,iostat=status) numerics%r(j)
        call log_read_check(m_name,s_name,13,status)
     end do
     call log_value("number of values read from r input file",i)

     close(unit=ninrz)

     call log_value("z input file",trim(file%zfile))
     call misc_getfileunit(ninrz)
     open(unit=ninrz,file=file%zfile,status='OLD',iostat=status)
     call log_open_check(m_name,s_name,14,status)

     i=0
     !! read into store values
     do j=1,numerics%npos
        i=i+1
        read(ninrz,*,iostat=status) numerics%z(j)
        call log_read_check(m_name,s_name,16,status)
     end do
     call log_value("number of values read from z input file",i)

  else

     call log_value("rz input file",trim(file%rzfile))
     call misc_getfileunit(ninrz)
     open(unit=ninrz,file=file%rzfile,status='OLD',iostat=status)
     call log_open_check(m_name,s_name,20,status)

     i=0
     do
        read(ninrz,*,iostat=status) zdum, zdum
        if(status/=0) exit
        i=i+1
     end do
     if (i<=1) call log_error(m_name,s_name,2,error_fatal,'Insufficient data in file')
     allocate(numerics%r(i),numerics%z(i),stat=status)
     call log_alloc_check(m_name,s_name,21,status)
     numerics%npos=i

     !! read into store values
     rewind(ninrz,iostat=status)
     call log_read_check(m_name,s_name,22,status)
     do j=1,i
        read(ninrz,*,iostat=status) numerics%r(j),numerics%z(j)
        call log_read_check(m_name,s_name,23,status)
     end do
     call log_value("number of values read from rz input file",i)

  end if

  close(unit=ninrz)

end  subroutine dcontrol_readatfile
!---------------------------------------------------------------------
!> wrapper for misc_lines2d
subroutine dcontrol_lines2d(numerics)
  !! arguments
  type(dnumerics_t), intent(inout) :: numerics !< module object
  !! local
  character(*), parameter :: s_name='dcontrol_lines2d' !< subroutine name

  call misc_lines2d(numerics%r,numerics%z,numerics%npos,numerics%ldiv)

end subroutine dcontrol_lines2d
!---------------------------------------------------------------------
!> get unit number of input
subroutine dcontrol_getunit(kunit)

  !! arguments
  integer, intent(out) :: kunit    !< log unit number

  kunit=nin

end subroutine dcontrol_getunit
!---------------------------------------------------------------------
!> delete object
subroutine dcontrol_delete(numerics)

  !! arguments
  type(dnumerics_t), intent(inout) :: numerics !< module object
  !! local
  character(*), parameter :: s_name='dcontrol_delete' !< subroutine name

  deallocate(numerics%r,numerics%z)

end subroutine dcontrol_delete
!---------------------------------------------------------------------
end module dcontrol_m
