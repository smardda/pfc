module bcontrol_m

  use const_kind_m
  use log_m
  use misc_m
  use position_h
  use skyl_h
  use dcontrol_h
  use dcontrol_m
  use skyl_m
  use fmesh_h
  use beq_h
  use beq_m

  implicit none
  private

! public subroutines
  public :: &
 &bcontrol_init, & !< open input control data file
 &bcontrol_closex, & !< close input control data file and flag
 &bcontrol_getunit,  & !< get unit number
 &bcontrol_read


! public types
  type, public :: bfiles_t
   !! file names
     character(len=80)  :: geoqout       !< output data
     character(len=80)  :: log            !< log file
     character(len=80)  :: eqdsk         !< equilibrium type eqtype file name
     character(len=80)  :: equil         !< generic equilibrium type file name
     character(len=80)  :: vtkdata         !< vtk input data file
     character(len=80)  :: laudata         !< launch data
     character(len=80)  :: cartv   !< DUPLICATE  \f$ \bf{B} \f$ in Cartesians
     character(len=80)  :: geoqx   !< \f$ \bf{B} \f$ in Cartesians
     character(len=80)  :: allcartv   !< DUPLICATE  all \f$ \bf{B} \f$ in Cartesians
     character(len=80)  :: geoqvolx   !< all \f$ \bf{B} \f$ in Cartesians
     character(len=80)  :: geofldxyz   !< INERT geometry in transformed Cartesians, field in geometry coordinates
     character(len=80)  :: geoqvolxyz   !< geometry and field in transformed Cartesians on Cartesian grid
     character(len=80)  :: fptv    !< DUPLICATE  vtk plot file of functions of \f$ (\psi,\theta) \f$
     character(len=80)  :: geoqvolm   !< all \f$ \bf{B} \f$ in mapped coordinates and components
     character(len=80)  :: geoqmm    !< DUPLICATE vtk plot file of functions of \f$ (\psi,\theta) \f$
     character(len=80)  :: allgeoq  !< DUPLICATE  plot file of all vtk info needed for SMITER-GEOQ
     character(len=80)  :: geoqm  !< plot file of all vtk info needed for SMITER-GEOQ
     character(len=80)  :: frzv    !< INERT vtk plot file of functions of \f$ (R,Z) \f$
     character(len=80)  :: frzzeta    !< INERT vtk plot file of functions of \f$ (R,Z,\zeta) \f$
     character(len=80)  :: frzxi    !< DUPLICATE  vtk plot file of functions of \f$ (R,Z,\xi) \f$
     character(len=80)  :: geofld  !< DUPLICATE  vtk plot file of geometry and field
     character(len=80)  :: geofldx  !< vtk plot file of geometry and field
     character(len=80)  :: geofldq  !< INERT vtk plot file of quantised geometry and field
     character(len=80)  :: gnu !< gnuplot file of \f$ \psi(R,Z) \f$
     character(len=80)  :: gnusil !< gnuplot file of silhouette as function of \f$ (R,Z) \f$
     character(len=80)  :: gnusilm !< gnuplot file of silhouette as function of \f$ (\psi,\theta) \f$
     character(len=80)  :: geoqfldxyz   !< \f$ \bf{B} \f$ in Cartesians on Cartesian grid
     character(len=80)  :: eqbdry   !< file containing boundary points from eqdsk
     character(len=80)  :: eqltr   !< file containing limiter points from eqdsk
     character(len=80)  :: fmesh   !< define special mesh on input file name
  end type bfiles_t


  type, public :: bplots_t
   !! vtk plot output selectors
     logical  :: cartv   !< DUPLICATE  \f$ \bf{B} \f$ in Cartesians
     logical  :: geoqx   !< \f$ \bf{B} \f$ in Cartesians
     logical  :: allcartv   !< DUPLICATE  all \f$ \bf{B} \f$ in Cartesians
     logical  :: geoqvolx   !< all \f$ \bf{B} \f$ in Cartesians
     logical  :: geoqfldxyz   !< \f$ \bf{B} \f$ in Cartesians on Cartesian grid
     logical  :: fptv    !< DUPLICATE  vtk plot file of functions of \f$ (\psi,\theta) \f$
     logical  :: geoqmm    !< DUPLICATE vtk plot file of functions of \f$ (\psi,\theta) \f$
     logical  :: geoqvolm    !< all \f$ \bf{B} \f$ in mapped coordinates and components
     logical  :: allgeoq    !< DUPLICATE  plot file of all vtk info needed for SMITER-GEOQ
     logical  :: geoqm    !< plot file of all vtk info needed for SMITER-GEOQ
     logical  :: frzv    !< INERT vtk plot file of functions of \f$ (R,Z) \f$
     logical  :: frzzeta    !< INERT vtk plot file of functions of \f$ (R,Z,\zeta) \f$
     logical  :: frzxi    !< DUPLICATE  vtk plot file of functions of \f$ (R,Z,\xi) \f$
     logical  :: geofld  !< DUPLICATE  vtk plot file of geometry and field
     logical  :: geofldx  !< vtk plot file of geometry and field
     logical  :: geofldq  !< INERT vtk plot file of quantised geometry and field
     logical  :: gnu !< gnuplot file of \f$ \psi(R,Z) \f$
     logical  :: gnusil !< gnuplot file of silhouette as function of \f$ (R,Z) \f$
     logical  :: gnusilm !< gnuplot file of silhouette as function of \f$ (\psi,\theta) \f$
     logical  :: gnuptz !< local variable
     logical  :: gnum !< gnuplot file of \f$ fns(\psi,\theta) \f$
     logical  :: geofldxyz   !< geometry and field in Cartesians
     logical  :: geoqvolxyz   !< geometry and field in Cartesians on Cartesian grid
     logical  :: eqbdry   !< produce files containing boundary and limiter points from eqdsk
  end type bplots_t


! private variables
  character(*), parameter :: m_name='bcontrol_m' !< module name
  integer  :: status   !< error status
  integer, save  :: nin=-1      !< control file unit number
  integer  :: ilog      !< for namelist dump after error
  logical :: iltest !< logical flag
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  character(len=80) :: root !< file root

  contains
!---------------------------------------------------------------------
!> open input control data file
subroutine bcontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='bcontrol_init' !< subroutine name
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

end  subroutine bcontrol_init
!---------------------------------------------------------------------
!> close input control data file and flag
subroutine bcontrol_closex

  !! local
  character(*), parameter :: s_name='bcontrol_closex' !< subroutine name

  !! close file unit
  close(unit=nin,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close file unit, ",i5)',nin
     call log_error(m_name,s_name,1,error_fatal,'Cannot close data file')
     stop
  end if
  nin=-2

end  subroutine bcontrol_closex
!---------------------------------------------------------------------
!> get unit number of input
subroutine bcontrol_getunit(kunit)

  !! arguments
  integer, intent(out) :: kunit    !< log unit number

  kunit=nin

end subroutine bcontrol_getunit
!---------------------------------------------------------------------
!> read data for this run
subroutine bcontrol_read(file,numerics,sknumerics,plot)

  !! arguments
  type(bfiles_t), intent(out) :: file !< file names
  type(bnumerics_t), intent(out) :: numerics !< input numerical parameters
  type(sknumerics_t), intent(out) :: sknumerics !< input numerical parameters
  type(bplots_t), intent(out) :: plot !< vtk plot selectors

  !!local
  character(*), parameter :: s_name='bcontrol_read' !< subroutine name
  character(len=80) :: eqdsk_input_file  !< EQDSK G format data input file
  character(len=80) :: equil_input_file  !< generic equilibrium input file
  character(len=80) :: mesh_input_file  !< generic field mesh input file
  character(len=80) :: vtk_input_file  !< vtk format geometry data input file
  character(len=80) :: lau_input_file  !< launch data input file
  character(len=80) :: icsuf  !< eqdsk file suffix
  logical :: filefound !< true of file exists
  integer(ki4) :: dummy_number !< dummy
  integer(ki4) :: indot !< position of dot in filename
  integer(ki4) :: ilent !< length of suffix string
  logical :: plot_cartv !< DUPLICATE  vtk plot selector
  logical :: plot_geoqx !< vtk plot selector
  logical :: plot_allcartv !< DUPLICATE  vtk plot selector
  logical :: plot_geoqvolx !< vtk plot selector
  logical :: plot_geofldxyz !< vtk plot selector
  logical :: plot_geoqvolxyz !< vtk plot selector
  logical :: plot_fptv !< DUPLICATE  vtk plot selector
  logical :: plot_geoqmm !< vtk plot selector
  logical :: plot_geoqvolm !< vtk plot selector
  logical :: plot_allgeoq !< DUPLICATE  vtk plot selector
  logical :: plot_geoqm !< vtk plot selector
  logical :: plot_frzv !< INERT vtk plot selector
  logical :: plot_frzzeta !< INERT vtk plot selector
  logical :: plot_frzxi !< DUPLICATE  vtk plot selector
  logical :: plot_geofld !< DUPLICATE  vtk plot selector
  logical :: plot_geofldx !< vtk plot selector
  logical :: plot_geofld_quantised !< vtk plot selector
  logical :: plot_gnu !< gnu plot selector
  logical :: plot_gnuptz !< local variable
  logical :: plot_gnum !< gnum plot selector
  logical :: plot_gnusil !< gnusil plot selector
  logical :: plot_gnusilm !< gnusilm plot selector
  logical :: plot_eqdsk_boundary !< eqbdry plot selector
  logical :: plot_geoqfldxyz !< save field in special format

  !! file names
  namelist /inputfiles/ &
 &eqdsk_input_file,vtk_input_file,lau_input_file, &
 &equil_input_file, &
 &mesh_input_file

  !! misc parameters
  namelist /miscparameters/ &
 &dummy_number

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_cartv, &
 &plot_geoqx, &
 &plot_allcartv, &
 &plot_geoqvolx, &
 &plot_geofldxyz, &
 &plot_geoqvolxyz, &
 &plot_fptv, &
 &plot_geoqmm, &
 &plot_geoqvolm, &
 &plot_allgeoq, &
 &plot_geoqm, &
 &plot_frzv, &
 &plot_frzzeta, &
 &plot_frzxi, &
 &plot_geofld, &
 &plot_geofldx, &
 &plot_geofld_quantised, &
 &plot_gnu, &
 &plot_gnuptz, &
 &plot_gnum, &
 &plot_gnusil, &
 &plot_gnusilm, &
 &plot_eqdsk_boundary, &
 &plot_geoqfldxyz

  !! read input file names
  vtk_input_file='null'
  lau_input_file='null'
  equil_input_file='null'
  mesh_input_file='null'
  eqdsk_input_file='null'
  read(nin,nml=inputfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_getunit(ilog)
     write(ilog,nml=inputfiles)
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
  end if

  if ( equil_input_file=='null') then
     equil_input_file   = eqdsk_input_file
  end if
  file%eqdsk   = equil_input_file
  file%equil   = equil_input_file
  file%fmesh   = mesh_input_file
  file%vtkdata   = vtk_input_file
  file%laudata   = lau_input_file

  !!check files exist

  call log_value("beq field data file, equil_input_file",trim(file%equil))
  inquire(file=equil_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Error : Unable to find Beq field data file, ",a)',equil_input_file
     call log_error(m_name,s_name,2,error_warning,'Beq field data file not found')
  end if
  !  call log_value("beq field data file, mesh_input_file",trim(file%fmesh))
  !  inquire(file=mesh_input_file,exist=filefound)
  !  if(.not.filefound) then
  !     !! error opening file
  !     print '("Fatal error: Unable to find Mesh field data file, ",a)',mesh_input_file
  !     call log_error(m_name,s_name,2,error_fatal,'Mesh field data file not found')
  !  end if

  ! set equilibrium type depending on suffix
  indot=index(equil_input_file,'.',.TRUE.)
  if (indot<2.OR.indot>80) then
     call log_error(m_name,s_name,2,error_fatal,'Beq field data file has no suffix')
  end if
  icsuf=equil_input_file(indot+1:)
  ilent=len_trim(icsuf)
  call lowor(icsuf,1,ilent)
  numerics%eqtype=icsuf

  call log_value("geometry data file, vtk_input_file",trim(file%vtkdata))
  if(file%vtkdata/='null') then
     inquire(file=vtk_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find Beq geometry data file, ",a)',vtk_input_file
        call log_error(m_name,s_name,3,error_fatal,'Beq geometry data file not found')
     end if
  end if

  if(file%laudata/='null') then
     inquire(file=lau_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find launch data file, ",a)',lau_input_file
        call log_error(m_name,s_name,4,error_fatal,'launch data file not found')
     end if
  end if

  if ( mesh_input_file/='null') then
     call log_error(m_name,s_name,22,error_warning,'Mesh input filename is not null')
     inquire(file=mesh_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find field mesh data file, ",a)',mesh_input_file
        call log_error(m_name,s_name,22,error_fatal,'Mesh field data file not found')
     end if
  end if

  !! create output file names from root

  !! output file
  file%geoqout = trim(root)//"_geoq.out"
  inquire(file=file%geoqout,exist=filefound)
  call log_value("geoq output file",trim(file%geoqout))
  if(filefound) then
     !! error opening file
     print '("Fatal error: Output file ",a, " already exists")',trim(file%geoqout)
     print '("Remove it and restart run")'
     call log_error(m_name,s_name,5,error_fatal,'Output file already exists')
  end if

  !!vtk file roots
  file%cartv     =trim(root)//"_cartv"
  file%geoqx     =trim(root)//"_geoqx"
  file%allcartv     =trim(root)//"_allcartv"
  file%geoqvolx     =trim(root)//"_geoqvolx"
  file%geofldxyz     =trim(root)//"_geofldxyz"
  file%geoqvolxyz     =trim(root)//"_geoqvolxyz"
  file%fptv    =trim(root)//"_fptv"
  file%geoqmm    =trim(root)//"_geoqmm"
  file%geoqvolm    =trim(root)//"_geoqvolm"
  file%allgeoq    =trim(root)//"_allgeoq"
  file%geoqm    =trim(root)//"_geoqm"
  file%frzv    =trim(root)//"_frzv"
  file%frzzeta    =trim(root)//"_frzzeta"
  file%frzxi    =trim(root)//"_frzxi"
  file%geofld =trim(root)//"_geofld"
  file%geofldx =trim(root)//"_geofldx"
  file%geofldq    =trim(root)//"_geofldq"
  !!gnu file roots
  file%gnu     =trim(root)//"_gnu"
  file%gnusil     =trim(root)//"_gnusil"
  file%gnusilm     =trim(root)//"_gnusilm"
  file%eqbdry     =file%equil(1:indot-1)//"_eqbdry"
  file%eqltr     =file%equil(1:indot-1)//"_eqltr"
  !! special field file format
  file%geoqfldxyz     =trim(root)//"_geoqfldxyz"
  file%geoqfldxyz     =trim(root)//"_fieldi"

  !! set default misc parameters
  dummy_number = 10

  !!read misc parameters
  read(nin,nml=miscparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading misc parameters")'
     call log_getunit(ilog)
     write(ilog,nml=miscparameters)
     call log_error(m_name,s_name,10,error_fatal,'Error reading misc parameters')
  end if

  !! check for valid data
  if(dummy_number<0) &
 &call log_error(m_name,s_name,11,error_fatal,'dummy_number must be >=0')

  !W     numerics%ndummy = dummy_number

  !! set default plot selections
  plot_cartv = .false.
  plot_geoqx = .false.
  plot_allcartv = .false.
  plot_geoqvolx = .false.
  plot_geofldxyz = .false.
  plot_geoqvolxyz = .false.
  plot_fptv = .false.
  plot_geoqmm = .false.
  plot_geoqvolm = .false.
  plot_allgeoq = .false.
  plot_geoqm = .false.
  plot_frzv = .false.
  plot_frzzeta = .false.
  plot_frzxi = .false.
  plot_geofld = .false.
  plot_geofldx = .false.
  plot_geofld_quantised = .false.
  plot_gnu = .false.
  plot_gnuptz = .false.
  plot_gnum = .false.
  plot_gnusil = .false.
  plot_gnusilm = .false.
  plot_eqdsk_boundary = .false.
  plot_geoqfldxyz = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     print '("Fatal error reading plot selections")'
     call log_getunit(ilog)
     write(ilog,nml=plotselections)
     call log_error(m_name,s_name,12,error_fatal,'Error reading plot selections')
  end if

  !! store values
  plot%geoqx = plot_cartv
  plot%geoqx   = plot_geoqx
  plot%geoqvolx = plot_allcartv
  plot%geoqvolx   = plot_geoqvolx
  plot%geofldxyz   = plot_geofldxyz
  plot%geoqvolxyz   = plot_geoqvolxyz
  plot%geoqfldxyz   = plot_geoqfldxyz
  plot%geoqmm    = plot_geoqmm.OR.plot_fptv
  plot%geoqvolm    = plot_geoqvolm.OR.plot_geoqmm
  plot%geoqm    = plot_geoqm
  plot%frzv     = plot_frzv
  plot%frzzeta     = plot_frzzeta
  plot%geofldx = plot_geofld
  plot%geofldx    = plot_geofldx
  plot%geofldq     = plot_geofld_quantised
  plot%gnu     = plot_gnu
  plot%gnuptz = plot_gnuptz
  plot%gnum     = plot_gnum
  plot%gnusil     = plot_gnusil
  plot%gnusilm     = plot_gnusilm
  plot%eqbdry     = plot_eqdsk_boundary
  if (plot_frzxi) then
     plot%geoqm = plot_frzxi
     file%geoqm = file%frzxi
     call log_error(m_name,s_name,20,error_warning,'Obsolete file handling feature frzxi activated')
  else if (plot_allgeoq) then
     plot%geoqm = plot_allgeoq
     file%geoqm = file%allgeoq
     call log_error(m_name,s_name,21,error_warning,'Obsolete file handling feature allgeoq activated')
  end if
  if (plot_gnuptz) then
     plot%gnum = plot_gnuptz
     call log_error(m_name,s_name,22,error_warning,'Obsolete file handling feature gnuptz activated')
  end if
  if (plot_geoqmm) then
     call log_error(m_name,s_name,23,error_warning,'Obsolete file handling switch geoqmm/fptv used')
  end if

  call beq_readcon(numerics,nin)

  numerics%eqbdry=plot%eqbdry
  numerics%eqbdryfile=file%eqbdry
  numerics%eqltrfile=file%eqltr

  iltest=numerics%skyl
  if (iltest) then
     call skyl_readcon(sknumerics,nin)
     numerics%skyl=(sknumerics%skyltyp>0)
  end if

end  subroutine bcontrol_read

end module bcontrol_m
