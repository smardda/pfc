program p_antest

  use const_kind_m
  use const_numphys_h
  use date_time_m
  use log_h
  use log_m
  use clock_m
  use misc_m
  use gfile_m
  use vfile_m
  use control_h
  use ray_m
  use spl2d_m
  use beq_h
  use beq_m
  use powelt_m
  use powcal_h
  use powres_h
  use skyl_h
  use prcontrol_h
  use geobjlist_h
  use geobjlist_m
  use radnsource_h
  use radnsource_m
  use atcontrol_h
  use atcontrol_m
  use antenna_h
  use antenna_m
  use antest_h
  use antest_m
  use smitermpi_h

  implicit none

! Local variables
  character(*), parameter :: m_name='antest' !< module name
  type(atparams_t)     :: param      !< run control parameters
  type(atfiles_t)     :: file      !< names of files
  type(atplots_t)     :: plot      !< diagnostic plot selectors
  type(antest_t)  :: antest      !< antest object
  type(antenna_t) :: antenna !< antenna data
  type(radnsource_t) :: radnsource !< radnsource data
  type(ray_t) :: ray !< ray data
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  type(numerics_t) :: gnumerics !< numerical control parameters for geobjlist
  type(beq_t)  :: beq      !< equilibrium B
  type(bnumerics_t)  :: numerics  !< numerical control parameters for equilibrium B
  type(fmesh_t)  :: fmesh      !< mesh (for duct problem)
  integer(ki4):: ninb !< unit for beq data
  integer(ki4) :: ifldspec !< field spec
  integer(ki4) :: ninv !< unit for radnsource data
  integer(ki4) :: ninr !< unit for ray data
  integer(ki4):: nin=0 !< unit for other data
  integer(ki4):: nplot=0 !< unit for vtk files
  integer(ki4):: nprint=0 !< unit for gnuplot files
  integer(ki4):: i !< loop variable
  integer(ki4):: j !< loop variable
  real(kr8) :: psi !< local variable
  real(kr8) :: normpsi !< local variable
  character(len=80) :: ibuf !< character workspace
  character(len=256) :: vtkdesc !< descriptor line for vtk files
  logical, parameter :: lpowcaltype=.true. !< version for subsequent powcal integration
  integer(ki4)  :: ishadow   !< if zero, no geometry shadowing
  integer(ki4)  :: status   !< error status

!! hemisphere diagnostics
  character(len=80) :: kchar !< local variable
  integer :: kplot !< local variable
!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(40)
  call clock_start(1,'run time')
!--------------------------------------------------------------------------
!! print header

  print *, '----------------------------------------------------'
  print *, 'antest: antenna test '
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long
!--------------------------------------------------------------------------
!! file root

!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run antest type at the command line:'
     print *, '   antest fileroot'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
  end if

!! start log
  call log_init(fileroot,timestamp)
!--------------------------------------------------------------------------
!! read control file

  call clock_start(2,'control initialisation time')
  call atcontrol_init(fileroot)
  call atcontrol_read(file,param,antest%n,antenna%n,radnsource%n,ray%n,numerics,plot)
  call clock_stop(2)
!--------------------------------------------------------------------------
!! other data for antest and/or read antest controls directly

  call clock_start(15,'initialisation time')
  call clock_stop(15)
!--------------------------------------------------------------------------
!! beq data

  call clock_start(3,'beq_init time')
!READFILE   file%equil='coarse.eqdsk'
  ifldspec=numerics%fldspec
  numerics%eqbdry=.FALSE.
!READFILE   write(*,*) 'file%gnu = ', file%gnu
  select case (numerics%eqtype)
  case ('ana','equ')
     print *, 'Fatal error: equilibrium file extension unsupported.'
     stop
  case ('geqdsk')
     numerics%fiesta=1
  case default
! default flag for external code generating file defining axisymmetric field
     numerics%fiesta=0
  end select
  call beq_readequil(beq,file%equil,numerics)
  call beq_move(beq,numerics)
  call beq_fixorigin(beq,numerics)
  call beq_init(beq,numerics,fmesh)
  call beq_sense(beq,0)
!write(*,*) 'numerics = ',numerics
  call clock_stop(3)

  call clock_start(6,'beq splines evaluation time')

  select case (ifldspec)
  case (2)
     call beq_centre(beq)
     call beq_bdryrb(beq)
     call beq_ipsi(beq)
  case (3)
     call beq_centre(beq)
     call beq_bdryrb(beq)
     call beq_ipsi(beq)
     call beq_rispld(beq)
  case default
     call beq_centre(beq)
     call beq_psilt(beq)
     call beq_thetalt(beq)
     call beq_rextrema(beq)
     call beq_rpsi(beq)
     call beq_rjpsi(beq)
     call beq_bdryrb(beq)
     call beq_ipsi(beq)
  end select
  call clock_stop(6)

!--------------------------------------------------------------------------
!! gnu plot to check input
  if(plot%gnu) then
     call clock_start(18,'gfile_gnu time')
     call gfile_init(trim(file%gnu),'psi sample in R-Z space',nprint)
     call spl2d_writeg(beq%psi,'sampl',nprint)
     call gfile_close
     call gfile_init(trim(file%gnu)//'_dR','dpsi/dR sample in R-Z space',nprint)
     call spl2d_writeg(beq%dpsidr,'sampl',nprint)
     call gfile_close
     call gfile_init(trim(file%gnu)//'_dZ','dpsi/dZ sample in R-Z space',nprint)
     call spl2d_writeg(beq%dpsidz,'sampl',nprint)
     call gfile_close
     if (beq%n%vactype=='nul') then
        call gfile_init(trim(file%gnu)//'_cart','all in cartesian',nprint)
        call beq_writeg(beq,'allcartesian',nprint)
        call gfile_close
     end if
     call clock_stop(18)
  end if
!--------------------------------------------------------------------------
!! initialise  geobjl data
  call clock_start(4,'geobjlist_init time')
!call atcontrol_getunit(ning)
!READFILE   file%antdata = 'ant11.vtk'
!READFILE   file%vtk = 'jdiv10.vtk'
!READFILE   file%radata = 'kprad_00224_6.txt'
!READFILE   write(*,*) 'file%antdata = ', file%antdata
  call geobjlist_read(antenna%hemisphere,file%antdata)
  print*, "antenna%hemisphere%ng=", antenna%hemisphere%ng
  antest%n%ntheta=antenna%hemisphere%nparam(1)
  antest%n%nphi=antenna%hemisphere%nparam(2)

!READFILE   write(*,*) 'file%vtk = ', file%vtk
  call geobjlist_read(antest%geobjl,file%vtk)
  print*, "antest%geobjl%ng=", antest%geobjl%ng

  call clock_stop(4)
!--------------------------------------------------------------------------
!! initialise  radnsource  data
  call clock_start(7,'radnsource_init time')

  ishadow=0
  call radnsource_init(radnsource, &
 &beq%ivac,beq%psi, &
 &beq%psiaxis,beq%psibdry, &
 &beq%n%rcen,beq%n%zcen, &
 &file%radata)

  call clock_stop(7)
!--------------------------------------------------------------------------
!! gnu plot to check input
  if(plot%gnu) then
     call clock_start(19,'gfile_gnu_radn time')
     call gfile_init(trim(file%gnu)//'_radn','radiation source in R-Z space',nprint)
     call spl2d_writeg(radnsource%srcspl2d,'sampl',nprint)
     call gfile_close
     call clock_stop(19)
  end if
!--------------------------------------------------------------------------
!! do the main work

  call clock_start(21,'evaluation time')

  select case(param%control)
  case default

     call radnsource_parameter(radnsource)
!call antest_solve(antest,beq,antest%geobjl,radnsource,ray)

     if (lpowcaltype) then
        call antenna_init(antenna)
        call ray_init(ray,beq%psi,radnsource%srcspl2d,radnsource%n%psicen,radnsource%n%psibdry)
        allocate(antest%pow(antest%geobjl%ng), stat=status)
        call log_alloc_check(m_name,m_name,15,status)
        call antenna_testelt(antenna,antest%geobjl,beq%psi,radnsource,ray,antest%pow)
     else
        call ray_allocate(ray)
        call antenna_walltest(antenna,antest%geobjl,beq,ray)
     end if

! call radnsource_normtesting(radnsource)

  end select

  call clock_stop(21)
!--------------------------------------------------------------------------
!! data diagnostics to log file - optional

  call clock_start(31,'diagnostics time')

  select case(param%control)
  case default
     call antest_dia(antest)
  end select

  call clock_stop(31)
!--------------------------------------------------------------------------
!! vtk plot file(s) - optional

  call clock_start(32,'vtk plot time')

  select case(param%control)
  case default
     if(plot%powx) then
        call geobjlist_makehedline(antest%geobjl,'radiation power in cartesians',vtkdesc)
        call vfile_init(file%powx,vtkdesc,nplot)
        call antest_writev(antest,'selector',nplot)
        call vfile_close
     end if
  end select

  call clock_stop(32)
!--------------------------------------------------------------------------
!! gnuplot file(s) - optional

  call clock_start(33,'gnuplot time')

  select case(param%control)
  case default
     if(plot%gnu) then
        call gfile_init(trim(file%gnu)//'_test','file header',nprint)
        call antest_writeg(antest,'selector',nprint)
        call gfile_close
     end if
  end select

  call clock_stop(33)
!--------------------------------------------------------------------------
!! output file

  call clock_start(34,'outfile_init time')
  call antest_initwrite(fileroot)
  call antest_write(antest)
  call antest_closewrite()
  call clock_stop(34)
!--------------------------------------------------------------------------
!! cleanup and closedown

  if (lpowcaltype) deallocate(antest%pow)
  call antest_delete(antest)
  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete
!--------------------------------------------------------------------------

end program p_antest
