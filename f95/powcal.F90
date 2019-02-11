program powcal_p

  use const_kind_m
  use const_numphys_h
  use date_time_m
  use log_m
  use misc_m
  use clock_m
  use position_h
  use fmesh_h
  use beq_h
  use posang_h
  use pcontrol_h
  use dcontrol_h
  use pcontrol_m
  use dcontrol_m
  use fmesh_m
  use beq_m
  use beqan_h
  use beqan_m
  use posang_m

  use bods_h
  use geobjlist_h
  use position_m
  use control_h
  use vcontrol_h
  use mcontrol_h
  use mcontrol_m
  use ls_m
  use btree_m
  use li_m
  use ld_m
  use dbtree_h
  use dbtree_m
  use geobj_m
  use query_m
  use datline_h
  use datline_m
  use geobjlist_m
  use indict_m
  use pcle_h
  use pcle_m

  use spl2d_m
  use spl3d_m
  use poutfile_m
  use moutfile_m
  use vfile_m
  use gfile_m
  use hdsfile_m

  use nrsolve_m
  use powelt_m
  use skyl_h
  use powres_h
  use powcal_h
  use termplane_h
  use edgprof_h
  use edgprof_m
  use termplane_m
!     use powres_m
  use powcal_m
  use odes_h
  use odes_m
  use stack_m
  use skyl_m

#ifdef WITH_MPI
  use mpi
#endif
  implicit none

! Local variables
  character(*), parameter :: m_name='powcal' !< module name
  type(pfiles_t)     :: file      !< names of files
  type(numerics_t)  :: gnumerics  !< numerical geobjlist control parameters
  type(pnumerics_t)  :: numerics  !< numerical control parameters
  type(onumerics_t)  :: onumerics  !< numerical ODE control parameters
  type(pplots_t)     :: plot      !< diagnostic plot selectors
  type(btree_t)  :: btree      !< binary tree
  type(powcal_t)  :: powcal      !< power and results geometry data structure
  type(geobjlist_t) :: gshadl !< shadowing geometrical objects
  type(date_time_t) :: timestamp !< timestamp of run
  type(date_time_t) :: timeold !< timestamp of preceding HDSGEN run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=256) :: vtkdesc !< descriptor line for vtk files

  integer:: nplot !< unit for vtk files
  integer:: nin=0 !< unit for other data
  integer(ki4):: ifldspec !< field specification
  real(kr8):: zivac !< value of I in field file
  real(kr8):: zfac !< ratio of I in different field files
  integer(ki4) :: iextra=0 !< extra info in file
!--------------------------------------------------------------------------
#ifdef WITH_MPI
  integer error, rank
  call MPI_Init ( error )
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, error)
  if (rank .eq. 0) then
#else
  integer, parameter :: rank = 0
#endif
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(30)
  call clock_start(1,'powcal run time')
!--------------------------------------------------------------------------
!! print header
  if ( rank .eq. 0) then
     print *, '----------------------------------------------------'
     print *, 'powcal: surface power deposition calculation '
     print *, '----------------------------------------------------'
     print '(a)', timestamp%long
  else
     print *, 'Process rank ', rank, ' started'
  end if
!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run powcal type at the command line:'
     print *, '   powcal fileroot'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
  end if

!! start log
  call log_init(fileroot,timestamp)
!--------------------------------------------------------------------------
!! read control file

  call clock_start(2,'pcontrol_init time')
  call pcontrol_init(fileroot)
  call pcontrol_read(file,numerics,onumerics,powcal%edgprof,plot)
  call clock_stop(2)
!--------------------------------------------------------------------------
!! beq part data output by geoq

  call clock_start(3,'geoq_init time')
! get fldspec from file
  call beq_readcheck(powcal%powres%beq,file%geoq,iextra)
  ifldspec=powcal%powres%beq%n%fldspec
  fld_specn: select case (ifldspec)
  case(1)
     call beq_readpart(powcal%powres%beq,file%geoq)
  case default
     call beq_readplus(powcal%powres%beq,file%geoq)
  end select fld_specn
  if (iextra==2) call skyl_read(powcal%powres%skyl,file%geoq)
  call pcontrol_fix(numerics,onumerics,ifldspec)
  call clock_stop(3)
!--------------------------------------------------------------------------
!! special vacuum field file
  call clock_start(9,'vacfld_init time')
  if (powcal%powres%beq%n%vacfile/='null') then
     call moutfile_read(powcal%powres%beq%n%vacfile,zivac,nin)
     call spl3d_read(powcal%powres%beq%vacfld,powcal%powres%beq%n%vacfile,nin)
     zfac=zivac/powcal%powres%beq%ivac
     call log_value("ripple field i to equilibrium field i ",zfac)
     if (abs(abs(zfac)-1.)>0.01) then
!print *, 'Fatal error: field I and BR values inconsistent.'
        call log_error(m_name,m_name,1,error_fatal,'eqm. currents inconsistent')
     else
        if (abs(zfac-1.)>0.01) then
           call log_error(m_name,m_name,2,error_warning,'eqm. currents have opposite sign')
        end if
     end if
  end if
  call clock_stop(9)

  if (numerics%nanalau==1) then
!--------------------------------------------------------------------------
!! special  launch description
!! read  geobjl data (shadowing geometry data)
!! must have associated X and B data
     call clock_start(5,'shadow geobjlist time')
     call powcal_readv(powcal,file%vtkres)
     gshadl=powcal%powres%geobjl
     call clock_stop(5)
!--------------------------------------------------------------------------
  else
!! launch from `results' geometry
!! read  geobjl data (`results' geometry data)
     call clock_start(4,'results geobjlist time')
     call powcal_readv(powcal,file%vtkres)
     call clock_stop(4)

!! read  geobjl data (shadowing geometry data)
     call clock_start(5,'shadow geobjlist time')
     call geobjlist_read(gshadl,file%vtkdata)
     call clock_stop(5)
  end if
!--------------------------------------------------------------------------
!! read  HDS  data

  call clock_start(6,'hds_init time')
!     line below needed for debugging version to run
!     write(*,*) 'debugging output',file%vtkdata
  call hdsfile_dinit(file%hdsdata,timeold)
  call hdsfile_read(gshadl,gnumerics,btree)
! rotate and quantise shadowing geometry
  call position_tfmlis(gshadl%posl,gshadl%tfmdata)
  call position_qtfmlis(gshadl%posl,gshadl%quantfm)
  call btree_limits(btree,powcal%powres%beq%n%xbb)
  call clock_stop(6)
!--------------------------------------------------------------------------
!! initialise powcal and odes

  call clock_start(7,'initialise powcal time')
  call odes_init(powcal%odes,onumerics)
  call powcal_init(powcal,numerics,plot,gnumerics)
  if (powcal%powres%flinends) then
     call gfile_init(trim(file%flinends),'ends',powcal%powres%nflends)
  end if

  call clock_stop(7)
!--------------------------------------------------------------------------
!! do the main work

  call clock_start(8,'powcal evaluation time')
  call powcal_refine(powcal,gshadl,btree)
  call clock_stop(8)
!--------------------------------------------------------------------------
! field line ends
if (rank .eq. 0) then
  if (powcal%powres%flinends) then
     call gfile_close
  end if

!! power diagnostics

!!plot cartesian power
  if(plot%powstatx) then
     call clock_start(10,'vfile_powstatx time')
     call geobjlist_makehedline(powcal%powres%geobjl,'power in cartesian space',vtkdesc)
     call vfile_init(file%powstatx,vtkdesc,nplot)
     call powcal_writev(powcal,'cartesian',nplot)
     call vfile_close
     call clock_stop(10)
  end if

!!plot cartesian power statistics
  if(plot%powx) then
     call clock_start(11,'vfile_powx time')
     call geobjlist_makehedline(powcal%powres%geobjl,'power statistics in cartesians',vtkdesc)
     call vfile_init(file%powx,vtkdesc,nplot)
     call powcal_writev(powcal,'allcartesian',nplot)
     call vfile_close
     call clock_stop(11)
  end if

!!plot wall power
! only valid for 'msum','middle' option
  if(plot%wall) then
     call clock_start(12,'vfile_wall time')
     call geobjlist_makehedline(gshadl,'wall power in cartesians',vtkdesc)
     call vfile_init(file%wall,vtkdesc,nplot)
     call powcal_xferpow(powcal,gshadl)
     call geobjlist_writev(gshadl,'geometry',nplot)
     call vfile_close
     call clock_stop(12)
  end if
!--------------------------------------------------------------------------
!! output file

  call clock_start(30,'outfile_init time')
  call poutfile_init(file,timestamp)
  call poutfile_write(powcal,timestamp)
  call poutfile_close
  call clock_stop(30)
end if
!--------------------------------------------------------------------------
!! cleanup and closedown
!      if (powcal%powres%n%nanalau==1) then
!      call powcal_delete(powcal)
!      call geobjlist_delete(gshadl)
!      else
  call powcal_delete(powcal)
  call geobjlist_delete(gshadl)
!      end if

  call clock_stop(1)

  if (rank .eq. 0) call clock_summary
  call log_close
  call clock_delete
!--------------------------------------------------------------------------

#ifdef WITH_MPI
  end if
  call MPI_Finalize ( error )
#endif
end program powcal_p
