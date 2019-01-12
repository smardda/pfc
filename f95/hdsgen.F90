program hdsgen_p

  use const_kind_m
  use const_numphys_h
  use control_h
  use dcontrol_h
  use vcontrol_h
  use date_time_m
  use control_m
  use log_m
  use misc_m
  use clock_m
  use position_m
  use geobj_m
  use posang_h
  use ls_m
  use btree_m
  use li_m
  use ld_m
  use dbtree_h
  use dbtree_m
  use termplane_h
  use pcle_h
  use pcle_m
  use bods_h
  use geobjlist_h
  use geobjlist_m
  use indict_m
  use vfile_m
  use outfile_m
  use hdsfile_m
  use query_m
  use position_h
  use fmesh_h
  use beq_h
  use posang_m
  use datline_h
  use datline_m
  use stack_m
  use spl2d_m
  use spl3d_m

  implicit none


! Local variables
  character(*), parameter :: m_name='hdsgen' !< module name
  type(files_t)     :: file      !< names of files
  type(numerics_t)  :: numerics  !< numerical parameters
  type(plots_t)     :: plot      !< diagnostic plot selectors
  type(btree_t)  :: btree      !< binary tree
  type(geobjlist_t) :: geobjl !< geometrical objects
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=256) :: vtkdesc !< descriptor line for vtk files

  integer:: nplot !< unit for vfiles
  integer:: nouth !< unit for hdsfiles

!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(20)
  call clock_start(1,'hdsgen run time')

!--------------------------------------------------------------------------
!! print header

  print *, '--------------------------------------------------'
  print *, 'hdsgen: map data to hybrid data structure'
  print *, '--------------------------------------------------'
  print '(a)', timestamp%long

!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run hdsgen type at the command line:'
     print *, '   hdsgen fileroot'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
  end if

!! start log
  call log_init(fileroot,timestamp)

!--------------------------------------------------------------------------
!! read control file

  call clock_start(2,'control_init time')
  call control_init(fileroot)
  call control_read(file,numerics,plot)
  call clock_stop(2)

!--------------------------------------------------------------------------
!! initialise  geobjl data

  call clock_start(3,'geobjlist_init time')
  call geobjlist_init(geobjl,file%vtkdata,numerics)
  call clock_stop(3)

!--------------------------------------------------------------------------
!! initialise btree

  call clock_start(4,'btree_init time')
  call btree_init(btree,numerics)
  call clock_stop(4)

!--------------------------------------------------------------------------
!! initial btree data diagnostics

! NOT NEEDED

!--------------------------------------------------------------------------
!!sort geobjl onto bins

  call clock_start(9,'geobjlist_sort time')
  call geobjlist_bin(geobjl,btree)
  numerics%geobj_coord_tfm=geobjl%quantfm
  call clock_stop(9)

!--------------------------------------------------------------------------
!! sorted btree and geobjl data diagnostics

  call btree_dia(btree)

!!plot mapped HDS
  if(plot%hdsm) then
     call clock_start(10,'vfile_hdsm time')
     call geobjlist_makehedline(geobjl,'hds lowest',vtkdesc)
     call vfile_init(file%hdsm,vtkdesc,nplot)
     call btree_writev(btree,numerics,'hds lowest',nplot)
     call vfile_close
     call clock_stop(10)
  end if

!!plot HDS quantised (default)
  if(plot%hdsq) then
     call clock_start(11,'vfile_hdsq time')
     call geobjlist_makehedline(geobjl,'hds quantised',vtkdesc)
     call vfile_init(file%hdsq,vtkdesc,nplot)
     call btree_writev(btree,numerics,'hds quantised',nplot)
     call vfile_close
     call clock_stop(11)
  end if

!!plot assigned geobj
  if(plot%geobjq) then
     call clock_start(14,'vfile_assignedgeobj time')
     call geobjlist_makehedline(geobjl,'assigned geobj',vtkdesc)
     call vfile_init(file%geobjq,vtkdesc,nplot)
! in quantised space
     call geobjlist_writev(geobjl,'geometry',nplot)
     call vfile_close
     call clock_stop(14)
  end if

!!plot unassigned geobjs
  if(plot%lostgeobj) then
     call clock_start(15,'vfile_unassignedgeobj time')
     call geobjlist_makehedline(geobjl,'unassigned geobj',vtkdesc)
     call vfile_init(file%lostgeobj,vtkdesc,nplot)
! in quantised space
     call geobjlist_writev(geobjl,'unassigned quantised',nplot)
     call vfile_close
     call clock_stop(15)
  end if

!!plot all geobj points
  if(plot%geoptq) then
     call clock_start(16,'vfile_allgeoptq time')
     call geobjlist_makehedline(geobjl,'all geoptq',vtkdesc)
     call vfile_init(file%geoptq,vtkdesc,nplot)
! in quantised space
     call geobjlist_writev(geobjl,'all quantised',nplot)
     call vfile_close
     call clock_stop(16)
  end if

!!plot density geobjs
  if(plot%densitygeobj) then
     call clock_start(17,'vfile_densitygeobj time')
     call geobjlist_makehedline(geobjl,'density geobj',vtkdesc)
     call vfile_init(file%densitygeobj,vtkdesc,nplot)
! in quantised space
     call geobjlist_writev(geobjl,'density quantised',nplot)
     call vfile_close
     call clock_stop(17)
  end if

!--------------------------------------------------------------------------
!! output file

  call clock_start(20,'outfile_init time')
  call outfile_init(file,timestamp)
  call outfile_write(geobjl,btree)
  call outfile_close
  call clock_stop(20)

!--------------------------------------------------------------------------
!! hds output file

  call clock_start(21,'hdsfile_init time')
  call hdsfile_init(file%hdslist,timestamp,nouth)
  call hdsfile_write(geobjl,numerics,btree)
  call hdsfile_close
  call clock_stop(21)

!--------------------------------------------------------------------------
!! cleanup and closedown
  call geobjlist_delete(geobjl)
  call btree_delete(btree)

  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete

!--------------------------------------------------------------------------

end program hdsgen_p
