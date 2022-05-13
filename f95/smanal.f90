!> @addtogroup groupname111
!> @{
!> @defgroup groupname9 SMARDDA-PFC-NEW
!> @{
program smanal_p
!> @}
!> @}
#ifdef WITH_MPI
  use mpi
#endif

  use const_kind_m
  use const_numphys_h
  use scontrol_h
  use scontrol_m
  use date_time_m
  use log_m
  use misc_m
  use clock_m
  use position_m
  use geobj_m
  use posang_h
!  use pcle_h
!  use pcle_m
  use bods_h
  use geobjlist_h
  use control_h
  use li_m
  use ld_m
  use ls_m
  use dbtree_h
  use dbtree_m
  use btree_m
  use stack_m
  use query_m
  use spl2d_m
  use spl3d_m
  use geobjlist_m
  use dcontrol_h
  use vcontrol_h
  use vfile_m
  use gfile_m
  use outfile_m
  use position_h
  use fmesh_h
  use beq_h
  use posang_m
  use datline_h
  use datline_m
  use smanal_h
  use smanal_m
  use soutfile_m
  use indict_m

  implicit none


! Local variables
  character(*), parameter :: m_name='smanal' !< module name
  type(sfiles_t)     :: file      !< names of files
  type(snumerics_t)  :: numerics  !< numerical parameters
  type(splots_t)     :: plot      !< diagnostic plot selectors
  type(smanal_t) :: smanal !< smanal object
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=80) :: cmesg !< message for output files
  character(len=256) :: vtkdesc !< descriptor line for vtk files

  integer(ki4):: iopt !< option
  integer:: nin !< unit for vtk files
  integer:: nplot !< unit for vtk file output
  integer:: nprint !< unit for gnu file output

#ifdef WITH_MPI
  integer error, rank
  call MPI_Init ( error )
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, error)
  if (rank .eq. 0) then
#endif
!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(20)
  call clock_start(1,'smanal run time')

!--------------------------------------------------------------------------
!! print header

  print *, '--------------------------------------------------'
  print *, 'smanal: analyse smardda output '
  print *, '--------------------------------------------------'
  print '(a)', timestamp%long

!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run smanal type at the command line:'
     print *, '   smanal fileroot'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
  end if

!! start log
  call log_init(fileroot,timestamp)

!--------------------------------------------------------------------------
!! read control file

  call clock_start(2,'scontrol_init time')
  call scontrol_init(fileroot)
  call scontrol_read(numerics,file,plot)
  call clock_stop(2)

!--------------------------------------------------------------------------
!! initialise geobjl data from powcal output file

  call clock_start(3,'geobjlist_init time')
  call smanal_read(smanal,file,numerics)
  call clock_stop(3)

!--------------------------------------------------------------------------
!! main work

  call clock_start(5,'analysis time')
  call smanal_rekey(smanal)
  call smanal_stats(smanal)
  call clock_stop(5)

!--------------------------------------------------------------------------
!! sorted gnuplot and geobjl data diagnostics


!! csv output of body statistics
  if(plot%gnu) then
     call clock_start(25,'vfile_gnuanal time')
     cmesg='Statistics of scalar field '//numerics%namescal
     call gfile_init(file%gnu,cmesg,nprint)
     call smanal_writeg(smanal,'gnuanal',nprint)
     call gfile_close
     if (smanal%n%rekey) then
        call gfile_init(file%gnure,cmesg,nprint)
        call smanal_writeg(smanal,'gnureanal',nprint)
     end if
     call clock_stop(25)
  end if

!! vtk output of body statistics (equivalent to a smanal_writev)
  if(plot%vtk) then
     call clock_start(11,'vfile_analysis time')
     cmesg='Statistics of '//numerics%namescal
     call geobjlist_makehedline(smanal%geobjl,cmesg,vtkdesc)
     call vfile_init(file%vtkfull,vtkdesc,nplot)
     call smanal_writev(smanal,'full',nplot)
     call clock_stop(11)
  end if

!! vtk output of body statistics as above but now small geobjl,
!! corresponding to one point per body
  if(plot%vtksmall) then
     call clock_start(12,'vfile_analysis time')
     cmesg='Small statistics of '//numerics%namescal
     call geobjlist_makehedline(smanal%geobjl,cmesg,vtkdesc)
     call vfile_init(file%vtksmall,vtkdesc,nplot)
     call smanal_writev(smanal,'small',nplot)
     call clock_stop(12)
  end if

!!--------------------------------------------------------------------------
!!! output file
!
  call clock_start(20,'outfile_init time')
!! just dictionaries for now
  call soutfile_init(file,timestamp)
  call soutfile_write(smanal)
  call soutfile_close
  call clock_stop(20)

!--------------------------------------------------------------------------
!! cleanup and closedown
  call smanal_delete(smanal)
  call clock_stop(1)
  call clock_summary

  call scontrol_close
  call log_close
  call clock_delete

!--------------------------------------------------------------------------

#ifdef WITH_MPI
  end if
  call MPI_Finalize ( error )
#endif
end program smanal_p
