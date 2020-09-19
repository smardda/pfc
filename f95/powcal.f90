module powcal_p

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
  use beqart_h
  use beqart_m
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

  use smitermpi_m
  use smitermpi_h

  implicit none

contains

  subroutine powcal_sub(arg,combined,decomp)

    ! Arguments
    integer(ki4) :: arg
    logical :: combined

    ! Local variables
    character(*), parameter :: m_name='powcal' !< module name
    type(pfiles_t)     :: file      !< names of files
    type(numerics_t)  :: gnumerics  !< numerical geobjlist control parameters
    type(pnumerics_t)  :: numerics  !< numerical control parameters
    type(onumerics_t)  :: onumerics  !< numerical ODE control parameters
    type(pplots_t)     :: plot      !< diagnostic plot selectors
    type(fmesh_t)  :: fmesh      !< mesh for Cartesian B
    type(btree_t)  :: btree      !< binary tree
    type(powcal_t)  :: powcal      !< power and results geometry data structure
    type(decomp_t)  :: decomp
    type(geobjlist_t) :: gshadl !< shadowing geometrical objects
    type(date_time_t) :: timestamp !< timestamp of run
    type(date_time_t) :: timeold !< timestamp of preceding HDSGEN run
    character(len=80) :: fileroot !< reference name for all files output by run
    character(len=256) :: vtkdesc !< descriptor line for vtk files

    integer(ki4):: nplot !< unit for vtk files
    integer(ki4):: nin=0 !< unit for other data
    integer(ki4):: ifldspec !< field specification
    real(kr8):: zivac !< value of I in field file
    real(kr8):: zfac !< ratio of I in different field files
    integer(ki4) :: iextra=0 !< extra info in file
    integer(ki4) :: ierr

    !--------------------------------------------------------------------------

    if(powcal%decomp%myrankpowcal == 0) call clock_init(36)

    !--------------------------------------------------------------------------

    powcal%decomp=decomp
    powcal%powres%geobjl%ngtype=0

    !! initialise timing
    if(powcal%decomp%myrankpowcal == 0) then
       call date_time_init(timestamp)
       call clock_start(1,'powcal run time')
    endif

    !--------------------------------------------------------------------------
    !! print header
    if(powcal%decomp%myrankpowcal .eq. 0) then
       print *, '----------------------------------------------------'
       print *, 'powcal: surface power deposition calculation '
       print *, '----------------------------------------------------'
       print '(a)', timestamp%long
    endif
    !--------------------------------------------------------------------------
    !! get file root from arg
    if(command_argument_count()<1) then
       !! no file root specified
       if(powcal%decomp%myrankpowcal .eq. 0) then
          print *, 'Fatal error: no file root name specified.'
          print *, 'To run powcal type at the command line:'
          print *, '   powcal fileroot'
       endif
       stop
    else
       !!get fileroot
       call get_command_argument(arg,value=fileroot)
    end if
    !! start log
    call log_init(fileroot,timestamp,decomp)
    !--------------------------------------------------------------------------
    !! read control file
    if(powcal%decomp%readpowcalcontrol) then
       call clock_start(2,'pcontrol_init time')
       call pcontrol_init(fileroot)
       ! Initialise fint to 1.0. It is used in userdefined formula
       powcal%edgprof%fint = 1.0
       call pcontrol_read(file,numerics,onumerics,powcal%edgprof,plot,combined)
       call clock_stop(2)
    endif
    if(powcal%decomp%bcastpowcalcontrol) then
       if(powcal%decomp%myrankpowcal == 0) call clock_start(29,'pcontrol broadcast')
       call bcast_pcontrol(file,numerics,onumerics,powcal%edgprof,&
            plot,powcal%decomp%myrankpowcal,powcal%decomp%commpowcal)
       if(powcal%decomp%myrankpowcal == 0) call clock_stop(29)
    endif
    !--------------------------------------------------------------------------
    !! beq part data output by geoq
    if(powcal%decomp%myrankpowcal == 0) call clock_start(3,'geoq_init time')
    ! Initialize values that are not read in beq_readpart
    powcal%powres%beq%nctrack = 0
    powcal%powres%beq%n%skyl = .false.
    if(powcal%decomp%readpowcaldata) then
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
    endif
    if(powcal%decomp%myrankpowcal == 0) call clock_stop(3)
    if(powcal%decomp%bcastpowcaldata) then
       if(powcal%decomp%myrankpowcal == 0) call clock_start(31,'beq broadcast')
       call bcast_beq(powcal%powres%beq,powcal%decomp%myrankgeoqshad,powcal%decomp%commgeoqshad)
       if(powcal%powres%beq%n%skyl) then
          call bcast_skyl(powcal%powres%skyl,powcal%decomp%myrankgeoqshad,powcal%decomp%commgeoqshad)
       endif
       if(powcal%decomp%myrankpowcal == 0) call clock_stop(31)
    endif

    ifldspec=powcal%powres%beq%n%fldspec
    call pcontrol_fix(numerics,onumerics,ifldspec)

    !--------------------------------------------------------------------------
    !! special vacuum field file
    if(powcal%decomp%myrankpowcal == 0) call clock_start(9,'vacfld_init time')
    vac_type: select case (powcal%powres%beq%n%vactype(1:3))
    case('sp3')

       call moutfile_read(powcal%powres%beq%n%vacfile,zivac,nin)
       call spl3d_read(powcal%powres%beq%vacfld,powcal%powres%beq%n%vacfile,nin)
       zfac=zivac/powcal%powres%beq%ivac
       call log_value("ripple field i to equilibrium field i ",zfac)
       if (abs(abs(zfac)-1.)>0.02) then
          !print *, 'Fatal error: field I and BR values inconsistent.'
          call log_error(m_name,m_name,1,error_fatal,'eqm. currents inconsistent')
       else
          if (abs(zfac-1.)>0.01) then
             call log_error(m_name,m_name,2,error_warning,'eqm. currents have opposite sign')
          end if
       end if
    case('txt')
  !! set up field mesh
       call fmesh_init(powcal%powres%beq%n%vacfile,nin)
  !! read field mesh
       call fmesh_readcon(fmesh,nin)
  !! read field
       call beqart_read(powcal%powres%beq%beqart,powcal%powres%beq%n%vacfile,'txt',nin)
       powcal%powres%beq%beqart%fmesh=fmesh
       call fmesh_close()
  ! input in cartesians convert to cyls
       call beqart_tfm(powcal%powres%beq%beqart,32,0)
       call beqart_sptfm(powcal%powres%beq%beqart,spl3d_rpower)
    end select vac_type
    if(powcal%decomp%myrankpowcal == 0) call clock_stop(9)
    if (numerics%nanalau==1) then
       !--------------------------------------------------------------------------
       !! special  launch description
       !! read  geobjl data (shadowing geometry data)
       !! must have associated X and B data
       if(powcal%decomp%myrankpowcal == 0) call clock_start(5,'shadow geobjlist time')

       call powcal_readv(powcal,file%vtkres)
       gshadl=powcal%powres%geobjl
       if(powcal%decomp%myrankpowcal == 0) call clock_stop(5)
       !--------------------------------------------------------------------------
    else
       !! read  geobjl data (shadowing geometry data)
       gshadl%ngtype=0 !LKUV
       gshadl%ngunassigned=0
       gshadl%coordbb = 0
       ! gshadl%nsampl = 0
       gshadl%binbb = 0
       if(powcal%decomp%readpowcaldata) then
          call clock_start(5,'shadow geobjlist time')
          call geobjlist_read(gshadl,file%vtkdata)
          call clock_stop(5)
       endif
       if(powcal%decomp%bcastpowcaldata) then
          if(powcal%decomp%myrankpowcal == 0) call clock_start(33,'shadow geobjlist broadcast time')
          call bcast_geobjlist(gshadl,powcal%decomp%myrankgeoqshad,powcal%decomp%commgeoqshad)
          ! Need to broadcast the converted positions
          if(combined) call bcast_position(gshadl%posl, &
               powcal%decomp%myrankgeoqshad,powcal%decomp%commgeoqshad)
          if(powcal%decomp%myrankpowcal == 0) call clock_stop(33)
       endif
       if(powcal%decomp%readpowcaldata) then
          !! launch from `results' geometry
          !! read  geobjl data (`results' geometry data)
          call clock_start(4,'results geobjlist time')
          call powcal_readv(powcal,file%vtkres)
          call clock_stop(4)
       end if
       if(powcal%decomp%bcastpowcaldata) then
          if(powcal%decomp%myrankpowcal == 0) call clock_start(32,'results geobjlist broadcast')
          call bcast_powcal(powcal%powres%geobjl,powcal%powres%vecx,powcal%powres%vecb,&
               powcal%decomp%myrankgeoqres,powcal%decomp%commgeoqres)
          if(combined) call bcast_position(powcal%powres%geobjl%posl,&
               powcal%decomp%myrankgeoqres,powcal%decomp%commgeoqres)
          if(powcal%decomp%myrankpowcal == 0) call clock_stop(32)
       endif
    end if
    !--------------------------------------------------------------------------
    !! read  HDS  data
    if(powcal%decomp%myrankpowcal == 0) call clock_start(6,'hds_init time')
    !     line below needed for debugging version to run
!    write(*,*) 'debugging output',file%vtkdata
    if(powcal%decomp%readpowcaldata) then
       call hdsfile_dinit(file%hdsdata,timeold)
       call hdsfile_read(gshadl,gnumerics,btree)
    endif

    if(powcal%decomp%bcastpowcaldata) then
       if(decomp%myrankpowcal == 0) call clock_start(34,"hds broadcast time")
       call bcast_hds(gshadl,gnumerics,btree,powcal%decomp%myrankhdsgen,powcal%decomp%commhdsgen)
       if(decomp%myrankpowcal == 0) call clock_stop(34)
    endif

    ! rotate and quantise shadowing geometry
    call position_tfmlis(gshadl%posl,gshadl%tfmdata)
    call position_qtfmlis(gshadl%posl,gshadl%quantfm)
    call btree_limits(btree,powcal%powres%beq%n%xbb)
    if(powcal%decomp%myrankpowcal == 0) call clock_stop(6)

    !--------------------------------------------------------------------------
    !! initialise powcal and odes
    if(powcal%decomp%myrankpowcal == 0) call clock_start(7,'initialise powcal time')
    call odes_init(powcal%odes,onumerics)
    call powcal_init(powcal,numerics,plot,gnumerics)

    if (powcal%powres%flinends .and. powcal%decomp%myrankpowcal .eq. 0) then
         call gfile_init(trim(file%flinends),'ends',powcal%powres%nflends)
    end if
    if(powcal%decomp%myrankpowcal == 0) call clock_stop(7)

    !--------------------------------------------------------------------------
    !! do the main work
    if(powcal%decomp%myrankpowcal == 0) call clock_start(8,'powcal evaluation time')
    call powcal_refine(powcal,gshadl,btree)
    if(powcal%decomp%myrankpowcal == 0) call clock_stop(8)

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !! collect all the evaluated paths
    if(powcal%decomp%myrankpowcal == 0) call clock_start(25,'MPI collect time')
    call powcalmpi_collect(powcal)
    if(powcal%decomp%myrankpowcal == 0) call clock_stop(25)

    !--------------------------------------------------------------------------

    if(powcal%decomp%myrankpowcal .eq. 0) then

       !! Flush the error buffer
       call log_flusherror

       ! field line ends - Need to check after merge
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
       !--------------------------------------------------------------------------

    endif

    !! cleanup and closedown
    call powcal_delete(powcal)
    call geobjlist_delete(gshadl)

    if(powcal%decomp%myrankpowcal == 0) then
       call clock_stop(1)
       call clock_summary
       call clock_delete
    endif

    call log_close

  end subroutine powcal_sub

end module powcal_p

#ifndef COMBINED
program powcal

  use powcal_p

  type(decomp_t) :: decomp

  !! Initialise MPI
  call powcalmpi_init(decomp)

  !--------------------------------------------------------------------------
  !! The actual program
  call powcal_sub(1,.False.,decomp)

  !--------------------------------------------------------------------------
  !! Finalize MPI
  call powcalmpi_finish(decomp)
  !--------------------------------------------------------------------------

end program powcal
#endif
