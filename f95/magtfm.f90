module magtfm_p

  use const_kind_m
  use const_numphys_h
  use date_time_m
  use mcontrol_h
  use dcontrol_h
  use apb_h
  use log_m
  use misc_m
  use clock_m
  use position_h
  use fmesh_h
  use beq_h
  use posang_h
  use posang_m
  use apb_m
  use gfile_m
  use spl3d_m
  use spl2d_m
  use smitermpi_h
  use smitermpi_m

  use position_m

  use mcontrol_m
  use moutfile_m
  use vfile_m
  use fmesh_m
  use beqart_h
  use beqart_m

  implicit none

contains
  
  subroutine magtfm_sub(arg,combined,decomp)

    ! Arguments
    integer(ki4), intent(in) :: arg !< The number of the command line argument conatining the input file
    logical, intent(in) :: combined !< Is this run as part of the combined smiter output
    type(decomp_t), intent(in) :: decomp !< is the code running as a single executable

    ! Local variables
    character(*), parameter :: m_name='magtfm' !< module name
    type(mfiles_t)     :: file      !< names of files
    type(apb_t)     :: apb      !< applied field input data structure
    type(apb_t)     :: apbsave      !< applied field saved data structure
    type(spl3d_t)     :: spl3d      !< applied field spl3d data structure
    type(mnumerics_t)  :: numerics  !< numerical control parameters
    type(date_time_t) :: timestamp !< timestamp of run
    character(len=80) :: fileroot !< reference name for all files output by run
    character(len=80), save :: iched='sp3' !< vacuum field file descriptor

    type(mplots_t)     :: plot      !< diagnostic plot selectors
    !     integer:: nplot !< unit for vtk files
    integer:: nout !< unit for mag files
    integer:: nprint !< unit for gnuplot files
    integer(ki4):: j !< loop variable
    real(kr8), dimension(3) :: zpt !< vector point
    real(kr8), dimension(3) :: zeval !< vector value
    real(kr8) :: scalereq  !< local variable
    integer(ki4), parameter :: iphchange=0 !< change phase of zeta/reflect if =3
    !--------------------------------------------------------------------------
    !! initialise timing

    call date_time_init(timestamp)
    call clock_init(30)
    call clock_start(1,'magtfm run time')
    !--------------------------------------------------------------------------
    !! print header

    print *, '----------------------------------------------------'
    print *, 'magtfm: transform applied field file'
    print *, '----------------------------------------------------'
    print '(a)', timestamp%long
    !--------------------------------------------------------------------------
    !! get file root from arg
    if(command_argument_count()<1) then
       !! no file root specified
       print *, 'Fatal error: no file root name specified.'
       print *, 'To run magtfm type at the command line:'
       print *, '   magtfm fileroot'
       stop
    else
       !!get fileroot
       call get_command_argument(1,value=fileroot)
    end if

    !! start log
    call log_init(fileroot,timestamp)
    !--------------------------------------------------------------------------
    !! read control file (nominally)

    call clock_start(2,'mcontrol_init time')
    call mcontrol_init(fileroot)
    call mcontrol_read(file,numerics,plot)
    call mcontrol_close
    call clock_stop(2)
    ! log parameters
    call log_value("Phase change parameter ",iphchange)
    call log_value("R factor in magnetic field ",spl3d_rpower)
    !--------------------------------------------------------------------------
    !! read  applied field data

    call clock_start(4,'field_data_read time')
    call apb_readmlab(apb,numerics,file%magdata)
    call apb_copy(apb,apbsave)
    call apb_checkperiod(apb)
    call clock_stop(4)

    !--------------------------------------------------------------------------
    !! output in cylindrical coordinates only

    if (trim(adjustl(file%magout))=='null') then

       if(plot%bcartx) then
          call clock_start(5,'writecart_init time')
          call moutfile_init(file%bcartx,numerics,'txt',nout)
          zpt(1)=numerics%rax
          zpt(2)=0.
          zpt(3)=0.
          call apb_eval(apb,zpt,zeval)
          scalereq=numerics%ireq/sqrt((zeval(1)**2+zeval(2)**2)*(zpt(1)**2+zpt(2)**2))
          call apb_scale(apb,numerics,scalereq)
          call apb_qartwrite(apb,numerics,nout)
          call moutfile_close()
          call clock_stop(5)
       end if
       !to tidy up and close down
       go to 10

    end if

    !dbg21   if(plot%bcartx) then    !dbg21
    !dbg21      call clock_start(5,'writecart_init time')    !dbg21
    !dbg21      call moutfile_init(file%bcartx,numerics,'txt',nout)    !dbg21
    !dbg21      call apb_qartwrite(apb,numerics,nout)    !dbg21
    !dbg21      call moutfile_close()    !dbg21
    !dbg21      call clock_stop(5)    !dbg21
    !dbg21      stop    !dbg21
    !dbg21   end if    !dbg21

    !--------------------------------------------------------------------------
    !! do the conversion
    call clock_start(6,'field transform time')
    !     line below needed for debugging version to run
    !  write(*,*) 'this output helps gfortran sometimes',(apb%pos3(j),j=1,2)
    !     write(*,*) 'second',((numerics%panbod(i,j),i=1,2),j=1,9)
    call apb_convert(apb,numerics,0,spl3d_rpower,iphchange)
    call spl3d_init(spl3d,apb%field,apb%ncpt,numerics%n0,&
         &apb%n1,apb%n2,apb%n3,&
         &apb%pos1,apb%pos2,apb%pos3,numerics%nord,&
         &numerics%kmin,numerics%kmax)
    call apb_copy(apbsave,apb)
    call clock_stop(6)
    !--------------------------------------------------------------------------
    !! scale

    call clock_start(7,'field scaling time')
    zpt(1)=numerics%rax
    zpt(2)=0.
    zpt(3)=0.
    call spl3d_eval(spl3d,zpt,zeval)
    scalereq=numerics%ireq/(zeval(3)*numerics%rax**(2-spl3d_rpower))
    call spl3d_scale(spl3d,scalereq)
    call clock_stop(7)
    !--------------------------------------------------------------------------
    !! set up the mode mask

    call clock_start(8,'mode masking time')
    call spl3d_mask(spl3d,numerics%cutoff,numerics%parity,numerics%ireq)
    call clock_stop(8)
    !--------------------------------------------------------------------------
    ! gnu plots
    if(plot%modes) then
       call clock_start(9,'gfile_modes time')
       call gfile_init(trim(file%modes),'maximum mode numbers',nprint)
       if (minval(numerics%kmax)<0) then
          call spl3d_writeg(spl3d,'countout',nprint)
       else
          call spl3d_writeg(spl3d,'modecount',nprint)
       end if
       call gfile_close
       call clock_stop(9)
    end if
    !--------------------------------------------------------------------------
    !! testing transform

    if(plot%testtfm) then
       call clock_start(30,'spl3d_testtfm time')
       call spl3d_testtfm(spl3d,apbsave%field,scalereq,iphchange)
       call clock_stop(30)
    end if
    !--------------------------------------------------------------------------
    !!plots
    !! diagnostic gnu plots
    if(plot%gnuv) then
       call clock_start(40,'gfile_gnuv time')
       call gfile_init(trim(file%gnuv),'field sample in R-Z space',nprint)
       call spl3d_writeg(spl3d,'sampl',nprint)
       call gfile_close
       call clock_stop(40)
    end if

    if(plot%gnu) then
       call clock_start(42,'gfile_gnu time')
       call gfile_init(trim(file%gnu),'field sample in k-space (Bzeta, BZ,BR)',nprint)
       call spl3d_writeg(spl3d,'strip',nprint)
       call gfile_close
       call clock_stop(42)
    end if
    !--------------------------------------------------------------------------
    !! output files
    call clock_start(32,'moutfile_init time')
    call moutfile_init(file%magout,numerics,iched,nout)
    !     write(*,*) 'iched=',iched
    call spl3d_write(spl3d,nout)
    call moutfile_close()
    call clock_stop(32)
    !--------------------------------------------------------------------------
    !! cleanup and closedown
    call spl3d_delete(spl3d)

10  continue
    call clock_stop(1)
    call clock_summary

    call log_close
    call clock_delete
    !--------------------------------------------------------------------------

  end subroutine magtfm_sub
  
end module magtfm_p

#ifndef COMBINED
!> @addtogroup groupname9
!> @{
program magtfm
!> @} 
  use magtfm_p
  
  type(decomp_t) :: decomp

  if(serialmpi_init()) call magtfm_sub(1,.False.,decomp)

  call serialmpi_finish()
  
end program magtfm
#endif
