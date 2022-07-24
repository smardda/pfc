!> @addtogroup groupname8
!> @{
module smitermpi_m
!> @}
  use smitermpi_h
  use log_m
  use powcal_h
  use powres_h
  use pcontrol_h
  use control_h
  use edgprof_h
  use odes_h
  use termplane_h
  use beq_h
  use spl2d_m
  use geobjlist_m
  use btree_m
  use ls_m
#ifdef WITH_MPI
  use mpi
#endif

  implicit none

  private

  public :: &
       smitermpi_init, &   !> initialise mpi
       smitermpi_finish, &    !> clean up and finalise mpi
       powcalmpi_init, &   !> initialise mpi
       powcalmpi_finish, &    !> clean up and finalise mpi
       serialmpi_init, &  !> initialises serial codes to run on 1 process
       serialmpi_finish, &  !> finished serial codes to running on 1 process
       powcalmpi_decomp, &  !> decompose data and tasks
       powcalmpi_collect, &   !> collect all the trace data to rank 0 for output
       bcast_pcontrol, &    !> broadcasts pcontrol_t values for powcal
       bcast_beq, &         !> broadcasts beq_t values for powcal
       bcast_powcal, &      !> broadcasts pcontrol_t values for powcal
       bcast_geobjlist, &   !> broadcasts geobjlist_t values for powcal and hdsgen
       bcast_hds, &        !> broadcasts hds values for powcal
       bcast_position, &   !> broadcasts position values for powcal and hdsgen
       bcast_skyl, &       !> broadcasts skylight_t values for powcal and hdsgen
       smitermpi_barrier !>

  character(*), parameter :: m_name='smitermpi_m' !< module name
  integer :: status !<

contains

  !!-----------------------------------------------------------------------------------
  !! Initialisationa and decomposition routines
  !!-----------------------------------------------------------------------------------

  !> calls mpi_init and defines decomp structure
  subroutine smitermpi_init(self)

    integer :: ierr !<
    logical :: init_flag !<

    type(decomp_t), intent(inout) :: self !< decomp data structure

#ifdef WITH_MPI
    call MPI_Initialized(init_flag,ierr)
    if(.not. init_flag) call MPI_Init(ierr)

    ! Duplicate MPI_COMM_WORLD to future proof the code
    call MPI_Comm_dup(MPI_COMM_WORLD,self%commsmiter,ierr)

    call MPI_Comm_rank(self%commsmiter,self%myrank,ierr)
    call MPI_Comm_size(self%commsmiter,self%nproc,ierr)
    if(self%myrank .eq. 0) print*, "Initialised MPI with", self%nproc ,"processes"
#else
    self%myrank=0
    self%nproc=1
#endif
    call smitermpi_control(self)

  end subroutine smitermpi_init

  !> assigns processes to each task and creates groups and communicators for the single exe smiter code
  !> this will become more complex once multiple time slices are used
  subroutine smitermpi_control(self)

    type(decomp_t), intent(inout) :: self !< decomp data structure

    ! Booleans to control running of individual steps
    if(self%nproc > 6) then
       self%procgeoqshad = 0
       self%procgeoqres = 1
       self%prochdsgen = 2
       if(self%myrank > 2) self%powcal = .True.
       self%bcastpowcal=.True.
       self%bcasthdsgen=.True.
       call smitermpi_decomp(self)
    else
       self%procgeoqshad = 0
       self%procgeoqres = 0
       self%prochdsgen = 0
       self%powcal = .True.
       self%commpowcal = self%commsmiter
       self%commgeoqshad = self%commsmiter
       self%commgeoqres = self%commsmiter
       self%commhdsgen = self%commsmiter
       self%nprocgeoqshad = 1
       self%nprocgeoqres = 1
       self%nprochdsgen = 1
       self%nprocpowcal = self%nproc
       self%myrankgeoqshad = self%myrank
       self%myrankgeoqres = self%myrank
       self%myrankhdsgen = self%myrank
       self%myrankpowcal = self%myrank
       if(self%nproc > 1) self%procgeoqres=1
       if(self%nproc > 2) then
          self%prochdsgen=2
          self%bcasthdsgen=.True.
          call smitermpi_decomp(self,.True.)
       endif
    endif

    ! now set the read and broadcast behaviour for powcal
    call powcalmpi_setio(self)

    ! Assign the global myrank_log and nproc_log
    if(self%myrank == self%procgeoqshad) myrank_log = self%myrankgeoqshad
    if(self%myrank == self%procgeoqres) myrank_log = self%myrankgeoqres
    if(self%myrank == self%prochdsgen) myrank_log = self%myrankhdsgen
    if(self%powcal) myrank_log = self%myrankpowcal

    ! This is for log ordering in the serial version
    nproc_log = self%nproc

  end subroutine smitermpi_control

  subroutine smitermpi_decomp(self,hdsonly)

    !! Arguments
    type(decomp_t), intent(inout) :: self !< decomp data structure
    logical, optional :: hdsonly !<

    integer(ki4) :: ierr !<
    integer(ki4) :: groupsmiter,groupgeoqshad,groupgeoqres,grouphdsgen,grouppowcal !<
    integer(ki4) :: grouppg,groupps,groupph,groupsh !<
    integer(ki4) :: excl(3),incl(1) !<

#ifdef WITH_MPI
    excl(1) = 0
    excl(2) = 1
    excl(3) = 2

    ! Create groups and communicators for geoq, hdsgen, powcal bcasts
    call MPI_Comm_group(self%commsmiter,groupsmiter,ierr)
    ! Create individual groups
    ! geoqshad
    incl = 0
    call MPI_Group_incl(groupsmiter,1,incl,groupgeoqshad,ierr)
    ! geoqres
    incl = 1
    call MPI_Group_incl(groupsmiter,1,incl,groupgeoqres,ierr)
    ! hdsgen
    incl = 2
    call MPI_Group_incl(groupsmiter,1,incl,grouphdsgen,ierr)

    ! Now do the geoqshad hdsgen union and return if hdsonly is true
    ! geoqshad and hdsgen
    call MPI_Group_union(groupgeoqshad,grouphdsgen,groupsh,ierr)
    call MPI_Comm_create(self%commsmiter,groupsh,self%commshadhds,ierr)
    call MPI_Group_rank(groupsh,self%myrankshadhds,ierr)
    call MPI_Group_size(groupsh,self%nprocshadhds,ierr)

    if(present(hdsonly)) return

    ! Now the rest of the unions
    ! powcal
    call MPI_Group_excl(groupsmiter,3,excl,grouppowcal,ierr)
    call MPI_Comm_create(self%commsmiter,grouppowcal,self%commpowcal,ierr)
    call MPI_Group_rank(grouppowcal,self%myrankpowcal,ierr)
    call MPI_Group_size(grouppowcal,self%nprocpowcal,ierr)

    ! geoqshad and powcal
    call MPI_Group_union(groupgeoqshad,grouppowcal,groupps,ierr)
    call MPI_Comm_create(self%commsmiter,groupps,self%commgeoqshad,ierr)
    call MPI_Group_rank(groupps,self%myrankgeoqshad,ierr)
    call MPI_Group_size(groupps,self%nprocgeoqshad,ierr)
    ! geoqres and powcal
    call MPI_Group_union(groupgeoqres,grouppowcal,grouppg,ierr)
    call MPI_Comm_create(self%commsmiter,grouppg,self%commgeoqres,ierr)
    call MPI_Group_rank(grouppg,self%myrankgeoqres,ierr)
    call MPI_Group_size(grouppg,self%nprocgeoqres,ierr)
    ! hdsgen and and powcal
    call MPI_Group_union(grouphdsgen,grouppowcal,groupph,ierr)
    call MPI_Comm_create(self%commsmiter,groupph,self%commhdsgen,ierr)
    call MPI_Group_rank(groupph,self%myrankhdsgen,ierr)
    call MPI_Group_size(groupph,self%nprochdsgen,ierr)
#endif

  end subroutine smitermpi_decomp

  !> calls mpi_init and defines decomp structure
  subroutine powcalmpi_init(self)

    integer :: ierr !<
    logical :: init_flag !<

    type(decomp_t), intent(inout) :: self !< decomp data structure

#ifdef WITH_MPI
    call MPI_Initialized(init_flag,ierr)
    if(.not. init_flag) call MPI_Init(ierr)

    ! Duplicate MPI_COMM_WORLD to future proof the code
    call MPI_Comm_dup(MPI_COMM_WORLD,self%commpowcal,ierr)
    ! Set the combined communicators to be equivalent to the commpowcal
    self%commgeoqshad = self%commpowcal
    self%commgeoqres = self%commpowcal
    self%commhdsgen = self%commpowcal

    call MPI_Comm_rank(self%commpowcal,self%myrank,ierr)
    call MPI_Comm_size(self%commpowcal,self%nproc,ierr)
    if(self%myrank .eq. 0) print*, "Initialised MPI with", self%nproc ,"processes"
#else
    self%myrank = 0
    self%nproc = 1
#endif
    self%nprocgeoqshad = self%nproc
    self%nprocgeoqres = self%nproc
    self%nprochdsgen = self%nproc
    self%nprocpowcal = self%nproc

    self%myrankgeoqshad = self%myrank
    self%myrankgeoqres = self%myrank
    self%myrankhdsgen = self%myrank
    self%myrankpowcal = self%myrank

    myrank_log = self%myrank

    ! now set the read and broadcast behaviour for powcal
    self%bcastpowcal = .False. ! This init is for standalone smiter
    call powcalmpi_setio(self)

  end subroutine powcalmpi_init

  !> calls mpi_finalize and cleans up decomp structure
  subroutine smitermpi_finish(self)

    integer :: ierr !<
    logical :: init_flag !<
    type(decomp_t), intent(inout) :: self !< decomp data structure

    if(allocated(self%localsize)) deallocate(self%localsize)
    if(allocated(self%localdisp)) deallocate(self%localdisp)

#ifdef WITH_MPI
    if(self%myrank .eq. 0) print*, "Finalizing MPI"

    call MPI_Initialized(init_flag,ierr)
    if(init_flag) call MPI_Finalize(ierr)
#endif

  end subroutine smitermpi_finish

  !> calls mpi_finalize and cleans up decomp structure
  subroutine powcalmpi_finish(self)

    integer :: ierr !<
    logical :: init_flag !<
    type(decomp_t), intent(inout) :: self !< decomp data structure

    if(allocated(self%localsize)) deallocate(self%localsize)
    if(allocated(self%localdisp)) deallocate(self%localdisp)
#ifdef WITH_MPI
    if(self%myrank .eq. 0) print*, "Finalizing MPI"

    call MPI_Initialized(init_flag,ierr)
    if(init_flag) call MPI_Finalize(ierr)
#endif

  end subroutine powcalmpi_finish

  !> initialises mpi for "serial" codes that run on single processes during mpi execution
  function serialmpi_init()

    integer :: ierr, myrank, commserial !<
    logical :: init_flag, serialmpi_init !<

    serialmpi_init = .False.
    myrank = 0

#ifdef WITH_MPI
    call MPI_Initialized(init_flag,ierr)
    if(.not. init_flag) call MPI_Init(ierr)
    ! Duplicate MPI_COMM_WORLD to future proof the code
    call MPI_Comm_dup(MPI_COMM_WORLD,commserial,ierr)

    call MPI_Comm_rank(commserial,myrank,ierr)
#endif

    if(myrank == 0) serialmpi_init = .True.

  end function serialmpi_init

  !> calls mpi_finalize for "serial" codes that run on single processes during mpi execution
  subroutine serialmpi_finish()

    integer :: ierr !<
    logical :: init_flag !<

#ifdef WITH_MPI
    call MPI_Initialized(init_flag,ierr)
    if(init_flag) call MPI_Finalize(ierr)
#endif

  end subroutine serialmpi_finish

  !> sets reading and broadcast flags for powcal
  subroutine powcalmpi_setio(self)

    integer :: ierr !<

    type(decomp_t), intent(inout) :: self !< decomp data structure

    self%bcastflag = .True.

    if(self%bcastflag) then  ! broadcasts are allowed in powcal
       if(self%myrank == 0) print*,"Powcal data will be broadcasted."
       if(self%nprocpowcal > 1) then  ! If there is more than one powcal process bcasts are required
          self%bcastpowcalcontrol = .True.
          self%bcastpowcaldata = .True.
       else                             ! Otherwise they aren't
          self%bcastpowcalcontrol = .False.
          self%bcastpowcaldata = .False.
       endif
       if(self%myrankpowcal == 0) then  ! powcal process 0 always reads the control data
          self%readpowcalcontrol = .True.
       else
          self%readpowcalcontrol = .False.
       endif
       if(self%myrankpowcal == 0 .and. .not. self%bcastpowcal) then  ! When non combined broadcasts are used only powcal process 0 reads data
          self%readpowcaldata = .True.
       else
          self%readpowcaldata = .False.
       endif
    else  ! broadcasts are turned off and every process must read all data
       if(self%myrank == 0) print*,"Powcal data will NOT be broadcasted - all processes will read from disk"
       self%bcastpowcalcontrol = .False.
       self%bcastpowcaldata = .False.
       self%readpowcalcontrol = .True.
       self%readpowcaldata = .True.
    endif

!    print*, "Proc ",self%myrank, self%bcastpowcalcontrol, self%bcastpowcaldata, self%readpowcalcontrol, self%readpowcaldata


  end subroutine powcalmpi_setio


  !> assigns a set of starting points to each process
  subroutine powcalmpi_decomp(self,npowe,ng)

    integer :: ierr !<
    integer, intent(in) :: npowe,ng !<
    type(decomp_t), intent(inout) :: self !< decomp data structure

!    if(self%myrankpowcal .eq. 0) print*, "Distributing", npowe ,"trajectories"
!    if(self%myrankpowcal .eq. 0) &
!         print*, "Powcal distributing trajectories over",self%nprocpowcal,"processes"

    ! Allocate the array containing the number of trajectories for each process
    allocate(self%localsize(0:self%nprocpowcal-1))
    allocate(self%localdisp(0:self%nprocpowcal-1))
    self%localsize = 0
    self%localdisp = 0

    call powcalmpi_localsize(self,npowe)
    self%firstpowe = self%localdisp(self%myrankpowcal) + 1
    self%mynpowe = self%localsize(self%myrankpowcal)
    self%lastpowe = self%firstpowe + self%mynpowe - 1

    ! This is temporary and should be dealt with properly
    self%firstng = self%firstpowe
    self%lastng = self%lastpowe


!    print*, "Process ",self%myrank,"calculating",self%mynpowe,"trajectories:",self%firstpowe,"to",self%lastpowe

  end subroutine powcalmpi_decomp

  !> Calculates the number of trajectories for each process and
  !> stores it in an array for use by the collect step
  subroutine powcalmpi_localsize(self,npowe)

    integer, intent(in) :: npowe !<
    type(decomp_t), intent(inout) :: self !< decomp data structure
    integer :: i, proc !<

    if(self%nprocpowcal .eq. 1) then
       self%localsize(0) = npowe
       return
    end if

    proc = 0
    i = 1
    do while(i .le. npowe)
       if(proc == self%nprocpowcal) proc = 0
       self%localsize(proc) = self%localsize(proc) + 1
       self%localdisp(proc+1:self%nprocpowcal-1) = self%localdisp(proc+1:self%nprocpowcal-1) + 1
       proc = proc + 1
       i = i + 1
    end do

  end subroutine powcalmpi_localsize


  subroutine smitermpi_barrier(comm)

    integer, intent(in) :: comm !<

    integer :: ierr !<
#ifdef WITH_MPI
    call MPI_Barrier(comm,ierr)
#endif

  end subroutine smitermpi_barrier

  !!----------------------------------------------------------------------------------
  !! Communication routines for broadcasting, scattering and collecting data
  !!----------------------------------------------------------------------------------


  !> Collects data from all processes for writing to disk
  subroutine powcalmpi_collect(self)

    integer :: ierr !<
    type(powcal_t), intent(inout) :: self   !< object data structure

#ifdef WITH_MPI
    ! MPI_IN_PLACE is used here, it makes the code a little ugly
    ! but saves a little memory and time. It only works if root and
    ! the other processes have different calls.

    if(self%decomp%myrankpowcal .eq. 0) then
       call MPI_Gatherv(MPI_IN_PLACE,self%decomp%mynpowe,MPI_REAL, &
            self%powres%pows,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
            0,self%decomp%commpowcal,ierr)

       call MPI_Gatherv(MPI_IN_PLACE,self%decomp%mynpowe,MPI_REAL, &
            self%powres%powa,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
            0,self%decomp%commpowcal,ierr)

       if (allocated(self%powres%psista)) then
          call MPI_Gatherv(MPI_IN_PLACE,self%decomp%mynpowe,MPI_REAL, &
               self%powres%psista,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
               0,self%decomp%commpowcal,ierr)
       endif

       if (allocated(self%powres%angle)) then
          call MPI_Gatherv(MPI_IN_PLACE,self%decomp%mynpowe,MPI_REAL, &
               self%powres%angle,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
               0,self%decomp%commpowcal,ierr)
       endif

       if (allocated(self%powres%lenpath)) then
          call MPI_Gatherv(MPI_IN_PLACE,self%decomp%mynpowe,MPI_REAL, &
               self%powres%lenpath,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
               0,self%decomp%commpowcal,ierr)
       endif

       if (allocated(self%powres%zbdotn)) then
          call MPI_Gatherv(MPI_IN_PLACE,self%decomp%mynpowe,MPI_REAL, &
               self%powres%zbdotn,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
               0,self%decomp%commpowcal,ierr)
       endif

       if (allocated(self%powres%bfield)) then
          call MPI_Gatherv(MPI_IN_PLACE,3*self%decomp%mynpowe,MPI_DOUBLE, &
               self%powres%bfield,3*self%decomp%localsize,3*self%decomp%localdisp,MPI_DOUBLE, &
               0,self%decomp%commpowcal,ierr)
       endif

       if (allocated(self%powres%norm)) then
          call MPI_Gatherv(MPI_IN_PLACE,3*self%decomp%mynpowe,MPI_DOUBLE, &
               self%powres%norm,3*self%decomp%localsize,3*self%decomp%localdisp,MPI_DOUBLE, &
               0,self%decomp%commpowcal,ierr)
       endif

       call MPI_Gatherv(MPI_IN_PLACE,self%decomp%mynpowe,MPI_REAL, &
            self%powres%pow,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
            0,self%decomp%commpowcal,ierr)

       call MPI_Gatherv(MPI_IN_PLACE,self%decomp%mynpowe,MPI_INTEGER, &
            self%powres%pmask,self%decomp%localsize,self%decomp%localdisp,MPI_INTEGER, &
            0,self%decomp%commpowcal,ierr)
    else
       call MPI_Gatherv(self%powres%pows,self%decomp%mynpowe, &
            MPI_REAL,self%powres%pows,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
            0,self%decomp%commpowcal,ierr)

       call MPI_Gatherv(self%powres%powa,self%decomp%mynpowe,MPI_REAL, &
            self%powres%powa,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
            0,self%decomp%commpowcal,ierr)

       if (allocated(self%powres%psista)) then
          call MPI_Gatherv(self%powres%psista,self%decomp%mynpowe,MPI_REAL, &
               self%powres%psista,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
               0,self%decomp%commpowcal,ierr)
       end if

       if (allocated(self%powres%angle)) then
          call MPI_Gatherv(self%powres%angle,self%decomp%mynpowe,MPI_REAL, &
               self%powres%angle,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
               0,self%decomp%commpowcal,ierr)
       end if

       if (allocated(self%powres%lenpath)) then
          call MPI_Gatherv(self%powres%lenpath,self%decomp%mynpowe,MPI_REAL, &
               self%powres%lenpath,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
               0,self%decomp%commpowcal,ierr)
       end if

       if (allocated(self%powres%zbdotn)) then
          call MPI_Gatherv(self%powres%zbdotn,self%decomp%mynpowe,MPI_REAL, &
               self%powres%zbdotn,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
               0,self%decomp%commpowcal,ierr)
       end if

       if (allocated(self%powres%bfield)) then
          call MPI_Gatherv(self%powres%bfield,3*self%decomp%mynpowe,MPI_DOUBLE, &
               self%powres%bfield,3*self%decomp%localsize,3*self%decomp%localdisp,MPI_DOUBLE, &
               0,self%decomp%commpowcal,ierr)
       end if

       if (allocated(self%powres%norm)) then
          call MPI_Gatherv(self%powres%norm,3*self%decomp%mynpowe,MPI_DOUBLE, &
               self%powres%norm,3*self%decomp%localsize,3*self%decomp%localdisp,MPI_DOUBLE, &
               0,self%decomp%commpowcal,ierr)
       end if

       call MPI_Gatherv(self%powres%pow,self%decomp%mynpowe,MPI_REAL, &
            self%powres%pow,self%decomp%localsize,self%decomp%localdisp,MPI_REAL, &
            0,self%decomp%commpowcal,ierr)

       call MPI_Gatherv(self%powres%pmask,self%decomp%mynpowe,MPI_INTEGER, &
            self%powres%pmask,self%decomp%localsize,self%decomp%localdisp,MPI_INTEGER, &
            0,self%decomp%commpowcal,ierr)

    endif
#endif
  end subroutine powcalmpi_collect


  !!----------------------------------------------------------------------------------
  !! Control file broadcast routines
  !!----------------------------------------------------------------------------------

  !> broadcasts info from the control file
  subroutine bcast_pcontrol(file,numerics,onumerics,edgprof,plot,rank,comm)

    !! arguments
    type(pfiles_t), intent(inout) :: file !< file names ! Temporary

    type(pnumerics_t), intent(inout) :: numerics !< input numerical parameters
    type(onumerics_t), intent(inout) :: onumerics !< input ODE numerical parameters
    type(edgprof_t), intent(inout) :: edgprof !< input ODE numerical parameters
    type(pplots_t), intent(out) :: plot !< vtk plot selectors
    integer(ki4) , intent(in) :: rank,comm !<

    call bcast_pnumerics(numerics,rank,comm)

    call bcast_edgprof(edgprof,rank,comm)

    call bcast_odes(onumerics,rank,comm)

    call bcast_file(file,rank,comm) ! Temporary

    call bcast_pplot(plot,rank,comm)

  end subroutine bcast_pcontrol

  !> Calls MPI Bcast to distribute the pnumerics structure to all processes
  subroutine bcast_pnumerics(self,rank,comm)

    !! Arguments
    type(pnumerics_t), intent(inout) :: self !< input numerical parameters
    integer(ki4) , intent(in) :: rank,comm !<

    integer :: ierr !<

#ifdef WITH_MPI
    call MPI_Bcast(self%f,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%lmid,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%ploss,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%sigma,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%nlevel,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%shadow,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%qpara0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%caltype,80,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(self%mtrack,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%ntrack,1,MPI_INTEGER,0,comm,ierr)

    ! Added during conversion to LIB/PFC structure
    call MPI_Bcast(self%usepsista,1,MPI_LOGICAL,0,comm,ierr)

    if(self%ntrack > 0) then
       if(rank .gt. 0) allocate(self%trackno(self%ntrack))
       call MPI_Bcast(self%trackno,self%ntrack,MPI_INTEGER,0,comm,ierr)
    endif

    call MPI_Bcast(self%nanalau,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nlaupow,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nshapow,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%ltermplane,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(self%ledgprof,1,MPI_LOGICAL,0,comm,ierr)

    ! Skylight
    call MPI_Bcast(self%lskyl,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(self%lflagtol,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(self%lpfpower,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(self%tolin,3,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%lskylpfr,1,MPI_LOGICAL,0,comm,ierr)

    if(self%ltermplane) then
       call bcast_termplanes(self%termp,rank,comm)
    else
       self%termp%ntermplane=0
    end if
#endif
  end subroutine bcast_pnumerics

  !> Calls MPI Bcast to distribute the termplane structure to all processes
  subroutine bcast_termplanes(self,rank,comm)

    !! Arguments
    type(termplane_t), intent(inout) :: self   !< termination planes numeric controls
    integer(ki4) , intent(in) :: rank,comm

    integer :: ierr
#ifdef WITH_MPI
    call MPI_Bcast(self%ntermplane,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%ntermactive,1,MPI_INTEGER,0,comm,ierr)
    if(self%ntermactive > 0) then
       if(rank .gt. 0) then
          allocate(self%termplanedir(self%ntermactive+1,5),stat=status)
          allocate(self%termplane(self%ntermactive+1,2),stat=status)
          allocate(self%termstore(self%ntermactive+1,3),stat=status)
       endif
       call MPI_Bcast(self%termplanedir,5*(self%ntermactive+1),MPI_INTEGER,0,comm,ierr)
       call MPI_Bcast(self%termplane,2*(self%ntermactive+1),MPI_REAL,0,comm,ierr)
       call MPI_Bcast(self%termstore,3*(self%ntermactive+1),MPI_REAL,0,comm,ierr)
    endif
#endif
  end subroutine bcast_termplanes

  !> Calls MPI Bcast to distribute the edgprof structure to all processes
  subroutine bcast_edgprof(self,rank,comm)

    !! Arguments
    type(edgprof_t), intent(inout) :: self !< input ODE numerical parameters
    integer(ki4) , intent(in) :: rank,comm !<

    integer :: ierr
#ifdef WITH_MPI
    call MPI_Bcast(self%formula,80,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(self%lmid,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%f,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%ploss,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%sigma,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%qpara0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%lmidnr,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%rqpara0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%postype,80,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(self%npos,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nrpams,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nipams,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%fint,1,MPI_DOUBLE_PRECISION,0,comm,ierr)

    formula_allocate: select case (self%formula)
    case('samples')
       allocate(self%pos(self%npos), self%prof(self%npos), stat=status)
       call MPI_Bcast(self%pos,self%npos,MPI_DOUBLE_PRECISION,0,comm,ierr)
       call MPI_Bcast(self%prof,self%npos,MPI_DOUBLE_PRECISION,0,comm,ierr)

    case('userdefined')
       if (self%nrpams>0) allocate(self%rpar(self%nrpams), stat=status)
       if (self%nipams>0) allocate(self%rpar(self%nipams), stat=status)
       call MPI_Bcast(self%rpar,self%nrpams,MPI_DOUBLE_PRECISION,0,comm,ierr)
       call MPI_Bcast(self%npar,self%nipams,MPI_INTEGER,0,comm,ierr)
       ! fint is nowhere initialised. Probably comes from user defined f90
       ! function
       call MPI_Bcast(self%fint,1,MPI_DOUBLE_PRECISION,0,comm,ierr)

    case default
    end select formula_allocate

    ! call MPI_Bcast(self%rblfac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    ! call MPI_Bcast(self%fpfac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    ! call MPI_Bcast(self%rblfacnr,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    ! call MPI_Bcast(self%fpfacnr,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    ! call MPI_Bcast(self%slfac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
#endif
  end subroutine bcast_edgprof

  !> Calls MPI Bcast to distribute the onumerics structure to all processes
  subroutine bcast_odes(self,rank,comm)

    !! Arguments
    type(onumerics_t), intent(inout) :: self !< input ODE numerical parameters
    integer(ki4) , intent(in) :: rank,comm !<

    integer :: ierr !<
#ifdef WITH_MPI
    call MPI_Bcast(self%npdim,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nstartcon,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%ntermcon,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%termcon,3,MPI_DOUBLE_PRECISION,0,comm,ierr)

    call MPI_Bcast(self%epsr,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%epsa,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%dt0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%stepmax,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%tmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
#endif
  end subroutine bcast_odes

  !> Calls MPI Bcast to distribute filenames to all processes - this will be removed during optmisation
  subroutine bcast_file(self,rank,comm)

    !! Arguments
    type(pfiles_t), intent(inout) :: self !< file names
    integer(ki4) , intent(in) :: rank,comm !<

    integer :: ierr
#ifdef WITH_MPI
    call MPI_Bcast(self%vtkdata,80,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(self%vtkres,80,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(self%hdsdata,80,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(self%geoq,80,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(self%laudata,80,MPI_CHARACTER,0,comm,ierr)
#endif
  end subroutine bcast_file


  subroutine bcast_pplot(self,rank,comm)

    !! Arguments
    type(pplots_t), intent(inout) :: self !< vtk plot selectors
    integer(ki4) , intent(in) :: rank,comm !<

    integer :: ierr !<
#ifdef WITH_MPI
    ! Bcast the plot structure which contains 15 logicals that are commented 
    call MPI_Bcast(self%cartv,15,MPI_LOGICAL,0,comm,ierr) 
!    call MPI_Bcast(self%cartv,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%powstatx,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%allcartv,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%powx,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%wall,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%ptzv,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%flincart,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%flinx,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%flinends,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%flinptz,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%flinm,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%allptzv,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%geofld,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%geofldq,1,MPI_LOGICAL,0,comm,ierr)
!    call MPI_Bcast(self%gnu,1,MPI_LOGICAL,0,comm,ierr)
#endif
  end subroutine bcast_pplot

  !!----------------------------------------------------------------------------------
  !! Equilibrium broadcast routines - beq
  !!----------------------------------------------------------------------------------

  subroutine bcast_beq(self,rank,comm)

    !! Arguments
    type(beq_t), intent(inout) :: self   !< object data structure
    integer(ki4) , intent(in) :: rank,comm !<
    integer :: ierr !<

#ifdef WITH_MPI
    ! Quite a few variables are missing from this routine however
    ! it currently works

    ! bnumerics broadcast - n in the beq structure
    call bcast_bnumerics(self%n,rank,comm)

    call MPI_Bcast(self%mr,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%mz,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%rmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%rmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%rmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%dr,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%zmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%zmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%dz,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%psiaxis,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%psibdry,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%rbdry(1),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%rbdry(2),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%psiqbdry,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%bpbdry,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%rxpt,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%zxpt,1,MPI_DOUBLE_PRECISION,0,comm,ierr)

    ! btotbdry - need to check if
    call MPI_Bcast(self%btotbdry,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_Bcast(self%nzets,1,MPI_INTEGER,0,comm,ierr)
    if(rank .gt. 0) allocate(self%f(self%mr+1), stat=status)
    call MPI_Bcast(self%f,self%mr+1,MPI_DOUBLE_PRECISION,0,comm,ierr)

    ! Splines
    fld_iter: select case (self%n%fldspec)
    case(1)         ! ITER case
       call bcast_spl2d(self%r,rank,comm)
       if(rank .gt. 0 ) call spl2d_initpart(self%r)
       call bcast_spl2d(self%z,rank,comm)
       if(rank .gt. 0 ) call spl2d_initpart(self%z)
       call bcast_spl2d(self%rjac,rank,comm)
       if(rank .gt. 0 ) call spl2d_initpart(self%rjac)
    case default
       call bcast_spl2d(self%psi,rank,comm)
       if(rank .gt. 0 ) call spl2d_initpart(self%psi)
       call bcast_spl2d(self%dpsidr,rank,comm)
       if(rank .gt. 0 ) call spl2d_initpart(self%dpsidr)
       call bcast_spl2d(self%dpsidz,rank,comm)
       if(rank .gt. 0 ) call spl2d_initpart(self%dpsidz)
    end select fld_iter

    fld_specn: select case (self%n%fldspec)
    case(3)
       call bcast_spl2d( self%rispldr,rank,comm)
       if(rank .gt. 0 ) call spl2d_initpart( self%rispldr )
       call bcast_spl2d( self%rispldz,rank,comm)
       if(rank .gt. 0 ) call spl2d_initpart( self%rispldz )
    end select fld_specn

    ! ripple field data
    call MPI_Bcast(self%ivac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)

    call MPI_Bcast(self%nctrack,1,MPI_INTEGER,0,comm,ierr)
    if(self%nctrack > 0) then
       if(rank .gt. 0) allocate(self%ctrackrz(self%nctrack,2), stat=status)
       call MPI_Bcast(self%ctrackrz,self%nctrack*2,MPI_DOUBLE_PRECISION,0,comm,ierr)
    endif

    call MPI_Bcast(self%rcenq,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%zcenq,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
#endif

  end subroutine bcast_beq

  subroutine bcast_bnumerics(self,rank,comm)

    !! Arguments
    type(bnumerics_t), intent(inout) :: self   !< object data structure
    integer(ki4) , intent(in) :: rank,comm !<

    integer :: ierr !<
#ifdef WITH_MPI
    ! Vac field type, fieldspec and centres
    call MPI_Bcast(self%vactype,3,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(self%fldspec,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%rcen,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%zcen,1,MPI_DOUBLE_PRECISION,0,comm,ierr)

    ! btotbbdry - need to check if
    call MPI_Bcast(self%zetamin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%zetamax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%ximin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%ximax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)

    ! ripple field data
    call MPI_Bcast(self%mrip,1,MPI_INTEGER,0,comm,ierr)
    if (self%mrip/=0) then
       call MPI_Bcast(self%arip,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    endif
    ! Won't need vacfile in futire as reading will only be done by proc 0
    call MPI_Bcast(self%vacfile,80,MPI_CHARACTER,0,comm,ierr)

    ! Mesh file would come here but I'm not sure this is necessary so will leave it for now

    ! Skylights
    call MPI_Bcast(self%objadd,11,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%skyl,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(self%skylpsi,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(self%skylcen,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(self%skyldbg,1,MPI_INTEGER,0,comm,ierr)
#endif

  end subroutine bcast_bnumerics

  subroutine bcast_spl2d(self,rank,comm)

    !! Arguments
    type(spl2d_t), intent(inout) :: self !< Spline struct
    integer(ki4) , intent(in) :: rank,comm !<

    !! Variables
    character(*), parameter :: s_name='bcast_spl2d' !< subroutine name
    integer :: splsize,ierr !<

#ifdef WITH_MPI
    call MPI_Bcast(self%lunif,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%h1,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%h2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%rh1,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%rh2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%n1,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%n2,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%n1p,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%n2p,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%pad1,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%pad2,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nord,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%noff,1,MPI_INTEGER,0,comm,ierr)
    !! org1 and org2 are set in spl2d_initpart

    ! Allocate and bcast 2D spline coeffs
    splsize = self%n1p*self%n2p
    if(rank .gt. 0) allocate(self%sampl(self%n1p,self%n2p), stat=status)
    call log_alloc_check(m_name,s_name,40,status)
    call MPI_Bcast(self%sampl,splsize,MPI_DOUBLE_PRECISION,0,comm,ierr)

    if(rank .gt. 0) allocate(self%coeff(self%n1p,self%n2p), stat=status)
    call MPI_Bcast(self%coeff,splsize,MPI_DOUBLE_PRECISION,0,comm,ierr)

    ! Allocate and bcast 1D spline data
    if(rank .gt. 0) allocate(self%pos1(self%n1p), stat=status)
    call MPI_Bcast(self%pos1,self%n1p,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(rank .gt. 0) allocate(self%pos2(self%n2p), stat=status)
    call MPI_Bcast(self%pos2,self%n2p,MPI_DOUBLE_PRECISION,0,comm,ierr)
    ! Knots
    if(rank .gt. 0) allocate(self%knot1(self%n1p+self%nord), stat=status)
    call MPI_Bcast(self%knot1,self%n1p+self%nord,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(rank .gt. 0) allocate(self%knot2(self%n2p+self%nord), stat=status)
    call MPI_Bcast(self%knot2,self%n2p+self%nord,MPI_DOUBLE_PRECISION,0,comm,ierr)
#endif
  end subroutine bcast_spl2d

  !!----------------------------------------------------------------------------------
  !! Geometry broadcast routines
  !!----------------------------------------------------------------------------------

  subroutine bcast_powcal(self,vecx,vecb,rank,comm)

    !! Arguments
    type(geobjlist_t), intent(inout) :: self   !< object data structure
    type(posveclis_t), intent(inout) :: vecx,vecb !<
    integer(ki4) , intent(in) :: rank,comm !<

    call bcast_geobjlist(self,rank,comm)

    call bcast_position(vecx,rank,comm)

    call bcast_position(vecb,rank,comm)

  end subroutine bcast_powcal


  !> Broadcast vtk geobjlist coordinate data
  subroutine bcast_geobjlist(self,rank,comm)

    !! Arguments
    type(geobjlist_t), intent(inout) :: self   !< data structure
    integer(ki4) , intent(in) :: rank,comm !<

    !! Variables
    integer(ki4) :: ierr, size, num !<
#ifdef WITH_MPI
    call MPI_Bcast(self%ngtype,1,MPI_INTEGER,0,comm,ierr)

    ! Position vectors
    call bcast_position(self%posl,rank,comm)

    ! list of nodes
    call MPI_Bcast(self%nnod,1,MPI_INTEGER,0,comm,ierr)
    if (self%nnod .gt. 0) then
       if(rank .gt.0) allocate(self%nodl(self%nnod), stat=status)
       call MPI_Bcast(self%nodl,self%nnod,MPI_INTEGER,0,comm,ierr)
    endif
    call MPI_Bcast(self%np,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nwset,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%ng,1,MPI_INTEGER,0,comm,ierr)

    ! obj
    if (self%ngtype==1.OR.self%nwset/=0) then
       if(rank .gt.0) allocate(self%obj(self%ng), stat=status)
       ! obj contains only 4 byte real weight
       call MPI_Bcast(self%obj%weight,self%ng,MPI_REAL,0,comm,ierr) ! alloc?
    endif
    ! nodl and obj2
    if (self%ngtype/=1) then
       if(rank .gt.0) allocate(self%obj2(self%ng), stat=status)
       ! In previous versions an MPI_TYPE was created here however this is simpler
       size = sizeof(self%obj2(1)) ! HJL !*WARNING*! NOT Fortran standard
       call MPI_Bcast(self%obj2,size*self%ng,MPI_CHARACTER,0,comm,ierr)
    endif

    call MPI_Bcast(self%ngunassigned,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nsampl,1,MPI_INTEGER,0,comm,ierr)
    ! These 3 broadcasts are declared ki2par which is a 2 byte integer
    call MPI_Bcast(self%numnparam,1,MPI_INTEGER2,0,comm,ierr)
    call MPI_Bcast(self%numrparam,1,MPI_INTEGER2,0,comm,ierr)
    num = self%numnparam
    call MPI_Bcast(self%nparam,num,MPI_INTEGER2,0,comm,ierr)
    num = self%numrparam
    call MPI_Bcast(self%rparam,num,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%coordbb,6,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%binbb,6,MPI_REAL,0,comm,ierr)

    ! Several broadcasts for the powcal code broadcasts were removed from here
    ! as they break the hdsgen broadcasts in the single exe version and appear unnecessary for powcal

    if(rank .gt. 0) call geobjlist_addcube(self)
#endif
  end subroutine bcast_geobjlist

  subroutine bcast_position(self,rank,comm)

    type(posveclis_t), intent(inout) :: self !< position list data
    integer(ki4) , intent(in) :: rank,comm !<

    integer(ki4) :: ierr, i, num !<
#ifdef WITH_MPI
    ! posvec is an array of 3 element 4 byte reals
    ! the position structure contains this an integer from position_h
    ! and since the addition of skylights some other info
    call MPI_Bcast(self%np,1,MPI_INTEGER,0,comm,ierr)
    if(rank .gt. 0 .and. .not. allocated(self%pos)) allocate(self%pos(self%np), stat=status)
    call MPI_Bcast(self%pos(1)%posvec(1),3*self%np,MPI_REAL,0,comm,ierr)
    ! Other bits
    ! These 3 broadcasts are declared ki2par which is a 2 byte intege
    call MPI_Bcast(self%numnparpos,1,MPI_INTEGER2,0,comm,ierr)
    num = self%numnparpos
    call MPI_Bcast(self%nparpos,num,MPI_INTEGER2,0,comm,ierr)
    call MPI_Bcast(self%numrparpos,1,MPI_INTEGER2,0,comm,ierr)
    num = self%numrparpos
    call MPI_Bcast(self%rparpos,num,MPI_REAL,0,comm,ierr)
#endif
  end subroutine bcast_position

  subroutine bcast_hds(geobjl,numerics,btree,rank,comm)

    !! Arguments
    type(geobjlist_t), intent(inout) :: geobjl !< point data
    type(numerics_t), intent(inout) :: numerics !< geobjlist controls
    type(btree_t), intent(inout) :: btree !< btree data
    integer(ki4) , intent(in) :: rank,comm !<

    !! Variables
    integer :: ierr !<

    call bcast_numerics(numerics,geobjl,rank,comm)

    call bcast_btree(btree,rank,comm)

  end subroutine bcast_hds

  subroutine bcast_btree(self,rank,comm)

    !! Arguments
    type(btree_t), intent(inout) :: self !< btree data
    integer(ki4), intent(in) :: rank,comm !<

    !! Variables
    integer :: ierr !<
#ifdef WITH_MPI
    call MPI_Bcast(self%nt,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%npter,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%ndepth,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nttype,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%mtype,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nttalg,1,MPI_INTEGER,0,comm,ierr)
    ! These three are declared ki2 which is a 4 byte integer
    call MPI_Bcast(self%nxyz,3,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%hxyz,3,MPI_REAL,0,comm,ierr)
    if(self%nt>0.and.self%npter>0.and.rank>0) then
       allocate(self%pter(self%npter,self%nt),stat=status)
    endif
    if(self%nt>0.and.rank>0) then
       allocate(self%desc(2,self%nt),self%corner(3,self%nt),stat=status)
    endif
    call MPI_Bcast(self%pter,self%nt * self%npter,MPI_INTEGER,0,comm,ierr)
    ! These three are declared ki2 which is a 4 byte integer
    call MPI_Bcast(self%desc,2*self%nt,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%corner,3*self%nt,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nexten,1,MPI_INTEGER,0,comm,ierr)
    if(self%nexten>0) then
       if(rank>0) then
          allocate(self%exten(3,self%nexten),stat=status)
       endif
       call MPI_Bcast(self%exten,3*self%nexten,MPI_INTEGER,0,comm,ierr)
    endif
    call MPI_Bcast(self%nroot,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%maxchildn,1,MPI_INTEGER,0,comm,ierr)

    call bcast_ls(self%objectls,rank,comm)
#endif
  end subroutine bcast_btree

  subroutine bcast_ls(self,rank,comm)

    !! arguments
    type(ls_t), intent(inout) :: self   !< ls data
    integer(ki4) , intent(in) :: rank,comm !<

    !! Variables
    integer :: ierr,ilist,isizel,isizeh !<
#ifdef WITH_MPI
    call MPI_Bcast(self%nlist,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%n2,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%listype,1,MPI_INTEGER,0,comm,ierr)

    if(rank > 0) then
       ilist = self%nlist
       isizel = ilist
       isizeh = self%n2
       call ls_init(self,isizeh,isizel)
       self%nlist = ilist
    endif
    if(self%nlist>0) then
       call MPI_Bcast(self%list(1,2),self%nlist,MPI_INTEGER,0,comm,ierr)
    endif
#endif

  end subroutine bcast_ls

  subroutine bcast_numerics(self,geobjl,rank,comm)

    !! Arguments
    type(numerics_t), intent(inout) :: self !< local variable
    type(geobjlist_t), intent(inout) :: geobjl !< point data
    integer(ki4) , intent(in) :: rank,comm !<

    integer :: ierr !<
    integer :: type !<
#ifdef WITH_MPI
    call MPI_Bcast(self%nsize,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nsizep,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nsizee,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nsizeh,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nsizel,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%ndepth,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nttype,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nttalg,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nxyz,3,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%hxyz,3,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%ngeobj,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%mtype,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%geomtype,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%mingeobjinbin,1,MPI_INTEGER,0,comm,ierr)
    ! nquante is declared ki2 which is a 4 byte integer
    call MPI_Bcast(self%nquante,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%maxtolerance,1,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%mintolerance,1,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%nbdcub,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%cornflag,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%lowcorner,3,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%upcorner,3,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%dilen,1,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%dolen,1,MPI_REAL,0,comm,ierr)

    ! tfmdata
    call MPI_Bcast(self%position_coord_tfm%scale,3,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%position_coord_tfm%offset,3,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%position_coord_tfm%matrix,9,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%position_coord_tfm%ntfm,1,MPI_INTEGER,0,comm,ierr)
    if(rank > 0) geobjl%tfmdata=self%position_coord_tfm

    ! quantfm data
    call MPI_Bcast(self%geobj_coord_tfm%hmin,3,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%geobj_coord_tfm%rhmin,3,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%geobj_coord_tfm%offvec,3,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(self%geobj_coord_tfm%nqtfm,1,MPI_INTEGER,0,comm,ierr)
    if(rank > 0) geobjl%quantfm=self%geobj_coord_tfm

    ! boubo
    call MPI_Bcast(geobjl%coordbb,6,MPI_REAL,0,comm,ierr)
    call MPI_Bcast(geobjl%binbb,6,MPI_REAL,0,comm,ierr)
#endif

  end subroutine bcast_numerics

  subroutine bcast_skyl(self,rank,comm)

    !! Arguments
    type(skyl_t), intent(inout) :: self !< local variable
    integer(ki4) , intent(in) :: rank,comm !<

    integer(ki4) :: ierr,ibdim,indbox(2) !<

#ifdef WITH_MPI
    ! numerics structure
    call bcast_sknumerics(self%n,rank,comm)

    ! skylight debug
    call MPI_Bcast(self%debug,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%ndskyl,10,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%ndskyln,1,MPI_INTEGER,0,comm,ierr)

    call MPI_Bcast(self%geoml,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%psilts,8,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%psidelta,4,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%geoml,1,MPI_DOUBLE_PRECISION,0,comm,ierr)

    call MPI_Bcast(self%dimbox,4,MPI_INTEGER,0,comm,ierr)

    ibdim = 2-mod(self%n%skyltyp,2)
    if (ibdim==2) then
       indbox(1)=max( self%dimbox(1,1),self%dimbox(1,2) )
       indbox(2)=max( self%dimbox(2,1),self%dimbox(2,2) )
    else
       indbox=self%dimbox(:,1)
    end if

    if(any(indbox < 1) .or. ibdim < 1) return

    if(rank > 0) then
       allocate(self%inboxr(indbox(1),ibdim),self%inboxz(indbox(1),ibdim), &
            self%ouboxr(indbox(2),ibdim),self%ouboxz(indbox(2),ibdim),&
            stat=status)
    endif

    call MPI_Bcast(self%inboxr,indbox(1)*ibdim,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%inboxz,indbox(1)*ibdim,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%ouboxr,indbox(2)*ibdim,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%ouboxz,indbox(2)*ibdim,MPI_DOUBLE_PRECISION,0,comm,ierr)
#endif
  end subroutine bcast_skyl

  subroutine bcast_sknumerics(self,rank,comm)

    !! Arguments
    type(sknumerics_t), intent(inout) :: self !< local variable
    integer(ki4), intent(in) :: rank,comm !<

    integer(ki4) :: ierr !<
#ifdef WITH_MPI
    call MPI_Bcast(self%psiexts,4,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%extsopt,4,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%nexts,4,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%binsize,4,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%extsdel,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%ngeoml,4,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%skyltyp,4,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%lrext,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(self%lzext,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(self%rextmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%rextmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%zextmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%zextmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%toli,3,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%div,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%control,3,MPI_INTEGER,0,comm,ierr)
#endif

  end subroutine bcast_sknumerics

  subroutine bcast_dnumerics(self,rank,comm)

    !! Arguments
    type(dnumerics_t), intent(inout) :: self !< local variable
    integer(ki4), intent(in) :: rank,comm !<

    integer :: ierr !<

#ifdef WITH_MPI
    call MPI_Bcast(self%tfm,80,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(self%npos,1,MPI_INTEGER,0,comm,ierr)
    if(rank .gt. 0) then
       allocate(self%r(self%npos),self%z(self%npos),stat=status)
    endif
    call MPI_Bcast(self%r,self%npos,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%z,self%npos,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%div,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(self%ldiv,1,MPI_INTEGER,0,comm,ierr)
    ! These 3 broadcasts are declared ki2par which is a 2 byte integer
    call MPI_Bcast(self%descode,1,MPI_INTEGER2,0,comm,ierr)
    call MPI_Bcast(self%csys,1,MPI_INTEGER2,0,comm,ierr)
    call MPI_Bcast(self%cunits,1,MPI_INTEGER2,0,comm,ierr)
    call MPI_Bcast(self%stang,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%finang,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%stpos,3,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(self%finpos,3,MPI_DOUBLE_PRECISION,0,comm,ierr)
#endif
  end subroutine bcast_dnumerics

end module smitermpi_m
