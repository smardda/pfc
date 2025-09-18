program smiter

#ifdef WITH_MPI
  use mpi
#endif

  use date_time_m
  use log_m
  use clock_m
  use geoq_p
  use hdsgen_p
  use powcal_p
  use smitermpi_m
  use smitermpi_h


  implicit none

  type(decomp_t) :: decomp
  integer(ki4) :: ierr
  double precision :: smiterstarttime,starttime,shadtime,restime,hdstime,powtime
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
#ifdef WITH_MPI
  integer(ki4) :: ready, status(MPI_STATUS_SIZE)
#else
  double precision :: MPI_Wtime, dummyarg
  !MPI_Wtime(dummyarg)=0
  MPI_Wtime()=0
#endif


  ! Initialize MPI
  call smitermpi_init(decomp)

  ! Can't use the inbuilt clock as it affects the clocks in each unit
  smiterstarttime = MPI_Wtime()

  if(decomp%myrank == decomp%procgeoqshad) then
     starttime = MPI_Wtime()
     call geoq_sub(1,.True.,decomp)
     shadtime = MPI_Wtime() - starttime
  endif

  if(decomp%myrank == decomp%procgeoqres) then
     starttime = MPI_Wtime()
     call geoq_sub(2,.True.,decomp)
     restime = MPI_Wtime() - starttime
  endif

  if(decomp%myrank == decomp%prochdsgen) then
     starttime = MPI_Wtime()
     call hdsgen_sub(3,.True.,decomp)
     hdstime = MPI_Wtime() - starttime
  endif

  if(decomp%powcal) then
     starttime = MPI_Wtime()
     call powcal_sub(4,.True.,decomp)
     powtime = MPI_Wtime() - starttime
  endif

  if(decomp%myrankpowcal == 1) print*,"powcal processes available after",MPI_Wtime()-smiterstarttime,"s"

  call smitermpi_barrier(decomp%commsmiter)
  if(decomp%myrank==0) print*, "smiter runtime = ", MPI_Wtime() - smiterstarttime
  if(decomp%myrank == decomp%procgeoqshad) print*, "geoq shadow finished in",shadtime,"s"
  if(decomp%myrank == decomp%procgeoqres) print*, "geoq results finished in",restime,"s"
  if(decomp%myrank == decomp%prochdsgen) print*, "hdsgen finished in",hdstime,"s"
  if(decomp%powcal .and. decomp%myrankpowcal == 0) print*, "powcal finished in",powtime,"s"
  call smitermpi_barrier(decomp%commsmiter)

  !--------------------------------------------------------------------------
  !! Finalize MPI
  call smitermpi_finish(decomp)
  !--------------------------------------------------------------------------

end program smiter
