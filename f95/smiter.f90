!> Runs the SMITER workflow and coordinates the MPI decomposition
program smiter

  use geoq_p,        only: geoq_sub
  use hdsgen_p,      only: hdsgen_sub
  use powcal_p,      only: powcal_sub
  use smitermpi_m,   only: smitermpi_init, smiter_wctime, &
                           smitermpi_barrier, smitermpi_finish
  use smitermpi_h,   only: decomp_t


  implicit none

  type(decomp_t) :: decomp
  double precision :: smiterstarttime,starttime
  double precision :: shadtime,restime,hdstime,powtime

  ! Initialize MPI environment and assign work to ranks
  call smitermpi_init(decomp)

  smiterstarttime = smiter_wctime()

  ! GEOQ shadow nad results are executed by their assigned rank
  if(decomp%myrank == decomp%procgeoqshad) then
     starttime = smiter_wctime()
     call geoq_sub(1,.true.,decomp)
     shadtime = smiter_wctime() - starttime
  end if

  if(decomp%myrank == decomp%procgeoqres) then
     starttime = smiter_wctime()
     call geoq_sub(2,.true.,decomp)
     restime = smiter_wctime() - starttime
  end if

  ! HDSGEN is executed by its dedicated rank
  if(decomp%myrank == decomp%prochdsgen) then
     starttime = smiter_wctime()
     call hdsgen_sub(3,.true.,decomp)
     hdstime = smiter_wctime() - starttime
  end if

  ! Al ranks in the POWCAL communicator participate in POWCAL
  if(decomp%powcal) then
     starttime = smiter_wctime()
     call powcal_sub(4,.true.,decomp)
     powtime = smiter_wctime() - starttime
  end if

  if(decomp%myrankpowcal == 1) &
     print '(a," = ",t50,g12.5,1x,a)', &
         "powcal processes available after WC-time", smiter_wctime()-smiterstarttime, "s"

  ! Synchronize all ranks before reporting timings
  call smitermpi_barrier(decomp%commsmiter)

  if(decomp%myrank==0) &
     print '(a," = ",t50,g12.5,1x,a)',  "smiter runtime ", smiter_wctime() - smiterstarttime, "s"
  if(decomp%myrank == decomp%procgeoqshad)  &
     print '(a," = ",t50,g12.5,1x,a)', "geoq shadow wall-clock time (WCT)", shadtime,"s"
  if(decomp%myrank == decomp%procgeoqres) &
     print '(a," = ",t50,g12.5,1x,a)',  "geoq results wall-clock time (WCT)", restime,"s"
  if(decomp%myrank == decomp%prochdsgen) &
     print '(a," = ",t50,g12.5,1x,a)', "hdsgen wall-clock time (WCT)", hdstime,"s"
  if(decomp%powcal .and. decomp%myrankpowcal == 0) &
     print '(a," = ",t50,g12.5,1x,a)', "powcal total wall-clock time (WCT)", powtime,"s"
  
  ! Ensure all ranks have completed reporting before finalizing MPI
  call smitermpi_barrier(decomp%commsmiter)

  !--------------------------------------------------------------------------
  !! Finalize MPI
  call smitermpi_finish(decomp)
  !--------------------------------------------------------------------------

end program smiter
