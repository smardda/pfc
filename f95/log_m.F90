module log_m

  use const_kind_m
  use date_time_m

#ifdef WITH_MPI
  use mpi
#endif
  implicit none
  private

! public subroutines
  public :: &
  log_init,   & !< create log file and initialise
  log_error,  & !< add error message to log
  log_value,  & !< add name = value to log
  log_alloc_check,  & !< check array allocation
  log_open_check,  & !< check file opening
  log_read_check,  & !< check result of read
  log_write_check,  & !< check result of write
  log_getunit,  & !< get unit number for logging (namelists)
  log_close     !< close log

  interface log_value
   module procedure log_value_ki2
   module procedure log_value_ki4
   module procedure log_value_kr4
   module procedure log_value_kr8
   module procedure log_value_char
   module procedure log_value_kl
  end interface

! public variables
  integer, public, parameter :: error_fatal=0 !< public variable
  integer, public, parameter :: error_serious=1 !< public variable
  integer, public, parameter :: error_warning=2 !< public variable
  integer, public, parameter :: error_info=3 !< public variable
  integer, public, parameter :: error_debug=4 !< public variable
  integer, public, parameter :: log_info=5 !< public variable


! private variables
  type(date_time_t) :: timedate !< timestamp
  character(*), parameter :: m_name='log_m' !< module name
  character(*), dimension(0:5), parameter :: errorname = & !< type of "error"
 &(/'Fatal  ','Serious','Warning', 'Info   ', 'Debug  ','Info   '/)
  character(128) :: logfile !< name of error logging file
  integer :: nlog    !< error logging unit number
  integer(ki4) :: i       !< loop counter
  integer(ki4) :: errorno       !<  number of error messages
  integer(ki4) :: seriouserrors !< number of serious errors
  integer(ki4),parameter :: maxseriouserrors=10 !< max no of serious errors

  contains
!---------------------------------------------------------------------
!> create log file and initialise
subroutine log_init(fileroot,timestamp)

  !! arguments
  character(*), intent(inout) :: fileroot !<root of log file name
  type(date_time_t), intent(in) :: timestamp !< timestamp
  !! local
  logical :: unitused !< flag to test unit is available
  integer(ki4) :: ilen  !< length of string

#ifdef WITH_MPI
  integer rank, error
  character(5)  :: log_id
#endif
  !! initialise counters !mpi5!
  errorno=0
  seriouserrors=0

  !! create filename
  timedate=timestamp
  ! strip ctl from end of name
  ilen=len_trim(fileroot)
  if (ilen>4) then
     if (fileroot(ilen-3:ilen)=='.ctl') ilen=ilen-4
  end if
  logfile=fileroot(1:ilen)//'.log'
  ! and reset file root
  fileroot=logfile(1:ilen)

#ifdef WITH_MPI
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, error)
  if (rank .ne. 0) return
  !write(log_id, '(I5.5)'), rank
  !logfile = trim(logfile)//'-'//trim(log_id)
#endif
  !! get file unit !mpi6!
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        nlog=i
        exit
     end if
  end do
  !! using this implies circular module dependency
  !! call misc_getfileunit(nlog)

  open(unit=nlog,file=logfile,status='REPLACE')

  !! write initial header

  write(nlog,'(/,a,a,/)') 'LOG FILE: ',trim(logfile)
  write(nlog, '(a,/)') timedate%long


end subroutine log_init
!---------------------------------------------------------------------
!>write message to error and info log
subroutine log_error(modname,subname,point,severity,errormessage)

  !! arguments
  character(*), intent(in) :: modname  !< module name
  character(*), intent(in) :: subname  !< subprogram name
  integer(ki4), intent(in) :: point    !< calling point
  integer(ki4), intent(in) :: severity !< error severity
  character(*), intent(in) :: errormessage  !< error message

#ifdef WITH_MPI
  integer rank, error
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, error)
  if (rank .ne. 0) return
#endif
  !! local !mpir!

  if(severity<log_info) then
     errorno=errorno+1
     write(nlog, '(i7.7,a,i2,a)') &
 &   errorno, ': '//errorname(severity)//': '// trim(modname)//':'//trim(subname)//':', point, ':- '//errormessage
  else
     write(nlog, '("Log  ",a,i2,a)') &
 &   ': '//errorname(severity)//': '// trim(modname)//':'//trim(subname)//':', point, ':- '//errormessage
  end if

  !! fatal error
  if(severity==error_fatal) then
     !!close down
     seriouserrors=seriouserrors+1
     write(nlog,'(a)') 'FATAL ERROR - run terminated'
     call log_close
     stop 1

     !!serious error
  else if(severity==error_serious) then
     seriouserrors=seriouserrors+1
     if(seriouserrors>=maxseriouserrors)then
        !!close down
        write(nlog,'(a)') 'TOO MANY SERIOUS ERRORS - run terminated'
        call log_close
        stop 1

     end if


     !!error warning
  else

     !!error info

     !!error debug

     !! log_info

  end if

end subroutine log_error
!---------------------------------------------------------------------
!> check for errors in array allocation
subroutine log_alloc_check(modname,subname,point,status)

  ! arguments
  character(len=*), intent(in) :: modname  !< module name
  character(len=*), intent(in) :: subname  !< subprogram name
  integer, intent(in) :: point    !< calling point
  integer, intent(in) :: status   !< error status flag

  if(status/=0)then
     call log_error(modname,subname,point,error_fatal,'allocation failed')
     !DBG else !DBG
     !DBG write(*,*) trim(modname), trim(subname), point, 'allocation OK' !DBG
  end if

end subroutine log_alloc_check
!---------------------------------------------------------------------
!> check for errors in file opening
subroutine log_open_check(modname,subname,point,status)

  ! arguments
  character(len=*), intent(in) :: modname  !< module name
  character(len=*), intent(in) :: subname  !< subprogram name
  integer, intent(in) :: point    !< calling point
  integer, intent(in) :: status   !< error status flag

  if(status/=0)then
     call log_error(modname,subname,point,error_fatal,'file opening failed')
  end if


end subroutine log_open_check
!---------------------------------------------------------------------
!> check for errors in read
subroutine log_read_check(modname,subname,point,status)

  ! arguments
  character(len=*), intent(in) :: modname  !< module name
  character(len=*), intent(in) :: subname  !< subprogram name
  integer, intent(in) :: point    !< calling point
  integer, intent(in) :: status   !< error status flag

  if(status/=0)then
     call log_error(modname,subname,point,error_fatal,'read failed')
  end if


end subroutine log_read_check
!---------------------------------------------------------------------
!> check for errors in write
subroutine log_write_check(modname,subname,point,status)

  ! arguments
  character(len=*), intent(in) :: modname  !< module name
  character(len=*), intent(in) :: subname  !< subprogram name
  integer, intent(in) :: point    !< calling point
  integer, intent(in) :: status   !< error status flag

  if(status/=0)then
     call log_error(modname,subname,point,error_fatal,'write failed')
  end if


end subroutine log_write_check
!---------------------------------------------------------------------
!> get unit number for logging (namelists)
subroutine log_getunit(kunit)

  !! arguments
  integer, intent(out) :: kunit    !< log unit number

  kunit=nlog

end subroutine log_getunit
!---------------------------------------------------------------------
!> close log file
subroutine log_close

#ifdef WITH_MPI
  integer rank, error
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, error)
  if (rank .ne. 0) return
#endif
  !! arguments !mpir!

  write(nlog, '(//,a,/,a)') ' Error Summary ',' ------------- '
  write(nlog, '(a,i7)') ' total number of errors   = ',errorno
  write(nlog, '(a,i7)') ' number of serious errors = ',seriouserrors
  write(nlog,'(a,/)') date_time_long()
  write(nlog,'(a,/)') 'END OF LOG FILE'
  close(unit=nlog)

end subroutine log_close
!---------------------------------------------------------------------
!> log name=value
subroutine log_value_ki2(varname,value,units)
  !! arguments
  character(*), intent(in) :: varname  !< variable name
  integer(ki2par), intent(in) :: value    !< variable value
  character(*), intent(in),optional :: units  !< units name

#ifdef WITH_MPI
  integer rank, error
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, error)
  if (rank .ne. 0) return
#endif
  !! output !mpir!
  if(present(units)) then
     write(nlog, '("Log  : ",a," = ",i10,1x,a)') varname,value,units
  else
     write(nlog, '("Log  : ",a," = ",i10)') varname,value
  end if
end subroutine log_value_ki2
subroutine log_value_ki4(varname,value,units)
  !! arguments
  character(*), intent(in) :: varname  !< variable name
  integer(ki4), intent(in) :: value    !< variable value
  character(*), intent(in),optional :: units  !< units name

#ifdef WITH_MPI
  integer rank, error
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, error)
  if (rank .ne. 0) return
#endif
  !! output !mpir!
  if(present(units)) then
     write(nlog, '("Log  : ",a," = ",i10,1x,a)') varname,value,units
  else
     write(nlog, '("Log  : ",a," = ",i10)') varname,value
  end if
end subroutine log_value_ki4
subroutine log_value_kr4(varname,value,units)
  !! arguments
  character(*), intent(in) :: varname  !< variable name
  real(ki4), intent(in) :: value    !< variable value
  character(*), intent(in), optional :: units  !< units name
#ifdef WITH_MPI
  integer rank, error
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, error)
  if (rank .ne. 0) return
#endif
  !! output !mpir!
  if(present(units)) then
     write(nlog, '("Log  : ",a32,"  = ",g12.5,1x,a)') varname,value,units
  else
     write(nlog, '("Log  : ",a," = ",g12.5)') varname,value
  end if
end subroutine log_value_kr4
subroutine log_value_kr8(varname,value,units)
  !! arguments
  character(*), intent(in) :: varname  !< variable name
  real(ki8), intent(in) :: value    !< variable value
  character(*), intent(in),optional :: units  !< units name
#ifdef WITH_MPI
  integer rank, error
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, error)
  if (rank .ne. 0) return
#endif
  !! output !mpir!

  if(present(units)) then
     write(nlog, '("Log  : ",a32,"  = ",1pe15.8,2x,a)') varname,value,units
  else
     write(nlog, '("Log  : ",a," = ",g20.8)') varname,value
  end if
end subroutine log_value_kr8
subroutine log_value_kl(varname,value)
  !! arguments
  character(*), intent(in) :: varname  !< variable name
  logical, intent(in) :: value    !< variable value
#ifdef WITH_MPI
  integer rank, error
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, error)
  if (rank .ne. 0) return
#endif
  !! output !mpir!
  write(nlog, '("Log  : ",a," = ",l2)') varname,value
end subroutine log_value_kl
subroutine log_value_char(varname,value)
  !! arguments
  character(*), intent(in) :: varname  !< variable name
  character(*), intent(in) :: value    !< variable value
#ifdef WITH_MPI
  integer rank, error
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, error)
  if (rank .ne. 0) return
#endif
  !! output !mpir!
  write(nlog, '("Log  : ",a," = ",a)') varname,value
end subroutine log_value_char


end module log_m
