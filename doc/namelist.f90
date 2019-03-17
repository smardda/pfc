program namelist_p

  use const_kind_m
  use const_numphys_h
  use date_time_m
  use namelist_d

  implicit none


! Local variables
  character(*), parameter :: m_name='namelist' !< module name
  type(dfiles_t)     :: file      !< names of files
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  integer(ki4):: nread !< unit for dat files
  integer(ki4):: nin !< unit for other data
  integer(ki4):: nscal !< number of scalars (body identifiers)
  integer(ki4):: i !< loop variable
  integer(ki4):: j !< loop variable
  integer(ki4) :: islen   !< length of input field filename
  integer(ki4), dimension(2):: idum !< dummy array
  logical :: iltest !< logical flag

!! start log
  call log_init(fileroot,timestamp)
!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(30)
  call clock_start(1,'namelist run time')
!--------------------------------------------------------------------------
!! print header

  print *, '----------------------------------------------------'
  print *, 'namelist: dummy program'
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long
!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run namelist type at the command line:'
     print *, '   namelist fileroot [opt]'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
!! strip any final '.dat' string
     islen=len_trim(fileroot)
     if (islen>4) then
        if (fileroot(islen-3:islen)=='.dat') then
           fileroot(islen-3:islen)='    '
        end if
     end if
!--------------------------------------------------------------------------
!! cleanup and closedown

  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete
!--------------------------------------------------------------------------

end program namelist_p
