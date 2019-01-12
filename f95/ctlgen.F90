program ctlgen

  use const_kind_m
  use date_time_m
  use log_m
  use misc_m
  use ccontrol_m

  character(len=2) :: arg1 !< 1st argument
  character(len=1) :: typoutfile !< type of file to produce
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=80) :: code !< first code input file intended for
  character(len=80) :: coden !< other code input file intended for
  integer   :: inarg   !< number of command arguments
!--------------------------------------------------------------------------
!! use arguments
  inarg=command_argument_count()
  if (inarg<2) then
!! not enough arguments specified
     stop 1
  end if
  call date_time_init(timestamp)
!--------------------------------------------------------------------------
!! print header
  print *, '----------------------------------------------------'
  print *, 'ctlgen: generate ctl file for vtktfm'
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long

!!get type of file to produce
  call get_command_argument(1,value=arg1)
  typoutfile=arg1(2:2)
!! get input file root
  call get_command_argument(2,value=fileroot)
  call log_init(fileroot,timestamp)
!! for future developments, -x input might select program x of SMITER suite
!! definitely x will be determined by STA statement in .ctlin file
!! ccontrol_read/write should be part of xcontrol_m module
  call ccontrol_init(fileroot,code)
  coden=code

!--------------------------------------------------------------------------
!! do the work
  target_code: select case(code)
  case('vtktfm')
     call ccontrol_readctlin(coden)

     call ccontrol_writectl
  end select target_code

  print *, 'ctlgen completed successfully'
  call log_close

end program ctlgen
