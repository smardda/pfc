program parse
  implicit none

  character(len=20) :: typetest !< type of test
  character(len=80) :: strtest !< string input variable
  character(len=80) :: strmin !< string input variable
  character(len=80) :: strmax !< string input variable
  integer   :: istatus   !< error status
  integer   :: inarg   !< number of command arguments
  integer   :: j     !< loop variable
  real   :: realno   !< input
  real   :: reallt   !< input range descriptor
!--------------------------------------------------------------------------
!! get arguments from arg
  inarg=command_argument_count()
  if (inarg<2) then
!! not enough arguments specified
     stop 1
  end if

!!get type of parser
  call get_command_argument(1,value=typetest)
!! get first argument to parse
  call get_command_argument(2,value=strtest)
  test_type: select case (typetest(1:6))
  case('number')
     read(strtest,*,iostat=istatus) realno
     if (istatus>0) then
        stop 1
     else
        test_value: select case (typetest(7:9))
        case('gt0')
           if (realno<=0) stop 1
        case('ge0')
           if (realno<0) stop 1
        end select test_value
     end if
  case('range_')
     read(strtest,*,iostat=istatus) realno
     if (istatus>0) then
        stop 1
     else
        if (inarg<3) stop 1
        call get_command_argument(3,value=strmin)
        read(strmin,*,iostat=istatus) reallt
        if (istatus>0) then
           stop 1
        else
           test_value2: select case (typetest(7:8))
           case('ge')
              if (realno<reallt) stop 1
           case('le')
              if (realno>reallt) stop 1
           case('gt')
              if (realno<=reallt) stop 1
           case('lt')
              if (realno>=reallt) stop 1
           end select test_value2
        end if
     end if
  case('vector')
     do j=2,inarg
        call get_command_argument(j,value=strtest)
        read(strtest,*,iostat=istatus) realno
        if (istatus>0) then
           stop 1
        else
           test_valuev: select case (typetest(7:9))
           case('gt0')
              if (realno<=0) stop 1
           case('ge0')
              if (realno<0) stop 1
           end select test_valuev
        end if
     end do
  end select test_type

end program parse
