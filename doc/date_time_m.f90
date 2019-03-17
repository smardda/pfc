!>
!!
!! Date and time in different formats
!!
!! \Bold{Type}
!!
!! Infrastructure;  Object oriented module.
!!
!! \Bold{Purpose}
!!
!! Provides a type containing variables for a short and long version of the
!! date and time and the subroutines to fill them.  The long type is
!! intended for anotating output files and the short type for generating unique
!! temporary file names.
!!
!! \Bold{Function}
!!
!! The module provides the routines \Tt{date_time_init} to initialise a
!! \Tt{date_time_t} type using the intrinsic date_and_time function and
!! \Tt{date_time_delete} to destroy the structure after use.
!<
module date_time_m

  implicit none
  private

! public parameters

  integer, public, parameter :: date_time_len_short = 15 !< local variable
  integer, public, parameter :: date_time_len_long = 32 !< local variable

! public types

!> Containing structure for run time strings.
  type, public :: date_time_t
     character(date_time_len_short) :: short !< local variable
     character(date_time_len_long)  :: long !< local variable
     integer, dimension(8) :: values !< local variable
  end type date_time_t

  public :: date_time_init
  public :: date_time_delete
  public :: date_time_long
  public :: date_time_short

!> Month array for long version of date
  character(9), dimension(12), parameter :: month = & !< local variable
 &(/ 'January  ', 'February ', &
 &'March    ', 'April    ', &
 &'May      ', 'June     ', &
 &'July     ', 'August   ', &
 &'September', 'October  ', &
 &'November ', 'December ' /)


  contains
subroutine date_time_init(self)

  type(date_time_t), intent(out) :: self !< local variable

  ! local variables
  character(8)  :: date !< local variable
  character(10) :: time !< local variable
  character(5)  :: zone !< local variable

  integer :: i !< loop variable

  call date_and_time(date, time, zone, self%values)

  write(self%short,'(A6,A6,A3)') date(3:8), time(1:6), time(8:10)

  write(self%long,'(2(I2.2,":"),I2.2," ",I0," ",A," ",I4)') &
 &self%values(5), self%values(6), self%values(7), self%values(3),&
 &trim(month(self%values(2))), self%values(1)

  do i = 1, 12
     if (self%short(i:i) == ' ') self%short(i:i) = '0'
  end do

end subroutine date_time_init
subroutine date_time_delete(self)

  type(date_time_t), intent(inout) :: self !< local variable

  self%values = 0

end subroutine date_time_delete


!> Convenience routine that returns a string containing the current
!! date and time (long format) without the need to initialise a type
function date_time_long() result(rv)

  character(date_time_len_long) :: rv !< local variable

  type(date_time_t) :: date_time !< local variable

  call date_time_init(date_time)
  rv = date_time%long
  call date_time_delete(date_time)

end function date_time_long


!> Convenience routine that returns a string containing the current
!! date and time (short format) without the need to initialise a type
function date_time_short() result(rv)

  character(date_time_len_short) :: rv !< local variable

  type(date_time_t) :: date_time !< local variable

  call date_time_init(date_time)
  rv = date_time%short
  call date_time_delete(date_time)

end function date_time_short

end module date_time_m
