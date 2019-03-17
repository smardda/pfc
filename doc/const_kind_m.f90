module const_kind_m

  implicit none
  private

!> 2 byte integer kind
  integer, parameter, public :: ki2 = selected_int_kind(4) 
!> number of bits in 2 byte integer kind
  integer, parameter, public :: ki2bits = 16 

!> 4 byte integer kind
  integer, parameter, public :: ki4 = selected_int_kind(9) 

!> 8 byte integer kind
  integer, parameter, public :: ki8 = selected_int_kind(18) 

!> 4 byte real kind
  integer, parameter, public :: kr4 = selected_real_kind(6) 

!> 8 byte real kind
  integer, parameter, public :: kr8 = selected_int_kind(15) 

!> real output format
  character(len=5), parameter, public :: cfmtr = 'G15.8'
!character(len=5), parameter, public :: cfmtr = 'G12.5'

!! derived parameters
  character(len=7), parameter, public :: cfmtbs = '('//cfmtr//')'
  character(len=12), parameter, public :: cfmtbs1 = '('//cfmtr//',1X)'
  character(len=13), parameter, public :: cfmtb1v = '(3(1X,'//cfmtr//'))'
  character(len=13), parameter, public :: cfmtbv1 = '(3('//cfmtr//',1X))'
  character(len=12), parameter, public :: cfmt1v = '3(1X,'//cfmtr//'))'
  character(len=12), parameter, public :: cfmtv1 = '3('//cfmtr//',1X))'
  character(len=12), parameter, public :: cfmt2v = '4(2X,'//cfmtr//'))'
  character(len=15), parameter, public :: cfmtbnv = '(999(1X,'//cfmtr//'))'
  character(len=14), parameter, public :: cfmtnv = '999(1X,'//cfmtr//'))'

end module const_kind_m
