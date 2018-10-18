module const_kind_m

  implicit none
  private

!> 2 byte integer kind
  integer, parameter, public :: ki2par = selected_int_kind(4) 
!> not 2 byte integer kind
  integer, parameter, public :: ki2 = selected_int_kind(9) 
!> number of bits in not 2 byte integer kind
  integer, parameter, public :: ki2bits = 32 

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
  character(len=5), parameter, public :: cfmte = 'E15.8'

!! derived parameters
  character(len=7), parameter, public :: cfmtbs = '('//cfmtr//')'
  character(len=12), parameter, public :: cfmtbs1 = '('//cfmtr//',1X)'
  character(len=13), parameter, public :: cfmtb1v = '(3(1X,'//cfmtr//'))'
  character(len=13), parameter, public :: cfmtbv1 = '(3('//cfmtr//',1X))'
  character(len=13), parameter, public :: cfmtbv2 = '(6('//cfmtr//',1X))'
  character(len=12), parameter, public :: cfmt1v = '3(1X,'//cfmtr//'))'
  character(len=12), parameter, public :: cfmtv1 = '3('//cfmtr//',1X))'
  character(len=12), parameter, public :: cfmt2v = '4(2X,'//cfmtr//'))'
  character(len=15), parameter, public :: cfmtbnv = '(999(1X,'//cfmtr//'))'
  character(len=14), parameter, public :: cfmtnv = '999(1X,'//cfmtr//'))'
  character(len=14), parameter, public :: cfmte1v = '3(1X,1P'//cfmte//'))'

end module const_kind_m
