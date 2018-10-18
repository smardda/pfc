module datline_m

  use const_kind_m
  use log_m
  use datline_h

  implicit none
  private

! public subroutines
  public :: &
  datline_init,  & !< initialise unit for reading
  datline_read,  & !< read logical line of dat file
  datline_fieldtype,  & !< establish field type in logical line of dat file
  datline_chkcomma, & !< check whether commas in line, effectively in the numeric data
  datline_split !< split logical line of dat file into chars

! public types

! private variables
  character(*), parameter :: m_name='datline_m' !< module name
  character(len=80) :: ibuff !< buffer for input/output
  character(len=8) :: blank8='        ' !< 8 blanks
  character(len=16) :: blank16='                ' !< 16 blanks
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer :: status  !< status flag
  logical :: iltest !< logical flag

  integer, save :: nread  !< file unit for input

  contains
!---------------------------------------------------------------------
!> initialise unit for reading
subroutine datline_init(self,kread)

  !! arguments
  type(datline_t), intent(inout) :: self   !< datline data
  integer, intent(inout) :: kread   !< unit number

  !! local
  character(*), parameter :: s_name='datline_init' !< subroutine name
  if (kread<=0) then
     call log_error(m_name,s_name,1,error_warning,'Negative or zero unit number')
  end if
  nread=kread

end subroutine datline_init
!---------------------------------------------------------------------
!> read logical line of dat file
subroutine datline_read(self,kbuff,kread)

  !! arguments
  type(datline_t), intent(out) :: self   !< datline data
  character(len=80), intent(in) :: kbuff   !< datline data
  integer, intent(inout), optional :: kread   !< unit number

  !! local
  character(*), parameter :: s_name='datline_read' !< subroutine name
  integer  :: iread !< actual unit number
  integer(ki4)  :: ilen !< length of character string
  integer(ki4)  :: ipos1!< positions in self
  integer(ki4)  :: ipos2 !< positions in self

  if(present(kread)) then
     iread=kread
  else
     iread=nread
  end if
  !! fill from buffered line
  self%line=' '
  self%line(1:80)=kbuff
  self%inline=1
  ilen=len_trim(kbuff)
  if (ilen>72) then
     !! read continuation lines from file
     ipos1=81
     do i=2, self%nlines
        read(iread,fmt='(a80)',iostat=status) ibuff
        call log_read_check(m_name,s_name,1,status)
        self%inline=self%inline+1
        ipos2=ipos1+79
        self%line(ipos1:ipos2)=ibuff
        ilen=len_trim(ibuff)
        if (ilen<=72) exit
        ipos1=ipos2+1
     end do
  end if

end subroutine datline_read
!---------------------------------------------------------------------
!> establish field type in logical line of dat file
!> subroutine datline_fieldtype(self,kctype)
subroutine datline_fieldtype(self)

  !! arguments
  type(datline_t), intent(inout) :: self   !< datline data
  !     character(len=*), intent(inout), optional :: kctype   !< type of field

  !! local
  character(len=*), parameter :: s_name='datline_fieldtype' !< subroutine name
  character(len=5)  :: ictype !< actual type string
  character(len=8)  :: field!< 1st field
  character(len=8)  :: fieldin !< 1st field
  integer(ki4)  :: ilen !< length of character string
  integer(ki4)  :: icomma !< is comma present

  fieldin=self%line(1:8)
  icomma=index(fieldin,',')
  if (icomma>0) then
     call log_error(m_name,s_name,1,error_fatal,'Free format dat file not processed')
  end if
  field=adjustl(fieldin)
  ilen=len_trim(field)
  if (field(ilen:ilen)=='*') then
     ictype='large'
  else
     ictype='small'
  end if

  !      if(present(kctype)) then
  !    !! output as argument
  !         kctype=ictype
  !      else
  !! save to self
  self%type=ictype
  !      end if

end subroutine datline_fieldtype
!---------------------------------------------------------------------
!> check whether commas in line, effectively in the numeric data
subroutine datline_chkcomma(self)

  !! arguments
  type(datline_t), intent(inout) :: self   !< datline data

  !! local
  character(len=*), parameter :: s_name='datline_chkcomma' !< subroutine name
  character(len=320)  :: fieldin  !<  field
  character(len=320)  :: field  !<  field
  integer(ki4)  :: ilen !< length of character string
  integer(ki4)  :: icomma !< is comma present

  fieldin=self%line
  field=adjustl(fieldin)
  ilen=len_trim(field)
  icomma=index(fieldin(9:ilen),',')
  if (icomma>0) then
     call log_error(m_name,s_name,1,error_warning,'Unexpected commas in input line')
  end if

end subroutine datline_chkcomma
!---------------------------------------------------------------------
!> split logical line of dat file
subroutine datline_split(self)

  !! arguments
  type(datline_t), intent(inout) :: self   !< datline data

  !! local
  character(len=*), parameter :: s_name='datline_split' !< subroutine name
  character(len=8)  :: field!< 1st field
  character(len=8)  :: fieldin !< 1st field
  integer(ki4)  :: ilen !< length of character string
  integer(ki4)  :: i8!< counters for 8 character fields
  integer(ki4)  :: ii8 !< counters for 8 character fields
  integer(ki4)  :: i16!< counters for 16 character fields
  integer(ki4)  :: ii16 !< counters for 16 character fields

  if (self%type=='large') then
     self%n8=1
     self%n16=4*self%inline
  else
     self%n8=1+8*self%inline
     self%n16=0
  end if

  allocate(self%flds8(self%n8), stat=status)
  if (self%n16>0) allocate(self%flds16(self%n16), stat=status)

  !! set entry name
  fieldin=self%line(1:8)
  field=adjustl(fieldin)
  ilen=len_trim(field)
  if (self%type=='large') then
     field(ilen:ilen)=' '
  end if
  self%flds8(1)=field

  if (self%type=='large') then
     i16=1
     ii16=1
     do j=1,self%inline
        do i=1,4
           self%flds16(i16)=self%line(16*ii16-7:16*ii16+8)
           if (self%flds16(i16)==blank16) then
              self%flds16(i16)(16:16)='0'
           end if
           i16=i16+1
           ii16=ii16+1
        end do
        ! skip 16 chars at end of line
        ii16=ii16+1
     end do
  else
     i8=2
     ii8=2
     do j=1,self%inline
        do i=1,8
           self%flds8(i8)=self%line(8*ii8-7:8*ii8+8)
           if (self%flds8(i8)==blank8) then
              self%flds8(i8)(8:8)='0'
           end if
           i8=i8+1
           ii8=ii8+1
        end do
        ! skip 16 chars at end of line
        ii8=ii8+2
     end do
  end if

end subroutine datline_split

end module datline_m
