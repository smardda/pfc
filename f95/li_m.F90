module li_m

  use const_kind_m
  use log_m

  implicit none
  private

! public subroutines
  public ::  &
  li_init,   & !< create simple linked list structure
  li_delete,   & !< delete simple linked list structure
  li_read,   & !< read simple linked list structure
  li_write, &  !< write simple linked list structure
  li_add,   & !< add to simple linked list structure
  li_addin,   & !< add into simple linked list structure
  li_count,   & !< count entries in simple linked list
  li_countpp,   & !< add one to number of entries in doubly linked list
  li_fetch, & !< return simple linked list
  li_write1, &  !< write simple linked list
  li_next1    !< return next entry in simple linked list

! public types
  type, public :: li_t
   !! linked list
     integer(ki4), dimension(:), allocatable  :: hoc   !< head of chain array
     integer(ki4), dimension(:), allocatable  :: last   !< last entry accessed of conten array
     integer(ki4), dimension(:), allocatable   :: conten   !< conten array
     integer(ki4) :: nhoc   !< number of entries in hoc array
     integer(ki4) :: nconten   !< number of entries in conten array
  end type li_t

!public variables

! private types

! private variables
  character(*), parameter :: m_name='li_m' !< module name
  integer   :: status   !< error status
  integer(ki4) :: inext !< next of simple linked list
  integer(ki4) :: i!< loop counters
  integer(ki4) :: j!< loop counters
  integer(ki4) :: k !< loop counters

  contains

!---------------------------------------------------------------------
! initialise simple linked list
subroutine li_init(self,ksizeh,ksize)

  !! arguments
  type(li_t), intent(out) :: self   !< linked list
  integer(ki4) , intent(in) :: ksizeh !< size of linked list hoc, set to unity
  integer(ki4) , intent(in) :: ksize !< size of linked list


  !! local
  character(*), parameter :: s_name='li_init' !< subroutine name
  integer(ki4) :: isize  !< local variable

  !! allocate storage
  if(ksize>0) then
     allocate(self%conten(ksize), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     self%nconten=0
  else
     call log_error(m_name,s_name,2,error_fatal,'No data')
  end if

  isize=ksizeh
  if(isize>0) then
     allocate(self%hoc(isize),self%last(isize), stat=status)
     call log_alloc_check(m_name,s_name,3,status)
  else
     call log_error(m_name,s_name,4,error_fatal,'No data')
  end if

  self%hoc=0
  self%last=0

end subroutine li_init
!---------------------------------------------------------------------
!> delete li_t
subroutine li_delete(self)

  !! arguments
  type(li_t), intent(inout) :: self !< particle coord data
  !! local

  deallocate(self%hoc)
  deallocate(self%last)
  deallocate(self%conten)

end subroutine li_delete
!---------------------------------------------------------------------
!> read in simple linked list data
subroutine li_read(self,kin)

  !! arguments
  type(li_t), intent(out) :: self   !< simple linked list data
  integer, intent(in) :: kin   !< input channel for simple linked list data


  !! local
  character(*), parameter :: s_name='li_read' !< subroutine name
  integer(ki4) :: inhoc  !< local variable
  integer(ki4) :: iconten  !< local variable

  ! read no in hoc array
  read(kin,*,iostat=status) inhoc
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading linked list nhoc')
  end if
  self%nhoc=inhoc

  ! read simple linked list hoc array

  read(kin,*,iostat=status) (self%hoc(j),j=1,inhoc)
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error reading linked list hoc')
  end if

  ! read no in conten array
  read(kin,*,iostat=status) iconten
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error reading linked list nconten')
  end if
  self%nconten=iconten

  ! read simple linked list conten array

  read(kin,*,iostat=status) (self%conten(j),j=1,iconten)
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error reading linked list conten')
  end if

end subroutine li_read
!---------------------------------------------------------------------
!> output simple linked list vector
subroutine li_write(self,kout)

  !! arguments
  type(li_t), intent(in) :: self   !< simple linked list data
  integer, intent(in) :: kout   !< output channel for simple linked list data


  !! local
  character(*), parameter :: s_name='li_write' !< subroutine name
  integer(ki4) :: inhoc  !< local variable
  integer(ki4) :: iconten  !< local variable

  ! output no in hoc array
  inhoc=self%nhoc
  write(kout,*,iostat=status) inhoc
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error writing linked list nhoc')
  end if

  ! output simple linked list hoc array

  write(kout,*,iostat=status) (self%hoc(j),j=1,inhoc)
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error writing linked list hoc')
  end if

  ! output no in conten array
  iconten=self%nconten
  write(kout,*,iostat=status) iconten
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error writing linked list nconten')
  end if

  ! output simple linked list conten array

  write(kout,*,iostat=status) (self%conten(j),j=1,iconten)
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error writing linked list conten')
  end if

end subroutine li_write
!---------------------------------------------------------------------
!> add to  linked list structure
subroutine li_add(self,khoc,knew)

  !! arguments
  type(li_t), intent(inout) :: self   !< linked list
  integer(ki4) , intent(inout) :: khoc !< hoc of linked list
  integer(ki4) , intent(in) :: knew !< new entry


  !! local
  character(*), parameter :: s_name='li_add' !< subroutine name
  integer(ki4) :: ihoc  !< local variable
  integer(ki4) :: icontsiz  !< local variable
  integer(ki4) :: iconten  !< local variable

  !! test storage
  ihoc=khoc
  iconten=self%nconten
  icontsiz=size(self%conten)
  if(ihoc<=icontsiz.AND.ihoc>=0) then
     if(knew<=icontsiz) then
        iconten=iconten+1
        if (ihoc==0) then
           !! new list
           self%conten(iconten)=0
           khoc=iconten
        else
           !! add to existing
           self%conten(knew)=khoc
           khoc=knew
        end if
        self%nconten=iconten
     else
        call log_error(m_name,s_name,1,error_fatal,'Contents array not large enough')
     end if
  else
     call log_error(m_name,s_name,2,error_fatal,'Hoc bound exceeded')
  end if

end subroutine li_add
!---------------------------------------------------------------------
!> add into linked list structure
subroutine li_addin(self,khoc,knew)

  !! arguments
  type(li_t), intent(inout) :: self   !< linked list
  integer(ki4) , intent(inout) :: khoc !< hoc of linked list
  integer(ki4) , intent(in) :: knew !< new entry

  !! local
  character(*), parameter :: s_name='li_addin' !< subroutine name
  integer(ki4) :: ihoc  !< local variable
  integer(ki4) :: icontsiz  !< local variable
  integer(ki4) :: iconten  !< local variable

  !! test storage
  ihoc=khoc
  iconten=self%nconten
  icontsiz=size(self%conten)
  if(ihoc<=icontsiz.AND.ihoc>=0) then
     if(knew<=icontsiz) then
        if (ihoc==0) then
           !! new list
           self%conten(knew)=0
           khoc=knew
        else
           !! add to existing
           self%conten(knew)=khoc
           khoc=knew
        end if
     else
        call log_error(m_name,s_name,1,error_fatal,'Contents array not large enough')
     end if
  else
     call log_error(m_name,s_name,2,error_fatal,'Hoc bound exceeded')
  end if

end subroutine li_addin
!---------------------------------------------------------------------
!> count entries in simple linked list
subroutine li_count(self,khoc,kn)

  !! arguments
  type(li_t), intent(in) :: self   !< linked list
  integer(ki4) , intent(in) :: khoc !< hoc of linked list
  integer(ki4) , intent(out) :: kn !< number in linked list

  !! local
  character(*), parameter :: s_name='li_count' !< subroutine name
  integer(ki4) :: iconten  !< local variable
  integer(ki4) :: inext  !< local variable
  integer(ki4) :: iold  !< local variable

  !! test storage
  iconten=self%nconten
  if(khoc<=iconten) then
     if (iconten==0.OR.khoc<=0) then
        !! empty list
        kn=0
     else
        !! count
        i=khoc
        kn=1
        do
           inext=self%conten(i)
           if (inext==0) exit
           kn=kn+1
           i=inext
        end do
     end if
  else
     call log_error(m_name,s_name,1,error_fatal,'Hoc bound exceeded')
  end if

end subroutine li_count
!---------------------------------------------------------------------
!> add one to number of entries in doubly linked list
subroutine li_countpp(self)

  !! arguments
  type(li_t), intent(inout) :: self   !< linked list

  !! local
  character(*), parameter :: s_name='li_countpp' !< subroutine name

  self%nconten=self%nconten+1

end subroutine li_countpp
!---------------------------------------------------------------------
!> return simple linked list
subroutine li_fetch(self,khoc,kout,knum)

  !! arguments
  type(li_t), intent(in) :: self   !< linked list
  integer(ki4), intent(in) :: khoc !< hoc of linked list
  integer(ki4), intent(inout) :: kout(*)   !< array large enough to contain list
  integer(ki4), intent(out) :: knum   !< number of proper entries in kout

  !! local
  character(*), parameter :: s_name='li_fetch' !< subroutine name
  integer(ki4) :: iconten  !< local variable
  integer(ki4) :: inext  !< local variable

  !! test storage
  iconten=self%nconten
  if(khoc<=iconten.AND.khoc>=0) then
     call li_count(self,khoc,knum)
     if (knum==0) then
        !! empty list
        return
     else
        !! add to array
        i=khoc
        j=knum
        kout(j)=i
     end if
     do
        inext=self%conten(i)
        if (inext==0) exit
        i=inext
        j=j-1
        kout(j)=i
     end do
  else
     call log_value('khoc',khoc)
     call log_error(m_name,s_name,2,error_fatal,'Hoc bound exceeded')
  end if

end subroutine li_fetch
!---------------------------------------------------------------------
!> write simple linked list
subroutine li_write1(self,khoc,kout)

  !! arguments
  type(li_t), intent(in) :: self   !< linked list
  integer(ki4) , intent(in) :: khoc !< hoc of linked list
  integer, intent(in) :: kout   !< output channel for simple linked list data

  !! local
  character(*), parameter :: s_name='li_write1' !< subroutine name
  integer(ki4) :: inum  !< local variable
  integer(ki4) :: iconten  !< local variable
  integer(ki4) :: inext  !< local variable
  integer(ki4), dimension(:), allocatable  :: iout  !< local variable

  !! test storage
  iconten=self%nconten
  if(khoc<=iconten.AND.khoc>=0) then
     call li_count(self,khoc,inum)
     if (inum==0) then
        !! empty list
        return
     else
        !! add to array
        allocate(iout(inum), stat=status)
        call log_alloc_check(m_name,s_name,1,status)
        i=khoc
        j=inum
        iout(j)=i
     end if
     do
        inext=self%conten(i)
        if (inext==0) exit
        i=inext
        j=j-1
        iout(j)=i
     end do
1        continue
  else
     call log_error(m_name,s_name,2,error_fatal,'Hoc bound exceeded')
  end if

  ! output simple linked list conten array
  write(kout,*,iostat=status) (iout(k),k=1,inum)
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error writing linked list conten')
  end if
  deallocate(iout)

end subroutine li_write1
!---------------------------------------------------------------------
!> return next entry in simple linked list
  integer(ki4) function li_next1(self,khoc)

!! arguments
  type(li_t), intent(inout) :: self   !< linked list
  integer(ki4) , intent(in) :: khoc !< hoc of linked list

!! local
  character(*), parameter :: s_name='li_next1' !< subroutine name
  integer(ki4) :: iconten  !< local variable
  integer(ki4) :: inext  !< local variable

!! test storage
  iconten=self%nconten
  if(khoc<=iconten.AND.khoc>=0) then
     i=khoc
     inext=0
     if (i<=0) then
!! empty list
     else
        j=self%last(1)
        if (j>0) then
! point to next
           inext=self%conten(j)
        else
! else point to head of chain
           inext=i
        end if
     end if
  else
     call log_error(m_name,s_name,1,error_fatal,'Hoc bound exceeded')
  end if

! output index and update
  li_next1=inext
  self%last(1)=inext

end function li_next1

end module li_m
