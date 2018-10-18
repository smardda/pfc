module ld_m

  use const_kind_m
  use log_m

  implicit none
  private

! public subroutines
  public ::  &
  ld_init,   & !< create doubly linked list structure
  ld_delete,   & !< delete doubly linked list structure
  ld_read,   & !< read doubly linked list structure
  ld_write, &  !< write doubly linked list structure
  ld_add,   & !< add to doubly linked list structure
  ld_addin,   & !< add into doubly linked list structure
  ld_remove,   & !< remove from doubly linked list structure
  ld_count,   & !< count entries in doubly linked list
  ld_countpp,   & !< add one to number of entries in doubly linked list
  ld_fetch, & !< return doubly linked list
  ld_write1, &  !< write doubly linked list
  ld_next1    !< return next entry in doubly linked list

! public types
  type, public :: ld_t
   !! linked list
     integer(ki4), dimension(:), allocatable  :: hoc   !< head of chain array
     integer(ki4), dimension(:), allocatable  :: last   !< last entry accessed of conten array
     integer(ki4), dimension(:,:), allocatable   :: conten   !< conten array
     integer(ki4) :: nhoc   !< number of entries in hoc array
     integer(ki4) :: nconten   !< number of entries in conten array (dynamic)
     integer(ki4) :: ndentry   !< number of data per entry including addressing
  end type ld_t

!public variables

! private types

! private variables
  character(*), parameter :: m_name='ld_m' !< module name
  integer   :: status   !< error status
  integer(ki4) :: inext !< next of doubly linked list
  integer(ki4) :: identry  !< size of one entry of list, including addresses
  integer(ki4) :: idata  !< number of data in one entry of list
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter

  contains

!---------------------------------------------------------------------
! initialise doubly linked list
subroutine ld_init(self,ksizeh,ksize,kdata)

  !! arguments
  type(ld_t), intent(out) :: self   !< linked list
  integer(ki4) , intent(in) :: ksizeh !< size of linked list hoc, set to unity
  integer(ki4) , intent(in) :: ksize !< size of linked list
  integer(ki4) , intent(in) :: kdata !< number of data per entry


  !! local
  character(*), parameter :: s_name='ld_init' !< subroutine name
  integer(ki4) :: isize  !< local variable

  !! allocate storage
  if(ksize>0.AND.kdata>0) then
     allocate(self%conten(2+kdata,ksize), stat=status)
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

  self%ndentry=2+kdata
  self%hoc=0
  self%last=0

end subroutine ld_init
!---------------------------------------------------------------------
!> delete ld_t
subroutine ld_delete(self)

  !! arguments
  type(ld_t), intent(inout) :: self !< particle coord data
  !! local

  deallocate(self%hoc)
  deallocate(self%last)
  deallocate(self%conten)

end subroutine ld_delete
!---------------------------------------------------------------------
!> read in doubly linked list data
subroutine ld_read(self,kin)

  !! arguments
  type(ld_t), intent(out) :: self   !< doubly linked list data
  integer, intent(in) :: kin   !< input channel for doubly linked list data


  !! local
  character(*), parameter :: s_name='ld_read' !< subroutine name
  integer(ki4) :: inhoc  !< local variable
  integer(ki4) :: iconten  !< local variable

  ! read no in hoc array
  read(kin,*,iostat=status) inhoc
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading linked list nhoc')
  end if
  self%nhoc=inhoc

  ! read doubly linked list hoc array

  read(kin,*,iostat=status) (self%hoc(j),j=1,inhoc)
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error reading linked list hoc')
  end if

  ! read no in conten array
  read(kin,*,iostat=status) iconten,identry
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error reading linked list nconten')
  end if
  self%nconten=iconten
  self%ndentry=identry

  ! read doubly linked list conten array

  read(kin,*,iostat=status) ((self%conten(k,j),k=1,identry),j=1,iconten)
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error reading linked list conten')
  end if

end subroutine ld_read
!---------------------------------------------------------------------
!> output doubly linked list vector
subroutine ld_write(self,kout)

  !! arguments
  type(ld_t), intent(in) :: self   !< doubly linked list data
  integer, intent(in) :: kout   !< output channel for doubly linked list data


  !! local
  character(*), parameter :: s_name='ld_write' !< subroutine name
  integer(ki4) :: inhoc  !< local variable
  integer(ki4) :: iconten  !< local variable

  ! output no in hoc array
  inhoc=self%nhoc
  write(kout,*,iostat=status) inhoc
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error writing linked list nhoc')
  end if

  ! output doubly linked list hoc array

  write(kout,*,iostat=status) (self%hoc(j),j=1,inhoc)
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error writing linked list hoc')
  end if

  ! output no in conten array
  iconten=self%nconten
  identry=self%ndentry
  write(kout,*,iostat=status) iconten,identry
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error writing linked list nconten')
  end if

  ! output doubly linked list conten array

  write(kout,*,iostat=status) ((self%conten(k,j),k=1,identry),j=1,iconten)
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error writing linked list conten')
  end if

end subroutine ld_write
!---------------------------------------------------------------------
!> add to  doubly linked list structure
subroutine ld_add(self,khoc,knew)

  !! arguments
  type(ld_t), intent(inout) :: self   !< linked list
  integer(ki4) , intent(inout) :: khoc !< hoc of linked list
  integer(ki4) , intent(in) :: knew(:) !< new entry


  !! local
  character(*), parameter :: s_name='ld_add' !< subroutine name
  integer(ki4) :: ihoc  !< local variable
  integer(ki4) :: icontsiz  !< local variable
  integer(ki4) :: iconten  !< local variable

  !! test storage
  ihoc=khoc
  iconten=self%nconten
  icontsiz=size(self%conten,2)
  if(khoc<=iconten) then
     if(iconten<icontsiz) then
        iconten=iconten+1
        if (ihoc==0) then
           !! new list
           self%conten(1,iconten)=0
           self%conten(2,iconten)=0
           self%conten(3:,iconten)=knew(:)
        else
           !! add to existing
           self%conten(1,khoc)=iconten
           self%conten(1,iconten)=0
           self%conten(2,iconten)=khoc
           self%conten(3:,iconten)=knew(:)
        end if
        khoc=iconten
        self%nconten=iconten
     else
        call log_error(m_name,s_name,1,error_fatal,'Contents array not large enough')
     end if
  else
     call log_error(m_name,s_name,2,error_fatal,'Hoc bound exceeded')
  end if

end subroutine ld_add
!---------------------------------------------------------------------
!> add into doubly linked list structure (TO DO) synonym for ld_add
subroutine ld_addin(self,khoc,knew)

  !! arguments
  type(ld_t), intent(inout) :: self   !< linked list
  integer(ki4) , intent(inout) :: khoc !< hoc of linked list
  integer(ki4) , dimension(:),intent(in) :: knew !< new entry

  !! local
  character(*), parameter :: s_name='ld_addin' !< subroutine name
  integer(ki4) :: ihoc  !< local variable
  integer(ki4) :: icontsiz  !< local variable
  integer(ki4) :: iconten  !< local variable

  !! test storage
  ihoc=khoc
  iconten=self%nconten
  icontsiz=size(self%conten)
  if(ihoc<=icontsiz.AND.ihoc>=0) then
     if(iconten<=icontsiz) then
        iconten=iconten+1
        if (ihoc==0) then
           !! new list
           self%conten(1,iconten)=0
           self%conten(2,iconten)=0
           self%conten(3:,iconten)=knew(:)
        else
           !! add to existing
           self%conten(1,khoc)=iconten
           self%conten(1,iconten)=0
           self%conten(2,iconten)=khoc
           self%conten(3:,iconten)=knew(:)
        end if
     else
        call log_error(m_name,s_name,1,error_fatal,'Contents array not large enough')
     end if
  else
     call log_error(m_name,s_name,2,error_fatal,'Hoc bound exceeded')
  end if

end subroutine ld_addin
!---------------------------------------------------------------------
!> remove from doubly linked list structure
subroutine ld_remove(self,khoc,key)

  !! arguments
  type(ld_t), intent(inout) :: self   !< linked list
  integer(ki4) , intent(inout) :: khoc !< hoc of linked list
  integer(ki4) , intent(in) :: key(:) !< key of entry

  !! local
  character(*), parameter :: s_name='ld_remove' !< subroutine name
  integer(ki4) :: ihoc  !< local variable
  integer(ki4) :: iconten  !< local variable
  integer(ki4) :: iprev  !< local variable
  integer(ki4) :: inext  !< local variable
  logical :: ilmatch !< local variable

  !! test storage
  ihoc=khoc
  iconten=self%nconten
  if(khoc<=iconten) then
     if (ihoc<=0) then
        !! empty list, nothing to do
     else
        !! find key
        i=ihoc
        do
           ilmatch=.TRUE.
           do k=1,self%ndentry-2
              ilmatch=ilmatch.AND.(self%conten(k+2,i)==key(k))
           end do
           if (ilmatch) exit
           inext=self%conten(2,i)
           if (inext==0) then
              call log_error(m_name,s_name,1,error_warning,'Key not found ')
              return
           end if
           i=inext
        end do
        iprev=self%conten(1,i)
        inext=self%conten(2,i)
        if (iprev>0) then
           self%conten(2,iprev)=inext
        else
           khoc=inext
        end if
        if (inext>0) self%conten(1,inext)=iprev
        self%conten(:,i)=0
     end if
  else
     call log_error(m_name,s_name,2,error_fatal,'Hoc bound exceeded')
  end if

end subroutine ld_remove
!---------------------------------------------------------------------
!> count entries in doubly linked list
subroutine ld_count(self,khoc,kn)

  !! arguments
  type(ld_t), intent(in) :: self   !< linked list
  integer(ki4) , intent(in) :: khoc !< hoc of linked list
  integer(ki4) , intent(out) :: kn !< number in linked list

  !! local
  character(*), parameter :: s_name='ld_count' !< subroutine name
  integer(ki4) :: iconten  !< local variable
  integer(ki4) :: inext  !< local variable

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
           inext=self%conten(2,i)
           if (inext==0) exit
           kn=kn+1
           i=inext
        end do
     end if
  else
     call log_error(m_name,s_name,1,error_fatal,'Hoc bound exceeded')
  end if

end subroutine ld_count
!---------------------------------------------------------------------
!> add one to number of entries in doubly linked list
subroutine ld_countpp(self)

  !! arguments
  type(ld_t), intent(inout) :: self   !< linked list

  !! local
  character(*), parameter :: s_name='ld_countpp' !< subroutine name

  self%nconten=self%nconten+1

end subroutine ld_countpp
!---------------------------------------------------------------------
!> return doubly linked list
subroutine ld_fetch(self,khoc,kout,knum)

  !! arguments
  type(ld_t), intent(in) :: self   !< linked list
  integer(ki4) , intent(in) :: khoc !< hoc of linked list
  integer(ki4), intent(inout) :: kout(:,:)   !< array large enough to contain doubly linked list
  integer(ki4), intent(out) :: knum   !< number of proper entries in kout

  !! local
  character(*), parameter :: s_name='ld_fetch' !< subroutine name
  integer(ki4) :: iconten  !< local variable
  integer(ki4) :: inext  !< local variable
  integer(ki4), dimension(:,:), allocatable  :: iout  !< local variable

  !! test storage
  iconten=self%nconten
  idata=self%ndentry-2
  if(khoc<=iconten) then
     call ld_count(self,khoc,knum)
     if (knum==0) then
        !! empty list
        return
     else
        !! add to array
        i=khoc
        j=knum
        kout(:,j)=self%conten(3:,j)
     end if
     do
        inext=self%conten(2,i)
        if (inext==0) exit
        i=inext
        j=j-1
        kout(:,j)=self%conten(3:,j)
     end do
  else
     call log_error(m_name,s_name,2,error_fatal,'Hoc bound exceeded')
  end if

end subroutine ld_fetch
!---------------------------------------------------------------------
!> write doubly linked list
subroutine ld_write1(self,khoc,kout)

  !! arguments
  type(ld_t), intent(in) :: self   !< linked list
  integer(ki4) , intent(in) :: khoc !< hoc of linked list
  integer, intent(in) :: kout   !< output channel for doubly linked list data

  !! local
  character(*), parameter :: s_name='ld_write1' !< subroutine name
  integer(ki4) :: inum  !< local variable
  integer(ki4) :: iconten  !< local variable
  integer(ki4) :: inext  !< local variable
  integer(ki4), dimension(:,:), allocatable  :: iout  !< local variable

  !! test storage
  iconten=self%nconten
  idata=self%ndentry-2
  if(khoc<=iconten) then
     call ld_count(self,khoc,inum)
     if (inum==0) then
        !! empty list
        return
     else
        !! add to array
        allocate(iout(inum,idata), stat=status)
        call log_alloc_check(m_name,s_name,1,status)
        i=khoc
        j=inum
        iout(:,j)=self%conten(3:,j)
     end if
     do
        inext=self%conten(2,i)
        if (inext==0) exit
        i=inext
        j=j-1
        iout(:,j)=self%conten(3:,j)
     end do
  else
     call log_error(m_name,s_name,2,error_fatal,'Hoc bound exceeded')
  end if

  ! output doubly linked list data array
  write(kout,*,iostat=status) ((iout(k,l),k=1,idata),l=1,inum)
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error writing linked list data')
  end if
  deallocate(iout)

end subroutine ld_write1
!---------------------------------------------------------------------
!> return next entry in doubly linked list
  integer(ki4) function ld_next1(self,khoc)

!! arguments
  type(ld_t), intent(inout) :: self   !< linked list
  integer(ki4) , intent(in) :: khoc !< hoc of linked list

!! local
  character(*), parameter :: s_name='ld_next1' !< subroutine name
  integer(ki4) :: iconten  !< local variable
  integer(ki4) :: inext  !< local variable

!! test storage
  iconten=self%nconten
  if(khoc<=iconten) then
     i=khoc
     inext=0
     if (i<=0) then
!! empty list
     else
        j=self%last(1)
        if (j>0) then
! point to next
           inext=self%conten(2,j)
        else
! else point to head of chain
           inext=i
        end if
     end if
  else
     call log_error(m_name,s_name,1,error_fatal,'Hoc bound exceeded')
  end if

! output index and update
  ld_next1=inext
  self%last(1)=inext

end function ld_next1

end module ld_m
