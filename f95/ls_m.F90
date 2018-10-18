module ls_m

  use const_kind_m
  use log_m
  use control_h

  implicit none
  private

! public subroutines
  public ::  &
  ls_init,   & ! create  list structure
  ls_delete,   & ! delete  list structure
  ls_read,   & ! read  list structure
  ls_write, &  ! write  list structure
  ls_write1, &  ! write partial list
  ls_add,   & ! add to  list structure
  ls_copy    ! copy entries in list structure

! public types
  type, public :: ls_t
   !! list
     integer(ki4), dimension(:,:), allocatable  :: list  !< list array
     integer(ki4) :: nlist   !< number of entries in `2' array
     integer(ki4) :: n2   !< second number of entries in `2' array
     integer(ki4) :: listype    !< Type of list
  end type ls_t

!public variables

! private types

! private variables
  character(*), parameter :: m_name='ls_m' !< module name
  integer   :: status   !< error status
  integer(ki4) :: first !< head of ls
  integer(ki4) :: inext !< next of ls
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: ji !< loop counter

  contains
!---------------------------------------------------------------------
!> initialise ls
subroutine ls_init(self,ksizeh,ksize,ktype)

  !! arguments
  type(ls_t), intent(out) :: self   !< list
  integer(ki4) , intent(in) :: ksize !< size of list
  integer(ki4) , intent(in) :: ksizeh !< second size of list
  integer(ki4) , intent(in), optional :: ktype !< type of list


  !! local
  character(*), parameter :: s_name='ls_init' !< subroutine name
  integer(ki4) :: isize  !< local variable

  isize=ksizeh
  !! check values
  if (present(ktype)) then
     tree_type: select case (ktype)
     case(1)
        isize=0
     end select tree_type
  end if

  !! allocate storage
  if(isize>=0) then
     if(ksize>0) then
        allocate(self%list(ksize,0:isize), stat=status)
        call log_alloc_check(m_name,s_name,1,status)
     else
        call log_error(m_name,s_name,2,error_fatal,'No data')
     end if
  else
     call log_error(m_name,s_name,3,error_fatal,'No data')
  end if

  self%nlist=0
  self%n2=isize

end subroutine ls_init
!---------------------------------------------------------------------
!! delete ls_t
subroutine ls_delete(self)

  !! arguments
  type(ls_t), intent(inout) :: self !< particle coord data
  !! local

  deallocate(self%list)

end subroutine ls_delete
!---------------------------------------------------------------------
!> read in ls data
subroutine ls_read(self,numerics,kin)

  !! arguments
  type(ls_t), intent(out) :: self   !< ls data
  type(numerics_t), intent(inout) :: numerics !< local variable
  integer, intent(in) :: kin   !< input channel for ls data


  !! local
  character(*), parameter :: s_name='ls_read' !< subroutine name
  integer(ki4) :: ilist !< local variable
  integer(ki4) :: in2  !< local variable
  integer(ki4) :: isizel !< local variable
  integer(ki4) :: isizeh  !< local variable

  ! read no in list array
  read(kin,*,iostat=status) ilist,in2
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading list nlist')
  end if
  isizel=ilist
  isizeh=in2
  numerics%nsizeh=isizeh
  numerics%nsizel=isizel
  call ls_init(self,isizeh,isizel)
  ! init sets nlist to zero
  self%nlist=ilist

  ! read ls list array (only last or 2 array signifies)
  read(kin,*,iostat=status) (self%list(j,2),j=1,ilist)
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error reading ls list')
  end if

end subroutine ls_read
!---------------------------------------------------------------------
!> output ls vector
subroutine ls_write(self,numerics,kout)

  !! arguments
  type(ls_t), intent(in) :: self   !< ls data
  type(numerics_t), intent(in) :: numerics !< local variable
  integer, intent(in) :: kout   !< output channel for ls data


  !! local
  character(*), parameter :: s_name='ls_write' !< subroutine name
  integer(ki4) :: ilist !< local variable
  integer(ki4) :: in2  !< local variable


  ! output no in list array
  ilist=self%nlist
  in2=numerics%nsizeh
  write(kout,*,iostat=status) ilist,in2
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error writing list nlist')
  end if

  ! output ls list array

  write(kout,*,iostat=status) (self%list(j,2),j=1,ilist)
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error writing list list')
  end if

end subroutine ls_write
!---------------------------------------------------------------------
subroutine ls_write1(self,kadr,kout)

  !! arguments
  type(ls_t), intent(in) :: self   !< list
  integer(ki4) , intent(in) :: kadr !< start number of list
  integer, intent(in) :: kout   !< output channel for ls data

  !! local
  character(*), parameter :: s_name='ls_write1' !< subroutine name
  integer(ki4) :: in  !< local variable
  integer(ki4), dimension(:), allocatable  :: iout  !< local variable
  integer(ki4) :: ilines !< number of output lines in block

  !! test storage
  if (kadr<=self%nlist) then
     in=self%list(kadr,2)
     if (in==0) then
        !! empty list
        return
     else
        !! add to array
        allocate(iout(in), stat=status)
        call log_alloc_check(m_name,s_name,1,status)
     end if
     loop_1 : do j=1,in
        iout(j)=self%list(kadr+j,2)
     end do loop_1
  else
     call log_error(m_name,s_name,2,error_fatal,'List array bound exceeded')
  end if

  ! output ls list array
  ilines=(in+9)/10
  !!write
  write(kout,'("layout ",2i8)') ilines,in
  write(kout,'(10(1x,i8))',iostat=status) (iout(k),k=1,in)
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error writing ls list')
  end if
  deallocate(iout)

end subroutine ls_write1
!---------------------------------------------------------------------
subroutine ls_add(self,kadr,k2,kitem)

  !! arguments
  type(ls_t), intent(inout) :: self   !< list
  integer(ki4) , intent(in) :: kadr !< first entry of list
  integer(ki4) , intent(in) :: k2 !< second entry of list
  integer(ki4) , intent(in) :: kitem !< item to add to list

  !! local
  character(*), parameter :: s_name='ls_add' !< subroutine name
  integer(ki4) :: isiz  !< local variable
  integer(ki4):: il  !< local variable
  integer(ki4):: iu  !< local variable

  !! test storage
  isiz=size(self%list,1)
  il=lbound(self%list,2)
  iu=ubound(self%list,2)
  if(il<=k2.AND.k2<=iu) then
     if(kadr<=isiz) then
        !! add to existing
        self%list(kadr,k2)=kitem
     else
        call log_error(m_name,s_name,1,error_fatal,'List array not large enough, increase btree_sizel')
     end if
  else
     call log_error(m_name,s_name,2,error_fatal,'List array not large enough, increase btree_sizel')
  end if

end subroutine ls_add
!---------------------------------------------------------------------
subroutine ls_copy(self,kadr1,k1,kadr2,k2,klen)

  !! arguments
  type(ls_t), intent(inout) :: self   !< list
  integer(ki4) , intent(in) :: kadr1 !< first entry of list subarray
  integer(ki4) , intent(in) :: k1 !< second entry of list
  integer(ki4) , intent(in) :: kadr2 !< first entry of list subarray
  integer(ki4) , intent(in) :: k2 !< second entry of list
  integer(ki4) , intent(in) :: klen !< length of list subarray

  !! local
  character(*), parameter :: s_name='ls_copy' !< subroutine name
  integer(ki4) :: isiz  !< local variable
  integer(ki4):: il  !< local variable
  integer(ki4):: iu  !< local variable
  integer(ki4):: ioff  !< local variable

  !! test storage
  ioff=klen-1
  isiz=size(self%list,1)
  il=lbound(self%list,2)
  iu=ubound(self%list,2)
  if(il<=k1.AND.k1<=iu.AND.il<=k2.AND.k2<=iu) then
     if(kadr1+ioff<=isiz.AND.kadr2+ioff<=isiz) then
        !! add to existing
        self%list(kadr2:kadr2+ioff,k2)=self%list(kadr1:kadr1+ioff,k1)
     else
        call log_error(m_name,s_name,1,error_fatal,'List array not large enough')
     end if
  else
     call log_error(m_name,s_name,2,error_fatal,'List array not large enough')
  end if

end subroutine ls_copy

end module ls_m
