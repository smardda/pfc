module indict_m

  use const_kind_m
  use const_numphys_h
  use log_m

  implicit none
  private

! public subroutines
  public ::  &
  indict_fromlist, &  !< set up dictionary
  indict_sort, &  !< sort keys into order
  indict_bycluster, &  !< dictionary from clustering of array
  indict, &  !< find entry in dictionary
  indict_write  !< write out dictionary

  interface indict
   module procedure indict_ki4
   module procedure indict_kr4  !kr4
   module procedure indict_kr8
  end interface

  interface indict_write
   module procedure indict_write_ki4
   module procedure indict_write_kr4  !kr4
   module procedure indict_write_kr8
  end interface

!public variables

! private types

! private variables
  character(*), parameter :: m_name='indict_m' !< module name
  integer   :: status   !< error status
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: ji !< loop counter

  contains
!---------------------------------------------------------------------
!> create dictionary from input list
subroutine indict_fromlist(self,list,kno)

  !! arguments
  integer(ki4), dimension(:), allocatable, intent(out) :: self !< dictionary
  integer(ki4), dimension(:), intent(in) :: list !< input list
  integer(ki4), intent(out) :: kno !< dictionary dimension

  !! local
  character(*), parameter :: s_name='indict_fromlist' !< subroutine name
  integer(ki4)  :: inbod !< number of entries in list
  integer(ki4), dimension(:), allocatable :: iwork !< scratch dictionary
  integer(ki4)  :: inew !< new entry in scratch dictionary

  !! allocate array
  inbod=ubound(list,1)
  if(inbod>0) then
     allocate(iwork(inbod),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No input list data')
  end if

  !! find different entries
  kno=1
  iwork(kno)=list(1)
  do j=2,inbod
     inew=1
     do k=1,kno
        if (list(j)==iwork(k)) then
           inew=0
           exit
        end if
     end do
     if (inew==1) then
        kno=kno+1
        iwork(kno)=list(j)
     end if
  end do
  !! copy to smaller output array
  allocate(self(kno),stat=status)
  call log_alloc_check(m_name,s_name,3,status)
  self=iwork(1:kno)

  deallocate(iwork)

end subroutine indict_fromlist
!---------------------------------------------------------------------
!> sort keys into order
subroutine indict_sort(self,selfo,kno,kopt)

  !! arguments
  integer(ki4), dimension(*), intent(inout) :: self !< dictionary
  integer(ki4), dimension(:), allocatable, intent(out) :: selfo !< dictionary
  integer(ki4), intent(out) :: kno !< dictionary dimension
  integer(ki4), intent(in) :: kopt !< isort flag option

  !! local
  character(*), parameter :: s_name='indict_sort' !< subroutine name

  !! allocate piggyback array if required
  if (abs(kopt)==2) then
     !! allocate
     if(kno>0) then
        allocate(selfo(kno),stat=status)
        call log_alloc_check(m_name,s_name,1,status)
     else
        call log_error(m_name,s_name,2,error_fatal,'No dictionary list data')
     end if
  else
     !! for ifort, must allocate even if not used in isort
     allocate(selfo(1),stat=status)
     call log_alloc_check(m_name,s_name,3,status)
  end if

  !! do sort
  call isort(self,selfo,kno,kopt)

end subroutine indict_sort
!---------------------------------------------------------------------
!> dictionary from clustering of array
subroutine indict_bycluster(self,p,kdim,pdel,kdict)

  !! arguments
  real(kr8), dimension(:), allocatable, intent(out) :: self !< dictionary
  real(kr8), dimension(kdim), intent(in) :: p !< array to be clustered
  integer(ki4), intent(in) :: kdim !< array size
  real(kr8), intent(in) :: pdel !< cluster distance
  integer(ki4), intent(out) :: kdict !< dictionary size
  !! local
  character(*), parameter :: s_name='indict_bycluster' !< subroutine name

  logical, dimension(:), allocatable :: marker !< local variable
  real(kr8) :: zmin !< minimum array value
  integer(ki4) :: ibin !< number of bins
  real(kr8), dimension(:,:), allocatable :: wdict !< working dictionary

  allocate(marker(kdim), wdict(2,kdim), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  marker=.TRUE.
  do j=1,kdim
     zmin=minval(p,marker)
     wdict(1,j)=zmin
     wdict(2,j)=zmin
     do k=1,kdim
        if (marker(k)) then
           if (p(k)-zmin<pdel) then
              marker(k)=.FALSE.
              wdict(2,j)=max(p(k),wdict(2,j))
           end if
        end if
     end do
     ibin=j
     if (.NOT.any(marker)) exit
  end do

  kdict=ibin+1
  allocate(self(kdict), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  do j=2,ibin
     self(j)=(wdict(2,j-1)+wdict(1,j))/2
  end do
  self(1)=-const_infty
  self(kdict)=const_infty

  deallocate(marker,wdict)

end subroutine indict_bycluster
!---------------------------------------------------------------------
function indict_ki4(dict,word,ndim)
  integer(ki4) :: indict_ki4 !< dictionary index
  integer(ki4), dimension(*), intent(in) :: dict !< dictionary
  integer(ki4), intent(in) :: word !< dictionary
  integer(ki4), intent(in) :: ndim !< dim of dict
  indict_ki4=0
  !write(*,*) 'ndim=',ndim
  !write(*,*) 'word=',word
  !write(*,*) 'dict=',dict
  do ji=1,ndim
     if (word==dict(ji)) then
        indict_ki4=ji
        return
     end if
  end do
end function indict_ki4

function indict_kr4(dict,word,ndim)
  integer(ki4) :: indict_kr4 !< dictionary index
  real(kr4), dimension(*), intent(in) :: dict !< dictionary
  real(kr4), intent(in) :: word !< dictionary
  integer(ki4), intent(in) :: ndim !< dictionary size
  indict_kr4=0
  !write(*,*) 'ndim=',ndim
  !write(*,*) 'word=',word
  !write(*,*) 'dict=',dict
  indict_kr4=0
  do ji=1,ndim-1
     if (word>dict(ji).AND.word<=dict(ji)) then
        indict_kr4=ji
        return
     end if
  end do
  if (word>dict(ndim)) indict_kr4=ndim
end function indict_kr4

function indict_kr8(dict,word,ndim)
  integer(ki4) :: indict_kr8 !< dictionary index
  real(kr8), dimension(*), intent(in) :: dict !< dictionary
  real(kr8), intent(in) :: word !< dictionary
  integer(ki4), intent(in) :: ndim !< dictionary size
  indict_kr8=0
  !write(*,*) 'ndim=',ndim
  !write(*,*) 'word=',word
  !write(*,*) 'dict=',dict
  indict_kr8=1
  do ji=1,ndim-1
     if (word>dict(ji).AND.word<=dict(ji+1)) then
        indict_kr8=ji
        return
     end if
  end do
  if (word>dict(ndim)) indict_kr8=ndim
end function indict_kr8
!---------------------------------------------------------------------
subroutine indict_write_ki4(self,ndim,nout)
  integer(ki4), dimension(*), intent(in) :: self !< dictionary
  integer(ki4), intent(in) :: ndim !< dictionary size
  integer, intent(in) :: nout !< output unit name

  !! local
  character(*), parameter :: s_name='indict_write_ki4' !< subroutine name

  write(nout,*,iostat=status) 'int_dictionary'
  call log_write_check(m_name,s_name,1,status)

  do ji=1,ndim
     write(nout,*,iostat=status) ji,self(ji)
     call log_write_check(m_name,s_name,2,status)
  end do

end subroutine indict_write_ki4

subroutine indict_write_kr4(self,ndim,nout)  !kr4
  real(kr4), dimension(*), intent(in) :: self !< dictionary  !kr4
  integer(ki4), intent(in) :: ndim !< dictionary size  !kr4
  integer, intent(in) :: nout !< output unit name  !kr4
  !kr4
  !! local  !kr4
  character(*), parameter :: s_name='indict_write_kr4' !< subroutine name  !kr4
  !kr4
  write(nout,*,iostat=status) 'real_dictionary'  !kr4
  call log_write_check(m_name,s_name,1,status)  !kr4
  !kr4
  do ji=1,ndim  !kr4
     write(nout,*,iostat=status) ji,self(ji)  !kr4
     call log_write_check(m_name,s_name,2,status)  !kr4
  end do  !kr4
  !kr4
end subroutine indict_write_kr4  !kr4

subroutine indict_write_kr8(self,ndim,nout)
  real(kr8), dimension(*), intent(in) :: self !< dictionary
  integer(ki4), intent(in) :: ndim !< dictionary size
  integer, intent(in) :: nout !< output unit name

  !! local
  character(*), parameter :: s_name='indict_write_kr8' !< subroutine name

  write(nout,*,iostat=status) 'real_kr8_dictionary'
  call log_write_check(m_name,s_name,1,status)

  do ji=1,ndim
     write(nout,*,iostat=status) ji,self(ji)
     call log_write_check(m_name,s_name,2,status)
  end do

end subroutine indict_write_kr8

end module indict_m
