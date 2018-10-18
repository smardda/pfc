module btree_m

  use const_kind_m
  use log_m
  use control_h
  use position_h
  use position_m
  use geobj_m
  use query_m
  use ls_m

  implicit none
  private

! public subroutines
  public ::  &
  btree_init,   & !< create  tree data structure
  btree_write, &  !< write  binary tree structure
  btree_read,   & !< read  binary tree structure
  btree_writev, & !< write (vtk) binary tree structure
  btree_add,    & !< split BSP node and assign to new nodes presorted array of start addresses
  btree_limits,    & !< limits of normalised btree
  btree_mdefn,   & !< define boxes for node split
  btree_madd,   & !< split node and assign to new nodes presorted array of start addresses
  btree_find,   & !< find node containing normalised point vector in BSP btree
  btree_mfind,   & !< locate node containing point vector in many-one btree
  btree_dia, &  !< object diagnostics to log file
  btree_newnode,   & !< return node containing point vector starting from given node
  btree_readnolab, & !< read unlabelled binary tree structure
  btree_floatvec, &  !< return float btree short integer vector
  btree_delete !< delete object
!<       btree_match,   & ! match in binary tree
!<       btree_next1, &  ! next at numbered level in binary tree
!<       btree_next, &  ! next in binary tree (bottom level)


! public types

!> type storing tree coord data
  type, public :: btree_t
     integer(ki4)  :: nt !< number of entries in tree
     integer(ki4), dimension(:,:), allocatable :: pter !< tree links
     integer(ki2), dimension(:,:), allocatable :: desc !< description
     integer(ki2), dimension(:,:), allocatable :: corner !< scaled position
     integer(ki2), dimension(:,:), allocatable :: exten !< array of exten descriptions
     integer(ki2), dimension(3) :: nxyz !< array of children
     integer(ki4) :: npter !< number of entries in pter
     integer(ki4) :: nexten !< number of entries in exten
     integer(ki4) :: nroot !< number of root record
     integer(ki2):: ndepth !< max depth of tree
     integer(ki4):: nttype !< type of tree
     integer(ki4):: mtype !< type of margin for quantising
     integer(ki4) :: nttalg !< type of binary tree algorithm for top node
     real(kr4), dimension(3)  :: hxyz !< specified mesh sizes
     type(ls_t) :: objectls !< general list object for children
     integer(ki4) :: maxchildn !< maximum number of children of a node
  end type btree_t

!public variables


! private types

! private variables
  character(*), parameter :: m_name='btree_m' !< module name
  integer   :: status   !< error status
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ii !< loop counter
  integer(ki4) :: jj !< loop counter
  integer(ki4) :: kk !< loop counter
  integer(ki2):: idepth !< max depth of tree
  integer :: ilog !< file unit for logging

  contains

!> create  tree data structure
subroutine btree_init(self,numerics)

  !! arguments
  type(btree_t), intent(out) :: self   !< binary tree
  type(numerics_t), intent(in) :: numerics   !< binary tree

  !! local
  character(*), parameter :: s_name='btree_init' !< subroutine name
  integer(ki4) :: isize !< size of binary tree
  integer(ki4) :: isizep !< size of binary tree pter
  integer(ki4) :: isizee !< size of exten array
  integer(ki4) :: isizeh !< second size of list array
  integer(ki4) :: isizel !< size of list array
  integer(ki4):: ittype !< type of tree
  integer(ki4):: ittalg !< type of top of tree algorithm
  integer(ki4) :: ibsize  !< btree size
  integer(ki2):: idepthe !< max depth of tree (provisional)
  integer(ki2):: iquante !< max scaled extent of tree (provisional, log)
  integer(ki4) :: igeobj  !< local variable
  integer(ki4) :: imtype  !< local variable
  integer(ki2), dimension(3) :: ixyz  !< local variable
  real(kr4), dimension(3) :: hxyz  !< local variable

  isize=numerics%nsize
  isizep=numerics%nsizep
  isizee=numerics%nsizee
  isizeh=numerics%nsizeh
  isizel=numerics%nsizel
  idepthe=numerics%ndepth
  ittype=numerics%nttype
  ixyz=numerics%nxyz
  ittalg=numerics%nttalg
  hxyz=numerics%hxyz
  imtype=numerics%mtype

  !! allocate storage
  if(isize>0.AND.isizep>0) then
     allocate(self%pter(isizep,isize), &
 &   stat=status)
     !! check successful allocation
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Unable to allocate memory')
     end if
  else
     call log_error(m_name,s_name,2,error_fatal,'No data')
  end if

  if(isize>0) then
     allocate(self%desc(2,isize),self%corner(3,isize), &
 &   stat=status)
     !! check successful allocation
     if(status/=0) then
        call log_error(m_name,s_name,3,error_fatal,'Unable to allocate memory')
     end if
  else
     call log_error(m_name,s_name,4,error_fatal,'No data')
  end if

  ibsize=isizee
  if(ibsize>0) then
     allocate(self%exten(3,ibsize), &
 &   stat=status)
     !! check successful allocation
     if(status/=0) then
        call log_error(m_name,s_name,5,error_fatal,'Unable to allocate memory')
     end if
  else
     call log_error(m_name,s_name,6,error_fatal,'No data')
  end if

  self%nt=1
  self%npter=isizep
  self%pter(1,1)=-1
  self%pter(2,1)=2
  self%pter(3,1)=0
  self%desc(1,1)=0
  self%desc(2,1)=1
  self%corner(:,1)=0

  self%ndepth=idepthe
  self%nttype=ittype
  self%nxyz=ixyz
  self%nttalg=ittalg
  self%hxyz=hxyz
  self%mtype=imtype

  if (ittype==2) then
     iquante=numerics%nquante-1
  else
     iquante=numerics%nquante
  end if
  self%nexten=1
  self%exten(:,1)=(/iquante,iquante,iquante/)

  call ls_init(self%objectls,isizeh,isizel)
  ! define dummy first record (0 objects) in list 0 and 2
  call ls_add(self%objectls,1,0,0)
  call ls_add(self%objectls,1,2,0)
  !INERT  call ls_add(self%objectls,2,0,0)
  !INERT  call ls_add(self%objectls,1,1,0)
  !INERT  call ls_add(self%objectls,2,1,1)
  ! define list of objects
  igeobj=numerics%ngeobj
  call ls_add(self%objectls,2,2,igeobj)
  do j=1,igeobj
     call ls_add(self%objectls,j+2,2,j)
  end do

end subroutine btree_init
!---------------------------------------------------------------------
!> write out binary tree data
subroutine btree_write(self,numerics,kout)

  !! arguments
  type(btree_t), intent(in) :: self   !< binary tree data
  type(numerics_t), intent(in) :: numerics    !< local variable
  integer, intent(in) :: kout   !< output channel for binary tree data

  !! local
  character(*), parameter :: s_name='btree_write' !< subroutine name
  integer(ki4) :: intr !< local variable
  integer(ki4) :: inpter  !< local variable
  integer(ki4) :: ittype  !< local variable
  integer(ki4) :: iexten  !< local variable
  integer(ki2), dimension(3) :: ixyz  !< local variable

  ! write btree label

  write(kout,'(a)',iostat=status) 'BTREE'
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error writing binary tree label')
  end if

  ! write no in binary tree array
  intr=self%nt
  inpter=self%npter
  idepth=self%ndepth
  ittype=self%nttype
  write(kout,*,iostat=status) intr,inpter,idepth,ittype
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error writing binary tree sizes')
  end if

  if (ittype/=1) then
     ixyz=self%nxyz
     write(kout,*,iostat=status) ixyz
     if(status/=0) then
        call log_error(m_name,s_name,9,error_fatal,'Error writing top node children')
     end if
  end if

  ! write binary tree arrays

  write(kout,*,iostat=status) &
  (self%pter(1:inpter,j),self%desc(:,j),self%corner(:,j),j=1,intr)
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error writing binary tree entries')
  end if

  ! write exten label
  write(kout,'(a)',iostat=status) 'EXTEN'
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error writing binary tree exten label')
  end if

  ! write no in exten array
  iexten=self%nexten
  write(kout,*,iostat=status) iexten
  if(status/=0) then
     call log_error(m_name,s_name,5,error_fatal,'Error writing binary tree nexten')
  end if

  ! write binary tree exten array
  write(kout,*,iostat=status) (self%exten(:,j),j=1,iexten)
  if(status/=0) then
     call log_error(m_name,s_name,6,error_fatal,'Error writing binary tree exten')
  end if

  ! write btree ls label
  write(kout,'(a)',iostat=status) 'LSLAB'
  if(status/=0) then
     call log_error(m_name,s_name,7,error_fatal,'Error writing binary tree ls label')
  end if

  ! write binary tree ls array
  call ls_write(self%objectls,numerics,kout)

end subroutine btree_write
!---------------------------------------------------------------------
!> read binary tree data
subroutine btree_read(self,numerics,kchar,kin)

  !! arguments
  type(btree_t), intent(inout) :: self   !< binary tree data
  type(numerics_t), intent(out) :: numerics    !< local variable
  character(len=5) :: kchar !< local variable
  integer, intent(in) :: kin   !< output channel for binary tree data


  !! local
  character(*), parameter :: s_name='btree_read' !< subroutine name
  integer(ki4) :: intr !< local variable
  integer(ki4) :: inpter  !< local variable
  integer(ki4) :: ittype  !< local variable
  integer(ki4) :: isize !< local variable
  integer(ki4) :: isizep  !< local variable
  integer(ki4) :: iexten  !< local variable
  integer(ki4) :: isizee  !< local variable
  integer(ki2), dimension(3) :: ixyz  !< local variable

  numerics%nbdcub=0
  if (kchar=='BTREE') then
     ! read no in binary tree array
     read(kin,*,iostat=status) intr,inpter,idepth,ittype
     if(status/=0) then
        call log_error(m_name,s_name,8,error_fatal,'Error reading binary tree sizes')
     end if
     self%nt=intr
     self%npter=inpter
     self%ndepth=idepth
     self%nttype=ittype

     if (ittype/=1) then
        read(kin,*,iostat=status) ixyz
        if(status/=0) then
           call log_error(m_name,s_name,9,error_fatal,'Error reading top node children')
        end if
        self%nxyz=ixyz
     end if

     ! allocate binary tree arrays
     isize=intr
     isizep=inpter
     !! allocate storage
     if(isize>0.AND.isizep>0) then
        allocate(self%pter(isizep,isize), &
 &      stat=status)
        !! check successful allocation
        if(status/=0) then
           call log_error(m_name,s_name,1,error_fatal,'Unable to allocate memory')
        end if
     else
        call log_error(m_name,s_name,2,error_fatal,'No data')
     end if

     if(isize>0) then
        allocate(self%desc(2,isize),self%corner(3,isize), &
 &      stat=status)
        !! check successful allocation
        if(status/=0) then
           call log_error(m_name,s_name,3,error_fatal,'Unable to allocate memory')
        end if
     else
        call log_error(m_name,s_name,4,error_fatal,'No data')
     end if
     ! read binary tree arrays
     read(kin,*,iostat=status) &
     (self%pter(1:inpter,j),self%desc(:,j),self%corner(:,j),j=1,intr)
     if(status/=0) then
        call log_error(m_name,s_name,3,error_fatal,'Error reading binary tree entries')
     end if
     numerics%nsize=isize
     numerics%nsizep=isizep

  else if (kchar=='EXTEN') then
     ! read no in exten array
     read(kin,*,iostat=status) iexten
     if(status/=0) then
        call log_error(m_name,s_name,5,error_fatal,'Error reading binary tree nexten')
     end if
     self%nexten=iexten

     ! allocate binary tree exten array
     isizee=iexten
     if(isizee>0) then
        allocate(self%exten(3,isizee), &
 &      stat=status)
        !! check successful allocation
        if(status/=0) then
           call log_error(m_name,s_name,5,error_fatal,'Unable to allocate memory')
        end if
     else
        call log_error(m_name,s_name,6,error_fatal,'No data')
     end if
     numerics%nsizee=isizee
     ! read binary tree exten array
     read(kin,*,iostat=status) (self%exten(:,j),j=1,iexten)
     if(status/=0) then
        call log_error(m_name,s_name,6,error_fatal,'Error reading binary tree exten')
     end if

  else if (kchar=='LSLAB') then
     ! read binary tree ls array
     call ls_read(self%objectls,numerics,kin)

  end if

end subroutine btree_read
!---------------------------------------------------------------------
!> write in visualisation format
subroutine btree_writev(self,numerics,kchar,kout)

  !! arguments
  type(btree_t), intent(in) :: self   !< binary tree data
  type(numerics_t), intent(inout) :: numerics !< local variable
  character(*),intent(in):: kchar !< control output
  integer, intent(in) :: kout   !< output channel for binary tree data

  !! local
  character(*), parameter :: s_name='btree_writev' !< subroutine name
  integer(ki4) :: intr  !< local variable
  integer(ki2), dimension(3) :: iorig  !< local variable
  integer(ki2), dimension(3) :: icell  !< local variable
  integer(ki2), dimension(3) :: ipos  !< local variable
  integer(ki2) :: iz  !< local variable
  integer(ki4) :: ioff  !< local variable
  integer(ki4) :: inode !< number of boxes in total
  integer(ki4) :: iptr2  !< local variable
  type(posvecl_t) :: zorig !< local variable
  type(posvecl_t) :: zorigt !< local variable
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zpost !< local variable

  !      type(posvecl_t) :: position_invqtfm,btree_floatvec
  !      external position_invqtfm,btree_floatvec
  iz=0
  ! call vfile_init must have been made
  !! write coordinates
  intr=self%nt
  write(kout,'(''DATASET UNSTRUCTURED_GRID'')')
  plot_type: select case (kchar)
  case('hds')
     inode=intr
     write(kout,'(''POINTS  '',i8, '' float'')') 8*inode

     do j=1,inode

        iorig=self%corner(:,j)
        icell=2**self%exten(:,self%desc(2,j))
        zorig=btree_floatvec(iorig)
        zorigt=position_invqtfm(zorig,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zorigt%posvec
        ! write corner positions
        ipos=iorig+(/icell(1),iz,iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,icell(2),iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/icell(1),icell(2),iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,iz,icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/icell(1),iz,icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,icell(2),icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+icell
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec

     end do

  case('hds lowest')
     ! only print lowest nodes in tree
     inode=0
     do j=1,intr
        iptr2=self%pter(3,j)
        if (iptr2>0) cycle
        inode=inode+1
     end do

     write(kout,'(''POINTS  '',i8, '' float'')') 8*inode

     do j=1,intr

        iptr2=self%pter(3,j)
        if (iptr2>0) cycle
        iorig=self%corner(:,j)
        icell=2**self%exten(:,self%desc(2,j))
        zorig=btree_floatvec(iorig)
        zorigt=position_invqtfm(zorig,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zorigt%posvec
        ! write corner positions
        ipos=iorig+(/icell(1),iz,iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,icell(2),iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/icell(1),icell(2),iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,iz,icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/icell(1),iz,icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,icell(2),icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+icell
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,numerics%geobj_coord_tfm)
        write(kout,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec

     end do

  case('hds quantised')
     inode=intr
     write(kout,'(''POINTS  '',i8, '' unsigned_int'')') 8*inode

     do j=1,intr

        iorig=self%corner(:,j)
        icell=2**self%exten(:,self%desc(2,j))
        ! write corner positions

        write(kout,'(i8,1x,i8,1x,i8)') iorig,iorig+(/icell(1),iz,iz/), &
 &      iorig+(/iz,icell(2),iz/),iorig+(/icell(1),icell(2),iz/), &
 &      iorig+(/iz,iz,icell(3)/),iorig+(/icell(1),iz,icell(3)/), &
 &      iorig+(/iz,icell(2),icell(3)/),iorig+icell

     end do

  case default
     call log_error(m_name,s_name,1,error_warning,'Plot type unrecognised')
  end select plot_type

  write(kout,'(''CELLS  '',i8,1x,i8)') inode, 9*inode

  ioff=0
  do j=1,inode

     ! write blocks
     write(kout,'("8 ",8i10)') (k,k=ioff,ioff+7)

     ioff=ioff+8

  end do

  !!write block types
  write(kout,'(''CELL_TYPES '', i8)') inode

  do j=1,inode
     write(kout,'("11")')
  end do

  !AB   !!format statements for vtk file
  !AB   write(kout,'(''CELL_DATA '', i8)') inode!number of leaves (from 'case(hds lowest)')
  !AB   write(kout,'(''SCALARS '', '' Number '',''int'')')
  !AB   write(kout,'(''LOOKUP_TABLE'','' default'')')
  !AB
  !AB   do j=1,self%nt !!loop over all tree entries
  !AB
  !AB     !!if entry is a leaf then print number of points in leaf
  !AB     !otherwise 'cycle'
  !AB     if (self%pter(3,j)>0) cycle
  !AB     write(kout,'(i4)') self%objectls%list(self%pter(2,j),2)
  !AB
  !AB   end do
  !AB

  close(kout)

end subroutine btree_writev
!---------------------------------------------------------------------
!> split BSP node and assign to new nodes presorted array of start addresses
subroutine btree_add(self,kd,knode,kadra)

  !! arguments
  type(btree_t), intent(inout) :: self   !< binary tree data
  integer(ki2), intent(in) :: kd   !< direction of split
  integer(ki4), intent(in) :: knode   !< node to be split
  integer(ki4), dimension(*), intent(in) :: kadra !< array of new list addresses

  !! local
  character(*), parameter :: s_name='btree_add' !< subroutine name
  integer(ki4) :: intr !< local variable
  integer(ki4) :: inpter  !< local variable
  integer(ki4) :: iexten  !< local variable
  integer(ki4) :: imatch  !< local variable
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: icorn  !< local variable
  integer(ki2), dimension(3) :: iextn  !< local variable
  integer(ki2), dimension(3) :: icornn !< new corner
  integer(ki2) :: iexth  !< local variable

  ! sizes of binary tree array
  intr=self%nt
  inpter=self%npter
  iexten=self%nexten
  ! check space
  if (intr+2>size(self%pter,2)) then
     call btree_dia(self)
     call log_error(m_name,s_name,1,error_fatal,'Binary tree size exceeded')
  end if
  ! complete description of original node
  ! record splitting direction
  self%desc(1,knode)=kd
  self%pter(2,knode)=intr+1
  self%pter(3,knode)=intr+2

  ! new extent array
  iext=self%exten(:,self%desc(2,knode))
  iextn=iext
  iexth=iext(kd)-1
  iextn(kd)=iexth

  ! test against existing entries
  imatch=0
  do j=iexten,1,-1
     if ( all(iextn.eq.self%exten(:,j)) ) then
        imatch=j
        exit
     end if
  end do

  if ( imatch==0 ) then
     ! add to exten array
     if (iexten+1>size(self%exten,2)) then
        call btree_dia(self)
        call log_error(m_name,s_name,2,error_fatal,'Binary tree exten size exceeded')
     end if
     iexten=iexten+1
     self%exten(:,iexten)=iextn
     imatch=iexten
  end if

  ! update corner
  icorn=self%corner(:,knode)
  icornn=icorn
  icornn(kd)=icornn(kd)+2**iexth

  ! add binary tree arrays

  ! first node
  intr=intr+1
  self%desc(1,intr)=0
  self%desc(2,intr)=imatch
  self%corner(:,intr)=icorn
  self%pter(1,intr)=knode
  if (kadra(1)>0) then
     self%pter(2,intr)=kadra(1)
     self%pter(3,intr)=0
  else
     self%pter(2,intr)=1
     self%pter(3,intr)=kadra(1)
  end if

  ! second node
  intr=intr+1
  self%desc(1,intr)=0
  self%desc(2,intr)=imatch
  self%corner(:,intr)=icornn
  self%pter(1,intr)=knode
  if (kadra(2)>0) then
     self%pter(2,intr)=kadra(2)
     self%pter(3,intr)=0
  else
     self%pter(2,intr)=1
     self%pter(3,intr)=kadra(2)
  end if

  ! sizes of binary tree array
  self%nt= intr
  self%nexten= iexten

end subroutine btree_add
!---------------------------------------------------------------------
!> limits of normalised btree
subroutine btree_limits(self,pbb)

  !! arguments
  type(btree_t), intent(in) :: self   !< octree data
  real(kr8), dimension(3,2), intent(out) :: pbb !< limits in normalised space

  !! local
  character(*), parameter :: s_name='btree_limits' !< subroutine name
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2) :: ittype !< local variable
  integer(ki2) :: ichil !< local variable
  integer(ki4), dimension(3) :: ibsiz !< actual size of largest cell
  integer(ki4), dimension(3,2) :: ibb !< actual size of tree

  !! number of children
  ! type of tree
  ittype=self%nttype
  if (ittype==1) then
     ! BSP
     ichil=1
  else if (all(self%nxyz==2)) then
     ! octree
     ichil=2
  else
     ! general nx x ny x nz tree
     ichil=3
  end if

  ! extent array
  iext=self%exten(:,1)
  ibsiz=2**iext

  ! limits
  ibb(:,1)=0
  ibb(:,2)=self%nxyz*ibsiz
  do jj=1,2
     do ii=1,3
        pbb(ii,jj)=real(ibb(ii,jj),kr8)
     end do
  end do

end subroutine btree_limits
!---------------------------------------------------------------------
!> define boxes for node split
subroutine btree_mdefn(self,kxyz,knode,kd,kextn,kcorna,kchildn,kbsiz)

  !! arguments
  type(btree_t), intent(in) :: self   !< octree data
  integer(ki2), dimension(3), intent(in) :: kxyz   !< defines split
  integer(ki4), intent(in) :: knode    !< local variable
  integer(ki2), dimension(3,*), intent(out) :: kcorna !< new corners
  integer(ki2), intent(out) :: kd   !< split direction/code
  integer(ki2), dimension(3), intent(out) :: kextn  !< local variable
  integer(ki4), dimension(3), intent(out) :: kbsiz !< actual size of box
  integer(ki2), intent(out) :: kchildn !< number of children

  !! local
  character(*), parameter :: s_name='btree_mdefn' !< subroutine name
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: icorn  !< local variable
  integer(ki2) :: iexth  !< local variable
  integer(ki2) :: ittype !< local variable
  integer(ki2) :: ichil !< local variable

  !! number of children
  ! type of tree
  ittype=self%nttype
  kchildn=kxyz(1)*kxyz(2)*kxyz(3)
  if (ittype==1) then
     ! BSP
     do jj=1,3
        if (kxyz(jj)==2) then
           kd=jj
           exit
        end if
     end do
     ichil=1
  else if (knode==1) then
     ! root of tree, size of octree
     ichil=3
     kd=3+kchildn
  else
     ! octree
     ichil=2
     kd=3+kchildn
  end if

  ! new extent array
  iext=self%exten(:,self%desc(2,knode))
  if (ichil==1) then
     kextn=iext
     iexth=iext(kd)-1
     kextn(kd)=iexth
     kbsiz=2**kextn
  else if (ichil==2) then
     kextn=iext-(/1,1,1/)
     kbsiz=2**kextn
  else
     ! header node has extent of one octree
     kextn=iext
     kbsiz=2**iext
  end if

  ! update corners
  icorn=self%corner(:,knode)
  do kk=1,kxyz(3)
     do jj=1,kxyz(2)
        do ii=1,kxyz(1)
           i=ii+kxyz(1)*(jj-1+(kk-1)*kxyz(2))
           kcorna(1,i)=icorn(1)+(ii-1)*kbsiz(1)
           kcorna(2,i)=icorn(2)+(jj-1)*kbsiz(2)
           kcorna(3,i)=icorn(3)+(kk-1)*kbsiz(3)
        end do
     end do
  end do

end subroutine btree_mdefn
!---------------------------------------------------------------------
!> split node and assign to new nodes presorted array of start addresses
subroutine btree_madd(self,kxyz,knode,kadra,kd,kextn,kcorna,kchildn,kbsiz)

  !! arguments
  type(btree_t), intent(inout) :: self   !< btree data
  integer(ki2), dimension(*), intent(in) :: kxyz   !< defines split
  integer(ki4), intent(in) :: knode   !< node to be split
  integer(ki4), dimension(*), intent(in) :: kadra !< array of new list addresses
  integer(ki2), intent(in) :: kd   !< split direction
  integer(ki2), dimension(3), intent(in) :: kextn  !< local variable
  integer(ki2), dimension(3,*), intent(in) :: kcorna !< new corners
  integer(ki2), intent(in) :: kchildn !< local variable
  integer(ki4), dimension(3), intent(in) :: kbsiz !< actual size of box


  !! local
  character(*), parameter :: s_name='btree_madd' !< subroutine name
  integer(ki4) :: intr !< local variable
  integer(ki4) :: inpter  !< local variable
  integer(ki4) :: iexten  !< local variable
  integer(ki4) :: imatch  !< local variable
  integer(ki2) :: ittype !< local variable


  ! entries of binary array
  ittype=self%nttype
  intr=self%nt
  inpter=self%npter
  iexten=self%nexten
  ! complete description of original node
  ! record splitting direction
  self%desc(1,knode)=kd
  self%pter(2,knode)=intr+1
  self%pter(3,knode)=intr+2

  ! check space
  if (intr+kchildn>size(self%pter,2)) then
     call btree_dia(self)
     call log_error(m_name,s_name,1,error_fatal,'Binary tree size exceeded')
  end if
  ! test against existing entries
  imatch=0
  do j=iexten,1,-1
     if ( all(kextn.eq.self%exten(:,j)) ) then
        imatch=j
        exit
     end if
  end do

  if ( imatch==0 ) then
     ! add to exten array
     if (iexten+1>size(self%exten,2)) then
        call btree_dia(self)
        call log_error(m_name,s_name,2,error_fatal,'Binary tree exten size exceeded')
     end if
     iexten=iexten+1
     self%exten(:,iexten)=kextn
     imatch=iexten
  end if

  ! add btree arrays

  ! first node
  intr=intr+1
  self%desc(1,intr)=0
  self%desc(2,intr)=imatch
  self%corner(:,intr)=kcorna(:,1)
  self%pter(1,intr)=knode
  if (kadra(1)>0) then
     self%pter(2,intr)=kadra(1)
     self%pter(3,intr)=0
  else
     self%pter(2,intr)=1
     self%pter(3,intr)=kadra(1)
  end if

  do jj=2,kchildn
     ! other nodes
     intr=intr+1
     self%desc(1,intr)=0
     self%desc(2,intr)=imatch
     self%corner(:,intr)=kcorna(:,jj)
     self%pter(1,intr)=knode
     if (kadra(jj)>0) then
        self%pter(2,intr)=kadra(jj)
        self%pter(3,intr)=0
     else
        self%pter(2,intr)=1
        self%pter(3,intr)=kadra(jj)
     end if
  end do

  ! sizes of btree array
  self%nt= intr
  self%nexten= iexten

end subroutine btree_madd
!---------------------------------------------------------------------
!> find node containing normalised point vector in BSP btree
subroutine btree_find(self,obj,query,knode)

  !! arguments
  type(btree_t), intent(in) :: self   !< binary tree data
  type(geobj_t), intent(in) :: obj !< geo object
  type(queryset_t), intent(in) :: query   !< query position data
  integer(ki4), intent(out) :: knode   !< node to be found

  !! local
  character(*), parameter :: s_name='btree_find' !< subroutine name
  integer(ki2) :: id  !< local variable
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: icorn  !< local variable
  integer(ki2), dimension(3) :: iextn  !< local variable

  !      logical geobj_innbox
  !      external geobj_innbox

  !! loop over depth
  k=1
  knode=0
  do i=1,self%ndepth+1
     if (self%pter(3,k) > 0 ) then
        ! node is not terminal, select child
        ! split direction
        id=self%desc(1,k)
        ! corner
        icorn=self%corner(:,k)
        ! extent
        iext=self%exten(:,self%desc(2,k))
        iextn=iext
        iextn(id)=iext(id)-1

        if ( geobj_innbox( obj,query%posl,icorn,iextn )  ) then
           ! in 0 position
           k=self%pter(2,k)
           cycle
        else
           ! in 1 position
           k=self%pter(3,k)
           cycle
        end if
     else if (self%pter(3,k) == -2 ) then
        ! empty node
        knode=k
        exit
     else if (self%pter(3,k) == -1 ) then
        ! leaf node
        knode=k
        exit
     else if (self%pter(3,k) == 0 ) then
        ! incomplete node
        knode=k
        exit
     end if
  end do

  ! check object really is in node
  if (knode<=0) then
     call btree_dia(self)
     call log_error(m_name,s_name,1,error_fatal,'Binary tree does not exist')
  else
     ! corner
     icorn=self%corner(:,knode)
     ! extent
     iext=self%exten(:,self%desc(2,knode))
     if ( .not.geobj_innbox( obj,query%posl,icorn,iext )  ) then
        call log_error(m_name,s_name,2,error_warning,'Object not found in binary tree')
        knode=-1
     end if
  end if

end subroutine btree_find
!---------------------------------------------------------------------
!> locate node containing point vector in many-one btree
subroutine btree_mfind(self,obj,poslis,knode)

  !! arguments
  type(btree_t), intent(in) :: self   !< binary tree data
  type(geobj_t), intent(in) :: obj !< geo object
  type(posveclis_t), intent(in) :: poslis   !< list of position data
  integer(ki4), intent(out) :: knode   !< node to be found

  !! local
  character(*), parameter :: s_name='btree_mfind' !< subroutine name
  type(posvecl_t) :: zpos !< local variable
  integer(ki4) :: inobj   !< geobj position
  integer(ki4) :: jj   !< loop counter
  integer(ki2), dimension(3) :: ivec !< integer geobj data
  integer(ki2), dimension(3) :: inxyz !< top level position

  integer(ki2) :: ichil  !< local variable
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: icorn  !< local variable
  integer(ki2), dimension(3) :: ixyz  !< local variable
  integer(ki2), dimension(3) :: iofftype  !< local variable
  integer(ki2):: iquant !< max scaled extent of tree (log)

  ! check for point
  if (obj%objtyp/=1) then
     call log_error(m_name,s_name,1,error_fatal,'Called with wrong object type')
  end if
  if (self%nttype==1) then
     ! btree_find more appropriate for BSP
     call log_error(m_name,s_name,2,error_fatal,'Called with wrong tree type')
  else if (self%nttype==2) then
     ! the magic factor or one for octree
     iofftype=0
  else
     iofftype=0
  end if

  ! top level analysis
  idepth=self%ndepth
  ixyz=self%nxyz
  inobj=obj%geobj
  zpos=poslis%pos(inobj)
  ivec=int(zpos%posvec)
  iquant=self%exten(2,self%desc(2,1))+iofftype(2)
  do jj=1,3
     inxyz(jj)=ishft(ivec(jj),-iquant)
     if (inxyz(jj).ge.ixyz(jj)) then
        call log_error(m_name,s_name,3,error_warning,'Object not found in binary tree')
        knode=-1
        return
     end if
  end do
  ichil=inxyz(1)+ixyz(1)*(inxyz(2)+inxyz(3)*ixyz(2))

  !D     write(*,*) zpos,inobj,ivec,inxyz
  !! top level node
  k=1
  knode=0
  if (self%pter(3,k) > 0 ) then
     ! node is not terminal, select appropriate child
     k=self%pter(2,k)+ichil
  end if

  if (k==1) then
     call log_error(m_name,s_name,4,error_fatal,'Binary tree has no entries')
  else
     !! loop over depth
     do i=1,idepth-1
        ichil=0
        do jj=1,3
           call mvbits(ivec(jj),iquant-i,1,ichil,jj-1)
        end do
        !D        write(*,*) i,k,ichil
        if (self%pter(3,k) > 0 ) then
           ! node is not terminal, select appropriate child
           k=self%pter(2,k)+ichil
           cycle
           ! split direction

        else if (self%pter(3,k) == -2 ) then
           ! empty node
           knode=k
           exit
        else if (self%pter(3,k) == -1 ) then
           ! leaf node
           knode=k
           exit
        else if (self%pter(3,k) == 0 ) then
           ! incomplete node
           knode=k
           exit
        end if
     end do
  end if

  ! check object really is in node
  if (knode<=0) then
     call btree_dia(self)
     call log_error(m_name,s_name,5,error_fatal,'Binary tree may be corrupt')
  else
     ! corner
     icorn=self%corner(:,knode)
     ! extent
     iext=self%exten(:,self%desc(2,knode))+iofftype
     if ( .not.geobj_innbox( obj,poslis,icorn,iext )  ) then
        call log_error(m_name,s_name,6,error_warning,'Object not found in binary tree')
        knode=-1
     end if
  end if

end subroutine btree_mfind
!---------------------------------------------------------------------
!> output to log file
subroutine btree_dia(self)

  !! arguments
  type(btree_t), intent(in) :: self !< module object
  !! local
  character(*), parameter :: s_name='btree_dia' !< subroutine name
  integer(ki2) :: iexmin !< minimum exten value
  integer(ki2) :: idescent !< binary depth (or descent) actually used
  integer(ki4) :: inleaf  !<  count number of leaves
  integer(ki4) :: inentry  !<  count number of entries
  integer(ki4) :: itlist  !<  local list type
  integer(ki4) :: jl  !< local loop

  call log_getunit(ilog)
  write(ilog,*) ' tree type ', self%nttype
  write(ilog,*) ' split at top of tree ', self%nxyz
  write(ilog,*) ' type of tree algorithm ', self%nttalg
  write(ilog,*) ' size of tree pointer array or numerics%nsize ', &
 &self%nt,' used out of  ',size(self%pter,2), ' (btree_size)'
  write(ilog,*) ' size of extents array or numerics%nsizee ', &
 &self%nexten,' used out of ',size(self%exten,2), ' (btree_sizee)'
  iexmin=minval( self%exten(1:3,1) )
  !dbg  write(*,*) ' 1 ', self%exten(1:3,1)  !dbg
  do jl=2,self%nexten
     iexmin=min( iexmin, minval( self%exten(1:3,jl) ) )
     !dbg     write(*,*) jl, self%exten(1:3,jl)  !dbg
  end do
  idescent=self%ndepth-iexmin
  write(ilog,*) ' descent of tree ', &
 &idescent, ' used out of ', self%ndepth, ' (log_2(quantising number))'
  ! count leaves and entries
  inleaf=0
  !! loop over all tree entries
  do jl=1,self%nt
     if (self%pter(3,jl)>0) cycle
     inleaf=inleaf+1
  end do
  write(ilog,*) ' number of leaves ', inleaf
  write(ilog,*) ' maximum number of entries in any node given by limit_geobj_in_bin'
  itlist=0 ! hardwired
  write(ilog,*) ' type of list structure ', itlist
  list_type: select case (itlist)
  case(0)
     write(ilog,*) 'list type  ls'
     write(ilog,*) ' number of entries in list array numerics%nsizel ', &
 &   self%objectls%nlist,' used out of ',size(self%objectls%list,1), ' (btree_sizel)'
  case(1)
     write(ilog,*) 'list type  li'
     !write(ilog,*) ' number of entries in hoc array ', &
     !&self%objectli%nhoc,' used out of ',self%n%nsizeh
     !     write(ilog,*) ' number of entries in conten array ', &
     ! &   self%objectli%nconten,' used out of ',self%n%nsizel,  ' (max_size_list_array)'
  case(2)
     write(ilog,*) 'list type  ld'
     !write(ilog,*) ' number of entries in hoc array ', &
     !&self%objectld%nhoc,' used out of ',self%n%nsizeh
     !     write(ilog,*) ' number of entries in conten array ', &
     ! &   self%objectld%nconten,' used out of ',self%n%nsizel, ' (max_size_list_array)'
  end select list_type

  !dbg  do jl=1,1000 ! dbg
  !dbg  write(ilog,*) jl ! dbg
  !dbg  end do ! dbg

end subroutine btree_dia
!---------------------------------------------------------------------
!> return node containing point vector starting from given node
function btree_newnode(self,kvec,k)

  ! This is only designed for octree.

  integer(ki4) :: btree_newnode !< returned function value
  !! arguments
  type(btree_t), intent(in) :: self   !< binary tree data
  integer(ki2), dimension(3), intent(in) :: kvec   !< position to be noded
  integer(ki4), intent(in) :: k   !< start node

  !! local
  character(*), parameter :: s_name='btree_newnode' !< subroutine name
  integer(ki4) :: jj   !< loop counter
  integer(ki2), dimension(3) :: ivec !< integer geobj data
  integer(ki2), dimension(3) :: inxyz !< top level position
  integer(ki2) :: ie !< integer element
  integer(ki2), dimension(3) :: iea !< array of ie
  integer(ki2) :: iesum !< integer sum
  integer(ki2) :: iemax !< integer max
  integer(ki4) :: inodn!< tree node number
  integer(ki4) :: inode !< tree node number
  integer(ki4) :: ik !< tree node number
  integer(ki4) :: ileaf !< tree index
  integer(ki4) :: ipow2  !< local variable
  integer(ki2) :: ichil  !< local variable
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: icorn  !< local variable
  integer(ki2) :: iext1  !< local variable
  integer(ki2) :: ishj !< local variable
  integer(ki2) :: idifj  !< local variable
  integer(ki2), dimension(3) :: ixyz  !< local variable
  integer(ki2):: iquant !< max scaled extent of tree (log)

  ivec=kvec
  !! Quick test for same or adjacent node
  ! check whether in input node
  iext=self%exten(:,self%desc(2,k))
  icorn=self%corner(:,k)
  iemax=0
  do jj=1,3
     iea(jj)=ishft( ieor(ivec(jj),icorn(jj)), -iext(jj) )
  end do
  iemax=max(iea(1),iea(2),iea(3))
  ! examine iemax
  if (iemax==0) then
     inode=k
     ! nodes at level iext(1)+1 (iext(1)=iext(2)=iext(3)) match
     ! bit x+1 corresponds to level x
  else if (iemax==1.AND.self%pter(1,k)>1) then
     ! in adjacent node
     inodn=k
     ipow2=1
     do jj=1,3
        idifj=ivec(jj)-icorn(jj)
        ishj=iea(jj)  * ipow2
        inodn=inodn+sign(ishj,idifj)
        ipow2=2*ipow2
     end do
     ileaf=self%pter(3,inodn)
     if (ileaf>0) then
        ! not lowest node, go down tree from here
        ik=inodn
        iext1=iext(1)

        ! loop over depth
        inode=0
        do i=1,iext1
           ichil=0
           do jj=1,3
              call mvbits(ivec(jj),iext1-i,1,ichil,jj-1)
           end do
           if (self%pter(3,ik) > 0 ) then
              ! node is not terminal, select appropriate child
              ik=self%pter(2,ik)+ichil
              cycle
              ! split direction

           else if (self%pter(3,ik) == -2 ) then
              ! empty node
              inode=ik
              exit
           else if (self%pter(3,ik) == -1 ) then
              ! leaf node
              inode=ik
              exit
           else if (self%pter(3,ik) == 0 ) then
              ! incomplete node
              inode=ik
              exit
           end if
        end do

     else
        inode=inodn
     end if
  else
     !! Simply start from top
     ! top level analysis
     ixyz=self%nxyz
     iquant=self%exten(2,self%desc(2,1))
     do jj=1,3
        inxyz(jj)=ishft(ivec(jj),-iquant)
        if (inxyz(jj).ge.ixyz(jj)) then
           call log_error(m_name,s_name,1,error_warning,'Object not found in binary tree')
           btree_newnode=-1
           return
        end if
     end do
     ichil=inxyz(1)+ixyz(1)*(inxyz(2)+inxyz(3)*ixyz(2))

     !! top level node
     ik=1
     inode=0
     if (self%pter(3,ik) > 0 ) then
        ! node is not terminal, select appropriate child
        ik=self%pter(2,ik)+ichil
     end if

     if (ik==1) then
        call log_error(m_name,s_name,2,error_fatal,'Binary tree has no entries')
     else

        !! loop over depth
        do i=1,self%ndepth-1
           ichil=0
           do jj=1,3
              call mvbits(ivec(jj),iquant-i,1,ichil,jj-1)
           end do
           if (self%pter(3,ik) > 0 ) then
              ! node is not terminal, select appropriate child
              ik=self%pter(2,ik)+ichil
              cycle
              ! split direction

           else if (self%pter(3,ik) == -2 ) then
              ! empty node
              inode=ik
              exit
           else if (self%pter(3,ik) == -1 ) then
              ! leaf node
              inode=ik
              exit
           else if (self%pter(3,ik) == 0 ) then
              ! incomplete node
              inode=ik
              exit
           end if
        end do

     end if
  end if

  ! check right node
  if (inode<=0) then
     call btree_dia(self)
     call log_error(m_name,s_name,5,error_fatal,'Binary tree does not exist')
  else
     ! corner
     icorn=self%corner(:,inode)
     ! extent
     iext=self%exten(:,self%desc(2,inode))
     iesum=0
     do jj=1,3
        ie=ishft( ieor(ivec(jj),icorn(jj)), -iext(jj) )
        iesum=iesum+ie
     end do
     if (iesum/=0) then

        call log_error(m_name,s_name,6,error_warning,'Object not found in binary tree')
        inode=-1
     end if
  end if
  btree_newnode=inode

end function btree_newnode
!---------------------------------------------------------------------
!> read unlabelled binary tree structure
subroutine btree_readnolab(self,numerics,kin)

  !! arguments
  type(btree_t), intent(out) :: self   !< binary tree data
  type(numerics_t), intent(inout) :: numerics !< local variable
  integer, intent(in) :: kin   !< input channel for binary tree data


  !! local
  character(*), parameter :: s_name='btree_readnolab' !< subroutine name
  integer(ki4) :: intr !< local variable
  integer(ki4) :: inpter  !< local variable
  integer(ki4) :: iexten  !< local variable

  ! read no in binary tree array
  read(kin,*,iostat=status) intr,inpter
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading binary tree sizes')
  end if
  self%nt=intr
  self%npter=inpter

  ! read binary tree arrays

  read(kin,*,iostat=status) &
  (self%pter(1:inpter,j),self%desc(:,j),self%corner(:,j),j=1,intr)
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error reading binary tree entries')
  end if

  ! read no in exten array
  read(kin,*,iostat=status) iexten
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error reading binary tree nexten')
  end if
  self%nexten=iexten

  ! read binary tree exten array
  read(kin,*,iostat=status) (self%exten(:,j),j=1,iexten)
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error reading binary tree exten')
  end if

  ! read binary tree ls array
  call ls_read(self%objectls,numerics,kin)

end subroutine btree_readnolab
!---------------------------------------------------------------------

  pure type(posvecl_t) function btree_floatvec(kvec)
  integer(ki2), dimension(3), intent(in) :: kvec  !< local variable
  type(posvecl_t) :: zvec !< local variable
  zvec%posvec(1)=float(kvec(1))
  zvec%posvec(2)=float(kvec(2))
  zvec%posvec(3)=float(kvec(3))
  btree_floatvec=zvec
  return
end function btree_floatvec
!---------------------------------------------------------------------
!> delete object
subroutine btree_delete(self)

  !! arguments
  type(btree_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='btree_delete' !< subroutine name
  integer(ki4) :: itlist  !<  local list type

  deallocate(self%pter)
  deallocate(self%desc)
  deallocate(self%corner)
  deallocate(self%exten)
  itlist=0 ! hardwired
  list_type: select case (itlist)
  case(0)
     call ls_delete(self%objectls)
  case(1)
     !  call li_delete(self%objectli)
  case(2)
     !  call ld_delete(self%objectld)
  end select list_type

end subroutine btree_delete

end module btree_m
