module bods_m

  use const_kind_m
  use log_m
  use const_numphys_h
  use bods_h
  use position_h
  use geobjlist_h
  use position_m
  use control_h
  use vcontrol_h
  use posang_h
  use ls_m
  use btree_m
  use geobj_m
  use query_m
  use posang_m
  use datline_h
  use datline_m
  use vfile_m
  use geobjlist_m

  implicit none
  private

! public subroutines
  public :: &
 &bods_init, & !< initialise vtktfm file bods
 &bods_initlist, & !< initialise datvtk file bods list only
 &bods_getlist, & !< get list of different bods in bods array
 &bods_write, & !< write out bods in separate files
 &bods_cumulate, & !< cumulate bods information
 &bods_cumulate2, & !< cumulate bods information and index
 &bods_delete !< delete bods information

! public types
!integer(ki4), parameter, public :: bods_made_up_data=1

! private variables
  character(*), parameter :: m_name='bods_m' !< module name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer :: status  !< status flag
  logical :: iltest !< logical flag
  character(len=256) :: vtkdesc !< descriptor line for vtk files

  contains
!---------------------------------------------------------------------
!> initialise vtktfm file bods
subroutine bods_init(self,geobjl,numerics)

  !! arguments
  type(bods_t), intent(inout) :: self !< bods
  type(geobjlist_t), intent(in) :: geobjl   !< geobj list data
  type(vnumerics_t)  :: numerics  !< numerical control parameters

  !! local
  character(*), parameter :: s_name='bods_init' !< subroutine name
  integer(ki4)  :: iobod !< number of objects in a bods
  integer(ki4)  :: mindx !< max number of bods in index

  !! allocate array
  iobod=geobjl%ng
  if( .NOT.allocated(self%list) ) then
     if( iobod>0 ) then
        allocate(self%list(iobod),stat=status)
        call log_alloc_check(m_name,s_name,1,status)
     else
        call log_error(m_name,s_name,2,error_fatal,'No bods data')
     end if
  end if
  self%nbod=iobod

  self%nindx=0
  self%maxbodsf=numerics%maxbodsf
  self%maxindx=numerics%maxindx

  if( .NOT.allocated(self%indx) ) then
     mindx=max(self%maxindx,1)
     allocate(self%indx(mindx),stat=status)
     call log_alloc_check(m_name,s_name,3,status)
  end if

end subroutine bods_init
!---------------------------------------------------------------------
!> initialise datvtk file bods list only
subroutine bods_initlist(self,geobjl,blockno)

  !! arguments
  type(bods_t), intent(out) :: self !< bods
  type(geobjlist_t), intent(in) :: geobjl   !< geobj list data
  integer(ki4), intent(in) :: blockno !< first block number

  !! local
  character(*), parameter :: s_name='bods_initlist' !< subroutine name
  integer(ki4)  :: iobod !< number of objects in a body

  !! allocate array
  iobod=geobjl%ng
  if(iobod>0) then
     allocate(self%list(iobod),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No bodies data')
  end if

  !! copy from geobjlist and integer convert
  if (allocated(geobjl%obj)) then
     do j=1,iobod
        self%list(j)=geobjl%obj(j)%weight+0.5
     end do
  else
     do j=1,iobod
        self%list(j)=blockno
     end do
  end if
  self%nbod=iobod

end subroutine bods_initlist
!---------------------------------------------------------------------
!> get list of different bods in bods array
subroutine bods_getlist(self,selfo)

  !! arguments
  type(bods_t), intent(in) :: self !< bods
  type(bods_t), intent(out) :: selfo !< bods short list

  !! local
  character(*), parameter :: s_name='bods_getlist' !< subroutine name
  integer(ki4)  :: inbod !< number of entries in self
  integer(ki4), dimension(:), allocatable :: iwork !< scratch bods list
  integer(ki4) :: ino !< bods short list dimension
  integer(ki4)  :: inew !< new entry in scratch list

  !! allocate array
  inbod=self%nbod
  if(inbod>0) then
     allocate(iwork(inbod),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No bods data')
  end if

  !! find different entries
  ino=1
  iwork(ino)=self%list(1)
  do j=2,inbod
     inew=1
     do k=1,ino
        if (self%list(j)==iwork(k)) then
           inew=0
           exit
        end if
     end do
     if (inew==1) then
        ino=ino+1
        iwork(ino)=self%list(j)
     end if
  end do
  !! copy to smaller output array
  allocate(selfo%list(ino),stat=status)
  call log_alloc_check(m_name,s_name,3,status)
  selfo%list=iwork(1:ino)
  selfo%nbod=ino

  deallocate(iwork)

end subroutine bods_getlist
!---------------------------------------------------------------------
!> write out bods in separate files
subroutine bods_write(self,geobjl,fileroot,kctyp,kcname,kheader)
  !! arguments
  type(bods_t), intent(in) :: self !< bods list
  type(geobjlist_t), intent(in) :: geobjl   !< geobj list data
  character(*),intent(in) :: fileroot !< file root
  character(*),intent(in) :: kctyp !< type of points compression (none)
  character(*),intent(in) :: kcname !< name of cell attribute
  integer(ki4), intent(in) :: kheader   !< header options

  !! local
  character(*), parameter :: s_name='bods_write' !< subroutine name
  type(geobjlist_t) :: glsmall   !< geobj list data small
  integer(ki4) :: infile   !< number of different bods files
  integer(ki4)  :: inpt !< number of points in generic geobjlist
  integer(ki4)  :: inobj !< number of objects in generic geobjlist
  integer(ki4)  :: innod !< number of nodal data in generic geobjlist
  integer(ki4)  :: inptsmall !< number of points in small geobjlist
  integer(ki4)  :: inobjsmall !< number of objects in small geobjlist
  integer(ki4)  :: innodsmall !< number of nodal data in small geobjlist
  integer(ki4)  :: iostart!< end points of object list if bods' objects contiguous
  integer(ki4)  :: ioend !< end points of object list if bods' objects contiguous
  integer :: iplot   !< output channel for bods list data
  integer(ki4) :: innd !< position of first entry for object in nodl
  integer(ki4) :: iobj !< local variable
  integer(ki4) :: inumpts !< length of object in nodl array
  integer(ki4) :: ityp !< local variable
  character(len=5) :: ic5  !< local variable
  character(len=80) :: icplot  !< local variable
  character(len=80) :: iched  !< local variable

  infile=max(1,maxval(self%list))

  ! write out all points
  inptsmall=geobjl%np

  do i=1,infile

     !! make up icplot from fileroot, iched from fileroot
     write(ic5,'(I5.5)') i
     icplot=trim(fileroot)//'_'//ic5
     iched='root was '//fileroot
     call geobjlist_makehedline(geobjl,iched,vtkdesc)
     call vfile_init(icplot,vtkdesc,iplot)

     !! construct small geobjl, first allocate storage

     inobjsmall=0
     innodsmall=0
     do j=1,geobjl%ng

        if (self%list(j)==i) then
           ioend=j
           inobjsmall=inobjsmall+1
           !     iobj=geobjl%obj2(j)%ptr
           ityp=geobjl%obj2(j)%typ
           inumpts=geobj_entry_table_fn(ityp)
           innodsmall=innodsmall+inumpts
        end if
     end do
     iostart=ioend+1-inobjsmall

     call geobjlist_iinit(glsmall,inptsmall,inobjsmall,innodsmall,2,1)

     !! fill small geobjl with data

     !! copy all points
     glsmall%posl=geobjl%posl

     inobj=0
     innd=1
     do j=1,geobjl%ng

        if (self%list(j)==i) then
           inobj=inobj+1
           iobj=geobjl%obj2(j)%ptr
           ityp=geobjl%obj2(j)%typ
           inumpts=geobj_entry_table_fn(ityp)
           do k=1,inumpts
              glsmall%nodl(innd-1+k)=geobjl%nodl(iobj-1+k)
           end do
           glsmall%obj2(inobj)%ptr=innd
           glsmall%obj2(inobj)%typ=ityp
           innd=innd+inumpts
        end if
     end do

     call geobjlist_ptcompress(glsmall)

     call geobjlist_writev(glsmall,'geometry',iplot)

     call vfile_iscalarwrite(self%list(iostart:ioend),inobjsmall,kcname,'CELL',iplot,kheader)

     call vfile_close

     call geobjlist_delete(glsmall)

  end do

end subroutine bods_write

!---------------------------------------------------------------------
!> cumulate bods information
subroutine bods_cumulate(self,geobjl,nfile,start,copy,kopt)

  !! arguments
  type(bods_t), intent(inout) :: self !< local variable
  type(geobjlist_t), intent(in) :: geobjl   !< geobj list data
  integer(ki4), intent(in) :: nfile   !< number of file
  integer(ki4), intent(in) :: start !< start number of copies required
  integer(ki4), intent(in) :: copy !< stop number of copies required
  integer(ki4), intent(in) :: kopt   !< unused option specifier

  !! local
  character(*), parameter :: s_name='bods_cumulate' !< subroutine name
  integer(ki4)  :: icopy !< actual number of copies required
  integer(ki4)  :: ibod !< number of objects in self
  integer(ki4)  :: inbod !< number of objects in geobjl
  integer(ki4), dimension(:), allocatable :: iwork !< integer work array
  integer(ki4) :: ing !< number of objects in geobjlist

  if (copy<=0) return
  icopy=copy+1-start
  !! allocate array
  ibod=self%nbod
  ing=geobjl%ng
  if(ing<0) call log_error(m_name,s_name,2,error_fatal,'No bods data')

  inbod=ibod+ing*icopy
  !set up replacement bods array and destroy old bods
  allocate(iwork(inbod), stat=status)
  call log_alloc_check(m_name,s_name,6,status)
  do j=1,ibod
     iwork(j)=self%list(j)
  end do
  i=ibod+1
  !! copy from geobjlist and integer convert
  do k=start,copy
     do j=1,ing
        if (allocated(geobjl%obj)) then
           iwork(i)=100*(self%maxbodsf*nfile+k)+geobjl%obj(j)%weight+0.5
        else
           iwork(i)=self%maxbodsf*nfile+k
        end if
        i=i+1
     end do
  end do
  deallocate(self%list)
  allocate(self%list(inbod), stat=status)
  call log_alloc_check(m_name,s_name,7,status)
  self%list=iwork
  self%nbod=inbod
  deallocate(iwork)

end subroutine bods_cumulate
!---------------------------------------------------------------------
!> cumulate bods information and index
subroutine bods_cumulate2(self,self2,geobjl,nfile,start,copy,kcall)

  !! arguments
  type(bods_t), intent(inout) :: self !< object variable
  integer(ki4), dimension(:), allocatable, intent(inout) :: self2 !< list of bods
  type(geobjlist_t), intent(in) :: geobjl   !< geobj list data
  integer(ki4), intent(in) :: nfile   !< number of file
  integer(ki4), intent(in) :: start !< start number of copies required
  integer(ki4), intent(in) :: copy !< stop number of copies required
  integer(ki4), intent(in) :: kcall   !< flag first call (negative suppresses index construction)

  !! local
  character(*), parameter :: s_name='bods_cumulate2' !< subroutine name
  integer(ki4)  :: icopy !< actual number of copies required
  integer(ki4)  :: ibod !< number of objects in self
  integer(ki4)  :: inbod !< number of objects in geobjl
  integer(ki4), dimension(:), allocatable :: iwork !< integer work array
  integer(ki4) :: ing !< number of objects in geobjlist
  integer(ki4) :: flbod !< flag change in bods number

  if (copy<=0) return
  icopy=copy+1-start
  !! allocate array
  ibod=self%nbod
  ing=geobjl%ng
  if(ing<0)call log_error(m_name,s_name,2,error_fatal,'No bods data')

  inbod=ibod+ing*icopy
  !set up replacement bods array and destroy old bods
  allocate(iwork(inbod), stat=status)
  call log_alloc_check(m_name,s_name,6,status)
  flbod=-1
  do j=1,ibod
     iwork(j)=self%list(j)
     if (kcall==0) then
        if (self%list(j)/=flbod) then
           flbod=self%list(j)
           self%nindx=self%nindx+1
           if (self%nindx<=self%maxindx) then
              self%indx(self%nindx)=self%maxbodsf+1
           else
              call log_error(m_name,s_name,3,error_fatal,'Bods index overflow, increase  max_bods_in_file')
           end if
        end if
     end if
  end do
  i=ibod+1
  !! update bods array
  do k=start,copy
     ! flags a change of file
     flbod=-1
     do j=1,ing
        if (kcall>=0) then
           if (self2(j)/=flbod) then
              flbod=self2(j)
              self%nindx=self%nindx+1
              if (self%nindx<=self%maxindx) then
                 self%indx(self%nindx)=self%maxbodsf*nfile+k
              else
                 call log_error(m_name,s_name,4,error_fatal,'Bods index overflow, increase  max_bods_in_file')
              end if
           end if
           iwork(i)=self%nindx
        else
           iwork(i)=self%maxbodsf*nfile+k
        end if
        i=i+1
     end do
  end do
  deallocate(self%list)
  allocate(self%list(inbod), stat=status)
  call log_alloc_check(m_name,s_name,7,status)
  self%list=iwork
  self%nbod=inbod
  deallocate(iwork)
  deallocate(self2)
  !dbgw write(*,*) 'nfile, self%nbod, self%indx = ',nfile,self%nbod,(self%indx(j),j=1,self%nbod) !dbgw

end subroutine bods_cumulate2
!---------------------------------------------------------------------
!> delete bods_t
subroutine bods_delete(self)

  !! arguments
  type(bods_t), intent(inout) :: self !< bods data

  deallocate(self%list)
  if (allocated(self%indx)) deallocate(self%indx)

end subroutine bods_delete


end module bods_m
