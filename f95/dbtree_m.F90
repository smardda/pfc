module dbtree_m

  use const_kind_m
  use const_numphys_h
  use ls_m
  use btree_m
  use li_m
  use ld_m
  use dbtree_h
  use log_m
  use misc_m
  use position_h
  use geobjlist_h
  use position_m
  use geobj_m

  implicit none
  private

! public subroutines
  public :: &
  dbtree_initfile,  & !< open file
  dbtree_readcon,  & !< read data from file
  dbtree_init,  & !< initialise object possibly using bbox to quantise
  dbtree_initfm, & !< initialise dbtree transformation (and nxyz)
  dbtree_addobj, & !< add object to nodal bin
  dbtree_box, & !< define boxes for node split
  dbtree_addnode, & !< update tree with split node and new nodes
  dbtree_find, & !< find node containing normalised point vector in BSP dbtree
  dbtree_mfind, & !< locate node containing point vector in tree (not BSP)
  dbtree_dia, &  !< object diagnostics to log file
  dbtree_initwrite, & !< open new file, making up name
  dbtree_write, &  !< write out object
  dbtree_writeg, &  !< write out object as gnuplot
  dbtree_writev, &  !< write out object as vtk
  dbtree_analyse, & !< analyse tree statistics
  dbtree_accum, & !< accumulate data in leaf, optionally just count leaves
  dbtree_geom, & !< calculate geometry of leaf
  dbtree_getreesca, & !< assign tree node values to array using mark indices
  dbtree_delete, & !< delete object
  dbtree_close, & !< close file
  dbtree_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='dbtree_m' !< module name
  integer  :: status   !< error status
  integer, save  :: nindbt=5     !< control file unit number
  integer, save  :: noutbo=6      !< output file unit number
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ii !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: kk !< loop counter
  integer(ki4) :: icall !< loop counter
  integer(ki4), save :: iwarn=0 !< count warnings for inability to split node
  integer :: ilog !< file unit for logging
  integer(ki4) :: idummy !< dummy

  contains
!---------------------------------------------------------------------
!> open file
subroutine dbtree_initfile(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer, intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='dbtree_initfile' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nindbt=i
  !! if (present(channel)) channel=i exit end if end do

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  call misc_getfileunit(nindbt)
  open(unit=nindbt,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open control data file')
     stop
  end if
  if (present(channel)) channel=nindbt

end subroutine dbtree_initfile
!---------------------------------------------------------------------
!> read data from file
subroutine dbtree_readcon(selfn,channel)

  !! arguments
  type(dbnumerics_t), intent(out) :: selfn !< type which data will be assigned to
  integer, intent(in),optional :: channel   !< input channel for object data structure

  integer(ki4) :: geometrical_type !< type of geometry
  real(kr4) :: min_tolerance !< numerical parameter
  real(kr4) :: max_tolerance !< numerical parameter
  integer(ki4) :: type_geobj_coord_scaling !< type of scaling of geobj x value
  integer(ki4) :: quantising_number !< local variable
  integer(ki4) :: no_geobj_records !< local variable
  integer(ki4) :: margin_type !< local variable
  integer(ki4) :: max_size_tree_array !<  max size of tree array
  integer(ki4) :: max_size_tree_pter_array !<  max size of tree pter array
  integer(ki4) :: max_size_exten_array !<  max size of exten array
  integer(ki4) :: max_size_list_array_hoc !<  max size of list array hoc
  integer(ki4) :: max_size_list_array !<  max size of list array
  integer(ki4) :: max_depth_tree !< max depth of tree
  !> type of binary tree
  !! for special top: 0 -default 1- use hxyz 2- use nxyz
  integer(ki4) :: type_tree_algorithm !< .
  integer(ki4) :: splitting_algorithm !< .
  integer(ki4) :: limit_geobj_in_bin !< maximum number of objects allowed any bin before splitting
  !> type of list structure used
  !! 0. numbered
  !! 1. simple, default
  !! 2. doubly-linked
  integer(ki4) :: type_list_structure !< .
  integer(ki4), dimension(3) :: top_tree_children !< top of tree children
  real(kr4), dimension(3) :: top_tree_spacings !< top of tree spacings
  integer(ki4) :: tree_type !< type of tree
  type(tfmdata_t) :: position_coord_tfm !< position \f$ x \f$ to \f$ x \f$

  !! local
  character(*), parameter :: s_name='dbtree_readcon' !< subroutine name

  !! dbtree parameters
  namelist /dbtreeparameters/ &
 &geometrical_type, &
 &min_tolerance, &
 &max_tolerance, &
 &quantising_number, &
 &no_geobj_records, &
 &margin_type, &
 &type_geobj_coord_scaling, &
 &max_size_tree_array, &
 &max_size_tree_pter_array, &
 &max_size_exten_array, &
 &max_size_list_array_hoc, &
 &max_size_list_array, &
 &max_depth_tree, &
 &type_tree_algorithm, &
 &splitting_algorithm, &
 &limit_geobj_in_bin, &
 &type_list_structure, &
 &top_tree_children, &
 &top_tree_spacings, &
 &tree_type

  !! set default dbtree parameters
  geometrical_type = 1
  min_tolerance = 1.0e-3
  max_tolerance = 1.0e-1
  quantising_number =4096
  no_geobj_records = 0
  margin_type = 0
  type_geobj_coord_scaling=2  ! allows for offset
  max_size_tree_array=2000000
  max_size_tree_pter_array=4
  max_size_exten_array=1000
  max_size_list_array_hoc=2
  max_size_list_array=2000000
  max_depth_tree=30
  !> Type of tree \n
  !! 1. BSP
  !! 2. octree
  !! 3. multi-octree (octree with special root node)
  tree_type=1
  type_tree_algorithm=1
  splitting_algorithm=0
  limit_geobj_in_bin=20
  type_list_structure=0
  top_tree_children=(/2,2,2/) !< for special top octree default
  top_tree_spacings=(/0.1,0.1,0.1/) !< for special top octree default

  if(present(channel).AND.channel/=0) then
     !! assume unit already open and reading infile
     nindbt=channel
  end if

  !!read dbtree parameters
  read(nindbt,nml=dbtreeparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading dbtree parameters")'
     call log_getunit(ilog)
     write(ilog,nml=dbtreeparameters)
     call log_error(m_name,s_name,5,error_fatal,'Error reading dbtree parameters')
  end if

  !! check for valid data
  if(geometrical_type<1) &
 &call log_error(m_name,s_name,16,error_fatal,'geometrical_type must be > 0')
  if(min_tolerance<0.0 .or. min_tolerance>max_tolerance) then
     call log_error(m_name,s_name,8,error_warning,'invalid min_tolerance, value reset')
     call log_value("min_tolerance",min_tolerance)
  end if
  if(max_tolerance<=0.0)then
     call log_value("max_tolerance",max_tolerance)
     call log_error(m_name,s_name,9,error_fatal,'Invalid max_tolerance, need value > 0')
  end if
  if(quantising_number<32) &
 &call log_error(m_name,s_name,10,error_fatal,'quantising_number must be > 31')

  if(type_geobj_coord_scaling<0) &
 &call log_error(m_name,s_name,11,error_fatal,'type_geobj_coord_scaling must be positive')

  !! check for valid data
  if(max_size_tree_array<1) &
 &call log_error(m_name,s_name,2,error_fatal,'max_size_tree_array must be > 0')
  if(max_size_tree_pter_array<1) &
 &call log_error(m_name,s_name,3,error_fatal,'max_size_tree_pter_array must be > 0')
  if(max_size_exten_array<1) &
 &call log_error(m_name,s_name,4,error_fatal,'max_size_exten_array must be > 0')
  if(max_size_list_array_hoc<1) &
 &call log_error(m_name,s_name,5,error_fatal,'max_size_list_array_hoc must be > 0')
  if(max_size_list_array<1) &
 &call log_error(m_name,s_name,6,error_fatal,'max_size_list_array must be > 0')
  if(max_depth_tree<1) &
 &call log_error(m_name,s_name,6,error_fatal,'max_depth_tree must be > 0')
  if(tree_type<1.OR.tree_type>3) &
 &call log_error(m_name,s_name,6,error_fatal,'tree_type must be > 0 and <=3')
  if(type_tree_algorithm<0.OR.type_tree_algorithm>2) &
 &call log_error(m_name,s_name,7,error_fatal,'type_tree_algorithm must be >= 0 and <=2')
  if(splitting_algorithm<0.OR.splitting_algorithm>3) &
 &call log_error(m_name,s_name,10,error_fatal,'splitting_algorithm must be >= 0 and <=3')
  if(limit_geobj_in_bin<1) &
 &call log_error(m_name,s_name,11,error_fatal,'limit_geobj_in_bin must be > 0')
  if(any(top_tree_children<1)) &
 &call log_error(m_name,s_name,12,error_fatal,'all top_tree_children must be > 0')
  if(any(top_tree_spacings<0.)) &
 &call log_error(m_name,s_name,14,error_fatal,'all top_tree_spacings must be > 0')

  call dbtree_quantise(selfn%nquante,quantising_number,type_geobj_coord_scaling)

  if (tree_type==1) then
     ! BSP
     if (max_depth_tree>3*selfn%nquante) then
        call log_error(m_name,s_name,21,error_warning,'depth too great for quantising_number,reset')
        !max_depth_tree=3*selfn%nquante
     end if
  else if (tree_type==2) then
     ! standard octree
     ! test depth not too great
     if (max_depth_tree>selfn%nquante) then
        call log_error(m_name,s_name,22,error_warning,'depth too great for quantising_number, reset')
        max_depth_tree=selfn%nquante
     end if
  end if

  call position_readcon(position_coord_tfm,nindbt)

  !! store values
  selfn%geomtype = geometrical_type
  selfn%ngeobj  = no_geobj_records
  selfn%mintolerance = min_tolerance
  selfn%mtype  = margin_type
  selfn%maxtolerance = max_tolerance
  selfn%nqtfm=type_geobj_coord_scaling
  selfn%nsize=max_size_tree_array
  selfn%nsizep=max_size_tree_pter_array
  selfn%nsizee=max_size_exten_array
  selfn%nsizeh=max_size_list_array_hoc
  selfn%nsizel=max_size_list_array
  selfn%ndepth=max_depth_tree
  selfn%nttype=tree_type
  selfn%nttalg=type_tree_algorithm
  selfn%splitalg=splitting_algorithm
  selfn%maxinbin=limit_geobj_in_bin
  selfn%ntlist=type_list_structure
  selfn%nxyz=top_tree_children
  selfn%hxyz=top_tree_spacings
  selfn%xtfm=position_coord_tfm

end subroutine dbtree_readcon
!---------------------------------------------------------------------
subroutine dbtree_quantise(nquante,kquante,kqtfm)

  !! arguments
  integer(ki4), intent(out) :: nquante   !< numeric control
  integer(ki4), intent(in):: kquante !< quantising number
  integer(ki4), intent(in):: kqtfm  !< not used variable

  !! local
  character(*), parameter :: s_name='dbtree_quantise' !< subroutine name
  real(kr4):: zrhmin  !< local variable
  integer(ki4):: i2  !< local variable
  integer(ki2):: iquante !< log base two of quantising_number

  ! quantisation
  i2=1
  do j=1,31
     i2=2*i2
     if (i2>=kquante) then
        ! log2 of quantising_number
        iquante=j
        zrhmin=float(i2)
        exit
     end if
  end do

  if (iquante.gt.ki2bits-1) then
     ! too large for type ki2 integers
     call log_error(m_name,s_name,1,error_warning,'quantising number too large, reset')
     iquante=ki2bits-2
  end if

  !! store values
  nquante  = iquante

end subroutine dbtree_quantise
!---------------------------------------------------------------------
!> initialise object possibly using bbox to quantise
subroutine dbtree_init(self,bbox,qtfmdata,numerics)

  !! arguments
  type(dbtree_t), intent(inout) :: self !< module object
  real(kr8), dimension(3,2), intent(inout) :: bbox !< bounding box corners
  type(quantfm_t), intent(in) :: qtfmdata !< geobj \f$ x \f$ to mesh units scaling
  type(dbnumerics_t), intent(in), optional :: numerics !< dbtree parameters

  !! local
  character(*), parameter :: s_name='dbtree_init' !< subroutine name
  integer(ki4) :: isize !< size of binary tree
  integer(ki4) :: isizep !< size of binary tree pter
  integer(ki4) :: isizee !< size of exten array
  integer(ki4) :: isizeh !< second size of list array
  integer(ki4) :: isizel !< size of list array
  integer(ki4) :: ingeobj  !< local variable

  if (present(numerics)) then
     self%n=numerics
  end if

  ! initialise tree storage
  isize=self%n%nsize
  isizep=self%n%nsizep
  isizee=self%n%nsizee
  isizeh=self%n%nsizeh
  isizel=self%n%nsizel

  !! allocate storage
  if(isize>0.AND.isizep>0) then
     allocate(self%pter(isizep,isize), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No data')
  end if

  if(isize>0) then
     allocate(self%desc(2,isize),self%corner(3,isize), stat=status)
     call log_alloc_check(m_name,s_name,3,status)
  else
     call log_error(m_name,s_name,4,error_fatal,'No data')
  end if

  if(isizee>0) then
     allocate(self%exten(3,isizee), stat=status)
     call log_alloc_check(m_name,s_name,5,status)
  else
     call log_error(m_name,s_name,6,error_fatal,'No data')
  end if

  ! initialise tfm and nxyz at top of tree
  call dbtree_initfm(self,qtfmdata,bbox)

  self%nt=1
  self%npter=isizep
  self%pter(1,1)=-1
  self%pter(2,1)=0
  self%pter(3,1)=-2
  self%pter(4,1)=0
  self%desc(1,1)=1
  self%desc(2,1)=1
  self%corner(:,1)=0
  self%minbdim=self%n%nquante
  self%maxallb=self%n%maxinbin

  list_type: select case (self%n%ntlist)
  case(0)
     call ls_init(self%objectls,isizeh,isizel)
     ! define dummy first record (0 objects) in list 0 and 2
     call ls_add(self%objectls,1,0,0)
     call ls_add(self%objectls,1,2,0)
     ! define list of objects
     ingeobj=self%n%ngeobj
     call ls_add(self%objectls,2,2,ingeobj)
     do j=1,ingeobj
        call ls_add(self%objectls,j+2,2,j)
     end do
  case(1)
     isizeh=1
     call li_init(self%objectli,isizeh,isizel)
  case(2)
     isizeh=2
     call ld_init(self%objectld,isizeh,isizel,1)
  end select list_type

  self%n%nsizeh=isizeh

end subroutine dbtree_init
!---------------------------------------------------------------------
!> initialise dbtree transformation using binbb
subroutine dbtree_initfm(self,qtfmdata,bbox)
  !! arguments
  type(dbtree_t), intent(inout) :: self !< dbtree data
  type(quantfm_t), intent(in) :: qtfmdata !< geobj \f$ x \f$ to mesh units scaling
  real(kr8), dimension(3,2), intent(inout) :: bbox !< bounding box corners

  !! local
  character(*), parameter :: s_name='dbtree_initfm' !< subroutine name
  integer(ki2):: iquante !< log base two of quantising_number
  real(kr4) :: margin !< margin added to bb
  real(kr4) :: minside !< min side of bb
  real(kr4) :: maxside !< max side of bb
  real(kr4) :: zhmin !< min side of bb for quantising
  real(kr4) :: zmind  !< local variable
  real(kr4) :: zpow2  !< local variable
  real(kr4), dimension(3) :: zlbb !< min corner of bb for quantising
  real(kr4), dimension(3) :: zubb !< max corner of bb for quantising
  real(kr4), dimension(3) :: zwork  !< local variable
  integer(ki2), dimension(3) :: ixyz  !< local variable
  integer(ki2) :: lmin  !< local variable
  integer(ki2) :: imxyz  !< local variable
  integer(ki2):: ittype !< type of tree
  integer(ki4) :: ipow2  !< local variable
  integer(ki4):: imtype  !< local variable
  real(kr4), dimension(3,2) :: zbinbb !< bb for geobj binning
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zposq !< local variable

  ! default quantisation
  ittype=self%n%nttype
  if (ittype==2) then
     iquante=self%n%nquante-1
  else
     iquante=self%n%nquante
  end if
  self%nexten=1
  self%exten(:,1)=(/iquante,iquante,iquante/)

  !  adjust bounding box
  zlbb=bbox(:,1)
  zubb=bbox(:,2)
  ! adjust depending on type of tree
  tree_types: select case (self%n%nttype)
  case(1)
     ! BSP
     ixyz=(/2,1,1/)
     self%maxchildn=2
  case(2)
     ! standard octree
     ixyz=(/2,2,2/)
     self%maxchildn=8
  case(3)
     !special top octree
     if (self%n%nttalg==0) then
        ! default, assumes hx=hy=hz unspecified
        !! define top level box sizes
        lmin=1
        minside=zubb(1)-zlbb(1)
        maxside=minside
        do i=2,3
           if(minside>zubb(i)-zlbb(i))then
              minside=zubb(i)-zlbb(i)
              lmin=i
           end if
           if(maxside<zubb(i)-zlbb(i))then
              maxside=zubb(i)-zlbb(i)
           end if

        end do

        margin=maxside*epsilon(maxside)
        do i=1,3
           ixyz(i)=int((margin+zubb(i)-zlbb(i))/minside,ki2)
        end do
        ! redefine upper bounds
        do i=1,3
           zubb(i)=zlbb(i)+float(ixyz(i))*minside
        end do
     else if (self%n%nttalg==1) then
        ! assumes hx,hy,hz specified
        zwork=(zubb-zlbb)/self%n%hxyz
        zmind=min(zwork(1),zwork(2),zwork(3))
        zpow2=2.
        do j=2,31
           ij=j
           if (zpow2.gt.zmind) exit
           zpow2=2.*zpow2
        end do
        self%n%ndepth=ij
        zwork=(zubb-zlbb+2.*margin)/(self%n%hxyz*zpow2)
        do i=1,3
           ixyz(i)=max(1,ceiling(zwork(i)))
           zubb(i)=zlbb(i)+float(ixyz(i))*self%n%hxyz(i)*zpow2
        end do
     else if (self%n%nttalg==2) then
        ! assumes nx,ny,nz specified
        ixyz=self%n%nxyz
     end if
     self%maxchildn=max(ixyz(1)*ixyz(2)*ixyz(3),8)
     ! test quantising not excessive
     imxyz=max(ixyz(1),ixyz(2),ixyz(3))
     ipow2=2
     do j=1,31
        ij=j
        if (ipow2.ge.imxyz) exit
        ipow2=2*ipow2
     end do
     if (self%n%nquante+ij.gt.ki2bits-2) then
        ! too large for type ki2 integers
        call log_error(m_name,s_name,3,error_warning,'quantising number too large, reset')
        iquante=ki2bits-2-ij
        self%n%nquante=iquante
        self%exten(:,1)=(/iquante,iquante,iquante/)
     end if
     ! test depth not too great
     if (self%n%ndepth>self%n%nquante) then
        call log_error(m_name,s_name,4,error_warning,'depth too great for quantising_number, reset')
        self%n%ndepth=self%n%nquante
     end if
  end select tree_types

  !! set bb (ignoring margins set above)
  self%binbb(:,1)=bbox(:,1)
  self%binbb(:,2)=bbox(:,2)

  !!set quantising to quantising transform
  self%quantfmq%nqtfm=2
  self%quantfmq%hmin=(/0.,0.,0./)
  zhmin=0.
  imtype=self%n%mtype

  ! set hmin vector
  if (self%n%nttype==1) then
     ! BSP
     self%quantfmq%hmin=(self%binbb(:,2)-self%binbb(:,1))/2**self%n%nquante
  else if (self%n%nttype==2) then
     ! octree
     self%quantfmq%hmin=(self%binbb(:,2)-self%binbb(:,1))/2**self%n%nquante
  else if (self%n%nttype==3) then
     ! multi-octree
     if (self%n%nttalg==0) then
        ! multi-octree, automatic quantisation
        zhmin=(self%binbb(lmin,2)-self%binbb(lmin,1))/2**self%n%nquante
        self%quantfmq%hmin=(/zhmin,zhmin,zhmin/)
     else if (self%n%nttalg==1) then
        ! multi-octree, using hxyz
        self%quantfmq%hmin=(self%binbb(:,2)-self%binbb(:,1))/(ixyz*2**self%n%nquante)
     else if (self%n%nttalg==2) then
        ! multi-octree, using nxyz
        self%quantfmq%hmin=(self%binbb(:,2)-self%binbb(:,1))/(ixyz*2**self%n%nquante)
     end if
  end if

  self%quantfmq%rhmin=1./self%quantfmq%hmin
  self%quantfmq%offvec=-self%binbb(:,1)*self%quantfmq%rhmin


  !!set quantising transform
  do j=1,2
     zpos%posvec=self%binbb(:,j)
     zposq=position_invqtfm(zpos,qtfmdata)
     zbinbb(:,j)=zposq%posvec
  end do

  self%quantfm%nqtfm=2
  self%quantfm%hmin=(/0.,0.,0./)
  zhmin=0.
  imtype=self%n%mtype

  ! set hmin vector
  if (self%n%nttype==1) then
     ! BSP
     self%quantfm%hmin=(zbinbb(:,2)-zbinbb(:,1))/2**self%n%nquante
  else if (self%n%nttype==2) then
     ! octree
     self%quantfm%hmin=(zbinbb(:,2)-zbinbb(:,1))/2**self%n%nquante
  else if (self%n%nttype==3) then
     ! multi-octree
     if (self%n%nttalg==0) then
        ! multi-octree, automatic quantisation
        zhmin=(zbinbb(lmin,2)-zbinbb(lmin,1))/2**self%n%nquante
        self%quantfm%hmin=(/zhmin,zhmin,zhmin/)
     else if (self%n%nttalg==1) then
        ! multi-octree, using hxyz
        self%quantfm%hmin=(zbinbb(:,2)-zbinbb(:,1))/(ixyz*2**self%n%nquante)
     else if (self%n%nttalg==2) then
        ! multi-octree, using nxyz
        self%quantfm%hmin=(zbinbb(:,2)-zbinbb(:,1))/(ixyz*2**self%n%nquante)
     end if
  end if

  self%quantfm%rhmin=1./self%quantfm%hmin
  self%quantfm%offvec=-zbinbb(:,1)*self%quantfm%rhmin

  self%n%nxyz=ixyz

end subroutine dbtree_initfm
!---------------------------------------------------------------------
!> initialise dbtree transformation (and nxyz)
subroutine dbtree_initfullm(self,bbox)
  !! arguments
  type(dbtree_t), intent(inout) :: self !< dbtree data
  real(kr8), dimension(3,2), intent(inout) :: bbox !< bounding box corners

  !! local
  character(*), parameter :: s_name='dbtree_initfullm' !< subroutine name
  integer(ki2):: iquante !< log base two of quantising_number
  real(kr4) :: margin !< margin added to bb
  real(kr4) :: minside !< min side of bb
  real(kr4) :: maxside !< max side of bb
  real(kr4) :: zhmin !< min side of bb for quantising
  real(kr4) :: zmind  !< local variable
  real(kr4) :: zpow2  !< local variable
  real(kr4), dimension(3) :: zlbb !< min corner of bb for quantising
  real(kr4), dimension(3) :: zubb !< max corner of bb for quantising
  real(kr4), dimension(3) :: zwork  !< local variable
  integer(ki2), dimension(3) :: ixyz  !< local variable
  integer(ki2) :: lmin  !< local variable
  integer(ki2) :: imxyz  !< local variable
  integer(ki2):: ittype !< type of tree
  integer(ki4) :: ipow2  !< local variable
  integer(ki4):: imtype  !< local variable

  ! default quantisation
  ittype=self%n%nttype
  if (ittype==2) then
     iquante=self%n%nquante-1
  else
     iquante=self%n%nquante
  end if
  self%nexten=1
  self%exten(:,1)=(/iquante,iquante,iquante/)

  !  adjust bounding box
  zlbb=bbox(:,1)
  zubb=bbox(:,2)
  ! adjust depending on type of tree
  tree_types: select case (self%n%nttype)
  case(1)
     ! BSP
     ixyz=(/2,1,1/)
     self%maxchildn=2
  case(2)
     ! standard octree
     ixyz=(/2,2,2/)
     self%maxchildn=8
  case(3)
     !special top octree
     if (self%n%nttalg==0) then
        ! default, assumes hx=hy=hz unspecified
        !! define top level box sizes
        lmin=1
        minside=zubb(1)-zlbb(1)
        maxside=minside
        do i=2,3
           if(minside>zubb(i)-zlbb(i))then
              minside=zubb(i)-zlbb(i)
              lmin=i
           end if
           if(maxside<zubb(i)-zlbb(i))then
              maxside=zubb(i)-zlbb(i)
           end if

        end do

        margin=maxside*epsilon(maxside)
        do i=1,3
           ixyz(i)=int((margin+zubb(i)-zlbb(i))/minside,ki2)
        end do
        ! redefine upper bounds
        do i=1,3
           zubb(i)=zlbb(i)+float(ixyz(i))*minside
        end do
     else if (self%n%nttalg==1) then
        ! assumes hx,hy,hz specified
        zwork=(zubb-zlbb)/self%n%hxyz
        zmind=min(zwork(1),zwork(2),zwork(3))
        zpow2=2.
        do j=2,31
           ij=j
           if (zpow2.gt.zmind) exit
           zpow2=2.*zpow2
        end do
        self%n%ndepth=ij
        zwork=(zubb-zlbb+2.*margin)/(self%n%hxyz*zpow2)
        do i=1,3
           ixyz(i)=max(1,ceiling(zwork(i)))
           zubb(i)=zlbb(i)+float(ixyz(i))*self%n%hxyz(i)*zpow2
        end do
     else if (self%n%nttalg==2) then
        ! assumes nx,ny,nz specified
        ixyz=self%n%nxyz
     end if
     self%maxchildn=max(ixyz(1)*ixyz(2)*ixyz(3),8)
     ! test quantising not excessive
     imxyz=max(ixyz(1),ixyz(2),ixyz(3))
     ipow2=2
     do j=1,31
        ij=j
        if (ipow2.ge.imxyz) exit
        ipow2=2*ipow2
     end do
     if (self%n%nquante+ij.gt.ki2bits-2) then
        ! too large for type ki2 integers
        call log_error(m_name,s_name,3,error_warning,'quantising number too large, reset')
        iquante=ki2bits-2-ij
        self%n%nquante=iquante
        self%exten(:,1)=(/iquante,iquante,iquante/)
     end if
     ! test depth not too great
     if (self%n%ndepth>self%n%nquante) then
        call log_error(m_name,s_name,4,error_warning,'depth too great for quantising_number, reset')
        self%n%ndepth=self%n%nquante
     end if
  end select tree_types


  !! set bb (with margin)
  self%binbb(:,1)=zlbb
  self%binbb(:,2)=zubb

  !!set quantising transform
  self%quantfm%nqtfm=2
  self%quantfm%hmin=(/0.,0.,0./)
  zhmin=0.
  imtype=self%n%mtype

  !! loop always executes at least once, even if imtype=margin_type=0
  do j=0,min(imtype,1)
     ! effect is to expand  box by 0, h or approx 3h+2h^2, depending on imtype
     if ( (self%n%nttype==3).AND.(self%n%nttalg==0) ) then
        ! multi-octree, automatic quantisation
        ! on second pass zhmin is non-zero,
        self%binbb(:,1)=self%binbb(:,1)-imtype*zhmin
        self%binbb(:,2)=self%binbb(:,2)+imtype*zhmin
     else
        self%binbb(:,1)=self%binbb(:,1)-imtype*self%quantfm%hmin
        self%binbb(:,2)=self%binbb(:,2)+imtype*self%quantfm%hmin
     end if
     !
     if (self%n%nttype==1) then
        ! BSP
        self%quantfm%hmin=(self%binbb(:,2)-self%binbb(:,1))/2**self%n%nquante
     else if (self%n%nttype==2) then
        ! octree
        self%quantfm%hmin=(self%binbb(:,2)-self%binbb(:,1))/2**self%n%nquante
     else if (self%n%nttype==3) then
        ! multi-octree
        if (self%n%nttalg==0) then
           ! multi-octree, automatic quantisation
           zhmin=(self%binbb(lmin,2)-self%binbb(lmin,1))/2**self%n%nquante
           self%quantfm%hmin=(/zhmin,zhmin,zhmin/)
        else if (self%n%nttalg==1) then
           ! multi-octree, using hxyz
           self%quantfm%hmin=(self%binbb(:,2)-self%binbb(:,1))/(ixyz*2**self%n%nquante)
        else if (self%n%nttalg==2) then
           ! multi-octree, using nxyz
           self%quantfm%hmin=(self%binbb(:,2)-self%binbb(:,1))/(ixyz*2**self%n%nquante)
        end if
     end if
  end do
  self%quantfm%rhmin=1./self%quantfm%hmin
  self%quantfm%offvec=-self%binbb(:,1)*self%quantfm%rhmin

  self%n%nxyz=ixyz

end subroutine dbtree_initfullm
!---------------------------------------------------------------------
!> add object to nodal bin
subroutine dbtree_addobj(self,obj,posl,ktry)
  !! arguments
  !> btree data \n
  !! For node k, pointer array entries are \n
  !! self%desc(1,k)   | Splitting direction (BSP) or meaningless (other) \n
  !! self%desc(2,k)   | Pointer to exten array of box sizes \n
  !! self%corner(:,k) | Smallest corner of box \n
  !! self%pter(1,k)   | Pointer to parent node \n
  !! self%pter(2,k)   | HOC for node list \n
  !! self%pter(3,k)   |  -2 (no entries), -1 (terminal)  or pointer to first child \n
  !! self%pter(4,k)   |  Number of objects in bin
  type(dbtree_t), intent(inout) :: self !< .
  type(geobj_t), intent(in) :: obj   !< obj definition
  type(posveclis_t), intent(in) :: posl !< posl vector
  integer(ki4), intent(inout) :: ktry !< try for node

  !! local
  character(*), parameter :: s_name='dbtree_addobj' !< subroutine name
  integer(ki4) :: inobj  !< object number
  integer(ki4),dimension(1) :: inobja  !< object number
  integer(ki4) :: ninode  !< number of entries in node
  integer(ki4)  :: iobj !< geo object
  integer(ki4), dimension(1)  :: iobja !< geo object
  type(geobj_t)  :: igeobj   !< local geobj definition
  integer(ki2), dimension(3) :: ixyz  !< local variable
  integer(ki2), dimension(3) :: iextn  !< local variable
  integer(ki2), dimension(3,8) :: icorna !< new corner array
  logical :: ilpsplit !< logical flag for splitting node (provisional)
  logical :: ilsplit !< logical flag for splitting node
  integer(ki4), dimension(8) :: iadra  !< local variable
  integer(ki4), dimension(8) :: ina  !< local variable
  integer(ki2), dimension(3) :: iord  !< local variable
  integer(ki4), dimension(:,:), save, allocatable :: inodea  !< local variable
  integer(ki4), dimension(:), save, allocatable :: inbox  !< local variable
  integer(ki4) :: innode !< number of proper entries in kout
  integer(ki4), dimension(3) :: ibsiz !< actual size of box
  integer(ki2) :: ibdim !< log_2 (smallest box dimension)
  integer(ki2) :: ichildn  !< number of children of node
  integer(ki2) :: id  !< direction in BSP case
  integer(ki2) :: ittype !< local variable
  integer(ki2) :: ittyp !< local variable
  integer(ki2), dimension(3) :: icorn  !< local variable
  integer(ki4) :: inode   !< local variable
  integer(ki2) :: i8   !< local variable
  integer(ki2) :: i2   !< local variable
  integer(ki2) :: i2p  !< local variable

  inobj=obj%geobj
  igeobj%objtyp=VTK_VERTEX

  ! locate node containing object
  inode=ktry
  if (self%n%nttype==1) then
     ! BSP
     call dbtree_find(self,obj,posl,inode)
  else
     call dbtree_mfind(self,obj,posl,inode)
  end if
  ! ignore any escapers
  if (inode==-1) then
     list_type: select case (self%n%ntlist)
     case(1)
        call li_countpp(self%objectli)
     case(2)
        call ld_countpp(self%objectld)
     end select list_type
     return
  end if

  ! add to linked list corresponding to node
  list_type2: select case (self%n%ntlist)
  case(1)
     call li_add(self%objectli,self%pter(2,inode),inobj)
  case(2)
     inobja(1)=inobj
     call ld_add(self%objectld,self%pter(2,inode),inobja)
  end select list_type2
  self%pter(4,inode)=self%pter(4,inode)+sign(1,self%pter(4,inode))
  ninode=self%pter(4,inode)
  self%maxallb=max(self%maxallb,abs(ninode))

  ilsplit=.FALSE.
  ! test for splitting
  ittype=self%n%nttype
  id=self%desc(1,inode)
  ilpsplit=(ninode>=self%n%maxinbin)
  if (ilpsplit) then
     ! provisionally split node
     ! octree split by default
     ixyz=(/2,2,2/)
     ittyp=ittype
     if (ittype==2.AND.inode==1) then
        ! head of tree
        ixyz=self%n%nxyz
     else if (ittype==1.AND.self%n%splitalg>0) then
        ittyp=12
     end if
     call dbtree_box(self,ittyp,ixyz,inode,idummy,iextn,icorna,ichildn,ibsiz,ibdim)

     if (ibdim<1) then
        ! run out of depth/descent, so do not in fact split
        iwarn=iwarn+1
        if (iwarn==1) then
           call log_error(m_name,s_name,2,error_warning,'Binary tree depth exceeded')
           call dbtree_dia(self)
        end if
        ! negate number of entries to mark box as unsplittable
        self%pter(4,inode)=-abs(self%pter(4,inode))
        self%minbdim=1
     else
        ilsplit=.TRUE.
        self%minbdim=min(self%minbdim,ibdim)
     end if
  end if

  if (ilsplit) then
     ! extract list of objects
     ! allocate array
     if (.NOT.allocated(inodea)) then
        allocate(inodea(1,self%n%maxinbin+1), inbox(self%n%maxinbin+1), stat=status)
        call log_alloc_check(m_name,s_name,3,status)
     end if
     ! fill with linked list
     list_type3: select case (self%n%ntlist)
     case(1)
        call li_fetch(self%objectli,self%pter(2,inode),inodea,innode)
     case(2)
        call ld_fetch(self%objectld,self%pter(2,inode),inodea,innode)
     end select list_type3
     !test     if (innode/=abs(self%pter(4,inode))) then !test
     !test        call log_error(m_name,s_name,3,error_warning,'mismatch in number of entries') !test
     !test     end if !test

     ! loop over children to define new nodes

     iadra=0
     ina=0
     ! determine which new box occupied by which objects
     ! first loop returns inbox array, second changes objectli
     ! after splitting direction  has been determined in BSP case
     loop_object: do j=1,self%n%maxinbin
        iobj=inodea(1,j)
        igeobj%geobj=iobj ! =ptr if not dealing with points
        loop_child: do i=1,ichildn
           ! loop over entries already in list
           icorn=icorna(:,i)
           if ( geobj_innbox( igeobj,posl,icorn,iextn )  ) then
              ina(i)=ina(i)+1
              inbox(j)=i
              exit ! remove if not dealing with points
           end if
        end do loop_child
     end do loop_object

     if (ittype==1) then
        ! BSP
        if (self%n%splitalg>0) then
           ! re-determine splitting direction
           call ordersplit(ina,iord)
           id=iord(self%n%splitalg)
        end if
        ixyz=(/1,1,1/)
        ixyz(id)=2
        call dbtree_box(self,ittype,ixyz,inode,id,iextn,icorna,ichildn,ibsiz,ibdim)
        ! now assign to 2 boxes indexed by i2p instead of 8
        ina=0
        loop_object1: do j=1,self%n%maxinbin
           i2=0
           i8=inbox(j)-1
           call mvbits(i8,id-1,1,i2,0)
           i2p=i2+1
           inbox(j)=i2p
           ina(i2p)=ina(i2p)+1
        end do loop_object1
     end if

     loop_object2: do j=1,self%n%maxinbin
        iobj=inodea(1,j)
        loop_child2: do i=1,ichildn
           ! loop over entries already in list
           if (inbox(j)==i) then
              list_type4: select case (self%n%ntlist)
              case(1)
                 call li_addin(self%objectli,iadra(i),iobj)
              case(2)
                 iobja(1)=iobj
                 call ld_addin(self%objectld,iadra(i),iobja)
              end select list_type4
           end if
        end do loop_child2
     end do loop_object2
     !dbg     write(*,*) 'icall',icall ! dbg
     !dbg     call dbtree_dia(self) ! dbg
     ! update tree with split node and new nodes using local arrays
     call dbtree_addnode(self,inode,iadra,ina,id,iextn,icorna,ichildn)
     !dbg1        continue ! dbg
  end if

  !dbg  icall=icall+1 ! dbg
  !dbg  write(*,*) icall ! dbg
  !if (allocated(inodea)) deallocate(inodea,inbox)

end subroutine dbtree_addobj
!---------------------------------------------------------------------
!> define boxes for node split
subroutine dbtree_box(self,kttype,kxyz,knode,kd,kextn,kcorna,kchildn,kbsiz,ibdim)

  !! arguments
  type(dbtree_t), intent(in) :: self   !< tree data
  integer(ki2), dimension(3), intent(in) :: kxyz   !< defines split
  integer(ki4), intent(in) :: kttype    !< local variable
  integer(ki4), intent(in) :: knode    !< local variable
  integer(ki2), dimension(3,8), intent(out) :: kcorna !< new corners
  integer(ki2), intent(out) :: kd   !< split direction/code
  integer(ki2), dimension(3), intent(out) :: kextn  !< local variable
  integer(ki4), dimension(3), intent(out) :: kbsiz !< actual size of box
  integer(ki2), intent(out) :: kchildn !< number of children
  integer(ki2), intent(out) :: ibdim !< log_2 (smallest box dimension)

  !! local
  character(*), parameter :: s_name='dbtree_box' !< subroutine name
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: icorn  !< local variable
  integer(ki2) :: iexth  !< local variable
  integer(ki2) :: ittype !< local variable
  integer(ki2) :: ichil !< local variable
  integer(ki4) :: jj   !< loop counter

  !! number of children
  ! type of tree
  ittype=kttype
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
  else if (ittype==12) then
     ! pseudo-octree
     ichil=2
     kd=3+kchildn
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
     ibdim=iexth
  else if (ichil==2) then
     kextn=iext-(/1,1,1/)
     kbsiz=2**kextn
     ibdim=kextn(1)
  else
     ! header node has extent of one octree (unless pseudo-octree)
     kextn=iext
     kbsiz=2**iext
     ibdim=kextn(1)
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

end subroutine dbtree_box
!---------------------------------------------------------------------
!> update tree with split node and new nodes
subroutine dbtree_addnode(self,knode,kadra,kna,kd,kextn,kcorna,kchildn)

  !! arguments
  type(dbtree_t), intent(inout) :: self   !< dbtree data
  integer(ki4), intent(in) :: knode   !< node to be split
  integer(ki4), dimension(8), intent(in) :: kadra !< array of new list hoc addresses
  integer(ki4), dimension(8), intent(in) :: kna !< number in each box
  integer(ki2), intent(in) :: kd   !< split direction
  integer(ki2), dimension(3), intent(in) :: kextn  !< new extents
  integer(ki2), dimension(3,8), intent(in) :: kcorna !< new corners
  integer(ki2), intent(in) :: kchildn !< local variable

  !! local
  character(*), parameter :: s_name='dbtree_addnode' !< subroutine name
  integer(ki4) :: intr !< count variable
  integer(ki4) :: iexten  !< local variable
  integer(ki4) :: imatch  !< local variable
  integer(ki2) :: ittype !< local variable
  integer(ki2) :: id !< local variable
  integer(ki4) :: jj !< loop counter

  ! entries of binary array
  ittype=self%n%nttype
  intr=self%nt
  iexten=self%nexten
  ! complete description of original node
  ! record splitting direction
  self%desc(1,knode)=kd
  self%pter(3,knode)=intr+1
  self%pter(4,knode)=0

  ! new nodes
  ! box extents, check space
  if (intr+kchildn>size(self%pter,2)) then
     call dbtree_dia(self)
     call log_error(m_name,s_name,1,error_fatal,'Binary tree size exceeded')
  end if
  ! test against existing entries
  imatch=0
  do kk=iexten,1,-1
     if ( all(kextn==self%exten(:,kk)) ) then
        imatch=kk
        exit
     end if
  end do

  if ( imatch==0 ) then
     ! add to exten array
     if (iexten+1>size(self%exten,2)) then
        call dbtree_dia(self)
        call log_error(m_name,s_name,2,error_fatal,'Binary tree exten size exceeded')
     end if
     iexten=iexten+1
     self%exten(:,iexten)=kextn
     imatch=iexten
  end if

  ! add dbtree arrays
  id=1+mod(kd,3) ! provisional only for BSP, else not used
  do jj=1,kchildn
     intr=intr+1
     self%desc(1,intr)=id
     self%desc(2,intr)=imatch
     self%corner(:,intr)=kcorna(:,jj)
     self%pter(1,intr)=knode
     self%pter(2,intr)=kadra(jj)
     if (kadra(jj)>0) then
        self%pter(3,intr)=-1
     else
        self%pter(3,intr)=-2
     end if
     self%pter(4,intr)=kna(jj)
  end do

  ! sizes of dbtree array
  self%nt=intr
  self%nexten=iexten

end subroutine dbtree_addnode
!---------------------------------------------------------------------
!> find node containing normalised point vector in BSP dbtree
subroutine dbtree_find(self,obj,poslis,knode)

  !! arguments
  type(dbtree_t), intent(in) :: self   !< binary tree data
  type(geobj_t), intent(in) :: obj !< geo object
  type(posveclis_t), intent(in) :: poslis   !< list of position data
  integer(ki4), intent(out) :: knode   !< node to be found

  !! local
  character(*), parameter :: s_name='dbtree_find' !< subroutine name
  integer(ki2) :: id !< local variable
  integer(ki2), dimension(3) :: ivec !< integer geobj data
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: icorn  !< local variable
  integer(ki2), dimension(3) :: iextn  !< local variable
  integer(ki4) :: jk   !< loop counter

  !! test against terminal input node
  if (knode>0.AND.self%pter(3,knode)<=0) then
     id=self%desc(1,knode)
     icorn=self%corner(:,knode)
     iext=self%exten(:,self%desc(2,knode))
     if ( geobj_innbox( obj,poslis,icorn,iext )  ) return
  end if

  !! loop over depth
  k=1
  knode=0
  do jk=1,self%n%ndepth+1
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

        if ( geobj_innbox( obj,poslis,icorn,iextn )  ) then
           ! in 0 position
           k=self%pter(3,k)
           cycle
        else
           ! in 1 position
           k=self%pter(3,k)+1
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
     call dbtree_dia(self)
     call log_error(m_name,s_name,1,error_fatal,'Binary tree too deep or may not exist')
  else
     ! corner
     icorn=self%corner(:,knode)
     ! extent
     iext=self%exten(:,self%desc(2,knode))
     if ( .not.geobj_innbox( obj,poslis,icorn,iext )  ) then
        call log_error(m_name,s_name,2,error_warning,'Object not found in binary tree')
        knode=-1
     end if
  end if

end subroutine dbtree_find
!---------------------------------------------------------------------
!> locate node containing point vector in tree (not BSP)
subroutine dbtree_mfind(self,obj,poslis,knode)

  !! arguments
  type(dbtree_t), intent(in) :: self   !< binary tree data
  type(geobj_t), intent(in) :: obj !< geo object
  type(posveclis_t), intent(in) :: poslis   !< list of position data
  integer(ki4), intent(inout) :: knode   !< node to be found

  !! local
  character(*), parameter :: s_name='dbtree_mfind' !< subroutine name
  type(posvecl_t) :: zpos !< local variable
  integer(ki4) :: inobj   !< geobj position
  integer(ki4) :: jj   !< loop counter
  integer(ki2), dimension(3) :: inxyz !< top level position
  integer(ki2) :: ichil  !< local variable
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: icorn  !< local variable
  integer(ki2), dimension(3) :: ixyz  !< local variable
  integer(ki2), dimension(3) :: ivec!< integer geobj data
  integer(ki2), dimension(3) :: iea !< array of ie
  integer(ki2) :: iesum !< integer sum
  integer(ki2) :: iquant !< max scaled extent of tree (log)

  ! check for point
  if (obj%objtyp/=1) then
     call log_error(m_name,s_name,1,error_fatal,'Called with wrong object type')
  end if

  if (knode>1.AND.self%pter(3,knode)<=0) then
     ! quick test if point is in a terminal input box (not top)
     icorn=self%corner(:,knode)
     iext=self%exten(:,self%desc(2,knode))
     if (geobj_innbox(obj,poslis,icorn,iext)) return
  end if

  ! top level analysis
  inobj=obj%geobj
  zpos=poslis%pos(inobj)
  ivec=int(zpos%posvec)
  ixyz=self%n%nxyz
  iquant=self%exten(2,self%desc(2,1))
  do jj=1,3
     inxyz(jj)=ishft(ivec(jj),-iquant)
     if (inxyz(jj)>=ixyz(jj)) then
        call log_error(m_name,s_name,3,error_warning,'Object not found in binary tree')
        call log_value('object number',obj%geobj)
        !call log_value('est node was',knode)
        !call log_value('iquant',iquant)
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
     k=self%pter(3,k)+ichil
  end if

  if (k==1) then
     knode=k
  else
     !! loop over depth
     do i=1,self%n%ndepth-1
        ichil=0
        do jj=1,3
           call mvbits(ivec(jj),iquant-i,1,ichil,jj-1)
        end do
        !D        write(*,*) i,k,ichil
        if (self%pter(3,k) > 0 ) then
           ! node is not terminal, select appropriate child
           k=self%pter(3,k)+ichil
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
     call log_error(m_name,s_name,5,error_warning,'Binary tree may be corrupt')
     call log_value('object number',obj%geobj)
     call log_value('node number',knode)
     call log_value('depth',i)
     !write(*,*) zpos !dbg
     call dbtree_dia(self)
     knode=-1
  else if (knode>1) then
     ! corner
     icorn=self%corner(:,knode)
     ! extent
     iext=self%exten(:,self%desc(2,knode))
     if ( .not.geobj_innbox( obj,poslis,icorn,iext )  ) then
        call log_error(m_name,s_name,6,error_warning,'Object not found in binary tree')
        call log_value('object number',obj%geobj)
        knode=-1
     end if
  end if

end subroutine dbtree_mfind
!---------------------------------------------------------------------
!> output to log file
subroutine dbtree_dia(self)

  !! arguments
  type(dbtree_t), intent(in) :: self !< module object
  !! local
  character(*), parameter :: s_name='dbtree_dia' !< subroutine name
  integer(ki2) :: iexmin !< minimum exten value
  integer(ki2) :: idescent !< binary depth (or descent) actually used
  integer(ki4) :: inleaf  !<  count number of leaves
  integer(ki4) :: inentry  !<  count number of entries
  integer(ki4) :: jl  !< local loop

  call log_getunit(ilog)
  write(ilog,*) ' tree type ', self%n%nttype
  write(ilog,*) ' split at top of tree ', self%n%nxyz
  write(ilog,*) ' type of tree algorithm ', self%n%nttalg
  write(ilog,*) ' size of tree pointer array or numerics%nsize ', &
 &self%nt,' used out of  ',self%n%nsize, ' (max_size_tree_array)'
  write(ilog,*) ' size of extents array or numerics%nsizee ', &
 &self%nexten,' used out of ',self%n%nsizee, ' (max_size_exten_array)'
  iexmin=minval( self%exten(1:3,1) )
  do jl=2,self%nexten
     iexmin=min( iexmin, minval( self%exten(1:3,jl) ) )
  end do
  idescent=self%n%nquante-iexmin
  write(ilog,*) ' descent of tree ', &
 &idescent, ' used out of ', self%n%nquante, ' (log_2(quantising number))'
  write(ilog,*) ' check descent of tree ', self%minbdim
  ! count leaves and entries
  inleaf=0
  inentry=0
  !! loop over all tree entries
  do jl=1,self%nt
     if (self%pter(3,jl)>0) cycle
     inleaf=inleaf+1
     inentry=inentry+abs(self%pter(4,jl))
  end do
  write(ilog,*) ' number of leaves ', inleaf
  write(ilog,*) ' check number of leaves ', self%nleaf
  write(ilog,*) ' number of entries ', inentry
  write(ilog,*) ' maximum number of entries in any node ', self%maxallb
  write(ilog,*) ' type of list structure ', self%n%ntlist
  list_type: select case (self%n%ntlist)
  case(0)
     write(ilog,*) 'list type  ls'
     write(ilog,*) ' number of entries in list array self%n%nsize or numerics%nsize ', &
 &   self%objectls%nlist,' used out of ',self%n%nsize, ' (max_size_tree_array)'
  case(1)
     write(ilog,*) 'list type  li'
     !write(ilog,*) ' number of entries in hoc array ', &
     !&self%objectli%nhoc,' used out of ',self%n%nsizeh
     write(ilog,*) ' number of entries in conten array ', &
 &   self%objectli%nconten,' used out of ',self%n%nsizel,  ' (max_size_list_array)'
  case(2)
     write(ilog,*) 'list type  ld'
     !write(ilog,*) ' number of entries in hoc array ', &
     !&self%objectld%nhoc,' used out of ',self%n%nsizeh
     write(ilog,*) ' number of entries in conten array ', &
 &   self%objectld%nconten,' used out of ',self%n%nsizel, ' (max_size_list_array)'
  end select list_type

  !dbg  do jl=1,1000 ! dbg
  !dbg  write(ilog,*) jl ! dbg
  !dbg  end do ! dbg

end subroutine dbtree_dia
!---------------------------------------------------------------------
!> open new file, making up name
subroutine dbtree_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer, intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='dbtree_initwrite' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  character(len=80) :: outputfile !< output file name

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then
  !! if (present(channel)) channel=i exit end if end do noutbo=i

  !! open file
  outputfile=trim(fileroot)//"_dbtree.out"
  call log_value("Control data file",trim(outputfile))
  call misc_getfileunit(noutbo)
  open(unit=noutbo,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutbo,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if
  if (present(channel)) channel=noutbo

end subroutine dbtree_initwrite
!---------------------------------------------------------------------
!> write dbtree data
subroutine dbtree_write(self,channel)

  !! arguments
  type(dbtree_t), intent(in) :: self   !< dbtree data structure
  integer, intent(in), optional :: channel   !< output channel for dbtree data structure

  !! local
  character(*), parameter :: s_name='dbtree_write' !< subroutine name
  integer :: iout   !< output channel for dbtree data structure

  !! sort out unit
  if(present(channel)) then
     iout=channel
  else
     iout=noutbo
  end if


end subroutine dbtree_write
!---------------------------------------------------------------------
!> write object data as gnuplot
subroutine dbtree_writeg(self,select,channel)

  !! arguments
  type(dbtree_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer, intent(in), optional :: channel   !< output channel for dbtree data structure

  !! local
  character(*), parameter :: s_name='dbtree_writeg' !< subroutine name
  integer :: iout   !< output channel for dbtree data structure

  call log_error(m_name,s_name,1,log_info,'gnuplot file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine dbtree_writeg
!---------------------------------------------------------------------
!> write object data as vtk
subroutine dbtree_writev(self,kchar,kclab,channel,geobjl)

  !! arguments
  type(dbtree_t), intent(in) :: self   !< binary tree data
  character(*),intent(in):: kchar !< control output geometry
  character(*),intent(in):: kclab !< control output data
  integer, intent(in) :: channel   !< output channel for binary tree data
  type(geobjlist_t), intent(in), optional :: geobjl   !< geobj list data

  !! local
  character(*), parameter :: s_name='dbtree_writev' !< subroutine name
  integer(ki4) :: intr  !< local variable
  integer(ki2), dimension(3) :: iorig  !< local variable
  integer(ki2), dimension(3) :: icell  !< local variable
  integer(ki2), dimension(3) :: ipos  !< local variable
  integer(ki2) :: iz  !< local variable
  integer(ki4) :: ioff  !< local variable
  integer(ki4) :: inode !< number of boxes in total
  integer(ki4) :: iptr2  !< local variable
  integer(ki4) :: islen1  !< local variable
  type(posvecl_t) :: zorig !< local variable
  type(posvecl_t) :: zorigt !< local variable
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zpost !< local variable
  character(len=80) :: iichar  !< local variable

  iichar=kchar
  if (present(geobjl)) then
     adjustplot_type: select case (kchar)
     case('hds')
        iichar='hds_lowest'
     case('hds_quantised')
        iichar='hds_quantised_lowest'
     end select adjustplot_type
  end if

  ! call to vfile_init must have been made
  iz=0
  !! write coordinates
  intr=self%nt
  write(channel,'(''DATASET UNSTRUCTURED_GRID'')')
  plot_geomtype: select case (iichar)
  case('hds')
     inode=intr
     write(channel,'(''POINTS  '',i8, '' float'')') 8*inode

     do j=1,inode

        iorig=self%corner(:,j)
        icell=2**self%exten(:,self%desc(2,j))
        zorig=btree_floatvec(iorig)
        zorigt=position_invqtfm(zorig,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zorigt%posvec
        ! write corner positions
        ipos=iorig+(/icell(1),iz,iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,icell(2),iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/icell(1),icell(2),iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,iz,icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/icell(1),iz,icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,icell(2),icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+icell
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec

     end do

  case('hds_lowest')
     ! only print lowest nodes in tree
     inode=0
     do j=1,intr
        iptr2=self%pter(3,j)
        if (iptr2>0) cycle
        inode=inode+1
     end do

     write(channel,'(''POINTS  '',i8, '' float'')') 8*inode

     do j=1,intr

        iptr2=self%pter(3,j)
        if (iptr2>0) cycle
        iorig=self%corner(:,j)
        icell=2**self%exten(:,self%desc(2,j))
        zorig=btree_floatvec(iorig)
        zorigt=position_invqtfm(zorig,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zorigt%posvec
        ! write corner positions
        ipos=iorig+(/icell(1),iz,iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,icell(2),iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/icell(1),icell(2),iz/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,iz,icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/icell(1),iz,icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+(/iz,icell(2),icell(3)/)
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec
        ipos=iorig+icell
        zpos=btree_floatvec(ipos)
        zpost=position_invqtfm(zpos,self%quantfm)
        write(channel,'(g14.7,1x,g14.7,1x,g14.7)') zpost%posvec

     end do

  case('hds_quantised')
     inode=intr
     write(channel,'(''POINTS  '',i8, '' unsigned_int'')') 8*inode

     do j=1,intr

        iorig=self%corner(:,j)
        icell=2**self%exten(:,self%desc(2,j))
        ! write corner positions

        write(channel,'(i8,1x,i8,1x,i8)') iorig,iorig+(/icell(1),iz,iz/), &
 &      iorig+(/iz,icell(2),iz/),iorig+(/icell(1),icell(2),iz/), &
 &      iorig+(/iz,iz,icell(3)/),iorig+(/icell(1),iz,icell(3)/), &
 &      iorig+(/iz,icell(2),icell(3)/),iorig+icell

     end do

  case default
     call log_error(m_name,s_name,1,error_warning,'Plot type unrecognised')
  end select plot_geomtype

  write(channel,'(''CELLS  '',i8,1x,i8)') inode, 9*inode

  ioff=0
  do j=1,inode
     ! write blocks
     write(channel,'("8 ",8i10)') (k,k=ioff,ioff+7)
     ioff=ioff+8
  end do

  !!write block types
  write(channel,'(''CELL_TYPES '', i8)') inode
  do j=1,inode
     write(channel,'("11")')
  end do

  if (kclab(1:1)/=' '.AND.present(geobjl)) then
     !!format statements for vtk file
     write(channel,'(''CELL_DATA '', i8)',iostat=status) inode
     call log_write_check(m_name,s_name,1,status)
     islen1=len_trim(kclab)
     write(channel,'(''SCALARS '',A,'' float 1'')',iostat=status)  kclab(1:islen1)
     call log_write_check(m_name,s_name,2,status)
     write(channel,'(''LOOKUP_TABLE default'')',iostat=status)
     call log_write_check(m_name,s_name,3,status)

     plot_datatype: select case (kclab)
     case('density')
        write(channel,cfmtbs1) (self%scal(j)/self%geom(j),j=1,self%nleaf)
     case('scalar')
        write(channel,cfmtbs1) (self%scal(j),j=1,self%nleaf)
     end select plot_datatype

  end if

  call log_error(m_name,s_name,1,log_info,'vtk file produced')

end subroutine dbtree_writev
!---------------------------------------------------------------------
!> analyse tree statistics
subroutine dbtree_analyse(self,kclab)

  !! arguments
  type(dbtree_t), intent(inout) :: self   !< binary tree data
  character(*),intent(in):: kclab !< control output data

  !! local
  character(*), parameter :: s_name='dbtree_analyse' !< subroutine name
  integer(ki4) :: inleaf  !<  count number of leaves

  count_leaves: select case (kclab)
  case('countleaf')
     ! count leaves
     inleaf=0
     !! first loop over all tree entries sets number of leaves
     do j=1,self%nt
        if (self%pter(3,j)>0) cycle
        inleaf=inleaf+1
     end do
     self%nleaf=inleaf
  end select count_leaves

end subroutine dbtree_analyse
!---------------------------------------------------------------------
!> accumulate data in leaf
subroutine dbtree_accum(self,kclab,geobjl)

  !! arguments
  type(dbtree_t), intent(inout) :: self   !< binary tree data
  character(*),intent(in):: kclab !< control output data
  type(geobjlist_t), intent(in), optional :: geobjl   !< geobj list data

  !! local
  character(*), parameter :: s_name='dbtree_accum' !< subroutine name
  integer(ki4) :: inleaf  !<  count number of leaves
  integer(ki4), dimension(:,:), allocatable :: inodea  !< local variable
  integer(ki4) :: innode !< number of proper entries in inodea
  integer(ki4) :: inodtest !< test
  integer(ki4) :: jk !< test
  type(geobj_t)  :: geobj   !< local geobj definition

  if (.NOT.allocated(self%scal)) then
     allocate(self%scal(1:self%nleaf), stat=status)
     !     else if (size(self%scal)/=self%nleaf) then
     !        deallocate(self%scal)
     !        allocate(self%scal(1:self%nleaf), stat=status)
  end if
  call log_alloc_check(m_name,s_name,1,status)
  self%scal=0

  if (.NOT.allocated(inodea)) then
     allocate(inodea(1,self%maxallb), stat=status)
     call log_alloc_check(m_name,s_name,2,status)
  end if

  !! second loop over all tree entries
  inleaf=0
  ! accumulate scalar
  do j=1,self%nt
     if (self%pter(3,j)>0) cycle
     inleaf=inleaf+1
     ! fill inodea with linked list
     list_type: select case (self%n%ntlist)
     case(1)
        call li_fetch(self%objectli,self%pter(2,j),inodea,innode)
     case(2)
        call ld_fetch(self%objectld,self%pter(2,j),inodea,innode)
     end select list_type
     !test        if (innode/=abs(self%pter(4,j))) then !test
     !test           call log_error(m_name,s_name,3,error_warning,'mismatch in number of entries') !test
     !test        end if !test
     ! test geobjl%posl(inodea(1,k)) is in node j
     !test        geobj%objtyp=VTK_VERTEX !test
     !test        inodtest=1 !test
     !test        do jk=1,innode !test
     !test            geobj%geobj=inodea(1,jk) !test
     !test            call dbtree_mfind(self,geobj,geobjl%posl,inodtest) !test
     !test           if (inodtest/=j) then !test
     !test              call log_error(m_name,s_name,4,error_warning,'mismatch in entry') !test
     !test           end if !test
     !test        end do !test
     ! now accumulate
     set_scalar: select case (kclab)
     case('scalar')
        do k=1,innode
           self%scal(inleaf)=self%scal(inleaf)+geobjl%obj(inodea(1,k))%weight
        end do
     case('number')
        self%scal(inleaf)=innode
     end select set_scalar
  end do

  deallocate(inodea)

end subroutine dbtree_accum
!---------------------------------------------------------------------
!> calculate area of leaf
subroutine dbtree_geom(self,kclab,geobjl)

  !! arguments
  type(dbtree_t), intent(inout) :: self   !< binary tree data
  character(*),intent(in):: kclab !< control output data
  type(geobjlist_t), intent(in), optional :: geobjl   !< geobj list data

  !! local
  character(*), parameter :: s_name='dbtree_geom' !< subroutine name
  integer(ki4) :: inleaf   !< local variable
  integer(ki2) :: iqtfm !< save transform type
  integer(ki2) :: iz !< local direction of normal \f$ 1,2,3 \f$
  integer(ki2) :: ix !< 1st direction in normal plane
  integer(ki2) :: iy !< 2nd direction in normal plane
  integer(ki2), dimension(3) :: icell  !< local variable
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zpost !< local variable


  inleaf=self%nleaf
  if (.NOT.allocated(self%geom)) allocate(self%geom(1:inleaf), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  self%geom=0

  if (kclab(3:6)=='area') then
     ! calculate area in direction normal to direction given by first entry
     ! (this is incorrect if there is only one leaf as first exten refers to
     ! size of virtual sub-boxes)
     read(kclab(1:1),'(i1)') iz
     ix=1+mod(iz,3)
     iy=1+mod(iz+1,3)
     ! save transform type, want only to scale
     iqtfm=self%quantfm%nqtfm
     self%quantfm%nqtfm=1
     inleaf=0
     do j=1,self%nt
        if (self%pter(3,j)>0) cycle
        inleaf=inleaf+1
        icell=2**self%exten(:,self%desc(2,j))
        zpos=btree_floatvec(icell)
        zpost=position_invqtfm(zpos,self%quantfm)
        self%geom(inleaf)=zpost%posvec(ix)*zpost%posvec(iy)
     end do
     ! restore
     self%quantfm%nqtfm=iqtfm
  end if

end subroutine dbtree_geom
!---------------------------------------------------------------------
!> assign tree node values to array using mark indices
subroutine dbtree_getreesca(self,kn,mark,kclab,psca)

  !! arguments
  type(dbtree_t), intent(inout) :: self   !< binary tree data
  integer(ki4), intent(in) :: kn !< arrays dimension
  integer(ki4), dimension(*), intent(in) :: mark !< array of indices
  character(*),intent(in):: kclab !< control output data
  real(kr4), dimension(:), allocatable, intent(out) :: psca !< scalar

  !! local
  character(*), parameter :: s_name='dbtree_getreesca' !< subroutine name

  if (kn>0) then
     allocate(psca(kn), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No data to plot')
  end if
  psca=0

  scalar_type: select case (kclab)
  case('scalar')
     do j=1,kn
        psca(j)=self%scal(mark(j))
     end do
  case('scalar2')
     do j=1,kn
        psca(j)=self%scal2(mark(j))
     end do
  case default
     do j=1,kn
        psca(j)=self%geom(mark(j))
     end do
  end select scalar_type

end subroutine dbtree_getreesca
!---------------------------------------------------------------------
!> close write file
subroutine dbtree_closewrite

  !! local
  character(*), parameter :: s_name='dbtree_closewrite' !< subroutine name

  !! close file
  close(unit=noutbo,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine dbtree_closewrite
!---------------------------------------------------------------------
!> delete object
subroutine dbtree_delete(self)

  !! arguments
  type(dbtree_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='dbtree_delete' !< subroutine name

  deallocate(self%pter)
  deallocate(self%desc)
  deallocate(self%corner)
  deallocate(self%exten)
  if (allocated(self%scal)) deallocate(self%scal)
  if (allocated(self%geom)) deallocate(self%geom)
  list_type: select case (self%n%ntlist)
  case(0)
     call ls_delete(self%objectls)
  case(1)
     call li_delete(self%objectli)
  case(2)
     call ld_delete(self%objectld)
  end select list_type

end subroutine dbtree_delete
!---------------------------------------------------------------------
!> close file
subroutine dbtree_close

  !! local
  character(*), parameter :: s_name='dbtree_close' !< subroutine name

  !! close file
  close(unit=nindbt,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine dbtree_close

!---------------------------------------------------------------------
!> return ordering of discrepancies in split
subroutine ordersplit(kna,kord)
  !! arguments
  integer(ki4), dimension(8), intent(in) :: kna !< number in each box
  integer(ki2), dimension(3), intent(out) :: kord  !< order of discrepancies, worst first
  !! local
  integer(ki4) :: isumh !< value of equal slit
  integer(ki4) :: idiscp  !< discrepancy in split
  integer(ki4) :: idiscpmax  !< max discrepancy in split
  integer(ki4) :: idiscpmin  !< min discrepancy in split
  integer(ki4), dimension(3) :: isuma  !< sum over lower boxes in direction of index
  integer(ki4) :: jl  !< local loop

  isumh=sum(kna(1:8))/2
  isuma=0
  ! x split
  do jl=1,4
     isuma(1)=isuma(1)+kna(2*jl-1)
  end do
  isuma(2)=kna(1)+kna(2)+kna(5)+kna(6)
  isuma(3)=sum(kna(1:4))

  idiscpmax=abs(isuma(1)-isumh)
  idiscpmin=idiscpmax
  kord=1
  do jl=2,3
     idiscp=abs(isuma(jl)-isumh)
     if (idiscp>idiscpmax) then
        idiscpmax=idiscp
        kord(1)=jl
     else if (idiscp<idiscpmin) then
        idiscpmin=idiscp
        kord(3)=jl
     end if
  end do

  if (kord(1)/=kord(3)) then
     kord(2)=6-kord(1)-kord(3)
  else
     ! all the same
     kord=(/1,2,3/)
  end if

end subroutine ordersplit

end module dbtree_m
