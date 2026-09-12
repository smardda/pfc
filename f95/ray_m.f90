module ray_m

  use log_m
  use const_numphys_h
  use const_kind_m
  use posang_h
  use posang_m
  use position_h
  use position_m
  use misc_m
  use beq_h
  use beq_m
  use radnsource_h
  use radnsource_m
  use geobjlist_h
  use geobjlist_m
  use btree_m
  use pcle_h
  use pcle_m

  implicit none
  private

! public subroutines
  public :: &
  ray_init, &  !< initialise object
  ray_initq, &  !< initialise object additionally with quantised information from beq
  ray_initfile,  & !< open file
  ray_readcon,  & !< read data from file
  ray_allocate, & !< allocate arrays associated with ray_t
  ray_launch,  & !< launch ray and sample points along ray
  ray_trackcoll, & !< loop along track of ray until collides with geometry
  ray_trackcount, & !< loop along track of ray counting steps to end
  ray_trackinlcfs, & !< loop along track of ray counting steps within LCFS
  ray_trackpower, &  !< loop along track of ray accumulating power
  ray_write, &  !< write out object
  ray_delete, & !< delete object
  ray_close !< close file

  type, public :: raynumerics_t
     integer(ki4) :: nsample = 0 !< user specified number of sample points along a ray
     real(kr4) :: termdistance = 0.0_kr4 !< maximum length of ray
     real(kr8) :: srcmin = 0.0_kr8 !< default 1, no source below this value, end track
     real(kr8) :: distfac = 0.0_kr8 !< default 1, multiplies maximum length based on mesh
     real(kr8) :: sepnfac = 0.0_kr8 !< default 1, multiplies minimum sample separation based on mesh
     real(kr8) :: abssepn = 0.0_kr8 !< minimum sample separation  (m)
     real(kr8) :: domeps = 0.0_kr8 !< start distance to ensure miss launch surface  (m)
     integer(ki4) :: distype = 0 !< default 1, use factor or absolute value if unity
     integer(ki4) :: quadtype = 0 !< default 1, default is trapezoidal rule
     integer(ki4) :: srctype = 0!< default 0, use flux values or source_threshold if unity
     character(len=80) :: formula = ''!< ray formula
     real(kr8) :: f = 0.0_kr8 !< power split (math variable name allowed)
  end type raynumerics_t

! type which defines/instantiates the object
  type, public :: ray_t
     type(raynumerics_t) :: n !< control  parameters
     real(kr8) :: power = 0.0_kr8 !< total power assigned to ray by flux function
     real(kr8) :: absdist = 0.0_kr8 !< maximum length (m) actually used
     real(kr8) :: dl = 0.0_kr8 !< step length between the sample points
     type(posveclis_t) :: cartposl !< cartesian list of sampled points along ray
     type(posveclis_t) :: polposl !< polar toroidal list of sampled points along ray
     real(kr8), dimension(:), allocatable :: fluxl !< list of flux values along ray
     integer(ki4) :: nquadrest = 0 !< estimated  number of sample (quadrature) points along a ray
     integer(ki4) :: nquadre  = 0!< actual  number of sample (quadrature) points along a ray
     real(kr8) :: psicen = 0.0_kr8 !< value of psi at centre of plasma
     real(kr8) :: psibdry = 0.0_kr8 !< value of psi at plasma boundary
     real(kr8) :: ximin = 0.0_kr8  !< minimum \f$ \xi \f$ value
     real(kr8) :: ximinp = 0.0_kr8 !< minimum \f$ \xi \f$ value \f$ +\epsilon \f$
     real(kr8) :: ximax = 0.0_kr8 !< maximum \f$ \xi \f$ value
     real(kr8) :: ximaxm = 0.0_kr8 !< maximum \f$ \xi \f$ value \f$ -\epsilon \f$
     real(kr8), dimension(3) :: domlen = 0.0_kr8 !< dimensions of domain (quantised)
     integer(ki4) :: nzets = 0 !< \f$ N_{\zeta s} \f$ such that \f$ \xi = N_{\zeta s} \zeta \f$
  end type ray_t

! private variables
  character(*), parameter :: m_name='ray_m' !< module name
  integer(ki4)  :: status  = 0 !< error status
  integer(ki4), save  :: ninray=5     !< control file unit number
  integer(ki4), save  :: noutray=6      !< output file unit number
  character(len=80), save :: controlfile = '' !< control file name
  character(len=80), save :: outputfile = '' !< output file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4)  :: ilog      !< for namelist dump after error
  character(len=29), parameter :: cfmtbv0c='(g14.7,", ",g14.7,", ",g14.7)'
  character(len=39), parameter :: cfmtbv1c='(g16.7,"",g16.7,"",g16.7)'
  real(kr8), parameter :: const_halfd=0.5_kr8 !< double precision half
  real(kr8), parameter :: const_twod=2.0_kr8 !< double precision two
  real(kr8), parameter :: const_zerod=0.0_kr8 !< double precision zero
  real(kr8), parameter :: const_oned=1.0_kr8 !< double precision unit
  real(kr8), parameter :: const_epsd=0.1_kr8 !< double precision tenth
  real(kr8), parameter :: const_eps2d=0.01_kr8 !< double precision hundredth
  real(kr8), parameter :: const_eps3d=0.001_kr8 !< double precision thousandth
  real(kr8), parameter :: const_eps4d=1.0e-4_kr8 !< double precision small
  real(kr8), parameter :: const_tend=10.0_kr8 !< double precision ten
  real(kr8), parameter :: const_ten8d=1.0e+8_kr8 !<  small infinity


  contains
!---------------------------------------------------------------------
!> initialise object
subroutine ray_init(self,psi,radn,psicen,psibdry)

  !! arguments
  type(ray_t), intent(inout) :: self !< module object
  type(spl2d_t), intent(inout) :: psi !< equilibrium flux data
  type(spl2d_t), intent(inout) :: radn !< radiation source array
  real(kr8), intent(in) :: psicen !< value of psi at centre of plasma
  real(kr8), intent(in) :: psibdry !< value of psi at boundary of plasma
  !! local
  character(*), parameter :: s_name='ray_init' !< subroutine name
  integer(ki4) :: nrayseg !< size of diagnostic arrays

  self%psicen=psicen
  self%psibdry=psibdry

  if (self%n%distype==0) then
     self%dl=self%n%distfac*min(psi%h1,psi%h2,radn%h1,radn%h2)
     self%absdist=self%n%sepnfac* &
 &   (  2*const_pid* max( psi%pos1(psi%n1p) , radn%pos1(radn%n1p) ) + &
 &   ( max( psi%pos2(psi%n2p)-psi%org2 , radn%pos2(radn%n2p)-radn%org2 ) )  )
     self%nquadrest=self%absdist/self%dl
  else
     self%dl=self%n%termdistance/self%n%nsample
     self%nquadrest=self%n%nsample
  end if
  nrayseg=self%nquadrest

  allocate(self%cartposl%pos(nrayseg), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  allocate(self%polposl%pos(nrayseg), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  allocate(self%fluxl(nrayseg), stat=status)
  call log_alloc_check(m_name,s_name,3,status)

end subroutine ray_init
!---------------------------------------------------------------------
!> initialise object additionally with quantised information from beq
subroutine ray_initq(self,ximin,ximinp,ximax,ximaxm,domlen,nzets)

  !! arguments
  type(ray_t), intent(inout) :: self !< module object
  real(kr8), intent(in) :: ximin  !< minimum \f$ \xi \f$ value
  real(kr8), intent(in) :: ximinp !< minimum \f$ \xi \f$ value \f$ +\epsilon \f$
  real(kr8), intent(in) :: ximax !< maximum \f$ \xi \f$ value
  real(kr8), intent(in) :: ximaxm !< maximum \f$ \xi \f$ value \f$ -\epsilon \f$
  real(kr8), dimension(3), intent(in) :: domlen !< dimensions of domain (quantised)
  integer(ki4), intent(in) :: nzets !< \f$ N_{\zeta s} \f$ such that \f$ \xi = N_{\zeta s} \zeta \f$

  !! local
  character(*), parameter :: s_name='ray_initq' !< subroutine name
  self%ximin=ximin
  self%ximinp=ximinp
  self%ximax=ximax
  self%ximaxm=ximaxm
  self%domlen=domlen
  self%nzets=nzets

end subroutine ray_initq
!---------------------------------------------------------------------
!> open file
subroutine ray_initfile(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer(ki4), intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='ray_initfile' !< subroutine name
  logical :: unitused !< flag to test unit is available

  if (trim(file)=='null') then
     call log_error(m_name,s_name,1,log_info,'null filename ignored')
     ninray=-1
     if (present(channel)) channel=ninray
     return
  end if

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  call misc_getfileunit(ninray)
  open(unit=ninray,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,2,error_fatal,'Cannot open control data file')
     stop
  end if
  if (present(channel)) then
     channel=ninray
  end if

end subroutine ray_initfile
!---------------------------------------------------------------------
!> read data from file
subroutine ray_readcon(selfn,channel)

  !! arguments
  type(raynumerics_t), intent(out) :: selfn !< type which data will be assigned to
  integer(ki4), intent(in),optional :: channel   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='ray_readcon' !< subroutine name
  integer(ki4) :: nsample !< default 100, user specified number of sample points along a ray
  real(kr8) :: termdistance !< default 10, user fixed maximum ray length (m)
  real(kr8) :: source_threshold !< default 1, no source below this value, end track
  real(kr8) :: distance_factor !< default 1, multiplies maximum length based on mesh
  real(kr8) :: separation_factor !< default 1, multiplies minimum sample separation based on mesh
  real(kr8) :: domain_epsilon !< default 10^{-4},  (m)
  real(kr8) :: absolute_distance !< maximum length (m) (NOT USED)
  real(kr8) :: absolute_separation !< minimum sample separation  (m)
  integer(ki4) :: distance_strategy !< default 1, use factor or absolute value if unity
  integer(ki4) :: integration_rule !< default 1, default is trapezoidal rule
  integer(ki4) :: source_threshold_strategy !< default 1, use flux values or source_threshold
  character(len=80) :: ray_formula !< local variable

  !! ray parameters
  namelist /rayparameters/ &
 &nsample, &
 &termdistance, &
 &source_threshold, &
 &distance_factor, &
 &separation_factor, &
 &domain_epsilon, &
 &   absolute_distance, &
 &   absolute_separation, &
 &   distance_strategy, &
 &   integration_rule, &
 &   source_threshold_strategy

     !! set default ray parameters
     nsample=100
     termdistance=10
     source_threshold=1
     distance_factor=1
     separation_factor=1
     domain_epsilon=const_eps4d
     absolute_distance=10
     absolute_separation=const_epsd
     distance_strategy=1
     integration_rule=1
     source_threshold_strategy=0

     if(present(channel).AND.channel/=0) then
        !! assume unit already open and reading infile
        ninray=channel
     end if

     !!read ray parameters
     read(ninray,nml=rayparameters,iostat=status)
     if(status/=0) then
        !!dump namelist contents to logfile to assist error location
        print '("Fatal error reading ray parameters")'
        call log_getunit(ilog)
        write(ilog,nml=rayparameters)
        call log_error(m_name,s_name,1,error_fatal,'Error reading ray parameters')
     end if

     call lowor(ray_formula,1,len_trim(ray_formula))
     !! check for valid data (TODO)

     !! store values
     selfn%nsample = nsample
     selfn%termdistance = termdistance
     selfn%srcmin=source_threshold
     selfn%distfac=distance_factor
     selfn%sepnfac=separation_factor
     selfn%domeps=domain_epsilon
     selfn%abssepn=absolute_separation
     selfn%distype=distance_strategy
     selfn%quadtype=integration_rule
     selfn%srctype=source_threshold_strategy

end  subroutine ray_readcon
!---------------------------------------------------------------------
!> allocate array associated with ray_t
subroutine ray_allocate(self)

     !! arguments
  type(ray_t), intent(inout) :: self !< module object
     !! local
  character(*), parameter :: s_name='ray_generic' !< subroutine name

     allocate(self%cartposl%pos(self%n%nsample), stat=status)
     call log_alloc_check(m_name,s_name,15,status)
     allocate(self%polposl%pos(self%n%nsample), stat=status)
     call log_alloc_check(m_name,s_name,15,status)
     allocate(self%fluxl(self%n%nsample), stat=status)
     call log_alloc_check(m_name,s_name,15,status)

end subroutine ray_allocate
!--------------------------------------------------------------------
!> loop over all launch directions (to be called when looping over triangles in antenna hemisphere
subroutine ray_launch(self,beq,centre,normal)

  !! arguments
  type(ray_t), intent(inout) :: self !< antenna object
  type(beq_t), intent(inout) :: beq !< local variable
  type(posvecl_t), intent(in) :: centre !< local variable
  type(posvecl_t), intent(in) :: normal !< local variable

  !! local
  character(*), parameter :: s_name='ray_launch' !< subroutine name
  integer(ki4) :: jpos !< loop counter over ray positions
  type(posang_t) :: posang !< posang data type to do coordinate conversions
  integer :: ucart !< unit for quick testing cartestian positions
  integer :: upol !< unit for quick testing polar toroidal positions
  integer :: upsi !< unit for quick testing flux (psi) positions
  real(kr8) :: z !< polar toroidal Z
  real(kr8) :: r !< polar toroidal R
  real(kr8) :: psi !< flux (psi) given by evaluation of psi(R,Z) spline

     !! test output .txts (just to see if the transform worked)
     !ucart = 33
     !open(unit = ucart, file = "testraycart2.txt")

     !  upol = 34
     !  open(unit = upol, file = "testraypolar2.txt")

     !  upsi = 35
     !  open(unit = upsi, file = "testrayflux2.txt")


     self%dl = self%n%termdistance/self%n%nsample ! set step length from user specified sample rate and distance

     do jpos = 1, self%n%nsample

        ! iterate to next sample point and store in posang
        posang%pos(1) = jpos*self%dl*normal%posvec(1) 
        posang%pos(2) = jpos*self%dl*normal%posvec(2)
        posang%pos(3) = jpos*self%dl*normal%posvec(3)

        ! select posang option
        posang%opt = 0  
        ! select posang units (0 for m)
        posang%units = 0 

        ! store cartesian sample points
        self%cartposl%pos(jpos)%posvec = centre%posvec + posang%pos 
        ! write out cartesian quick test
        !     write(ucart,*) self%cartposl%pos(jpos)%posvec 

        ! transform cartestion to polar toroidal
        call posang_invtfm(posang, 0) 

        ! store polar toroidal sample points
        self%polposl%pos(jpos)%posvec = posang%pos 
        ! write out polar quick test
        !     write(upol,cfmtbv0c) self%polposl%pos(jpos)%posvec 

        ! spline evaluation psi(R,Z) to determine value of psi of point on ray
        ! 

        ! set value of z for spline eval
        z = self%polposl%pos(jpos)%posvec(1) 
        ! set value of r for spline eval
        r = self%polposl%pos(jpos)%posvec(2) 
        ! carry out spline eval
        call spl2d_eval(beq%psi,z,r,psi) 
        ! store value of psi in list
        self%fluxl(jpos) = psi 
        !     write(upsi,cfmtbv0c) self%fluxl(jpos)

     end do

end subroutine ray_launch
!--------------------------------------------------------------------
!> loop along track of ray until collides with geometry
subroutine ray_trackcoll(self,psi,radn,gshadl,btree,centre,normal,kplot)

     !! arguments
  type(ray_t), intent(inout) :: self !< antenna object
  type(spl2d_t), intent(inout) :: psi !< equilibrium flux data
  type(spl2d_t), intent(inout) :: radn !< radiation source array
  type(geobjlist_t), intent(inout) :: gshadl !< shadow geobj list data
  type(btree_t), intent(in) :: btree   !< binary tree data
  type(posvecl_t), intent(in) :: centre  !< triangle centre
  type(posvecl_t), intent(in) :: normal !< triangle normal
  integer(ki4), intent(in) :: kplot   !< output channel for ray data structure

     !! local
  character(*), parameter :: s_name='ray_trackcoll' !< subroutine name
  integer(ki4) :: jpos !< loop counter over ray track
  type(posang_t) :: zposang !< posang data type to do coordinate conversions
  integer :: ucart !< unit for quick testing cartestian positions
  integer :: upol !< unit for quick testing polar toroidal positions
  integer :: upsi !< unit for quick testing flux (psi) positions
  real(kr8) :: zr0    !<   \f$ R_0(x,y,z) \f$
  real(kr8) :: zz0    !<   \f$ Z_0(x,y,z) \f$
  real(kr8) :: zxi0    !<   \f$ \xi_0(x,y,z) \f$
  real(kr8) :: zxi0q   !<  quantised  \f$ \xi_0(x,y,z) \f$
  real(kr8) :: zxi   !<  \f$ \xi(x,y,z) \f$
  real(kr8) :: zxiq   !<  quantised  \f$ \xi(x,y,z) \f$
  type(posnode_t) :: xo !< start pcle posn/node
  type(posnode_t) :: xn !< end pcle posn/node
  type(posnode_t) :: xm !< intermediate pcle posn/node
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  real(kr8) :: zrq    !< quantised \f$ R(x,y,z) \f$
  real(kr8) :: zzq    !< quantised \f$ Z(x,y,z) \f$
  real(kr8) :: zzeta    !<   \f$ \zeta \f$ coordinate of point on track
  real(kr8) :: zzeta0    !<   \f$ \zeta_0 \f$
  real(kr8) :: zzmax    !<   \f$ Z \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zpsimin    !<   \f$ \psi \f$ at minimum including starting value
  real(kr8) :: zpsimax    !<   \f$ \psi \f$ at maximum including starting value
  real(kr8) :: zpsi !< flux (psi) given by evaluation of psi(R,Z) spline
  type(geobj_t) :: iobj !< geo object
  !> Non-zero at fieldline termination
  !! - (-2) Line hit terminating plane
  !! - (-1) Line left octree
  !! - Positive - number of triangle hit
  integer(ki4) :: nobjhit=0  !< .
  real(kr8) :: zf1 !< fraction
  real(kr8) :: lenpath !< length of path in m
  type(posvecl_t) :: zpos  !<  position data
  type(posvecl_t) :: zposo  !<  position data
  type(posvecl_t) :: zposq   !< quantised position data
  type(posvecl_t) :: ztypec   !< use to do implicit type conversion
  type(posveclis_t) :: zposl   !< list of position data
  integer(ki4) :: ipos !< last value of loop counter for ray track
  integer(ki4) :: ik !< sector number
  integer(ki4) :: inode !< node number
  logical :: lcoll   !< collision on path

     ! convert start point to (R,Z) and m
     zposang%pos  = centre%posvec
     zposang%opt = 0
     zposang%units = 0
     call posang_invtfm(zposang, 0)
     zr0=zposang%pos(1)
     zz0=zposang%pos(2)
     xo%posvec=zposang%pos
     ! ensure all node components initialised
     xo%node=0
     xn%node=0
     xm%node=0

     ! find start point in tree
     if (.not.allocated(zposl%pos)) then
        allocate(zposl%pos(1), stat=status)
        call log_alloc_check(m_name,s_name,1,status)
     end if
     zzeta=zposang%pos(3)
     zzeta0=zzeta
     zxi=self%nzets*zzeta
     zposang%pos(3)=zxi
     !gshadl%quantfm%nqtfm=2
     ztypec%posvec=zposang%pos
     zposq=position_qtfm(ztypec,gshadl%quantfm)
     iobj%objtyp=1
     iobj%geobj=1
     zposl%pos(1)=zposq
     call btree_mfind(btree,iobj,zposl,inode)
     if (inode<0) then
        call log_error(m_name,s_name,2,error_warning,'cannot locate start point in tree')
        return
     end if

     ! loop to collision

     zxi0q=zposq%posvec(3)
     ipos=0
     lcoll=.FALSE.
     ik=0
     do jpos = 1, self%nquadrest
        nobjhit=0

        lenpath=(ipos+self%n%domeps)*self%dl
        ! next sample point and store in polars
        zposang%pos  = centre%posvec + lenpath*normal%posvec
        ! select posang option
        zposang%opt = 0
        ! select posang units
        zposang%units = 0
        ! store cartesian sample points
        self%cartposl%pos(jpos)%posvec = zposang%pos
        ! write out cartesian quick test
        !    write(ucart,*) self%cartposl%pos(jpos)%posvec
        ! transform cartesian to polar toroidal
        call posang_invtfm(zposang, 0)
        ! store polar toroidal sample points
        self%polposl%pos(jpos)%posvec = zposang%pos
        ! write out polar quick test
        !    write(upol,cfmtbv0c) self%polposl%pos(jpos)%posvec

        ! quantised coordinates
        zzeta = self%polposl%pos(jpos)%posvec(3)
        ! convert to xi
        zxi=self%nzets*zzeta
        zposang%pos(3)=zxi
        !gshadl%quantfm%nqtfm=2
        ztypec%posvec=zposang%pos
        zposq=position_qtfm(ztypec,gshadl%quantfm)
        xn%posvec=zposq%posvec
        ! store value of psi in list
        zrq = xn%posvec(1)
        zzq = xn%posvec(2)
        call spl2d_eval(psi,zrq,zzq,zpsi)
        self%fluxl(jpos) = zpsi
        xm=xn
        zxiq = xn%posvec(3)
        ! check for leaving periodic cell
        if (zxiq<self%ximin) then
           ! interpolate for xm=x(0) and check versus objects
           zf1=(self%ximin-zxi0q)/(zxiq-zxi0q)
           zf1=max(min(zf1,const_oned),const_zerod)
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           ! only terminate if trajectory is long enough
           if (jpos>1) then
              call pcle_move(xo,xm,jpos,gshadl,btree,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) exit
           xo=xm
           xo%posvec(3)=self%ximax
           xo%node=0
           xn%posvec(3)=min(self%ximaxm, xn%posvec(3)+self%domlen(3))
           ik=ik-1
        else if (zxiq>self%ximax) then
           ! interpolate for xm=x(2pi) and check versus objects
           zf1=(self%ximax-zxi0q)/(zxiq-zxi0q)
           zf1=max(min(zf1,const_oned),const_zerod)
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           ! only terminate if trajectory is long enough
           if (jpos>1) then
              call pcle_move(xo,xm,jpos,gshadl,btree,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) exit
           xo=xm
           xo%posvec(3)=self%ximin
           xo%node=0
           xn%posvec(3)=min(self%ximinp, xn%posvec(3)-self%domlen(3))
           ik=ik+1
        end if

        if (jpos>1) then
           call pcle_move(xo,xn,jpos,gshadl,btree,nobjhit)
           lcoll=(nobjhit/=0)
        end if
        if (lcoll) then
           xm=xn
           exit
        end if

        xo=xn
        zxi0q=zxiq

        ipos=ipos+1
     end do

     if (ipos==0) then
        call log_error(m_name,s_name,1,log_info,'track has no samples')
     else if (lcoll) then
        ! change entry to xm in self%polposl%pos,self%cartposl%pos (TO DO)
        ! back to polars
        ztypec%posvec=xm%posvec
        zpos=position_invqtfm(ztypec,gshadl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos%posvec(3)+2*const_pid*ik)/self%nzets
        zpos%posvec(3)=zzeta
        ! transform polar toroidal to cartesian and m
        zposang%pos=zpos%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,0)
        self%cartposl%pos(ipos)%posvec = zposang%pos
     else
        call log_error(m_name,s_name,2,log_info,'track terminated unexpectedly')
     end if
     self%nquadre=ipos-1
     !WA1 write(*,*) 'r,z,psi,length', zr,zz,zpsi,self%nquadre !WA1
     if (kplot>0.and.ipos>0) then
        zr=self%polposl%pos(ipos)%posvec(1)
        zz=self%polposl%pos(ipos)%posvec(2)
        zzeta=self%polposl%pos(ipos)%posvec(3)
        write(kplot,*)  zr0,zz0,zzeta0,zr,zz,zzeta
     end if
     !WA2 write(*,*) 'r0,z0,r,z', zr0,zz0,zr,zz !WA2

end subroutine ray_trackcoll
!--------------------------------------------------------------------
!> loop along track of ray counting steps to end
subroutine ray_trackcount(self,psi,radn,centre,normal)

     !! arguments
  type(ray_t), intent(inout) :: self !< antenna object
  type(spl2d_t), intent(inout) :: psi !< equilibrium flux data
  type(spl2d_t), intent(inout) :: radn !< radiation source array
  type(posvecl_t), intent(in) :: centre  !< triangle centre
  type(posvecl_t), intent(in) :: normal !< triangle normal

     !! local
  character(*), parameter :: s_name='ray_trackcount' !< subroutine name
  integer(ki4) :: jpos !< loop counter over ray track
  type(posang_t) :: posang !< posang data type to do coordinate conversions
  integer :: ucart !< unit for quick testing cartestian positions
  integer :: upol !< unit for quick testing polar toroidal positions
  integer :: upsi !< unit for quick testing flux (psi) positions
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zr0    !<   \f$ R_0(x,y,z) \f$
  real(kr8) :: zrmin    !<   \f$ R \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zrmax    !<   \f$ R \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  real(kr8) :: zz0    !<   \f$ Z_0(x,y,z) \f$
  real(kr8) :: zzmin    !<   \f$ Z \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zzmax    !<   \f$ Z \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zpsimin    !<   \f$ \psi \f$ at minimum including starting value
  real(kr8) :: zpsimax    !<   \f$ \psi \f$ at maximum including starting value
  real(kr8) :: zpsi !< flux (psi) given by evaluation of psi(R,Z) spline
  !real(kr8) :: zradn !< source (radn) given by evaluation of radn(R,Z) spline
  real(kr8) :: tracklen !< length of track
  integer(ki4) :: ipos !< last value of loop counter for ray track
  !integer(ki4) :: isrcmin !< above source threshold

     ! first point extends effective domain if necessary
     ! convert to (R,Z)
     posang%pos  = centre%posvec
     posang%opt = 0
     posang%units = 0
     call posang_invtfm(posang, 0)
     zr0=posang%pos(1)
     zz0=posang%pos(2)
     !! extend if necessary
     zrmin=min((1-const_eps2d)*zr0,radn%org1)
     zrmax=max((1+const_eps2d)*zr0,radn%pos1(radn%n1p))
     zzmin=min( (1-const_eps2d*sign(const_oned,zz0))*zz0,radn%org2 )
     zzmax=max( (1+const_eps2d*sign(const_oned,zz0))*zz0,radn%pos2(radn%n2p) )
     call spl2d_eval(psi,zr0,zz0,zpsi)
     if (self%psicen>self%psibdry) then
        zpsimin=min( (1-const_eps2d*sign(const_oned,zpsi))*zpsi , self%psibdry )
        zpsimax=max( (1+const_eps2d*sign(const_oned,zpsi))*zpsi , self%psicen )
     else
        zpsimin=min( (1-const_eps2d*sign(const_oned,zpsi))*zpsi , self%psicen )
        zpsimax=max( (1+const_eps2d*sign(const_oned,zpsi))*zpsi , self%psibdry )
     end if

     ! loop to determine where to stop tracking (source is an array, not a function)

     do jpos = 1, self%nquadrest
        ipos=jpos
        tracklen=jpos*self%dl
        ! next sample point and store in polars
        !WA posang%pos(1) = tracklen*normal%posvec(1)
        !WA posang%pos(2) = tracklen*normal%posvec(2)
        !WA posang%pos(3) = tracklen*normal%posvec(3)
        posang%pos  = centre%posvec + tracklen*normal%posvec

        ! select posang option
        posang%opt = 0
        ! select posang units (0 for m?)
        posang%units = 0

        ! store cartesian sample points
        self%cartposl%pos(jpos)%posvec = posang%pos
        ! write out cartesian quick test
        !    write(ucart,*) self%cartposl%pos(jpos)%posvec

        ! transform cartesian to polar toroidal
        call posang_invtfm(posang, 0)

        ! store polar toroidal sample points
        self%polposl%pos(jpos)%posvec = posang%pos
        ! write out polar quick test
        !    write(upol,cfmtbv0c) self%polposl%pos(jpos)%posvec

        !! spline evaluation psi(R,Z) to determine value of psi of point on ray
        ! set value of r for spline eval
        zr = self%polposl%pos(jpos)%posvec(1)
        ! set value of z for spline eval
        zz = self%polposl%pos(jpos)%posvec(2)

        ! Check leaving box given by limits of radn mesh and starting value
        if  ( ( (zr-zrmin) * (zrmax-zr) < 0 ) .OR. &
 &      ( (zz-zzmin) * (zzmax-zz) < 0 ) ) exit
        ! Check leaving range of psi (srctype=0 assumed)
        call spl2d_eval(psi,zr,zz,zpsi)
        if   ( (zpsi-zpsimin) * (zpsimax-zpsi) < 0 ) exit
        !call spl2d_eval(radn,r,z,zradn)
        ! store value of psi in list
        self%fluxl(jpos) = zpsi

     end do

     if (ipos==1) then
        call log_error(m_name,s_name,1,log_info,'track has no samples')
     end if
     self%nquadre=ipos-1
     !WA1 write(*,*) 'r,z,psi,length', zr,zz,zpsi,self%nquadre !WA1
     write(*,*) 'r0,z0,r,z', zr0,zz0,zr,zz !WA2

end subroutine ray_trackcount
!--------------------------------------------------------------------
!> loop along track of ray counting steps within LCFS
subroutine ray_trackinlcfs(self,psi,radn,centre,normal)

     !! arguments
  type(ray_t), intent(inout) :: self !< antenna object
  type(spl2d_t), intent(inout) :: psi !< equilibrium flux data
  type(spl2d_t), intent(inout) :: radn !< radiation source array
  type(posvecl_t), intent(in) :: centre  !< triangle centre
  type(posvecl_t), intent(in) :: normal !< triangle normal

     !! local
  character(*), parameter :: s_name='ray_trackinlcfs' !< subroutine name
  integer(ki4) :: jpos !< loop counter over ray track
  type(posang_t) :: posang !< posang data type to do coordinate conversions
  integer :: ucart !< unit for quick testing cartestian positions
  integer :: upol !< unit for quick testing polar toroidal positions
  integer :: upsi !< unit for quick testing flux (psi) positions
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zr0    !<   \f$ R_0(x,y,z) \f$
  real(kr8) :: zrmin    !<   \f$ R \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zrmax    !<   \f$ R \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  real(kr8) :: zz0    !<   \f$ Z_0(x,y,z) \f$
  real(kr8) :: zzmin    !<   \f$ Z \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zzmax    !<   \f$ Z \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zpsibdry    !<   \f$ \psi \f$ at boundary
  real(kr8) :: zsig !< +1 if psi increasing from centre (cf. beq_sigset), else (-1)
  real(kr8) :: zpsi !< flux (psi) given by evaluation of psi(R,Z) spline
  real(kr8) :: zpsin !< flux (psi) at next step
  real(kr8) :: zpsinorm !< normalisation of \f$ \psi \f$, i.e. a representative magnitude
  !real(kr8) :: zradn !< source (radn) given by evaluation of radn(R,Z) spline
  real(kr8) :: tracklen !< length of track
  integer(ki4) :: ipos !< last value of loop counter for ray track
  logical :: linlcfs  !< true when track inside LCFS

  ! psi at boundary adjusted to antenna always outside LCFS
  linlcfs=.FALSE.
  ! convert to (R,Z)
  posang%pos  = centre%posvec
  posang%opt = 0
  posang%units = 0
  call posang_invtfm(posang, 0)
  zr0=posang%pos(1)
  zz0=posang%pos(2)
  !! extend effective domain if necessary
  zrmin=min((const_oned-const_eps2d)*zr0,radn%org1)
  zrmax=max((const_oned+const_eps2d)*zr0,radn%pos1(radn%n1p))
  zzmin=min( (1-const_eps2d*sign(const_oned,zz0))*zz0,radn%org2 )
  zzmax=max( (1+const_eps2d*sign(const_oned,zz0))*zz0,radn%pos2(radn%n2p) )
  zsig=sign(const_oned,self%psibdry-self%psicen)
  zpsinorm=abs(self%psibdry-self%psicen)/2
  call spl2d_eval(psi,zr0,zz0,zpsi)
  zpsibdry=self%psibdry-const_eps2d*zsig*zpsinorm

  ! loop to determine where to stop tracking (source is an array, not a function)

  do jpos = 1, self%nquadrest
     ipos=jpos
     tracklen=jpos*self%dl
     ! next sample point and store in polars
     !WA posang%pos(1) = tracklen*normal%posvec(1) !WA
     !WA posang%pos(2) = tracklen*normal%posvec(2) !WA
     !WA posang%pos(3) = tracklen*normal%posvec(3) !WA
     posang%pos  = centre%posvec + tracklen*normal%posvec

     ! select posang option
     posang%opt = 0
     ! select posang units (0 for m)
     posang%units = 0

     ! store cartesian sample points
     self%cartposl%pos(jpos)%posvec = posang%pos
     !WA3 write out cartesian quick test
     !WA3    write(ucart,*) self%cartposl%pos(jpos)%posvec !WA3

     ! transform cartesian to polar toroidal
     call posang_invtfm(posang, 0)

     ! store polar toroidal sample points
     self%polposl%pos(jpos)%posvec = posang%pos
     !WA3 write out polar quick test
     !WA3    write(upol,cfmtbv0c) self%polposl%pos(jpos)%posvec !WA3

     !! spline evaluation psi(R,Z) to determine value of psi of point on ray
     ! set value of r for spline eval
     zr = self%polposl%pos(jpos)%posvec(1)
     ! set value of z for spline eval
     zz = self%polposl%pos(jpos)%posvec(2)

     ! Check leaving box given by limits of radn mesh and starting value
     if  ( (zr-zrmin)*(zrmax-zr)<0 .OR. (zz-zzmin)*(zzmax-zz)<0 ) exit
     ! Check entering and leaving LCFS
     call spl2d_eval(psi,zr,zz,zpsin)
     if  ( (zpsi-zpsibdry)*(zpsin-zpsibdry) < 0 ) then
        if ( (zpsi-zpsibdry)*zsig>0 ) then
           ! entry into LCFS
           zr0=zr
           zz0=zz
           linlcfs=.TRUE.
        else
           ! exit from both LCFS and loop
           if (linlcfs) exit
        end if
     end if
     !WA3 call spl2d_eval(radn,r,z,zradn) !WA3
     ! store value of psi in list
     zpsi=zpsin
     self%fluxl(jpos) = zpsi

  end do

  if (ipos==1 .OR. .NOT.linlcfs ) then
  ! never intersected LCFS
     self%nquadre=0
     call log_error(m_name,s_name,1,log_info,'track has no samples')
   else
  !WA1 write(*,*) 'r,z,psi,length', zr,zz,zpsi,self%nquadre !WA1
     write(*,*) 'r0,z0,r,z', zr0,zz0,zr,zz !WA2
   end if
     self%nquadre=ipos-1

end subroutine ray_trackinlcfs
!--------------------------------------------------------------------
!> loop along track of ray accumulating power
subroutine ray_trackpower(self,psi,radn,powtot)

     !! arguments
  type(ray_t), intent(inout) :: self !< antenna object
  type(spl2d_t), intent(inout) :: psi !< equilibrium flux data
  type(spl2d_t), intent(inout) :: radn !< radiation source array
  real(kr8) :: powtot !< integral of power along ray

     !! local
  character(*), parameter :: s_name='ray_trackpower' !< subroutine name
  integer(ki4) :: jpos !< loop counter over ray track
  type(posang_t) :: posang !< posang data type to do coordinate conversions
  real(kr8) :: z !< polar toroidal Z
  real(kr8) :: r !< polar toroidal R
  real(kr8) :: zpsi !< flux (psi) given by evaluation of psi(R,Z) spline
  real(kr8) :: zradn !< radiation power
  real(kr8) :: ztotsum !< integrated power along ray
  real(kr8) :: zlen   !< physical length (working)
  real(kr8), dimension(3)  :: zlenvec !< physical length vector (working)

     ! loop over track
     ztotsum=0
     do jpos = 1, self%nquadre-1
        ! set value of z for spline eval
        r = (self%polposl%pos(jpos+1)%posvec(1)+self%polposl%pos(jpos)%posvec(1))/2
        ! set value of r for spline eval
        z = (self%polposl%pos(jpos+1)%posvec(2)+self%polposl%pos(jpos)%posvec(2))/2
        ! carry out spline eval for radiation
        call spl2d_eval(radn,r,z,zradn)
        zradn = max(zradn, const_zerod)
        !! segment length
        zlenvec=self%cartposl%pos(jpos+1)%posvec-self%cartposl%pos(jpos)%posvec
        zlen=sqrt( max(const_zerod,zlenvec(1)**2+zlenvec(2)**2+zlenvec(3)**2) )
        ! integrate
        ztotsum=ztotsum+zradn*zlen
     end do
     powtot=ztotsum

end subroutine ray_trackpower
!---------------------------------------------------------------------
!> write ray data
subroutine ray_write(self,channel)

     !! arguments
  type(ray_t), intent(in) :: self   !< ray data structure
  integer(ki4), intent(in), optional :: channel   !< output channel for ray data structure

     !! local
  character(*), parameter :: s_name='ray_write' !< subroutine name
  integer(ki4) :: iout   !< output channel for ray data structure

     !! sort out unit
     if(present(channel)) then
        iout=channel
     else
        iout=noutray
     end if

     write(iout,*,iostat=status) 'ray_formula'
     call log_write_check(m_name,s_name,18,status)
     write(iout,*,iostat=status) self%n%formula
     call log_write_check(m_name,s_name,19,status)
     write(iout,*,iostat=status) 'f'
     call log_write_check(m_name,s_name,20,status)
     write(iout,*,iostat=status) self%n%f
     call log_write_check(m_name,s_name,21,status)

end subroutine ray_write
!---------------------------------------------------------------------
!> delete object
subroutine ray_delete(self)

     !! arguments
  type(ray_t), intent(inout) :: self !< module object
     !! local
  character(*), parameter :: s_name='ray_delete' !< subroutine name

     deallocate(self%cartposl%pos)
     deallocate(self%polposl%pos)
     deallocate(self%fluxl)

end subroutine ray_delete
!---------------------------------------------------------------------
!> close file
subroutine ray_close

     !! local
  character(*), parameter :: s_name='ray_close' !< subroutine name

     !! close file
     close(unit=ninray,iostat=status)
     if(status/=0)then
        !! error closing file
        print '("Fatal error: Unable to close control file, ",a)',controlfile
        call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
        stop
     end if

end subroutine ray_close

end module ray_m
