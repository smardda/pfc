module geoq_m

  use const_kind_m
  use log_m
  use const_numphys_h
  use geobjlist_h
  use position_m
  use control_h
  use dcontrol_h
  use ls_m
  use bcontrol_m
  use dcontrol_m
  use btree_m
  use geobj_m
  use query_m
  use position_h
  use fmesh_h
  use beq_h
  use posang_h
  use posang_m
  use geobjlist_m
  use spl2d_m
  use skyl_h
  use skyl_m
  use gfile_m
  use beq_m
  use geoq_h

  implicit none
  private

! public subroutines
  public :: &
  geoq_delete,   & !< delete  geoq data structure
  geoq_beqscale,   & !< scale geobjlist according to beq data
  geoq_read,   & !< read (vtk)  geoq data structure
  geoq_init,   & !< initialise geometry+field quantities
  geoq_objaddcon, & !< control addition of object(s) into model
  geoq_psilimiter,   & !< calculate limits of limiter object
  geoq_psisilh,   & !< calculate \f$ \psi \f$ of silhouette object
  geoq_dsilhcont,   & !< distance between silhouette and flux contour
  geoq_writev, &    !< write (vtk)  geoq data structure
  geoq_writeg    !< write (gnuplot)  geoq data structure

  private :: &
  geoq_skylspec, & !< special addition of skylight object(s) into model
  geoq_skyladd,   & !< skylight(s) into geobjlist for faster tracing
  geoq_objadd,   & !< object(s) into geobjlist for faster tracing
  geoq_skylpsi, & !< skylight defined using flux values
  geoq_skylpsi1,   & !< assist set up of skylight based on point flux values
  geoq_skylpsi2,   & !< assist set up of skylight based on centroid flux values
  geoq_skylcen, & !< skylight(s) defined using plasma centre line in PFR
  geoq_skylext !< skylight extent in flux terms


! public types
!public variables

! private types

! private variables
  character(*), parameter :: m_name='geoq_m' !< module name
  real(kr8), dimension(:), allocatable :: work1 !< 1D work array
  real(kr8), dimension(:), allocatable :: work1a !< 1D work array
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  integer   :: status   !< error status
  integer :: nin   !< input channel for geobj data
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  type(posveclis_t) :: rposl   !< list of position data
  integer(ki4) :: inpos !< number of position vectors to be read
  integer(ki4) :: inobj !< number of geobj records to be read
  integer(ki4) :: inls !< number of entries in list
  integer(ki4) :: idum !< dummy integer
  real(kr8), parameter :: deltal=1 !< silhoutte sampling length (mm)
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> delete geoq_t
subroutine geoq_delete(self)

  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  call geobjlist_delete(self%objl)
  call beq_delete(self%beq)

end subroutine geoq_delete
!---------------------------------------------------------------------
!> scale geobjlist according to beq data
subroutine geoq_beqscale(self,kt)

  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  integer(ki4), intent(in) :: kt    !<  type of transform

  !! local
  character(*), parameter :: s_name='geoq_beqscale' !< subroutine name
  real(kr8), dimension(3)   :: zpar !< scaling parameters
  type(geobjlist_t) :: zobjl !< geobj coord data and useful bits
  type(posang_t) :: zposang !< position and vector involving angles

  zpar=0
  if (self%beq%n%xiopt==2) then
     ! copy to workspace
     zobjl%np=self%objl%np
     !! allocate storage
     if(zobjl%np>0) then
        allocate(zobjl%posl%pos(zobjl%np),stat=status)
        call log_alloc_check(m_name,s_name,1,status)
     else
        call log_error(m_name,s_name,2,error_fatal,'No data')
     end if
     ! inert?2/8/18 zobjl%posl%np=self%objl%posl%np
     !! copy positions
     do j=1,zobjl%np
        zobjl%posl%pos(j)%posvec=self%objl%posl%pos(j)%posvec
     end do
     zobjl%nparam=self%objl%nparam
     zobjl%posl%nparpos=self%objl%posl%nparpos
     print '("number of geobj coordinates copied = ",i10)',zobjl%np
     call log_value("number of geobj coordinates copied ",zobjl%np)

     do j=1,zobjl%np
        ! transform positions to R-Z-zeta space
        zposang%pos=self%objl%posl%pos(j)%posvec
        zposang%opt=0 ; zposang%units=-3
        call posang_invtfm(zposang,0)
        !?7/8/18 self%objl%posl%pos(j)%posvec=zposang%pos
        zobjl%posl%pos(j)%posvec=zposang%pos
     end do

     ! scale in zeta
     zpar(2)=0
     zpar(3)=real(self%beq%nzets)
     call geobjlist_spectfm(zobjl,kt,zpar,3)
     zobjl%posl%nparpos(2)=2

     ! put back into geobjlist
     do j=1,zobjl%np
        self%objl%posl%pos(j)%posvec=zobjl%posl%pos(j)%posvec
     end do
     self%objl%nparam=zobjl%nparam
     self%objl%posl%nparpos=zobjl%posl%nparpos

     deallocate(zobjl%posl%pos)
  end if

end subroutine geoq_beqscale
!---------------------------------------------------------------------
!> read geobj coordinates
subroutine geoq_read(self,infileo,infileb)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  character(*),intent(in) :: infileo !< name of geobjlist input file
  character(*),intent(in) :: infileb !< name of beq input file

  !! local
  character(*), parameter :: s_name='geoq_read' !< subroutine name
  integer(ki4):: ipsibig !< flag whether psi overlarge by 2pi

  call geobjlist_read(self%objl,infileo)
  ipsibig=abs(self%beq%n%psibig)
  call log_value("psi 2pi too big if unity ",ipsibig)
  call beq_readequil(self%beq,infileb,self%beq%n)
  call log_error(m_name,s_name,70,log_info,'geoq input data files read')

end subroutine geoq_read
!---------------------------------------------------------------------
!> initialise geometry+field quantities
subroutine geoq_init(self)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  !! local
  character(*), parameter :: s_name='geoq_init' !< subroutine name

  if(self%beq%n%bdryopt/=3.AND.self%beq%n%bdryopt/=6) then
     ! calculate limits based on geometry, unless are silhouette options (3,6)
     call geoq_psilimiter(self)
  end if
  ! default no X-point, position R=0
  self%beq%rxpt=0
  self%beq%zxpt=0
  self%beq%psixpt=0

  ! psi boundary definition
  boundary_type: select case (self%beq%n%bdryopt)
  case (1,5,9)
     ! user defined
     self%beq%psibdry=self%beq%n%psiref
     self%beq%psiltr=self%beq%n%psiref
  case (2)
     ! boundary based on nodes of geometry only unless value in eqdsk over-rides
     if (beq_rsig()>0) then
        self%beq%psibdry=min(self%beq%psiqbdry,self%beq%psiltr)
     else
        self%beq%psibdry=max(self%beq%psiqbdry,self%beq%psiltr)
     end if
  case (3)
     ! boundary based on geometry sampled at rate given by parameter deltal
     call geoq_psisilh(self)
     !B     if (beq_rsig()>0) then
     !B        self%beq%psibdry=min(self%beq%psiqbdry,self%beq%psiltr)
     !B     else
     !B        self%beq%psibdry=max(self%beq%psiqbdry,self%beq%psiltr)
     !B     end if
     !B     ! geoq_override for specific project
     !B     if (BEQ_OVERRIDE_ITER) then
     !B        call log_error(m_name,s_name,49,log_info,'override for ITER')
     self%beq%psibdry=self%beq%psiltr
     !B     end if
  case (6)
     ! boundary based on geometry sampled at rate given by parameter deltal
     call geoq_psisilh(self)
     if (beq_rsig()>0) then
        self%beq%psibdry=min(self%beq%psiqbdry,self%beq%psiltr)
     else
        self%beq%psibdry=max(self%beq%psiqbdry,self%beq%psiltr)
     end if
  case (4,8)
     ! boundary from X-point
     call beq_psix(self%beq)
     !B     if (beq_rsig()>0) then
     !B        self%beq%psibdry=min(self%beq%psiqbdry,self%beq%psixpt)
     !B     else
     !B        self%beq%psibdry=max(self%beq%psiqbdry,self%beq%psixpt)
     !B     end if
     !B     ! geoq_override for specific project
     !B     if (BEQ_OVERRIDE_ITER) then
     !B     call log_error(m_name,s_name,50,log_info,'override for ITER')
     self%beq%psibdry=self%beq%psixpt
     self%beq%psiltr=self%beq%psixpt
     !B     end if
  case (7,10)
     ! boundary from X-point
     call beq_psix(self%beq)
     if (beq_rsig()>0) then
        self%beq%psibdry=min(self%beq%psiqbdry,self%beq%psixpt)
     else
        self%beq%psibdry=max(self%beq%psiqbdry,self%beq%psixpt)
     end if
  case (11,12)
     ! boundary from EQDSK regardless
     self%beq%psibdry=self%beq%psiqbdry
  case (13,14,15)
     ! boundary based on nodes of geometry only
     self%beq%psibdry=self%beq%psiltr
  end select boundary_type

  if (beq_rsig()*(self%beq%psiqbdry-self%beq%psibdry)<0) then
     call log_error(m_name,s_name,1,error_warning,'boundary psi outside separatrix')
  end if

  call log_error(m_name,s_name,80,log_info,'boundary psi value set')
  call log_value("psi boundary reference",self%beq%psibdry)
  call log_value("psi limiter or X-point",self%beq%psiltr)
  call log_value("psi plasma boundary",self%beq%psiqbdry)

end subroutine geoq_init
!---------------------------------------------------------------------
!> control addition of objects into model
subroutine geoq_objaddcon(self)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  !! local
  character(*), parameter :: s_name='geoq_objaddcon' !< subroutine name
  integer :: iin      !< local control file unit number
  type(dfiles_t) :: file !< file names specifying \f$ (R,Z) \f$ geometry
  type(dnumerics_t) :: numerics !< control numerics
  integer(ki4) :: jpla !<  number of object planes to add to geobjlist
  integer(ki2par) :: igcode !< integer scalar geometry code
  real(kr4) :: zetamin   !<  minimum \f$ \zeta \f$ of any point
  real(kr4) :: zetamax   !<  maximum \f$ \zeta \f$ of any point
  logical :: filedata !< data file specifies geometry

  filedata=.FALSE.

  if (sum(self%beq%n%objadd)>0.OR.self%beq%n%skylcen) then
     ! find angular extent of geometry
     call geobjlist_angext(self%objl,zetamin,zetamax)
     ! get unit for input
     call bcontrol_getunit(iin)
  end if

  ! flux based skylight before add in more geometry
  call geoq_skylspec(self,zetamin,zetamax)

  ! objects based on datvtkparameters input
  do j=1,size(self%beq%n%objadd)
     igcode=j-1
     desc_type: select case (igcode)
     case(GEOBJ_ABSORB, GEOBJ_INVISI, GEOBJ_ERRLOS, GEOBJ_CUTOUT)
        do jpla=1,self%beq%n%objadd(igcode)
           call dcontrol_readnum(numerics,iin,filedata)
           ! test
           if (numerics%descode-igcode/=0) then
              call log_error(m_name,s_name,1,error_warning,'Object description does not match')
              call log_value("Requested description code",igcode)
              call log_value("Found description code",numerics%descode)
           end if
           if (filedata) then
              call dcontrol_readprogfiles(file,iin)
              call dcontrol_readatfile(file,numerics)
              call dcontrol_lines2d(numerics)
           end if
           numerics%stang=zetamin
           numerics%finang=zetamax
           call geoq_objadd(self,numerics,igcode)
           call dcontrol_delete(numerics)
        end do
     case(GEOBJ_SKYLIT)
        do jpla=1,self%beq%n%objadd(igcode)
           call dcontrol_readnum(self%skyl%dn,iin,filedata)
           if (filedata) then
              call dcontrol_readprogfiles(file,iin)
              call dcontrol_readatfile(file,numerics)
              call dcontrol_lines2d(numerics)
           end if
           self%skyl%dn%stang=zetamin
           self%skyl%dn%finang=zetamax
           call geoq_skyladd(self,jpla)
        end do
     end select desc_type
  end do

end subroutine geoq_objaddcon
!---------------------------------------------------------------------
!> control special addition of skylight(s) into model
subroutine geoq_skylspec(self,zetamin,zetamax)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  real(kr4), intent(in) :: zetamin   !<  minimum \f$ \zeta \f$ of any point
  real(kr4), intent(in) :: zetamax   !<  maximum \f$ \zeta \f$ of any point

  !! local
  character(*), parameter :: s_name='geoq_skylspec' !< subroutine name
  integer :: iin      !< local control file unit number
  integer(ki4) :: jpla !<  number of skylight planes to add to geobjlist

  if (self%beq%n%skylpsi.OR.self%beq%n%skylcen) then
     ! define plasma centre line
     call beq_ctrax(self%beq)
  end if

  ! skylights based on flux
  if (self%beq%n%skylpsi) then
     call geoq_skylpsi(self)
  else
     self%skyl%n%skyltyp=0
  end if

  ! skylights based on plasma centre line
  if (self%beq%n%skylcen) then
     self%skyl%dn%stang=zetamin
     self%skyl%dn%finang=zetamax
     call geoq_skylcen(self)
  end if

end subroutine geoq_skylspec
!---------------------------------------------------------------------
!> object(s) into geobjlist for faster tracing
subroutine geoq_objadd(self,numerics,kgcode)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  type(dnumerics_t), intent(in) :: numerics !< control numerics
  integer(ki2par), intent(in) :: kgcode !<  object description

  !! local
  character(*), parameter :: s_name='geoq_objadd' !< subroutine name
  type(geobjlist_t) :: addobj !< skylight geobjlist
  integer(ki4) :: iopt=0 !< option for cumulate (no weights)

  call geobjlist_create3d(addobj,numerics,kgcode)

  call geobjlist_cumulate(self%objl,addobj,1,1,iopt,-1_ki2par)

  call geobjlist_delete(addobj)

end subroutine geoq_objadd
!---------------------------------------------------------------------
!> skylight(s) into geobjlist for faster tracing
subroutine geoq_skyladd(self,kpla)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  integer(ki4), intent(in) :: kpla !<  number of skylight plane to add to geobjlist

  !! local
  character(*), parameter :: s_name='geoq_skyladd' !< subroutine name
  type(geobjlist_t) :: skylobj !< skylight geobjlist
  integer :: iunit      !< local debug file unit number
  integer(ki4) :: iopt=0 !< option for cumulate (no weights)

  call geobjlist_create3d(skylobj,self%skyl%dn,GEOBJ_SKYLIT)
  if (self%beq%n%skyldbg>0) then
     iunit=8+kpla
     call geobjlist_writev(skylobj,'geometry',self%skyl%ndskyl(iunit))
  end if

  call geobjlist_cumulate(self%objl,skylobj,1,1,iopt,-1_ki2par)

  call geobjlist_delete(skylobj)

end subroutine geoq_skyladd
!---------------------------------------------------------------------
!> skylight defined using flux values
subroutine geoq_skylpsi(self)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  !! local
  character(*), parameter :: s_name='geoq_skylpsi' !< subroutine name
  integer(ki4) :: ilower !<  lower skylight type ( 0 or 1 )
  integer(ki4) :: iupper !<  upper skylight type ( 0 or 2 )
  integer(ki4) :: ityps !<  upper or lower skylight type
  integer(ki4) :: idoub !<  number of call to skylpsiN
  real(kr8), dimension(:,:), allocatable :: ctrackrz !< central track
  integer(ki4) :: nctrack !<  size of ctrackrz array
  integer(ki4), dimension(2) :: indbox !<  size of boxr,z arrays in skylpsiN
  integer(ki4) :: ibdim !<  second dimension of boxr,z arrays
  integer(ki4) :: jdoub !<  loop for possible lower and upper divertor
  real(kr8), dimension(2,4) :: zlts !<  limits for boxr,z arrays in skylpsiN
  real(kr8), dimension(4) :: zdeltas !<  delta for boxr,z arrays in skylpsiN
  real(kr8), dimension(2,2) :: zlt !<  limits for boxr,z arrays in skylpsiN
  real(kr8), dimension(2) :: zdelta !<  delta for boxr,z arrays in skylpsiN

  plot_skylight_type: select case (self%skyl%n%skyltyp)
  case default
     !     lower only
     ilower=1
     iupper=0
     indbox=self%skyl%n%nexts(1:2)
     ibdim=1
  case(2)
     !     both
     ilower=1
     iupper=2
     indbox(1)=max(self%skyl%n%nexts(1),self%skyl%n%nexts(3))
     indbox(2)=max(self%skyl%n%nexts(2),self%skyl%n%nexts(4))
     ibdim=2
  case(3)
     !     upper only
     ilower=0
     iupper=2
     indbox=self%skyl%n%nexts(3:4)
     ibdim=1
  end select plot_skylight_type
  ! initialise main (bin) arrays
  allocate(self%skyl%inboxr(indbox(1),ibdim),self%skyl%inboxz(indbox(1),ibdim),stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  self%skyl%inboxr=1.1*const_pushinf
  self%skyl%inboxz=0
  allocate(self%skyl%ouboxr(indbox(2),ibdim),self%skyl%ouboxz(indbox(2),ibdim),stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  self%skyl%ouboxr=1.1*const_pushinf
  self%skyl%ouboxz=0
  idoub=0
  do jdoub=1,2
     if (jdoub==1) then
        ityps=ilower
     else if (jdoub==2) then
        ityps=iupper
     end if
     if (ityps==0) cycle
     idoub=idoub+1
     call beq_ctrack1r(self%beq,ctrackrz,nctrack,ityps)
     if (self%beq%n%skyldbg>0) then
        call gfile_rwrite(ctrackrz(1,1),ctrackrz(1,2),nctrack,&
 &      'Track of extremum through plasma centre',self%skyl%ndskyl(1))
     end if
     call geoq_skylext(self,ctrackrz,nctrack,zlts,zdeltas,ityps)
     if (ityps==1) then
        indbox=self%skyl%n%nexts(1:2)
        zlt=zlts(:,1:2)
        zdelta=zdeltas(1:2)
     else if (ityps==2) then
        indbox=self%skyl%n%nexts(3:4)
        zlt=zlts(:,3:4)
        zdelta=zdeltas(3:4)
     end if
     if (self%skyl%n%control(1)==0) then
        call geoq_skylpsi1(self,ctrackrz,nctrack,indbox,zlt,zdelta,ityps,idoub)
     else
        call geoq_skylpsi2(self,ctrackrz,nctrack,indbox,zlt,zdelta,ityps,idoub)
     end if
     deallocate(ctrackrz)
  end do

end subroutine geoq_skylpsi
!---------------------------------------------------------------------
!> skylight(s) defined using plasma centre line in PFR
subroutine geoq_skylcen(self)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  !! local
  character(*), parameter :: s_name='geoq_skylcen' !< subroutine name
  integer(ki4) :: ilower !<  lower skylight type ( 0 or 1 )
  integer(ki4) :: iupper !<  upper skylight type ( 0 or 2 )
  integer(ki4) :: ityps !<  upper or lower skylight type
  integer(ki4) :: idoub !<  number of call to skyladd
  integer(ki4) :: jdoub !<  loop for possible lower and upper divertor
  real(kr8), dimension(2) :: zrz    !<  End point position \f$ (R,Z) \f

  plot_skylight_type: select case (self%skyl%n%skyltyp)
  case default
     !     lower only
     ilower=1
     iupper=0
  case(2)
     !     both
     ilower=1
     iupper=2
  case(3)
     !     upper only
     ilower=0
     iupper=2
  end select plot_skylight_type
  ! loop over possible plasma flux null configurations
  idoub=0
  do jdoub=1,2
     if (jdoub==1) then
        ityps=ilower
     else if (jdoub==2) then
        ityps=iupper
     end if
     if (ityps==0) cycle
     idoub=idoub+1
     call beq_ctrackpt(self%beq,self%beq%ctrackrz,self%beq%nctrack,&
 &   zrz,self%beq%psibdry,ityps)
     call skyl_dnumerics(self%skyl,self%beq%ctrackrz,self%beq%nctrack,&
 &   zrz(2),ityps,idoub)
     call geoq_skyladd(self,0)
  end do

end subroutine geoq_skylcen
!---------------------------------------------------------------------
!> skylight extent in flux terms
subroutine geoq_skylext(self,ctrackrz,nctrack,plts,pdeltas,ktyps)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  real(kr8), dimension(:,:), allocatable, intent(inout) :: ctrackrz !< central track
  integer(ki4), intent(in) :: nctrack !<  size of track array
  real(kr8), dimension(2,4), intent(out) :: plts !< psi limits
  real(kr8), dimension(4), intent(out) :: pdeltas !< psi bin size
  integer(ki4), intent(in) :: ktyps !<  lower (1) or upper (2) skylight type

  !! local
  character(*), parameter :: s_name='geoq_skylext' !< subroutine name
  type(posang_t) :: posang !< position and vector involving angle
  integer(ki4) :: idir !< ! direction of travel in \f$ Z \f$, (-1) in lower, (+1) in upper
  real(kr8) :: zdz    !<  bin spacing in \f$ Z \f$
  integer(ki4) :: icenz   !< local variable
  integer(ki4) :: intotz   !< local variable
  integer(ki4) :: iz   !< local variable
  integer(ki4) :: inz !<  dimension of bin arrays
  integer(ki4) :: izlti   !< local variable
  integer(ki4) :: izlto   !< local variable
  real(kr8) :: rlentol    !<  smaller than this length, assume on flux centreline
  real(kr8) :: distrack     !< local variable
  real(kr8) :: zpsi    !<  \f$ \psi \f$
  real(kr8), dimension(2) :: zrz    !<  Query point position \f$ (R,Z) \f$
  real(kr8) :: zr    !<  Query point \f$ R \f$
  real(kr8) :: zz    !<  Query point \f$ Z \f$
  real(kr8) :: psiinr   !<  inmost \f$ \psi \f$ of any point
  real(kr8) :: psioutr   !<  outmost \f$ \psi \f$ of any point
  real(kr8) :: zsign    !<  sign factor
  logical :: ilnwin !<  within eqdsk window
  integer(ki4) :: iupper   !< local variable
  real(kr8), dimension(:), allocatable :: inwallr !< inner wall limits in \f$ R \f$
  real(kr8), dimension(:), allocatable :: inwallz !< inner wall limit in \f$ Z \f$
  real(kr8), dimension(:), allocatable :: ouwallr !< outer wall limits in \f$ R \f$
  real(kr8), dimension(:), allocatable :: ouwallz !< outer wall limit in \f$ Z \f$

  ! direction of travel in z
  idir=2*ktyps-3

  if (self%skyl%n%extsopt(2+idir)==2.OR.self%skyl%n%extsopt(3+idir)==2) then
     zsign=beq_rsig()
     ! initialise main (bin) arrays
     ! dimensions
     zdz=self%skyl%geoml*self%skyl%n%ngeoml
     icenz=self%beq%n%zcen/zdz
     intotz=(self%beq%zmax-self%beq%zmin)/zdz
     rlentol=(self%beq%rmax-self%beq%rmin)*self%skyl%n%extsdel
     inz=intotz/2-idir*icenz+1+1
     allocate(inwallr(inz),inwallz(inz),ouwallr(inz),ouwallz(inz),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     inwallr=self%beq%rmin
     inwallz=1.1*const_pushinf
     ouwallr=self%beq%rmax
     ouwallz=1.1*const_pushinf

     izlti=inz
     izlto=inz
     do i=1,self%objl%np
        ! transform positions to R-Z-zeta space
        posang%pos=self%objl%posl%pos(i)%posvec
        posang%opt=0 ; posang%units=-3
        call posang_invtfm(posang,0)
        zr=posang%pos(1)
        zz=posang%pos(2)
        ! test within eqdsk window
        ilnwin= ( (zr-self%beq%rmin)*(zr-self%beq%rmax)<0 .AND. &
 &      (zz-self%beq%zmin)*(zz-self%beq%zmax)<0 )
        if (.NOT.ilnwin) cycle
        zrz=(/zr,zz/)
        call beq_ctrackcq(self%beq,ctrackrz,nctrack,zrz,distrack)
        iz=1+idir*(zz-self%beq%n%zcen)/zdz
        if (iz<1) cycle
        if (distrack<0) then
           if ( zr-inwallr(iz)>0 ) then
              inwallr(iz)=zr
              inwallz(iz)=zz
           end if
           if ( abs(distrack)<rlentol ) izlti=min(izlti,iz)
        else
           if ( ouwallr(iz)-zr>0 ) then
              ouwallr(iz)=zr
              ouwallz(iz)=zz
           end if
           if ( abs(distrack)<rlentol ) izlto=min(izlto,iz)
        end if
     end do
  end if

  if (self%beq%n%skyldbg>0) then
     call gfile_rwrite(inwallr,inwallz,izlti,'Estimate of inner wall silhouette',self%skyl%ndskyl(2))
     call gfile_rwrite(ouwallr,ouwallz,izlto,'Estimate of outer wall silhouette',self%skyl%ndskyl(3))
  end if

  if (self%skyl%n%extsopt(2+idir)==1) then
     plts(1,2+idir)=self%beq%psixpt-zsign*abs(self%skyl%n%psiexts(2+idir)/2)
     plts(2,2+idir)=self%beq%psixpt+zsign*abs(self%skyl%n%psiexts(2+idir)/2)
  else
     psiinr=zsign*const_pushinf
     psioutr=-zsign*const_pushinf
     do i=1,izlti
        zz=inwallz(i)
        if (zz>const_pushinf) cycle
        if ( (zz-self%skyl%n%zextmin)*(zz-self%skyl%n%zextmax)>0 ) exit
        zr=inwallr(i)
        if ( (zr-self%skyl%n%rextmin)*(zr-self%skyl%n%rextmax)>0 ) exit
        call spl2d_evaln(self%beq%psi,zr,zz,1,zpsi)
        if ( (zpsi-psiinr)*zsign<0 ) psiinr=zpsi
        if ( (zpsi-psioutr)*zsign>0 ) psioutr=zpsi
     end do
     plts(1,2+idir)=psiinr
     plts(2,2+idir)=psioutr
  end if

  if (self%skyl%n%extsopt(3+idir)==1) then
     plts(1,3+idir)=self%beq%psixpt-zsign*abs(self%skyl%n%psiexts(3+idir)/2)
     plts(2,3+idir)=self%beq%psixpt+zsign*abs(self%skyl%n%psiexts(3+idir)/2)
  else
     psiinr=zsign*const_pushinf
     psioutr=-zsign*const_pushinf
     do i=1,izlto
        zz=ouwallz(i)
        if (zz>const_pushinf) cycle
        if ( (zz-self%skyl%n%zextmin)*(zz-self%skyl%n%zextmax)>0 ) exit
        zr=ouwallr(i)
        if ( (zr-self%skyl%n%rextmin)*(zr-self%skyl%n%rextmax)>0 ) exit
        call spl2d_evaln(self%beq%psi,zr,zz,1,zpsi)
        if ( (zpsi-psiinr)*zsign<0 ) psiinr=zpsi
        if ( (zpsi-psioutr)*zsign>0 ) psioutr=zpsi
     end do
     plts(1,3+idir)=psiinr
     plts(2,3+idir)=psioutr
  end if

  do j=2,3
     pdeltas(j+idir)=(plts(2,j+idir)-plts(1,j+idir))/self%skyl%n%nexts(j+idir)
  end do

  deallocate(inwallr,inwallz,ouwallr,ouwallz)

end subroutine geoq_skylext
!---------------------------------------------------------------------
!> assist set up of skylight based on flux values
subroutine geoq_skylpsi1(self,ctrackrz,nctrack,kndbox,plt,pdelta,ktyps,kcall)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  real(kr8), dimension(:,:), allocatable, intent(in) :: ctrackrz !< central track
  integer(ki4), intent(in) :: nctrack !<  size of track array
  integer(ki4), dimension(2), intent(in) :: kndbox !<  size of box arrays
  real(kr8), dimension(2,2), intent(in) :: plt !<  limits for bin arrays
  real(kr8), dimension(2), intent(in) :: pdelta !<  delta for bin arrays
  integer(ki4), intent(in) :: ktyps !<  lower (1) or upper (2) skylight type
  integer(ki4), intent(in) :: kcall !<  number of call

  !! local
  character(*), parameter :: s_name='geoq_skylpsi1' !< subroutine name
  type(posang_t) :: posang !< position and vector involving angle
  integer(ki4) :: irid !< ! direction of travel in \f$ Z \f$ (+1) in lower, (-1) in upper
  integer(ki4) :: inou     !< whether on HFS (in, 1) or LFS (ou, 2) of centre
  integer(ki4) :: iupper   !< local variable
  integer(ki4) :: itp   !< local variable
  real(kr8) :: zpsi    !<  \f$ \psi \f$
  real(kr8), dimension(2) :: zrz    !<  Query point position \f$ (R,Z) \f$
  real(kr8) :: zr    !<  Query point \f$ R \f$
  real(kr8) :: zz    !<  Query point \f$ Z \f$
  real(kr8) :: zeta    !<  \f$ \zeta \f$
  logical :: ilnwin !<  within skylight window
  real(kr8) :: distrack     !< signed distance from central track

  ! direction of travel in z (negative of idir used elsewhere)
  irid=-(2*ktyps-3)

  ! complete initialisation
  self%skyl%inboxz(:,kcall)=(2-ktyps)*self%skyl%n%zextmin+(ktyps-1)*self%skyl%n%zextmax
  self%skyl%ouboxz(:,kcall)=(2-ktyps)*self%skyl%n%zextmin+(ktyps-1)*self%skyl%n%zextmax

  self%skyl%dimbox(:,kcall)=kndbox
  self%skyl%psilts(:,:,kcall)=plt
  self%skyl%psidelta(:,kcall)=pdelta

  do i=1,self%objl%np
     ! transform positions to R-Z-zeta space
     posang%pos=self%objl%posl%pos(i)%posvec
     posang%opt=0 ; posang%units=-3
     call posang_invtfm(posang,0)
     zr=posang%pos(1)
     zz=posang%pos(2)
     !
     ! test within skylight window
     ilnwin= ( (zr-self%skyl%n%rextmin)*(zr-self%skyl%n%rextmax)<0 .AND. &
 &   (zz-self%skyl%n%zextmin)*(zz-self%skyl%n%zextmax)<0 )
     if (.NOT.ilnwin) cycle
     if (irid*(zz-self%beq%n%zcen)>0) cycle
     call spl2d_evaln(self%beq%psi,zr,zz,1,zpsi)
     zrz=(/zr,zz/)
     call beq_ctrackcq(self%beq,ctrackrz,nctrack,zrz,distrack)
     inou=2+sign(0.5_kr8,distrack) ! 1 if distrack<0, 2 if distrack>0
     ! check within allowable range of psi
     if ( (zpsi-plt(1,inou))*(zpsi-plt(2,inou)) > 0 ) cycle
     itp=min( 1+int((zpsi-plt(1,inou))/pdelta(inou)), kndbox(inou) )
     if (inou==1) then
        if ( (zz-self%skyl%inboxz(itp,kcall))*irid>0 ) then
           self%skyl%inboxr(itp,kcall)=zr
           self%skyl%inboxz(itp,kcall)=zz
        end if
     else
        if ( (zz-self%skyl%ouboxz(itp,kcall))*irid>0 ) then
           self%skyl%ouboxr(itp,kcall)=zr
           self%skyl%ouboxz(itp,kcall)=zz
        end if
     end if
  end do

  call skyl_fixup1(self%skyl,ktyps,kcall,self%skyl%n%control)
  call skyl_fixup2(self%skyl,ktyps,kcall,self%skyl%n%control)

  if (self%beq%n%skyldbg>0) then
     call gfile_rwrite(self%skyl%inboxr(1,kcall),self%skyl%inboxz(1,kcall),kndbox(1),&
 &   'Preliminary inner skylight',self%skyl%ndskyl(4))
     call gfile_rwrite(self%skyl%ouboxr(1,kcall),self%skyl%ouboxz(1,kcall),kndbox(2),&
 &   'Preliminary outer skylight',self%skyl%ndskyl(5))

     allocate(work1(kndbox(1)),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     do j=1,kndbox(1)
        work1(j)=self%skyl%psilts(1,1,kcall)+(j-1)*self%skyl%psidelta(1,kcall)
     end do
     call gfile_rwrite(work1(1),self%skyl%inboxz(1,kcall),kndbox(1),&
 &   'Inner skylight Z vs flux',self%skyl%ndskyl(6))
     allocate(work1a(kndbox(2)),stat=status)
     call log_alloc_check(m_name,s_name,2,status)
     do j=1,kndbox(2)
        work1a(j)=self%skyl%psilts(1,2,kcall)+(j-1)*self%skyl%psidelta(2,kcall)
     end do
     call gfile_rwrite(work1a(1),self%skyl%ouboxz(1,kcall),kndbox(2),&
 &   'Outer skylight Z vs flux',self%skyl%ndskyl(7))
     deallocate(work1,work1a)
  end if

end subroutine geoq_skylpsi1
!---------------------------------------------------------------------
!> assist set up of skylight based on flux values
subroutine geoq_skylpsi2(self,ctrackrz,nctrack,kndbox,plt,pdelta,ktyps,kcall)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  real(kr8), dimension(:,:), allocatable, intent(in) :: ctrackrz !< central track
  integer(ki4), intent(in) :: nctrack !<  size of track array
  integer(ki4), dimension(2), intent(in) :: kndbox !<  size of box arrays
  real(kr8), dimension(2,2), intent(in) :: plt !<  limits for bin arrays
  real(kr8), dimension(2), intent(in) :: pdelta !<  delta for bin arrays
  integer(ki4), intent(in) :: ktyps !<  lower (1) or upper (2) skylight type
  integer(ki4), intent(in) :: kcall !<  number of call

  !! local
  character(*), parameter :: s_name='geoq_skylpsi2' !< subroutine name
  type(posang_t) :: posang !< position and vector involving angle
  integer(ki4) :: irid !< ! direction of travel in \f$ Z \f$ (+1) in lower, (-1) in upper
  integer(ki4) :: inou     !< whether on HFS (in, 1) or LFS (ou, 2) of centre
  integer(ki4) :: iupper   !< local variable
  integer(ki4) :: itp   !< local variable
  real(kr4), dimension(3) :: zbary !< vector of barycentre of geobj
  type(geobj_t) :: iobj !< object
  real(kr8) :: zpsi    !<  \f$ \psi \f$
  real(kr8), dimension(2) :: zrz    !<  Query point position \f$ (R,Z) \f$
  real(kr8) :: zr    !<  Query point \f$ R \f$
  real(kr8) :: zz    !<  Query point \f$ Z \f$
  real(kr8) :: zeta    !<  \f$ \zeta \f$
  logical :: ilnwin !<  within skylight window
  real(kr8) :: distrack     !< signed distance from central track

  ! direction of travel in z (negative of idir used elsewhere)
  irid=-(2*ktyps-3)

  ! complete initialisation
  self%skyl%inboxz(:,kcall)=(2-ktyps)*self%skyl%n%zextmin+(ktyps-1)*self%skyl%n%zextmax
  self%skyl%ouboxz(:,kcall)=(2-ktyps)*self%skyl%n%zextmin+(ktyps-1)*self%skyl%n%zextmax

  self%skyl%dimbox(:,kcall)=kndbox
  self%skyl%psilts(:,:,kcall)=plt
  self%skyl%psidelta(:,kcall)=pdelta

  do i=1,self%objl%ng
     ! calculate barycentre
     iobj%geobj=self%objl%obj2(i)%ptr
     iobj%objtyp=self%objl%obj2(i)%typ
     call geobj_centre(iobj,self%objl%posl,self%objl%nodl,zbary)
     ! transform positions to R-Z-zeta space
     posang%pos=zbary
     posang%opt=0 ; posang%units=-3
     call posang_invtfm(posang,0)
     zr=posang%pos(1)
     zz=posang%pos(2)
     !
     ! test within skylight window
     ilnwin= ( (zr-self%skyl%n%rextmin)*(zr-self%skyl%n%rextmax)<0 .AND. &
 &   (zz-self%skyl%n%zextmin)*(zz-self%skyl%n%zextmax)<0 )
     if (.NOT.ilnwin) cycle
     if (irid*(zz-self%beq%n%zcen)>0) cycle
     call spl2d_evaln(self%beq%psi,zr,zz,1,zpsi)
     zrz=(/zr,zz/)
     call beq_ctrackcq(self%beq,ctrackrz,nctrack,zrz,distrack)
     inou=2+sign(0.5_kr8,distrack) ! 1 if distrack<0, 2 if distrack>0
     ! check within allowable range of psi
     if ( (zpsi-plt(1,inou))*(zpsi-plt(2,inou)) > 0 ) cycle
     itp=min( 1+int((zpsi-plt(1,inou))/pdelta(inou)), kndbox(inou) )
     if (inou==1) then
        if ( (zz-self%skyl%inboxz(itp,kcall))*irid>0 ) then
           self%skyl%inboxr(itp,kcall)=zr
           self%skyl%inboxz(itp,kcall)=zz
        end if
     else
        if ( (zz-self%skyl%ouboxz(itp,kcall))*irid>0 ) then
           self%skyl%ouboxr(itp,kcall)=zr
           self%skyl%ouboxz(itp,kcall)=zz
        end if
     end if
  end do

  call skyl_fixup1(self%skyl,ktyps,kcall,self%skyl%n%control)
  call skyl_fixup2(self%skyl,ktyps,kcall,self%skyl%n%control)

  if (self%beq%n%skyldbg>0) then
     call gfile_rwrite(self%skyl%inboxr(1,kcall),self%skyl%inboxz(1,kcall),kndbox(1),&
 &   'Preliminary inner skylight',self%skyl%ndskyl(4))
     call gfile_rwrite(self%skyl%ouboxr(1,kcall),self%skyl%ouboxz(1,kcall),kndbox(2),&
 &   'Preliminary outer skylight',self%skyl%ndskyl(5))

     allocate(work1(kndbox(1)),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     do j=1,kndbox(1)
        work1(j)=self%skyl%psilts(1,1,kcall)+(j-1)*self%skyl%psidelta(1,kcall)
     end do
     call gfile_rwrite(work1(1),self%skyl%inboxz(1,kcall),kndbox(1),&
 &   'Inner skylight Z vs flux',self%skyl%ndskyl(6))
     allocate(work1a(kndbox(2)),stat=status)
     call log_alloc_check(m_name,s_name,2,status)
     do j=1,kndbox(2)
        work1a(j)=self%skyl%psilts(1,2,kcall)+(j-1)*self%skyl%psidelta(2,kcall)
     end do
     call gfile_rwrite(work1a(1),self%skyl%ouboxz(1,kcall),kndbox(2),&
 &   'Outer skylight Z vs flux',self%skyl%ndskyl(7))
     deallocate(work1,work1a)
  end if

end subroutine geoq_skylpsi2
!---------------------------------------------------------------------
!> calculate limits of limiter object
subroutine geoq_psilimiter(self)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  !! local
  character(*), parameter :: s_name='geoq_psilimiter' !< subroutine name
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zrmin    !<   \f$ R \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zrmax    !<   \f$ R \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  real(kr8) :: zzmin    !<   \f$ Z \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zzmax    !<   \f$ Z \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zpsi    !<  \f$ \psi(R,Z) \f$
  real(kr8) :: zf    !<  \f$ \f(\psi) \f$
  type(posang_t) :: posang !< position and vector involving anglea
  real(kr8) :: zpsimin    !<  least value of \f$ \psi \f$ on geometry
  real(kr8) :: zpsimax    !<  largest value of \f$ \psi \f$ on geometry
  real(kr8) :: ztheta    !<  \f$ \theta(R,Z) \f$
  real(kr8) :: zthetamin    !<  least value of \f$ \theta \f$ on geometry
  real(kr8) :: zthetamax    !<  largest value of \f$ \theta \f$ on geometry
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zsr    !< \f$  r_i \f$
  real(kr8) :: cylj    !<  current estimated in cylindrical approx.
  real(kr8) :: zbr    !<  radial field component
  real(kr8) :: zbz    !<  vertical field component
  real(kr8) :: zbt    !<  toroidal field component

  zpsimin=1.E+8
  zpsimax=-1.E+8
  zthetamin=1.E+8
  zthetamax=-1.E+8
  do j=1,self%objl%np
     ! transform positions to R-Z-zeta space
     posang%pos=self%objl%posl%pos(j)%posvec
     posang%opt=0 ; posang%units=-3
     call posang_invtfm(posang,0)
     zr=posang%pos(1)
     zz=posang%pos(2)
     if (self%beq%n%search==1) then
        ! must lie within search box
        if (zr<self%beq%n%lkrsta.OR.zr>self%beq%n%lkrend.OR. &
        zz<self%beq%n%lkzsta.OR.zz>self%beq%n%lkzend) cycle
     end if
     ! then to flux coordinates
     !! evaluate psi
     !      call spl2d_eval(self%beq%psi,zr,zz,zpsi)
     !   zpsimin=min(zpsi,zpsimin)
     call posang_psitfm(posang,self%beq)
     zpsi=posang%pos(1)
     ztheta=posang%pos(2)
     if (zpsi<zpsimin) then
        zpsimin=zpsi
        zrmin=zr ; zzmin=zz
     end if
     !   zpsimax=max(zpsi,zpsimax)
     if (zpsi>zpsimax) then
        zpsimax=zpsi
        zrmax=zr ; zzmax=zz
     end if
     zthetamin=min(zthetamin,ztheta)
     zthetamax=max(zthetamax,ztheta)
  end do

  if (beq_rsig()>0) then
     ! equiv. (self%psiaxis<self%psiqbdry)
     ! psi increasing outward
     self%beq%psiltr=zpsimin
     zr=zrmin ; zz=zzmin
     self%beq%psiotr=zpsimax
  else
     self%beq%psiltr=zpsimax
     zr=zrmax ; zz=zzmax
     self%beq%psiotr=zpsimin
  end if

  self%beq%thetagmax=zthetamax
  self%beq%thetagmin=zthetamin

  call spl2d_evaln(self%beq%dpsidr,zr,zz,1,zdpdr)
  call spl2d_evaln(self%beq%dpsidz,zr,zz,2,zdpdz)
  call spl2d_evaln(self%beq%psi,zr,zz,2,zpsi)
  self%beq%rbdry=zr
  self%beq%bpbdry=(1/zr)*sqrt( max(0.,(zdpdr**2+zdpdz**2)) )

  ! evaluate I aka f at psi
  call spleval(self%beq%f,self%beq%mr,self%beq%psiaxis,self%beq%psiqbdry,zpsi,zf,1)
  self%beq%btotbdry=sqrt( max(0.,(self%beq%bpbdry**2+(zf/zr)**2)) )

  call log_error(m_name,s_name,1,log_info,'Reference boundary values')
  call log_value("SMITER-GEOQ psiltr ",self%beq%psiltr)
  call log_value("SMITER-GEOQ rbdry ",self%beq%rbdry)
  call log_value("SMITER-GEOQ zbdry ",zz)
  call log_value("SMITER-GEOQ bpbdry ",self%beq%bpbdry)
  call log_value("SMITER-GEOQ btotbdry ",self%beq%btotbdry)
  zsr=sqrt( max(0.,(zr-self%beq%n%rcen)**2+(zz-self%beq%n%zcen)**2) )
  cylj=abs(zsr*self%beq%bpbdry/2.e-7)
  call log_value("Estimated cylindrical current ",cylj)

  zbr=-(1/zr)*zdpdz
  zbz=(1/zr)*zdpdr
  zbt=zf/zr
  call log_value("radial field component",zbr)
  call log_value("vertical field component",zbz)
  call log_value("toroidal field component",zbt)

  call log_error(m_name,s_name,2,log_info,'Object limits ')
  call log_value("OBJECT LIMIT psimax ",zpsimax)
  call log_value("OBJECT LIMIT psimin ",zpsimin)
  call log_value("OBJECT LIMIT thetamax ",zthetamax)
  call log_value("OBJECT LIMIT thetamin ",zthetamin)

end subroutine geoq_psilimiter
!---------------------------------------------------------------------
!> calculate psi of silhouette object
subroutine geoq_psisilh(self)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  !! local
  character(*), parameter :: s_name='geoq_psisilh' !< subroutine name
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zrmin    !<   \f$ R \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zrmax    !<   \f$ R \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  real(kr8) :: zzmin    !<   \f$ Z \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zzmax    !<   \f$ Z \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zpsi    !<  \f$ \psi(R,Z) \f$
  real(kr8) :: zf    !<  \f$ \f(\psi) \f$
  type(posang_t) :: posang !< position and vector involving anglea
  real(kr8) :: zpsimin    !<  least value of \f$ \psi \f$ on geometry
  real(kr8) :: zpsimax    !<  largest value of \f$ \psi \f$ on geometry
  real(kr8), dimension(3)  :: zv1 !<   3-D Cartesian vector
  real(kr8), dimension(3)  :: zvd !<   3-D Cartesian vector difference
  real(kr8) :: zline    !<  length of line
  integer(ki4) :: inline    !<  length of line in segments
  integer(ki4) :: iobj   !< object pointer
  integer(ki4) :: inode   !< point (node) pointer
  integer(ki4) :: ityp   !< object type
  integer(ki4) :: inumpts   !< number of points in object
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zsr    !< \f$  r_i \f$
  real(kr8) :: cylj    !<  current estimated in cylindrical approx.
  real(kr8) :: zbr    !<  radial field component
  real(kr8) :: zbz    !<  vertical field component
  real(kr8) :: zbt    !<  toroidal field component

  zpsimin=1.E+8
  zpsimax=-1.E+8
  loop_object: do j=1,self%objl%ng
     iobj=self%objl%obj2(j)%ptr
     ityp=self%objl%obj2(j)%typ
     inumpts=geobj_entry_table_fn(ityp)
     do k=1,inumpts-1
        ! length of line
        inode=self%objl%nodl(iobj)
        zv1=self%objl%posl%pos(inode)%posvec
        inode=self%objl%nodl(iobj+1)
        zvd=self%objl%posl%pos(inode)%posvec-zv1
        zline=sqrt( max(0.,zvd(1)**2+zvd(2)**2+zvd(3)**2) )
        inline=zline/deltal
        !dbg   if (j<=10) write(*,*) inode,zvd,inline   !dbg
        loop_segment: do l=1,max(inline,1)
           ! transform positions to R-Z-zeta space
           posang%pos=zv1+(l-1)*(zvd/inline)
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           zr=posang%pos(1)
           zz=posang%pos(2)
           if (self%beq%n%search==1) then
              ! must lie within search box
              if (zr<self%beq%n%lkrsta.OR.zr>self%beq%n%lkrend.OR. &
              zz<self%beq%n%lkzsta.OR.zz>self%beq%n%lkzend) cycle
           end if
           ! evaluate psi
           call spl2d_eval(self%beq%psi,zr,zz,zpsi)
           !dbg   if (j<=10) write(*,*) zr,zz,zpsi  !dbg
           !   zpsimin=min(zpsi,zpsimin)
           if (zpsi<zpsimin) then
              zpsimin=zpsi
              zrmin=zr ; zzmin=zz
           end if
           !   zpsimax=max(zpsi,zpsimax)
           if (zpsi>zpsimax) then
              zpsimax=zpsi
              zrmax=zr ; zzmax=zz
           end if
        end do loop_segment
     end do
  end do loop_object

  !dbg   write(*,*) zrmin,zrmax,zpsimin,zpsimax !dbg
  if (beq_rsig()>0) then
     ! equiv. (self%psiaxis<self%psiqbdry)
     ! psi increasing outward
     self%beq%psiltr=zpsimin
     zr=zrmin ; zz=zzmin
     self%beq%psiotr=zpsimax
  else
     self%beq%psiltr=zpsimax
     zr=zrmax ; zz=zzmax
     self%beq%psiotr=zpsimin
  end if
  !dbg write(*,*) 'psiltr,psiotr', self%beq%psiltr,self%beq%psiotr !dbg
  call spl2d_evaln(self%beq%dpsidr,zr,zz,1,zdpdr)
  call spl2d_evaln(self%beq%dpsidz,zr,zz,2,zdpdz)
  self%beq%rbdry=zr
  self%beq%bpbdry=(1/zr)*sqrt( max(0.,(zdpdr**2+zdpdz**2)) )

  ! evaluate I aka f at psi
  zpsi=self%beq%psiotr
  call spleval(self%beq%f,self%beq%mr,self%beq%psiaxis,self%beq%psiqbdry,zpsi,zf,1)
  self%beq%btotbdry=sqrt( max(0.,(self%beq%bpbdry**2+(zf/zr)**2)) )
  !dbg   write(*,*) self%beq%psiltr,self%beq%psiotr !dbg
  call log_error(m_name,s_name,1,log_info,'Reference boundary values')
  call log_value("SMITER-GEOQ psiltr ",self%beq%psiltr)
  call log_value("SMITER-GEOQ rbdry ",self%beq%rbdry)
  call log_value("SMITER-GEOQ zbdry ",zz)
  call log_value("SMITER-GEOQ bpbdry ",self%beq%bpbdry)
  call log_value("SMITER-GEOQ btotbdry ",self%beq%btotbdry)
  zsr=sqrt( max(0.,(zr-self%beq%n%rcen)**2+(zz-self%beq%n%zcen)**2) )
  cylj=abs(zsr*self%beq%bpbdry/2.e-7)
  call log_value("Estimated cylindrical current ",cylj)

  zbr=-(1/zr)*zdpdz
  zbz=(1/zr)*zdpdr
  zbt=zf/zr
  call log_value("radial field component",zbr)
  call log_value("vertical field component",zbz)
  call log_value("toroidal field component",zbt)

end subroutine geoq_psisilh
!---------------------------------------------------------------------
!> calculate distance between silhouette and flux contour as function of silhouette position
subroutine geoq_dsilhcont(self,prc,pzc,prs,pzs,knear,pdist)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  real(kr8), dimension(:), intent(in) :: prc !< \f$ R \f$ values of contour
  real(kr8), dimension(:), intent(in) :: pzc !< \f$ Z \f$ values of contour
  real(kr8), dimension(:), allocatable, intent(out) :: prs !< \f$ R \f$ values of silhouette
  real(kr8), dimension(:), allocatable, intent(out) :: pzs !< \f$ Z \f$ values of silhouette
  integer(ki4), dimension(:), allocatable, intent(out) :: knear !< nearest segment of silhouette
  real(kr8), dimension(:), allocatable, intent(out) :: pdist !< distances of contour from silhouette

  !! local
  character(*), parameter :: s_name='geoq_dsilhcont' !< subroutine name
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  type(posang_t) :: posang !< position and vector involving anglea
  real(kr8) :: zdist    !<  initial least value of distance
  real(kr8) :: zdistsq    !<  initial least value of distance squared
  real(kr8) :: zdsq    !<  distance squared
  real(kr8), dimension(3)  :: zv1 !<   3-D Cartesian vector
  real(kr8), dimension(3)  :: zvd !<   3-D Cartesian vector difference
  real(kr8) :: zline    !<  length of line
  integer(ki4) :: inc    !<  length of contour in segments
  integer(ki4) :: inseg    !<  length of silhouette in segments
  integer(ki4) :: inline    !<  length of line in segments
  integer(ki4) :: il    !<  loope counter for line in segments
  integer(ki4) :: iobj   !< object pointer
  integer(ki4) :: inode   !< point (node) pointer
  integer(ki4) :: ityp   !< object type
  integer(ki4) :: inumpts   !< number of points in object

  inc=size(prc)
  ! work out size of distance array
  inseg=0
  loop_size: do j=1,self%objl%ng
     iobj=self%objl%obj2(j)%ptr
     ityp=self%objl%obj2(j)%typ
     inumpts=geobj_entry_table_fn(ityp)
     do k=1,inumpts-1
        ! length of line
        inode=self%objl%nodl(iobj)
        zv1=self%objl%posl%pos(inode)%posvec
        inode=self%objl%nodl(iobj+1)
        zvd=self%objl%posl%pos(inode)%posvec-zv1
        zline=sqrt( max(0.,zvd(1)**2+zvd(2)**2+zvd(3)**2) )
        inline=zline/deltal
        inseg=inseg+inline
     end do
  end do loop_size

  !! allocate distance/silhouette position storage
  if(inseg>0) then
     allocate(prs(inseg),pzs(inseg),knear(inseg),pdist(inseg),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No silhouette position data')
  end if

  ! calculate positions for distance array
  il=0
  loop_object: do j=1,self%objl%ng
     iobj=self%objl%obj2(j)%ptr
     ityp=self%objl%obj2(j)%typ
     inumpts=geobj_entry_table_fn(ityp)
     do k=1,inumpts-1
        ! length of line
        inode=self%objl%nodl(iobj)
        zv1=self%objl%posl%pos(inode)%posvec
        inode=self%objl%nodl(iobj+1)
        zvd=self%objl%posl%pos(inode)%posvec-zv1
        zline=sqrt( max(0.,zvd(1)**2+zvd(2)**2+zvd(3)**2) )
        inline=zline/deltal
        loop_segment: do l=1,inline
           ! transform positions to R-Z-zeta space
           posang%pos=zv1+(l-1)*(zvd/inline)
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           il=il+1
           ! save
           prs(il)=posang%pos(1)
           pzs(il)=posang%pos(2)
        end do loop_segment
     end do
  end do loop_object

  ! now calculate distance
  loop_segments: do j=1,inseg
     zr=prs(j)
     zz=pzs(j)
     ! evaluate distance to contour
     zdist=1.E+8
     zdistsq=zdist**2
     do k=1,inc
        ! squared distance between points
        zdsq=(zr-prc(k))**2+(zz-pzc(k))**2
        if (zdsq<zdistsq) then
           zdistsq=zdsq
           knear(j)=k
        end if
     end do
     pdist(j)=sqrt( max(0._kr8,zdistsq) )
  end do loop_segments

  call log_error(m_name,s_name,1,log_info,'Nearest contour segments calculated')

end subroutine geoq_dsilhcont
!---------------------------------------------------------------------
!> write (vtk)  geoq data structure
subroutine geoq_writev(self,kchar,kplot)

  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  character(*), intent(in) :: kchar  !< case
  integer :: kplot   !< output channel for vis. data


  !! local
  character(*), parameter :: s_name='geoq_writev' !< subroutine name
  type(posang_t) :: posang !< position and vector involving angles
  integer(ki4) :: isum !< sum workspace
  integer(ki2) :: inn !< number of nodes defining geobj
  integer(ki4) :: iplall !< flag plot 'all' option
  !  integer(ki4) :: irzzeta !< flag plot 'all' option
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4), dimension(3) :: zbary !< vector of barycentre of geobj
  real(kr4) :: zmag !< magnitude of normal vector
  type(geobj_t) :: iobj !< object
  logical, parameter :: output_cartv=.TRUE. !< local variable
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zpostfm !< local variable
  type(tfmdata_t) :: ztfmdata !< local variable
  integer(ki4) :: ityp   !< object type

  ! preamble
  ibuf1=kchar
  !  irzzeta=0
  iplall=0
  plot_init_type: select case (kchar)
  case('geometry')
  case('all')
     iplall=1
  case('frzzeta')
     !     irzzeta=1
     iplall=1
  case('frzxi')
     !     irzzeta=2
     iplall=1
  case('fxyz')
     !     irzzeta=3
     iplall=1
  case('allcart')
  end select plot_init_type

  plot_type: select case (kchar)
  case('geometry','all','frzzeta','frzxi','fxyz')
     write(kplot,'(''DATASET UNSTRUCTURED_GRID'')')
     write(kplot,'(''POINTS '',I8, '' float'')') self%objl%np
     ! transform positions
     geom_type: select case (ibuf1)
     case('geometry','all')
        do j=1,self%objl%np
           ! transform positions to psi-theta-zeta space and write
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           call posang_psitfm(posang,self%beq)
           call posang_writev(posang,kplot)
        end do
     case('frzzeta')
        do j=1,self%objl%np
           ! transform positions to R-Z-zeta space and write
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           call posang_writev(posang,kplot)
        end do
     case('frzxi')
        ! default identity
        ztfmdata%ntfm=2
        ztfmdata%matrix(1,:)=(/1.,  0.,  0. /)
        ztfmdata%matrix(2,:)=(/0.,  1.,  0. /)
        ztfmdata%matrix(3,:)=(/0.,  0.,  1. /)
        ztfmdata%scale=(/1.,  1.,  1. /)
        ztfmdata%offset=(/0.,  0.,  0. /)
        ztfmdata%scale(3)=real(self%beq%nzets)

        do j=1,self%objl%np
           ! transform positions to R-Z-xi space and write
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           zpos%posvec=posang%pos
           zpostfm=position_tfm(zpos,ztfmdata)
           posang%pos=zpostfm%posvec
           call posang_writev(posang,kplot)
        end do
     case('fxyz')
        ztfmdata=self%beq%fmesh%tfmdata

        do j=1,self%objl%np
           ! transform positions from CAD (geometry in mm) coordinates to duct  (m)
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0 ; posang%units=-3
           call posang_invcartfm(posang,ztfmdata,0)
           call posang_writev(posang,kplot)
        end do
     end select geom_type
     write(kplot, '('' '')')

     ! count entries in CELLS
     isum=0
     do j=1,self%objl%ng
        inn=geobj_entry_table_fn(self%objl%obj2(j)%typ)
        isum=isum+inn
     end do

     ! output CELL data
     write(kplot,'(''CELLS '',I8,1X,I8)') self%objl%ng,self%objl%ng+isum
     i=1
     do j=1,self%objl%ng
        inn=geobj_entry_table_fn(self%objl%obj2(j)%typ)
        write(kplot,'(8(1X,I8))') inn,(self%objl%nodl(i+ij-1)-1,ij=1,inn)
        i=i+inn
     end do
     write(kplot, '('' '')')

     ! output CELL types
     write(kplot,'(''CELL_TYPES '',I8)') self%objl%ng
     do j=1,self%objl%ng
        ityp=self%objl%obj2(j)%typ
        write(kplot,'(1X,I8)') geobj_type_fn(ityp)
     end do
     write(kplot, '('' '')')

     if (self%objl%nparam(2)==1) then
        write(kplot,'(''CELL_DATA '',I8)') self%objl%ng
        ! output geometry codes
        write(kplot,'(''SCALARS Code int'')')
        write(kplot,'(''LOOKUP_TABLE default'')')
        do j=1,self%objl%ng
           ityp=self%objl%obj2(j)%typ
           write(kplot,'(1X,I8)') geobj_code_fn(ityp)
        end do

        write(kplot, '('' '')')
     end if

     write(kplot,'(''POINT_DATA '',I8)') self%objl%np

     ! output cartesian vectors (mm) usually, no transform

     if (output_cartv) then
        write(kplot,'(''VECTORS Xcart float'')')
        do j=1,self%objl%np
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0
           call posang_writev(posang,kplot)
        end do
     end if

     if(iplall==1)then

        geom_type1:select case (ibuf1)
        case('fxyz')
           write(kplot, '('' '')')
           ! output cartesian B vectors in duct coordinates
           write(kplot,'(''VECTORS Bcart float'')')
           ztfmdata=self%beq%fmesh%tfmdata
           do j=1,self%objl%np
              ! B evaluated at position in CAD coordinates
              posang%pos=self%objl%posl%pos(j)%posvec
              posang%opt=0 ; posang%units=-3
              call beq_b(self%beq,posang,0)
              ! now convert B to duct components (and coordinates)
              call posang_invcartfm(posang,ztfmdata,0)
              call posang_writev(posang,kplot,2)
           end do
        case default
           write(kplot, '('' '')')
           ! output cartesian B vectors
           write(kplot,'(''VECTORS Bcart float'')')
           do j=1,self%objl%np
              ! position
              posang%pos=self%objl%posl%pos(j)%posvec
              posang%opt=0 ; posang%units=-3
              !BB! first position to cyl polars to get psi derivs
              !BB         call posang_invtfm(posang,0)
              !BB         zr=posang%pos(1)
              !BB         zz=posang%pos(2)
              !BB         call spl2d_eval(self%beq%dpsidr,zr,zz,zdpdr)
              !BB         call spl2d_eval(self%beq%dpsidz,zr,zz,zdpdz)
              !BB         call spl2d_eval(self%beq%psi,zr,zz,zpsi)
              !BB! evaluate I aka f at psi
              !BB         call spleval(self%beq%f,self%beq%mr,self%beq%psiaxis,self%beq%psiqbdry,zpsi,zf,1)
              !BB! B in toroidal-cyl polars
              !BB         posang%vec(1)=-zdpdz/zr ; posang%vec(2)=zdpdr/zr ; posang%vec(3)=zf/zr
              !BB         posang%opt=17
              !BB! vec to cartesians and output
              !BB         call posang_tfm(posang,-3)
              call beq_b(self%beq,posang,0)
              call posang_writev(posang,kplot,2)
           end do
        end select geom_type1

     end if

  case('allcart')
     write(kplot,'(''DATASET UNSTRUCTURED_GRID'')')
     write(kplot,'(''POINTS '',I8, '' float'')') self%objl%np
     do j=1,self%objl%np
        ! write positions in Cartesians (mm)
        posang%pos=self%objl%posl%pos(j)%posvec
        posang%opt=0 ; posang%units=-3
        call posang_writev(posang,kplot)
     end do
     write(kplot, '('' '')')

     ! count entries in CELLS
     isum=0
     do j=1,self%objl%ng
        inn=geobj_entry_table_fn(self%objl%obj2(j)%typ)
        isum=isum+inn
     end do

     ! output CELL data
     write(kplot,'(''CELLS '',I8,1X,I8)') self%objl%ng,self%objl%ng+isum
     i=1
     do j=1,self%objl%ng
        inn=geobj_entry_table_fn(self%objl%obj2(j)%typ)
        write(kplot,'(8(1X,I8))') inn,(self%objl%nodl(i+ij-1)-1,ij=1,inn)
        i=i+inn
     end do
     write(kplot, '('' '')')

     ! output CELL types B

     write(kplot,'(''CELL_TYPES '',I8)') self%objl%ng
     do j=1,self%objl%ng
        ityp=self%objl%obj2(j)%typ
        write(kplot,'(1X,I8)') geobj_type_fn(ityp)
     end do

     write(kplot, '('' '')')

     write(kplot,'(''CELL_DATA '',I8)') self%objl%ng
     if (self%objl%nparam(2)==1) then
        ! output geometry codes
        write(kplot,'(''SCALARS Code int'')')
        write(kplot,'(''LOOKUP_TABLE default'')')
        do j=1,self%objl%ng
           ityp=self%objl%obj2(j)%typ
           write(kplot,'(1X,I8)') geobj_code_fn(ityp)
        end do
     end if
     ! output normal vectors
     write(kplot,'(''VECTORS Normal float'')')
     do j=1,self%objl%ng
        iobj%geobj=self%objl%obj2(j)%ptr
        iobj%objtyp=self%objl%obj2(j)%typ
        call geobj_normal(iobj,self%objl%posl,self%objl%nodl,znormal,zmag)
        write(kplot, cfmtb1v ) znormal
     end do
     write(kplot, '('' '')')
     ! output cartesian B vectors as cell data
     write(kplot,'(''VECTORS Bcell float'')')
     do j=1,self%objl%ng
        ! calculate barycentre
        iobj%geobj=self%objl%obj2(j)%ptr
        iobj%objtyp=self%objl%obj2(j)%typ
        call geobj_centre(iobj,self%objl%posl,self%objl%nodl,zbary)
        ! evaluate field at barycentre
        posang%pos=zbary
        posang%opt=0 ; posang%units=-3
        call beq_b(self%beq,posang,0)
        call posang_writev(posang,kplot,2)
     end do
     ! output cartesian B vectors as cell data
     write(kplot,'(''VECTORS B-RZzeta float'')')
     do j=1,self%objl%ng
        ! calculate barycentre
        iobj%geobj=self%objl%obj2(j)%ptr
        iobj%objtyp=self%objl%obj2(j)%typ
        call geobj_centre(iobj,self%objl%posl,self%objl%nodl,zbary)
        ! evaluate field at barycentre
        posang%pos=zbary
        posang%opt=0 ; posang%units=-3
        call beq_b(self%beq,posang,0)
        ! convert to polar-toroidal
        posang%opt=16
        call posang_invtfm(posang,-3)
        call posang_writev(posang,kplot,2)
     end do

     write(kplot, '('' '')')

     ! output psi scalar
     write(kplot,'(''POINT_DATA '',I8)') self%objl%np
     write(kplot,'(''SCALARS Psi float'')')
     write(kplot,'(''LOOKUP_TABLE default'')')
     do j=1,self%objl%np
        ! transform positions to psi-theta-zeta space and write psi
        posang%pos=self%objl%posl%pos(j)%posvec
        posang%opt=0 ; posang%units=-3
        call posang_invtfm(posang,0)
        call posang_psitfm(posang,self%beq)
        write(kplot, cfmtbs ) posang%pos(1)
     end do

     write(kplot, '('' '')')
     ! output cartesian B vectors
     write(kplot,'(''VECTORS Bcart float'')')
     do j=1,self%objl%np
        ! field at nodes
        posang%pos=self%objl%posl%pos(j)%posvec
        posang%opt=0 ; posang%units=-3
        call beq_b(self%beq,posang,0)
        call posang_writev(posang,kplot,2)
     end do
  end select plot_type

end subroutine geoq_writev

!---------------------------------------------------------------------
!> write (gnu)  geoq data structure
subroutine geoq_writeg(self,kchar,kout,kopt)

  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  character(*), intent(in) :: kchar  !< case
  integer, intent(in) :: kout   !< output channel for gnuplot data
  !> which points plotted
  !! 0=original geometry, 1=geometry and extra, 2=extra geometry only
  integer(ki4), intent(in) :: kopt   !< .

  !! local
  character(*), parameter :: s_name='geoq_writeg' !< subroutine name
  type(posang_t) :: posang !< position and vector involving angles
  integer(ki4) :: is   !< start loop
  integer(ki4) :: ie   !< end loop

  if (kopt==0) then
     is=1
     ie=self%objl%posl%np
  else if (kopt==1) then
     is=1
     ie=self%objl%np
  else if (kopt==2) then
     is=self%objl%posl%np+1
     ie=self%objl%np
  end if

  plot_type: select case (kchar)
  case('gnusil')
     ! positions in R-Z-zeta space
     do j=is,ie
        posang%pos=self%objl%posl%pos(j)%posvec
        posang%opt=0 ; posang%units=-3
        call posang_invtfm(posang,0)
        write(kout,'(1x,i9,'//cfmt2v,iostat=status) &
 &      j,posang%pos(1),posang%pos(2),posang%pos(3)
        call log_write_check(m_name,s_name,1,status)
     end do
  case('gnusilm')
     ! positions in psi-theta-zeta space
     do j=is,ie
        posang%pos=self%objl%posl%pos(j)%posvec
        posang%opt=0 ; posang%units=-3
        call posang_invtfm(posang,0)
        call posang_psitfm(posang,self%beq)
        write(kout,'(1x,i9,'//cfmt2v,iostat=status) &
 &      j,posang%pos(1),posang%pos(2),posang%pos(3)
        call log_write_check(m_name,s_name,2,status)
     end do
  end select plot_type

end subroutine geoq_writeg

end module geoq_m
