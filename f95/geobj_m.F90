module geobj_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use position_h
  use position_m

  implicit none
  private

! public subroutines
  public :: geobj_readv, & !< read in visualisation format
  geobj_writev, & !< write in visualisation format
  geobj_writestl, & !< write object in stl format
  geobj_inbox, & !< test whether within real box
  geobj_innbox, &  !< test whether within normalised box
  geobj_linehits, &  !< test whether line hits
  geobj_trihitsbox, &  !< test whether triangle intersects real box
  geobj_normal, &  !< calculate normal to geobj
  geobj_area, &  !< calculate area of geobj
  geobj_nodes, &  !< return geobj nodal values
  geobj_centre, &  !< calculate barycentre of geobj
  geobj_sample, & !< calculate sample point of geobj
  geobj_entry_table_fn, &  !< return number of entries in geobj_entry_table
  geobj_type_fn, &  !< extract last 8 bits in geobj_entry_table
  geobj_code_fn  !< get geometry code


! public types
  type, public :: geobj1_t
   !! geobj1 data type ! type 1 - weighted particles
   ! virtual identifiers (just 1...n)
   !     integer(ki4)  :: lab !< geobj1 identifier (index of point)
     real(kr4) :: weight !< object weight
  end type geobj1_t

  type, public :: geobj2_t
   !! geobj2 data type ! type 2 - list of labels
   ! identifiers
     integer(ki4) :: ptr !< geobj2 identifier (pointer to label array)
   !> geobj2
   !! geometry type coded using vtk type (5 for triangles etc.) plus
   !! code describing function of each object (type+256*integer in list)
   !! 0. geometry
   !! 1. periodic boundary (high)
   !! 2. periodic boundary (low)
   !! 3. transparent geometry, totally ignored
   !! 4. skylight
   !! 5. terminating geometry, no termplane instructions
   !! 6. beancan
   !! 7. cutout
   !! 32+ptr. termination plane
   !! 64+ptr. boundary interface
     integer(ki2) :: typ !< -
  end type geobj2_t

  type, public :: geobj_t
     integer(ki4) :: geobj !< geobj index
     integer(ki2) :: objtyp !< geobj type
  end type geobj_t

!public variables
!> table from ParaView manual
  integer(ki4), dimension(27), parameter, public :: geobj_entry_table = & !< vtk variable
 &(/1,0,2,0,3,0,0,4,4,4,8,8,6,5,10,12,0,0,0,0,3,6,8,10,20,15,13/)
!     1 2 3 4 5 6 7 8 9 0 1 2 3 4 5  6  7 8 9 0 1 2 3 4  5  6  7
  integer(ki4), parameter, public :: geobj_max_entry_table = 20 !< vtk variable
!> define most of cell types likely to be needed
  integer(ki2par),  parameter, public :: VTK_VERTEX=1 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_POLY_VERTEX=2 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_LINE=3 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_POLY_LINE=4 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_TRIANGLE=5 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_TRIANGLE_STRIP=6 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_POLYGON=7 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_PIXEL=8 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_QUAD=9 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_TETRA=10 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_VOXEL=11 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_HEXAHEDRON=12 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_WEDGE=13 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_PYRAMID=14 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_QUADRATIC_EDGE=21 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_QUADRATIC_TRIANGLE=22 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_QUADRATIC_QUAD=23 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_QUADRATIC_TETRA=24 !< vtk variable
  integer(ki2par),  parameter, public :: VTK_QUADRATIC_HEXAHEDRON=25 !< vtk variable

  integer(ki2par),  parameter, public :: GEOBJ_BIT_SHIFT=8 !< for encoding geometry
  integer(ki2par),  parameter, public :: GEOBJ_POW=2**GEOBJ_BIT_SHIFT !< for encoding geometry
  integer(ki2par),  parameter, public :: GEOBJ_ABSORB=0!< First wall (absorbing boundary)
  integer(ki2par),  parameter, public :: GEOBJ_PERMAX=1!< Periodic boundary (Max third angle)
  integer(ki2par),  parameter, public :: GEOBJ_PERMIN=2!< Periodic boundary (Min third angle)
  integer(ki2par),  parameter, public :: GEOBJ_INVISI=3!< Transparent geometry, totally ignored
  integer(ki2par),  parameter, public :: GEOBJ_SKYLIT=4!< Skylight
  integer(ki2par),  parameter, public :: GEOBJ_ESCAPE=5!< Escape boundary (Loss)
  integer(ki2par),  parameter, public :: GEOBJ_ERRLOS=6!< Beancan (Error Loss)
  integer(ki2par),  parameter, public :: GEOBJ_CUTOUT=7!< Cutout (Intended Loss)
  integer(ki2par),  parameter, public :: GEOBJ_TERMPL=32!< Terminating plane, pointer to %termplane array
  integer(ki2par),  parameter, public :: GEOBJ_INTERF=64!< Boundary interface, pointer to boundary array TO BE DONE
!> Offset between geometry and handling codes
!!  Plan is to remove use of nobjhit codes -1 and -2 so this can be set to zero
  integer(ki2par),  parameter, public :: GEOBJ_OFFSET=-2!< -
! private types

! private variables
  character(*), parameter :: m_name='geobj_m' !< module name
  integer   :: status   !< error status
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki2) :: inn !< number of nodes


  contains
!---------------------------------------------------------------------
!> read in geobj data
subroutine geobj_readv(self,k,kin)

  !! arguments
  type(geobj_t), intent(out) :: self   !< geobj data
  integer(ki4), intent(in) :: k   !< part of geobj data
  integer, intent(in) :: kin   !< input channel for geobj data


  !! local
  character(*), parameter :: s_name='geobj_readv' !< subroutine name

  if (k==1) then
     read(kin,*,iostat=status) self%geobj
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Error reading geobj')
     end if
  else if (k==2) then
     read(kin,*,iostat=status) self%objtyp
     if(status/=0) then
        call log_error(m_name,s_name,2,error_fatal,'Error reading objtyp')
     end if
  end if

end subroutine geobj_readv
!---------------------------------------------------------------------
!> output geobj vector
subroutine geobj_writev(self,k,kout)

  !! arguments
  type(geobj_t), intent(in) :: self   !< geobj data
  integer(ki4), intent(in) :: k   !< part of geobj data
  integer, intent(in) :: kout   !< output channel for geobj data


  !! local
  character(*), parameter :: s_name='geobj_writev' !< subroutine name

  if (k==1) then
     write(kout,'(1x,i8)',iostat=status) self%geobj
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Error writing geobj index')
     end if
  else if (k==2) then
     write(kout,'(1x,i8)',iostat=status) self%objtyp
     if(status/=0) then
        call log_error(m_name,s_name,2,error_fatal,'Error writing objtyp')
     end if
  end if

end subroutine geobj_writev
!---------------------------------------------------------------------
!> output geobj vectors in stl format
subroutine geobj_writestl(self,posl,nodl,k,kout)

  !! arguments
  type(geobj_t), intent(in) :: self   !< geobj data
  type(posveclis_t), intent(in) :: posl   !< list of positions
  integer(ki4), dimension(*), intent(in) :: nodl   !< list of nodes
  integer(ki4), intent(in) :: k   !< part of geobj data (inert, but keep for consistency)
  integer, intent(in) :: kout   !< output channel for geobj data

  !! local
  character(*), parameter :: s_name='geobj_writestl' !< subroutine name
  integer(ki4) :: ii   !< geobj position
  integer(ki4) :: jj   !< loop counter
  real(kr4), dimension(3,8) :: xnodes !< x(compt,node) of obj
  integer(ki4), dimension(8) :: inod !< nodes of obj
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4) :: zmag !< magnitude of normal vector
  integer(ki2) :: ityp !< object type

  ityp=ibits(self%objtyp,0,GEOBJ_BIT_SHIFT)
  if (ityp==VTK_TRIANGLE) then
     ! triangle type
     inn=geobj_entry_table(ityp)
     ii=self%geobj
     do jj=1,inn
        inod(jj)=nodl(ii+jj-1)
        xnodes(:,jj)=posl%pos(inod(jj))%posvec
     end do
     call geobj_normal(self,posl,nodl,znormal,zmag)
     write(kout,'(''facet normal '''//cfmte1v) (znormal(i),i=1,3)
     write(kout,'(''    outer loop'')')
     write(kout,'(''        vertex '''//cfmte1v) (xnodes(i,1),i=1,3)
     write(kout,'(''        vertex '''//cfmte1v) (xnodes(i,2),i=1,3)
     write(kout,'(''        vertex '''//cfmte1v) (xnodes(i,3),i=1,3)
     write(kout,'(''    endloop'')')
     write(kout,'(''endfacet'')')

  else
     ! not triangle type
     call log_error(m_name,s_name,1,error_fatal,'Can only write triangle geobj')
  end if

end subroutine geobj_writestl
!---------------------------------------------------------------------
!> is geobj in real box
function geobj_inbox(self,posl,nodl,box)

  !! arguments
  type(geobj_t), intent(in) :: self   !< geobj
  type(posveclis_t), intent(in) :: posl   !< list of positions
  integer(ki4), dimension(*), intent(in) :: nodl   !< list of nodes
  real(kr4), dimension(3,2), intent(in) :: box !< bounding box corners

  !! local
  character(*), parameter :: s_name='geobj_inbox' !< subroutine name
  logical :: geobj_inbox !< function variable type
  type(posvecl_t) :: zpos   !< position
  logical :: linint   !< geobj in interval
  logical :: linbox   !< geobj in interval
  integer(ki4) :: ii   !< geobj position
  integer(ki4) :: jj   !< loop counter
  real(kr4), dimension(3,8) :: xnodes !< x(compt,node) of obj
  integer(ki4), dimension(8) :: inod !< nodes of obj
  real(kr4), dimension(3) :: xc !< centre of  box
  real(kr4), dimension(3) :: hh !< half side of box
  integer(ki2) :: ityp !< object type

  ityp=ibits(self%objtyp,0,GEOBJ_BIT_SHIFT)

  if (self%objtyp==1) then
     ! particle type
     ii=self%geobj
     zpos=posl%pos(ii)
     linbox=.true.
     do jj=1,3
        linint=( box(jj,1)<=zpos%posvec(jj) )
        linint=linint.and.( zpos%posvec(jj)<box(jj,2) )
        linbox=linbox.and.linint
     end do

  else

     ityp=ibits(self%objtyp,0,GEOBJ_BIT_SHIFT)
     if (ityp==VTK_TRIANGLE) then
        ! triangle type
        inn=geobj_entry_table(ityp)
        ii=self%geobj
        do jj=1,inn
           inod(jj)=nodl(ii+jj-1)
           xnodes(:,jj)=posl%pos(inod(jj))%posvec
        end do
        xc=0.5*(box(:,1)+box(:,2))
        hh=0.500005*abs(box(:,2)-box(:,1))
        linbox=geobj_trihitsbox(xc,hh,xnodes)
     end if
  end if

  geobj_inbox=linbox

end function geobj_inbox
!---------------------------------------------------------------------
!> is geobj in logical box
function geobj_innbox(self,posl,kcorn,kext)

  !! arguments
  type(geobj_t), intent(in) :: self   !< geobj
  type(posveclis_t), intent(in) :: posl   !< list of positions
  integer(ki2), dimension(3), intent(in) :: kcorn !< bounding box corner
  integer(ki2), dimension(3), intent(in) :: kext !< bounding box extents


  !! local
  character(*), parameter :: s_name='geobj_innbox' !< subroutine name
  logical :: geobj_innbox !< function variable type
  type(posvecl_t) :: zpos   !< position
  logical :: linbox   !< geobj in interval
  integer(ki4) :: ii   !< geobj position
  integer(ki4) :: jj   !< loop counter
  integer(ki2), dimension(3) :: ivec !< integer geobj data
  integer(ki2) :: ie !< integer element
  integer(ki2) :: iesum !< integer sum

  if (self%objtyp==1) then
     ! particle type
     ii=self%geobj
     zpos=posl%pos(ii)
     ivec=int(zpos%posvec)
     !         if (any(ivec<0)) then
     !            linbox=.false.
     !         else
     iesum=0
     do jj=1,3
        ie=ishft( ieor(kcorn(jj),ivec(jj)), -kext(jj) )
        iesum=iesum+ie
     end do
     linbox=(iesum==0)
     !         end if
     !      else if (self%objtyp==VTK_TRIANGLE) then
     !         inn=geobj_entry_table(self%objtyp)
     !! triangle type
     !         ii=self%geobj
     !         do jj=1,inn
     !         inod(jj)=nodl(ii+jj-1)
     !         xnodes(:,jj)=posl%pos(inod(jj))%posvec
     !         end do
     !         zh=0.5*float(2**kext)
     !         xc=float(kcorn)+zh
     !         hh=1.00001*zh
     !         linbox=geobj_trihitsbox(xc,hh,xnodes)
  else
     call log_error(m_name,s_name,1,error_fatal,'Unrecognised object type')
  end if

  geobj_innbox=linbox

end function geobj_innbox
!---------------------------------------------------------------------
!> does line hit geobj
subroutine geobj_linehits(self,poso,posn,posl,nodl,posh,ps,psmin,lhit)

  !! arguments
  type(geobj_t), intent(in) :: self   !< geobj definition
  type(posvecl_t), intent(in) :: poso   !< line start position
  type(posvecl_t), intent(in) :: posn   !< line end position
  type(posveclis_t), intent(in) :: posl   !< list of positions
  integer(ki4), dimension(*), intent(in) :: nodl   !< list of nodes
  type(posvecl_t), intent(out) :: posh   !< position where object hit
  real(kr4), intent(out) :: ps !< path length to object or zero if not hit
  real(kr4), intent(in) :: psmin !< minimum path length to object so far found
  logical, intent(out) :: lhit !< true if object hit

  !! local
  character(*), parameter :: s_name='geobj_linehits' !< subroutine name
  integer(ki4) :: ii   !< geobj position
  integer(ki4) :: jj   !< loop counter
  integer(ki2) :: inn !< number of nodes defining geobj
  integer(ki4), dimension(8) :: inod !< nodes of obj
  real(kr4), dimension(3,8) :: xnodes !< x(compt,node) of obj
  real(kr4), dimension(3) :: z01 !< side of triangle
  real(kr4), dimension(3) :: z02 !< side of triangle
  real(kr4), dimension(3) :: znormal !< normal to plane of polygon
  real(kr4) :: zdn !< dot product
  real(kr4) :: zd !< dot product
  real(kr4) :: zs !< local path length to object
  integer(ki2) :: i0 !< coordinate relative to dominant
  integer(ki2) :: i1 !< coordinate relative to dominant
  integer(ki2) :: i2 !< coordinate relative to dominant
  real(kr4), dimension(3) :: zp !< intersection point
  real(kr4), dimension(3) :: z0p !< intersection point relative to triangle
  real(kr4) :: zalpha !< local coordinate of intersection point
  real(kr4) :: zbeta !< local coordinate of intersection point
  integer(ki2) :: ityp !< object type

  ! defaults
  lhit=.false.
  ps=0.

  ityp=ibits(self%objtyp,0,GEOBJ_BIT_SHIFT)
  if (ityp==VTK_TRIANGLE) then
     ! triangle type
     inn=geobj_entry_table(ityp)
     ii=self%geobj
     do jj=1,inn
        inod(jj)=nodl(ii+jj-1)
        xnodes(:,jj)=posl%pos(inod(jj))%posvec
     end do
     z01=xnodes(:,2)-xnodes(:,1)
     z02=xnodes(:,3)-xnodes(:,1)
     znormal(1)=z01(2)*z02(3)-z01(3)*z02(2)
     znormal(2)=z01(3)*z02(1)-z01(1)*z02(3)
     znormal(3)=z01(1)*z02(2)-z01(2)*z02(1)
     zdn=dot_product(posn%posvec-poso%posvec,znormal)
     if (abs(zdn)<const_pusheps) then
        return
     else
        zd=-dot_product(xnodes(:,1),znormal)
        zs=-(zd+dot_product(poso%posvec,znormal))/zdn
        !           if ( zs<0. .OR. zs>psmin ) return
        if ( zs<0. .OR. zs>1.00001*psmin ) return
     end if

     ! check intersection point inside triangle
     i0=maxloc( (/abs(znormal(1)),abs(znormal(2)),abs(znormal(3))/) , dim=1  )
     i1=1+mod(i0,3)
     i2=1+mod(i1,3)
     zp=poso%posvec*(1.-zs)+posn%posvec*zs
     z0p=zp-xnodes(:,1)
     if (abs(z01(i1))<const_pusheps) then
        zbeta=z0p(i1)/z02(i1)
        if (zbeta>=0. .AND. zbeta<=1. ) then
           zalpha=(z0p(i2)-zbeta*z02(i2))/z01(i2)
           lhit= ( zalpha>=0. .AND. zalpha+zbeta<=1. )
        end if
     else
        zbeta=(z01(i1)*z0p(i2)-z0p(i1)*z01(i2)) &
 &      /(z01(i1)*z02(i2)-z02(i1)*z01(i2))
        zalpha=(z0p(i1)-zbeta*z02(i1))/z01(i1)
        !           lhit= ( zalpha>=0. .AND. zbeta>=0. .AND. zalpha+zbeta<=1. )
        lhit= ( zalpha>-const_pusheps .AND. zbeta>-const_pusheps .AND. zalpha+zbeta<1.00001 )
     end if
  end if
  ! change defaults only if hit
  if (lhit) then
     ps=zs
     posh%posvec=zp
  end if

end subroutine geobj_linehits
!---------------------------------------------------------------------
!> return .true. if face intersects binning box
!! unrolled loop version for general brick
function geobj_trihitsbox(xc,hh,xnodes)

  implicit none

  !! arguments
  real(kr4), dimension(3), intent(in) :: xc !< centre of  box
  real(kr4), dimension(3), intent(in) :: hh !< half side of box
  real(kr4), dimension(3,3), intent(in) :: xnodes !< x(compt,node) of face

  !! local
  character(*), parameter :: s_name='geobj_trihitsbox' !< subroutine name
  logical :: geobj_trihitsbox !< function variable type
  real(kr4), dimension(3) :: fnormal !< normal to face element
  real(kr4), dimension(3) :: s1 !< edge vector 1 of face element
  real(kr4), dimension(3) :: s2 !< edge vector 2 of face element
  real(kr4), dimension(3) :: s3 !< edge vector 3 of face element
  real(kr4), dimension(3) :: r1 !< box radii 1
  real(kr4), dimension(3) :: r2 !< box radii 2
  real(kr4), dimension(3) :: r3 !< box radii 3
  real(kr4) :: ndotd !< normal dot diag node
  real(kr4), dimension(3) :: vxv1 !< node 1 vertex-vertex cross products
  real(kr4), dimension(3) :: vxv2 !< node 2 vertex-vertex cross products
  real(kr4), dimension(3) :: vxv3 !< node 3 vertex-vertex cross products
  real(kr4), dimension(3) :: v1 !< tria vertex 1 relative to box centre
  real(kr4), dimension(3) :: v2 !< tria vertex 1 relative to box centre
  real(kr4), dimension(3) :: v3 !< tria vertex 1 relative to box centre
  logical :: cross !< true if triangle crosses box

  cross=.true.

  !! get relative vertex components v
  v1=xnodes(:,1)-xc
  v2=xnodes(:,2)-xc
  v3=xnodes(:,3)-xc
  !!-------------------------------------------------------
  !> mesh aligned bb test
  cross=.not.( &
 &min(v1(1),v2(1),v3(1))> hh(1) .or. &
 &max(v1(1),v2(1),v3(1))<-hh(1) .or. &
 &min(v1(2),v2(2),v3(2))> hh(2) .or. &
 &max(v1(2),v2(2),v3(2))<-hh(2) .or. &
 &min(v1(3),v2(3),v3(3))> hh(3) .or. &
 &max(v1(3),v2(3),v3(3))<-hh(3))

  !!if cross is true,  bounding boxes overlap
  if(cross) then
     !!-------------------------------------------------------
     !> now do the diagonal intersection test
     !! compute edge vectors
     s1=v3-v2
     s2=v1-v3
     s3=v2-v1
     !! compute normal (s1 x s2) and n.diag
     fnormal(1)=s1(2)*s2(3)-s1(3)*s2(2)
     fnormal(2)=s1(3)*s2(1)-s1(1)*s2(3)
     fnormal(3)=s1(1)*s2(2)-s1(2)*s2(1)
     ndotd=  abs(fnormal(1))*hh(1) &
 &   +abs(fnormal(2))*hh(2) &
 &   +abs(fnormal(3))*hh(3)

     if(abs(dot_product(fnormal,v3))>ndotd) then
        cross=.false.
     else
        cross=.true.
        !!-------------------------------------------------------
        !> finally do the 9 edge cross product tests
        !! side vectors are s-edge(compt)
        s1=abs(s1)
        s2=abs(s2)
        s3=abs(s3)

        vxv1(1)=v3(2)*v2(3)-v3(3)*v2(2)
        vxv2(1)=v1(2)*v3(3)-v1(3)*v3(2)
        vxv3(1)=v2(2)*v1(3)-v2(3)*v1(2)
        r1(1)=s1(2)*hh(3)+s1(3)*hh(2)
        r2(1)=s2(2)*hh(3)+s2(3)*hh(2)
        r3(1)=s3(2)*hh(3)+s3(3)*hh(2)
        cross=.not.( &
 &      min(vxv1(1),-vxv2(1)-vxv3(1))> r1(1) .or. &
 &      max(vxv1(1),-vxv2(1)-vxv3(1))<-r1(1) .or. &
 &      min(vxv2(1),-vxv3(1)-vxv1(1))> r2(1) .or. &
 &      max(vxv2(1),-vxv3(1)-vxv1(1))<-r2(1) .or. &
 &      min(vxv3(1),-vxv1(1)-vxv2(1))> r3(1) .or. &
 &      max(vxv3(1),-vxv1(1)-vxv2(1))<-r3(1))

        if(cross) then
           vxv1(2)=v3(3)*v2(1)-v3(1)*v2(3)
           vxv2(2)=v1(3)*v3(1)-v1(1)*v3(3)
           vxv3(2)=v2(3)*v1(1)-v2(1)*v1(3)
           r1(2)=s1(3)*hh(1)+s1(1)*hh(3)
           r2(2)=s2(3)*hh(1)+s2(1)*hh(3)
           r3(2)=s3(3)*hh(1)+s3(1)*hh(3)
           cross=.not.( &
 &         min(vxv1(2),-vxv2(2)-vxv3(2))> r1(2) .or. &
 &         max(vxv1(2),-vxv2(2)-vxv3(2))<-r1(2) .or. &
 &         min(vxv2(2),-vxv3(2)-vxv1(2))> r2(2) .or. &
 &         max(vxv2(2),-vxv3(2)-vxv1(2))<-r2(2) .or. &
 &         min(vxv3(2),-vxv1(2)-vxv2(2))> r3(2) .or. &
 &         max(vxv3(2),-vxv1(2)-vxv2(2))<-r3(2))

           if(cross)then
              vxv1(3)=v3(1)*v2(2)-v3(2)*v2(1)
              vxv2(3)=v1(1)*v3(2)-v1(2)*v3(1)
              vxv3(3)=v2(1)*v1(2)-v2(2)*v1(1)
              r1(3)=s1(1)*hh(2)+s1(2)*hh(1)
              r2(3)=s2(1)*hh(2)+s2(2)*hh(1)
              r3(3)=s3(1)*hh(2)+s3(2)*hh(1)
              cross=.not.( &
 &            min(vxv1(3),-vxv2(3)-vxv3(3))> r1(3) .or. &
 &            max(vxv1(3),-vxv2(3)-vxv3(3))<-r1(3) .or. &
 &            min(vxv2(3),-vxv3(3)-vxv1(3))> r2(3) .or. &
 &            max(vxv2(3),-vxv3(3)-vxv1(3))<-r2(3) .or. &
 &            min(vxv3(3),-vxv1(3)-vxv2(3))> r3(3) .or. &
 &            max(vxv3(3),-vxv1(3)-vxv2(3))<-r3(3))
           end if
        end if
     end if
  end if
  geobj_trihitsbox=cross
end function geobj_trihitsbox
!---------------------------------------------------------------------
!> calculate normal to geobj
subroutine geobj_normal(self,posl,nodl,pnormal,pmag)

  !! arguments
  type(geobj_t), intent(in) :: self   !< geobj definition
  type(posveclis_t), intent(in) :: posl   !< list of positions
  integer(ki4), dimension(*), intent(in) :: nodl   !< list of nodes
  real(kr4), dimension(3), intent(out) :: pnormal !< unit normal vector
  real(kr4), intent(out) :: pmag !< magnitude of normal vector

  !! local
  character(*), parameter :: s_name='geobj_normal' !< subroutine name
  integer(ki4) :: ii   !< geobj position
  integer(ki4) :: jj   !< loop counter
  integer(ki2) :: inn !< number of nodes defining geobj
  integer(ki4), dimension(8) :: inod !< nodes of obj
  real(kr4), dimension(3,8) :: xnodes !< x(compt,node) of obj
  real(kr4), dimension(3) :: z01 !< side of triangle
  real(kr4), dimension(3) :: z02 !< side of triangle
  real(kr4) :: zdsq !< dot product squared
  real(kr4) :: zd !< dot product
  real(kr4) :: zref !< reference length based on maximum component length
  integer(ki2) :: ityp !< object type

  ityp=ibits(self%objtyp,0,GEOBJ_BIT_SHIFT)
  if (ityp==VTK_TRIANGLE) then
     ! triangle type
     inn=geobj_entry_table(ityp)
     ii=self%geobj
     do jj=1,inn
        inod(jj)=nodl(ii+jj-1)
        xnodes(:,jj)=posl%pos(inod(jj))%posvec
     end do
     z01=xnodes(:,2)-xnodes(:,1)
     z02=xnodes(:,3)-xnodes(:,1)
     pnormal(1)=z01(2)*z02(3)-z01(3)*z02(2)
     pnormal(2)=z01(3)*z02(1)-z01(1)*z02(3)
     pnormal(3)=z01(1)*z02(2)-z01(2)*z02(1)
     zdsq=dot_product(pnormal,pnormal)
     zd=sqrt( max(0.,zdsq) )
     zref=max( abs(maxval(z01)), abs(minval(z01) ), abs(maxval(z02)), abs(minval(z02)) )
     if (zd<const_pusheps*zref) then
        call log_error(m_name,s_name,1,error_warning,'triangle is degenerate')
     else
        pnormal=pnormal/zd
     end if
  else
     call log_error(m_name,s_name,10,error_warning,'object is not a triangle')
     pnormal=0
     zd=1
  end if

  pmag=zd

end subroutine geobj_normal
!---------------------------------------------------------------------
!> calculate area of geobj
subroutine geobj_area(self,posl,nodl,parea)

  !! arguments
  type(geobj_t), intent(in) :: self   !< geobj definition
  type(posveclis_t), intent(in) :: posl   !< list of positions
  integer(ki4), dimension(*), intent(in) :: nodl   !< list of nodes
  real(kr4), intent(out) :: parea !< area

  !! local
  character(*), parameter :: s_name='geobj_area' !< subroutine name
  integer(ki4) :: ii   !< geobj position
  integer(ki4) :: jj   !< loop counter
  integer(ki2) :: inn !< number of nodes defining geobj
  integer(ki4), dimension(8) :: inod !< nodes of obj
  real(kr4), dimension(3,8) :: xnodes !< x(compt,node) of obj
  real(kr4), dimension(3) :: z01 !< side of triangle
  real(kr4), dimension(3) :: z02 !< side of triangle
  real(kr4), dimension(3) :: zcx1 !< cross-product
  real(kr4), dimension(3) :: zcx2 !< cross-product
  real(kr4) :: zarea !< cumulative area
  real(kr4) :: zd !< dot product
  real(kr4) :: zref !< reference length based on maximum component length
  integer(ki2) :: ityp !< object type

  ityp=ibits(self%objtyp,0,GEOBJ_BIT_SHIFT)
  if (ityp==VTK_TRIANGLE.OR.ityp==VTK_QUAD) then
     ! triangle or quad type
     inn=geobj_entry_table(ityp)
     ii=self%geobj
     do jj=1,inn
        inod(jj)=nodl(ii+jj-1)
        xnodes(:,jj)=posl%pos(inod(jj))%posvec
     end do
     z01=xnodes(:,2)-xnodes(:,1)
     z02=xnodes(:,3)-xnodes(:,1)
     call cross_product(z01,z02,zcx1)
     zarea=sqrt(dot_product(zcx1,zcx1))
     if (inn==4) then
        z01=xnodes(:,2)-xnodes(:,4)
        z02=xnodes(:,3)-xnodes(:,4)
        call cross_product(z01,z02,zcx1)
        zarea=zarea+sqrt(dot_product(zcx1,zcx1))
     end if
     zref=max( abs(maxval(z01)), abs(minval(z01) ), abs(maxval(z02)), abs(minval(z02)) )
     if (zarea<const_pusheps*zref) then
        call log_error(m_name,s_name,1,error_warning,'object is degenerate')
        zarea=const_pusheps*zref
     end if
  else
     call log_error(m_name,s_name,10,error_warning,'object is not a triangle or quadrilateral')
     zarea=0
  end if

  parea=zarea/2

end subroutine geobj_area
!---------------------------------------------------------------------
!> calculate nodes to geobj
subroutine geobj_nodes(self,posl,nodl,pnodes)

  !! arguments
  type(geobj_t), intent(in) :: self   !< geobj definition
  type(posveclis_t), intent(in) :: posl   !< list of positions
  integer(ki4), dimension(*), intent(in) :: nodl   !< list of nodes
  real(kr4), dimension(3,8), intent(out) :: pnodes !< nodes vector

  !! local
  character(*), parameter :: s_name='geobj_nodes' !< subroutine name
  integer(ki4) :: ii   !< geobj position
  integer(ki4) :: jj   !< loop counter
  integer(ki2) :: inn !< number of nodes defining geobj
  integer(ki4), dimension(8) :: inod !< nodes of obj
  integer(ki2) :: ityp !< object type

  ityp=ibits(self%objtyp,0,GEOBJ_BIT_SHIFT)
  inn=geobj_entry_table(ityp)
  if (inn>=1.AND.inn<=8) then
     ii=self%geobj
     do jj=1,inn
        inod(jj)=nodl(ii+jj-1)
        pnodes(:,jj)=posl%pos(inod(jj))%posvec
     end do
  else
     call log_error(m_name,s_name,1,error_warning,'no or too many nodes')
     pnodes=0
  end if

end subroutine geobj_nodes
!---------------------------------------------------------------------
!> calculate barycentre of geobj
subroutine geobj_centre(self,posl,nodl,pcentr)

  !! arguments
  type(geobj_t), intent(in) :: self   !< geobj definition
  type(posveclis_t), intent(in) :: posl   !< list of positions
  integer(ki4), dimension(*), intent(in) :: nodl   !< list of nodes
  real(kr4), dimension(3), intent(out) :: pcentr !< centre

  !! local
  character(*), parameter :: s_name='geobj_centre' !< subroutine name
  integer(ki4) :: ii   !< geobj position
  integer(ki4) :: jj   !< loop counter
  integer(ki2) :: inn !< number of nodes defining geobj
  integer(ki4), dimension(8) :: inod !< nodes of obj
  real(kr4), dimension(3,8) :: xnodes !< x(compt,node) of obj
  integer(ki2) :: ityp !< object type

  ityp=ibits(self%objtyp,0,GEOBJ_BIT_SHIFT)
  if (ityp==VTK_TRIANGLE) then
     ! triangle type
     inn=geobj_entry_table(ityp)
     ii=self%geobj
     do jj=1,inn
        inod(jj)=nodl(ii+jj-1)
        xnodes(:,jj)=posl%pos(inod(jj))%posvec
     end do
     pcentr=(xnodes(:,1)+xnodes(:,2)+xnodes(:,3))/3
  else
     call log_error(m_name,s_name,10,error_warning,'object is not a triangle')
  end if

end subroutine geobj_centre
!---------------------------------------------------------------------
!> calculate sample point of geobj
subroutine geobj_sample(self,posl,nodl,pt,psampl)

  !! arguments
  type(geobj_t), intent(in) :: self   !< geobj definition
  type(posveclis_t), intent(in) :: posl   !< list of positions
  integer(ki4), dimension(*), intent(in) :: nodl   !< list of nodes
  real(kr4), dimension(3), intent(in) :: pt !< sample point
  real(kr4), dimension(3), intent(out) :: psampl !< sample value

  !! local
  character(*), parameter :: s_name='geobj_sample' !< subroutine name
  integer(ki4) :: ii   !< geobj position
  integer(ki4) :: jj   !< loop counter
  integer(ki2) :: inn !< number of nodes defining geobj
  integer(ki4), dimension(8) :: inod !< nodes of obj
  real(kr4), dimension(3,8) :: xnodes !< x(compt,node) of obj
  integer(ki2) :: ityp !< object type

  ityp=ibits(self%objtyp,0,GEOBJ_BIT_SHIFT)
  if (ityp==VTK_TRIANGLE) then
     ! triangle type
     inn=geobj_entry_table(ityp)
     ii=self%geobj
     do jj=1,inn
        inod(jj)=nodl(ii+jj-1)
        xnodes(:,jj)=posl%pos(inod(jj))%posvec
     end do
     psampl=pt(1)*xnodes(:,1)+pt(2)*xnodes(:,2)+pt(3)*xnodes(:,3)
  else
     call log_error(m_name,s_name,10,error_warning,'object is not a triangle')
  end if

end subroutine geobj_sample
!---------------------------------------------------------------------
!> return number of entries in geobj_entry_table
function geobj_entry_table_fn(ktyp)

  !! arguments
  integer(ki2), intent(in) :: ktyp   !< object type

  !! local
  character(*), parameter :: s_name='geobj_entry_table_fn' !< subroutine name
  integer(ki2) :: geobj_entry_table_fn !< function variable type
  integer(ki2) :: ityp   !< geobj type unscrambled

  ityp=ibits(ktyp,0,GEOBJ_BIT_SHIFT)
  geobj_entry_table_fn=geobj_entry_table(ityp)

end function geobj_entry_table_fn
!---------------------------------------------------------------------
!> extract last 8 bits in geobj_entry_table
function geobj_type_fn(ktyp)

  !! arguments
  integer(ki2), intent(in) :: ktyp   !< object type

  !! local
  character(*), parameter :: s_name='geobj_type_fn' !< subroutine name
  integer(ki2par) :: geobj_type_fn !< function variable type
  integer(ki2par) :: ityp   !< geobj type unscrambled

  ityp=ibits(ktyp,0,GEOBJ_BIT_SHIFT)
  geobj_type_fn=ityp

end function geobj_type_fn
!---------------------------------------------------------------------
!> get geometry code
function geobj_code_fn(ktyp)

  !! arguments
  integer(ki2), intent(in) :: ktyp   !< object type

  !! local
  character(*), parameter :: s_name='geobj_code_fn' !< subroutine name
  integer(ki2par) :: geobj_code_fn !< function variable type
  integer(ki2par) :: igcode   !< code for geometry

  igcode=ishft(ktyp,-GEOBJ_BIT_SHIFT)
  geobj_code_fn=igcode

end function geobj_code_fn
!---------------------------------------------------------------------
!> vector cross product (kr4) of the input vectors a and b
subroutine cross_product(a,b,crossout)

  !! arguments
  real(kr4), dimension(3), intent(in) :: a !< input vector 1
  real(kr4), dimension(3), intent(in) :: b !< input vector 2
  real(kr4), dimension(3), intent(out) :: crossout !< output vector

  crossout(1) = a(2)*b(3)-a(3)*b(2)
  crossout(2) = a(3)*b(1)-a(1)*b(3)
  crossout(3) = a(1)*b(2)-a(2)*b(1)

end subroutine cross_product

end module geobj_m
