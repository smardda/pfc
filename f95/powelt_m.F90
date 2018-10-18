module powelt_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use vfile_m
  use btree_m
  use powcal_h
  use powres_h
  use edgprof_h
  use edgprof_m
  use position_h
  use position_m
  use geobj_m
  use spl2d_m
  use spl3d_m
  use posang_h
  use posang_m
  use odes_h
  use odes_m
  use beq_m
  use pcle_h
  use pcle_m
!     use termplane_h

  implicit none
  private

! public subroutines
  public :: &
  powelt_init, & !< create tables for element centre addressing
  powelt_pos, & !< position of element centre in tracking coordinates
  powelt_vec,   & !< magnetic field in cartesians at centre of element
  powelt_normal, & !< normal in cartesians at centre of element
  powelt_normalm, & !< normal in tracking coordinates at centre of element
  powelt_dep, & !< power deposition formula of element
  powelt_move, & !< motion along fieldline local, axisymm, flux coords
  powelt_move0, & !< motion along fieldline local, axisymm, cylindricals
  powelt_move1, & !< motion along fieldline global, axisymm, cylindricals
  powelt_move2, & !< motion along fieldline global, analytic ripple
  powelt_move3, & !< motion along fieldline global, ripple field OBSOLETE
  powelt_move4, & !< motion along fieldline global, ripple field
  powelt_move5, & !< motion along fieldline 5th scheme for testing
  powelt_avg, & !< average value of power on element
  powelt_devn, & !< normalised deviation of power on element
  powelt_stat, & !< average deviation of power on element
  powelt_write, & !< write data associated with element
  powelt_writev !< write (vtk) data associated with element
! private subroutines
  private :: &
  powelt_addr, & !< address in powcal structure
  powelt_setpow !< save data needed for power deposition calculation

! public types
!> data structure describing powelt
  type, public :: powelt_t
     integer(ki4)  :: ie !< geometrical object
     integer(ki4)  :: je !< sub-element object
  end type powelt_t

! public variables
  integer(ki4), parameter, public :: powelt_max_baryc_table = 16 !< local variable
  integer(ki4), parameter, public :: powelt_max_powelt_table = 3 !< local variable
  real(kr8), dimension(3,powelt_max_baryc_table), public :: powelt_weights !< local variable
  real(kr8), dimension(powelt_max_baryc_table), public :: powelt_level_table !< local variable
  integer(ki4), dimension(3,3),public :: powelt_table !< local variable

! private types

! private variables
  character(*), parameter :: m_name='powelt_m' !< module name
  logical, parameter :: lpmidplane=.FALSE. !< test for midplane intersection
  logical :: ilpmidplane   !< is 'midplane' test active
  logical, parameter :: loutdt=.FALSE. !< output for timestep analysis
  real(kr8), dimension(:), allocatable :: work1 !< 1D work array
  type(posveclis_t) :: wposl   !< list of position data
  integer   :: status   !< error status
  integer :: nplot   !< output channel for powelt data
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: idum !< dummy integer
  logical, parameter :: debug=.FALSE. !< logical flag for debugging
  integer(ki4) :: ilevel   !< refinement level for address
  integer(ki4) :: iinlevel   !< number of subelements at refinement level
  integer(ki4) :: inpow   !< address of element
  integer(ki4) :: iiobj   !< geobj position
  integer(ki4), save :: nskyl1=0   !< number of starting points "skylighted" 1
  integer(ki4), save :: nskyl2=0   !< number of starting points "skylighted" 2
  integer(ki2) :: inn !< number of nodes defining geobj
  integer(ki4), dimension(8) :: inod !< nodes of obj
  real(kr4), dimension(3,8) :: xnodes !< x(compt,node) of obj
  character(len=80) :: ibuff !< buffer for input/output
  character(len=80) :: icfile !< local file-name

  contains
!---------------------------------------------------------------------
!> position of element centre in tracking coordinates
subroutine powelt_init

  !! local
  character(*), parameter :: s_name='powelt_init' !< subroutine name
  integer(ki4), dimension(powelt_max_baryc_table), parameter :: baryc_table1 = & !< local variable
 &(/4,8,2,2,2,5,5,10,7,7,1,1,4,4,1,1/)
  integer(ki4), dimension(powelt_max_baryc_table), parameter :: baryc_table2 = & !< local variable
 &(/4,2,8,2,5,2,5,1,4,1,7,10,7,1,4,1/)
  real(kr8), parameter :: baryc_norm=12_kr8 !< local variable
  integer(ki4), dimension(powelt_max_powelt_table*powelt_max_powelt_table), parameter :: table = & !< local variable
 &(/0,1,4,1,4,16,1,3,12/)
  integer(ki4), dimension(powelt_max_baryc_table), parameter :: level_table = & !< local variable
 &(/1,(2,j=1,3),(3,j=1,12)/)
  !       1 2 3 4 5 6 7 8  9 0 1 2  3 4 5 6  7 8 9 0 1 2 3 4  5  6  7

  powelt_table=reshape(table, (/powelt_max_powelt_table,powelt_max_powelt_table/) )
  ! note in version where all triangles pre-computed, only i=1 entry is used ever
  do i=1,powelt_max_baryc_table
     powelt_weights(1,i)=baryc_table1(i)/baryc_norm
     powelt_weights(2,i)=baryc_table2(i)/baryc_norm
     powelt_weights(3,i)=1-powelt_weights(1,i)-powelt_weights(2,i)
  end do
  powelt_level_table=level_table

end subroutine powelt_init
!---------------------------------------------------------------------
!> position of element centre in tracking coordinates
subroutine powelt_pos(self,powcal,pos)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(in) :: powcal !< powcal data structure
  real(kr4), dimension(3), intent(out) :: pos !< output position vector

  !! local
  character(*), parameter :: s_name='powelt_pos' !< subroutine name
  type(geobj_t) :: iobj !< object

  inpow=powelt_addr(self,powcal%powres%npowe)
  iobj%geobj=powcal%powres%geobjl%obj2(inpow)%ptr
  iobj%objtyp=powcal%powres%geobjl%obj2(inpow)%typ
  call geobj_nodes(iobj,powcal%powres%geobjl%posl,powcal%powres%geobjl%nodl,xnodes)
  pos=powelt_weights(1,1)*xnodes(:,1)+&
 &powelt_weights(2,1)*xnodes(:,2)+powelt_weights(3,1)*xnodes(:,3)

end subroutine powelt_pos
!---------------------------------------------------------------------
!> magnetic field in cartesians at centre of element
subroutine powelt_vec(self,powcal,pb)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(in) :: powcal !< powcal data structure
  real(kr4), dimension(3), intent(out) :: pb !< output vector

  !! local
  character(*), parameter :: s_name='powelt_vec' !< subroutine name

  inpow=powelt_addr(self,powcal%powres%npowe)
  iiobj=powcal%powres%geobjl%obj2(inpow)%ptr
  inn=geobj_entry_table_fn(powcal%powres%geobjl%obj2(inpow)%typ)
  if (inn/=3) then
     call log_error(m_name,s_name,1,error_warning,'object is not a triangle')
     pb=0
  else
     do j=1,inn
        inod(j)=powcal%powres%geobjl%nodl(iiobj+j-1)
        xnodes(:,j)=powcal%powres%vecb%pos(inod(j))%posvec
     end do
     pb=powelt_weights(1,1)*xnodes(:,1)+&
 &   powelt_weights(2,1)*xnodes(:,2)+powelt_weights(3,1)*xnodes(:,3)
  end if

end subroutine powelt_vec
!---------------------------------------------------------------------
!> normal in cartesians at centre of element
subroutine powelt_normal(self,powcal,pnormal)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(in) :: powcal !< powcal data structure
  real(kr4), dimension(3), intent(out) :: pnormal !< normal vector

  !! local
  character(*), parameter :: s_name='powelt_normal' !< subroutine name
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4) :: zmag !< magnitude of normal vector
  type(geobj_t) :: iobj !< object

  inpow=powelt_addr(self,powcal%powres%npowe)
  iobj%geobj=powcal%powres%geobjl%obj2(inpow)%ptr
  iobj%objtyp=powcal%powres%geobjl%obj2(inpow)%typ
  call geobj_normal(iobj,powcal%powres%vecx,powcal%powres%geobjl%nodl,znormal,zmag)

  pnormal=znormal

end subroutine powelt_normal
!---------------------------------------------------------------------
!> normal in tracking coordinates at centre of element
subroutine powelt_normalm(self,powcal,pnormal)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(in) :: powcal !< powcal data structure
  real(kr4), dimension(3), intent(out) :: pnormal !< normal vector

  !! local
  character(*), parameter :: s_name='powelt_normalm' !< subroutine name
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4) :: zmag !< magnitude of normal vector
  type(geobj_t) :: iobj !< object

  inpow=powelt_addr(self,powcal%powres%npowe)
  iobj%geobj=powcal%powres%geobjl%obj2(inpow)%ptr
  iobj%objtyp=powcal%powres%geobjl%obj2(inpow)%typ
  call geobj_normal(iobj,powcal%powres%geobjl%posl,powcal%powres%geobjl%nodl,znormal,zmag)

  pnormal=znormal

end subroutine powelt_normalm
!---------------------------------------------------------------------
!> power deposition formula of element
subroutine powelt_dep(self,powcal,gshadl)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure
  type(geobjlist_t), intent(in), optional :: gshadl   !< shadow data structure

  !! local
  character(*), parameter :: s_name='powelt_dep' !< subroutine name
  real(kr4) :: zpow !< power deposited
  real(kr4), dimension(3) :: zpos !< position vector
  type(posvecl_t) :: zpos1   !< one position data
  type(posvecl_t) :: zpos2   !< position data dequantised
  real(kr4), dimension(3) :: zb !< magnetic field vector
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4) :: zbdotn !< \f$ \bf{B} .\bf{n} \f$
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  real(kr8) :: zpsi    !<  \f$ \psi(R,Z) \f$
  real(kr8) :: zzeta    !<  \f$ \zeta(R,Z) \f$
  real(kr8) :: zpsid    !<  \f$ \psi(R,Z)-\psi_b \f$
  real(kr4) :: zmag     !< local variable
  integer(ki4) :: iflag !< status of element
  type(geobj_t) :: iobj !< object
  type(posang_t) :: zposang   !< posang data structure

  inpow=powelt_addr(self,powcal%powres%npowe)

  iflag=powcal%powres%pow(inpow)+0.5
  if (iflag==2) then
     ! fieldline was not followed to termination
     calcn_type0: select case (powcal%n%caltype)
     case('afws','local','msus','global')
        ! assume unshadowed
     case('msum','middle')
        ! for mid-lane launch, assumed no power deposited  and exit
        powcal%powres%pow(inpow)=0
        return
     end select calcn_type0
  end if

  if (iflag/=0) then
     ! element is not shadowed
     call powelt_pos(self,powcal,zpos)
     zpos1%posvec=zpos

     calcn_type: select case (powcal%n%caltype)
     case('afws','local','msus','global')
        ! evaluate B.n from element
        call powelt_vec(self,powcal,zb)
        call powelt_normal(self,powcal,znormal)
        zbdotn=dot_product(zb,znormal)
     case('msum','middle')
        field_type: select case (powcal%powres%beq%n%fldspec)
        case(1)
           call log_error(m_name,s_name,6,error_fatal,'No code for case')
        case default
           ! use saved object number to get normal and position
           ! only OK if ripple defined via F
           iobj%geobj=powcal%powres%geobjl%obj(inpow)%weight+0.5
           !        write(*,*) 'iog=',iobj%geobj
           iobj%objtyp=gshadl%obj2(iobj%geobj)%typ
           call geobj_normal(iobj,gshadl%posl,gshadl%nodl,znormal,zmag)
           call geobj_centre(iobj,gshadl%posl,gshadl%nodl,zpos)
           ! calculate F at position
           zposang%pos=zpos
           call beq_b(powcal%powres%beq,zposang,1)
           ! F vector to cartesians
           ! dequantise position
           zpos1%posvec=zpos
           zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
           ! convert xi to zeta
           zzeta=(zpos2%posvec(3))/powcal%powres%beq%nzets
           zpos2%posvec(3)=zzeta
           zposang%pos=zpos2%posvec
           zposang%opt=1 ; zposang%units=0
           ! convert F to field vector
           zr=zpos2%posvec(1)
           zposang%vec(1)=zposang%vec(1)/zr ; zposang%vec(2)=zposang%vec(2)/zr
           call posang_tfm(zposang,-3)
           zb=zposang%vec
           zbdotn=dot_product(zb,znormal)
        end select field_type
     end select calcn_type


     calcn_type2: select case (powcal%n%caltype)
     case('afws','local','msum','middle')
        field_type2: select case (powcal%powres%beq%n%fldspec)
        case(1)
           ! psi is the position (unquantised)
           zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
           zpsi=zpos2%posvec(1)
        case default
           zr=zpos1%posvec(1)
           zz=zpos1%posvec(2)
           call spl2d_eval(powcal%powres%beq%psi,zr,zz,zpsi)
        end select field_type2
     case('msus','global')
        ! make saved psi available
        if (iflag==1) then
           zpsi=powcal%powres%geobjl%obj(inpow)%weight
        end if
        field_type3: select case (powcal%powres%beq%n%fldspec)
        case(1,2)
           if (iflag/=1) then
              zpsi=powcal%powres%geobjl%obj(inpow)%weight
              ! psi is the position (unquantised)
              zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
              zpsi=zpos2%posvec(1)
           end if
        case default
           ! use elemental psi
           zr=zpos1%posvec(1)
           zz=zpos1%posvec(2)
           call spl2d_eval(powcal%powres%beq%psi,zr,zz,zpsi)
        end select field_type3
        powcal%powres%geobjl%obj(inpow)%weight=iflag
     end select calcn_type2

     ! psibdry is not quantised
     zpsid=zpsi-powcal%powres%beq%psibdry
     zpow=zbdotn*edgprof_fn(powcal%edgprof,zpsid)
     powcal%powres%pow(inpow)=zpow
     powcal%powres%pmask(inpow)=nint(sign(1._kr4,zpow))
     powcal%powres%psista(inpow)=zpsi
     !        write(*,*) inpow, zpow, zbdotn, zpsi
  end if

end subroutine powelt_dep
!---------------------------------------------------------------------
!> motion along fieldline local, axisymm, flux coords
subroutine powelt_move(self,powcal,gshadl,btree)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure
  type(btree_t), intent(inout),optional :: btree !< btree data

  !! local
  character(*), parameter :: s_name='powelt_move' !< subroutine name
  integer(ki4), parameter :: ipback=0 !< back-track test if unity
  real(kr4), dimension(3) :: zpos !< start position vector
  real(kr8), dimension(3) :: zposd !< start position vector
  real(kr8), dimension(3) :: xpath !< start path position vector Cartesians
  real(kr8) :: zpsi !< value of \f$ \psi \f$
  real(kr8) :: ztheta !< value of \f$ \theta \f$
  real(kr8) :: zf !< value of \f$ f \f$
  real(kr8) :: zt0 !< start value of \f$ t \f$
  type(posang_t) :: zposang   !< posang data structure
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4), dimension(3) :: znormalm !< unit normal vector in mapped coords
  real(kr8), dimension(3) :: znormald !< unit normal vector
  real(kr8) :: lenpath !< length of path in \f$ t \f$
  real(kr8) :: zdfac !< direction factor
  real(kr8) :: zs !< \f$ \pm 1 \f$
  real(kr8) :: zpdotn !< path dotted with normal vector
  real(kr8) :: zpdotnm !< path dotted with normal vector in mapped coords
  integer(ki4) :: ierr !< error flag
  integer(ki4) :: ip !< entries in track
  type(geobj_t) :: iobj !< geo object
  integer(ki4) :: inode  !< local variable
  !> Non-zero at fieldline termination
  !! - (-2) Line hit terminating plane
  !! - (-1) Line left octree
  !! - Positive - number of triangle hit
  integer(ki4) :: nobjhit=0  !< .
  integer(ki4) :: ibacktr !< flag if back-tracking
  type(posnode_t) :: xo !< position variable
  type(posnode_t) :: xn !< position variable
  type(posvecl_t) :: zpos1   !< one position data
  type(posvecl_t) :: zpos2   !< one position data
  type(posveclis_t) :: rposl   !< list of position data
  real(kr8), dimension(2) :: zk !< sector number
  type(posnode_t) :: xm !< local variable
  real(kr8) :: zzeta !< value of \f$ \zeta \f$
  real(kr8) :: zpsim !< value of \f$ \psi \f$ at field line end (diagnostic)
  logical :: lcoll   !< collision on path

  ! find element number
  inpow=powelt_addr(self,powcal%powres%npowe)
  !dbgwr write(*,*) "inpow=",inpow !dbgwr
  !  write(*,*) powcal%n%ntrack !dbgwr
  !  write(*,*) powcal%n%trackno !dbgwr
  ! check whether looking at small number of tracks
  if (powcal%n%ntrack>0) then
     ! is this element track in the test set
     do j=1,powcal%n%ntrack
        if (powcal%n%trackno(j)==inpow) go to 1
     end do
     return
  end if
1     continue
  ! mark powelt as unshadowed by default
  powcal%powres%pow(inpow)=1
  ! needed in shadowed case
  if (powcal%n%shadow>0.AND..NOT.allocated(rposl%pos)) then
     allocate(rposl%pos(1),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else if (powcal%n%shadow<-1) then
     return
  end if

  ! calculate, save and test initial position
  call powelt_pos(self,powcal,zpos)
  zposd=zpos
  zt0=zposd(3)
  zpsi=zposd(1)
  ztheta=zposd(2)
  ! check psi lies within range
  if ((zpsi-powcal%powres%psimin<0).OR.(zpsi-powcal%powres%psimax>0)) then
     call log_error(m_name,s_name,1,error_warning,'Starting psi out of range')
     return
  end if
  ! check theta lies within range
  if ((ztheta-powcal%powres%thetamin<0).OR.(ztheta-powcal%powres%thetamax>0)) then
     call log_error(m_name,s_name,2,error_warning,'Starting theta out of range')
     return
  end if
  ! Examine start trajectory behaviour in mapped coords
  call powelt_normalm(self,powcal,znormalm)
  znormald=znormalm
  ! start trajectory
  powcal%odes%ndt=1
  powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=zposd
  ! evaluate I aka f at psi
  call spleval(powcal%powres%beq%f,powcal%powres%beq%mr,powcal%powres%beq%psiaxis,&
  powcal%powres%beq%psiqbdry,zpsi,zf,1)
  ! set initial time step and guess first new point of track
  call odes_rjstep1(powcal%odes,ierr,powcal%powres%qfac*zf,powcal%powres%beq%rjac)
  xpath=powcal%odes%vecp%pos(powcal%odes%ndt+1)%posvec-zposd
  zpdotnm=dot_product(xpath,znormald)

  ! Check behaviour in Cartesians
  call powelt_normal(self,powcal,znormal)
  znormald=znormal
  ! ensure path is outward from surface
  xpath=0
  zs=1
  do j=0,1
     ! dequantise
     zpos1%posvec=powcal%odes%vecp%pos(powcal%odes%ndt+j)%posvec
     zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
     zposang%pos=zpos2%posvec
     zposang%opt=2 ; zposang%units=0
     call posang_invpsitfm(zposang,powcal%powres%beq)
     call posang_tfm(zposang,-3)
     xpath=zposang%pos+zs*xpath
     zs=-zs
  end do
  zpdotn=dot_product(xpath,znormald)
  if (zpdotn*zpdotnm*beq_rsig()<0) then
     call log_error(m_name,s_name,3,error_warning,'Starting direction uncertain')
     !        zpdotn=sign(1._kr8,zpdotnm)
     write(*,*) inpow, zpdotn, zpdotnm !DPRT
  end if
  !     zdfac=sign(1._kr8,zpdotn)
  zdfac=sign(powcal%powres%qfac,zpdotn)
  if (powcal%n%shadow>0) then
     ! find start point in tree
     rposl%pos(1)=powcal%odes%vecp%pos(powcal%odes%ndt)
     rposl%np=1
     iobj%geobj=1
     iobj%objtyp=1 ! for point
     call btree_mfind(btree,iobj,rposl,inode)
     if (inode<0) then
        call log_error(m_name,s_name,4,error_warning,'cannot locate start point in tree')
        return
     end if
     xo%posvec=rposl%pos(1)%posvec
     xo%node=inode
  end if
  ibacktr=0

100   continue
  ! time step loop for path
  lenpath=0
  ierr=0
  lcoll=.FALSE.
  !     powcal%odes%near=0
  loop_path: do

     lenpath=lenpath+abs(powcal%odes%dt)

     call odes_rjstep(powcal%odes,ierr,zf,zdfac,powcal%powres%beq%rjac)
     ! check for termination due to solution failure on this field-line
     if (ierr>0) exit

     if (powcal%n%shadow>0) then
        ! check for collision or exit from computational domain
        xn%posvec=powcal%odes%vecp%pos(powcal%odes%ndt)%posvec
        xn%node=0
        if (abs(lenpath)>powcal%odes%n%termcon(1)) then
           ! find new node and position if hits boundary
           call pcle_move(xo,xn,1,gshadl,btree,nobjhit)
           lcoll=(nobjhit/=0)
        end if
        if (lcoll) then
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
           ! mark as really a collision, not just due to leaving volume
           if (nobjhit>0)  powcal%powres%pow(inpow)=0
           !D           write(*,*) 'Final posn',xn%posvec! writediagn
           exit
        end if
        ! save new data
        xo=xn
        ! check for near-collision and repeat step (once) if necessary TO DO???
     end if

     ! check for termination due to overlong trajectory (number of steps)
     if (powcal%odes%ndt>=powcal%odes%n%stepmax) exit
     ! this one covered when test exit from domain
     ! if (abs(powcal%odes%vecp%pos(powcal%odes%ndt)%posvec(3)-zt0)>abs(powcal%odes%n%tmax-zt0)) exit
     if (abs(lenpath)>abs(powcal%odes%n%tmax-zt0)) exit

  end do loop_path


  ip=powcal%odes%ndt
  powcal%odes%vecp%np=ip
  !
  if (ibacktr==ipback.AND.powcal%powres%flinm) then
     ! open file to record mapped track
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' trackm'')') self%ie,self%je
     write(icfile,'(''trackm'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     ! write track in mapped coordinates
     ! de-quantise
     call position_invqtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     call position_writelis(powcal%odes%vecp,'track',nplot)
     ! re-quantise
     call position_qtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     ! close file
     close(nplot)
  end if

  if (ibacktr==ipback.AND.powcal%powres%flinx) then
     ! open file to record Cartesian track
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' trackx'')') self%ie,self%je
     write(icfile,'(''trackx'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     ! write track in Cartesian coordinates???
     allocate(wposl%pos(ip), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     do l=1,ip
        ! dequantise
        zpos1%posvec=powcal%odes%vecp%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        zposang%pos=zpos2%posvec
        zposang%opt=2 ; zposang%units=0
        call posang_invpsitfm(zposang,powcal%powres%beq)
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=ip
     ! write track
     call position_writelis(wposl,'track',nplot)
     ! close file and deallocate
     deallocate(wposl%pos)
     close(nplot)
  end if

  if (nobjhit<=GEOBJ_OFFSET.AND.ibacktr==ipback.AND.powcal%powres%flinends) then
     ! write first and last points in Cartesian coordinates
     ! first convert track points
     allocate(wposl%pos(3), stat=status)
     call log_alloc_check(m_name,s_name,11,status)
     wposl%pos(1)%posvec=powcal%odes%vecp%pos(1)%posvec
     zk(1)=powcal%odes%posk(1)
     wposl%pos(2)%posvec=xm%posvec
     zk(2)=powcal%odes%posk(ip)
     do l=1,2
        ! dequantise
        zpos1%posvec=wposl%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*zk(l))/powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=3
     ! open file to record Cartesian track and write
     wposl%pos(3)%posvec(1)=self%ie
     wposl%pos(3)%posvec(2)=self%je
     wposl%pos(3)%posvec(3)=zpsim
     do l=1,3
        call position_writev(wposl%pos(l),powcal%powres%nflends)
     end do
     ! deallocate
     deallocate(wposl%pos)
  end if

  if (nobjhit<0.AND.ibacktr<ipback) then
     ibacktr=ipback
     zdfac=-zdfac
     powcal%odes%ndt=powcal%odes%ndt-1
     goto 100
  end if

end subroutine powelt_move
!---------------------------------------------------------------------
!> motion along fieldline local, axisymm, cylindricals
subroutine powelt_move0(self,powcal,gshadl,btree)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure
  type(btree_t), intent(inout),optional :: btree !< btree data

  !! local
  character(*), parameter :: s_name='powelt_move0' !< subroutine name
  integer(ki4), parameter :: ipback=0 !< back-track test if unity
  real(kr4), dimension(3) :: zpos !< start position vector
  real(kr8), dimension(3) :: zposd !< start position vector
  real(kr8), dimension(3) :: xpath !< start path position vector Cartesians
  real(kr8) :: zpsi !< value of \f$ \psi \f$
  real(kr8) :: zzeta !< value of \f$ \zeta \f$
  real(kr8) :: zt0 !< start value of \f$ t \f$
  type(posang_t) :: zposang   !< posang data structure
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4), dimension(3) :: znormalm !< unit normal vector in mapped coords
  real(kr8), dimension(3) :: znormald !< unit normal vector
  real(kr8), dimension(6) :: zdfaca !< direction factor and offsets
  real(kr8) :: zs !< \f$ \pm 1 \f$
  real(kr8) :: zpdotn !< path dotted with normal vector
  real(kr8) :: zpdotnm !< path dotted with normal vector in mapped coords
  integer(ki4) :: ierr !< error flag
  integer(ki4) :: ip !< entries in track
  type(geobj_t) :: iobj !< geo object
  integer(ki4) :: inode  !< local variable
  integer(ki4) :: indt !< shorthand for current timestep number
  !> Non-zero at fieldline termination
  !! - (-2) Line hit terminating plane
  !! - (-1) Line left octree
  !! - Positive - number of triangle hit
  integer(ki4) :: nobjhit=0  !< .
  integer(ki4) :: ibacktr !< flag if back-tracking
  type(posnode_t) :: xo !< local variable
  type(posnode_t) :: xn !< local variable
  type(posnode_t) :: xm !< local variable
  real(kr8) :: domlen !< size of (periodic) domain in \f$ \xi \f$
  real(kr8) :: lenpath !< length of path in \f$ t \f$
  real(kr8) :: ztime !< time \f$ t \f$
  real(kr8) :: zdt !< time increment \f$ \Delta t \f$
  type(posvecl_t) :: zpos1   !< one position data
  type(posvecl_t) :: zpos2   !< one position data
  type(posveclis_t) :: rposl   !< list of position data
  real(kr4) :: phylenpath !< physical length of path
  logical :: lcoll   !< collision on path
  logical :: lexit   !< have left domain covered by HDS
  real(kr8) :: zpsim !< value of \f$ \psi \f$ at field line end (diagnostic)
  real(kr8) :: zmin!< min and max
  real(kr8) :: zmax !< min and max
  real(kr8), dimension(2) :: zk !< sector number
  integer(ki4), save :: firstcall  !< flag up first step of new trajectory

  inpow=powelt_addr(self,powcal%powres%npowe)
  !dbgwr write(*,*) "inpow=",inpow !dbgwr
  ! check whether looking at small number of tracks
  if (powcal%n%ntrack>0) then
     ! is this element track in the test set
     do j=1,powcal%n%ntrack
        if (powcal%n%trackno(j)==inpow) go to 1
     end do
     return
  end if
1     continue
  domlen=powcal%powres%beq%domlen(3)
  ! mark powelt as unknown (2) by default
  powcal%powres%pow(inpow)=2
  ! needed in shadowed case
  if (powcal%n%shadow>0.AND..NOT.allocated(rposl%pos)) then
     allocate(rposl%pos(1),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else if (powcal%n%shadow<-1) then
     return
  end if
  !DVEC      allocate(work1(powcal%odes%n%stepmax),stat=status) !DVEC
  !DVEC      call log_alloc_check(m_name,s_name,0,status) !DVEC

  ! calculate, save and test initial position
  call powelt_pos(self,powcal,zpos)
  ! Examine start trajectory behaviour in RZxi
  call powelt_normalm(self,powcal,znormalm)
  znormald=znormalm
  zposd=zpos
  zt0=zposd(3)
  ! start trajectory
  powcal%odes%ndt=1
  powcal%odes%t=zt0
  powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=zposd
  powcal%odes%posk(powcal%odes%ndt)=0
  !psi      call spl2d_eval(powcal%powres%beq%psi,zposd(1),zposd(2),zpsi)
  !psi! evaluate I aka f at psi
  !psi      call spleval(powcal%powres%beq%f,powcal%powres%beq%mr,powcal%powres%beq%psiaxis,&
  !psi      powcal%powres%beq%psiqbdry,zpsi,zf,1)
  ! set initial time step and guess first new point of track
  zdfaca(1:3)=powcal%powres%qfaca
  zdfaca(4:6)=powcal%powres%geobjl%quantfm%offvec
  call odes_1ststep1(powcal%odes,ierr,zdfaca,&
 &powcal%powres%beq%psi,powcal%powres%beq%rispldr,powcal%powres%beq%rispldz,0)
  if (ierr>0) return
  xpath=powcal%odes%vecp%pos(powcal%odes%ndt+1)%posvec-zposd
  zpdotnm=dot_product(xpath,znormald)

  ! Check behaviour in Cartesians
  call powelt_normal(self,powcal,znormal)
  znormald=znormal
  ! ensure path is outward from surface
  xpath=0
  zs=1
  do j=0,1
     ! dequantise
     zpos1%posvec=powcal%odes%vecp%pos(powcal%odes%ndt+j)%posvec
     zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
     ! convert xi to zeta
     zzeta=zpos2%posvec(3)/powcal%powres%beq%nzets
     zpos2%posvec(3)=zzeta
     !
     zposang%pos=zpos2%posvec
     zposang%opt=1 ; zposang%units=0
     call posang_tfm(zposang,-3)
     xpath=zposang%pos+zs*xpath
     zs=-zs
  end do
  zpdotn=dot_product(xpath,znormald)
  if (zpdotn*zpdotnm<0) then
     call log_error(m_name,s_name,3,error_warning,'Starting direction uncertain')
     write(*,*) inpow, zpdotn, zpdotnm
     ! default to shadowed
     powcal%powres%pow(inpow)=0
  end if
  ! adjust sign of timestep
  powcal%odes%dt=sign(1._kr8,zpdotn)*powcal%odes%dt
  ! scaling factor no longer carries timestep sign
  zdfaca(1:3)=powcal%powres%qfaca
  zdfaca(4:6)=powcal%powres%geobjl%quantfm%offvec
  if (powcal%n%shadow>0) then
     ! find start point in tree
     rposl%pos(1)=powcal%odes%vecp%pos(powcal%odes%ndt)
     rposl%np=1
     iobj%geobj=1
     iobj%objtyp=1 ! for point
     call btree_mfind(btree,iobj,rposl,inode)
     if (inode<0) then
        call log_error(m_name,s_name,4,error_warning,'cannot locate start point in tree')
        return
     end if
     xo%posvec=rposl%pos(1)%posvec
     xo%node=inode
  end if
  ibacktr=0
  if (loutdt) write(170,*)
  if (loutdt) write(170,*)
  if (loutdt) write(170,*) '#',self%ie,self%je

100   continue
  ! initialise loop for path
  firstcall=1
  lenpath=0
  ierr=0
  ! time step loop for path
  loop_path: do
     lcoll=.FALSE.
     nobjhit=0
     zpsim=0
     !     write(*,*) 'step=', powcal%odes%ndt

     ztime=powcal%odes%t

     zposd=xo%posvec
     call spl2d_eval(powcal%powres%beq%psi,zposd(1),zposd(2),zpsi)
     !psi! evaluate I aka f at psi
     !psi      call spleval(powcal%powres%beq%f,powcal%powres%beq%mr,powcal%powres%beq%psiaxis,&
     !psi      powcal%powres%beq%psiqbdry,zpsi,zf,1)
     call odes_1ststepcon(powcal%odes,ierr,zdfaca,firstcall,&
 &   powcal%powres%beq%psi,powcal%powres%beq%rispldr,powcal%powres%beq%rispldz)
     ! check for termination due to solution failure on this field-line
     if (ierr>0) exit

     zdt=powcal%odes%t-ztime
     !DVEC      work1(powcal%odes%ndt-1)=zdt  !DVEC
     !DVEC      write(*,*) 'tstep'  !DVEC
     lenpath=lenpath+abs(zdt)
     if (loutdt) write(170,*) zdt

     if (powcal%n%shadow>0) then
        ! check for collision or exit from computational domain
        xn%posvec=powcal%odes%vecp%pos(powcal%odes%ndt)%posvec
        zmin=minval(xn%posvec-powcal%powres%beq%n%xbb(:,1))
        zmax=maxval(xn%posvec-powcal%powres%beq%n%xbb(:,2))
        lexit=(zmin<0.OR.zmax>0)
        if (lexit) then
           nobjhit=-1
           call powelt_setpow(self,powcal,nobjhit,xo,inpow,gshadl)
           exit
        end if
        !xi         zxi=xn%posvec(3)
        !xi! check for leaving periodic cell
        !xi         if (zxi<powcal%powres%beq%n%ximin) then
        !xi! interpolate for xm=x(0) and check versus objects
        !xi            zf1=(powcal%powres%beq%n%ximin-xo%posvec(3))/(zxi-xo%posvec(3))
        !xi            xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
        !xi            xm%node=0
        !xi! only terminate if trajectory is long enough
        !xi            if (abs(lenpath)>powcal%odes%n%termcon(1)) then
        !xi               call pcle_movet(xo,xm,xo%node,gshadl,btree, &
        !xi     &         powcal%n%termp,nobjhit)
        !xi               lcoll=(nobjhit/=0)
        !xi            end if
        !xi            if (lcoll) then
        !xi               powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
        !xi               call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
        !xi               exit
        !xi            end if
        !xi            powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)-1
        !xi            xo=xm
        !xi            xo%posvec(3)=powcal%powres%beq%n%ximax
        !xi            xo%node=0
        !xi            xn%posvec(3)=min(powcal%powres%beq%ximaxm, xn%posvec(3)+domlen)
        !xi            powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
        !xi
        !xi         else if (zxi>powcal%powres%beq%n%ximax) then
        !xi! interpolate for xm=x(2pi) and check versus objects
        !xi            zf1=(powcal%powres%beq%n%ximax-xo%posvec(3))/(zxi-xo%posvec(3))
        !xi            xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
        !xi            xm%node=0
        !xi! only terminate if trajectory is long enough
        !xi            if (abs(lenpath)>powcal%odes%n%termcon(1)) then
        !xi               call pcle_movet(xo,xm,xo%node,gshadl,btree, &
        !xi     &         powcal%n%termp,nobjhit)
        !xi               lcoll=(nobjhit/=0)
        !xi            end if
        !xi            if (lcoll) then
        !xi               powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
        !xi               call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
        !xi! record psi at boundary if hit termplane boundary
        !xi               if (nobjhit<=GEOBJ_OFFSET) zpsim=powcal%powres%geobjl%obj(inpow)%weight
        !xi               exit
        !xi            end if
        !xi            powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)+1
        !xi            xo=xm
        !xi            xo%posvec(3)=powcal%powres%beq%n%ximin
        !xi            xo%node=0
        !xi            xn%posvec(3)=max(powcal%powres%beq%ximinp, xn%posvec(3)-domlen)
        !xi            powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
        !xi         end if

        ! find new node and position if hits boundary
        xn%node=0
        ! only terminate if trajectory is long enough
        if (abs(lenpath)>powcal%odes%n%termcon(1)) then
           call pcle_movet(xo,xn,xo%node,gshadl,btree, &
 &         powcal%n%termp,nobjhit)
           lcoll=(nobjhit/=0)
        end if
        if (lcoll) then
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
           call powelt_setpow(self,powcal,nobjhit,xo,inpow,gshadl)
           !! record psi at boundary if hit termplane boundary
           !           if (nobjhit<=GEOBJ_OFFSET) zpsim=powcal%powres%geobjl%obj(inpow)%weight
           !D           write(*,*) 'Final posn',xn%posvec! writediagn
           exit
        end if
        ! check whether hit midplane
        !mid         if (ilpmidplane) then
        !mid            if ((xn%posvec(2)-powcal%odes%n%termcon(3))*powcal%odes%n%ntermcon>0) then
        !mid               lcoll=.TRUE.
        !mid! interpolate for xm=x(midplane)
        !mid               zf1=(powcal%odes%n%termcon(3)-xo%posvec(2))/(xn%posvec(2)-xo%posvec(2))
        !mid! fix up in case zf1 garbage caused by xi boundary interaction
        !mid               zf1=max(min(zf1,1.),0.)
        !mid               xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
        !mid               xm%node=0
        !mid               nobjhit=-2
        !mid               call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
        !mid               zpsim=powcal%powres%geobjl%obj(inpow)%weight
        !mid               xn=xm
        !mid               exit
        !mid            end if
        !mid         end if
        ! save new data
        xo=xn
        ! check for near-collision and repeat step (once) if necessary TO DO???
     end if

     !DP !! psi diagnostic !DP
     !DP      zr=xo%posvec(1) !DP
     !DP      zz=xo%posvec(2) !DP
     !DP      call spl2d_eval(powcal%powres%beq%psi,zr,zz,zpsi) !DP
     !DP! dequantise xo !DP
     !DP      zpos1%posvec=xo%posvec !DP
     !DP      zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm) !DP
     !DP     write(*,*) lenpath,zpos2%posvec(1),zpos2%posvec(2),zpsi !DP
     ! check for termination due to overlong trajectory (number of steps)
     indt=powcal%odes%ndt
     if (indt>=powcal%odes%n%stepmax) exit
     ! or length of path in t
     ! (often covered when test exit from domain)
     if (abs(lenpath)>abs(powcal%odes%n%tmax-zt0)) exit
  end do loop_path


  ip=powcal%odes%ndt
  powcal%odes%vecp%np=ip
  !
  if (ibacktr==ipback.AND.powcal%powres%flinm) then
     ! open file to record RZxi
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' lenpath= '',1pg12.5,'' objhit= '',I8)') &
 &   self%ie,self%je,lenpath,nobjhit
     write(icfile,'(''trackm'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     ! write track in RZxi coordinates
     ! de-quantise
     call position_invqtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     do l=1,ip
        powcal%odes%vecp%pos(l)%posvec(3)=&
 &      powcal%odes%vecp%pos(l)%posvec(3)+2*const_pid*powcal%odes%posk(l)
     end do
     call position_writelis(powcal%odes%vecp,'track',nplot)
     do l=1,ip
        powcal%odes%vecp%pos(l)%posvec(3)=&
 &      powcal%odes%vecp%pos(l)%posvec(3)-2*const_pid*powcal%odes%posk(l)
     end do
     ! re-quantise
     call position_qtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     ! close file
     close(nplot)
  end if

  if (ibacktr==ipback.AND.powcal%powres%flinx) then
     ! write track in Cartesian coordinates
     ! first convert track
     allocate(wposl%pos(ip), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     do l=1,ip
        ! dequantise
        zpos1%posvec=powcal%odes%vecp%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*powcal%odes%posk(l))&
 &      /powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=ip
     call position_lenlis(wposl,phylenpath)
     !DVEC         do l=1,ip-1 !DVEC
     !DVEC         zr=sqrt(wposl%pos(l)%posvec(1)**2+wposl%pos(l)%posvec(2)**2)/1000 !DVEC
     !DVEC         wposl%pos(l)%posvec=(wposl%pos(l+1)%posvec-wposl%pos(l)%posvec)/& !DVEC
     !DVEC      &  (zr*work1(l)) !DVEC
     !DVEC         end do !DVEC
     !DVEC         wposl%pos(ip)=wposl%pos(ip-1) !DVEC
     !DVEC! write track (or field vectors if !DVEC lines executed)

     ! open file to record Cartesian track and write
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' lenpath= '',1pg12.5,'' objhit= '',I8)') &
 &   self%ie,self%je,phylenpath,nobjhit
     !        write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' trackx'')') self%ie,self%je
     write(icfile,'(''trackx'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     call position_writelis(wposl,'track',nplot)
     !DIAG!   dump end point !DIAG
     !DIAG!        call position_writev(wposl%pos(ip),88) !DIAG
     ! close file and deallocate
     deallocate(wposl%pos)
     !DVEC         deallocate(work1)  !DVEC
     close(nplot)
  end if

  if (nobjhit<=GEOBJ_OFFSET.AND.ibacktr==ipback.AND.powcal%powres%flinends) then
     ! write first and last points in Cartesian coordinates
     ! first convert track points
     allocate(wposl%pos(3), stat=status)
     call log_alloc_check(m_name,s_name,11,status)
     wposl%pos(1)%posvec=powcal%odes%vecp%pos(1)%posvec
     zk(1)=powcal%odes%posk(1)
     wposl%pos(2)%posvec=xm%posvec
     zk(2)=powcal%odes%posk(ip)
     do l=1,2
        ! dequantise
        zpos1%posvec=wposl%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*zk(l))/powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=3
     ! open file to record Cartesian track and write
     wposl%pos(3)%posvec(1)=self%ie
     wposl%pos(3)%posvec(2)=self%je
     wposl%pos(3)%posvec(3)=zpsim
     do l=1,3
        call position_writev(wposl%pos(l),powcal%powres%nflends)
     end do
     ! deallocate
     deallocate(wposl%pos)
  end if

  if (ibacktr<ipback) then
     ibacktr=ipback
     !        zdfaca(1:3)=-zdfaca(1:3)
     powcal%odes%ndt=powcal%odes%ndt-1
     powcal%odes%n%stepmax=2*powcal%odes%n%stepmax
     powcal%odes%dt=-powcal%odes%dt
     goto 100
  end if

end subroutine powelt_move0
!---------------------------------------------------------------------
!> motion along fieldline global, axisymm, cylindricals
subroutine powelt_move1(self,powcal,gshadl,btree)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure
  type(btree_t), intent(inout),optional :: btree !< btree data

  !! local
  character(*), parameter :: s_name='powelt_move1' !< subroutine name
  integer(ki4), parameter :: ipback=0 !< back-track test if unity
  real(kr4), dimension(3) :: zpos !< start position vector
  real(kr8), dimension(3) :: zposd !< start position vector
  real(kr8), dimension(3) :: xpath !< start path position vector Cartesians
  !psi real(kr8) :: zpsi !< value of \f$ \psi \f$
  real(kr8) :: zzeta !< value of \f$ \zeta \f$
  real(kr8) :: zt0 !< start value of \f$ t \f$
  type(posang_t) :: zposang   !< posang data structure
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4), dimension(3) :: znormalm !< unit normal vector in mapped coords
  real(kr8), dimension(3) :: znormald !< unit normal vector
  real(kr8), dimension(6) :: zdfaca !< direction factor and offsets
  real(kr8) :: zs !< \f$ \pm 1 \f$
  real(kr8) :: zpdotn !< path dotted with normal vector
  real(kr8) :: zpdotnm !< path dotted with normal vector in mapped coords
  integer(ki4) :: ierr !< error flag
  integer(ki4) :: ip !< entries in track
  type(geobj_t) :: iobj !< geo object
  integer(ki4) :: inode  !< local variable
  integer(ki4) :: indt !< shorthand for current timestep number
  !> Non-zero at fieldline termination
  !! - (-2) Line hit terminating plane
  !! - (-1) Line left octree
  !! - Positive - number of triangle hit
  integer(ki4) :: nobjhit=0  !< .
  integer(ki4) :: ibacktr !< flag if back-tracking
  type(posnode_t) :: xo !< local variable
  type(posnode_t) :: xn !< local variable
  type(posnode_t) :: xm !< local variable
  real(kr8) :: zf1 !< fraction of path length
  real(kr8) :: zxi !< value of \f$ \xi \f$
  real(kr8) :: domlen !< size of (periodic) domain in \f$ \xi \f$
  real(kr8) :: lenpath !< length of path in \f$ t \f$
  real(kr8) :: ztime !< time \f$ t \f$
  real(kr8) :: zdt !< time increment \f$ \Delta t \f$
  type(posvecl_t) :: zpos1   !< one position data
  type(posvecl_t) :: zpos2   !< one position data
  type(posveclis_t), save :: rposl   !< list of position data
  real(kr4) :: phylenpath !< physical length of path
  logical :: lpset   !< true if flux at fieldline start-point has been set
  logical :: linpfr   !< true if fieldline start-point in private flux region
  logical :: lcoll   !< collision on path
  real(kr8) :: zpsim !< value of \f$ \psi \f$ at field line end (diagnostic)
  real(kr8) :: zpsi !<  \f$ \psi \f$ debugging diagnostic !DP
  real(kr8) :: zr !<  \f$ \r \f$ debugging diagnostic !DP
  real(kr8) :: zz !<  \f$ \z \f$ debugging diagnostic !DP
  real(kr8), dimension(2) :: zk !< sector number
  integer(ki4), save :: firstcall  !< flag up first step of new trajectory
  integer(ki4) :: inou     !< whether on HFS (in, 1) or LFS (ou, 2) of centre
  real(kr8) :: distrack     !< signed distance from central track
  integer(ki4) :: itp   !< indexes psi limits array
  integer(ki4) :: ityps   !< test versus lower (1) or upper (2) skylight type
  integer(ki4) :: irid   !< \f$ \pm 1 \f$ depending whether below (+) or above (-) centre
  integer(ki4) :: interm   !< index of dynamical entry in termplane
  real(kr8), dimension(2) :: zrz    !<  Query point position \f$ (R,Z) \f$
  real(kr8), dimension(2) :: zlts !< limits
  real(kr8) :: zskyl !<  \f$ \z \f$ value for skylight

  inpow=powelt_addr(self,powcal%powres%npowe)
  !dbgwr write(*,*) "inpow=",inpow !dbgwr
  ! check whether looking at small number of tracks
  if (powcal%n%ntrack>0) then
     ! is this element track in the test set
     do j=1,powcal%n%ntrack
        if (powcal%n%trackno(j)==inpow) go to 1
     end do
     return
  end if
1     continue
  domlen=powcal%powres%beq%domlen(3)
  ! ignore midplane if termination planes present
  if (powcal%n%ltermplane) then
     ilpmidplane=.FALSE.
     powcal%n%termp%termstore(1,:)=(/1.e30,1.e30,1.e30/)
  else
     ilpmidplane=lpmidplane
  end if
  ! mark powelt as unknown (2) by default
  powcal%powres%pow(inpow)=2
  ! needed in shadowed case
  if (powcal%n%shadow>0.AND..NOT.allocated(rposl%pos)) then
     allocate(rposl%pos(1),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else if (powcal%n%shadow<-1) then
     return
  end if
  !DVEC      allocate(work1(powcal%odes%n%stepmax),stat=status) !DVEC
  !DVEC      call log_alloc_check(m_name,s_name,0,status) !DVEC

  ! calculate, save and test initial position
  call powelt_pos(self,powcal,zpos)
  ! Examine start trajectory behaviour in RZxi
  call powelt_normalm(self,powcal,znormalm)
  znormald=znormalm
  zposd=zpos
  zt0=zposd(3)
  lpset=.FALSE.
  if (.NOT.powcal%n%lpfpower) then
     ! if no power in private flux region, return if in private flux region
     call spl2d_eval(powcal%powres%beq%psi,zposd(1),zposd(2),zpsi)
     linpfr= ( beq_rsig()*(powcal%powres%beq%psibdry-zpsi)>0 )
     if (linpfr) then
        powcal%powres%pow(inpow)=0
        return
     end if
     lpset=.TRUE.
  end if
  ! start trajectory
  powcal%odes%ndt=1
  powcal%odes%t=zt0
  powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=zposd
  powcal%odes%posk(powcal%odes%ndt)=0
  !psi      call spl2d_eval(powcal%powres%beq%psi,zposd(1),zposd(2),zpsi)
  !psi! evaluate I aka f at psi
  !psi      call spleval(powcal%powres%beq%f,powcal%powres%beq%mr,powcal%powres%beq%psiaxis,&
  !psi      powcal%powres%beq%psiqbdry,zpsi,zf,1)
  ! set initial time step and guess first new point of track
  zdfaca(1:3)=powcal%powres%qfaca
  zdfaca(4:6)=powcal%powres%geobjl%quantfm%offvec
  call odes_1ststep1(powcal%odes,ierr,zdfaca,&
 &powcal%powres%beq%psi,powcal%powres%beq%rispldr,powcal%powres%beq%rispldz,0)
  if (ierr>0) return
  xpath=powcal%odes%vecp%pos(powcal%odes%ndt+1)%posvec-zposd
  zpdotnm=dot_product(xpath,znormald)

  ! Check behaviour in Cartesians
  call powelt_normal(self,powcal,znormal)
  znormald=znormal
  ! ensure path is outward from surface
  xpath=0
  zs=1
  do j=0,1
     ! dequantise
     zpos1%posvec=powcal%odes%vecp%pos(powcal%odes%ndt+j)%posvec
     zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
     ! convert xi to zeta
     zzeta=zpos2%posvec(3)/powcal%powres%beq%nzets
     zpos2%posvec(3)=zzeta
     !
     zposang%pos=zpos2%posvec
     zposang%opt=1 ; zposang%units=0
     call posang_tfm(zposang,-3)
     xpath=zposang%pos+zs*xpath
     zs=-zs
  end do
  zpdotn=dot_product(xpath,znormald)
  if (zpdotn*zpdotnm<0) then
     call log_error(m_name,s_name,3,error_warning,'Starting direction uncertain')
     write(*,*) inpow, zpdotn, zpdotnm
     ! default to shadowed
     powcal%powres%pow(inpow)=0
  end if
  ! default number of active termplanes
  interm=powcal%n%termp%ntermplane
  if (powcal%powres%beq%n%skylpsi) then
     if (.NOT.lpset) then
        call spl2d_eval(powcal%powres%beq%psi,zposd(1),zposd(2),zpsi)
     end if
     linpfr=( beq_rsig()*(powcal%powres%beq%psibdry-zpsi)>0 )
     !! optionally skip flux skylight test if in private flux region
     !if (linpfr.AND..NOT.powcal%n%lskylpfr) go to 2
     ! always skip skylight test if in private flux region
     if (.NOT.linpfr) then
        ityps=0
        if (zposd(2)<powcal%powres%beq%zcenq) then
           ! lower skylight active
           if (powcal%powres%skyl%n%skyltyp<3) ityps=1
        else
           ! upper skylight only active
           ityps=4-powcal%powres%skyl%n%skyltyp
        end if
        if (ityps==1.OR.ityps==2) then
           zrz=(/zposd(1),zposd(2)/)
           call beq_ctrackcq(powcal%powres%beq,powcal%powres%beq%ctrackrz,&
           powcal%powres%beq%nctrack,zrz,distrack)
           inou=2+sign(0.5_kr8,distrack) ! 1 if distrack<0, 2 if distrack>0
           ! check within allowable range of psi
           zlts=powcal%powres%skyl%psilts(1:2,inou,ityps)
           if ( (zpsi-zlts(1))*(zpsi-zlts(2)) <= 0 ) then
              itp=1+int((zpsi-zlts(1))/powcal%powres%skyl%psidelta(inou,ityps))
              itp=min( itp, powcal%powres%skyl%dimbox(inou,ityps) )
              if (inou==1) zskyl=powcal%powres%skyl%inboxz(itp,ityps)
              if (inou==2) zskyl=powcal%powres%skyl%ouboxz(itp,ityps)
              if(debug)write(800,*) "psi, Z values", zpsi, zskyl,zposd(2)
              if (abs(zskyl-zposd(2))<powcal%powres%skyl%n%toli(2)) then
                 if (debug) then
                    nskyl1=nskyl1+1
                    write(800,*) "coincident skylight", nskyl1
                 end if
                 return
              else
                 irid=3-2*ityps
                 if ((zposd(2)-zskyl)*irid>0) then
                    if (debug) then
                       nskyl2=nskyl2+1
                       write(800,*) "above skylight", nskyl2
                    end if
                    return
                 end if
              end if
              ! set dynamic termplane based on skylight
              interm=powcal%n%termp%ntermplane+1
              powcal%n%termp%termplane(interm,1)=zskyl
              powcal%n%termp%termplanedir(interm,2)=irid
           end if
        end if
     end if
  end if
  !2     continue
  ! set number of active termplanes
  powcal%n%termp%ntermactive=interm
  ! adjust sign of timestep
  powcal%odes%dt=sign(1._kr8,zpdotn)*powcal%odes%dt
  ! scaling factor no longer carries timestep sign
  zdfaca(1:3)=powcal%powres%qfaca
  zdfaca(4:6)=powcal%powres%geobjl%quantfm%offvec
  if (powcal%n%shadow>0) then
     ! find start point in tree
     rposl%pos(1)=powcal%odes%vecp%pos(powcal%odes%ndt)
     rposl%np=1
     iobj%geobj=1
     iobj%objtyp=1 ! for point
     call btree_mfind(btree,iobj,rposl,inode)
     if (inode<0) then
        call log_error(m_name,s_name,4,error_warning,'cannot locate start point in tree')
        return
     end if
     xo%posvec=rposl%pos(1)%posvec
     xo%node=inode
  end if
  ibacktr=0
  if (loutdt) write(170,*)
  if (loutdt) write(170,*)
  if (loutdt) write(170,*) '#',self%ie,self%je

100   continue
  ! initialise loop for path
  lenpath=0
  firstcall=1
  ierr=0
  !     powcal%odes%near=0
  ! time step loop for path
  loop_path: do
     lcoll=.FALSE.
     nobjhit=0
     zpsim=0
     if (loutdt) write(*,*) 'step=', powcal%odes%ndt

     ztime=powcal%odes%t

     !psi zposd=xo%posvec
     !psi call spl2d_eval(powcal%powres%beq%psi,zposd(1),zposd(2),zpsi)
     !psi! evaluate I aka f at psi
     !psi      call spleval(powcal%powres%beq%f,powcal%powres%beq%mr,powcal%powres%beq%psiaxis,&
     !psi      powcal%powres%beq%psiqbdry,zpsi,zf,1)
     call odes_1ststepcon(powcal%odes,ierr,zdfaca,firstcall,&
 &   powcal%powres%beq%psi,powcal%powres%beq%rispldr,powcal%powres%beq%rispldz)
     ! check for termination due to solution failure on this field-line
     if (ierr>0) exit

     zdt=powcal%odes%t-ztime
     !DVEC      work1(powcal%odes%ndt-1)=zdt  !DVEC
     !DVEC      write(*,*) 'tstep'  !DVEC
     lenpath=lenpath+abs(zdt)
     if (loutdt) write(170,*) 'tstep',zdt

     if (powcal%n%shadow>0) then
        ! check for collision or exit from computational domain
        xn%posvec=powcal%odes%vecp%pos(powcal%odes%ndt)%posvec
        zxi=powcal%odes%t
        ! check for leaving periodic cell
        if (zxi<powcal%powres%beq%n%ximin) then
           ! interpolate for xm=x(0) and check versus objects
           zf1=(powcal%powres%beq%n%ximin-ztime)/(zxi-ztime)
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           xm%node=0
           ! only terminate if trajectory is long enough
           if (abs(lenpath)>powcal%odes%n%termcon(1)) then
              call pcle_movet(xo,xm,xo%node,gshadl,btree, &
 &            powcal%n%termp,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) then
              powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              exit
           end if
           powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)-1
           xo=xm
           xo%posvec(3)=powcal%powres%beq%n%ximax
           xo%node=0
           powcal%odes%t=min(powcal%powres%beq%ximaxm, powcal%odes%t+domlen)
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec(3)=powcal%odes%t
           xn%posvec(3)=powcal%odes%t

        else if (zxi>powcal%powres%beq%n%ximax) then
           ! interpolate for xm=x(2pi) and check versus objects
           zf1=(powcal%powres%beq%n%ximax-ztime)/(zxi-ztime)
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           xm%node=0
           ! only terminate if trajectory is long enough
           if (abs(lenpath)>powcal%odes%n%termcon(1)) then
              call pcle_movet(xo,xm,xo%node,gshadl,btree, &
 &            powcal%n%termp,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) then
              powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              ! record psi at boundary if hit termplane boundary
              if (nobjhit<=GEOBJ_OFFSET) zpsim=powcal%powres%geobjl%obj(inpow)%weight
              exit
           end if
           powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)+1
           xo=xm
           xo%posvec(3)=powcal%powres%beq%n%ximin
           xo%node=0
           powcal%odes%t=max(powcal%powres%beq%ximinp, powcal%odes%t-domlen)
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec(3)=powcal%odes%t
           xn%posvec(3)=powcal%odes%t
        end if

        ! find new node and position if hits boundary
        xn%node=0
        ! only terminate if trajectory is long enough
        if (abs(lenpath)>powcal%odes%n%termcon(1)) then
           call pcle_movet(xo,xn,xo%node,gshadl,btree, &
 &         powcal%n%termp,nobjhit)
           lcoll=(nobjhit/=0)
        end if
        if (lcoll) then
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
           call powelt_setpow(self,powcal,nobjhit,xn,inpow,gshadl)
           ! record psi at boundary if hit termplane boundary
           if (nobjhit<=GEOBJ_OFFSET) zpsim=powcal%powres%geobjl%obj(inpow)%weight
           !D           write(*,*) 'Final posn',xn%posvec! writediagn
           exit
        end if
        ! save new data
        xo=xn
        ! check for near-collision and repeat step (once) if necessary TO DO???
     end if

     !! psi diagnostic !DP
     if (loutdt) zr=xo%posvec(1) !DP
     if (loutdt) zz=xo%posvec(2) !DP
     if (loutdt) call spl2d_eval(powcal%powres%beq%psi,zr,zz,zpsi) !DP
     !DP! dequantise xo !DP
     if (loutdt) zpos1%posvec=xo%posvec !DP
     if (loutdt) zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm) !DP
     if (loutdt) write(*,*) lenpath,zpos2%posvec(1),zpos2%posvec(2),zpsi !DP
     ! check for termination due to overlong trajectory (number of steps)
     indt=powcal%odes%ndt
     if (indt>=powcal%odes%n%stepmax) exit
     ! or length of path in t
     ! (often covered when test exit from domain)
     if (abs(lenpath)>abs(powcal%odes%n%tmax-zt0)) exit
  end do loop_path


  ip=powcal%odes%ndt
  powcal%odes%vecp%np=ip
  !
  if (ibacktr==ipback.AND.powcal%powres%flinm) then
     ! open file to record RZxi
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' lenpath= '',1pg12.5,'' objhit= '',I8)') &
 &   self%ie,self%je,lenpath,nobjhit
     write(icfile,'(''trackm'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     ! write track in RZxi coordinates
     ! de-quantise
     call position_invqtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     do l=1,ip
        powcal%odes%vecp%pos(l)%posvec(3)=&
 &      powcal%odes%vecp%pos(l)%posvec(3)+2*const_pid*powcal%odes%posk(l)
     end do
     call position_writelis(powcal%odes%vecp,'track',nplot)
     do l=1,ip
        powcal%odes%vecp%pos(l)%posvec(3)=&
 &      powcal%odes%vecp%pos(l)%posvec(3)-2*const_pid*powcal%odes%posk(l)
     end do
     ! re-quantise
     call position_qtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     ! close file
     close(nplot)
  end if

  if (ibacktr==ipback.AND.powcal%powres%flinx) then
     ! write track in Cartesian coordinates
     ! first convert track
     allocate(wposl%pos(ip), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     do l=1,ip
        ! dequantise
        zpos1%posvec=powcal%odes%vecp%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*powcal%odes%posk(l))&
 &      /powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=ip
     call position_lenlis(wposl,phylenpath)
     !DVEC         do l=1,ip-1 !DVEC
     !DVEC         zr=sqrt(wposl%pos(l)%posvec(1)**2+wposl%pos(l)%posvec(2)**2)/1000 !DVEC
     !DVEC         wposl%pos(l)%posvec=(wposl%pos(l+1)%posvec-wposl%pos(l)%posvec)/& !DVEC
     !DVEC      &  (zr*work1(l)) !DVEC
     !DVEC         end do !DVEC
     !DVEC         wposl%pos(ip)=wposl%pos(ip-1) !DVEC
     !DVEC! write track (or field vectors if !DVEC lines executed)

     ! open file to record Cartesian track and write
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' lenpath= '',1pg12.5,'' objhit= '',I8)') &
 &   self%ie,self%je,phylenpath,nobjhit
     !        write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' trackx'')') self%ie,self%je
     write(icfile,'(''trackx'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     call position_writelis(wposl,'track',nplot)
     !DIAG!   dump end point !DIAG
     !DIAG!        call position_writev(wposl%pos(ip),88) !DIAG
     ! close file and deallocate
     deallocate(wposl%pos)
     !DVEC         deallocate(work1)  !DVEC
     close(nplot)
  end if

  if (nobjhit<=GEOBJ_OFFSET.AND.ibacktr==ipback.AND.powcal%powres%flinends) then
     ! write first and last points in Cartesian coordinates
     ! first convert track points
     allocate(wposl%pos(3), stat=status)
     call log_alloc_check(m_name,s_name,11,status)
     wposl%pos(1)%posvec=powcal%odes%vecp%pos(1)%posvec
     zk(1)=powcal%odes%posk(1)
     wposl%pos(2)%posvec=xm%posvec
     zk(2)=powcal%odes%posk(ip)
     do l=1,2
        ! dequantise
        zpos1%posvec=wposl%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*zk(l))/powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=3
     ! open file to record Cartesian track and write
     wposl%pos(3)%posvec(1)=self%ie
     wposl%pos(3)%posvec(2)=self%je
     wposl%pos(3)%posvec(3)=zpsim
     do l=1,3
        call position_writev(wposl%pos(l),powcal%powres%nflends)
     end do
     ! deallocate
     deallocate(wposl%pos)
  end if

  if (ibacktr<ipback) then
     ibacktr=ipback
     !        zdfaca(1:3)=-zdfaca(1:3)
     powcal%odes%ndt=powcal%odes%ndt-1
     powcal%odes%n%stepmax=2*powcal%odes%n%stepmax
     powcal%odes%dt=-powcal%odes%dt
     goto 100
  end if

end subroutine powelt_move1
!---------------------------------------------------------------------
!> motion along fieldline global, analytic ripple
subroutine powelt_move2(self,powcal,gshadl,btree)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure
  type(btree_t), intent(inout),optional :: btree !< btree data

  !! local
  character(*), parameter :: s_name='powelt_move2' !< subroutine name
  integer(ki4), parameter :: ipback=0 !< back-track test if unity
  real(kr4), dimension(3) :: zpos !< start position vector
  real(kr8), dimension(3) :: zposd !< start position vector
  real(kr8), dimension(3) :: xpath !< start path position vector Cartesians
  real(kr8) :: zzeta !< value of \f$ \zeta \f$
  real(kr8) :: zt0 !< start value of \f$ t \f$
  type(posang_t) :: zposang   !< posang data structure
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4), dimension(3) :: znormalm !< unit normal vector in mapped coords
  real(kr8), dimension(3) :: znormald !< unit normal vector
  real(kr8), dimension(3) :: zdfaca !< direction factor
  real(kr8) :: zs !< \f$ \pm 1 \f$
  real(kr8) :: zpdotn !< path dotted with normal vector
  real(kr8) :: zpdotnm !< path dotted with normal vector in mapped coords
  integer(ki4) :: ierr !< error flag
  integer(ki4) :: ip !< entries in track
  type(geobj_t) :: iobj !< geo object
  integer(ki4) :: inode  !< local variable
  integer(ki4) :: indt !< shorthand for current timestep number
  !> Non-zero at fieldline termination
  !! - (-2) Line hit terminating plane
  !! - (-1) Line left octree
  !! - Positive - number of triangle hit
  integer(ki4) :: nobjhit=0  !< .
  integer(ki4) :: ibacktr !< flag if back-tracking
  type(posnode_t) :: xo !< local variable
  type(posnode_t) :: xn !< local variable
  type(posnode_t) :: xm !< local variable
  real(kr8) :: zf1 !< fraction
  real(kr8) :: zt !< value of \f$ t \f$
  real(kr8) :: domlen !< size of (periodic) domain in \f$ t \f$
  real(kr8) :: lenpath !< length of path in \f$ t \f$
  type(posvecl_t) :: zpos1   !< one position data
  type(posvecl_t) :: zpos2   !< one position data
  type(posveclis_t) :: rposl   !< list of position data
  logical :: lcoll   !< collision on path

  inpow=powelt_addr(self,powcal%powres%npowe)
  ! check whether looking at small number of tracks
  if (powcal%n%ntrack>0) then
     ! is this element track in the test set
     do j=1,powcal%n%ntrack
        if (powcal%n%trackno(j)==inpow) go to 1
     end do
     return
  end if
1     continue
  domlen=powcal%powres%beq%n%ximax-powcal%powres%beq%n%ximin
  ! mark powelt as unshadowed by default
  powcal%powres%pow(inpow)=1
  ! needed in shadowed case
  if (powcal%n%shadow>0.AND..NOT.allocated(rposl%pos)) then
     allocate(rposl%pos(1),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else if (powcal%n%shadow<-1) then
     return
  end if

  ! calculate, save and test initial position
  call powelt_pos(self,powcal,zpos)
  ! Examine start trajectory behaviour in RZxi
  call powelt_normalm(self,powcal,znormalm)
  znormald=znormalm
  zposd=zpos
  zt0=zposd(3)
  ! start trajectory
  powcal%odes%ndt=1
  powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=zposd
  powcal%odes%posk(powcal%odes%ndt)=0
  ! set initial time step and guess first new point of track
  zdfaca=powcal%powres%qfaca
  call odes_2ndstep1(powcal%odes,ierr,zdfaca,beq_ripple_h1,&
 &powcal%powres%beq%psi,powcal%powres%beq%dpsidr,powcal%powres%beq%dpsidz,0)
  if (ierr>0) return
  xpath=powcal%odes%vecp%pos(powcal%odes%ndt+1)%posvec-zposd
  zpdotnm=dot_product(xpath,znormald)

  ! Check behaviour in Cartesians
  call powelt_normal(self,powcal,znormal)
  znormald=znormal
  ! ensure path is outward from surface
  xpath=0
  zs=1
  do j=0,1
     ! dequantise
     zpos1%posvec=powcal%odes%vecp%pos(powcal%odes%ndt+j)%posvec
     zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
     ! convert xi to zeta
     zzeta=zpos2%posvec(3)/powcal%powres%beq%nzets
     zpos2%posvec(3)=zzeta
     !
     zposang%pos=zpos2%posvec
     zposang%opt=1 ; zposang%units=0
     call posang_tfm(zposang,-3)
     xpath=zposang%pos+zs*xpath
     zs=-zs
  end do
  zpdotn=dot_product(xpath,znormald)
  if (zpdotn*zpdotnm<0) then
     call log_error(m_name,s_name,3,error_warning,'Starting direction uncertain')
     write(*,*) inpow, zpdotn, zpdotnm
  end if
  ! adjust sign of timestep
  powcal%odes%dt=sign(1._kr8,zpdotn)*powcal%odes%dt
  ! scaling factor no longer carries timestep sign
  zdfaca=powcal%powres%qfaca
  if (powcal%n%shadow>0) then
     ! find start point in tree
     rposl%pos(1)=powcal%odes%vecp%pos(powcal%odes%ndt)
     rposl%np=1
     iobj%geobj=1
     iobj%objtyp=1 ! for point
     call btree_mfind(btree,iobj,rposl,inode)
     if (inode<0) then
        call log_error(m_name,s_name,4,error_warning,'cannot locate start point in tree')
        return
     end if
     xo%posvec=rposl%pos(1)%posvec
     xo%node=inode
  end if
  ibacktr=0

100   continue
  ! time step loop for path
  lenpath=0
  ierr=0
  !     powcal%odes%near=0
  loop_path: do
     lcoll=.FALSE.
     nobjhit=0
     !     write(*,*) 'step=', powcal%odes%ndt

     lenpath=lenpath+abs(powcal%odes%dt)

     call odes_2ndstep(powcal%odes,ierr,zdfaca,beq_ripple_h1,&
 &   powcal%powres%beq%psi,powcal%powres%beq%dpsidr,powcal%powres%beq%dpsidz)
     ! check for termination due to solution failure on this field-line
     if (ierr>0) exit

     if (powcal%n%shadow>0) then
        ! check for collision or exit from computational domain
        xn%posvec=powcal%odes%vecp%pos(powcal%odes%ndt)%posvec
        zt=xn%posvec(3)
        ! check for leaving periodic cell
        if (zt<powcal%powres%beq%n%ximin) then
           ! interpolate for xm=x(0) and check versus objects
           zf1=(powcal%powres%beq%n%ximin-xo%posvec(3))/(zt-xo%posvec(3))
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           xm%node=0
           ! only terminate if trajectory is long enough
           if (abs(lenpath)>powcal%odes%n%termcon(1)) then
              call pcle_movet(xo,xm,xo%node,gshadl,btree, &
 &            powcal%n%termp,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) then
              powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
              ! mark as really a collision, not just due to leaving volume
              if (nobjhit>0)  powcal%powres%pow(inpow)=0
              exit
           end if
           powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)-1
           xo=xm
           xo%posvec(3)=powcal%powres%beq%n%ximax
           xo%node=0
           xn%posvec(3)=min(powcal%powres%beq%ximaxm, xn%posvec(3)+domlen)
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec

        else if (zt>powcal%powres%beq%n%ximax) then
           ! interpolate for xm=x(2pi) and check versus objects
           zf1=(powcal%powres%beq%n%ximax-xo%posvec(3))/(zt-xo%posvec(3))
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           xm%node=0
           ! only terminate if trajectory is long enough
           if (abs(lenpath)>powcal%odes%n%termcon(1)) then
              call pcle_movet(xo,xm,xo%node,gshadl,btree, &
 &            powcal%n%termp,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) then
              powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
              ! mark as really a collision, not just due to leaving volume
              if (nobjhit>0)  powcal%powres%pow(inpow)=0
              exit
           end if
           powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)+1
           xo=xm
           xo%posvec(3)=powcal%powres%beq%n%ximin
           xo%node=0
           xn%posvec(3)=max(powcal%powres%beq%ximinp, xn%posvec(3)-domlen)
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
        end if

        ! find new node and position if hits boundary
        xn%node=0
        ! only terminate if trajectory is long enough
        if (abs(lenpath)>powcal%odes%n%termcon(1)) then
           call pcle_movet(xo,xn,xo%node,gshadl,btree, &
 &         powcal%n%termp,nobjhit)
           lcoll=(nobjhit/=0)
        end if
        if (lcoll) then
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
           ! mark as really a collision, not just due to leaving volume
           if (nobjhit>0)  powcal%powres%pow(inpow)=0
           !D           write(*,*) 'Final posn',xn%posvec! writediagn
           exit
        end if
        ! save new data
        xo=xn
        ! check for near-collision and repeat step (once) if necessary TO DO???
     end if

     ! check for termination due to overlong trajectory (number of steps)
     indt=powcal%odes%ndt
     if (indt>=powcal%odes%n%stepmax) exit
     ! or length of path in t
     ! (often covered when test exit from domain)
     if (abs(lenpath)>abs(powcal%odes%n%tmax-zt0)) exit
  end do loop_path


  ip=powcal%odes%ndt
  powcal%odes%vecp%np=ip
  !
  if (ibacktr==ipback.AND.powcal%powres%flinm) then
     ! open file to record RZxi
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' trackm'')') self%ie,self%je
     write(icfile,'(''trackm'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     ! write track in RZxi coordinates
     ! de-quantise
     call position_invqtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     call position_writelis(powcal%odes%vecp,'track',nplot)
     ! re-quantise
     call position_qtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     ! close file
     close(nplot)
  end if

  if (ibacktr==ipback.AND.powcal%powres%flinx) then
     ! open file to record Cartesian track
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' trackx'')') self%ie,self%je
     write(icfile,'(''trackx'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     ! write track in Cartesian coordinates
     allocate(wposl%pos(ip), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     do l=1,ip
        ! dequantise
        zpos1%posvec=powcal%odes%vecp%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*powcal%odes%posk(l))&
 &      /powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=ip
     ! write track
     call position_writelis(wposl,'track',nplot)
     ! close file and deallocate
     deallocate(wposl%pos)
     close(nplot)
  end if

  if (ibacktr<ipback) then
     ibacktr=ipback
     !        zdfaca=-zdfaca
     powcal%odes%ndt=powcal%odes%ndt-1
     powcal%odes%n%stepmax=2*powcal%odes%n%stepmax
     powcal%odes%dt=-powcal%odes%dt
     goto 100
  end if

end subroutine powelt_move2
!---------------------------------------------------------------------
!> motion along fieldline global, ripple field OBSOLETE
subroutine powelt_move3(self,powcal,gshadl,btree)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure
  type(btree_t), intent(inout),optional :: btree !< btree data

  !! local
  character(*), parameter :: s_name='powelt_move3' !< subroutine name
  integer(ki4), parameter :: ipback=0 !< back-track test if unity
  real(kr4), dimension(3) :: zpos !< start position vector
  real(kr8), dimension(3) :: zposd !< start position vector
  real(kr8), dimension(3) :: xpath !< start path position vector Cartesians
  real(kr8) :: zzeta !< value of \f$ \zeta \f$
  real(kr8) :: zt0 !< start value of \f$ t \f$
  type(posang_t) :: zposang   !< posang data structure
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4), dimension(3) :: znormalm !< unit normal vector in mapped coords
  real(kr8), dimension(3) :: znormald !< unit normal vector
  real(kr8), dimension(3) :: zdfaca !< direction factor
  real(kr8) :: zs !< \f$ \pm 1 \f$
  real(kr8) :: zpdotn !< path dotted with normal vector
  real(kr8) :: zpdotnm !< path dotted with normal vector in mapped coords
  integer(ki4) :: ierr !< error flag
  integer(ki4) :: ip !< entries in track
  type(geobj_t) :: iobj !< geo object
  integer(ki4) :: inode  !< local variable
  integer(ki4) :: indt !< shorthand for current timestep number
  !> Non-zero at fieldline termination
  !! - (-2) Line hit terminating plane
  !! - (-1) Line left octree
  !! - Positive - number of triangle hit
  integer(ki4) :: nobjhit=0  !< .
  integer(ki4) :: ibacktr !< flag if back-tracking
  type(posnode_t) :: xo !< local variable
  type(posnode_t) :: xn !< local variable
  type(posnode_t) :: xm !< local variable
  real(kr8) :: zf1 !< fraction
  real(kr8) :: zxi !< value of \f$ \xi \f$
  real(kr8) :: domlen !< size of (periodic) domain in \f$ \xi \f$
  real(kr8) :: lenpath !< length of path in \f$ t \f$
  real(kr8) :: ztime !< time \f$ t \f$
  real(kr8) :: zdt !< time increment \f$ \Delta t \f$
  type(posvecl_t) :: zpos1   !< one position data
  type(posvecl_t) :: zpos2   !< one position data
  type(posveclis_t) :: rposl   !< list of position data
  real(kr4) :: phylenpath !< physical length of path
  logical :: lcoll   !< collision on path
  real(kr8) :: zpsim !< value of \f$ \psi \f$ at field line end (diagnostic)
  real(kr8), dimension(2) :: zk !< sector number

  inpow=powelt_addr(self,powcal%powres%npowe)
  ! check whether looking at small number of tracks
  if (powcal%n%ntrack>0) then
     ! is this element track in the test set
     do j=1,powcal%n%ntrack
        if (powcal%n%trackno(j)==inpow) go to 1
     end do
     return
  end if
1     continue
  domlen=powcal%powres%beq%domlen(3)
  ! mark powelt as unknown (2) by default
  powcal%powres%pow(inpow)=2
  ! needed in shadowed case
  if (powcal%n%shadow>0.AND..NOT.allocated(rposl%pos)) then
     allocate(rposl%pos(1),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else if (powcal%n%shadow<-1) then
     return
  end if
  !DVEC      allocate(work1(powcal%odes%n%stepmax),stat=status) !DVEC
  !DVEC      call log_alloc_check(m_name,s_name,0,status) !DVEC

  ! calculate, save and test initial position
  call powelt_pos(self,powcal,zpos)
  ! Examine start trajectory behaviour in RZxi
  call powelt_normalm(self,powcal,znormalm)
  znormald=znormalm
  zposd=zpos
  zt0=0
  ! start trajectory
  powcal%odes%ndt=1
  powcal%odes%t=0
  powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=zposd
  powcal%odes%posk(powcal%odes%ndt)=0
  !! evaluate I aka f at psi
  !      call spleval(powcal%powres%beq%f,powcal%powres%beq%mr,powcal%powres%beq%psiaxis,&
  !      powcal%powres%beq%psiqbdry,zpsi,zf,1)
  ! set initial time step and guess first new point of track
  zdfaca=powcal%powres%qfaca
  call odes_3rdstep1(powcal%odes,ierr,zdfaca,powcal%powres%beq%vacfld,&
 &powcal%powres%beq%psi,powcal%powres%beq%dpsidr,powcal%powres%beq%dpsidz,0)
  if (ierr>0) return
  xpath=powcal%odes%vecp%pos(powcal%odes%ndt+1)%posvec-zposd
  zpdotnm=dot_product(xpath,znormald)

  ! Check behaviour in Cartesians
  call powelt_normal(self,powcal,znormal)
  znormald=znormal
  ! ensure path is outward from surface
  xpath=0
  zs=1
  do j=0,1
     ! dequantise
     zpos1%posvec=powcal%odes%vecp%pos(powcal%odes%ndt+j)%posvec
     zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
     ! convert xi to zeta
     zzeta=zpos2%posvec(3)/powcal%powres%beq%nzets
     zpos2%posvec(3)=zzeta
     !
     zposang%pos=zpos2%posvec
     zposang%opt=1 ; zposang%units=0
     call posang_tfm(zposang,-3)
     xpath=zposang%pos+zs*xpath
     zs=-zs
  end do
  zpdotn=dot_product(xpath,znormald)
  if (zpdotn*zpdotnm<0) then
     call log_error(m_name,s_name,3,error_warning,'Starting direction uncertain')
     write(*,*) inpow, zpdotn, zpdotnm
     ! default to shadowed
     powcal%powres%pow(inpow)=0
  end if
  ! adjust sign of timestep
  powcal%odes%dt=sign(1._kr8,zpdotn)*powcal%odes%dt
  ! scaling factor no longer carries timestep sign
  zdfaca=powcal%powres%qfaca
  if (powcal%n%shadow>0) then
     ! find start point in tree
     rposl%pos(1)=powcal%odes%vecp%pos(powcal%odes%ndt)
     rposl%np=1
     iobj%geobj=1
     iobj%objtyp=1 ! for point
     call btree_mfind(btree,iobj,rposl,inode)
     if (inode<0) then
        call log_error(m_name,s_name,4,error_warning,'cannot locate start point in tree')
        return
     end if
     xo%posvec=rposl%pos(1)%posvec
     xo%node=inode
  end if
  ibacktr=0
  if (loutdt) write(170,*)
  if (loutdt) write(170,*)
  if (loutdt) write(170,*) '#',self%ie,self%je

100   continue
  ! time step loop for path
  lenpath=0
  ierr=0
  !     powcal%odes%near=0
  loop_path: do
     lcoll=.FALSE.
     nobjhit=0
     !     write(*,*) 'step=', powcal%odes%ndt

     ztime=powcal%odes%t

     call odes_3rdstep(powcal%odes,ierr,zdfaca,powcal%powres%beq%vacfld,&
 &   powcal%powres%beq%psi,powcal%powres%beq%dpsidr,powcal%powres%beq%dpsidz)
     ! check for termination due to solution failure on this field-line
     if (ierr>0) exit

     zdt=powcal%odes%t-ztime
     !DVEC      work1(powcal%odes%ndt-1)=zdt  !DVEC
     !DVEC      write(*,*) 'tstep'  !DVEC
     lenpath=lenpath+abs(zdt)
     if (loutdt) write(170,*) zdt

     if (powcal%n%shadow>0) then
        ! check for collision or exit from computational domain
        xn%posvec=powcal%odes%vecp%pos(powcal%odes%ndt)%posvec
        zxi=xn%posvec(3)
        ! check for leaving periodic cell
        if (zxi<powcal%powres%beq%n%ximin) then
           ! interpolate for xm=x(0) and check versus objects
           zf1=(powcal%powres%beq%n%ximin-xo%posvec(3))/(zxi-xo%posvec(3))
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           xm%node=0
           ! only terminate if trajectory is long enough
           if (abs(lenpath)>powcal%odes%n%termcon(1)) then
              call pcle_move(xo,xm,xo%node,gshadl,btree,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) then
              powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              exit
           end if
           powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)-1
           xo=xm
           xo%posvec(3)=powcal%powres%beq%n%ximax
           xo%node=0
           xn%posvec(3)=min(powcal%powres%beq%ximaxm, xn%posvec(3)+domlen)
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec

        else if (zxi>powcal%powres%beq%n%ximax) then
           ! interpolate for xm=x(2pi) and check versus objects
           zf1=(powcal%powres%beq%n%ximax-xo%posvec(3))/(zxi-xo%posvec(3))
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           xm%node=0
           ! only terminate if trajectory is long enough
           if (abs(lenpath)>powcal%odes%n%termcon(1)) then
              call pcle_move(xo,xm,xo%node,gshadl,btree,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) then
              powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              exit
           end if
           powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)+1
           xo=xm
           xo%posvec(3)=powcal%powres%beq%n%ximin
           xo%node=0
           xn%posvec(3)=max(powcal%powres%beq%ximinp, xn%posvec(3)-domlen)
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
        end if

        ! find new node and position if hits boundary
        xn%node=0
        ! only terminate if trajectory is long enough
        if (abs(lenpath)>powcal%odes%n%termcon(1)) then
           call pcle_move(xo,xn,xo%node,gshadl,btree,nobjhit)
           lcoll=(nobjhit/=0)
        end if
        if (lcoll) then
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
           call powelt_setpow(self,powcal,nobjhit,xn,inpow,gshadl)
           !D           write(*,*) 'Final posn',xn%posvec! writediagn
           exit
        end if
        ! check whether hit midplane
        if (lpmidplane) then
           if ((xn%posvec(2)-powcal%odes%n%termcon(3))*powcal%odes%n%ntermcon>0) then
              lcoll=.TRUE.
              ! interpolate for xm=x(midplane)
              zf1=(powcal%odes%n%termcon(3)-xo%posvec(2))/(xn%posvec(2)-xo%posvec(2))
              xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
              xm%node=0
              nobjhit=-2
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              zpsim=powcal%powres%geobjl%obj(inpow)%weight
              xn=xm
              exit
           end if
        end if
        ! save new data
        xo=xn
        ! check for near-collision and repeat step (once) if necessary TO DO???
     end if

     !DP !! psi diagnostic !DP
     !DP      zr=xo%posvec(1) !DP
     !DP      zz=xo%posvec(2) !DP
     !DP      call spl2d_eval(powcal%powres%beq%psi,zr,zz,zpsi) !DP
     !DP! dequantise xo !DP
     !DP      zpos1%posvec=xo%posvec !DP
     !DP      zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm) !DP
     !DP     write(*,*) lenpath,zpos2%posvec(1),zpos2%posvec(2),zpsi !DP
     ! check for termination due to overlong trajectory (number of steps)
     indt=powcal%odes%ndt
     if (indt>=powcal%odes%n%stepmax) exit
     ! or length of path in t
     ! (often covered when test exit from domain)
     if (abs(lenpath)>abs(powcal%odes%n%tmax-zt0)) exit
  end do loop_path


  ip=powcal%odes%ndt
  powcal%odes%vecp%np=ip
  !
  if (ibacktr==ipback.AND.powcal%powres%flinm) then
     ! open file to record RZxi
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' lenpath= '',1pg12.5,'' objhit= '',I8)') &
 &   self%ie,self%je,lenpath,nobjhit
     write(icfile,'(''trackm'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     ! write track in RZxi coordinates
     ! de-quantise
     call position_invqtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     do l=1,ip
        powcal%odes%vecp%pos(l)%posvec(3)=&
 &      powcal%odes%vecp%pos(l)%posvec(3)+2*const_pid*powcal%odes%posk(l)
     end do
     call position_writelis(powcal%odes%vecp,'track',nplot)
     do l=1,ip
        powcal%odes%vecp%pos(l)%posvec(3)=&
 &      powcal%odes%vecp%pos(l)%posvec(3)-2*const_pid*powcal%odes%posk(l)
     end do
     ! re-quantise
     call position_qtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     ! close file
     close(nplot)
  end if

  if (ibacktr==ipback.AND.powcal%powres%flinx) then
     ! write track in Cartesian coordinates
     ! first convert track
     allocate(wposl%pos(ip), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     do l=1,ip
        ! dequantise
        zpos1%posvec=powcal%odes%vecp%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*powcal%odes%posk(l))&
 &      /powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=ip
     call position_lenlis(wposl,phylenpath)
     !DVEC         do l=1,ip-1 !DVEC
     !DVEC         zr=sqrt(wposl%pos(l)%posvec(1)**2+wposl%pos(l)%posvec(2)**2)/1000 !DVEC
     !DVEC         wposl%pos(l)%posvec=(wposl%pos(l+1)%posvec-wposl%pos(l)%posvec)/& !DVEC
     !DVEC      &  (zr*work1(l)) !DVEC
     !DVEC         end do !DVEC
     !DVEC         wposl%pos(ip)=wposl%pos(ip-1) !DVEC
     !DVEC! write track (or field vectors if !DVEC lines executed)

     ! open file to record Cartesian track and write
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' lenpath= '',1pg12.5,'' objhit= '',I8)') &
 &   self%ie,self%je,phylenpath,nobjhit
     !        write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' trackx'')') self%ie,self%je
     write(icfile,'(''trackx'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     call position_writelis(wposl,'track',nplot)
     !DIAG!   dump end point !DIAG
     !DIAG!        call position_writev(wposl%pos(ip),88) !DIAG
     ! close file and deallocate
     deallocate(wposl%pos)
     !DVEC         deallocate(work1)  !DVEC
     close(nplot)
  end if

  if (nobjhit<=GEOBJ_OFFSET.AND.ibacktr==ipback.AND.powcal%powres%flinends) then
     ! write first and last points in Cartesian coordinates
     ! first convert track points
     allocate(wposl%pos(3), stat=status)
     call log_alloc_check(m_name,s_name,11,status)
     wposl%pos(1)%posvec=powcal%odes%vecp%pos(1)%posvec
     zk(1)=powcal%odes%posk(1)
     wposl%pos(2)%posvec=xm%posvec
     zk(2)=powcal%odes%posk(ip)
     do l=1,2
        ! dequantise
        zpos1%posvec=wposl%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*zk(l))/powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=3
     ! open file to record Cartesian track and write
     wposl%pos(3)%posvec(1)=self%ie
     wposl%pos(3)%posvec(2)=self%je
     wposl%pos(3)%posvec(3)=zpsim
     do l=1,3
        call position_writev(wposl%pos(l),powcal%powres%nflends)
     end do
     ! deallocate
     deallocate(wposl%pos)
  end if

  if (ibacktr<ipback) then
     ibacktr=ipback
     !        zdfaca=-zdfaca
     powcal%odes%ndt=powcal%odes%ndt-1
     powcal%odes%n%stepmax=2*powcal%odes%n%stepmax
     powcal%odes%dt=-powcal%odes%dt
     goto 100
  end if

end subroutine powelt_move3
!---------------------------------------------------------------------
!> motion along fieldline global, ripple field
subroutine powelt_move4(self,powcal,gshadl,btree)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure
  type(btree_t), intent(inout),optional :: btree !< btree data

  !! local
  character(*), parameter :: s_name='powelt_move4' !< subroutine name
  integer(ki4), parameter :: ipback=0 !< back-track test if unity
  real(kr4), dimension(3) :: zpos !< start position vector
  real(kr8), dimension(3) :: zposd !< start position vector
  real(kr8), dimension(3) :: xpath !< start path position vector Cartesians
  real(kr8) :: zzeta !< value of \f$ \zeta \f$
  real(kr8) :: zt0 !< start value of \f$ t \f$
  type(posang_t) :: zposang   !< posang data structure
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4), dimension(3) :: znormalm !< unit normal vector in mapped coords
  real(kr8), dimension(3) :: znormald !< unit normal vector
  real(kr8), dimension(3) :: zdfaca !< direction factor
  real(kr8) :: zs !< \f$ \pm 1 \f$
  real(kr8) :: zpdotn !< path dotted with normal vector
  real(kr8) :: zpdotnm !< path dotted with normal vector in mapped coords
  integer(ki4) :: ierr !< error flag
  integer(ki4) :: ip !< entries in track
  type(geobj_t) :: iobj !< geo object
  integer(ki4) :: inode  !< local variable
  integer(ki4) :: indt !< shorthand for current timestep number
  !> Non-zero at fieldline termination
  !! - (-2) Line hit terminating plane
  !! - (-1) Line left octree
  !! - Positive - number of triangle hit
  integer(ki4) :: nobjhit=0  !< .
  integer(ki4) :: ibacktr !< flag if back-tracking
  type(posnode_t) :: xo !< local variable
  type(posnode_t) :: xn !< local variable
  type(posnode_t) :: xm !< local variable
  real(kr8) :: zf1 !< fraction of path length
  real(kr8) :: zxi !< value of \f$ \xi \f$
  real(kr8) :: domlen !< size of (periodic) domain in \f$ \xi \f$
  real(kr8) :: lenpath !< length of path in \f$ t \f$
  real(kr8) :: ztime !< time \f$ t \f$
  real(kr8) :: zdt !< time increment \f$ \Delta t \f$
  type(posvecl_t) :: zpos1   !< one position data
  type(posvecl_t) :: zpos2   !< one position data
  type(posveclis_t) :: rposl   !< list of position data
  real(kr4) :: phylenpath !< physical length of path
  logical :: lcoll   !< collision on path
  real(kr8) :: zpsim !< value of \f$ \psi \f$ at field line end (diagnostic)
  real(kr8), dimension(2) :: zk !< sector number

  inpow=powelt_addr(self,powcal%powres%npowe)
  ! check whether looking at small number of tracks
  if (powcal%n%ntrack>0) then
     ! is this element track in the test set
     do j=1,powcal%n%ntrack
        if (powcal%n%trackno(j)==inpow) go to 1
     end do
     return
  end if
1     continue
  domlen=powcal%powres%beq%domlen(3)
  ! mark powelt as unknown (2) by default
  powcal%powres%pow(inpow)=2
  ! needed in shadowed case
  if (powcal%n%shadow>0.AND..NOT.allocated(rposl%pos)) then
     allocate(rposl%pos(1),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else if (powcal%n%shadow<-1) then
     return
  end if
  !DVEC      allocate(work1(powcal%odes%n%stepmax),stat=status) !DVEC
  !DVEC      call log_alloc_check(m_name,s_name,0,status) !DVEC

  ! calculate, save and test initial position
  call powelt_pos(self,powcal,zpos)
  ! Examine start trajectory behaviour in RZxi
  call powelt_normalm(self,powcal,znormalm)
  znormald=znormalm
  zposd=zpos
  zt0=0
  ! start trajectory
  powcal%odes%ndt=1
  powcal%odes%t=0
  powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=zposd
  powcal%odes%posk(powcal%odes%ndt)=0
  !! evaluate I aka f at psi
  !      call spleval(powcal%powres%beq%f,powcal%powres%beq%mr,powcal%powres%beq%psiaxis,&
  !      powcal%powres%beq%psiqbdry,zpsi,zf,1)
  ! set initial time step and guess first new point of track
  zdfaca=powcal%powres%qfaca
  call odes_3rdstep1(powcal%odes,ierr,zdfaca,powcal%powres%beq%vacfld,&
 &powcal%powres%beq%psi,powcal%powres%beq%dpsidr,powcal%powres%beq%dpsidz,0)
  if (ierr>0) return
  xpath=powcal%odes%vecp%pos(powcal%odes%ndt+1)%posvec-zposd
  zpdotnm=dot_product(xpath,znormald)

  ! Check behaviour in Cartesians
  call powelt_normal(self,powcal,znormal)
  znormald=znormal
  ! ensure path is outward from surface
  xpath=0
  zs=1
  do j=0,1
     ! dequantise
     zpos1%posvec=powcal%odes%vecp%pos(powcal%odes%ndt+j)%posvec
     zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
     ! convert xi to zeta
     zzeta=zpos2%posvec(3)/powcal%powres%beq%nzets
     zpos2%posvec(3)=zzeta
     !
     zposang%pos=zpos2%posvec
     zposang%opt=1 ; zposang%units=0
     call posang_tfm(zposang,-3)
     xpath=zposang%pos+zs*xpath
     zs=-zs
  end do
  zpdotn=dot_product(xpath,znormald)
  if (zpdotn*zpdotnm<0) then
     call log_error(m_name,s_name,3,error_warning,'Starting direction uncertain')
     write(*,*) inpow, zpdotn, zpdotnm
     ! default to shadowed
     powcal%powres%pow(inpow)=0
  end if
  ! adjust sign of timestep
  powcal%odes%dt=sign(1._kr8,zpdotn)*powcal%odes%dt
  ! scaling factor no longer carries timestep sign
  zdfaca=powcal%powres%qfaca
  if (powcal%n%shadow>0) then
     ! find start point in tree
     rposl%pos(1)=powcal%odes%vecp%pos(powcal%odes%ndt)
     rposl%np=1
     iobj%geobj=1
     iobj%objtyp=1 ! for point
     call btree_mfind(btree,iobj,rposl,inode)
     if (inode<0) then
        call log_error(m_name,s_name,4,error_warning,'cannot locate start point in tree')
        return
     end if
     xo%posvec=rposl%pos(1)%posvec
     xo%node=inode
  end if
  ibacktr=0
  if (loutdt) write(170,*)
  if (loutdt) write(170,*)
  if (loutdt) write(170,*) '#',self%ie,self%je

100   continue
  ! initialise loop for path
  lenpath=0
  ierr=0
  ! ignore midplane if termination planes present
  if (powcal%n%ltermplane) then
     ilpmidplane=.FALSE.
     powcal%n%termp%termstore(1,:)=(/1.e30,1.e30,1.e30/)
  else
     ilpmidplane=lpmidplane
  end if
  !     powcal%odes%near=0
  ! time step loop for path
  loop_path: do
     lcoll=.FALSE.
     nobjhit=0
     zpsim=0
     !     write(*,*) 'step=', powcal%odes%ndt

     ztime=powcal%odes%t

     call odes_3rdstep(powcal%odes,ierr,zdfaca,powcal%powres%beq%vacfld,&
 &   powcal%powres%beq%psi,powcal%powres%beq%dpsidr,powcal%powres%beq%dpsidz)
     ! check for termination due to solution failure on this field-line
     if (ierr>0) exit

     zdt=powcal%odes%t-ztime
     !DVEC      work1(powcal%odes%ndt-1)=zdt  !DVEC
     !DVEC      write(*,*) 'tstep'  !DVEC
     lenpath=lenpath+abs(zdt)
     if (loutdt) write(170,*) zdt

     if (powcal%n%shadow>0) then
        ! check for collision or exit from computational domain
        xn%posvec=powcal%odes%vecp%pos(powcal%odes%ndt)%posvec
        zxi=xn%posvec(3)
        ! check for leaving periodic cell
        if (zxi<powcal%powres%beq%n%ximin) then
           ! interpolate for xm=x(0) and check versus objects
           zf1=(powcal%powres%beq%n%ximin-xo%posvec(3))/(zxi-xo%posvec(3))
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           xm%node=0
           ! only terminate if trajectory is long enough
           if (abs(lenpath)>powcal%odes%n%termcon(1)) then
              call pcle_movet(xo,xm,xo%node,gshadl,btree, &
 &            powcal%n%termp,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) then
              powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              exit
           end if
           powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)-1
           xo=xm
           xo%posvec(3)=powcal%powres%beq%n%ximax
           xo%node=0
           xn%posvec(3)=min(powcal%powres%beq%ximaxm, xn%posvec(3)+domlen)
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec

        else if (zxi>powcal%powres%beq%n%ximax) then
           ! interpolate for xm=x(2pi) and check versus objects
           zf1=(powcal%powres%beq%n%ximax-xo%posvec(3))/(zxi-xo%posvec(3))
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           xm%node=0
           ! only terminate if trajectory is long enough
           if (abs(lenpath)>powcal%odes%n%termcon(1)) then
              call pcle_movet(xo,xm,xo%node,gshadl,btree, &
 &            powcal%n%termp,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) then
              powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              ! record psi at boundary if hit termplane boundary
              if (nobjhit<=GEOBJ_OFFSET) zpsim=powcal%powres%geobjl%obj(inpow)%weight
              exit
           end if
           powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)+1
           xo=xm
           xo%posvec(3)=powcal%powres%beq%n%ximin
           xo%node=0
           xn%posvec(3)=max(powcal%powres%beq%ximinp, xn%posvec(3)-domlen)
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
        end if

        ! find new node and position if hits boundary
        xn%node=0
        ! only terminate if trajectory is long enough
        if (abs(lenpath)>powcal%odes%n%termcon(1)) then
           call pcle_movet(xo,xn,xo%node,gshadl,btree, &
 &         powcal%n%termp,nobjhit)
           lcoll=(nobjhit/=0)
        end if
        if (lcoll) then
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
           call powelt_setpow(self,powcal,nobjhit,xn,inpow,gshadl)
           ! record psi at boundary if hit termplane boundary
           if (nobjhit<=GEOBJ_OFFSET) zpsim=powcal%powres%geobjl%obj(inpow)%weight
           !D           write(*,*) 'Final posn',xn%posvec! writediagn
           exit
        end if
        ! check whether hit midplane
        if (ilpmidplane) then
           if ((xn%posvec(2)-powcal%odes%n%termcon(3))*powcal%odes%n%ntermcon>0) then
              lcoll=.TRUE.
              ! interpolate for xm=x(midplane)
              zf1=(powcal%odes%n%termcon(3)-xo%posvec(2))/(xn%posvec(2)-xo%posvec(2))
              ! fix up in case zf1 garbage caused by xi boundary interaction
              zf1=max(min(zf1,1.),0.)
              xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
              xm%node=0
              nobjhit=-2
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              zpsim=powcal%powres%geobjl%obj(inpow)%weight
              xn=xm
              exit
           end if
        end if
        ! save new data
        xo=xn
        ! check for near-collision and repeat step (once) if necessary TO DO???
     end if

     !DP !! psi diagnostic !DP
     !DP      zr=xo%posvec(1) !DP
     !DP      zz=xo%posvec(2) !DP
     !DP      call spl2d_eval(powcal%powres%beq%psi,zr,zz,zpsi) !DP
     !DP! dequantise xo !DP
     !DP      zpos1%posvec=xo%posvec !DP
     !DP      zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm) !DP
     !DP     write(*,*) lenpath,zpos2%posvec(1),zpos2%posvec(2),zpsi !DP
     ! check for termination due to overlong trajectory (number of steps)
     indt=powcal%odes%ndt
     if (indt>=powcal%odes%n%stepmax) exit
     ! or length of path in t
     ! (often covered when test exit from domain)
     if (abs(lenpath)>abs(powcal%odes%n%tmax-zt0)) exit
  end do loop_path


  ip=powcal%odes%ndt
  powcal%odes%vecp%np=ip
  !
  if (ibacktr==ipback.AND.powcal%powres%flinm) then
     ! open file to record RZxi
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' lenpath= '',1pg12.5,'' objhit= '',I8)') &
 &   self%ie,self%je,lenpath,nobjhit
     write(icfile,'(''trackm'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     ! write track in RZxi coordinates
     ! de-quantise
     call position_invqtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     do l=1,ip
        powcal%odes%vecp%pos(l)%posvec(3)=&
 &      powcal%odes%vecp%pos(l)%posvec(3)+2*const_pid*powcal%odes%posk(l)
     end do
     call position_writelis(powcal%odes%vecp,'track',nplot)
     do l=1,ip
        powcal%odes%vecp%pos(l)%posvec(3)=&
 &      powcal%odes%vecp%pos(l)%posvec(3)-2*const_pid*powcal%odes%posk(l)
     end do
     ! re-quantise
     call position_qtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     ! close file
     close(nplot)
  end if

  if (ibacktr==ipback.AND.powcal%powres%flinx) then
     ! write track in Cartesian coordinates
     ! first convert track
     allocate(wposl%pos(ip), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     do l=1,ip
        ! dequantise
        zpos1%posvec=powcal%odes%vecp%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*powcal%odes%posk(l))&
 &      /powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=ip
     call position_lenlis(wposl,phylenpath)
     !DVEC         do l=1,ip-1 !DVEC
     !DVEC         zr=sqrt(wposl%pos(l)%posvec(1)**2+wposl%pos(l)%posvec(2)**2)/1000 !DVEC
     !DVEC         wposl%pos(l)%posvec=(wposl%pos(l+1)%posvec-wposl%pos(l)%posvec)/& !DVEC
     !DVEC      &  (zr*work1(l)) !DVEC
     !DVEC         end do !DVEC
     !DVEC         wposl%pos(ip)=wposl%pos(ip-1) !DVEC
     !DVEC! write track (or field vectors if !DVEC lines executed)

     ! open file to record Cartesian track and write
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' lenpath= '',1pg12.5,'' objhit= '',I8)') &
 &   self%ie,self%je,phylenpath,nobjhit
     !        write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' trackx'')') self%ie,self%je
     write(icfile,'(''trackx'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     call position_writelis(wposl,'track',nplot)
     !DIAG!   dump end point !DIAG
     !DIAG!        call position_writev(wposl%pos(ip),88) !DIAG
     ! close file and deallocate
     deallocate(wposl%pos)
     !DVEC         deallocate(work1)  !DVEC
     close(nplot)
  end if

  if (nobjhit<=GEOBJ_OFFSET.AND.ibacktr==ipback.AND.powcal%powres%flinends) then
     ! write first and last points in Cartesian coordinates
     ! first convert track points
     allocate(wposl%pos(3), stat=status)
     call log_alloc_check(m_name,s_name,11,status)
     wposl%pos(1)%posvec=powcal%odes%vecp%pos(1)%posvec
     zk(1)=powcal%odes%posk(1)
     wposl%pos(2)%posvec=xm%posvec
     zk(2)=powcal%odes%posk(ip)
     do l=1,2
        ! dequantise
        zpos1%posvec=wposl%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*zk(l))/powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=3
     ! open file to record Cartesian track and write
     wposl%pos(3)%posvec(1)=self%ie
     wposl%pos(3)%posvec(2)=self%je
     wposl%pos(3)%posvec(3)=zpsim
     do l=1,3
        call position_writev(wposl%pos(l),powcal%powres%nflends)
     end do
     ! deallocate
     deallocate(wposl%pos)
  end if

  if (ibacktr<ipback) then
     ibacktr=ipback
     !        zdfaca=-zdfaca
     powcal%odes%ndt=powcal%odes%ndt-1
     powcal%odes%n%stepmax=2*powcal%odes%n%stepmax
     powcal%odes%dt=-powcal%odes%dt
     goto 100
  end if

end subroutine powelt_move4
!---------------------------------------------------------------------
!> motion along fieldline 5th scheme for testing
subroutine powelt_move5(self,powcal,gshadl,btree)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure
  type(btree_t), intent(inout),optional :: btree !< btree data

  !! local
  character(*), parameter :: s_name='powelt_move5' !< subroutine name
  integer(ki4), parameter :: ipback=0 !< back-track test if unity
  real(kr4), dimension(3) :: zpos !< start position vector
  real(kr8), dimension(3) :: zposd !< start position vector
  real(kr8), dimension(3) :: xpath !< start path position vector Cartesians
  real(kr8) :: zpsi !< value of \f$ \psi \f$
  real(kr8) :: zzeta !< value of \f$ \zeta \f$
  real(kr8) :: zf !< value of \f$ f \f$
  real(kr8) :: zt0 !< start value of \f$ t \f$
  type(posang_t) :: zposang   !< posang data structure
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4), dimension(3) :: znormalm !< unit normal vector in mapped coords
  real(kr8), dimension(3) :: znormald !< unit normal vector
  real(kr8), dimension(6) :: zdfaca !< direction factor and offsets
  real(kr8) :: zs !< \f$ \pm 1 \f$
  real(kr8) :: zpdotn !< path dotted with normal vector
  real(kr8) :: zpdotnm !< path dotted with normal vector in mapped coords
  integer(ki4) :: ierr !< error flag
  integer(ki4) :: ip !< entries in track
  type(geobj_t) :: iobj !< geo object
  integer(ki4) :: inode  !< local variable
  integer(ki4) :: indt !< shorthand for current timestep number
  !> Non-zero at fieldline termination
  !! - (-2) Line hit terminating plane
  !! - (-1) Line left octree
  !! - Positive - number of triangle hit
  integer(ki4) :: nobjhit=0  !< .
  integer(ki4) :: ibacktr !< flag if back-tracking
  type(posnode_t) :: xo !< local variable
  type(posnode_t) :: xn !< local variable
  type(posnode_t) :: xm !< local variable
  real(kr8) :: zf1 !< fraction of path length
  real(kr8) :: zxi !< value of \f$ \xi \f$
  real(kr8) :: domlen !< size of (periodic) domain in \f$ \xi \f$
  real(kr8) :: lenpath !< length of path in \f$ t \f$
  real(kr8) :: ztime !< time \f$ t \f$
  real(kr8) :: zdt !< time increment \f$ \Delta t \f$
  type(posvecl_t) :: zpos1   !< one position data
  type(posvecl_t) :: zpos2   !< one position data
  type(posveclis_t) :: rposl   !< list of position data
  real(kr4) :: phylenpath !< physical length of path
  logical :: lcoll   !< collision on path
  real(kr8) :: zpsim !< value of \f$ \psi \f$ at field line end (diagnostic)
  real(kr8), dimension(2) :: zk !< sector number
  integer(ki4), save :: firstcall  !< flag up first step of new trajectory

  inpow=powelt_addr(self,powcal%powres%npowe)
  ! check whether looking at small number of tracks
  if (powcal%n%ntrack>0) then
     ! is this element track in the test set
     do j=1,powcal%n%ntrack
        if (powcal%n%trackno(j)==inpow) go to 1
     end do
     return
  end if
1     continue
  domlen=powcal%powres%beq%domlen(3)
  ! mark powelt as unknown (2) by default
  powcal%powres%pow(inpow)=2
  ! needed in shadowed case
  if (powcal%n%shadow>0.AND..NOT.allocated(rposl%pos)) then
     allocate(rposl%pos(1),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else if (powcal%n%shadow<-1) then
     return
  end if
  !DVEC      allocate(work1(powcal%odes%n%stepmax),stat=status) !DVEC
  !DVEC      call log_alloc_check(m_name,s_name,0,status) !DVEC

  ! calculate, save and test initial position
  call powelt_pos(self,powcal,zpos)
  ! Examine start trajectory behaviour in RZxi
  call powelt_normalm(self,powcal,znormalm)
  znormald=znormalm
  zposd=zpos
  zt0=0
  ! start trajectory
  powcal%odes%ndt=1
  powcal%odes%t=0
  powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=zposd
  powcal%odes%posk(powcal%odes%ndt)=0
  call spl2d_eval(powcal%powres%beq%psi,zposd(1),zposd(2),zpsi)
  ! evaluate I aka f at psi
  call spleval(powcal%powres%beq%f,powcal%powres%beq%mr,powcal%powres%beq%psiaxis,&
  powcal%powres%beq%psiqbdry,zpsi,zf,1)
  ! set initial time step and guess first new point of track
  zdfaca(1:3)=powcal%powres%qfaca
  zdfaca(4:6)=powcal%powres%geobjl%quantfm%offvec
  call odes_4thstep1(powcal%odes,ierr,zdfaca,zf,&
 &powcal%powres%beq%psi,powcal%powres%beq%dpsidr,powcal%powres%beq%dpsidz,0)
  if (ierr>0) return
  xpath=powcal%odes%vecp%pos(powcal%odes%ndt+1)%posvec-zposd
  zpdotnm=dot_product(xpath,znormald)

  ! Check behaviour in Cartesians
  call powelt_normal(self,powcal,znormal)
  znormald=znormal
  ! ensure path is outward from surface
  xpath=0
  zs=1
  do j=0,1
     ! dequantise
     zpos1%posvec=powcal%odes%vecp%pos(powcal%odes%ndt+j)%posvec
     zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
     ! convert xi to zeta
     zzeta=zpos2%posvec(3)/powcal%powres%beq%nzets
     zpos2%posvec(3)=zzeta
     !
     zposang%pos=zpos2%posvec
     zposang%opt=1 ; zposang%units=0
     call posang_tfm(zposang,-3)
     xpath=zposang%pos+zs*xpath
     zs=-zs
  end do
  zpdotn=dot_product(xpath,znormald)
  if (zpdotn*zpdotnm<0) then
     call log_error(m_name,s_name,3,error_warning,'Starting direction uncertain')
     write(*,*) inpow, zpdotn, zpdotnm
     ! default to shadowed
     powcal%powres%pow(inpow)=0
  end if
  ! adjust sign of timestep
  powcal%odes%dt=sign(1._kr8,zpdotn)*powcal%odes%dt
  ! scaling factor no longer carries timestep sign
  zdfaca(1:3)=powcal%powres%qfaca
  zdfaca(4:6)=powcal%powres%geobjl%quantfm%offvec
  if (powcal%n%shadow>0) then
     ! find start point in tree
     rposl%pos(1)=powcal%odes%vecp%pos(powcal%odes%ndt)
     rposl%np=1
     iobj%geobj=1
     iobj%objtyp=1 ! for point
     call btree_mfind(btree,iobj,rposl,inode)
     if (inode<0) then
        call log_error(m_name,s_name,4,error_warning,'cannot locate start point in tree')
        return
     end if
     xo%posvec=rposl%pos(1)%posvec
     xo%node=inode
  end if
  ibacktr=0
  if (loutdt) write(170,*)
  if (loutdt) write(170,*)
  if (loutdt) write(170,*) '#',self%ie,self%je

100   continue
  ! initialise loop for path
  firstcall=1
  lenpath=0
  ierr=0
  ! ignore midplane if termination planes present
  if (powcal%n%ltermplane) then
     ilpmidplane=.FALSE.
     powcal%n%termp%termstore(1,:)=(/1.e30,1.e30,1.e30/)
  else
     ilpmidplane=lpmidplane
  end if
  !     powcal%odes%near=0
  ! time step loop for path
  loop_path: do
     lcoll=.FALSE.
     nobjhit=0
     zpsim=0
     !     write(*,*) 'step=', powcal%odes%ndt

     ztime=powcal%odes%t

     zposd=xo%posvec
     call spl2d_eval(powcal%powres%beq%psi,zposd(1),zposd(2),zpsi)
     ! evaluate I aka f at psi
     call spleval(powcal%powres%beq%f,powcal%powres%beq%mr,powcal%powres%beq%psiaxis,&
     powcal%powres%beq%psiqbdry,zpsi,zf,1)
     call odes_4thstepcon(powcal%odes,ierr,zdfaca,firstcall,zf,&
 &   powcal%powres%beq%psi,powcal%powres%beq%dpsidr,powcal%powres%beq%dpsidz)
     ! check for termination due to solution failure on this field-line
     if (ierr>0) exit

     zdt=powcal%odes%t-ztime
     !DVEC      work1(powcal%odes%ndt-1)=zdt  !DVEC
     !DVEC      write(*,*) 'tstep'  !DVEC
     lenpath=lenpath+abs(zdt)
     if (loutdt) write(170,*) zdt

     if (powcal%n%shadow>0) then
        ! check for collision or exit from computational domain
        xn%posvec=powcal%odes%vecp%pos(powcal%odes%ndt)%posvec
        zxi=xn%posvec(3)
        ! check for leaving periodic cell
        if (zxi<powcal%powres%beq%n%ximin) then
           ! interpolate for xm=x(0) and check versus objects
           zf1=(powcal%powres%beq%n%ximin-xo%posvec(3))/(zxi-xo%posvec(3))
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           xm%node=0
           ! only terminate if trajectory is long enough
           if (abs(lenpath)>powcal%odes%n%termcon(1)) then
              call pcle_movet(xo,xm,xo%node,gshadl,btree, &
 &            powcal%n%termp,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) then
              powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              exit
           end if
           powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)-1
           xo=xm
           xo%posvec(3)=powcal%powres%beq%n%ximax
           xo%node=0
           xn%posvec(3)=min(powcal%powres%beq%ximaxm, xn%posvec(3)+domlen)
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec

        else if (zxi>powcal%powres%beq%n%ximax) then
           ! interpolate for xm=x(2pi) and check versus objects
           zf1=(powcal%powres%beq%n%ximax-xo%posvec(3))/(zxi-xo%posvec(3))
           xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
           xm%node=0
           ! only terminate if trajectory is long enough
           if (abs(lenpath)>powcal%odes%n%termcon(1)) then
              call pcle_movet(xo,xm,xo%node,gshadl,btree, &
 &            powcal%n%termp,nobjhit)
              lcoll=(nobjhit/=0)
           end if
           if (lcoll) then
              powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xm%posvec
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              ! record psi at boundary if hit termplane boundary
              if (nobjhit<=GEOBJ_OFFSET) zpsim=powcal%powres%geobjl%obj(inpow)%weight
              exit
           end if
           powcal%odes%posk(powcal%odes%ndt)=powcal%odes%posk(powcal%odes%ndt)+1
           xo=xm
           xo%posvec(3)=powcal%powres%beq%n%ximin
           xo%node=0
           xn%posvec(3)=max(powcal%powres%beq%ximinp, xn%posvec(3)-domlen)
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
        end if

        ! find new node and position if hits boundary
        xn%node=0
        ! only terminate if trajectory is long enough
        if (abs(lenpath)>powcal%odes%n%termcon(1)) then
           call pcle_movet(xo,xn,xo%node,gshadl,btree, &
 &         powcal%n%termp,nobjhit)
           lcoll=(nobjhit/=0)
        end if
        if (lcoll) then
           powcal%odes%vecp%pos(powcal%odes%ndt)%posvec=xn%posvec
           call powelt_setpow(self,powcal,nobjhit,xn,inpow,gshadl)
           ! record psi at boundary if hit termplane boundary
           if (nobjhit<=GEOBJ_OFFSET) zpsim=powcal%powres%geobjl%obj(inpow)%weight
           !D           write(*,*) 'Final posn',xn%posvec! writediagn
           exit
        end if
        ! check whether hit midplane
        if (ilpmidplane) then
           if ((xn%posvec(2)-powcal%odes%n%termcon(3))*powcal%odes%n%ntermcon>0) then
              lcoll=.TRUE.
              ! interpolate for xm=x(midplane)
              zf1=(powcal%odes%n%termcon(3)-xo%posvec(2))/(xn%posvec(2)-xo%posvec(2))
              ! fix up in case zf1 garbage caused by xi boundary interaction
              zf1=max(min(zf1,1.),0.)
              xm%posvec=(1-zf1)*xo%posvec+zf1*xn%posvec
              xm%node=0
              nobjhit=-2
              call powelt_setpow(self,powcal,nobjhit,xm,inpow,gshadl)
              zpsim=powcal%powres%geobjl%obj(inpow)%weight
              xn=xm
              exit
           end if
        end if
        ! save new data
        xo=xn
        ! check for near-collision and repeat step (once) if necessary TO DO???
     end if

     !DP !! psi diagnostic !DP
     !DP      zr=xo%posvec(1) !DP
     !DP      zz=xo%posvec(2) !DP
     !DP      call spl2d_eval(powcal%powres%beq%psi,zr,zz,zpsi) !DP
     !DP! dequantise xo !DP
     !DP      zpos1%posvec=xo%posvec !DP
     !DP      zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm) !DP
     !DP     write(*,*) lenpath,zpos2%posvec(1),zpos2%posvec(2),zpsi !DP
     ! check for termination due to overlong trajectory (number of steps)
     indt=powcal%odes%ndt
     if (indt>=powcal%odes%n%stepmax) exit
     ! or length of path in t
     ! (often covered when test exit from domain)
     if (abs(lenpath)>abs(powcal%odes%n%tmax-zt0)) exit
  end do loop_path


  ip=powcal%odes%ndt
  powcal%odes%vecp%np=ip
  !
  if (ibacktr==ipback.AND.powcal%powres%flinm) then
     ! open file to record RZxi
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' lenpath= '',1pg12.5,'' objhit= '',I8)') &
 &   self%ie,self%je,lenpath,nobjhit
     write(icfile,'(''trackm'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     ! write track in RZxi coordinates
     ! de-quantise
     call position_invqtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     do l=1,ip
        powcal%odes%vecp%pos(l)%posvec(3)=&
 &      powcal%odes%vecp%pos(l)%posvec(3)+2*const_pid*powcal%odes%posk(l)
     end do
     call position_writelis(powcal%odes%vecp,'track',nplot)
     do l=1,ip
        powcal%odes%vecp%pos(l)%posvec(3)=&
 &      powcal%odes%vecp%pos(l)%posvec(3)-2*const_pid*powcal%odes%posk(l)
     end do
     ! re-quantise
     call position_qtfmlis(powcal%odes%vecp,powcal%powres%geobjl%quantfm)
     ! close file
     close(nplot)
  end if

  if (ibacktr==ipback.AND.powcal%powres%flinx) then
     ! write track in Cartesian coordinates
     ! first convert track
     allocate(wposl%pos(ip), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     do l=1,ip
        ! dequantise
        zpos1%posvec=powcal%odes%vecp%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*powcal%odes%posk(l))&
 &      /powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=ip
     call position_lenlis(wposl,phylenpath)
     !DVEC         do l=1,ip-1 !DVEC
     !DVEC         zr=sqrt(wposl%pos(l)%posvec(1)**2+wposl%pos(l)%posvec(2)**2)/1000 !DVEC
     !DVEC         wposl%pos(l)%posvec=(wposl%pos(l+1)%posvec-wposl%pos(l)%posvec)/& !DVEC
     !DVEC      &  (zr*work1(l)) !DVEC
     !DVEC         end do !DVEC
     !DVEC         wposl%pos(ip)=wposl%pos(ip-1) !DVEC
     !DVEC! write track (or field vectors if !DVEC lines executed)

     ! open file to record Cartesian track and write
     write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' lenpath= '',1pg12.5,'' objhit= '',I8)') &
 &   self%ie,self%je,phylenpath,nobjhit
     !        write(ibuff,'(''elt= '',I6,'' sub= '',I2,'' trackx'')') self%ie,self%je
     write(icfile,'(''trackx'',I7.7,I2.2)') self%ie,self%je
     call vfile_init(icfile,ibuff,nplot)
     call position_writelis(wposl,'track',nplot)
     !DIAG!   dump end point !DIAG
     !DIAG!        call position_writev(wposl%pos(ip),88) !DIAG
     ! close file and deallocate
     deallocate(wposl%pos)
     !DVEC         deallocate(work1)  !DVEC
     close(nplot)
  end if

  if (nobjhit<=GEOBJ_OFFSET.AND.ibacktr==ipback.AND.powcal%powres%flinends) then
     ! write first and last points in Cartesian coordinates
     ! first convert track points
     allocate(wposl%pos(3), stat=status)
     call log_alloc_check(m_name,s_name,11,status)
     wposl%pos(1)%posvec=powcal%odes%vecp%pos(1)%posvec
     zk(1)=powcal%odes%posk(1)
     wposl%pos(2)%posvec=xm%posvec
     zk(2)=powcal%odes%posk(ip)
     do l=1,2
        ! dequantise
        zpos1%posvec=wposl%pos(l)%posvec
        zpos2=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
        ! convert xi to zeta
        zzeta=(zpos2%posvec(3)+2*const_pid*zk(l))/powcal%powres%beq%nzets
        zpos2%posvec(3)=zzeta
        zposang%pos=zpos2%posvec
        zposang%opt=1 ; zposang%units=0
        call posang_tfm(zposang,-3)
        wposl%pos(l)%posvec=zposang%pos
     end do
     wposl%np=3
     ! open file to record Cartesian track and write
     wposl%pos(3)%posvec(1)=self%ie
     wposl%pos(3)%posvec(2)=self%je
     wposl%pos(3)%posvec(3)=zpsim
     do l=1,3
        call position_writev(wposl%pos(l),powcal%powres%nflends)
     end do
     ! deallocate
     deallocate(wposl%pos)
  end if

  if (ibacktr<ipback) then
     ibacktr=ipback
     !        zdfaca(1:3)=-zdfaca(1:3)
     powcal%odes%ndt=powcal%odes%ndt-1
     powcal%odes%n%stepmax=2*powcal%odes%n%stepmax
     powcal%odes%dt=-powcal%odes%dt
     goto 100
  end if

end subroutine powelt_move5
!---------------------------------------------------------------------
!> average value of power on element
subroutine powelt_avg(self,powcal)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure

  !! local
  character(*), parameter :: s_name='powelt_avg' !< subroutine name
  type(powelt_t) :: zself   !< object data structure
  real(kr4) :: zsum !< sum of powers

  zself%ie=self%ie
  zsum=0
  ! find extremal values of pow on element
  ! find number in level from je, although je usually set to iinlevel
  ilevel=powelt_level_table(self%je)
  iinlevel=powelt_table(ilevel,2)
  do j=1,iinlevel
     zself%je=j
     inpow=powelt_addr(zself,powcal%powres%npowe)
     zsum=zsum+powcal%powres%pow(inpow)
  end do

  powcal%powres%powa(self%ie)=zsum/iinlevel

end subroutine powelt_avg
!---------------------------------------------------------------------
!> normalised deviation of power on element
subroutine powelt_devn(self,powcal)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure

  !! local
  character(*), parameter :: s_name='powelt_devn' !< subroutine name
  type(powelt_t) :: zself   !< object data structure
  real(kr4) :: zpmax  !< local variable
  real(kr4) :: zpmin  !< local variable
  real(kr4) :: zpest  !< local variable
  real(kr4) :: zdev  !< local variable

  zself%ie=self%ie
  ! find extremal values of pow on element
  ! find number in level from je, although je usually set to iinlevel
  ilevel=powelt_level_table(self%je)
  iinlevel=powelt_table(ilevel,2)
  do j=1,iinlevel
     zself%je=j
     inpow=powelt_addr(zself,powcal%powres%npowe)
     if (j==1) then
        zpmax=powcal%powres%pow(inpow)
        zpmin=powcal%powres%pow(inpow)
     else
        zpmax=max(zpmax,powcal%powres%pow(inpow))
        zpmin=min(zpmin,powcal%powres%pow(inpow))
     end if
  end do

  zpest=(abs(zpmax)+abs(zpmin))/2
  !??? need better test using eps
  if (zpest>0) then
     zdev=(zpmax-zpmin)*50/zpest
  else
     zdev=0
  end if
  powcal%powres%pows(self%ie)=zdev

end subroutine powelt_devn
!---------------------------------------------------------------------
!> average deviation of power on element
subroutine powelt_stat(self,powcal)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure

  !! local
  character(*), parameter :: s_name='powelt_stat' !< subroutine name
  type(powelt_t) :: zself   !< object data structure
  real(kr4) :: zsumn !< sum over higher refinement level
  real(kr4) :: zsumm !< sum over minus one refinement level
  integer(ki4) :: ilevel   !< higher refinement level
  integer(ki4) :: imlevel   !< minus one refinement level
  integer(ki4) :: iinlevel   !< number of subelements at higher level
  integer(ki4) :: iimlevel   !< number of subelements at minus one level

  real(kr4) :: zdev  !< local variable

  zself%ie=self%ie
  ! find extremal values of pow on element
  ! find number in level from je, although je usually set to iinlevel
  ilevel=powelt_level_table(self%je)
  imlevel=max(ilevel-1,1)
  iinlevel=powelt_table(ilevel,2)
  iimlevel=powelt_table(imlevel,2)
  zsumm=0
  do j=1,iimlevel
     zself%je=j
     inpow=powelt_addr(zself,powcal%powres%npowe)
     zsumm=zsumm+powcal%powres%pow(inpow)
  end do
  zsumn=zsumm
  do j=iimlevel+1,iinlevel
     zself%je=j
     inpow=powelt_addr(zself,powcal%powres%npowe)
     zsumn=zsumn+powcal%powres%pow(inpow)
  end do

  zdev=zsumm/iimlevel-zsumn/iinlevel
  powcal%powres%pows(self%ie)=zdev

end subroutine powelt_stat
!---------------------------------------------------------------------
!> write data associated with element
subroutine powelt_write(self,kout)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  integer, intent(in) :: kout   !< input channel for object data structure


  !! local
  character(*), parameter :: s_name='powelt_write' !< subroutine name

  write(kout,*,iostat=status) self
  call log_write_check(m_name,s_name,1,status)

end subroutine powelt_write
!---------------------------------------------------------------------
!> write (vtk) data associated with element
subroutine powelt_writev(self,kplot)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  integer, intent(in) :: kplot   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='powelt_writev' !< subroutine name

  write(kplot,*,iostat=status) self
  call log_write_check(m_name,s_name,1,status)

end subroutine powelt_writev
!---------------------------------------------------------------------
!> address in powcal structure
function powelt_addr(self,kn)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  integer(ki4), intent(in) :: kn   !< number of objects in powcal data structure

  !! local
  character(*), parameter :: s_name='powelt_addr' !< subroutine name
  integer(ki4) :: powelt_addr !< local variable
  integer(ki4) :: iaddr   !< address

  ilevel=powelt_level_table(self%je)
  iaddr=powelt_table(ilevel,1)*kn-powelt_table(ilevel,2)&
 &+powelt_table(ilevel,3)*self%ie+self%je
  powelt_addr=iaddr

end function powelt_addr
!---------------------------------------------------------------------
!> save data needed for power deposition calculation
subroutine powelt_setpow(self,powcal,nobjhit,zm,inpow,gshadl)

  !! arguments
  type(powelt_t), intent(in) :: self   !< object data structure
  type(powcal_t), intent(inout) :: powcal !< powcal data structure
  integer(ki4), intent(in) :: nobjhit   !< object hit
  type(posnode_t), intent(in) :: zm !< final position
  integer(ki4), intent(in) :: inpow   !< index of element
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure

  !! local
  character(*), parameter :: s_name='powelt_setpow' !< subroutine name
  real(kr8) :: zr !< value of \f$ R \f$
  real(kr8) :: zz !< value of \f$ Z \f$
  type(posvecl_t) :: zpos !< final position
  type(posvecl_t) :: zpos1   !< one position data
  real(kr8) :: zpsi !< value of \f$ \psi \f$
  integer(ki4) :: ihit   !< corrected object hit

  ihit=nobjhit
  if (lpmidplane) then
     ! valid exit only if trapped by midplane test (case 'msus','global')
     if (nobjhit==-1.OR.nobjhit==0) ihit=1
  end if

  if (ihit>0.OR. ihit==-10)  then
     ! really a collision
     calcn_type: select case (powcal%n%caltype)
     case('afws','local','msus','global')
        !! shadowed, no power deposited
        powcal%powres%pow(inpow)=0
     case('msum','middle')
        !! hit something, power deposited upon it
        powcal%powres%pow(inpow)=1
        powcal%powres%geobjl%obj(inpow)%weight=ihit
     end select calcn_type
  else if (ihit<0) then
     ! leaving volume
     calcn_type2: select case (powcal%n%caltype)
     case('afws','local')
        !! unshadowed, power to be deposited in this case
        powcal%powres%pow(inpow)=1
     case('msum','middle')
        !! no  power to be deposited in this case (may be error)
        powcal%powres%pow(inpow)=0
     case('msus','global')
        !! power calculated from psi at exit
        field_type: select case (powcal%powres%beq%n%fldspec)
        case(1)
           ! psi is the position (unquantised)
           zpos1%posvec=zm%posvec
           zpos=position_invqtfm(zpos1,powcal%powres%geobjl%quantfm)
           zpsi=zpos%posvec(1)
        case default
           zr=zm%posvec(1)
           zz=zm%posvec(2)
           call spl2d_eval(powcal%powres%beq%psi,zr,zz,zpsi)
        end select field_type
        powcal%powres%pow(inpow)=1
        powcal%powres%geobjl%obj(inpow)%weight=zpsi
     end select calcn_type2
  end if

end subroutine powelt_setpow

end module powelt_m
