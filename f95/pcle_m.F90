module pcle_m

  use const_kind_m
  use log_m
  use const_numphys_h
  use position_h
  use pcle_h
  use position_m
  use control_h
  use ls_m
  use btree_m
  use geobj_m
  use query_m
  use geobjlist_h
  use termplane_h

  implicit none
  private

! public subroutines
  public :: &
  pcle_movet, &    !< move pcle with terminating plane test
  pcle_move    !< move pcle

! public types

! private types

! private variables
  character(*), parameter :: m_name='pcle_m' !< module name
  integer   :: status   !< error status
  integer(ki4) :: inadr   !< address in position list
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: hitobj !< geobj hit
  integer(ki4) :: inls !< number of entries in list
  integer(ki2) :: ileaf !< leaf test
  logical, parameter :: debug=.FALSE. !< logical flag for debugging

  contains
!---------------------------------------------------------------------
!> move particle with terminating plane test
subroutine pcle_movet(selfo,selfn,kstep,geol,btree,termp,kobj)
  !! arguments
  type(posnode_t), intent(in) :: selfo !< input pcle posn/node
  type(posnode_t), intent(inout) :: selfn !< output pcle posn/node
  integer(ki4), intent(in) :: kstep   !< step number
  type(geobjlist_t), intent(in) :: geol !< geobj list data
  type(btree_t), intent(in) :: btree   !< binary tree data
  type(termplane_t), intent(inout) :: termp   !< termination plane data
  integer(ki4), intent(out) :: kobj   !< hit object

  !! local
  character(*), parameter :: s_name='pcle_movet' !< subroutine name
  integer(ki2par) :: igcode   !< code for geometry
  real(kr8) :: zf1 !< fraction of path length
  real(kr8) :: zf1t !< fraction of path length (numerator)
  real(kr8) :: zf1b !< fraction of path length (denominator)
  real(kr8) :: zfmin !< smallest fraction of path length
  integer(ki4) :: idir !< coordinate direction
  integer(ki4) :: idirc !< coordinate direction to apply condition on coordinate
  integer(ki4) :: jt !< loop over termination planes
  integer(ki4) :: idirs !< sign of coordinate direction
  integer(ki4) :: idircs !< sign of coordinate direction to apply condition on coordinate
  integer(ki4) :: iop !< code describing termination type

  call pcle_move(selfo,selfn,kstep,geol,btree,kobj)

  igcode=0
  if (kobj>0) igcode=geobj_code_fn(geol%obj2(kobj)%typ)
  if (igcode>0) then
     kobj=GEOBJ_OFFSET-1-igcode
     if (igcode>=GEOBJ_SKYLIT) return
  end if

  if (termp%ntermplane==0) return

  ! test on path which has been reduced if there is a collision
  zfmin=1.1
  do jt=1,termp%ntermactive
     idir=termp%termplanedir(jt,1)
     idirs=termp%termplanedir(jt,2)
     iop=termp%termplanedir(jt,3)
     idirc=termp%termplanedir(jt,4)
     idircs=termp%termplanedir(jt,5)
     terminate_type: select case (iop)
     case (0)
        ! terminate if goes past a plane
        if ((selfn%posvec(idir)-termp%termplane(jt,1))*idirs>0) then
           !rewritten to avoid occasional division by zero
           !zf1=(termp%termplane(jt,1)-selfo%posvec(idir))&
           !&         /(selfn%posvec(idir)-selfo%posvec(idir))
           !if (zf1<=zfmin) then
           !zfmin=zf1
           zf1t=termp%termplane(jt,1)-selfo%posvec(idir)
           zf1b=selfn%posvec(idir)-selfo%posvec(idir)
           if (abs(zf1t)<=zfmin*abs(zf1b)) then
              zfmin=zf1t/zf1b
              kobj=-2
              if (debug) then
                 if ( jt==termp%ntermactive ) write(801,*) 'termplane n'
              end if
           end if
        end if
     case (1)
        ! terminate if intersects a plane, idirs irrelevant
        if ((selfn%posvec(idir)-termp%termplane(jt,1))*(selfo%posvec(idir)-termp%termplane(jt,1))<=0) then
           !rewritten to avoid occasional division by zero
           !zf1=(termp%termplane(jt,1)-selfo%posvec(idir))/(selfn%posvec(idir)-selfo%posvec(idir))
           !if (zf1<=zfmin) then
           !zfmin=zf1
           zf1t=termp%termplane(jt,1)-selfo%posvec(idir)
           zf1b=selfn%posvec(idir)-selfo%posvec(idir)
           if (abs(zf1t)<=zfmin*abs(zf1b)) then
              zfmin=zf1t/zf1b
              kobj=-2
           end if
        end if
     case (2)
        if (termp%termstore(1,idir)<1.e29) then
           ! terminate if reaches an extremum
           if ((selfo%posvec(idir)-termp%termstore(1,idir))*(selfo%posvec(idir)-selfn%posvec(idir))>=0) then
              if (idirs*(selfo%posvec(idir)-selfn%posvec(idir))>=0) then
                 ! must be specific type
                 zf1=0
                 if (zf1<=zfmin) then
                    zfmin=zf1
                    kobj=-2
                 end if
              end if
           end if
        else
           ! save first point
           termp%termstore(1,:)=selfo%posvec
        end if
     case (3)
        idirc=termp%termplanedir(jt,4)
        idircs=termp%termplanedir(jt,5)
        ! terminate if intersects a plane, satisfying auxilliary condition
        if ((selfn%posvec(idirc)-termp%termplane(jt,2))*idircs>0) then
           if ((selfn%posvec(idir)-termp%termplane(jt,1))*(selfo%posvec(idir)-termp%termplane(jt,1))<=0) then
              !rewritten to avoid occasional division by zero
              !zf1=(termp%termplane(jt,1)-selfo%posvec(idir))/(selfn%posvec(idir)-selfo%posvec(idir))
              !if (zf1<=zfmin) then
              !zfmin=zf1
              zf1t=termp%termplane(jt,1)-selfo%posvec(idir)
              zf1b=selfn%posvec(idir)-selfo%posvec(idir)
              if (abs(zf1t)<=zfmin*abs(zf1b)) then
                 zfmin=zf1t/zf1b
                 kobj=-2
              end if
           end if
        end if
     end select terminate_type
  end do

  if (kobj==-2) then
     selfn%posvec=(1-zfmin)*selfo%posvec+zfmin*selfn%posvec
     selfn%node=0
  end if

end subroutine pcle_movet
!---------------------------------------------------------------------
!> move particle
subroutine pcle_move(selfo,selfn,kstep,geol,btree,kobj)

  !! arguments
  type(posnode_t), intent(in) :: selfo !< input pcle posn/node
  type(posnode_t), intent(inout) :: selfn !< output pcle posn/node
  integer(ki4), intent(in) :: kstep   !< step number
  type(geobjlist_t), intent(in) :: geol !< geobj list data
  type(btree_t), intent(in) :: btree   !< binary tree data
  integer(ki4), intent(out) :: kobj   !< hit object

  !! local
  character(*), parameter :: s_name='pcle_move' !< subroutine name
  logical :: lcoll   !< collision on path
  type(geobj_t) :: iobj !< local variable
  type(posveclis_t) :: zposl   !< list of position data
  integer(ki4) :: jj !< loop counter
  real(kr4) :: zd1  !< local variable
  real(kr4) :: zdel !< local variable
  real(kr4) :: zdela !< local variable
  real(kr4) :: zdeln  !< local variable
  real(kr4) :: zalf1 !< local variable
  real(kr4) :: zalf2 !< local variable
  real(kr4) :: zalf3  !< local variable
  real(kr4), dimension(3) :: zalf  !< local variable
  type(posvecl_t) :: zposo !< local variable
  type(posvecl_t) :: zposn !< local variable
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zzpos !< local variable
  type(posnode_t) :: zpcle !< local variable
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: ialf !<  integer direction vector
  integer(ki2) :: ialf1 !< local variable
  integer(ki2) :: iexten1  !< local variable
  integer(ki2) :: ishj !< local variable
  integer(ki2) :: idifj  !< local variable
  integer(ki2) :: isgalf  !< local variable
  integer(ki2) :: d1 !< local variable
  integer(ki2) :: d2 !< local variable
  integer(ki2) :: d3  !< local variable
  integer(ki2) :: i1 !< local variable
  integer(ki2) :: i2 !< local variable
  integer(ki2) :: i3  !< local variable
  integer(ki2) :: il2 !< local variable
  integer(ki2) :: il3  !< local variable
  integer(ki2) :: i1o !< loop 1 variables
  integer(ki2) :: i1r !< loop 1 variables
  integer(ki2) :: i1j !< loop 1 variables
  integer(ki2) :: i1n !< 2 and 3 not needed
  integer(ki2) :: i2n !< 2 and 3 not needed
  integer(ki2) :: i3n !< 2 and 3 not needed
  integer(ki2) :: i1d !< integer coordinates (d=departure)
  integer(ki2) :: i2d !< integer coordinates (d=departure)
  integer(ki2) :: i3d !< integer coordinates (d=departure)
  integer(ki4) :: ipow2  !< local variable
  real(kr4), dimension(3) :: zvec !< local variable
  real(kr4), dimension(3) :: zveco !< local variable
  real(kr4), dimension(3) :: zvecn !< local variable
  real(kr4), dimension(3) :: zvdif !< local variable
  real(kr4), dimension(3) :: dela  !< local variable
  integer(ki2), dimension(3) :: ivec!< integer geobj data
  integer(ki2), dimension(3) :: iveco!< integer geobj data
  integer(ki2), dimension(3) :: ivecn!< integer geobj data
  integer(ki2), dimension(3) :: ivecc !< integer geobj data
  integer(ki2), dimension(3) :: icorn !< new corner
  integer(ki2), dimension(3) :: icornxt !< next corner
  integer(ki2), dimension(3) :: iea !< array of ie
  integer(ki2) :: iesum !< integer sum
  integer(ki4) :: iupnode !< tree node number
  integer(ki4) :: inodn !< tree node number
  integer(ki4) :: iinode !< tree node number
  integer(ki4) :: inode !< tree node number
  integer(ki4) :: inod2 !< tree node number
  integer(ki4) :: inod3 !< tree node number

  logical, parameter :: iljump=.TRUE. !< allow jumps in 1 direction
  logical, parameter :: ilrawpos=.false. !< positions to be quantised

  ! set defaults
  hitobj=0
  inodn=0
  lcoll=.false.

  ! local copies of positions
  zposo%posvec=selfo%posvec
  zposn%posvec=selfn%posvec
  if (ilrawpos) then
     ! quantise particle position
     zzpos=position_qtfm(zposo,geol%quantfm)
     zposo=zzpos
     zzpos=position_qtfm(zposn,geol%quantfm)
     zposn=zzpos
  end if

  inode=selfo%node
  if (kstep==0.OR.inode==0) then
     ! find inode, starting node, if not defined or on first step
     ! calculate inode from position
     status=0
     if (.not.allocated(zposl%pos)) allocate(zposl%pos(1), stat=status)
     !! check successful allocation
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Unable to allocate memory')
     end if
     iobj%objtyp=1
     iobj%geobj=1
     zposl%pos(1)=zposo
     call btree_mfind(btree,iobj,zposl,inode)
  end if
  ! check
  if (inode<=0) then
     call log_error(m_name,s_name,2,error_warning,'Bad particle node ')
     return
  end if

  !! main working
  ! test for collision in current node
  zpcle%node=inode
  call pcle_collinnode(zposo,zposn,geol,btree,zpcle,hitobj)
  lcoll=(hitobj/=0)
  if (lcoll) goto 200

  ! Evaluate iesum test which compares old and new positions versus node
  zveco=zposo%posvec
  iveco=int(zveco)
  zvecn=zposn%posvec
  ivecn=int(zvecn)
  iext=btree%exten(:,btree%desc(2,inode))
  do jj=1,3
     iea(jj)=ishft( ieor(iveco(jj),ivecn(jj)), -iext(jj) )
  end do
  iesum=iea(1)+iea(2)+iea(3)

  ! three way split depending on iesum
  if (iesum==0) then
     !! first branch
     goto 200

     !! second branch
  else if (iesum==1) then
     ! test adjacent node
     iupnode=btree%pter(1,inode)
     if (iupnode>1) then
        inodn=inode
        ipow2=1
        do jj=1,3
           idifj=ivecn(jj)-iveco(jj)
           ishj=iea(jj)  * ipow2
           inodn=inodn+sign(ishj,idifj)
           ipow2=2*ipow2
        end do
        ileaf=btree%pter(3,inodn)
        if (ileaf>0) then
           ! not lowest node, fall through and track particle
        else if (ileaf==-1) then
           ! test for collision in adjacent node
           zpcle%node=inodn
           call pcle_collinnode(zposo,zposn,geol,btree,zpcle,hitobj)
           lcoll=(hitobj/=0)
           goto 200
        else
           ! empty node
           zpcle%node=inodn
           goto 200
        end if
     end if
  end if

  !! Logically the third branch of iesum
  ! follow particle path
  ! find direction in which path is longest and call it d1
  zvdif=zvecn-zveco
  d1=maxloc( (/abs(zvdif(1)),abs(zvdif(2)),abs(zvdif(3))/) , dim=1)
  d2=mod(d1,3)+1
  d3=mod(d2,3)+1
  ! initialise loop over d1
  zalf(d1)=sign( 1._kr4,zvecn(d1)-zveco(d1) )
  zalf(d2)=(zvecn(d2)-zveco(d2))/abs(zvecn(d1)-zveco(d1))
  zalf(d3)=(zvecn(d3)-zveco(d3))/abs(zvecn(d1)-zveco(d1))
  zalf(d2)=sign( max(const_pusheps,abs(zalf(d2))), zalf(d2) )
  zalf(d3)=sign( max(const_pusheps,abs(zalf(d3))), zalf(d3) )
  zalf1=zalf(d1)
  zalf2=zalf(d2)
  zalf3=zalf(d3)
  ! ialf entries are zero for negative i1 travel
  do jj=1,3
     ialf(jj)=max( min( ceiling(zalf(jj)) ,1 ) , 0 )
  end do
  ialf1=ialf(d1)
  isgalf=2*ialf1-1
  ! initial and final integer values of particle position
  ! also gives start and end values of i1
  ivec=iveco
  i1d=1-ialf1+iveco(d1)
  ivec(d1)=i1d
  i2d=iveco(d2)
  i3d=iveco(d3)
  i1n=ivecn(d1)+1-ialf1
  i2n=ivecn(d2)
  i3n=ivecn(d3)
  i1o=i1d
  ! virtual start position for path, at cell wall in 1 direction
  zvec(d1)=real(i1d)
  zd1=abs(zveco(d1)-zvec(d1))
  zvec(d2)=zveco(d2)-zalf2*zd1
  zvec(d3)=zveco(d3)-zalf3*zd1

  !! optimised loop, track along 1 direction, watch for motion in 2 and 3
  ! using il2 and il3
  !            loop_path: do i1=i1o,i1n,isgalf
  i1=i1o

199   continue

  iea=(/0,0,0/)
  icorn=btree%corner(:,inode)
  iext=btree%exten(:,btree%desc(2,inode))

  if (iljump) then
     ! option to advance in 1 direction as far as possible without leaving box
     ! ie. jump more than one cell in 1 per step
     ! Estimate max length in d1 direction (less than distance to end)
     iexten1=ishft( 1,iext(d1) )
     zdeln=abs(zvecn(d1)-zvec(d1))
     do jj=1,3
        ! corners in direction of travel
        icornxt(jj)=icorn(jj)+ishft(ialf(jj),iext(jj))
        dela(jj)=(real(icornxt(jj))-zvec(jj))/zalf(jj)
     end do
     ! dela should always be positive
     zdela=min(dela(1),dela(2),dela(3))
     !D          write(*,*) 'dela',dela(1),dela(2),dela(3),zdeln! writediagn
     !D        if (zdeln<120.) then! writediagn
     !D          write(*,*) 'getting near'! writediagn
     !D        end if! writediagn
     zdel=min(zdela,zdeln)
     i1r=int(zdel)
     i1j=max( 0, min( i1r,iexten1-1 ) )
     i1=i1+isgalf*i1j
     zvec(d2)=zvec(d2)+zalf2*i1j
     i2d=int(zvec(d2))
     zvec(d3)=zvec(d3)+zalf3*i1j
     i3d=int(zvec(d3))
     ivec=int(zvec)
     ivec(d1)=i1
     ! next i1 after this should now generate point(s) in new box
  end if
  !D       write(*,*) icorn,iext! writediagn
  !D       if (iljump) write(*,*) i1,zvec(d2),zvec(d3)! writediagn
  !D       if (.not.iljump) write(*,*) ivec,zvec! writediagn


  ! check in 2 direction
  i2=i2d
  zvec(d2)=zvec(d2)+zalf2
  i2d=int(zvec(d2))
  il2=abs(i2-i2d)
  inod2=inode
  ivecc=ivec
  inodn=inod2
  if (il2>=1) then
     ! test i1,i2d,i3 still in same box
     icorn=btree%corner(:,inod2)
     iext=btree%exten(:,btree%desc(2,inod2))
     ivecc(d2)=i2d
     iea(d2)=ishft( ieor(ivecc(d2),icorn(d2)), -iext(d2) )
     if (iea(d2)/=0) then
        ! find new node
        inod2=btree_newnode(btree,ivecc,inode)
        zpcle%node=inod2
        call pcle_collinnode(zposo,zposn,geol,btree,zpcle,hitobj)
        lcoll=(hitobj/=0)
        if (lcoll) goto 200
        inodn=inod2
     end if
  end if

  ! check in 3 direction
  i3=i3d
  zvec(d3)=zvec(d3)+zalf3
  i3d=int(zvec(d3))
  il3=abs(i3-i3d)
  inod3=inode
  if (il3>=1) then
     ! test i1,i2,i3d still in same box
     icorn=btree%corner(:,inod3)
     iext=btree%exten(:,btree%desc(2,inod3))
     ivecc=ivec
     ivecc(d3)=i3d
     iea(d3)=ishft( ieor(ivecc(d3),icorn(d3)), -iext(d3) )
     if (iea(d3)/=0) then
        ! find new node
        inod3=btree_newnode(btree,ivecc,inode)
        zpcle%node=inod3
        call pcle_collinnode(zposo,zposn,geol,btree,zpcle,hitobj)
        lcoll=(hitobj/=0)
        if (lcoll) goto 200
        inodn=inod3

        if (iea(d2)/=0) then
           ! check in both 2 and 3 directions
           ! test i1,i2d,i3d versus box i1,i2d,i3
           ivecc=ivec
           ivecc(d2)=i2d
           ivecc(d3)=i3d
           ! find new node
           inodn=btree_newnode(btree,ivecc,inod2)
           zpcle%node=inodn
           call pcle_collinnode(zposo,zposn,geol,btree,zpcle,hitobj)
           lcoll=(hitobj/=0)
           if (lcoll) goto 200
        end if
     else
        ivecc(d2)=i2d
     end if
  end if

  ! now move in 1 direction
  ivec=ivecc
  i1=i1+isgalf
  ivec(d1)=i1
  zvec(d1)=real(i1)
  ! Only now are ivec and zvec at a position on the track
  iea(d1)=ishft( ieor(ivec(d1),icorn(d1)), -iext(d1) )
  if (iea(d1)/=0) then
     ! find new node
     iinode=inodn
     inodn=btree_newnode(btree,ivec,iinode)
     zpcle%node=inodn
     call pcle_collinnode(zposo,zposn,geol,btree,zpcle,hitobj)
     lcoll=(hitobj/=0)
     if (lcoll) goto 200
  end if
  inode=inodn

  if ( (i1-i1n)*isgalf >= 0 ) then
     ! test whether in same box as final position
     iext=btree%exten(:,btree%desc(2,inode))
     iesum=0
     do jj=1,3
        iea(jj)=ishft( ieor(ivec(jj),ivecn(jj)), -iext(jj) )
        iesum=iesum+iea(jj)
     end do

     if (iesum==0) then
        ! find new node and quit
        inodn=btree_newnode(btree,ivec,inode)
        zpcle%node=inodn
        call pcle_collinnode(zposo,zposn,geol,btree,zpcle,hitobj)
        lcoll=(hitobj/=0)
        goto 200
     else if ( (i1-i1n)*isgalf > 0 ) then
        ! have overshot
        ! leave in right node
        inodn=btree_newnode(btree,ivecn,inode)
        zpcle%node=inodn
        lcoll=.FALSE.
        goto 200
     end if

  end if
  !           end do loop_path
  goto 199

200   continue

  ! tidy up
  if (lcoll) then
     zpos%posvec=zpcle%posvec
     zzpos=zpos
     if (ilrawpos) then
        ! return new particle position to real space
        zzpos=position_invqtfm(zpos,geol%quantfm)
     end if
     zpcle%posvec=zzpos%posvec
  else
     zpcle%posvec=selfn%posvec
  end if

  ! copy locals to output vars
  selfn=zpcle
  kobj=hitobj

end subroutine pcle_move
!---------------------------------------------------------------------
!> test for collision in current node
subroutine pcle_collinnode(poso,posn,geol,btree,pcle,kobj)
  !! arguments
  type(posvecl_t), intent(in) :: poso !< input pcle posn
  type(posvecl_t), intent(inout) :: posn !< output pcle posn
  type(geobjlist_t), intent(in) :: geol !< geobj list data
  type(btree_t), intent(in) :: btree   !< binary tree data
  type(posnode_t), intent(inout) :: pcle !< new particle
  integer(ki4), intent(out) :: kobj   !< hit object

  !! local
  character(*), parameter :: s_name='pcle_collinnode' !< subroutine name
  type(geobj_t) :: iobj !< local variable
  type(posvecl_t) :: zzpos !< local variable
  real(kr4) :: zs !< path length
  real(kr4) :: zsmin !< path length
  integer(ki4) :: inode  !< local variable
  integer(ki4) :: indx  !< local variable
  integer(ki2), dimension(3) :: ivecn !< integer geobj data
  integer(ki2), dimension(3) :: icorn !< new corner
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2) :: ie !< integer element
  integer(ki2) :: iesum !< integer sum
  real(kr4), dimension(3) :: zvdif !< local variable
  integer(ki2) :: ib !< component with largest path difference
  logical :: lhit !< flag hit
  logical :: ilhit !< flag hit

  kobj=0
  inode=pcle%node
  ! special start
  if (inode<0) then
     kobj=-1
     pcle%posvec=posn%posvec
     return
  end if
  ! special end
  ileaf=btree%pter(3,inode)
  if (ileaf==-2) then
     ! empty node, no collision
  else
     inadr=btree%pter(2,inode)
     inls=btree%objectls%list(inadr,2)
     zsmin=const_pushinf
     lhit=.false.
     ! test against all objects
     do l=1,inls
        !        if (l==13) then! writediagn
        !            write(*,*) 'getting near'! writediagn
        !        end if! writediagn
        indx=btree%objectls%list(inadr+l,2)
        iobj%objtyp=geol%obj2(indx)%typ
        iobj%geobj=geol%obj2(indx)%ptr
        call geobj_linehits(iobj,poso,posn,geol%posl,geol%nodl,zzpos,zs,zsmin,ilhit)
        if (ilhit) then
           zsmin=zs
           pcle%posvec=zzpos%posvec
           kobj=indx
           lhit=.true.
        end if
     end do

     !          write(*,*) 'node ',inode! writediagn
     !          write(*,*) btree%corner(:,inode)! writediagn
     !          write(*,*) (btree%objectls%list(inadr+l,2),l=1,inls)! writediagn

     if (lhit) then
        ! check nearest collision takes place in range
        zvdif=posn%posvec-poso%posvec
        ib=maxloc( (/abs(zvdif(1)),abs(zvdif(2)),abs(zvdif(3))/) , dim=1)
        if ((posn%posvec(ib)-pcle%posvec(ib))*(pcle%posvec(ib)-poso%posvec(ib))>0) then
           ! check nearest collision takes place in node
           ! iesum test
           icorn=btree%corner(:,inode)
           ivecn=int(pcle%posvec)
           iext=btree%exten(:,btree%desc(2,inode))
           iesum=0
           do j=1,3
              ie=ishft( ieor(icorn(j),ivecn(j)), -iext(j) )
              iesum=iesum+ie
           end do
           !     kobj=0 if iesum/=0
           kobj=max(1-iesum,0)*kobj
        else
           kobj=0
        end if
     end if
  end if

end subroutine pcle_collinnode

end module pcle_m
