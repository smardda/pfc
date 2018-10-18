module posang_m

  use const_kind_m
  use const_numphys_h
  use position_m
  use posang_h
  use log_m
  use position_h
  use fmesh_h
  use beq_h
  use spl2d_m

  implicit none
  private

! public subroutines
  public ::posang_tfm,& !< convert polar-toroidal to cartesian coordinates
  posang_invtfm,& !< convert cartesian back to polar-toroidal coordinates
  posang_psitfm,&   !< convert polar to flux coordinates
  posang_invpsitfm,& !< convert flux back to polar coordinates
  posang_cartfm,&   !< rotate and translate to geometry cartesian coordinates
  posang_invcartfm,& !< rotate and translate from geometry cartesian coordinates
  posang_units,&   !< convert units of position only
  posang_tfmlis, & !> transform list of positions as posangs
  posang_invtfmlis, & !> inverse transform list of positions as posangs
  posang_writev !< output posang vector

! private variables
  character(*), parameter :: m_name='posang_m' !< module name
  real(kr8), dimension(:), allocatable :: work1 !< 1D work array
  real(kr8) :: zpsi    !<  \f$ \psi \f$
  real(kr8) :: ztheta    !<  \f$ \theta \f$
  real(kr8) :: zeta    !<  \f$ \zeta \f$
  real(kr8) :: zcos    !<  \f$ \cos(\zeta) \f$
  real(kr8) :: zsin    !<  \f$ \sin(\zeta) \f$
  real(kr8) :: zdpdr    !< \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !< \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zr    !<   \f$ R \f$
  real(kr8) :: zz    !<   \f$ Z \f$
  real(kr8) :: zv1    !<  1-component of  vector
  real(kr8) :: zv2    !<  2-component of  vector
  real(kr8) :: zv3    !<  3-component of  vector
  real(kr8) :: zvr    !<  \f$ R \f$-component of  vector
  real(kr8) :: zvzeta    !<  \f$ \zeta \f$-component of  vector
  real(kr8) :: zvz    !<  \f$ Z \f$-component of  vector
  real(kr8) :: zx1    !<  1-component of position vector
  real(kr8) :: zx2    !<  2-component of position vector
  real(kr8) :: zx3    !<  3-component of position vector
  integer   :: status   !< error status
  real(kr8) :: zufac !< units conversion factor, power of ten
  real(kr8) :: ztfmfac !< units conversion to metres (i.e. factor 10**0)
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: iopt !< copy of format option
  character(len=80) :: ibuff !< buffer for input/output
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> convert polar-toroidal to cartesian coordinates
subroutine posang_tfm(self,kunits)

  !! arguments
  type(posang_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(in) :: kunits !< units of result (-3 for millimetres)

  !! local
  character(*), parameter :: s_name='posang_tfm' !< subroutine name

  iopt=self%opt

  zufac=10.**(self%units-kunits)
  position_format: select case (iopt)
  case(1,17) ! position
     zeta=self%pos(3) ; zsin=sin(zeta) ; zcos=cos(zeta)
     zr=self%pos(1)
     zz=self%pos(2)
     zx1=zr*zcos
     zx2=-zr*zsin
     zx3=zz
     self%pos(1)=zufac*zx1 ; self%pos(2)=zufac*zx2 ; self%pos(3)=zufac*zx3
     self%opt=0 ; self%units=kunits

     vector_format: select case (iopt)
     case(17) ! position and vector
        zvr=self%vec(1) ; zvz=self%vec(2) ; zvzeta=self%vec(3)
        zv1=(1/zr)*(zx1*zvr+zx2*zvzeta)
        zv2=(1/zr)*(zx2*zvr-zx1*zvzeta)
        zv3=zvz
        self%vec(1)=zv1 ; self%vec(2)=zv2 ; self%vec(3)=zv3
        self%opt=16
     end select vector_format

  case default
     call log_error(m_name,s_name,1,error_fatal,'Vector in wrong format')
  end select position_format

end subroutine posang_tfm
!---------------------------------------------------------------------
!> convert cartesian back to polar-toroidal coordinates
subroutine posang_invtfm(self,kunits)

  !! arguments
  type(posang_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(in) :: kunits !< units of result (0 for metres)

  !! local
  character(*), parameter :: s_name='posang_invtfm' !< subroutine name

  iopt=self%opt

  zufac=10.**(self%units-kunits)
  position_format: select case (iopt)
  case(0,16) ! position
     zx1=self%pos(1)*zufac ; zx2=self%pos(2)*zufac ; zx3=self%pos(3)*zufac
     zr=sqrt( max(0.,zx1**2+zx2**2) )
     zeta=atan2(-zx2,zx1)
     zz=zx3
     self%pos(1)=zr ; self%pos(2)=zz ; self%pos(3)=zeta
     self%opt=1 ; self%units=kunits

     vector_format: select case (iopt)
     case(16) ! position and vector
        zv1=self%vec(1) ; zv2=self%vec(2) ; zv3=self%vec(3)
        zvr=(-1/zr)*(zx1*zv1+zx2*zv2)
        zvzeta=(1/zr)*(-zx2*zv1+zx1*zv2)
        zvz=zv3
        self%vec(1)=zvr ; self%vec(2)=zvz ; self%vec(3)=zvzeta
        self%opt=17
     end select vector_format

  case default
     call log_error(m_name,s_name,1,error_fatal,'Vector in wrong format')
  end select position_format

end subroutine posang_invtfm
!---------------------------------------------------------------------
!> convert polar to flux coordinates
subroutine posang_psitfm(self,pbeq)

  !! arguments
  type(posang_t), intent(inout) :: self   !< object data structure
  type(beq_t), intent(inout) :: pbeq   !< beq data structure


  !! local
  character(*), parameter :: s_name='posang_psitfm' !< subroutine name

  iopt=self%opt

  position_format: select case (iopt)
  case(1,17) ! position
     zr=self%pos(1) ; zz=self%pos(2) ; zeta=self%pos(3)
     call spl2d_eval(pbeq%psi,zr,zz,zpsi)
     ztheta=atan2(zz-pbeq%n%zcen,zr-pbeq%n%rcen)
     if (ztheta<-const_pid/2) ztheta=2*const_pid+ztheta
     zx1=zpsi
     zx2=ztheta
     zx3=zeta
     self%pos(1)=zx1 ; self%pos(2)=zx2 ; self%pos(3)=zx3
     self%opt=2

     vector_format: select case (iopt)
     case(17) ! position and vector
        call log_error(m_name,s_name,2,error_fatal,'Case not implemented')
        self%vec(1)=zv1 ; self%vec(2)=zv2 ; self%vec(3)=zv3
        self%opt=18
     end select vector_format

  case default
     call log_error(m_name,s_name,1,error_fatal,'Vector in wrong format')
  end select position_format

end subroutine posang_psitfm
!---------------------------------------------------------------------
!> convert flux back to polar coordinates
subroutine posang_invpsitfm(self,pbeq)

  !! arguments
  type(posang_t), intent(inout) :: self   !< object data structure
  type(beq_t), intent(inout) :: pbeq   !< beq data structure


  !! local
  character(*), parameter :: s_name='posang_invpsitfm' !< subroutine name

  iopt=self%opt

  position_format: select case (iopt)
  case(2,18) ! position
     zpsi=self%pos(1) ; ztheta=self%pos(2) ; zeta=self%pos(3)
     call spl2d_evaln(pbeq%r,zpsi,ztheta,1,zr)
     call spl2d_evaln(pbeq%z,zpsi,ztheta,2,zz)
     zx1=zr
     zx2=zz
     zx3=zeta
     self%pos(1)=zx1 ; self%pos(2)=zx2 ; self%pos(3)=zx3
     self%opt=1

     vector_format: select case (iopt)
     case(18) ! position and vector
        call log_error(m_name,s_name,2,error_fatal,'Case not implemented')
        self%vec(1)=zv1 ; self%vec(2)=zv2 ; self%vec(3)=zv3
        self%opt=17
     end select vector_format

  case default
     call log_error(m_name,s_name,1,error_fatal,'Vector in wrong format')
  end select position_format

end subroutine posang_invpsitfm
!---------------------------------------------------------------------
!> rotate and translate to geometry cartesian coordinates
subroutine posang_cartfm(self,tfmdata,kunits)

  !! arguments
  type(posang_t), intent(inout) :: self   !< object data structure
  type(tfmdata_t), intent(in) :: tfmdata   !< data defining transform
  integer(ki4), intent(in) :: kunits !< units of result (-3 for millimetres)

  !! local
  character(*), parameter :: s_name='posang_cartfm' !< subroutine name
  type(posvecl_t) :: zpos !< vector
  type(posvecl_t) :: zpost !< transformed vector
  type(tfmdata_t) :: ztfmdata   !< scratch transform

  iopt=self%opt

  zufac=10.**(self%units-kunits)
  ! convert to metres (units=0) before applying transform
  ztfmfac=10.**self%units
  zufac=zufac/ztfmfac
  position_format: select case (iopt)
  case(0,16) ! position
     zpos%posvec=ztfmfac*self%pos
     zpost=position_tfm(zpos,tfmdata)
     self%pos(1)=zufac*zpost%posvec(1)
     self%pos(2)=zufac*zpost%posvec(2)
     self%pos(3)=zufac*zpost%posvec(3)
     self%opt=0 ; self%units=kunits

     vector_format: select case (iopt)
     case(16) ! position and vector
        zpos%posvec=self%vec
        ztfmdata=tfmdata
        ztfmdata%offset=0.
        zpost=position_tfm(zpos,ztfmdata)
        self%vec(1)=zpost%posvec(1)
        self%vec(2)=zpost%posvec(2)
        self%vec(3)=zpost%posvec(3)
        self%opt=16
     end select vector_format

  case default
     call log_error(m_name,s_name,1,error_fatal,'Vector in wrong format')
  end select position_format

end subroutine posang_cartfm
!---------------------------------------------------------------------
!> rotate and translate from geometry cartesian coordinates
subroutine posang_invcartfm(self,tfmdata,kunits)

  !! arguments
  type(posang_t), intent(inout) :: self   !< object data structure
  type(tfmdata_t), intent(in) :: tfmdata   !< data defining transform
  integer(ki4), intent(in) :: kunits !< units of result (-3 for millimetres)

  !! local
  character(*), parameter :: s_name='posang_invcartfm' !< subroutine name
  type(posvecl_t) :: zpos !< vector
  type(posvecl_t) :: zpost !< transformed vector
  type(tfmdata_t) :: ztfmdata   !< scratch transform

  iopt=self%opt

  zufac=10.**(self%units-kunits)
  ! convert to metres (units=0) before applying transform
  ztfmfac=10.**self%units
  zufac=zufac/ztfmfac
  position_format: select case (iopt)
  case(0,16) ! position
     zpos%posvec=ztfmfac*self%pos
     zpost=position_invtfm(zpos,tfmdata)
     self%pos(1)=zufac*zpost%posvec(1)
     self%pos(2)=zufac*zpost%posvec(2)
     self%pos(3)=zufac*zpost%posvec(3)
     self%opt=0 ; self%units=kunits

     vector_format: select case (iopt)
     case(16) ! position and vector
        zpos%posvec=self%vec
        ztfmdata=tfmdata
        ztfmdata%offset=0.
        zpost=position_invtfm(zpos,ztfmdata)
        self%vec(1)=zpost%posvec(1)
        self%vec(2)=zpost%posvec(2)
        self%vec(3)=zpost%posvec(3)
        self%opt=16
     end select vector_format

  case default
     call log_error(m_name,s_name,1,error_fatal,'Vector in wrong format')
  end select position_format

end subroutine posang_invcartfm
!---------------------------------------------------------------------
!> convert units of position only
subroutine posang_units(self,kunits)

  !! arguments
  type(posang_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(in) :: kunits !< units of result (-3 for millimetres)

  !! local
  character(*), parameter :: s_name='posang_units' !< subroutine name
  type(posvecl_t) :: zpos !< vector

  iopt=self%opt

  zufac=10.**(self%units-kunits)
  position_format: select case (iopt)
  case(0,16) ! cartesians
     zpos%posvec=self%pos
     self%pos(1)=zufac*zpos%posvec(1)
     self%pos(2)=zufac*zpos%posvec(2)
     self%pos(3)=zufac*zpos%posvec(3)
     self%opt=0 ; self%units=kunits
  case(1,17) ! polar-toroidal
     zpos%posvec=self%pos
     self%pos(1)=zufac*zpos%posvec(1)
     self%pos(2)=zufac*zpos%posvec(2)
     self%opt=0 ; self%units=kunits
  case default
     call log_error(m_name,s_name,1,error_fatal,'Position in wrong format')
  end select position_format

end subroutine posang_units
!---------------------------------------------------------------------
!> transform list of positions as posangs
subroutine posang_tfmlis(self,kunits,kzetp)

  !! arguments
  type(posveclis_t), intent(inout) :: self !< posang list data
  integer(ki4), intent(in) :: kunits !< units of result (-3 for millimetres)
  integer(ki4), intent(in) :: kzetp !< reciprocal factor to convert to zeta from xi

  !! local
  character(*), parameter :: s_name='posang_tfmlis' !< subroutine name
  type(posvecl_t) :: zpos !< local variable
  type(posang_t) :: zposang !<local

  !! transform list of positions as posangs
  do j=1,self%np
     zposang%pos = self%pos(j)%posvec ; zposang%opt = 1 ; zposang%units = 0 ! metres, polars
     zposang%pos(3) = zposang%pos(3)/kzetp
     call posang_tfm(zposang,kunits)
     self%pos(j)%posvec=zposang%pos
  end do
  self%nparpos(1)=kunits
  self%nparpos(2)=zposang%opt

end subroutine posang_tfmlis
!---------------------------------------------------------------------
!> inverse transform list of positions as posangs
subroutine posang_invtfmlis(self,kunits,kzetp)

  !! arguments
  type(posveclis_t), intent(inout) :: self !< posang list data
  integer(ki4), intent(in) :: kunits !< units of result (-3 for millimetres)
  integer(ki4), intent(in) :: kzetp !< factor to convert to zeta from xi

  !! local
  character(*), parameter :: s_name='posang_invtfmlis' !< subroutine name
  type(posvecl_t) :: zpos !< local variable
  type(posang_t) :: zposang !<local

  !! transform list of positions as posangs
  do j=1,self%np
     zposang%pos = self%pos(j)%posvec ; zposang%opt = 0 ; zposang%units = -3 ! millimetres, cartesians
     zposang%pos(3) = zposang%pos(3)*kzetp
     call posang_invtfm(zposang,kunits)
     self%pos(j)%posvec=zposang%pos
  end do
  self%nparpos(1)=kunits
  self%nparpos(2)=zposang%opt

end subroutine posang_invtfmlis
!---------------------------------------------------------------------
!> output posang vector
subroutine posang_writev(self,kplot,kopt)

  !! arguments
  type(posang_t), intent(in) :: self   !< posang data
  integer, intent(in) :: kplot   !< output channel for posang data
  integer(ki4), intent(in), optional :: kopt   !< select vector if=2


  !! local
  character(*), parameter :: s_name='posang_writev' !< subroutine name

  if(.NOT.present(kopt)) then

     write(kplot,cfmtbv1,iostat=status) self%pos
     if(status/=0) then
        call log_error(m_name,s_name,1,error_fatal,'Error writing position')
     end if

  else

     if(kopt==2) then
        write(kplot,cfmtbv1,iostat=status) self%vec
        if(status/=0) then
           call log_error(m_name,s_name,2,error_fatal,'Error writing vector')
        end if
     else if(kopt==3) then
        write(kplot,cfmtbv2,iostat=status) self%pos,self%vec
        if(status/=0) then
           call log_error(m_name,s_name,3,error_fatal,'Error writing vector')
        end if
     end if

  end if

end subroutine posang_writev

end module posang_m
