module nrsolve_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use misc_m
  use date_time_m
  use spl2d_m
  use pcontrol_h
!     use beq_m, only : beq_ripple_h1

  implicit none
  private

! public subroutines
  public :: &
 &nrsolve_init, &
 &nrsolve_close, &
 &nrsolve_spl2dfn, &
 &nrsolve_funct

! public types
  integer :: nout !< output file unit

! private variables
  character(*), parameter :: m_name='nrsolve_m' !< module name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: idum !< dummy integer
  integer(ki4), parameter :: maxit = 100 !< maximum number of iterations allowed
  real(kr8) , parameter :: cvgnce = 1.e-8_kr8 !<  \f$ \epsilon \f$ for convergence
  real(kr8) , parameter :: dereps = 1.e-8_kr8 !<  \f$ \epsilon \f$ for vanishing derivative
  real(kr8) , parameter :: derepssq = dereps*dereps !<  \f$ \epsilon \f$ squared for vanishing derivative

  contains
!---------------------------------------------------------------------
!> open output file
subroutine nrsolve_init(file,timestamp)

  !! argument
  type(pfiles_t), intent(in) :: file !< file names
  type(date_time_t), intent(in) :: timestamp !< timestamp of run

  !! local
  character(*), parameter :: s_name='nrsolve_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  !! get file do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nout=i exit end if end do

  call misc_getfileunit(nout)
  open(unit=nout,file='nrsolve.log',status='replace')

  !! header information
  write(nout,'(a)') trim(timestamp%long)

end subroutine nrsolve_init
!---------------------------------------------------------------------
!> close output files
subroutine nrsolve_close

  close(unit=nout)

end  subroutine nrsolve_close
!---------------------------------------------------------------------
!> solve for a 2-D spline function
subroutine nrsolve_spl2dfn(spl2d,dspldr,dspldz,pval,pr,pz,kerr,kopt)

  !! arguments
  type(spl2d_t), intent(inout) :: spl2d !< spl2d object data structure (\f$ \psi(R,Z) \f$)
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data structure
  real(kr8), intent(inout) :: pval !< requested value of spl2d
  real(kr8), intent(inout) :: pr !< root in 1-coordinate
  real(kr8), intent(inout) :: pz !< root in 2-coordinate
  integer(ki4), intent(out) :: kerr !< error flag
  integer(ki4), intent(in), optional :: kopt   !< options

  ! local variables
  character(*), parameter :: s_name='nrsolve_spl2dfn' !< subroutine name
  real(kr8) :: zs   !< parameter value
  real(kr8) :: zr   !< 1-coordinate value
  real(kr8) :: zz   !< 2-coordinate value
  real(kr8) :: zp    !<  \f$ \psi(R,Z) \f$
  real(kr8) :: zpold    !<  old \f$ \psi(R,Z) \f$
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zdp    !<  \f$ |{\nabla\psi}| \f$
  real(kr8) :: zdpsq  !<  \f$ |{\nabla\psi}|^2 \f$
  real(kr8), dimension(maxit) :: zra !< local variable
  real(kr8), dimension(maxit) :: zza !< local variable
  real(kr8), dimension(maxit) :: zpa !< local variable
  real(kr8), dimension(maxit) :: zdpa !< local variable
  real(kr8), dimension(maxit) :: zdpna !< local variable
  real(kr8),save :: zh1    !<  normalisation of derivative
  real(kr8),save :: zh2    !<  normalisation of derivative
  integer(ki4),save :: icall=0    !<  call number

  kerr=0
  icall=icall+1
  if(present(kopt)) then
     zh1=pr
     zh2=pz
     return
  end if
  ! initialise
  zs=0
  zr=pr
  zz=pz
  zpold=pval
  !D       write(*,*) 'spl2dfn',zr,zz,zpold    !D

  ! loop
  do i=1,maxit
     call spl2d_evaln(spl2d,zr,zz,1,zp)
     zra(i)=zr
     zza(i)=zz
     zpa(i)=zp
     if (abs(zp-pval)<cvgnce) then
        pr=zr
        pz=zz
        return
     end if

     call spl2d_evaln(dspldr,zr,zz,2,zdpdr)
     call spl2d_evaln(dspldz,zr,zz,2,zdpdz)
     zdpsq=zdpdr**2+zdpdz**2
     !D      if (i==1) write(*,*) icall,zdpsq    !D
     if (zdpsq<derepssq) then
        ! zero derivative condition
        kerr=10
        exit
     end if
     zdp=sqrt(zdpsq)
     zdpa(i)=zdp
     zdpna(i)=zdp
     zs=zs-(zp-pval)/zdp
     zr=pr+zs*(zdpdr/(zdp*zh1))
     zz=pz+zs*(zdpdz/(zdp*zh2))
     zpold=zp
     !D       write(*,*) i,zr,zz,zpold    !D
  end do

  ! error conditions and warn
  if (kerr==0) kerr=1 ! exceeded max number of its.
  !D      write(*,*) 'kerr=',kerr    !D
  call log_error(m_name,s_name,kerr,error_warning,'Convergence failure')

end subroutine nrsolve_spl2dfn
!---------------------------------------------------------------------
!> solve for a 1-D function with parameters
subroutine nrsolve_funct(pfunct,pval,px,kparam,param,kerr)

  !! arguments
  real(kr8), intent(in) :: pval !< requested value of pfunct
  real(kr8), intent(inout) :: px !< root
  integer(ki4), intent(in) :: kparam !< number of parameters
  real(kr8), intent(in), dimension(kparam) :: param !< parameters
  integer(ki4), intent(out) :: kerr !< error flag
  real(kr8), external :: pfunct !< local variable

  ! local variables
  character(*), parameter :: s_name='nrsolve_funct' !< subroutine name
  real(kr8) :: zs   !< parameter value
  real(kr8) :: zx   !< coordinate value
  real(kr8) :: zf    !<  function value
  real(kr8) :: zfold   !<  old function value
  real(kr8) :: zdfdx    !<  \f$ df/dx \f$
  real(kr8) :: z    !< local variable

  ! initialise
  kerr=0
  zs=param(1)
  zx=px
  z=0
  zfold=pval
  !D       write(*,*) 'funct',zs,zx,zfold    !D

  ! loop
  do i=1,maxit
     zf=pfunct(zx,z,zs,0)
     if (abs(zf-zfold)<cvgnce) then
        px=zx
        return
     end if

     zdfdx=pfunct(zx,z,zs,3)
     if (abs(zdfdx)<dereps) then
        ! zero derivative condition
        kerr=10
        exit
     end if
     zx=zx-(zf-pval)/zdfdx
     zfold=zf
     !D       write(*,*) i,zx,zfold    !D
  end do

  ! error conditions and warn
  if (kerr==0) kerr=1 ! exceeded max number of its.
  call log_error(m_name,s_name,kerr,error_warning,'Convergence failure')

end subroutine nrsolve_funct

end module nrsolve_m
