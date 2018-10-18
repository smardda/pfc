module odes_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use spl2d_m
  use spl3d_m
  use odes_h
  use nrsolve_m

  implicit none
  private

! public subroutines
  public :: &
  odes_init, &   !< initialise  odes data structure
  odes_readcon, &   !< namelist input of controls
  odes_rjstep, &   !< special timestepping for \f$ R/J \f$
  odes_1ststepcon, & !< phi-stepping for axisymm 3-D vector field
  odes_2ndstep, &   !< special timestepping for 2nd scheme
  odes_3rdstep, & !< special timestepping for 3-D vector field
  odes_4thstep, & !< special timestepping for axisymm 3-D vector field
  odes_4thstepcon, & !< extraspecial timestepping for axisymm 3-D vector field
  odes_rjstep1, &   !< special first timestep for \f$ R/J \f$
  odes_1ststep1, &   !< special first timestep for phi-stepping axisymm 3-D vector field
  odes_2ndstep1, &   !< special first timestep for 2nd scheme
  odes_3rdstep1, & !< special first timestep for 3-D vector field
  odes_4thstep1, & !< special first timestep for axisymm 3-D vector field
  odes_writev, &   !< write (vtk) odes data (dummy)
  odes_write, &  !< output odes data (dummy)
  odes_delete  !< delete odes data

! public types

! public variables

! private types

! private variables
  character(*), parameter :: m_name='odes_m' !< module name
  real(kr8), dimension(:), allocatable :: work1 !< 1D work array
  integer   :: status   !< error status
!integer :: nin   !< input channel for ZZZZ data
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: idum !< dummy integer
  integer(ki4), parameter  :: mpdim=2 !< dimension of system
  integer(ki4), parameter  :: mpvdim=3 !< dimension of 3rd system
  integer(ki4), parameter  :: ipmaxhalf=10 !< maximum number of halvings
  integer(ki4), parameter  :: ipsiconst=0 !< enforce psi=const. (if non-zero)
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> initialise  odes data structure
subroutine odes_init(self,numerics)

  !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  type(onumerics_t), intent(inout) :: numerics   !< object control data structure

  !! local
  character(*), parameter :: s_name='odes_init' !< subroutine name
  real(kr8), parameter :: zepsmach=1.e-20 !< machine  \f$ \epsilon \f$
  real(kr8), parameter :: zepsrmin=1.e-12 + 2*zepsmach !< min permitted value of epsr
  real(kr8), parameter :: zepsbig=33*zepsmach !< no. significantly bigger than epsmach
  integer(ki4), parameter  :: ifacup=5 !< magic factor for timestep control
  integer(ki4), parameter  :: ifacdn=10 !< magic factor for timestep control
  real(kr8), parameter :: zffac=0.9 !< magic factor for timestep control
  integer(ki4), parameter  :: ipdim=1 !< dimension of system
  !ZAD     integer(ki4), parameter  :: isord=3 !< max scheme order !ZAD
  integer(ki4), parameter  :: isord=2 !< max scheme order !DOD

  self%epsmach=zepsmach
  self%epsrmin=zepsrmin
  self%epsbig=zepsbig
  self%nfacup=ifacup
  self%nfacdn=ifacdn
  self%ffac=zffac
  ! is overwritten below
  self%npdim=ipdim
  self%nsord=isord

  self%facu=ifacup
  self%facut=(zffac/self%facu)**isord
  self%facd=1._kr8/ifacdn
  self%facdt=(zffac/self%facd)**isord
  self%rsord=1._kr8/isord
  self%rsordb=1._kr8/(6*isord)

  ! initialise control structure
  self%n=numerics

  !! set up dependent variables
  self%n%epsr=max(self%epsrmin,self%n%epsr)
  self%reps=2/self%n%epsr
  self%epsar=2*self%n%epsa/self%n%epsr
  self%dtmin=self%epsbig*self%n%dt0
  self%npdim=max(self%n%npdim,ipdim)
  !! array for old difference
  allocate(self%g3do(self%npdim), stat=status)
  call log_alloc_check(m_name,s_name,9,status)
  !! allocate ODE position dimension
  allocate(self%poso(self%npdim), stat=status)
  call log_alloc_check(m_name,s_name,10,status)
  allocate(self%posa(self%npdim), stat=status)
  call log_alloc_check(m_name,s_name,11,status)
  allocate(self%posb(self%npdim), stat=status)
  call log_alloc_check(m_name,s_name,12,status)
  !! allocate position vector storage
  allocate(self%vecp%pos(self%n%stepmax), stat=status)
  call log_alloc_check(m_name,s_name,20,status)
  allocate(self%posk(self%n%stepmax), stat=status)
  call log_alloc_check(m_name,s_name,21,status)
  self%vecp%np=0

end subroutine odes_init
!---------------------------------------------------------------------
!> namelist input of controls
subroutine odes_readcon(selfn,kin)

  !! arguments
  type(onumerics_t), intent(out) :: selfn   !< object numerics data structure
  integer, intent(in),optional :: kin   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='odes_readcon' !< subroutine name
  real(kr8) :: rel_error !< relative error
  real(kr8) :: abs_error !< absolute error
  real(kr8) :: initial_dt !< initial timestep
  integer(ki4) :: max_numsteps   !< max number of steps
  real(kr8) :: max_zeta !< maximum allowed change in \f$ \zeta \f$ (radians)
  integer(ki4) :: system_order !< system order
  integer(ki4) :: start_control !< start control
  integer(ki4) :: termination_control !< termination control
  real(kr8), dimension(3)  :: termination_parameters !< termination parameters

  !! odes parameters
  namelist /odesparameters/ &
 &system_order, &
 &start_control, termination_control, termination_parameters, &
 &rel_error, abs_error, initial_dt, max_numsteps, max_zeta

  !! set default odes parameters
  system_order=0
  start_control=0
  termination_control=0
  termination_parameters=0._kr8
  termination_parameters(1)=1
  termination_parameters(2)=1000000
  termination_parameters(3)=1000000
  rel_error=.0001_kr8
  abs_error=.0001_kr8
  initial_dt=0.1_kr8
  max_numsteps=100
  max_zeta=1.0_kr8

  !!read odes parameters
  read(kin,nml=odesparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_warning,'Error reading odes parameters')
     call log_getunit(ilog)
     write(ilog,nml=odesparameters)
     print '("Error reading odes parameters")'
  end if

  !! check for valid data
  if(system_order<0.OR.system_order>10) &
 &call log_error(m_name,s_name,2,error_fatal,'system_order must be >=0 and <=10')
  if(termination_control<-1.OR.termination_control>10) &
 &call log_error(m_name,s_name,3,error_fatal,'termination_control must be >=-1 and <=10')
  if(start_control<0.OR.start_control>10) &
 &call log_error(m_name,s_name,4,error_fatal,'start_control must be >=0 and <=10')
  if(rel_error<0.OR.rel_error>1) &
 &call log_error(m_name,s_name,11,error_fatal,'rel_error must be >=0 and <=1')
  if(abs_error<0.OR.abs_error>1) &
 &call log_error(m_name,s_name,12,error_fatal,'abs_error must be >=0 and <=1')
  if(initial_dt<=0) &
 &call log_error(m_name,s_name,13,error_fatal,'initial_dt must be >0')
  if(max_numsteps<1) &
 &call log_error(m_name,s_name,14,error_fatal,'max_numsteps must be positive')
  if(max_zeta<0) &
 &call log_error(m_name,s_name,15,error_fatal,'max_zeta must be >=0')
  if(max_zeta>3*const_pid) &
 &call log_error(m_name,s_name,16,error_warning,'max_zeta >3pi')

  !! store values
  selfn%npdim=system_order
  selfn%nstartcon=start_control
  selfn%ntermcon=termination_control
  selfn%termcon=termination_parameters
  selfn%epsr=rel_error
  selfn%epsa=abs_error
  selfn%dt0=initial_dt
  selfn%stepmax=max_numsteps
  selfn%tmax=max_zeta

end subroutine odes_readcon
!---------------------------------------------------------------------
!> special timestepping for \f$ R/J \f$
subroutine odes_rjstep(self,kerr,pf,pdfac,rjspl2d)

  !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  real(kr8), intent(inout) :: pf   !< \f$ f \f$
  real(kr8), intent(in) :: pdfac   !< direction factor
  type(spl2d_t), intent(inout) :: rjspl2d   !< \f$ R/J \f$
  integer(ki4), intent(inout) :: kerr   !< error flag

  !! local
  character(*), parameter :: s_name='odes_rjstep' !< subroutine name
  real(kr8) :: g2 !< local variable
  real(kr8) :: g3 !< local variable
  real(kr8), dimension(3) :: gdia !< local variable
  real(kr8) :: g3dn !< local variable
  real(kr8) :: g23n !< local variable
  real(kr8) :: rs !< local variable
  real(kr8) :: zhf !< local variable
  real(kr8) :: zdt !< local variable
  real(kr8) :: psi !< local variable
  real(kr8) :: g !< local variable
  real(kr8) :: t !< local variable

  psi=self%vecp%pos(self%ndt)%posvec(1)
  g=self%vecp%pos(self%ndt)%posvec(2)
  t=self%vecp%pos(self%ndt)%posvec(3)
201   continue
  call odes_rkf23n(odes_rjfunct,t,self%dt,g,g2,g3,gdia,psi,pf*pdfac,rjspl2d)
  g3dn=abs(g3-g2)
  g23n=abs(g)+abs(g3)+self%epsar
  rs=self%reps*g3dn/g23n

  if (rs<=1) then
     !     step succeeded, but if possible increase t-step for next step
     zhf=self%facu
     if (rs>self%facut) zhf=self%ffac/rs**self%rsord
     if (kerr==-1) zhf=1
     zdt=max(zhf*self%dt,self%dtmin)
     t=t+sign(1._kr8,pdfac)*self%dt
     self%ndt=self%ndt+1
     !     save position
     self%vecp%pos(self%ndt)%posvec(1)=psi
     self%vecp%pos(self%ndt)%posvec(2)=g3
     self%vecp%pos(self%ndt)%posvec(3)=t
     !D       write(*,*) i,t,self%dt,abs(g) !D
     !D       write(20,*)i,t,self%dt,g,g2,g3 !D

     self%posa(1)=min(g2,g3,gdia(1),gdia(2),gdia(3))
     self%posb(1)=max(g2,g3,gdia(1),gdia(2),gdia(3))
     self%dt=zdt
     ! cancel error flag (-1)
     kerr=0
     !        g=g3
     !        if (zg3>self%glt) g=g3-2.*self%glt
     return
  end if
  !
  !     step failed, reduce t-step , equivalent to setting lfail=t
  kerr=-1

  zhf=self%facd
  if (rs<self%facdt) zhf=self%ffac/rs**self%rsord
  zdt=zhf*self%dt
  self%dt=zdt
  !D        write(*,*) 'failure',i,self%dt,g,g2,g3 !D
  !     try again unless t-step too small
  if (self%dt>=self%dtmin) goto 201

  call log_error(m_name,s_name,1,error_warning,'t-step too small')
  kerr=1

end subroutine odes_rjstep
!---------------------------------------------------------------------
!> phi-stepping for axisymm 3-D vector field
subroutine odes_1ststepcon(self,kerr,pdfaca,firstcall,psi,rispldr,rispldz)

  !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(inout) :: kerr   !< error flag
  real(kr8), dimension(*), intent(in) :: pdfaca   !< direction factor array
  integer(ki4), intent(inout) :: firstcall  !< flag up first step of new trajectory
  type(spl2d_t), intent(inout) :: psi !< spl2d object data structure
  type(spl2d_t), intent(inout) :: rispldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: rispldz !< derivative spl2d object data structure

  !! local
  character(*), parameter :: s_name='odes_1ststepcon' !< subroutine name
  integer(ki4), parameter :: ifixdt=0 !< timestep control (1=fix,0=SW,2=Soderlind)
  real(kr8), dimension(mpdim) :: g !< local variable
  real(kr8), dimension(mpdim) :: g2 !< local variable
  real(kr8), dimension(mpdim) :: g3 !< local variable
  real(kr8), dimension(mpdim,3) :: gdia !< local variable
  real(kr8), dimension(mpdim) :: g3d !< local variable
  real(kr8), dimension(mpdim) :: g23d !< local variable
  real(kr8), dimension(mpdim) :: g3do !< local variable
  real(kr8) :: rs !< local variable
  real(kr8) :: zhf !< local variable
  real(kr8) :: zdt !< local variable
  real(kr8) :: zdir !< local variable
  real(kr8) :: zgincr !< local variable
  real(kr8) :: rso !< local variable
  real(kr8) :: t !< local variable
  integer(ki4) :: ik  !< winding number

  g=self%vecp%pos(self%ndt)%posvec(1:2)
  !     t=self%vecp%pos(self%ndt)%posvec(3)
  t=self%t
  ik=self%posk(self%ndt)
  ! save direction of advance
  zdir=sign(1._kr8,self%dt)
  !DIA  write(*,*) 'Entering 1ststepcon' !DIA
  !DIA  write(*,*) 't,self%dt,g,g2,g3,gdia,pdfaca',t,self%dt,g,g2,g3,gdia,pdfaca(1) !DIA
201   continue
  call odes_rkf23a(odes_1stfunct,t,self%dt,g,g2,g3,gdia,pdfaca,rispldr,rispldz)
  if (ifixdt==1) then
     self%dt=zdir*self%n%dt0
     self%t=self%t+self%dt
     self%ndt=self%ndt+1
     self%posk(self%ndt)=ik
     !     save position
     self%vecp%pos(self%ndt)%posvec(1:2)=g3
     self%vecp%pos(self%ndt)%posvec(3)=t
     kerr=0
     return
  end if
  g3d=abs(g3-g2)
  ! Soderlind control requires value from older timestep
  if (firstcall==1) then
     firstcall=0
     g3do=g3d
  else
     g3do=self%g3do(1:2)
  end if
  self%g3do(1:2)=g3d
  g23d=abs(g)+abs(g3)+self%epsar
  !DES rs=self%reps*maxval(g3d/g23d)
  rs=self%reps*maxval(g3d)
  !DOD        write(*,*) 'g3d=',g3d,'g23d=',g23d !DOD
  !DOD        write(*,*) 'rs=',rs !DOD

  if (rs<=1) then
     !     step succeeded, but if possible increase t-step for next step
     zhf=self%facu
     if (rs>self%facut) then
        if (ifixdt==0) then
           zhf=self%ffac/rs**self%rsord
        else
           rso=self%reps*maxval(g3do)
           zhf=self%ffac/(rs*rso/4)**self%rsordb
        end if
     end if
     zgincr=maxval(abs(g3-g))
     if (kerr==-1.OR.zgincr>self%n%termcon(2)) zhf=1
     !DOD              write(*,*) 'zhf=',zhf !DOD
     ! set new timestep, allow negative timestep
     zdt=zdir*max(zhf*abs(self%dt),self%dtmin)
     !        t=t+zdt
     self%t=self%t+self%dt
     self%ndt=self%ndt+1
     self%posk(self%ndt)=ik
     !     save position
     self%vecp%pos(self%ndt)%posvec(1:2)=g3
     self%vecp%pos(self%ndt)%posvec(3)=self%t
     !DOD              write(*,*) 'pass,t,dt,g3=',t,self%dt,g3 !DOD
     !D            write(20,*)i,t,self%dt,g,g2,g3 !D

     !        self%posa(1)=min(g2,g3,gdia(1),gdia(2),gdia(3))
     !        self%posb(1)=max(g2,g3,gdia(1),gdia(2),gdia(3))
     self%dt=zdt
     kerr=0
     !        g=g3
     !        if (zg3>self%glt) g=g3-2.*self%glt
     return
  end if
  !
  !     step failed, reduce t-step , equivalent to setting lfail=t
  kerr=-1

  zhf=self%facd
  if (rs<self%facdt) then
     if (ifixdt==0) then
        zhf=self%ffac/rs**self%rsord
     else
        rso=self%reps*maxval(g3do)
        zhf=self%ffac/(rs*rso/4)**self%rsordb
     end if
  end if
  zdt=zhf*self%dt
  self%dt=zdt
  !DOD        write(*,*) 'fail,t,dt,g3=',t,self%dt,g3 !DOD
  !DOD        write(*,*) 'failextra g,g2=',g,g2 !DOD
  !     try again unless t-step too small
  if (abs(self%dt)>=self%dtmin) goto 201

  call log_error(m_name,s_name,1,error_warning,'t-step too small')
  kerr=1

end subroutine odes_1ststepcon
!---------------------------------------------------------------------
!> special timestepping for 2nd scheme
subroutine odes_2ndstep(self,kerr,pdfaca,pf,psi,dspldr,dspldz)

  !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(inout) :: kerr   !< error flag
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  type(spl2d_t), intent(inout) :: psi !< spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data structure
  real(kr8), external :: pf !< local variable

  !! local
  character(*), parameter :: s_name='odes_2ndstep' !< subroutine name
  real(kr8), dimension(mpdim) :: g !< dependent variable at start time
  real(kr8), dimension(mpdim) :: g2 !< 2nd order accurate dependent variable at new time
  real(kr8), dimension(mpdim) :: g3 !< 3rd order accurate dependent variable at new time
  real(kr8), dimension(mpdim) :: g3dn !< local variable
  real(kr8), dimension(mpdim) :: g23n !< local variable
  real(kr8), dimension(mpdim,3) :: gdia !< diagnostics
  real(kr8) :: rs !< timestep control variable
  real(kr8) :: zhf !< timestep control variable
  real(kr8) :: rsta !<  \f$ R \f$  start value
  real(kr8) :: zsta !<  \f$ Z \f$  start value
  real(kr8) :: tsta !<  \f$ \xi \f$  start value
  real(kr8) :: dtsta !<  \f$ d\xi \f$  start value
  real(kr8) :: zval !<  required function value
  real(kr8) :: zdtn !< timestep for next timstep
  real(kr8) :: zdt !< timestep
  real(kr8) :: t !< \f$ \xi \f$ coordinate value
  real(kr8) :: zdth !< \f$ d\xi/2 \f$ for timestep
  real(kr8) :: zr   !<  \f$ R \f$ value
  real(kr8) :: zz   !<  \f$ Z \f$ value
  real(kr8) :: zt   !<  \f$ \xi \f$ value
  real(kr8) :: zpsi   !<  \f$ \psi \f$ value
  integer(ki4) :: ihalf  !< number of steps to half timestep
  integer(ki4) :: immax  !< number of steps to timestep
  integer(ki4) :: ik  !< winding number
  real(kr8) :: zfac   !<  \f$ d\xi \f$ factor
  real(kr8) :: zdir   !<  sign of \f$ d\xi \f$ advance
  real(kr8) :: thalf   !<  time at half timestep
  real(kr8) :: zrdiff!<  changes in \f$ R \f$ at half timestep
  real(kr8) :: zrdiffo   !<  changes in \f$ R \f$ at half timestep
  integer(ki4), parameter :: iparam=1 !< number of parameters
  real(kr8), dimension(iparam) :: zparam !< parameters

  ! starting values
  rsta=self%vecp%pos(self%ndt)%posvec(1)
  zsta=self%vecp%pos(self%ndt)%posvec(2)
  tsta=self%vecp%pos(self%ndt)%posvec(3)
  ik=self%posk(self%ndt)
  dtsta=self%dt
  zval=pf(rsta,zsta,tsta,0)
  zdt=dtsta
  ! save direction of advance
  zdir=sign(1._kr8,dtsta)

  ! solve step

201   continue

  ! calculate position at half-timestep (as function solution)
  zdth=zdt/2
  zfac=1
  immax=1
  do i=1, ipmaxhalf
     ihalf=immax
     immax=2*immax
     zr=rsta
     zt=tsta
     do j=1, ihalf
        zt=zt+zfac*zdth
        zparam(1)=zt
        call nrsolve_funct(pf,zval,zr,iparam,zparam,kerr)
        if (kerr>0) exit
     end do
     if (kerr==0) exit
     zfac=zfac/2
  end do
  if (kerr>0) return ! catastrophic error, give up this integration
  thalf=zt
  zrdiff=abs(zr-rsta)

  call spl2d_eval(psi,zr,zsta,zpsi)
  ! Runge-Kutta step, full timestep
  g(1)=zr
  g(2)=zsta
  t=tsta
  call odes_rkf23m(odes_2ndfunct,t,zdt,g,g2,g3,gdia,pdfaca,pf,dspldr,dspldz)
  g3dn=abs(g3-g2)
  g23n=abs(g)+abs(g3)+self%epsar
  rs=self%reps*maxval(g3dn/g23n)

  if (rs<=1) then
     !     step succeeded, but if possible increase t-step for next step
     zrdiffo=abs(g3(1)-zr)
     ! may not help to ensure psi conserved, commented out
     zr=g3(1)
     zz=g3(2)
     if (ipsiconst/=0) then
        call nrsolve_spl2dfn(psi,dspldr,dspldz,zpsi,zr,zz,kerr)
        !E          zrzdiff=sqrt((g3(1)-zr)**2+(g3(2)-zz)**2) !E
        !E          call spl2d_eval(psi,zr,zz,zpsicheck) !E
        !E          write(*,*) zrdiff,zrzdiff,abs(zpsi-zpsicheck) !E
        if (kerr>0) kerr=0 ! ignore failure to converge
     end if
     zhf=self%facu
     if (rs>self%facut) zhf=self%ffac/rs**self%rsord
     if (kerr==-1) zhf=1
     ! do not increase timestep if ripple is causing greater travel
     if (zrdiff>zrdiffo) zhf=1
     ! set new timestep, change to allow negative timestep
     zdtn=zdir*max(zhf*abs(zdt),self%dtmin)
     ! update t
     t=t+zdt
     ! cancel error flag (-1)
     kerr=0
     !D       write(*,*) i,t,zdt,abs(g)  !D
     !D       write(20,*)i,t,zdt,g,g2,g3  !D

  else
     !     step failed, reduce t-step , equivalent to setting lfail=t
     kerr=-1
     zhf=self%facd
     if (rs<self%facdt) zhf=self%ffac/rs**self%rsord
     zdtn=zhf*zdt
     !D        write(*,*) 'failure',i,self%dt,g,g2,g3 !D
     !     try again unless t-step too small
     if (abs(self%dt)>=self%dtmin) then
        zdt=zdtn
        goto 201
     else

        call log_error(m_name,s_name,1,error_warning,'t-step too small')
        kerr=1
        return
     end if
  end if

  rsta=zr
  zsta=zz
  tsta=thalf
  zval=pf(rsta,zsta,tsta,0)
  ! calculate position at second half-time step (as function solution)
  zfac=1
  immax=1
  do i=1, ipmaxhalf
     ihalf=immax
     immax=2*immax
     zr=rsta
     zt=tsta
     do j=1, ihalf
        zt=zt+zfac*zdth
        zparam(1)=zt
        call nrsolve_funct(pf,zval,zr,iparam,zparam,kerr)
        if (kerr>0) exit
     end do
     if (kerr==0) exit
     zfac=zfac/2
  end do
  if (kerr>0) return ! catastrophic error, give up this integration

  self%ndt=self%ndt+1
  self%vecp%pos(self%ndt)%posvec(1)=zr
  self%vecp%pos(self%ndt)%posvec(2)=zsta
  self%vecp%pos(self%ndt)%posvec(3)=zt
  self%posk(self%ndt)=ik
  ! timstep for next time
  self%dt=zdtn

end subroutine odes_2ndstep
!---------------------------------------------------------------------
!> special timestepping for 3-D vector field
subroutine odes_3rdstep(self,kerr,pdfaca,spl3d,psi,dspldr,dspldz)

  !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(inout) :: kerr   !< error flag
  real(kr8), dimension(mpvdim), intent(in) :: pdfaca   !< direction factor array
  type(spl3d_t), intent(inout) :: spl3d !< spl2d object data structure
  type(spl2d_t), intent(inout) :: psi !< spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data str`

  !! local
  character(*), parameter :: s_name='odes_3rdstep' !< subroutine name
  integer(ki4), parameter :: ifixdt=0 !< fix timestep
  real(kr8), dimension(mpvdim) :: g !< local variable
  real(kr8), dimension(mpvdim) :: g2 !< local variable
  real(kr8), dimension(mpvdim) :: g3 !< local variable
  real(kr8), dimension(mpvdim,3) :: gdia !< local variable
  real(kr8), dimension(mpvdim) :: g3d !< local variable
  real(kr8), dimension(mpvdim) :: g23d !< local variable
  real(kr8) :: rs !< local variable
  real(kr8) :: zhf !< local variable
  real(kr8) :: zdt !< local variable
  real(kr8) :: zdir !< local variable
  real(kr8) :: zgincr !< local variable
  real(kr8) :: t !< local variable
  integer(ki4) :: ik  !< winding number

  g=self%vecp%pos(self%ndt)%posvec
  t=self%t
  ik=self%posk(self%ndt)
  ! save direction of advance
  zdir=sign(1._kr8,self%dt)
201   continue
  call odes_rkf23d(odes_3rdfunct,t,self%dt,g,g2,g3,gdia,pdfaca,spl3d,dspldr,dspldz)
  if (ifixdt==1) then
     self%dt=zdir*self%n%dt0
     self%t=self%t+self%dt
     self%ndt=self%ndt+1
     self%posk(self%ndt)=ik
     !     save position
     self%vecp%pos(self%ndt)%posvec=g3
     kerr=0
     return
  end if
  g3d=abs(g3-g2)
  g23d=abs(g)+abs(g3)+self%epsar
  !DES     rs=self%reps*maxval(g3d/g23d)
  rs=self%reps*maxval(g3d)
  !DOD      write(*,*) 'g3d=',g3d,'g23d=',g23d !DOD
  !DOD      write(*,*) 'rs=',rs !DOD

  if (rs<=1) then
     !     step succeeded, but if possible increase t-step for next step
     zhf=self%facu
     if (rs>self%facut) zhf=self%ffac/rs**self%rsord
     zgincr=maxval(abs(g3-g))
     if (kerr==-1.OR.zgincr>self%n%termcon(2)) zhf=1
     !DOD         write(*,*) 'zhf=',zhf !DOD
     ! set new timestep, allow negative timestep
     zdt=zdir*max(zhf*abs(self%dt),self%dtmin)
     !        t=t+zdt
     self%t=self%t+self%dt
     self%ndt=self%ndt+1
     self%posk(self%ndt)=ik
     !     save position
     self%vecp%pos(self%ndt)%posvec=g3
     !DOD         write(*,*) 'pass,t,dt,g3=',t,self%dt,g3 !DOD
     !D       write(20,*)i,t,self%dt,g,g2,g3 !D

     !        self%posa(1)=min(g2,g3,gdia(1),gdia(2),gdia(3))
     !        self%posb(1)=max(g2,g3,gdia(1),gdia(2),gdia(3))
     self%dt=zdt
     kerr=0
     !        g=g3
     !        if (zg3>self%glt) g=g3-2.*self%glt
     return
  end if
  !
  !     step failed, reduce t-step , equivalent to setting lfail=t
  kerr=-1

  zhf=self%facd
  if (rs<self%facdt) zhf=self%ffac/rs**self%rsord
  zdt=zhf*self%dt
  self%dt=zdt
  !DOD      write(*,*) 'fail,t,dt,g3=',t,self%dt,g3 !DOD
  !DOD      write(*,*) 'failextra g,g2=',g,g2 !DOD
  !     try again unless t-step too small
  if (abs(self%dt)>=self%dtmin) goto 201

  call log_error(m_name,s_name,1,error_warning,'t-step too small')
  kerr=1

end subroutine odes_3rdstep
!---------------------------------------------------------------------
!> special timestepping for axisymm 3-D vector field
subroutine odes_4thstep(self,kerr,pdfaca,pf,psi,dspldr,dspldz)

  !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(inout) :: kerr   !< error flag
  real(kr8), dimension(*), intent(in) :: pdfaca   !< direction factor array
  real(kr8), intent(inout) :: pf   !< \f$ f \f$
  type(spl2d_t), intent(inout) :: psi !< spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data str`

  !! local
  character(*), parameter :: s_name='odes_4thstep' !< subroutine name
  integer(ki4), parameter :: ifixdt=0 !< fix timestep
  real(kr8), dimension(mpvdim) :: g !< local variable
  real(kr8), dimension(mpvdim) :: g2 !< local variable
  real(kr8), dimension(mpvdim) :: g3 !< local variable
  real(kr8), dimension(mpvdim,3) :: gdia !< local variable
  real(kr8), dimension(mpvdim) :: g3d !< local variable
  real(kr8), dimension(mpvdim) :: g23d !< local variable
  real(kr8) :: rs !< local variable
  real(kr8) :: zhf !< local variable
  real(kr8) :: zdt !< local variable
  real(kr8) :: zdir !< local variable
  real(kr8) :: zgincr !< local variable
  real(kr8) :: t !< local variable
  integer(ki4) :: ik  !< winding number

  g=self%vecp%pos(self%ndt)%posvec
  t=self%t
  ik=self%posk(self%ndt)
  ! save direction of advance
  zdir=sign(1._kr8,self%dt)
201   continue
  call odes_rkf23p(odes_4thfunct,t,self%dt,g,g2,g3,gdia,pdfaca,pf,dspldr,dspldz)
  if (ifixdt==1) then
     self%dt=zdir*self%n%dt0
     self%t=self%t+self%dt
     self%ndt=self%ndt+1
     self%posk(self%ndt)=ik
     !     save position
     self%vecp%pos(self%ndt)%posvec=g3
     kerr=0
     return
  end if
  g3d=abs(g3-g2)
  g23d=abs(g)+abs(g3)+self%epsar
  !DES     rs=self%reps*maxval(g3d/g23d)
  rs=self%reps*maxval(g3d)
  !DOD      write(*,*) 'g3d=',g3d,'g23d=',g23d !DOD
  !DOD      write(*,*) 'rs=',rs !DOD

  if (rs<=1) then
     !     step succeeded, but if possible increase t-step for next step
     zhf=self%facu
     if (rs>self%facut) zhf=self%ffac/rs**self%rsord
     zgincr=maxval(abs(g3-g))
     if (kerr==-1.OR.zgincr>self%n%termcon(2)) zhf=1
     !DOD         write(*,*) 'zhf=',zhf !DOD
     ! set new timestep, allow negative timestep
     zdt=zdir*max(zhf*abs(self%dt),self%dtmin)
     !        t=t+zdt
     self%t=self%t+self%dt
     self%ndt=self%ndt+1
     self%posk(self%ndt)=ik
     !     save position
     self%vecp%pos(self%ndt)%posvec=g3
     !DOD         write(*,*) 'pass,t,dt,g3=',t,self%dt,g3 !DOD
     !D       write(20,*)i,t,self%dt,g,g2,g3 !D

     !        self%posa(1)=min(g2,g3,gdia(1),gdia(2),gdia(3))
     !        self%posb(1)=max(g2,g3,gdia(1),gdia(2),gdia(3))
     self%dt=zdt
     kerr=0
     !        g=g3
     !        if (zg3>self%glt) g=g3-2.*self%glt
     return
  end if
  !
  !     step failed, reduce t-step , equivalent to setting lfail=t
  kerr=-1

  zhf=self%facd
  if (rs<self%facdt) zhf=self%ffac/rs**self%rsord
  zdt=zhf*self%dt
  self%dt=zdt
  !DOD      write(*,*) 'fail,t,dt,g3=',t,self%dt,g3 !DOD
  !DOD      write(*,*) 'failextra g,g2=',g,g2 !DOD
  !     try again unless t-step too small
  if (abs(self%dt)>=self%dtmin) goto 201

  call log_error(m_name,s_name,1,error_warning,'t-step too small')
  kerr=1

end subroutine odes_4thstep
!---------------------------------------------------------------------
!> extra-special timestepping for axisymm 3-D vector field
subroutine odes_4thstepcon(self,kerr,pdfaca,firstcall,pf,psi,dspldr,dspldz)

  !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(inout) :: kerr   !< error flag
  real(kr8), dimension(*), intent(in) :: pdfaca   !< direction factor array
  real(kr8), intent(inout) :: pf   !< \f$ f \f$
  type(spl2d_t), intent(inout) :: psi !< spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data str`
  integer(ki4), intent(inout) :: firstcall  !< flag up first step of new trajectory

  !! local
  character(*), parameter :: s_name='odes_4thstepcon' !< subroutine name
  integer(ki4), parameter :: ifixdt=0 !< timestep control (1=fix,0=SW,2=Soderlind)
  integer(ki4), parameter :: izadu=1 !< timestep error by Zadunaisky
  real(kr8), dimension(mpvdim) :: g !< local variable
  real(kr8), dimension(mpvdim) :: g2 !< local variable
  real(kr8), dimension(mpvdim) :: g3 !< local variable
  real(kr8), dimension(mpvdim,3) :: gdia !< local variable
  real(kr8), dimension(mpvdim) :: g3d !< local variable
  real(kr8), dimension(mpvdim) :: g23d !< local variable
  real(kr8), dimension(mpvdim) :: g3do !< local variable
  real(kr8) :: rs !< local variable
  real(kr8) :: zhf !< local variable
  real(kr8) :: zdt !< local variable
  real(kr8) :: zdir !< local variable
  real(kr8) :: zgincr !< local variable
  real(kr8) :: rso !< local variable
  real(kr8) :: t !< local variable
  integer(ki4) :: ik  !< winding number

  g=self%vecp%pos(self%ndt)%posvec
  t=self%t
  ik=self%posk(self%ndt)
  ! save direction of advance
  zdir=sign(1._kr8,self%dt)
201   continue
  call odes_rkf23p(odes_4thfunct,t,self%dt,g,g2,g3,gdia,pdfaca,pf,dspldr,dspldz)
  if (ifixdt==1) then
     self%dt=zdir*self%n%dt0
     self%t=self%t+self%dt
     self%ndt=self%ndt+1
     self%posk(self%ndt)=ik
     !     save position
     self%vecp%pos(self%ndt)%posvec=g3
     kerr=0
     return
  end if
  if (izadu==1) then
     g3d=(g3-g)/self%dt
  else
     g3d=abs(g3-g2)
  end if
  ! Soderlind control (and Zadunaisky) requires value from older timestep
  if (firstcall==1) then
     firstcall=0
     g3do=g3d
  else
     g3do=self%g3do
     end if
     self%g3do=g3d
     g23d=abs(g)+abs(g3)+self%epsar
     !DES     rs=self%reps*maxval(g3d/g23d)
     if (izadu==1) then
        rs=self%reps*0.5_kr8*maxval(abs(g3d-g3do))*abs(self%dt)
     else
        rs=self%reps*maxval(g3d)
     end if
     !DOD      write(*,*) 'g3d=',g3d,'g23d=',g23d !DOD
     !DOD      write(*,*) 'g3do=',g3do !DOD
     !DOD      write(*,*) 'rs=',rs !DOD

     if (rs<=1) then
        !     step succeeded, but if possible increase t-step for next step
        zhf=self%facu
        if (rs>self%facut) then
           if (ifixdt==0) then
              zhf=self%ffac/rs**self%rsord
           else
              rso=self%reps*maxval(g3do)
              zhf=self%ffac/(rs*rso/4)**self%rsordb
           end if
        end if
        zgincr=maxval(abs(g3-g))
        if (kerr==-1.OR.zgincr>self%n%termcon(2)) zhf=1
        !DOD         write(*,*) 'zhf=',zhf !DOD
        ! set new timestep, allow negative timestep
        zdt=zdir*max(zhf*abs(self%dt),self%dtmin)
        !        t=t+zdt
        self%t=self%t+self%dt
        self%ndt=self%ndt+1
        self%posk(self%ndt)=ik
        !     save position
        self%vecp%pos(self%ndt)%posvec=g3
        !DOD         write(*,*) 'pass,t,dt,g3=',t,self%dt,g3 !DOD
        !D       write(20,*)i,t,self%dt,g,g2,g3 !D

        !        self%posa(1)=min(g2,g3,gdia(1),gdia(2),gdia(3))
        !        self%posb(1)=max(g2,g3,gdia(1),gdia(2),gdia(3))
        self%dt=zdt
        kerr=0
        !        g=g3
        !        if (zg3>self%glt) g=g3-2.*self%glt
        return
     end if
     !
     !     step failed, reduce t-step , equivalent to setting lfail=t
     kerr=-1

     zhf=self%facd
     if (rs<self%facdt) then
        if (ifixdt==0) then
           zhf=self%ffac/rs**self%rsord
        else
           rso=self%reps*maxval(g3do)
           zhf=self%ffac/(rs*rso/4)**self%rsordb
        end if
     end if
     zdt=zhf*self%dt
     self%dt=zdt
     !DOD      write(*,*) 'fail,t,dt,g3=',t,self%dt,g3 !DOD
     !DOD      write(*,*) 'failextra g,g2=',g,g2 !DOD
     !     try again unless t-step too small
     if (abs(self%dt)>=self%dtmin) goto 201

     call log_error(m_name,s_name,1,error_warning,'t-step too small')
     kerr=1

end subroutine odes_4thstepcon
!---------------------------------------------------------------------
!> special first timestep for \f$ R/J \f$
subroutine odes_rjstep1(self,kerr,pf,rjspl2d)

     !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(out) :: kerr   !< error flag
  real(kr8), intent(in) :: pf   !< \f$ f \f$
  type(spl2d_t), intent(inout) :: rjspl2d   !< \f$ R/J \f$

     !! local
  character(*), parameter :: s_name='odes_rjstep1' !< subroutine name
  integer(ki4), parameter :: ipdtfac=1 !< multiplies dt on first step
  real(kr8) :: g !< function value
  real(kr8) :: gdot !< function derivative value
  real(kr8) :: gn !< function norm
  real(kr8) :: gdotn !< function derivative norm
  real(kr8) :: zepst  !< local variable
  real(kr8) :: zpsi  !< local variable
  real(kr8) :: t  !< local variable
  real(kr8), dimension(mpdim,3) :: gdia !< diagnostic

     kerr=0
     zpsi=self%vecp%pos(self%ndt)%posvec(1)
     g=self%vecp%pos(self%ndt)%posvec(2)
     t=self%vecp%pos(self%ndt)%posvec(3)
     call odes_rjfunct(t,g,gdot,zpsi,gdia,pf,rjspl2d)
     gn=abs(g)
     gdotn=abs(gdot)
     zepst =self%n%epsr*gn+self%n%epsa
     if (zepst .lt. gdotn**self%nsord) then
        self%dt=(zepst/gdotn)**self%rsord
     else
        self%dt=self%n%dt0
     end if
     self%dt=ipdtfac*max(self%dt,self%dtmin)
     ! estimate first step of trajectory with positive dt
     self%vecp%pos(self%ndt+1)%posvec(1)=zpsi
     self%vecp%pos(self%ndt+1)%posvec(2)=g+gdot*self%dt
     self%vecp%pos(self%ndt+1)%posvec(3)=t+self%dt

end subroutine odes_rjstep1
!---------------------------------------------------------------------
!> special first timestep for phi-stepping axisymm 3-D vector field
subroutine odes_1ststep1(self,kerr,pdfaca,psi,rispldr,rispldz,kflag)

     !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(out) :: kerr   !< error flag
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  type(spl2d_t), intent(inout) :: psi !< spl2d object data st ructure
  type(spl2d_t), intent(inout) :: rispldr !< derivative spl2d object data st ructure
  type(spl2d_t), intent(inout) :: rispldz !< derivative spl2d object data st ructure
  integer(ki4), intent(in) :: kflag   !< zero for full half-step, 1 for half half-step

     !! local
  character(*), parameter :: s_name='odes_1ststep1' !< subroutine name
  integer(ki4), parameter :: ipdtfac=1 !< multiplies dt on first step
  real(kr8), dimension(mpdim) :: g !< function value
  real(kr8), dimension(mpdim) :: gdot !< function derivative value
  real(kr8) :: gn !< function norm
  real(kr8) :: gdotn !< function derivative norm
  real(kr8), dimension(mpdim,3) :: gdia !< diagnostics
  real(kr8) :: zepst  !< local variable
  real(kr8) :: t !< time
  real(kr8) :: zr   !<  \f$ R \f$ value
  real(kr8) :: zz   !<  \f$ Z \f$ value
     !     real(kr8) :: zt   !<  \f$ \xi \f$ value
  real(kr8) :: zdt !< timestep

     kerr=0
     g(1)=self%vecp%pos(self%ndt)%posvec(1)
     g(2)=self%vecp%pos(self%ndt)%posvec(2)
     t=self%vecp%pos(self%ndt)%posvec(3)
     call odes_1stfunct(t,g,gdot,gdia,pdfaca,rispldr,rispldz)
     gn=maxval(abs(g))
     gdotn=maxval(abs(gdot))
     zepst=self%n%epsr*gn+self%n%epsa
     if (zepst < gdotn**self%nsord) then
        self%dt=(zepst/gdotn)**self%rsord
     else
        self%dt=self%n%dt0
     end if
     !DIA     write(*,*) 'g,gn,gdotn,self%n%dt0',g,gn,gdotn,self%n%dt0 !DIA
     !DIA     write(*,*) 'ipdtfac,self%dt,self%dtmin',ipdtfac,self%dt,self%dtmin !DIA
     self%dt=ipdtfac*max(self%dt,self%dtmin)
     zdt=self%dt
     ! first step of trajectory with positive dt
     zr=g(1)+gdot(1)*zdt
     zz=g(2)+gdot(2)*zdt
     t=t+zdt
     self%vecp%pos(self%ndt+1)%posvec(1)=zr
     self%vecp%pos(self%ndt+1)%posvec(2)=zz
     self%vecp%pos(self%ndt+1)%posvec(3)=t

     !  self%dt already contains next time-step

end subroutine odes_1ststep1
!---------------------------------------------------------------------
!> special first timestep for 2nd scheme
subroutine odes_2ndstep1(self,kerr,pdfaca,pf,psi,dspldr,dspldz,kflag)

     !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(out) :: kerr   !< error flag
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  type(spl2d_t), intent(inout) :: psi !< spl2d object data st ructure
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data st ructure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data st ructure
  integer(ki4), intent(in) :: kflag   !< zero for full half-step, 1 for half half-step
  real(kr8), external :: pf !< local variable

     !! local
  character(*), parameter :: s_name='odes_2ndstep1' !< subroutine name
  integer(ki4), parameter :: ipdtfac=1 !< multiplies dt on first step
  real(kr8), dimension(mpdim) :: g !< function value
  real(kr8), dimension(mpdim) :: gdot !< function derivative value
  real(kr8) :: gn !< function norm
  real(kr8) :: gdotn !< function derivative norm
  real(kr8), dimension(mpdim,3) :: gdia !< diagnostics
  real(kr8) :: zepst  !< local variable
  real(kr8) :: t !< time
  real(kr8) :: zr   !<  \f$ R \f$ value
  real(kr8) :: zz   !<  \f$ Z \f$ value
  real(kr8) :: zt   !<  \f$ \xi \f$ value
  real(kr8) :: zpsi   !<  \f$ \psi \f$ value
  real(kr8) :: zdth !< half timestep
  real(kr8) :: rsta !<  \f$ R \f$  start value
  real(kr8) :: zsta !<  \f$ Z \f$  start value
  real(kr8) :: tsta !<  \f$ \xi \f$  start value
  real(kr8) :: dtsta !<  \f$ d\xi \f$  start value
  real(kr8) :: zval !<  required function value
  real(kr8) :: zdt1h !< first half timestep
  integer(ki4) :: ihalf  !< number of steps to time + half timestep
  integer(ki4) :: immax  !< number of steps to time + timestep
  real(kr8) :: zfac   !<  \f$ d\xi \f$ factor
  integer(ki4), parameter :: iparam=1 !< number of parameters
  real(kr8), dimension(iparam) :: zparam !< parameters

     ! first half-step by ODE solution
     kerr=0
     g(1)=self%vecp%pos(self%ndt)%posvec(1)
     g(2)=self%vecp%pos(self%ndt)%posvec(2)
     tsta=self%vecp%pos(self%ndt)%posvec(3)
     t=tsta
     call spl2d_eval(psi,g(1),g(2),zpsi)
     call odes_2ndfunct(t,0.D0,g,gdot,gdia,pdfaca,pf,dspldr,dspldz)
     gn=maxval(abs(g))
     gdotn=maxval(abs(gdot))
     zepst=self%n%epsr*gn+self%n%epsa
     if (zepst < gdotn**self%nsord) then
        self%dt=(zepst/gdotn)**self%rsord
     else
        self%dt=self%n%dt0
     end if
     self%dt=ipdtfac*max(self%dt,self%dtmin)
     zdth=self%dt/2
     ! first half-step of trajectory with positive dt
     zr=g(1)+gdot(1)*zdth
     zz=g(2)+gdot(2)*zdth
     t=t+zdth
     if (ipsiconst/=0) then
        ! ensure psi conserved
        call nrsolve_spl2dfn(psi,dspldr,dspldz,zpsi,zr,zz,kerr)
        if (kerr>0) kerr=0 ! ignore failure to converge
     end if

     ! second half-step by function solution
     ! starting values (tsta unchanged)
     rsta=zr
     zsta=zz
     dtsta=self%dt
     zval=pf(rsta,zsta,tsta,0)
     zdt1h=zdth

     ! solve step
     zfac=1
     immax=1
     do i=1, ipmaxhalf
        ihalf=immax
        immax=2*immax
        zr=rsta
        zt=tsta
        do j=1, ihalf
           zt=zt+zfac*zdt1h
           zparam(1)=zt
           call nrsolve_funct(pf,zval,zr,iparam,zparam,kerr)
           if (kerr>0) exit
        end do
        if (kerr==0) exit
        zfac=zfac/2
     end do
     if (kerr>0) return ! catastrophic error, give up this integration

     self%vecp%pos(self%ndt+1)%posvec(1)=zr
     self%vecp%pos(self%ndt+1)%posvec(2)=zsta
     self%vecp%pos(self%ndt+1)%posvec(3)=zt
     !  self%dt already contains next time-step

end subroutine odes_2ndstep1
!---------------------------------------------------------------------
!> special first timestep for 3-D vector field
subroutine odes_3rdstep1(self,kerr,pdfaca,spl3d,psi,dspldr,dspldz,kflag)

     !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(out) :: kerr   !< error flag (dummy)
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  type(spl3d_t), intent(inout) :: spl3d !< spl2d object data structure
  type(spl2d_t), intent(inout) :: psi !< spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data str`
  integer(ki4), intent(in) :: kflag   !< dummy


     !! local
  character(*), parameter :: s_name='odes_3rdstep1' !< subroutine name
  integer(ki4), parameter :: ipdtfac=1 !< multiplies dt on first step
  real(kr8), dimension(mpvdim) :: g !< function value
  real(kr8), dimension(mpvdim) :: gdot !< function derivative value
  real(kr8), dimension(mpvdim,3) :: gdia !< diagnostics
  real(kr8) :: gn !< function norm
  real(kr8) :: gdotn !< function derivative norm
  real(kr8) :: zepst  !< local variable
  real(kr8) :: t !< pseudo-time

     kerr=0
     g=self%vecp%pos(self%ndt)%posvec
     t=self%t
     !DA      do j=1,80 !DA
     !DA      g(3)=g(3)+1. !DA
     call odes_3rdfunct(t,g,gdot,gdia,pdfaca,spl3d,dspldr,dspldz)
     !DA      end do !DA
     !DA      stop !DA
     gn=maxval(abs(g))
     gdotn=maxval(abs(gdot))
     !DES     zepst=self%n%epsr*gn+self%n%epsa
     zepst=self%n%epsr+self%n%epsa
     if (zepst .lt. gdotn**self%nsord) then
        self%dt=(zepst/gdotn)**self%rsord
     else
        self%dt=self%n%dt0
     end if
     self%dt=ipdtfac*max(self%dt,self%dtmin)
     ! estimate first step of trajectory with positive dt
     self%vecp%pos(self%ndt+1)%posvec=g+gdot*self%dt

end subroutine odes_3rdstep1
!---------------------------------------------------------------------
!> special first timestep for axisymm 3-D vector field
subroutine odes_4thstep1(self,kerr,pdfaca,pf,psi,dspldr,dspldz,kflag)

     !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(out) :: kerr   !< error flag (dummy)
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  real(kr8), intent(inout) :: pf   !< \f$ f \f$
  type(spl2d_t), intent(inout) :: psi !< spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data str`
  integer(ki4), intent(in) :: kflag   !< dummy


     !! local
  character(*), parameter :: s_name='odes_4thstep1' !< subroutine name
  integer(ki4), parameter :: ipdtfac=1 !< multiplies dt on first step
  real(kr8), dimension(mpvdim) :: g !< function value
  real(kr8), dimension(mpvdim) :: gdot !< function derivative value
  real(kr8), dimension(mpvdim,3) :: gdia !< diagnostics
  real(kr8) :: gn !< function norm
  real(kr8) :: gdotn !< function derivative norm
  real(kr8) :: zepst  !< local variable
  real(kr8) :: t !< pseudo-time

     kerr=0
     g=self%vecp%pos(self%ndt)%posvec
     t=self%t
     !DA      do j=1,80 !DA
     !DA      g(3)=g(3)+1. !DA
     call odes_4thfunct(t,g,gdot,gdia,pdfaca,pf,dspldr,dspldz)
     !DA      end do !DA
     !DA      stop !DA
     gn=maxval(abs(g))
     gdotn=maxval(abs(gdot))
     !DES     zepst=self%n%epsr*gn+self%n%epsa
     zepst=self%n%epsr+self%n%epsa
     if (zepst .lt. gdotn**self%nsord) then
        self%dt=(zepst/gdotn)**self%rsord
     else
        self%dt=self%n%dt0
     end if
     self%dt=ipdtfac*max(self%dt,self%dtmin)
     ! estimate first step of trajectory with positive dt
     self%vecp%pos(self%ndt+1)%posvec=g+gdot*self%dt

end subroutine odes_4thstep1
!---------------------------------------------------------------------
!> RKF 23 step
subroutine odes_rkf23n(pfunct,pt,pdt,py,py2,py3,pdia,&
     psi,pf,rjspl2d)
     !     rkf23 scheme for ITER work
  real(kr8), intent(in) :: pt !< starting time
  real(kr8), intent(in) :: pdt !< timestep
  real(kr8), intent(inout) :: py !< starting value
  real(kr8), intent(inout) :: py2 !< 2nd order estimate
  real(kr8), intent(inout) :: py3 !< 3rd order estimate
  real(kr8), intent(out), dimension(*) :: pdia !< diagnostic
  real(kr8), intent(in) :: psi !< \f$ \psi \f$
  real(kr8), intent(in) :: pf   !< \f$ f \f$
  type(spl2d_t), intent(inout) :: rjspl2d   !< \f$ R/J \f$ spline

     external pfunct

     !! local
  character(*), parameter :: s_name='odes_rkf23n' !< subroutine name
  real(kr8) :: zy1  !< 1st order estimate
  real(kr8) :: zy1b  !< 1st order estimate
  real(kr8) :: zy0  !< 0th order estimate
  real(kr8) :: zy0b  !< 0th order estimate
  real(kr8) :: zy1bb  !< intermediate 0th order estimate
  real(kr8) :: zy2b  !< 2nd order estimate
  real(kr8) :: zydot!< derivative estimates
  real(kr8) :: zyd1!< derivative estimates
  real(kr8) :: zyd1bb   !< derivative estimates

     call pfunct(pt,py,zydot,psi,pdia,pf,rjspl2d)
     zy1=py+pdt*zydot
     call pfunct(pt,zy1,zyd1,psi,pdia,pf,rjspl2d)
     zy1b=py+pdt*zyd1
     py2=(zy1+zy1b)/2

     zy0=py+(pdt/3)*zydot
     zy0b=py+(2*pdt/3)*zyd1
     zy1bb=(zy0+zy0b)/2
     call pfunct(pt,zy1bb,zyd1bb,psi,pdia,pf,rjspl2d)
     zy2b=py+pdt*zyd1bb
     py3=(py2+zy2b)/2
     pdia(1)=zy1
     pdia(2)=zy1b
     pdia(3)=zy2b

end subroutine odes_rkf23n
!---------------------------------------------------------------------
!> RKF 23 step for R and Z
subroutine odes_rkf23a(pfunct,pt,pdt,py,py2,py3,pdia,pdfaca,rispldr,rispldz)
     !     rkf23 scheme for 'msus','global' work
  real(kr8), intent(in) :: pt !< starting time
  real(kr8), intent(in) :: pdt !< timestep
  real(kr8),  intent(inout), dimension(mpdim) :: py !< starting value
  real(kr8),  intent(inout), dimension(mpdim) :: py2 !< 2nd order estimate
  real(kr8),  intent(inout), dimension(mpdim) :: py3 !< 3rd order estimate
  real(kr8), intent(out), dimension(mpdim,*) :: pdia !< diagnostic
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  type(spl2d_t), intent(inout) :: rispldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: rispldz !< derivative spl2d object data structure

     external pfunct

     !! local
  character(*), parameter :: s_name='odes_rkf23a' !< subroutine name
  real(kr8), dimension(mpdim) :: zy1 !< 1st order estimate
  real(kr8), dimension(mpdim) :: zy1b !< 1st order estimate
  real(kr8), dimension(mpdim) :: zy1bb !< intermediate 0th order estimate
  real(kr8), dimension(mpdim) :: zy2b !< 2nd order estimate
  real(kr8), dimension(mpdim) :: zydot !< derivative estimate
  real(kr8), dimension(mpdim) :: zyd1 !< derivative estimate
  real(kr8), dimension(mpdim) :: zyd1bb !< derivative estimate

     call pfunct(pt,py,zydot,pdia,pdfaca,rispldr,rispldz)
     zy1=py+pdt*zydot
     call pfunct(pt,zy1,zyd1,pdia,pdfaca,rispldr,rispldz)
     zy1b=py+pdt*zyd1
     py2=(zy1+zy1b)/2

     zy1bb=(py+py2)/2
     call pfunct(pt,zy1bb,zyd1bb,pdia,pdfaca,rispldr,rispldz)
     zy2b=py+pdt*zyd1bb
     py3=(py2+2*zy2b)/3
     pdia(:,1)=zy1
     pdia(:,2)=zy1b
     pdia(:,3)=zy2b

end subroutine odes_rkf23a
!---------------------------------------------------------------------
!> RKF 23 step for explicit t-dependence
subroutine odes_rkf23m(pfunct,pt,pdt,py,py2,py3,pdia,&
     pdfaca,pf,dspldr,dspldz)
     !     rkf23 scheme for 'msus','global' work
  real(kr8), intent(in) :: pt !< starting time
  real(kr8), intent(in) :: pdt !< timestep
  real(kr8),  intent(inout), dimension(mpdim) :: py !< starting value
  real(kr8),  intent(inout), dimension(mpdim) :: py2 !< 2nd order estimate
  real(kr8),  intent(inout), dimension(mpdim) :: py3 !< 3rd order estimate
  real(kr8), intent(out), dimension(mpdim,*) :: pdia !< diagnostic
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data structure

     external pfunct
  real(kr8), external :: pf !< local variable

     !! local
  character(*), parameter :: s_name='odes_rkf23m' !< subroutine name
  real(kr8), dimension(mpdim) :: zy1 !< 1st order estimate
  real(kr8), dimension(mpdim) :: zy1b !< 1st order estimate
  real(kr8), dimension(mpdim) :: zy1bb !< intermediate 0th order estimate
  real(kr8), dimension(mpdim) :: zy2b !< 2nd order estimate
  real(kr8), dimension(mpdim) :: zydot !< derivative estimate
  real(kr8), dimension(mpdim) :: zyd1 !< derivative estimate
  real(kr8), dimension(mpdim) :: zyd1bb !< derivative estimate

     call pfunct(pt,pdt,py,zydot,pdia,pdfaca,pf,dspldr,dspldz)
     !     call pfunct(pt+pdt/2,py,zydot,pdia,pdfaca,pf,dspldr,dspldz)
     zy1=py+pdt*zydot
     call pfunct(pt,pdt,zy1,zyd1,pdia,pdfaca,pf,dspldr,dspldz)
     !     call pfunct(pt+pdt/2,zy1,zyd1,pdia,pdfaca,pf,dspldr,dspldz)
     zy1b=py+pdt*zyd1
     py2=(zy1+zy1b)/2

     zy1bb=(py+py2)/2
     call pfunct(pt,pdt,zy1bb,zyd1bb,pdia,pdfaca,pf,dspldr,dspldz)
     zy2b=py+pdt*zyd1bb
     py3=(py2+2*zy2b)/3
     pdia(:,1)=zy1
     pdia(:,2)=zy1b
     pdia(:,3)=zy2b

end subroutine odes_rkf23m
!---------------------------------------------------------------------
!> RKF 23 step for 3-D vector, t-dependence
subroutine odes_rkf23d(pfunct,pt,pdt,py,py2,py3,pdia,&
     pdfaca,spl3d,dspldr,dspldz)
     !     rkf23 scheme for 'msus','global' work
  real(kr8), intent(in) :: pt !< starting time
  real(kr8), intent(in) :: pdt !< timestep
  real(kr8),  intent(inout), dimension(mpvdim) :: py !< starting value
  real(kr8),  intent(inout), dimension(mpvdim) :: py2 !< 2nd order estimate
  real(kr8),  intent(inout), dimension(mpvdim) :: py3 !< 3rd order estimate
  real(kr8), intent(out), dimension(mpvdim,*) :: pdia !< diagnostic
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  type(spl3d_t), intent(inout) :: spl3d !< spl3d object data structure
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data structure

     external pfunct

     !! local
  character(*), parameter :: s_name='odes_rkf23d' !< subroutine name
  real(kr8), dimension(mpvdim) :: zy1 !< 1st order estimate
  real(kr8), dimension(mpvdim) :: zy1b !< 1st order estimate
  real(kr8), dimension(mpvdim) :: zy1bb !< intermediate 0th order estimate
  real(kr8), dimension(mpvdim) :: zy2b !< 2nd order estimate
  real(kr8), dimension(mpvdim) :: zydot !< derivative estimate
  real(kr8), dimension(mpvdim) :: zyd1 !< derivative estimate
  real(kr8), dimension(mpvdim) :: zyd1bb !< derivative estimate

     call pfunct(pt,py,zydot,pdia,pdfaca,spl3d,dspldr,dspldz)
     !     call pfunct(pt+pdt/2,py,zydot,pdia,pdfaca,spl3d,dspldr,dspldz)
     zy1=py+pdt*zydot
     call pfunct(pt,zy1,zyd1,pdia,pdfaca,spl3d,dspldr,dspldz)
     !     call pfunct(pt+pdt/2,zy1,zyd1,pdia,pdfaca,spl3d,dspldr,dspldz)
     zy1b=py+pdt*zyd1
     py2=(zy1+zy1b)/2

     zy1bb=(py+py2)/2
     call pfunct(pt,zy1bb,zyd1bb,pdia,pdfaca,spl3d,dspldr,dspldz)
     zy2b=py+pdt*zyd1bb
     py3=(py2+2*zy2b)/3
     pdia(:,1)=zy1
     pdia(:,2)=zy1b
     pdia(:,3)=zy2b

end subroutine odes_rkf23d
!---------------------------------------------------------------------
!> RKF 23 step for 3-D vector, t-dependence
subroutine odes_rkf23p(pfunct,pt,pdt,py,py2,py3,pdia,&
     pdfaca,pf,dspldr,dspldz)
     !     rkf23 scheme for 'msus','global' work
  real(kr8), intent(in) :: pt !< starting time
  real(kr8), intent(in) :: pdt !< timestep
  real(kr8),  intent(inout), dimension(mpvdim) :: py !< starting value
  real(kr8),  intent(inout), dimension(mpvdim) :: py2 !< 2nd order estimate
  real(kr8),  intent(inout), dimension(mpvdim) :: py3 !< 3rd order estimate
  real(kr8), intent(out), dimension(mpvdim,*) :: pdia !< diagnostic
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  real(kr8), intent(inout) :: pf   !< \f$ f \f$
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data structure

     external pfunct

     !! local
  character(*), parameter :: s_name='odes_rkf23p' !< subroutine name
  real(kr8), dimension(mpvdim) :: zy1 !< 1st order estimate
  real(kr8), dimension(mpvdim) :: zy1b !< 1st order estimate
  real(kr8), dimension(mpvdim) :: zy1bb !< intermediate 0th order estimate
  real(kr8), dimension(mpvdim) :: zy2b !< 2nd order estimate
  real(kr8), dimension(mpvdim) :: zydot !< derivative estimate
  real(kr8), dimension(mpvdim) :: zyd1 !< derivative estimate
  real(kr8), dimension(mpvdim) :: zyd1bb !< derivative estimate

     call pfunct(pt,py,zydot,pdia,pdfaca,pf,dspldr,dspldz)
     !     call pfunct(pt+pdt/2,py,zydot,pdia,pdfaca,pf,dspldr,dspldz)
     zy1=py+pdt*zydot
     call pfunct(pt,zy1,zyd1,pdia,pdfaca,pf,dspldr,dspldz)
     !     call pfunct(pt+pdt/2,zy1,zyd1,pdia,pdfaca,pf,dspldr,dspldz)
     zy1b=py+pdt*zyd1
     py2=(zy1+zy1b)/2

     zy1bb=(py+py2)/2
     call pfunct(pt,zy1bb,zyd1bb,pdia,pdfaca,pf,dspldr,dspldz)
     zy2b=py+pdt*zyd1bb
     py3=(py2+2*zy2b)/3
     pdia(:,1)=zy1
     pdia(:,2)=zy1b
     pdia(:,3)=zy2b

end subroutine odes_rkf23p
!---------------------------------------------------------------------
!> special function \f$ R/J \f$
subroutine odes_rjfunct(pt,py,pydot,psi,pdia,pf,rjspl2d)

     !! arguments
  real(kr8), intent(inout) :: pt   !< \f$ \zeta \f$ (not needed)
  real(kr8), intent(inout) :: py   !< \f$ \theta \f$
  real(kr8), intent(inout) :: pydot  !< \f$ d\theta d\zeta \f$
  real(kr8), intent(inout) :: psi   !< \f$ \psi \f$ fixed for trajectory
  real(kr8), intent(out), dimension(mpdim,*) :: pdia !< diagnostic
  real(kr8), intent(in) :: pf   !< \f$ f \f$
  type(spl2d_t), intent(inout) :: rjspl2d   !< \f$ R/J \f$ spline

     !! local
  character(*), parameter :: s_name='odes_rjfunct' !< subroutine name
  real(kr8) :: zrj   !< \f$ R/J \f$ value

     call spl2d_eval(rjspl2d,psi,py,zrj)
     pydot=zrj/pf !< sign change kludge
     !D     write(*,*) pt,py,pydot,psi,pf !D

end subroutine odes_rjfunct
!---------------------------------------------------------------------
!> advance in 2-D autonomous
subroutine odes_1stfunct(pt,py,pydot,pdia,pdfaca,rispldr,rispldz)

     !! arguments
  real(kr8), intent(in) :: pt !< time
  real(kr8),  intent(inout), dimension(mpdim) :: py !< starting value
  real(kr8), intent(inout), dimension(mpdim) :: pydot !< \f$ d\theta d\zeta \f$
  real(kr8), intent(out), dimension(mpdim,*) :: pdia !< diagnostic
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  type(spl2d_t), intent(inout) :: rispldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: rispldz !< derivative spl2d object data structure

     !! local
  character(*), parameter :: s_name='odes_1stfunct' !< subroutine name
  real(kr8) :: zridpdr    !<  \f$ \frac{R}{I}\frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zridpdz    !<  \f$ \frac{R}{I}\frac{\partial\psi}{\partial Z} \f$

     call spl2d_evaln(rispldr,py(1),py(2),1,zridpdr)
     call spl2d_evaln(rispldz,py(1),py(2),2,zridpdz)
     pydot(1)=pdfaca(1)*zridpdz
     pydot(2)=pdfaca(2)*zridpdr
     !D     write(*,*) pt,py,pydot !D
     pdia(1,1)=zridpdr
     pdia(2,1)=zridpdz

end subroutine odes_1stfunct
!---------------------------------------------------------------------
!> advance in 2-D with time dependent RHS
subroutine odes_2ndfunct(pt,pdt,py,pydot,pdia,pdfaca,pf,dspldr,dspldz)

     !! arguments
  real(kr8), intent(in) :: pt !< time
  real(kr8), intent(in) :: pdt !< timestep
  real(kr8),  intent(inout), dimension(mpdim) :: py !< starting value
  real(kr8), intent(inout), dimension(mpdim) :: pydot !< \f$ d\theta d\zeta \f$
  real(kr8), intent(out), dimension(mpdim,*) :: pdia !< diagnostic
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data structure

  real(kr8), external :: pf !< local variable

     !! local
  character(*), parameter :: s_name='odes_2ndfunct' !< subroutine name
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zdhdr    !<  \f$ \frac{\partial H}{\partial R} \f$
  real(kr8),parameter :: pgauss=0.788675135D0 !<  Gauss point on [0,1] \f$ 1/2+\sqrt{3}/6 \f$
  real(kr8),parameter :: pgaussm=1-pgauss    !<  Gauss point on [0,1]

     call spl2d_evaln(dspldr,py(1),py(2),1,zdpdr)
     call spl2d_evaln(dspldz,py(1),py(2),2,zdpdz)
     zdhdr=(pf(py(1),py(2),pt+pgaussm*pdt,1)+ &
 &   pf(py(1),py(2),pt+pgauss*pdt,1))/2
     pydot(1)=-pdfaca(1)*zdpdz/zdhdr
     pydot(2)=pdfaca(2)*zdpdr/zdhdr
     !D     write(*,*) pt,py,pydot !D
     pdia(1,1)=zdpdr
     pdia(2,1)=zdpdz

end subroutine odes_2ndfunct
!---------------------------------------------------------------------
!> advance 3-D vector field
subroutine odes_3rdfunct(pt,py,pydot,pdia,pdfaca,spl3d,dspldr,dspldz)

     !! arguments
  real(kr8), intent(in) :: pt !< time
  real(kr8),  intent(inout), dimension(mpvdim) :: py !< starting value
  real(kr8), intent(inout), dimension(mpvdim) :: pydot !< \f$ d\theta d\zeta \f$
  real(kr8), intent(out), dimension(mpvdim,*) :: pdia !< diagnostic
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  type(spl3d_t), intent(inout) :: spl3d !< spl3d object data structure
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data structure

     !! local
  character(*), parameter :: s_name='odes_3rdfunct' !< subroutine name
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8), dimension(3)  :: zvecf    !<  Vacuum field

     call spl2d_evaln(dspldr,py(1),py(2),1,zdpdr)
     call spl2d_evaln(dspldz,py(1),py(2),2,zdpdz)
     call spl3d_evalm(spl3d,py,zvecf)
     pydot(1)=pdfaca(1)*(-zdpdz+zvecf(1))
     pydot(2)=pdfaca(2)*(zdpdr+zvecf(2))
     pydot(3)=pdfaca(3)*zvecf(3)
     !D     write(*,*) pt,py,pydot !D
     pdia(1,1)=-zdpdz+zvecf(1)
     pdia(2,1)=zdpdr+zvecf(2)
     pdia(3,1)=zvecf(3)
     !D     write(*,*) 'pt,py,pydot,zvecf=',pt,py,pydot,zvecf !D

end subroutine odes_3rdfunct
!---------------------------------------------------------------------
!> advance 3-D vector field
subroutine odes_4thfunct(pt,py,pydot,pdia,pdfaca,pf,dspldr,dspldz)

     !! arguments
  real(kr8), intent(in) :: pt !< time
  real(kr8),  intent(inout), dimension(mpvdim) :: py !< starting value
  real(kr8), intent(inout), dimension(mpvdim) :: pydot !< \f$ d\theta d\zeta \f$
  real(kr8), intent(out), dimension(mpvdim,*) :: pdia !< diagnostic
  real(kr8), intent(in), dimension(*) :: pdfaca !< direction factor
  real(kr8), intent(inout) :: pf   !< \f$ f \f$
  type(spl2d_t), intent(inout) :: dspldr !< derivative spl2d object data structure
  type(spl2d_t), intent(inout) :: dspldz !< derivative spl2d object data structure

     !! local
  character(*), parameter :: s_name='odes_4thfunct' !< subroutine name
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z} \f$

     call spl2d_evaln(dspldr,py(1),py(2),1,zdpdr)
     call spl2d_evaln(dspldz,py(1),py(2),2,zdpdz)
     pydot(1)=-pdfaca(1)*zdpdz
     pydot(2)=pdfaca(2)*zdpdr
     pydot(3)=pdfaca(3)*pf/(py(1)-pdfaca(4))
     !D     write(*,*) 'pt,py,pydot,zvecf=',pt,py,pydot,zvecf !D
     !D     write(*,*) pt,py,pydot !D
     pdia(1,1)=-zdpdz
     pdia(2,1)=zdpdr
     pdia(3,1)=pf

end subroutine odes_4thfunct
!---------------------------------------------------------------------
!> write (vtk) odes data (dummy)
subroutine odes_writev(self)

     !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure

     !! local
  character(*), parameter :: s_name='odes_writev' !< subroutine name


end subroutine odes_writev
!---------------------------------------------------------------------
!> output odes data (dummy)
subroutine odes_write(self)

     !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure

     !! local
  character(*), parameter :: s_name='odes_write' !< subroutine name


end subroutine odes_write
!---------------------------------------------------------------------
!> delete odes data (dummy)
subroutine odes_delete(self)

     !! arguments
  type(odes_t), intent(inout) :: self   !< object data structure

     !! local
  character(*), parameter :: s_name='odes_delete' !< subroutine name

     deallocate(self%g3do)
     deallocate(self%poso)
     deallocate(self%posa)
     deallocate(self%posb)
     deallocate(self%vecp%pos)
     deallocate(self%posk)

end subroutine odes_delete

end module odes_m
