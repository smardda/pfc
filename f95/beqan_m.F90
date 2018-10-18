module beqan_m

  use beqan_h
  use log_m
  use misc_m
  use const_kind_m

  implicit none
  private

! public subroutines
  public :: &
  beqan_init,  & !< open file
  beqan_readcon,  & !< read data from file
  beqan_solovevmove,  & !< calculate field, velocity and flux
  beqan_solovev,  & !< calculate analytic field after Solovev
  beqan_userdefined,  & !< userdefined analytic field
  beqan_generic,  & !< generic analytic call
  beqan_initwrite, & !< open new output file
  beqan_write, &  !< write out object
  beqan_delete, & !< delete object
  beqan_close, & !< close file
  beqan_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='beqan_m' !< module name
  integer  :: status   !< error status
  integer,save  :: ninba=5      !< control file unit number
  integer,save  :: noutba=6      !< output file unit number
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter

  contains
!---------------------------------------------------------------------
!> open file
subroutine beqan_init(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer, intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='beqan_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then channel=i exit end if end do ninba=i

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  call misc_getfileunit(ninba)
  open(unit=ninba,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open control data file')
     stop
  end if
  if (present(channel)) channel=ninba

end  subroutine beqan_init
!---------------------------------------------------------------------
!> read data from file
subroutine beqan_readcon(self,kin)

  !! arguments
  type(beqan_t), intent(out) :: self !< type which data will be assigned to
  integer, intent(in),optional :: kin   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='beqan_readcon' !< subroutine name
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  character(len=80) :: equilibrium_formula !< equilibrium equilibrium_formula
  real(kr8) :: BPhi !< Toroidal field
  real(kr8) :: Psi0 !< Flux through plasma
  real(kr8) :: R1 !< toroidal inner radius
  real(kr8) :: R2 !< toroidal outer radius
  real(kr8) :: Rm !< \f$ R \f$ of maximum height
  real(kr8) :: Zm !< \f$ Z \f$ of maximum height
  real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters !< general real parameters
  integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters !< general integer_parameters
  integer(ki4) :: number_of_real_parameters !< number of real parameters
  integer(ki4) :: number_of_integer_parameters !< number of integer parameters


  real(kr8) :: rmin !< \f$ R_{\min} \f$
  real(kr8) :: rmax !< \f$ R_{\max} \f$
  integer(ki4) :: nw !< number of grid points in R-direction
  real(kr8) :: zmin !< \f$ Z_{\min} \f$
  real(kr8) :: zmax !< \f$ Z_{\max} \f$
  integer(ki4) :: nh !< number of grid points in Z-direction

  !! equil parameters
  namelist /equilparameters/ equilibrium_formula, BPhi, Psi0, R1, R2, Rm, Zm, &
 &general_real_parameters, number_of_real_parameters, &
 &general_integer_parameters, number_of_integer_parameters

  namelist /meshparameters/ rmin, rmax, nw, zmin, zmax, nh

  !! set default equil parameters
  equilibrium_formula='solovev'
  !! set default field parameters
  BPhi=1.0
  Psi0=0.76225
  ! JET geometry
  R1=2.0
  R2=4.0
  Rm=2.828427125
  Zm=1.75
  ! doubled JET for approx ITER
  BPhi=2.0
  Psi0=1.5245
  ! doubled JET for approx ITER geometry
  R1=4.0
  R2=8.0
  Rm=5.65685425
  Zm=3.5

  !! userdefined formula
  general_real_parameters=0
  general_integer_parameters=0
  number_of_real_parameters=0
  number_of_integer_parameters=0

  if(present(kin)) then
     !! assume unit already open and reading infile
     ninba=kin
  end if

  !!read equil parameters
  read(ninba,nml=equilparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_warning,'Error reading equil parameters')
     call log_getunit(ilog)
     write(ilog,nml=equilparameters)
     print '("Error reading equil parameters")'
  end if

  if (R1<=0.OR.R2<=0) then
     call log_error(m_name,s_name,2,error_fatal,'Solovev R1, R2 should be positive')
  end if
  if (R1>=R2) then
     call log_error(m_name,s_name,3,error_fatal,'Solovev R1 should be less than R2')
  end if
  if (Rm>=R2.OR.Rm<=R1) then
     call log_error(m_name,s_name,4,error_fatal,'Solovev Rm should lie between R1 and R2')
  end if

  !! store equil values
  call lowor(equilibrium_formula,1,len_trim(equilibrium_formula))
  self%eqparms%formula=equilibrium_formula
  !! store field values
  self%eqparms%bphi=BPhi
  self%eqparms%psi0=Psi0
  !! store geometry values
  self%eqparms%r1=R1
  self%eqparms%r2=R2
  self%eqparms%rm=Rm
  self%eqparms%zm=Zm
  self%psiaxis=0
  self%psibdry=Psi0

  equil_formula: select case (equilibrium_formula)
  case ('solovev')
     self%a=(R2-R1)/2
     self%r0sq=(R1**2+R2**2)/2
     self%rxsq=((R1**2)*(R2**2)-Rm**4)/(R1**2+R2**2-2*(Rm**2))
     self%esq=(Zm**2)/(R1**2+R2**2-2*(Rm**2))
     self%rc=0.5*(R1+R2)

  case ('userdefined')
     if(number_of_real_parameters<0) &
 &   call log_error(m_name,s_name,14,error_fatal,'number of real parameters must be >=0')
     if(number_of_real_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of real_p rameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,15,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters<0) &
 &   call log_error(m_name,s_name,16,error_fatal,'number of integer parameters must be >=0')
     if(number_of_integer_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of integer parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,17,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters==0.AND.number_of_real_parameters==0) &
 &   call log_error(m_name,s_name,18,error_fatal,'no parameters set')

  case default
     call log_error(m_name,s_name,20,error_fatal,'Equilibrium with this formula not implemented')
  end select equil_formula

  self%eqparms%nrpams=number_of_real_parameters
  self%eqparms%nipams=number_of_integer_parameters

  !! allocate arrays and assign

  formula_allocate: select case (equilibrium_formula)
  case('userdefined')
     if (number_of_real_parameters>0) allocate(self%eqparms%rpar(number_of_real_parameters), stat=status)
     call log_alloc_check(m_name,s_name,65,status)
     if (number_of_integer_parameters>0) allocate(self%eqparms%npar(number_of_integer_parameters), stat=status)
     call log_alloc_check(m_name,s_name,66,status)
     self%eqparms%rpar=general_real_parameters(:number_of_real_parameters)
     self%eqparms%npar=general_integer_parameters(:number_of_integer_parameters)
  case default
  end select formula_allocate


  !! set default mesh parameters
  rmin=3
  rmax=9
  nw=100
  zmin=-5
  zmax=5
  nh=100
  !!read and store mesh parameters
  read(ninba,nml=meshparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,5,error_warning,'Error reading mesh parameters')
     call log_getunit(ilog)
     write(ilog,nml=meshparameters)
     print '("Error reading mesh parameters")'
  end if

  if (nw<=0) then
     call log_error(m_name,s_name,6,error_fatal,'need positive number of radial mesh points ')
  end if
  self%mr=nw-1
  if (nh<=0) then
     call log_error(m_name,s_name,6,error_fatal,'need positive number of vertical mesh points ')
  end if
  self%mz=nh-1

  if (rmin>rmax) then
     call log_error(m_name,s_name,12,error_warning,'rmin should be less than rmax')
     self%rmin=rmax
     self%rmax=rmin
  else if (rmin==rmax) then
     call log_error(m_name,s_name,13,error_fatal,'rmin should not be equal rmax')
  else
     self%rmin=rmin
     self%rmax=rmax
  end if

  if (zmin>zmax) then
     call log_error(m_name,s_name,15,error_warning,'zmin should be less than zmax')
     self%zmin=zmax
     self%zmax=zmin
  else if (zmin==zmax) then
     call log_error(m_name,s_name,16,error_fatal,'zmin should not be equal zmax')
  else
     self%zmin=zmin
     self%zmax=zmax
  end if

end  subroutine beqan_readcon
!---------------------------------------------------------------------
!> calculate field, velocity and flux
subroutine beqan_solovevmove(self,y,ydot,B,psi)

  !! arguments
  type(beqan_t), intent(in) :: self !< type containing field definition parameters
  real(kr8), dimension(6), intent(inout) :: y !< position/velocity array
  real(kr8), dimension(6), intent(inout) :: ydot !< velocity/accel array
  real(kr8), dimension(3),intent(out) :: B !< field at position
  real(kr8), intent(out) :: psi !< flux position

  !! local variables
  character(*), parameter :: s_name='beqan_solovevmove' !< subroutine name
  real(kr8) :: term1 !<   \f$ (R_2^2 - R_1^2)^2 \f$
  real(kr8) :: term2 !<   \f$ R^2 - R_x^2 \f$
  real(kr8) :: term3 !<   \f$ R^2 - R_0^2 \f$
  real(kr8) :: chi !< local variable
  real(kr8), parameter :: q=1.602E-19 !< charge
  real(kr8), parameter :: m=9.11E-31 !< mass
  real(kr8) :: qmratio !< charge over mass

  !!calculate intermediate terms
  term1=(self%eqparms%r2**2-self%eqparms%r1**2)**2
  term2=y(1)**2-self%rxsq
  term3=y(1)**2-self%r0sq
  chi=(term3**2+(y(3)**2/self%esq)*term2)*(4/term1)
  qmratio=q/m

  !! calculate components of field
  B(1)=((-8*self%eqparms%psi0*y(3))/self%esq)*(term2/term1) !< \f$ B_R \f$
  B(2)=((8*y(1)*self%eqparms%psi0)/(term1))*(2*term3+(y(3)**2/self%esq)) !< \f$ B_Z \f$
  B(3)=(self%eqparms%bphi*self%rc)/y(1) !< \f$ B_\zeta \f$
  B=qmratio*B !< set field term to include constant terms of Lorentz eqn

  !! calculate flux
  psi=self%eqparms%psi0*chi

  !! calculate ydot components
  ydot(1)=y(4)
  ydot(2)=y(5)
  ydot(3)=y(6)
  ydot(4)=-B(2)*y(5)+B(3)*y(6)+y(1)*y(5)**2
  ydot(5)=(B(2)*y(4))/(y(1)**2)+(-B(1)*y(6))/(y(1)**2)-(2*y(4)*y(5))/y(1)
  ydot(6)=-B(3)*y(4)+B(1)*y(5)

end subroutine beqan_solovevmove
!---------------------------------------------------------------------
!> calculate analytic field after Solovev
subroutine beqan_solovev(self,y,B,psi)

  !! arguments
  type(beqan_t), intent(in) :: self !< type containing field definition parameters
  real(kr8), dimension(3), intent(in) :: y !< position in cyls (R,.,Z)
  real(kr8), dimension(3),intent(out) :: B !< field at position
  real(kr8), intent(out) :: psi !< flux position

  !! local variables
  character(*), parameter :: s_name='beqan_solovev' !< subroutine name
  real(kr8) :: term1 !<   \f$ (R_2^2 - R_1^2)^2 \f$
  real(kr8) :: term2 !<   \f$ R^2 - R_x^2 \f$
  real(kr8) :: term3 !<   \f$ R^2 - R_0^2 \f$
  real(kr8) :: chi !< local variable

  !!calculate intermediate terms
  term1=(self%eqparms%r2**2-self%eqparms%r1**2)**2
  term2=y(1)**2-self%rxsq
  term3=y(1)**2-self%r0sq
  chi=(term3**2+(y(3)**2/self%esq)*term2)*(4/term1)

  !! calculate components of field
  B(1)=((-8*self%eqparms%psi0*y(3))/self%esq)*(term2/term1) !< \f$ B_R  \f$
  B(2)=((8*y(1)*self%eqparms%psi0)/(term1))*(2*term3+(y(3)**2/self%esq)) !< \f$ B_Z \f$
  B(3)=(self%eqparms%bphi*self%rc) !< \f$ B_\zeta \f$ note multiplied by \f$ R \f$ hence \f$ (=I=RB_T) \f$

  !! calculate flux
  psi=self%eqparms%psi0*chi

end subroutine beqan_solovev
!---------------------------------------------------------------------
!> user defined analytic field
subroutine beqan_userdefined(self,y,B,psi)

  !! arguments
  type(beqan_t), intent(in) :: self !< type containing field definition parameters
  real(kr8), dimension(3), intent(in) :: y !< position in cyls (R,.,Z)
  real(kr8), dimension(3),intent(out) :: B !< field at position
  real(kr8), intent(out) :: psi !< flux position

  !! local variables
  character(*), parameter :: s_name='beqan_userdefined' !< subroutine name

  !! user defines field here
  B=0
  psi=0

end subroutine beqan_userdefined
!>---------------------------------------------------------------------
subroutine beqan_generic(self,y,B,psi)

  !! arguments
  type(beqan_t), intent(in) :: self !< type containing field definition parameters
  real(kr8), dimension(3), intent(in) :: y !< position in cyls (R,.,Z)
  real(kr8), dimension(3),intent(out) :: B !< field at position
  real(kr8), intent(out) :: psi !< flux position

  !! local variables
  character(*), parameter :: s_name='beqan_generic' !< subroutine name

  !! select analytic formula
  formula_chosen: select case (self%eqparms%formula)
  case('solovev')
     call beqan_solovev(self,y,B,psi)
  case('userdefined')
     call beqan_userdefined(self,y,B,psi)
  end select formula_chosen

end subroutine beqan_generic
!---------------------------------------------------------------------
!> open new file
subroutine beqan_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer, intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='beqan_initwrite' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  character(len=80) :: outputfile !< output file name

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then channel=i exit end if end do noutba=i

  !! open file
  outputfile=trim(fileroot)//"_beqan.out"
  call log_value("Control data file",trim(outputfile))
  call misc_getfileunit(noutba)
  open(unit=noutba,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutba,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if
  if (present(channel)) channel=noutba

end subroutine beqan_initwrite
!---------------------------------------------------------------------
!> write beqan data
subroutine beqan_write(self,kout)

  !! arguments
  type(beqan_t), intent(in) :: self   !< beqan data structure
  integer, intent(in), optional :: kout   !< output channel for beqan data structure

  !! local
  character(*), parameter :: s_name='beqan_write' !< subroutine name
  integer :: iout   !< output channel for beqan data structure

  !! sort out unit
  if(present(kout)) then
     iout=kout
  else
     iout=noutba
  end if

  write(iout,*,iostat=status) 'equilibrium_formula'
  call log_write_check(m_name,s_name,1,status)
  write(iout,*,iostat=status) self%eqparms%formula
  call log_write_check(m_name,s_name,2,status)
  write(iout,*,iostat=status) 'psi0'
  call log_write_check(m_name,s_name,3,status)
  write(iout,*,iostat=status) self%eqparms%psi0
  call log_write_check(m_name,s_name,4,status)
  write(iout,*,iostat=status) 'bphi'
  call log_write_check(m_name,s_name,5,status)
  write(iout,*,iostat=status) self%eqparms%bphi
  call log_write_check(m_name,s_name,6,status)
  write(iout,*,iostat=status) 'r1'
  call log_write_check(m_name,s_name,7,status)
  write(iout,*,iostat=status) self%eqparms%r1
  call log_write_check(m_name,s_name,8,status)
  write(iout,*,iostat=status) 'r2'
  call log_write_check(m_name,s_name,9,status)
  write(iout,*,iostat=status) self%eqparms%r2
  call log_write_check(m_name,s_name,10,status)
  write(iout,*,iostat=status) 'rm'
  call log_write_check(m_name,s_name,18,status)
  write(iout,*,iostat=status) self%eqparms%rm
  call log_write_check(m_name,s_name,19,status)
  write(iout,*,iostat=status) 'zm'
  call log_write_check(m_name,s_name,20,status)
  write(iout,*,iostat=status) self%eqparms%zm
  call log_write_check(m_name,s_name,21,status)


  write(iout,*,iostat=status) 'nrpams'
  call log_write_check(m_name,s_name,46,status)
  write(iout,*,iostat=status) self%eqparms%nrpams
  call log_write_check(m_name,s_name,47,status)
  if (self%eqparms%nrpams>0) then
     write(iout,*,iostat=status) 'real_parameters'
     call log_write_check(m_name,s_name,48,status)
     write(iout,*,iostat=status) self%eqparms%rpar
     call log_write_check(m_name,s_name,49,status)
  end if

  write(iout,*,iostat=status) 'nipams'
  call log_write_check(m_name,s_name,50,status)
  write(iout,*,iostat=status) self%eqparms%nipams
  call log_write_check(m_name,s_name,51,status)
  if (self%eqparms%nipams>0) then
     write(iout,*,iostat=status) 'integer_parameters'
     call log_write_check(m_name,s_name,52,status)
     write(iout,*,iostat=status) self%eqparms%npar
     call log_write_check(m_name,s_name,53,status)
  end if

end subroutine beqan_write
!---------------------------------------------------------------------
!> delete object
subroutine beqan_delete(self)

  !! arguments
  type(beqan_t), intent(inout) :: self !< type containing analytic equil parameters
  !! local
  character(*), parameter :: s_name='beqan_delete' !< subroutine name

  formula_deallocate: select case (self%eqparms%formula)
  case('userdefined')
     if (self%eqparms%nrpams>0) deallocate(self%eqparms%rpar)
     if (self%eqparms%nipams>0) deallocate(self%eqparms%npar)
  case default
  end select formula_deallocate

end subroutine beqan_delete
!---------------------------------------------------------------------
!> close file
subroutine beqan_close

  !! local
  character(*), parameter :: s_name='beqan_close' !< subroutine name

  !! close file
  close(unit=ninba,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine beqan_close
!---------------------------------------------------------------------
!> close write file
subroutine beqan_closewrite

  !! local
  character(*), parameter :: s_name='beqan_closewrite' !< subroutine name

  !! close file
  close(unit=noutba,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine beqan_closewrite

end module beqan_m
