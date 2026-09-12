module antest_m

  use antest_h
  use log_m
  use const_numphys_h
  use const_kind_m
  use vfile_m
  use radnsource_h
  use radnsource_m
  use ray_m
  use antenna_h
  use antenna_m
  use beq_h
  use beq_m
  use geobjlist_h
  use geobjlist_m

  implicit none
  private

! public subroutines
  public :: &
  antest_initfile,  & !< open file
  antest_readcon,  & !< read data from file
  antest_solve,  & !< generic subroutine
  antest_userdefined,  & !< user-defined function
  antest_fn, &  !< general external function call
  antest_dia, &  !< object diagnostics to log file
  antest_initwrite, & !< open new file, making up name
  antest_write, &  !< write out object
  antest_writeg, &  !< write out object as gnuplot
  antest_writev, &  !< write out object as vtk
  antest_delete, & !< delete object
  antest_close, & !< close file
  antest_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='antest_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4), save  :: ninat=5     !< control file unit number
  integer(ki4), save  :: noutat=6      !< output file unit number
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4)  :: ilog      !< for namelist dump after error

  contains
!---------------------------------------------------------------------
!> open file
subroutine antest_initfile(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer(ki4), intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='antest_initfile' !< subroutine name
  logical :: unitused !< flag to test unit is available

  if (trim(file)=='null') then
     call log_error(m_name,s_name,1,log_info,'null filename ignored')
     return
  end if

  !! get file unit
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        ninat=i
        if (present(channel)) channel=i
        exit
     end if
  end do

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  open(unit=ninat,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,2,error_fatal,'Cannot open control data file')
     stop
  end if

end subroutine antest_initfile
!---------------------------------------------------------------------
!> read data from file
subroutine antest_readcon(selfn,channel)

  !! arguments
  type(atnumerics_t), intent(out) :: selfn !< type which data will be assigned to
  integer(ki4), intent(in),optional :: channel   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='antest_readcon' !< subroutine name
  character(len=80) :: antest_formula !< formula to be used
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  real(kr8) :: power_split !< variable with meaningful name
  real(kr8) :: number_of_samples_in_theta !< self-explanatory local
  real(kr8) :: number_of_samples_in_phi !< self-explanatory local
  real(kr8) :: case_select_for_sampling_pattern !< self-explanatory local
  real(kr8) :: number_of_samples_along_ray !< self-explanatory local
  real(kr8) :: terminating_distance_of_ray !< self-explanatory local

  real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters  !< self-explanatory variable
  integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters  !< self-explanatory variable
  integer(ki4) :: number_of_real_parameters  !< self-explanatory variable
  integer(ki4) :: number_of_integer_parameters  !< self-explanatory variable

  !! antest parameters
  namelist /antestparameters/ &
 &number_of_samples_in_theta , &
 &number_of_samples_in_phi , &
 &case_select_for_sampling_pattern , &
 &number_of_samples_along_ray, &
 &terminating_distance_of_ray, &
 &power_split, antest_formula, &
 &general_real_parameters, number_of_real_parameters, &
 &general_integer_parameters, number_of_integer_parameters

  !! set default antest parameters
  power_split=0.5_kr8

  antest_formula='unset'
  general_real_parameters=0
  general_integer_parameters=0
  number_of_real_parameters=0
  number_of_integer_parameters=0
  number_of_samples_in_theta = 25
  number_of_samples_in_phi = 25
  case_select_for_sampling_pattern = 1
  number_of_samples_along_ray = 30
  terminating_distance_of_ray = 100

  if(present(channel).AND.channel/=0) then
     !! assume unit already open and reading infile
     ninat=channel
  end if

  !!read antest parameters
  read(ninat,nml=antestparameters,iostat=status)
  if(status/=0) then
     !!dump namelist contents to logfile to assist error location
     print '("Fatal error reading antest parameters")'
     call log_getunit(ilog)
     write(ilog,nml=antestparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading antest parameters')
  end if

  call lowor(antest_formula,1,len_trim(antest_formula))
  !! check for valid data

  formula_chosen: select case (antest_formula)
  case('unset','exp')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,11,error_fatal,'power_split must be >=0 and <=1')

  case('expdouble')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,21,error_fatal,'power_split must be >=0 and <=1')

  case('userdefined')
     if(number_of_real_parameters<0) &
 &   call log_error(m_name,s_name,44,error_fatal,'number of real parameters must be >=0')
     if(number_of_real_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of real parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,45,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters<0) &
 &   call log_error(m_name,s_name,46,error_fatal,'number of integer parameters must be >=0')
     if(number_of_integer_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of integer parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,47,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters==0.AND.number_of_real_parameters==0) &
 &   call log_error(m_name,s_name,48,error_fatal,'no parameters set')

  end select formula_chosen

  !! store values
  selfn%formula=antest_formula

  selfn%f=power_split

  !! allocate arrays and assign

  selfn%nrpams=number_of_real_parameters
  selfn%nipams=number_of_integer_parameters

  formula_allocate: select case (antest_formula)

  case('userdefined')
     if (number_of_real_parameters>0) then
        allocate(selfn%rpar(number_of_real_parameters), stat=status)
        call log_alloc_check(m_name,s_name,65,status)
        selfn%rpar=general_real_parameters(:number_of_real_parameters)
     end if
     if (number_of_integer_parameters>0) then
        allocate(selfn%npar(number_of_integer_parameters), stat=status)
        call log_alloc_check(m_name,s_name,66,status)
        selfn%npar=general_integer_parameters(:number_of_integer_parameters)
     end if
  case default
  end select formula_allocate
  selfn%ntheta = number_of_samples_in_theta
  selfn%nphi = number_of_samples_in_phi
  selfn%option = case_select_for_sampling_pattern
  selfn%nraysample = number_of_samples_along_ray
  selfn%termdistance = terminating_distance_of_ray

end  subroutine antest_readcon
!---------------------------------------------------------------------
!> generic subroutine
subroutine antest_solve(self,beq,wall,radnsource,ray)

  !! arguments
  type(antest_t), intent(inout) :: self !< module object
  type(beq_t), intent(inout) :: beq !< beq data
  type(geobjlist_t), intent(inout) :: wall !< triangulated wall segment
  type(radnsource_t), intent(inout) :: radnsource !< radnsource data
  type(ray_t), intent(inout) :: ray !< ray data
  !! local
  character(*), parameter :: s_name='antest_solve' !< subroutine name


  !call antenna_align(self%antenna,hemisphere)
  !call antenna_normal(self%antenna,hemisphere)
  call antenna_walltest(self%antenna,wall,beq,ray)
end subroutine antest_solve
!---------------------------------------------------------------------
!> output to log file
subroutine antest_dia(self)

  !! arguments
  type(antest_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='antest_dia' !< subroutine name

  call log_value("power ",self%totpow)
  call log_value("ntheta ", self%n%ntheta)
  call log_value("nphi ", self%n%nphi)
  call log_value("nphi ", self%n%option)

end subroutine antest_dia
!---------------------------------------------------------------------
!> userdefined function
function antest_userdefined(self,psi)

  !! arguments
  type(antest_t), intent(in) :: self !< module object
  real(kr8) :: antest_userdefined !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='antest_userdefined' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  integer(ki4) :: ilocal !< local integer variable

  zpos=psi
  pow=0._kr8
  !> user defines \Tt{pow} here
  !! .....
  !! return antest
  antest_userdefined=pow

end function antest_userdefined
!---------------------------------------------------------------------
!> general external function call
function antest_fn(self,psi)

  !! arguments
  type(antest_t), intent(in) :: self !< module object
  real(kr8) :: antest_fn !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='antest_fn' !< subroutine name
  real(kr8) :: pow !< local variable

  pow=0._kr8
  !! select antest
  formula_chosen: select case (self%n%formula)
  case('userdefined')
     pow=antest_userdefined(self,psi)
  end select formula_chosen

  !! return antest
  antest_fn=pow

end function antest_fn
!---------------------------------------------------------------------
!> open new file, making up name
subroutine antest_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer(ki4), intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='antest_initwrite' !< subroutine name
  logical :: unitused !< flag to test unit is available
  character(len=80) :: outputfile !< output file name

  !! get file unit
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        if (present(channel)) channel=i
        exit
     end if
  end do
  noutat=i

  !! open file
  outputfile=trim(fileroot)//"_antest.out"
  call log_value("Control data file",trim(outputfile))
  open(unit=noutat,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutat,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if

end subroutine antest_initwrite
!---------------------------------------------------------------------
!> write antest data
subroutine antest_write(self,channel)

  !! arguments
  type(antest_t), intent(in) :: self   !< antest data structure
  integer(ki4), intent(in), optional :: channel   !< output channel for antest data structure

  !! local
  character(*), parameter :: s_name='antest_write' !< subroutine name
  integer(ki4) :: iout   !< output channel for antest data structure

  !! sort out unit
  if(present(channel)) then
     iout=channel
  else
     iout=noutat
  end if

  !! write to .out file

  write(iout,*,iostat=status) 'antest_formula'
  call log_write_check(m_name,s_name,18,status)
  write(iout,*,iostat=status) self%n%formula
  call log_write_check(m_name,s_name,19,status)

  write(iout,*,iostat=status) 'ntheta'
  call log_write_check(m_name,s_name,20,status)
  write(iout,*,iostat=status)  self%n%ntheta                         !! ntheta
  call log_write_check(m_name,s_name,21,status)

  write(iout,*,iostat=status) 'nphi'
  call log_write_check(m_name,s_name,22,status)
  write(iout,*,iostat=status)  self%n%nphi                           !! nphi
  call log_write_check(m_name,s_name,23,status)


  !  write(iout,*,iostat=status) 'f'
  !  call log_write_check(m_name,s_name,26,status)
  !  write(iout,*,iostat=status) self%n%f
  !  call log_write_check(m_name,s_name,27,status)

  write(iout,*,iostat=status) 'nrpams'
  call log_write_check(m_name,s_name,46,status)
  write(iout,*,iostat=status) self%n%nrpams
  call log_write_check(m_name,s_name,47,status)

  if (self%n%nrpams>0) then
     write(iout,*,iostat=status) 'real_parameters'
     call log_write_check(m_name,s_name,48,status)
     write(iout,*,iostat=status) self%n%rpar
     call log_write_check(m_name,s_name,49,status)
  end if
  write(iout,*,iostat=status) 'nipams'
  call log_write_check(m_name,s_name,50,status)
  write(iout,*,iostat=status) self%n%nipams
  call log_write_check(m_name,s_name,51,status)
  if (self%n%nipams>0) then
     write(iout,*,iostat=status) 'integer_parameters'
     call log_write_check(m_name,s_name,52,status)
     write(iout,*,iostat=status) self%n%npar
     call log_write_check(m_name,s_name,53,status)
  end if

end subroutine antest_write
!---------------------------------------------------------------------
!> write object data as gnuplot
subroutine antest_writeg(self,select,channel)

  !! arguments
  type(antest_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for antest data structure

  !! local
  character(*), parameter :: s_name='antest_writeg' !< subroutine name
  integer(ki4) :: iout   !< output channel for antest data structure

  call log_error(m_name,s_name,1,log_info,'gnuplot file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine antest_writeg
!---------------------------------------------------------------------
!> write object data as vtk
subroutine antest_writev(self,select,kplot)

  !! arguments
  type(antest_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in) :: kplot   !< output channel for antest data structure

  !! local
  character(*), parameter :: s_name='antest_writev' !< subroutine name
  integer(ki4) :: iout   !< output channel for antest data structure
  integer(ki4) :: ifirst=1   !< unity for first call to vfile_rscalarwrite

  plot_type: select case(select)
  case('cartesian')

  case default

     call geobjlist_writev(self%geobjl,'geometry',kplot)
     call vfile_rscalarwrite(self%pow,self%geobjl%ng,'Q-avg','CELL',kplot,ifirst)

  end select plot_type

  call log_error(m_name,s_name,1,log_info,'vtk file produced')

end subroutine antest_writev
!---------------------------------------------------------------------
!> close write file
subroutine antest_closewrite

  !! local
  character(*), parameter :: s_name='antest_closewrite' !< subroutine name

  !! close file
  close(unit=noutat,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine antest_closewrite
!---------------------------------------------------------------------
!> delete object
subroutine antest_delete(self)

  !! arguments
  type(antest_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='antest_delete' !< subroutine name

  formula_deallocate: select case (self%n%formula)
  case('userdefined')
     if (self%n%nrpams>0) deallocate(self%n%rpar)
     if (self%n%nipams>0) deallocate(self%n%npar)
  case default
  end select formula_deallocate

  call antenna_delete(self%antenna)
  call geobjlist_delete(self%geobjl)

end subroutine antest_delete
!---------------------------------------------------------------------
!> close file
subroutine antest_close

  !! local
  character(*), parameter :: s_name='antest_close' !< subroutine name

  !! close file
  close(unit=ninat,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine antest_close

end module antest_m
