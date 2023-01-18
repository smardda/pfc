!> @addtogroup groupname8
!> @{
module powcal_m
!> @}
  use const_kind_m
  use const_numphys_h
  use control_h
  use log_m
  use spl2d_m
  use spl3d_m
  use powcal_h
  use powres_h
!     use powres_m
  use geobjlist_h
  use position_h
  use fmesh_h
  use fmesh_m
  use beq_h
  use vfile_m
  use position_m
  use beq_m
  use btree_m
  use powelt_m
  use geobjlist_m
  use nrsolve_m
  use pcontrol_h
  use termplane_h
  use termplane_m
  use edgprof_m
  use edgprof_h
  use smitermpi_h
  use smitermpi_m
  
  implicit none
  private

! public subroutines
  public :: &
  powcal_readv, &   !< read (vtk) powcal data
  powcal_init, &   !< initialise  powcal data structure
  powcal_readcon, &   !< namelist input of controls
  powcal_move, &   !< loop over entries in results
  powcal_calc, &   !< coordinate calculation of power deposition
  powcal_writev, &   !< write (vtk) powcal data
  powcal_refine, &   !< set refinement level
  powcal_write, &  !< output powcal data
  powcal_xferpow, & !< transfer power to shadow wall geometry
  powcal_delete  !< delete powcal data

!public variables

! private types

! private variables
  logical :: double_null=.FALSE. !<  if .TRUE., read in values required for double null
  character(*), parameter :: m_name='powcal_m' !< module name
  integer   :: status   !< error status
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  integer(ki4) :: level=1   !< refinement level, dynamic if inlevel>1
  integer(ki4) :: inlevel=1   !< number of subelements at refinement level
  integer(ki4) :: imlevel=1   !< number of subelements at previous refinement level +1
  integer :: nin   !< input channel for geobj data
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: idum !< dummy integer
  real(kr8) :: zdum !< dummy real
  logical :: iltest !< logical flag
  integer(ki4) :: infilelevel   !< max refinement level from input
  integer(ki4), dimension(:), allocatable :: work !< save objlist variable
  real(kr4), dimension(:), allocatable :: rwork !< real workspace

  contains
!---------------------------------------------------------------------
!> read (vtk) powcal data
subroutine powcal_readv(self,vtkfile)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  character(*), intent(in)       :: vtkfile !< name of input file


  !! local
  character(*), parameter :: s_name='powcal_readv' !< subroutine name
  integer(ki4) :: ifmt !< numeric data format default

  call geobjlist_read(self%powres%geobjl,vtkfile)
  idum=0
  ifmt=1
!HJL  call position_readveclis(self%powres%vecx,vtkfile,'Xcart',idum,ifmt,1)
  self%powres%vecx%np=0 
  call position_readveclis(self%powres%vecx,vtkfile,'Xcart',idum,ifmt)
  self%powres%vecb%np=self%powres%vecx%np
  idum=0
  call position_readveclis(self%powres%vecb,vtkfile,'Bcart',idum,ifmt,1)

end subroutine powcal_readv
!---------------------------------------------------------------------
!> initialise  powcal data structure
subroutine powcal_init(self,numerics,plot,gnumerics,rxpt,zxpt,psixpt)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< object data structure
  type(pnumerics_t), intent(in) :: numerics   !< object control data structure
  type(pplots_t), intent(in) :: plot   !< vtk plot selectors
  type(numerics_t), intent(inout), optional :: gnumerics   !< object control data structure
  real(kr8) :: rxpt(2)   !<  position x point \f$ R \f$
  real(kr8) :: zxpt(2)   !<  position x point \f$ Z \f$
  real(kr8) :: psixpt(2) !<  \f$ \psi \f$ at x point

  !! local
  character(*), parameter :: s_name='powcal_init' !< subroutine name
  real(kr8) :: zsign !< sign factor
  real(kr8) :: r0 !< inner wall R position at z=zcen
  real(kr8) :: r5 !< outer wall R position at z=zcen
  real(kr4), dimension(3) :: pos !< output position vector
  type(powelt_t) :: elt   !< element

  ! initialise control structure
  self%n=numerics

  ! set sign of psi coordinate transform
  call beq_rsigset(self%powres%beq)

  ! initialise powelt data
  call powelt_init

  ! check refinement level
  infilelevel=max(1,self%powres%geobjl%nparam(1))
  call log_value("data file has refinement level ",infilelevel)
  call log_value("input demands refinement level ",self%n%nlevel)
  if (self%n%nlevel>infilelevel) then
     call log_error(m_name,s_name,1,error_warning,'Requested refinement level too high, using file value')
     self%n%nlevel=infilelevel
  end if

  ! calculate number of level 1 triangles
  self%powres%npowe=self%powres%geobjl%ng/powelt_table(infilelevel,2)

  ! this would ideally be called from the powcal program unit but needs to be here for now
  call powcalmpi_decomp(self%decomp,self%powres%npowe,self%powres%geobjl%ng)

  ! set size of pow array
  if(self%decomp%myrankpowcal.eq.0) then
     allocate(self%powres%pow(self%powres%geobjl%ng), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     allocate(self%powres%pow(self%decomp%firstng:self%decomp%lastng), stat=status)
  endif
  !calculate R0 and R5, needed for 4 region case
  do i=self%decomp%firstpowe,self%decomp%lastpowe
     elt%ie=i
     do j=imlevel,inlevel
        elt%je=j
        call powelt_pos(elt,self,pos)
        !look for inner wall R position at z=zcen
        if(pos(1)<self%powres%beq%n%rcen) then
           if(ABS(pos(2)-self%powres%beq%n%zcen)<1E-2) R0=pos(1)
        !look for outer wall R position at z=zcen
        else  
           if(ABS(pos(2)-self%powres%beq%n%zcen)<1E-2) R5=pos(1)    
        end if
     end do
  end do
  WRITE(*,*) 'R0 = ',R0,' R5 = ',R5
  ! set size of power statistics arrays
  if(self%decomp%myrankpowcal.eq.0) then
     allocate(self%powres%pows(self%powres%npowe), stat=status)
     call log_alloc_check(m_name,s_name,2,status)
     allocate(self%powres%powa(self%powres%npowe), stat=status)
     call log_alloc_check(m_name,s_name,3,status)
  else
     allocate(self%powres%pows(self%decomp%firstpowe:self%decomp%lastpowe), stat=status)
     allocate(self%powres%powa(self%decomp%firstpowe:self%decomp%lastpowe), stat=status)
  endif
  psixpt=self%powres%beq%psixptarr
  zsign=beq_rsig()
  if(self%powres%beq%number_regions==1) then
    call edgprof_factors(self%edgprof,&
       &self%powres%beq%rbdry,self%powres%beq%bpbdry,self%powres%beq%btotbdry,zsign)
  else if(self%powres%beq%number_regions==2) then
     call edgprof_factors_2_region(self%edgprof,&
          &self%powres%beq%rbdry,self%powres%beq%bpbdry,self%powres%beq%btotbdry,zsign,&
          self%powres%beq%btotbdryarr,self%powres%beq%rbdryarr,&
          self%powres%beq%bpbdryarr,self%powres%beq%n%rcen)
  else if(self%powres%beq%number_regions==4) then
     call edgprof_factors_4_region(self%edgprof,&
          &self%powres%beq%rbdry,self%powres%beq%bpbdry,self%powres%beq%btotbdry,zsign,&
          self%powres%beq%btotbdryarr,self%powres%beq%rbdryarr,&
          self%powres%beq%bpbdryarr,r0,r5)
  end if
  ! track controls
  self%powres%flinx=plot%flinx
  self%powres%flinends=plot%flinends
  self%powres%flinm=plot%flinm

  calcn_type: select case (self%n%caltype)
  case('afws','local','msus','global')
     ! set size of psista, angle and other diagnostics arrays
     if(self%decomp%myrankpowcal.eq.0) then
        allocate(self%powres%psista(self%powres%geobjl%ng), stat=status)
        call log_alloc_check(m_name,s_name,4,status)
        allocate(self%powres%angle(self%powres%geobjl%ng), stat=status)
        call log_alloc_check(m_name,s_name,5,status)
        allocate(self%powres%zbdotn(self%powres%geobjl%ng), stat=status)
        call log_alloc_check(m_name, s_name,6,status)
        allocate(self%powres%bfield(3, self%powres%geobjl%ng), stat=status)
        call log_alloc_check(m_name, s_name,7,status)
        allocate(self%powres%norm(3, self%powres%geobjl%ng), stat=status)
        call log_alloc_check(m_name, s_name,8,status)
        allocate(self%powres%lenpath(self%powres%geobjl%ng), stat=status)
        call log_alloc_check(m_name, s_name,9,status)
     else
        allocate(self%powres%psista(self%decomp%firstng:self%decomp%lastng), stat=status)
        allocate(self%powres%angle(self%decomp%firstng:self%decomp%lastng), stat=status)
        allocate(self%powres%zbdotn(self%decomp%firstng:self%decomp%lastng), stat=status)
        allocate(self%powres%lenpath(self%decomp%firstng:self%decomp%lastng), stat=status)
        allocate(self%powres%bfield(3, self%decomp%firstng:self%decomp%lastng), stat=status)
        allocate(self%powres%norm(3, self%decomp%firstng:self%decomp%lastng), stat=status)
     endif
     self%powres%psista=0
     self%powres%angle=0
     self%powres%zbdotn=0
     self%powres%bfield=0
     self%powres%norm=0
  case default
  end select calcn_type
  ! set size and zero of mask array
  if(self%decomp%myrankpowcal.eq.0) then
     allocate(self%powres%pmask(self%powres%geobjl%ng), stat=status)
     call log_alloc_check(m_name,s_name,9,status)
  else
     allocate(self%powres%pmask(self%decomp%firstng:self%decomp%lastng), stat=status)
  endif
  self%powres%pmask=0

  field_type: select case (self%powres%beq%n%fldspec)
  case(1)
     ! limits on psi-theta interpolation
     call spl2d_ptlimits(self%powres%beq%rjac,self%powres%psimin,self%powres%psimax,&
     self%powres%thetamin,self%powres%thetamax)
  case(2,3)
  end select field_type

  if (self%n%shadow/=0) then
     if (self%n%shadow==-1) then
        ! testing, make quantising transform the identity
        gnumerics%geobj_coord_tfm%hmin=1
        gnumerics%geobj_coord_tfm%rhmin=1
        gnumerics%geobj_coord_tfm%offvec=0
        gnumerics%geobj_coord_tfm%nqtfm=1
        ! and position transform
        gnumerics%position_coord_tfm%scale=1
        gnumerics%position_coord_tfm%offset=0
        gnumerics%position_coord_tfm%matrix(1,:)=(/1,0,0/)
        gnumerics%position_coord_tfm%matrix(2,:)=(/0,1,0/)
        gnumerics%position_coord_tfm%matrix(3,:)=(/0,0,1/)
        gnumerics%position_coord_tfm%ntfm=1
     end if
     ! quantise relevant quantities
     call powcal_quant(self,gnumerics,rxpt,zxpt,psixpt)
  end if

end subroutine powcal_init
!---------------------------------------------------------------------
!> quantise relevant quantities
subroutine powcal_quant(self,gnumerics,rxpt,zxpt,psixpt)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< object data structure
  type(numerics_t), intent(in) :: gnumerics   !< object control data structure
  real(kr8) :: rxpt(2)   !<  position x point \f$ R \f$
  real(kr8) :: zxpt(2)   !<  position x point \f$ Z \f$
  real(kr8) :: psixpt(2) !<  \f$ \psi \f$ at x point
    
  !! local
  character(*), parameter :: s_name='powcal_quant' !< subroutine name
  type(posveclis_t) :: zposl !< local variable
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zposq !< local variable
  integer(ki4) :: inqtfm   !< save transform type
  integer(ki4) :: inkn1   !< number of knots in 1-direction
  integer(ki4) :: inkn2   !< number of knots in 2-direction
  integer(ki4) :: inkn   !< maximum number of knots in either direction
  real(kr8), dimension(3)  :: zpar !< scaling parameters
  real(kr8) :: zh1 !< scaling parameter
  real(kr8) :: zh2 !< scaling parameter
  real(kr8) :: zh3 !< scaling parameter
  real(kr8) :: zo1 !< offset parameter
  real(kr8) :: zo3 !< offset parameter
  real(kr8) :: zepsdif !<  \f$ \epsilon \f$ for \f$ \zeta \f$ boundary crossing
  real(kr8) :: zrnns !< \f$ \zeta \f$ scaling
  real(kr8) :: zarnns  !< \f$ \zeta \f$ component scaling
  integer(ki4) :: itermp !< total number of termination planes
  integer(ki4) :: idir !< coordinate direction \f$ (R,Z,\xi) \f$ in order \f$ 1,2,3 \f$
  integer(ki4) :: idirc !< coordinate direction of auxilliary condition
  integer(ki4), dimension(2) :: indbox   !<  extents of box arrays
  integer(ki4) :: ibdim !<  second dimension of boxr,z arrays
  integer(ki4) :: ntyps !<  upper or lower skylight type
  integer(ki4) :: interm   !< index of dynamical entry in termplane
  integer(ki4) :: lc   !< used in line count of double null case
    
  rxpt=self%powres%beq%rxptarr
  zxpt=self%powres%beq%zxptarr
     
  zpar=0
  self%powres%beq%domlen=0

  self%powres%geobjl%tfmdata=gnumerics%position_coord_tfm
  self%powres%geobjl%quantfm=gnumerics%geobj_coord_tfm
  ! rotate and quantise results geometry positions
  call position_tfmlis(self%powres%geobjl%posl,self%powres%geobjl%tfmdata)
  call position_qtfmlis(self%powres%geobjl%posl,self%powres%geobjl%quantfm)

     allocate(zposl%pos(2),stat=status)
     zposl%pos(1)%posvec=(/rxpt(1), zxpt(1), 0._kr8/)
     zposl%pos(2)%posvec=(/rxpt(2), zxpt(2), 0._kr8/)
     zposl%np=2
     call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
     rxpt(1)=zposl%pos(1)%posvec(1)
     rxpt(2)=zposl%pos(2)%posvec(1)
     zxpt(1)=zposl%pos(1)%posvec(2)
     zxpt(2)=zposl%pos(2)%posvec(2)
     ! domain control
     deallocate(zposl%pos)

  field_type: select case (self%powres%beq%n%fldspec)
  case (1)
     ! mapped (psi-theta-zeta) case
     ! multiplies RHS of ode
     self%powres%qfac=self%powres%geobjl%quantfm%hmin(2)/self%powres%geobjl%quantfm%hmin(3)

     ! limits on psi-theta interpolation, also timestep control
     ! beware that psibdry must not be quantised
     allocate(zposl%pos(6),stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     zposl%pos(1)%posvec=0 ! scaled below
     zposl%pos(2)%posvec=(/self%powres%beq%rjac%org1, self%powres%beq%rjac%org2, 0._kr8/)
     zposl%pos(3)%posvec=(/self%powres%psimin, self%powres%thetamin, 0._kr8/)
     zposl%pos(4)%posvec=(/self%powres%psimax, self%powres%thetamax, self%odes%n%tmax/)
     zposl%pos(5)%posvec=(/self%powres%beq%psiaxis, 0._kr8, 0._kr8/)
     zposl%pos(6)%posvec=(/self%powres%beq%psiqbdry, 0._kr8, 0._kr8/)
     zposl%np=6
     call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
     self%powres%beq%rjac%org1=zposl%pos(2)%posvec(1)
     self%powres%beq%rjac%org2=zposl%pos(2)%posvec(2)
     self%powres%psimin=zposl%pos(3)%posvec(1)
     self%powres%thetamin=zposl%pos(3)%posvec(2)
     self%powres%psimax=zposl%pos(4)%posvec(1)
     self%powres%thetamax=zposl%pos(4)%posvec(2)
     self%powres%beq%domlen(1)=self%powres%psimax-self%powres%psimin
     self%powres%beq%domlen(2)=self%powres%thetamax-self%powres%thetamin
     self%odes%n%tmax=zposl%pos(4)%posvec(3)
     self%powres%beq%psiaxis=zposl%pos(5)%posvec(1)
     self%powres%beq%psiqbdry=zposl%pos(6)%posvec(1)

     zpos%posvec=(/self%powres%beq%rjac%h1, self%powres%beq%rjac%h2 , self%odes%n%dt0/)
     inqtfm=self%powres%geobjl%quantfm%nqtfm
     ! case 1 just scales
     self%powres%geobjl%quantfm%nqtfm=1
     zposq=position_qtfm(zpos,self%powres%geobjl%quantfm)
     self%powres%beq%rjac%h1=zposq%posvec(1)
     self%powres%beq%rjac%rh1=1/zposq%posvec(1)
     self%powres%beq%rjac%rh2=1/zposq%posvec(2)
     self%powres%beq%rjac%h2=zposq%posvec(2)
     self%odes%n%dt0=zposq%posvec(3)
     self%powres%geobjl%quantfm%nqtfm=inqtfm
     deallocate(zposl%pos)

     ! now knots
     inkn1=self%powres%beq%rjac%n1p+self%powres%beq%rjac%nord
     inkn2=self%powres%beq%rjac%n2p+self%powres%beq%rjac%nord

     inkn=max(inkn1,inkn2)
     allocate(zposl%pos(inkn),stat=status)
     call log_alloc_check(m_name,s_name,20,status)
     do j=1,inkn
        zposl%pos(j)%posvec=0
     end do
     do j=1,inkn1
        zposl%pos(j)%posvec(1)=self%powres%beq%rjac%knot1(j)
     end do
     do j=1,inkn2
        zposl%pos(j)%posvec(2)=self%powres%beq%rjac%knot2(j)
     end do
     zposl%np=inkn
     call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
     do j=1,inkn1
        self%powres%beq%rjac%knot1(j)=zposl%pos(j)%posvec(1)
     end do
     do j=1,inkn2
        self%powres%beq%rjac%knot2(j)=zposl%pos(j)%posvec(2)
     end do
     deallocate(zposl%pos)

  case default
     ! R-Z-xi or R-Z-zeta case

     ! limits on R-Z interpolation
     call spl2d_ptscale(self%powres%beq%psi,self%powres%geobjl%quantfm)
     call spl2d_ptscale(self%powres%beq%dpsidr,self%powres%geobjl%quantfm)
     call spl2d_ptscale(self%powres%beq%dpsidz,self%powres%geobjl%quantfm)
     field_typeinr: select case (self%powres%beq%n%fldspec)
     case (3)
        call spl2d_ptscale(self%powres%beq%rispldr,self%powres%geobjl%quantfm)
        call spl2d_ptscale(self%powres%beq%rispldz,self%powres%geobjl%quantfm)
     end select field_typeinr
     
     ! domain control
     allocate(zposl%pos(6),stat=status)
     call log_alloc_check(m_name,s_name,30,status)
     zposl%pos(1)%posvec=(/0._kr8, 0._kr8, self%powres%beq%n%ximin/)
     zposl%pos(2)%posvec=(/0._kr8, 0._kr8, self%powres%beq%n%ximax/)
     zepsdif=const_epsbdry*(self%powres%beq%n%ximax-self%powres%beq%n%ximin)
     zposl%pos(3)%posvec=(/0._kr8, 0._kr8, self%powres%beq%n%ximin+zepsdif/)
     zposl%pos(4)%posvec=(/0._kr8, 0._kr8, self%powres%beq%n%ximax-zepsdif/)
     zposl%pos(5)%posvec=(/0._kr8, self%odes%n%termcon(3), 0._kr8/)
     zposl%pos(6)%posvec=(/self%powres%beq%n%rcen, self%powres%beq%n%zcen, 0._kr8/)
     zposl%np=6
     call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
     self%powres%beq%n%ximin=zposl%pos(1)%posvec(3)
     self%powres%beq%n%ximax=zposl%pos(2)%posvec(3)
     self%powres%beq%domlen(3)=self%powres%beq%n%ximax-self%powres%beq%n%ximin
     self%powres%beq%ximinp=zposl%pos(3)%posvec(3)
     self%powres%beq%ximaxm=zposl%pos(4)%posvec(3)
     self%odes%n%termcon(3)=zposl%pos(5)%posvec(2)
     self%powres%beq%rcenq=zposl%pos(6)%posvec(1)
     self%powres%beq%zcenq=zposl%pos(6)%posvec(2)
     !
     inqtfm=self%powres%geobjl%quantfm%nqtfm
     ! case 1 only scales, as required here
     self%powres%geobjl%quantfm%nqtfm=1

     zpar(3)=self%powres%beq%nzets
     zposl%pos(1)%posvec=(/0._kr8, 0._kr8, zpar(3)*self%odes%n%tmax/)
     zposl%pos(2)%posvec=(/0._kr8, 0._kr8, zpar(3)*self%odes%n%dt0/)
     zposl%np=2
     call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
     self%odes%n%tmax=zposl%pos(1)%posvec(3)
     self%odes%n%dt0=zposl%pos(2)%posvec(3)
     ! restore
     self%powres%geobjl%quantfm%nqtfm=inqtfm
     deallocate(zposl%pos)

     ! assuming termplane only in non-mapped field type
     if (self%n%ltermplane) then
        ! scale termination planes
        itermp=self%n%termp%ntermplane
        allocate(zposl%pos(2*itermp),stat=status)
        call log_alloc_check(m_name,s_name,31,status)
        do j=1,itermp
           zposl%pos(j)%posvec=0
           idir=self%n%termp%termplanedir(j,1)
           zposl%pos(j)%posvec(idir)=self%n%termp%termplane(j,1)
           zposl%pos(itermp+j)%posvec=0
           idirc=self%n%termp%termplanedir(j,4)
           zposl%pos(itermp+j)%posvec(idirc)=self%n%termp%termplane(j,2)
        end do
        zposl%np=2*itermp
        call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
        do j=1,itermp
           idir=self%n%termp%termplanedir(j,1)
           self%n%termp%termplane(j,1)=zposl%pos(j)%posvec(idir)
           idirc=self%n%termp%termplanedir(j,4)
           self%n%termp%termplane(j,2)=zposl%pos(itermp+j)%posvec(idirc)
        end do
        deallocate(zposl%pos)
     end if

     ! multiplies RHS of ode
     self%powres%qfaca=self%powres%geobjl%quantfm%hmin(3)/self%powres%geobjl%quantfm%hmin

     vac_type: select case (self%powres%beq%n%vactype(1:3))
     case('nul')
        if (self%powres%beq%n%mrip/=0) then
           !! special analytic model, has to be fldspec=2
           ! set ripple parameters (before scaling)
           zdum=beq_ripple_h1(1._kr8,0._kr8,0._kr8,-1,&
 &         self%powres%beq%n%mrip,self%powres%beq%ivac,self%powres%beq%n%arip)
           ! set scale function arguments
           zh1=self%powres%geobjl%quantfm%hmin(1)
           zh2=self%powres%geobjl%quantfm%hmin(2)
           zh3=self%powres%geobjl%quantfm%hmin(3)
           zdum=beq_ripple_h1(1._kr8,0._kr8,0._kr8,-1,-1,zh1,zh3)
           zo1=self%powres%geobjl%quantfm%offvec(1)
           zo3=self%powres%geobjl%quantfm%offvec(3)
           zdum=beq_ripple_h1(1._kr8,0._kr8,0._kr8,-1,-2,zo1,zo3)
           zrnns=real(self%powres%beq%n%mrip,kr8)/self%powres%beq%nzets
           zarnns=self%powres%beq%n%arip*zrnns
           zdum=beq_ripple_h1(1._kr8,0._kr8,0._kr8,-1,-3,zrnns,zarnns)
           call nrsolve_spl2dfn(self%powres%beq%psi,self%powres%beq%dpsidr,&
 &         self%powres%beq%dpsidz,zdum,&
 &         zh1,zh2,idum,0)
        else
           field_typeinr2: select case (self%powres%beq%n%fldspec)
           case (3)
              ! no 3-cpt so value of qfaca(3) does not matter
              self%powres%qfaca=self%powres%qfaca/self%powres%beq%nzets
           case default ! should be only fldspec=2
              ! working with f=I in 3-component, so need 1/h_1 factor for extra R
              self%powres%qfaca(3)=1._kr8/self%powres%geobjl%quantfm%hmin(1)
           end select field_typeinr2
        end if
     case('sp3')
        !! has to be fldspec=2 here
        ! ripple defined by spl3d vacuum field
        call spl3d_ptscale(self%powres%beq%vacfld,self%powres%geobjl%quantfm)
        ! here scale consistent with xi update
        self%powres%qfaca(3)=self%powres%beq%nzets
     case('txt')
        !! has to be fldspec=2 here
        ! ripple defined by beqart vacuum field
        call fmesh_qtfm(self%powres%beq%beqart%fmesh,self%powres%geobjl%quantfm)
        ! here scale consistent with xi update
        self%powres%qfaca(3)=self%powres%beq%nzets
     end select vac_type

     ! normalise termination controls
     self%odes%n%termcon(2)=self%odes%n%termcon(2)*self%powres%beq%domlen(3)

     !! skylight
     if (self%n%lskyl) then
        if (self%powres%beq%n%skylpsi) then
           ! overwrite with changed powcalparameters value
           if (self%n%lflagtol) self%powres%skyl%n%toli=self%n%tolin
           ! ctrackrz array
           allocate(zposl%pos(self%powres%beq%nctrack),stat=status)
           call log_alloc_check(m_name,s_name,51,status)
           do j=1,self%powres%beq%nctrack
              zposl%pos(j)%posvec=0
              zposl%pos(j)%posvec(1)=self%powres%beq%ctrackrz(j,1)
              zposl%pos(j)%posvec(2)=self%powres%beq%ctrackrz(j,2)
           end do
           zposl%np=self%powres%beq%nctrack
           call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
           do j=1,self%powres%beq%nctrack
              self%powres%beq%ctrackrz(j,1)=zposl%pos(j)%posvec(1)
              self%powres%beq%ctrackrz(j,2)=zposl%pos(j)%posvec(2)
           end do
           deallocate(zposl%pos)
           ibdim=mod(self%powres%skyl%n%skyltyp,2)
           if (ibdim==2) then
              indbox(1)=max( self%powres%skyl%dimbox(1,1),self%powres%skyl%dimbox(1,2) )
              indbox(2)=max( self%powres%skyl%dimbox(2,1),self%powres%skyl%dimbox(2,2) )
           else
              indbox=self%powres%skyl%dimbox(:,1)
           end if
           !! inboxr,z arrays
           allocate(zposl%pos(indbox(1)),stat=status)
           call log_alloc_check(m_name,s_name,52,status)
           do k=1,ibdim
              do j=1,indbox(1)
                 zposl%pos(j)%posvec=0
                 zposl%pos(j)%posvec(1)=self%powres%skyl%inboxr(j,k)
                 zposl%pos(j)%posvec(2)=self%powres%skyl%inboxz(j,k)
              end do
              zposl%np=indbox(1)
              call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
              do j=1,indbox(1)
                 self%powres%skyl%inboxr(j,k)=zposl%pos(j)%posvec(1)
                 self%powres%skyl%inboxz(j,k)=zposl%pos(j)%posvec(2)
              end do
           end do
           deallocate(zposl%pos)
           !! ouboxr,z arrays
           allocate(zposl%pos(indbox(2)),stat=status)
           call log_alloc_check(m_name,s_name,53,status)
           do k=1,ibdim
              do j=1,indbox(2)
                 zposl%pos(j)%posvec=0
                 zposl%pos(j)%posvec(1)=self%powres%skyl%ouboxr(j,k)
                 zposl%pos(j)%posvec(2)=self%powres%skyl%ouboxz(j,k)
              end do
              zposl%np=indbox(2)
              call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
              do j=1,indbox(2)
                 self%powres%skyl%ouboxr(j,k)=zposl%pos(j)%posvec(1)
                 self%powres%skyl%ouboxz(j,k)=zposl%pos(j)%posvec(2)
              end do
           end do
           deallocate(zposl%pos)
           ntyps=self%powres%skyl%n%skyltyp
           ! allow for divertor skylight
           if (self%n%termp%ntermplane<=0) then
              call log_error(m_name,s_name,60,error_fatal,'termplane(s) must be defined for flux-based  skylight(s)')
           end if
           ! quantise positions
           allocate(zposl%pos(2),stat=status)
           call log_alloc_check(m_name,s_name,30,status)
           zposl%pos(1)%posvec=(/0._kr8, self%powres%beq%zmin, 0._kr8/)
           zposl%pos(2)%posvec=(/0._kr8, self%powres%beq%zmax, 0._kr8/)
           zposl%np=2
           call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
           self%powres%beq%zmin=zposl%pos(1)%posvec(2)
           self%powres%beq%zmax=zposl%pos(2)%posvec(2)
           deallocate(zposl%pos)
           allocate(zposl%pos(1),stat=status)
           ! only scale
           self%powres%geobjl%quantfm%nqtfm=1
           zposl%pos(1)%posvec=(/0._kr8, self%powres%beq%dz, 0._kr8/)
           zposl%np=1
           call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
           self%powres%beq%dz=zposl%pos(1)%posvec(2)
           ! restore
           self%powres%geobjl%quantfm%nqtfm=inqtfm
           deallocate(zposl%pos)
        else
           call log_error(m_name,s_name,61,error_warning,'skylight requested is not defined by flux')
        end if
     else
        ! prevent application of flux skylight
        self%powres%beq%n%skylpsi=.FALSE.
     end if
  end select field_type

end subroutine powcal_quant
!---------------------------------------------------------------------
!> namelist input of controls
subroutine powcal_readcon(selfn,kin,number_of_regions)

  !! arguments
  type(pnumerics_t), intent(inout) :: selfn   !< object control data structure
  integer, intent(in) :: kin   !< input channel for object data structure
  integer(ki4) :: number_of_regions !< number of regions 

  !! local
  character(*), parameter :: s_name='powcal_readcon' !< subroutine name
  real(kr8) :: power_split !< power split ion to electron direction
  real(kr8), dimension(:), allocatable :: decay_length !< power decay length at outer midplane (m)
  real(kr8) , dimension(:), allocatable:: power_loss !< power crossing the separatrix (W)
  real(kr8) , dimension(:), allocatable:: diffusion_length !< (m)
  real(kr8) :: q_parallel0 !< parallel reference flux (W/sqm)
  integer(ki4) :: refine_level   !< level of refinement to be used
  logical :: use_pfc_flux   !< default .FALSE., if true, use value of psi on PFC to evaluate Q
  logical :: flux_from_pfc   !< default .FALSE., if true, use value of psi on PFC to evaluate Q
  integer(ki4) :: shadow_control   !< control shadowing
  character(len=80) :: calculation_type !< calculation type
  integer(ki4) :: number_of_tracks !< nonzero, only compute selected tracks
  integer(ki4), parameter  :: maximum_number_of_tracks=100 !< maximum number of tracks allowed
  integer(ki4), dimension(maximum_number_of_tracks) :: element_numbers !< tracks selected by element number
  logical :: termination_planes !<  termination planes present
  logical :: more_profiles !<  use edgprofparameters
  logical :: power_in_pfr !<  if .TRUE., power in private flux region
  logical :: employ_skyl_pfr !<  if .TRUE., employ skylights in private flux region
  integer(ki4) :: analytic_launch_type !< analytic launch type
  integer(ki4) :: power_on_launch_geo !< power on launch geometry
  integer(ki4) :: power_on_shadow_geo !< power on shadow geometry
  real(kr8), dimension(3) :: ztola !< default tolerance
  !> Distance that assumes intersection with skylight, in quantised units, in each direction
  !!  typically 1/1024 of geometry, default 1.
  real(kr8), dimension(3) :: skylight_tolerance !< .
  logical :: skylight !< activate skylights if present, default .TRUE.

  !! powcal parameters
  namelist /powcalparameters/ &
 &power_split, decay_length, power_loss, refine_level, shadow_control, &
 &use_pfc_flux, &
 &flux_from_pfc, &
 &q_parallel0, &
 &diffusion_length, &
 &calculation_type, &
 &number_of_tracks, element_numbers, &
 &termination_planes, &
 &more_profiles, &
 &power_in_pfr, &
 &employ_skyl_pfr, &
 &skylight, &
 &skylight_tolerance, double_null, &
 &analytic_launch_type, power_on_launch_geo, power_on_shadow_geo


  allocate(decay_length(number_of_regions),power_loss(number_of_regions),diffusion_length(number_of_regions))
  allocate(selfn%lmid(number_of_regions),selfn%ploss(number_of_regions),selfn%sigma(number_of_regions))
  !! set default powcal parameters
  power_split=0.5_kr8
  decay_length=.012_kr8
  diffusion_length=0
  power_loss=7.5e+06_kr8
  refine_level=1
  use_pfc_flux=.FALSE.
  flux_from_pfc=.FALSE.
  shadow_control=0
  !! ignore q_parallel0 if negative
  q_parallel0=-1_kr8

  calculation_type='unset'
  number_of_tracks=0
  if (selfn%mtrack>maximum_number_of_tracks) then
     call log_error(m_name,s_name,1,error_warning,'maximum_number_of_tracks parameter too small')
  end if
  element_numbers=0
  termination_planes=.FALSE.
  more_profiles=.FALSE.
  power_in_pfr=.TRUE.
  employ_skyl_pfr=.TRUE.
  skylight=.TRUE.
  ztola=1._kr8
  skylight_tolerance=ztola
  analytic_launch_type=0
  power_on_launch_geo=0
  power_on_shadow_geo=0

  !!read powcal parameters
  read(kin,nml=powcalparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading powcal parameters")'
     call log_getunit(ilog)
     write(ilog,nml=powcalparameters)
     call log_error(m_name,s_name,2,error_fatal,'Error reading powcal parameters')
  end if

  !! check for valid data
  calcn_type: select case (calculation_type)
  case('afws','local','msus','global','msum','middle')
     ! valid options
     call log_value('calculation type', calculation_type)
  case('unset')
     call log_error(m_name,s_name,5,error_warning,'calculation_type not set, default to local')
     calculation_type='local'
     call log_value('calculation type', calculation_type)
  case default
     call log_error(m_name,s_name,6,error_fatal,'calculation_type not recognised')
  end select calcn_type
  if(power_split<0.OR.power_split>1) &
 &call log_error(m_name,s_name,11,error_fatal,'power_split must be >=0 and <=1')
 do l=1,number_of_regions
  if(decay_length(l)<=0) &
 &call log_error(m_name,s_name,12,error_fatal,'decay_length (m) must be >0')
  if(diffusion_length(l)<0) &
 &call log_error(m_name,s_name,13,error_fatal,'diffusion_length (m) must be >=0')
  if(q_parallel0<=0.AND.power_loss(l)<=0) &
 &call log_error(m_name,s_name,14,error_fatal,'power_loss (W) must be >0')
 end do
  ! need to check against value in input file
  if(refine_level<1.OR.refine_level>3) &
 &call log_error(m_name,s_name,15,error_fatal,'refine_level must be 1, 2 or 3')

  if(analytic_launch_type<0.OR.analytic_launch_type>10) &
 &call log_error(m_name,s_name,23,error_fatal,'analytic_launch_type must be >=0 and <=10')
  if(power_on_launch_geo<0.OR.power_on_launch_geo>10) &
 &call log_error(m_name,s_name,24,error_fatal,'power_on_launch_geo must be >=0 and <=10')
  if(power_on_shadow_geo<0.OR.power_on_shadow_geo>10) &
 &call log_error(m_name,s_name,25,error_fatal,'power_on_shadow_geo must be >=0 and <=10')
  if(number_of_tracks<0) &
 &call log_error(m_name,s_name,26,error_fatal,'number_of_tracks must be non-negative')
  if(number_of_tracks>0.AND.minval(element_numbers)<0) &
 &call log_error(m_name,s_name,27,error_fatal,'element numbers must be non-negative')

  !! store values
  selfn%f=power_split
  selfn%lmid=decay_length
  selfn%ploss=power_loss
  selfn%sigma=diffusion_length
  selfn%nlevel=refine_level
  selfn%usepsista=.NOT.use_pfc_flux.AND..NOT.flux_from_pfc
  selfn%shadow=shadow_control
  selfn%qpara0=q_parallel0

  selfn%caltype=calculation_type
  selfn%ntrack=number_of_tracks
  if (number_of_tracks>0) then
     allocate(selfn%trackno(selfn%ntrack), stat=status)
     call log_alloc_check(m_name,s_name,30,status)
     do j=1,number_of_tracks
        selfn%trackno(j)=element_numbers(j)+1
     end do
  end if
  selfn%ltermplane=termination_planes
  selfn%ledgprof=more_profiles
  selfn%lpfpower=power_in_pfr
  selfn%lskylpfr=employ_skyl_pfr
  selfn%lskyl=skylight
  selfn%lflagtol=ANY(skylight_tolerance/=ztola)
  selfn%tolin=skylight_tolerance
  selfn%nanalau=analytic_launch_type
  selfn%nlaupow=power_on_launch_geo
  selfn%nshapow=power_on_shadow_geo
  if (selfn%ltermplane) then
     call termplane_readcon(selfn%termp,kin)
  else
     ! shut off more detailed tests
     selfn%termp%ntermplane=0
  end if

end subroutine powcal_readcon
!---------------------------------------------------------------------
!> loop over entries in results
subroutine powcal_move(self,zelt,gshadl,btree)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  type(powelt_t) :: zelt   !< power element
  type(geobjlist_t), intent(inout) :: gshadl   !< shadow data structure
  type(btree_t), intent(inout) :: btree !< btree data

  !! local
  character(*), parameter :: s_name='powcal_move' !< subroutine name

  ! check for axisymmetric !mpi1!
  vac_type: select case (self%powres%beq%n%vactype(1:3))
  case('nul')
     if (self%powres%beq%n%mrip/=0) then
        ! not axisymmetric
        ! analytic ripple model, has to be fldspec=2
        ! no termplane criteria
        do i=self%decomp%firstpowe,self%decomp%lastpowe
           zelt%ie=i
           do j=imlevel,inlevel
              zelt%je=j
              call powelt_move2(zelt,self,gshadl,btree)
           end do
           if(self%decomp%myrankpowcal==0 .and. mod(i,10000)==0) then
              write(*,"(i3,a1,f6.3,a1)",advance="no") 100*(i)/self%decomp%mynpowe,"%"
           endif           
        end do
     else
        ! axisymmetric (no ripple field)
        field_type: select case (self%powres%beq%n%fldspec)
        case (1) ! should only be for caltype='local', no termplane
           do i=self%decomp%firstpowe,self%decomp%lastpowe
              zelt%ie=i
              do j=imlevel,inlevel
                 zelt%je=j
                 call powelt_move(zelt,self,gshadl,btree)
              end do
              if(self%decomp%myrankpowcal==0 .and. mod(i,10000)==0) then
                 write(*,"(i3,a1,f6.3,a1)",advance="no") 100*(i)/self%decomp%mynpowe,"%"
              endif
           end do
        case default ! should only be fldspec=3, all caltypes with termplane
           calcn_type: select case (self%n%caltype)
           case ('afws','local')
              do i=self%decomp%firstpowe,self%decomp%lastpowe
                 zelt%ie=i
                 do j=imlevel,inlevel
                    zelt%je=j
                    call powelt_move0(zelt,self,gshadl,btree)
                 end do
                 if(self%decomp%myrankpowcal==0 .and. mod(i,10000)==0) then
                    write(*,"(i3,a1,f6.3,a1)",advance="no") 100*(i)/self%decomp%mynpowe,"%"
                 endif
              end do
           case ('msus','global','msum','middle')
              do i=self%decomp%firstpowe,self%decomp%lastpowe
                 zelt%ie=i
                 do j=imlevel,inlevel
                    zelt%je=j
                    call powelt_move1(zelt,self,gshadl,btree)
                 end do
                 if(self%decomp%myrankpowcal==0 .and. mod(i,10000)==0) then
                    write(*,"(i3,a1,f6.3,a1)",advance="no") 100*(i)/self%decomp%mynpowe,"%"
                 endif
              end do
              if(self%decomp%myrankpowcal==0) write(*,*)
           end select calcn_type
        end select field_type
     end if
     ! remaining cases non-axisymmetric, ripple defined by vacuum file
  case('txt')
     if (self%n%ltermplane) then
        ! more sophisticated termination criteria
        do i=self%decomp%firstpowe,self%decomp%lastpowe
           zelt%ie=i
           do j=imlevel,inlevel
              zelt%je=j
              call powelt_move4qart(zelt,self,gshadl,btree)
           end do
           if(self%decomp%myrankpowcal==0 .and. mod(i,10000)==0) then
              write(*,"(i3,a1,f6.3,a1)",advance="no") 100*(i)/self%decomp%mynpowe,"%"
           endif

        end do
     else
        ! older midplane-based termination criteria
     call log_error(m_name,s_name,1,error_warning,'Invalid option')
     end if
  case default
     if (self%n%ltermplane) then
        ! more sophisticated termination criteria
        do i=self%decomp%firstpowe,self%decomp%lastpowe
           zelt%ie=i
           do j=imlevel,inlevel
              zelt%je=j
              call powelt_move4(zelt,self,gshadl,btree)
           end do
           if(self%decomp%myrankpowcal==0 .and. mod(i,10000)==0) then
              write(*,"(i3,a1,f6.3,a1)",advance="no") 100*(i)/self%decomp%mynpowe,"%"
           endif
        end do
     else
        ! older midplane-based termination criteria
        do i=self%decomp%firstpowe,self%decomp%lastpowe
           zelt%ie=i
           do j=imlevel,inlevel
              zelt%je=j
              call powelt_move3(zelt,self,gshadl,btree)
           end do
           if(self%decomp%myrankpowcal==0 .and. mod(i,10000)==0) then
              write(*,"(i3,a1,f6.3,a1)",advance="no") 100*(i)/self%decomp%mynpowe,"%"
           endif
        end do
     end if
  end select vac_type

end subroutine powcal_move
!---------------------------------------------------------------------
!> coordinate calculation of power deposition (=powres_XX)
subroutine powcal_calc(self,post,gshadl,rxpt,zxpt,psixpt)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  type(powelt_t), intent(inout) :: post   !< powelt data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure
  real(kr8) :: rxpt(2)   !<  position x point \f$ R \f$
  real(kr8) :: zxpt(2)   !<  position x point \f$ Z \f$
  real(kr8) :: psixpt(2) !<  \f$ \psi \f$ at x point
    
  !! local
  character(*), parameter :: s_name='powcal_calc' !< subroutine name

  ! calculate power deposition
  do i=self%decomp%firstpowe,self%decomp%lastpowe
     post%ie=i
     do j=imlevel,inlevel
        post%je=j
        call powelt_dep(post,self,gshadl,self%powres%beq%rcenq,self%powres%beq%zcenq,rxpt,zxpt,psixpt)
     end do
  end do

end subroutine powcal_calc
!---------------------------------------------------------------------
!> coordinate calculation of power statistics (=powres_XX)
subroutine powcal_stat(self,gshadl)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure

  !! local
  character(*), parameter :: s_name='powcal_stat' !< subroutine name
  type(powelt_t) :: zelt   !< powelt data instance

  ! calculate power fluctuation statistic
  zelt%je=inlevel
  do i=self%decomp%firstpowe,self%decomp%lastpowe
     zelt%ie=i
     call powelt_stat(zelt,self)
  end do
  ! calculate power average
  zelt%je=inlevel
  do i=self%decomp%firstpowe,self%decomp%lastpowe
     zelt%ie=i
     call powelt_avg(zelt,self)
  end do

end subroutine powcal_stat
!---------------------------------------------------------------------
!> write (vtk) powcal data
subroutine powcal_writev(self,kchar,kplot)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  character(*), intent(in) :: kchar !< case
  integer, intent(inout) :: kplot   !< output channel for vis. data

  !! local
  character(*), parameter :: s_name='powcal_writev' !< subroutine name
  type(posvecl_t), dimension(:), allocatable :: workpos !< save objlist variable
  integer(ki4) :: ing   !< save objlist variable
  integer(ki4) :: iwset   !< record whether weights present in geobjl
  integer(ki4) :: ifirst=1   !< unity for first call to vfile_rscalarwrite

  ing=self%powres%geobjl%ng
  iwset=self%powres%geobjl%nwset
  self%powres%geobjl%nwset=0
  plot_type: select case (kchar)
  case('psi-theta-zeta')
     ! rearrange nodl, keep copy of structure in work and ing
     ! ??reset number of objects to be output according to refine_level statistics input
     self%powres%geobjl%ng=ing*powelt_table(self%powres%n%nlevel,2)/powelt_table(infilelevel,2)
     self%powres%geobjl%ng=self%powres%npowe
     allocate(work(3*self%powres%geobjl%ng), stat=status)
     call log_alloc_check(m_name,s_name,11,status)
     work=self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)
     call geobjlist_nodlmv(self%powres%geobjl,infilelevel,1,self%powres%npowe)
     call geobjlist_writev(self%powres%geobjl,'geometry',kplot)
     ! restore nodl
     self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)=work
     deallocate(work)
     call vfile_rscalarwrite(self%powres%powa,self%powres%geobjl%ng,'Q-avg','CELL',kplot,ifirst)
     call vfile_rscalarwrite(self%powres%pows,self%powres%geobjl%ng,'Q-dev','CELL',kplot,0)

  case('allpsi-theta-zeta')
     ! reset number of objects to be output according to refine_level input
     self%powres%geobjl%ng=ing*powelt_table(self%powres%n%nlevel,2)/powelt_table(infilelevel,2)
     if (infilelevel>self%powres%n%nlevel) then
        allocate(work(3*self%powres%geobjl%ng), stat=status)
        call log_alloc_check(m_name,s_name,31,status)
        work=self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)
        call geobjlist_nodlmv(self%powres%geobjl,infilelevel,self%powres%n%nlevel,self%powres%npowe)
     end if
     call geobjlist_writev(self%powres%geobjl,'geometry',kplot)
     call vfile_rscalarwrite(abs(self%powres%pow),self%powres%geobjl%ng,'Q','CELL',kplot,ifirst)
     if (infilelevel>self%powres%n%nlevel) then
        ! restore nodl
        self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)=work
     end if

  case('cartesian')
     ! replace R-Z-xi positions with Cartesian
     allocate(workpos(self%powres%geobjl%np), stat=status)
     call log_alloc_check(m_name,s_name,21,status)
     workpos=self%powres%geobjl%posl%pos
     ! have to explicitly code as loop
     ! self%powres%geobjl%posl%pos=self%powres%vecx%pos
     ! to avoid ifort 12.0.0 segfault
     do l=1,self%powres%geobjl%np
        self%powres%geobjl%posl%pos(l)=self%powres%vecx%pos(l)
     end do
     ! rearrange nodl, keep copy of structure in work and ing
     ! reset number of objects to be output according to refine_level statistics input
     self%powres%geobjl%ng=ing*powelt_table(self%n%nlevel,2)/powelt_table(infilelevel,2)
     self%powres%geobjl%ng=self%powres%npowe
     allocate(work(3*self%powres%geobjl%ng), stat=status)
     call log_alloc_check(m_name,s_name,22,status)
     work=self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)
     call geobjlist_nodlmv(self%powres%geobjl,infilelevel,1,self%powres%npowe)
     call geobjlist_writev(self%powres%geobjl,'geometry',kplot)
     self%powres%geobjl%posl%pos=workpos
     deallocate(workpos)
     ! restore nodl
     self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)=work
     deallocate(work)
     call vfile_rscalarwrite(self%powres%powa,self%powres%geobjl%ng,'Q-avg','CELL',kplot,ifirst)
     call vfile_rscalarwrite(self%powres%pows,self%powres%geobjl%ng,'Q-dev','CELL',kplot,0)

  case('allcartesian')
     ! replace R-Z-xi positions with Cartesian
     allocate(workpos(self%powres%geobjl%np), stat=status)
     call log_alloc_check(m_name,s_name,40,status)
     !DBG write(*,*) "np=",self%powres%geobjl%np !DBG
     !DBG write(*,*) "vecx size=",size(self%powres%vecx%pos) !DBG
     !DBG write(*,*) "geobjl size=",size(self%powres%geobjl%posl%pos) !DBG
     workpos=self%powres%geobjl%posl%pos
     !DBG write(*,*) '01' !DBG
     ! have to explicitly code as loop
     ! self%powres%geobjl%posl%pos=self%powres%vecx%pos
     ! to avoid ifort 12.0.0 segfault
     do l=1,self%powres%geobjl%np
        self%powres%geobjl%posl%pos(l)=self%powres%vecx%pos(l)
     end do
     !DBG write(*,*) '02' !DBG
     ! reset number of objects to be output according to refine_level input
     self%powres%geobjl%ng=ing*powelt_table(self%n%nlevel,2)/powelt_table(infilelevel,2)
     !DBG write(*,*) '03' !DBG
     if (infilelevel>self%n%nlevel) then
        !DBG write(*,*) '04' !DBG
        !DBG write(*,*) infilelevel,self%n%nlevel !DBG
        allocate(work(3*self%powres%geobjl%ng), stat=status)
        call log_alloc_check(m_name,s_name,41,status)
        work=self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)
        call geobjlist_nodlmv(self%powres%geobjl,infilelevel,self%n%nlevel,self%powres%npowe)
     end if
     call geobjlist_writev(self%powres%geobjl,'geometry',kplot)
     !DBG write(*,*) '05' !DBG
     allocate(rwork(self%powres%geobjl%ng), stat=status)
     call log_alloc_check(m_name,s_name,42,status)
     rwork=abs(self%powres%pow(:self%powres%geobjl%ng))
     call vfile_rscalarwrite(rwork,self%powres%geobjl%ng,'Q','CELL',kplot,ifirst)
     deallocate(rwork)
     if (allocated(self%powres%psista)) then
        call vfile_rscalarwrite(self%powres%psista,self%powres%geobjl%ng,'psista','CELL',kplot,0)
     end if
     if (allocated(self%powres%angle)) then
        call vfile_rscalarwrite(self%powres%angle,self%powres%geobjl%ng,'angle','CELL',kplot,0)
     end if
     if (allocated(self%powres%lenpath)) then
        call vfile_rscalarwrite(self%powres%lenpath, self%powres%geobjl%ng, 'lenpath', 'CELL', kplot, 0)
     end if
     if (allocated(self%powres%zbdotn)) then
        call vfile_rscalarwrite(self%powres%zbdotn,self%powres%geobjl%ng,'zbdotn','CELL',kplot,0)
     end if
     if (allocated(self%powres%bfield)) then
        call vfile_dvectorwrite(self%powres%bfield,self%powres%geobjl%ng,'bfield','CELL',kplot,0)
     end if
     if (allocated(self%powres%norm)) then
        call vfile_dvectorwrite(self%powres%norm,self%powres%geobjl%ng,'normal','CELL',kplot,0)
     end if
     call vfile_iscalarwrite(self%powres%pmask,self%powres%geobjl%ng,'Msign','CELL',kplot,0)
     if (infilelevel>self%n%nlevel) then
        ! restore nodl
        self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)=work
     end if
     self%powres%geobjl%posl%pos=workpos
     deallocate(workpos)

  end select plot_type

  ! restore weight marker, number of objects
  self%powres%geobjl%nwset=iwset
  self%powres%geobjl%ng=ing

end subroutine powcal_writev
!---------------------------------------------------------------------
!> set refinement level
subroutine powcal_refine(self,gshadl,btree,rxpt,zxpt,psixpt)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure
  type(btree_t), intent(inout), optional :: btree !< btree data
  real(kr8) :: rxpt(2)   !<  position x point \f$ R \f$
  real(kr8) :: zxpt(2)   !<  position x point \f$ Z \f$
  real(kr8) :: psixpt(2) !<  \f$ \psi \f$ at x point
    
  !! local
  character(*), parameter :: s_name='powcal_refine' !< subroutine name
  type(powelt_t) :: post   !< powelt data structure
  
  calcn_type: select case (self%n%caltype)
  case('afws','local')
  case('msum','middle')
     ! need launch and shadow obj arrays  (=weight)s
     allocate(self%powres%geobjl%obj(self%powres%geobjl%ng), stat=status)
     self%powres%geobjl%nwset=2
     allocate(gshadl%obj(gshadl%ng), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     gshadl%nwset=2
  case('msus','global')
     ! need launch obj array (=weight)
     allocate(self%powres%geobjl%obj(self%powres%geobjl%ng), stat=status)
     call log_alloc_check(m_name,s_name,2,status)
     self%powres%geobjl%nwset=2
  end select calcn_type
  imlevel=1
  do level=1,self%n%nlevel
     inlevel=powelt_table(level,2)
     if (self%n%shadow/=0) then
        call powcal_move(self,post,gshadl,btree)
     end if
     call powcal_calc(self,post,gshadl,rxpt,zxpt,psixpt)
     imlevel=inlevel+1
  end do

  call powcal_stat(self,gshadl)

end subroutine powcal_refine
!---------------------------------------------------------------------
!> write powcal data
subroutine powcal_write(self,kout)

  !! arguments
  type(powcal_t), intent(in) :: self   !< powcal data structure
  integer, intent(in) :: kout   !< output channel for powcal data structure

  !! local
  character(*), parameter :: s_name='powcal_write' !< subroutine name

  write(kout,*,iostat=status) 'npowe'
  call log_write_check(m_name,s_name,1,status)
  write(kout,*,iostat=status) self%powres%npowe
  call log_write_check(m_name,s_name,2,status)
  write(kout,*,iostat=status) 'npow'
  call log_write_check(m_name,s_name,3,status)
  write(kout,*,iostat=status) self%powres%geobjl%ng
  call log_write_check(m_name,s_name,4,status)
  !X      write(kout,*,iostat=status) 'rblfac'
  !X      call log_write_check(m_name,s_name,5,status)
  !X      write(kout,*,iostat=status) self%powres%rblfac
  !X      call log_write_check(m_name,s_name,5,status)
  !X      write(kout,*,iostat=status) 'fpfac'
  !X      call log_write_check(m_name,s_name,6,status)
  !X      write(kout,*,iostat=status) self%powres%fpfac
  !X      call log_write_check(m_name,s_name,7,status)
  !
  write(kout,*,iostat=status) 'pows'
  call log_write_check(m_name,s_name,10,status)
  write(kout,*,iostat=status) self%powres%pows
  call log_write_check(m_name,s_name,11,status)
  write(kout,*,iostat=status) 'pow'
  call log_write_check(m_name,s_name,12,status)
  write(kout,*,iostat=status) self%powres%pow
  call log_write_check(m_name,s_name,13,status)

  !X      write(kout,*,iostat=status) 'f'
  !X      call log_write_check(m_name,s_name,20,status)
  !X      write(kout,*,iostat=status) self%n%f
  !X      call log_write_check(m_name,s_name,21,status)
  !X      write(kout,*,iostat=status) 'lmid'
  !X      call log_write_check(m_name,s_name,22,status)
  !X      write(kout,*,iostat=status) self%n%lmid
  !X      call log_write_check(m_name,s_name,23,status)
  !X      write(kout,*,iostat=status) 'ploss'
  !X      call log_write_check(m_name,s_name,24,status)
  !X      write(kout,*,iostat=status) self%n%ploss
  !X      call log_write_check(m_name,s_name,25,status)
  write(kout,*,iostat=status) 'nlevel'
  call log_write_check(m_name,s_name,26,status)
  write(kout,*,iostat=status) self%n%nlevel
  call log_write_check(m_name,s_name,27,status)
  write(kout,*,iostat=status) 'shadow'
  call log_write_check(m_name,s_name,28,status)
  write(kout,*,iostat=status) self%n%shadow
  call log_write_check(m_name,s_name,29,status)
  write(kout,*,iostat=status) 'start_control'
  call log_write_check(m_name,s_name,30,status)
  write(kout,*,iostat=status) self%odes%n%nstartcon
  call log_write_check(m_name,s_name,31,status)
  write(kout,*,iostat=status) 'termination_control'
  call log_write_check(m_name,s_name,32,status)
  write(kout,*,iostat=status) self%odes%n%ntermcon
  call log_write_check(m_name,s_name,33,status)
  write(kout,*,iostat=status) 'termination_parameters'
  call log_write_check(m_name,s_name,34,status)
  write(kout,*,iostat=status) self%odes%n%termcon
  call log_write_check(m_name,s_name,35,status)
  write(kout,*,iostat=status) 'analytic_launch_type'
  call log_write_check(m_name,s_name,36,status)
  write(kout,*,iostat=status) self%n%nanalau
  call log_write_check(m_name,s_name,37,status)
  write(kout,*,iostat=status) 'power_on_launch_geo'
  call log_write_check(m_name,s_name,38,status)
  write(kout,*,iostat=status) self%n%nlaupow
  call log_write_check(m_name,s_name,39,status)
  write(kout,*,iostat=status) 'power_on_shadow_geo'
  call log_write_check(m_name,s_name,40,status)
  write(kout,*,iostat=status) self%n%nshapow
  call log_write_check(m_name,s_name,41,status)
  !X      write(kout,*,iostat=status) 'sigma'
  !X      call log_write_check(m_name,s_name,42,status)
  !X      write(kout,*,iostat=status) self%n%sigma
  !X      call log_write_check(m_name,s_name,43,status)

  call edgprof_write(self%edgprof,kout)

  write(kout,*,iostat=status) 'pmask'
  call log_write_check(m_name,s_name,12,status)
  write(kout,*,iostat=status) self%powres%pmask
  call log_write_check(m_name,s_name,13,status)

end subroutine powcal_write
!---------------------------------------------------------------------
!> transfer power to shadow wall geometry
subroutine powcal_xferpow(self,gshadl)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  type(geobjlist_t), intent(inout) :: gshadl   !< shadow data structure

  !! local
  character(*), parameter :: s_name='powcal_xferpow' !< subroutine name
  integer(ki4) :: iobj   !< object number

  if (gshadl%nwset==0) then
     allocate(gshadl%obj(gshadl%ng), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     gshadl%nwset=2
  end if
  ! zero assigned power
  do j=1,gshadl%ng
     gshadl%obj(j)%weight=0._kr4
  end do
  ! assign power to saved object number
  do j=1,self%powres%geobjl%ng
     iobj=self%powres%geobjl%obj(j)%weight+0.5
     if (iobj>0) then
        gshadl%obj(iobj)%weight=self%powres%pow(j)
     end if
  end do

end subroutine powcal_xferpow
!---------------------------------------------------------------------
!> delete powcal_t
subroutine powcal_delete(self)

  !! arguments
  type(powcal_t), intent(inout) :: self !<  powcal data
  !! local

  call position_deleteveclis(self%powres%vecb)
  call position_deleteveclis(self%powres%vecx)
  call geobjlist_delete(self%powres%geobjl)
  deallocate(self%powres%pows)
  deallocate(self%powres%pow)
  deallocate(self%powres%pmask)
  calcn_type: select case (self%n%caltype)
  case('afws','local','msus','global')
     deallocate(self%powres%psista)
     deallocate(self%powres%angle)
     deallocate(self%powres%zbdotn)
     deallocate(self%powres%lenpath)
     deallocate(self%powres%bfield)
     deallocate(self%powres%norm)
  case default
  end select calcn_type

  field_type: select case (self%powres%beq%n%fldspec)
  case(1)
     call beq_deletepart(self%powres%beq)
  case default
     call beq_deleteplus(self%powres%beq)
  end select field_type
  call edgprof_delete(self%edgprof)

end subroutine powcal_delete

end module powcal_m
