module skyl_m

  use log_m
  use misc_m
  use dcontrol_h
  use const_numphys_h
  use const_kind_m
  use fmesh_h
  use spl2d_m
  use spl3d_m
  use geobjlist_h
  use geobjlist_m
  use vfile_m
  use skyl_h
  use beq_h
  use dcontrol_m

  implicit none
  private

! public subroutines
  public :: &
  skyl_initfile,  & !< open file
  skyl_init,  & !< initialise object
  skyl_readcon,  & !< read data from file
  skyl_fixupn, & !< fixup skyl numerics data structure
  skyl_fixup1, & !< fix up for skylight, replace missing data
  skyl_fixup2, & !< fix up for skylight, running 3-point extremum
  skyl_provis, &  !< produce one array for each skylight by compressing in/ouboxrz
  skyl_read, &  !< read in object
  skyl_write, &  !< write out object
  skyl_writeg, & !< write gnu skyl data structure
  skyl_delete, & !< delete object
  skyl_close !< close file

! private variables
  character(*), parameter :: m_name='skyl_m' !< module name
  integer  :: status   !< error status
  integer, save  :: ninso=5     !< control file unit number
  integer, save  :: noutso=6      !< output file unit number
  character(len=80), dimension(10), save :: dbgfile !< debug file names
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer  :: ilog      !< for namelist dump after error
  character(len=80) :: ibuff !< buffer for input/output
  character(len=256) :: vtkdesc !< descriptor line for vtk files
  real(kr8), dimension(:), allocatable :: work1 !< 1D work array
  real(kr8), dimension(:,:), allocatable :: work2 !< 2D work array

  contains
!---------------------------------------------------------------------
!> open file
subroutine skyl_initfile(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer, intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='skyl_initfile' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  if (trim(file)=='null') then
     call log_error(m_name,s_name,1,log_info,'null filename ignored')
     return
  end if

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then
  !! ninso=i if (present(channel)) channel=i exit end if end do

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  call misc_getfileunit(ninso)
  open(unit=ninso,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,2,error_fatal,'Cannot open control data file')
     stop
  end if
  if (present(channel)) channel=ninso

end subroutine skyl_initfile
!---------------------------------------------------------------------
!> initialise object
subroutine skyl_init(self,beq,geobjl,fileroot)

  !! arguments
  type(skyl_t), intent(inout) :: self !< module object
  type(beq_t), intent(in) :: beq !< beq object
  type(geobjlist_t), intent(in) :: geobjl !< geobjlist object
  character(*), intent(in) :: fileroot !< file root

  !! local
  character(*), parameter :: s_name='skyl_init' !< subroutine name
  integer(ki4) :: jpla !<  number of skylight planes to add to geobjlist
  real(kr8) :: rmin !< \f$ R_{\min} \f$ for equilibrium data
  real(kr8) :: rmax !< \f$ R_{\max} \f$ for equilibrium data
  real(kr8) :: zmin !< \f$ Z_{\min} \f$ for equilibrium data
  real(kr8) :: zmax !< \f$ Z_{\max} \f$ for equilibrium data
  real(kr4) :: totarea !< total area of geobjlist objects

  rmin=beq%rmin
  rmax=beq%rmax
  zmin=beq%zmin
  zmax=beq%zmax

  self%eps=min( const_epsbdry*(rmax-rmin), const_epsbdry*(zmax-zmin) )

  if ( .NOT.self%n%lrext.OR.(self%n%rextmax-self%n%rextmin)<const_epsbdry*(rmax-rmin) ) then
     self%n%rextmin=rmin
     self%n%rextmax=rmax
  else
     self%n%rextmin=max(self%n%rextmin,rmin)
     self%n%rextmax=min(self%n%rextmax,rmax)
  end if

  if ( .NOT.self%n%lzext.OR.(self%n%zextmax-self%n%zextmin)<const_epsbdry*(zmax-zmin) ) then
     self%n%zextmin=zmin
     self%n%zextmax=zmax
  else
     self%n%zextmin=max(self%n%zextmin,zmin)
     self%n%zextmax=min(self%n%zextmax,zmax)
  end if

  call geobjlist_totarea(geobjl,totarea)

  self%geoml=sqrt(totarea/(2*geobjl%ng))

  if (beq%n%skyldbg>0) then
     self%ndskyln=10
     dbgfile(1)=trim(fileroot)//"_centr_skyl.gnu"
     dbgfile(2)=trim(fileroot)//"_insil_skyl.gnu"
     dbgfile(3)=trim(fileroot)//"_ousil_skyl.gnu"
     dbgfile(4)=trim(fileroot)//"_inest_skyl.gnu"
     dbgfile(5)=trim(fileroot)//"_ouest_skyl.gnu"
     dbgfile(6)=trim(fileroot)//"_inpsi_skyl.gnu"
     dbgfile(7)=trim(fileroot)//"_oupsi_skyl.gnu"
     dbgfile(8)=trim(fileroot)//"_skyl"
     dbgfile(9)=trim(fileroot)//"_skyl1"
     dbgfile(10)=trim(fileroot)//"_skyl2"
     do j=1,7
        call misc_getfileunit(self%ndskyl(j))
        open(unit=self%ndskyl(j), file=dbgfile(j),form='FORMATTED',iostat=status)
        if(status/=0)then
           !! error opening file
           call log_error(m_name,s_name,j,error_fatal,'Error opening skyl debug file')
        else
           call log_error(m_name,s_name,j,log_info,'skyl debug file opened')
        end if
     end do
     do j=8,10
        call geobjlist_makehedline(geobjl,'skylight debug',vtkdesc)
        call vfile_init(dbgfile(j),vtkdesc,self%ndskyl(j))
     end do
  else
     self%ndskyln=0
     self%ndskyl=0
  end if

end subroutine skyl_init
!---------------------------------------------------------------------
!> read data from file
subroutine skyl_readcon(selfn,kin)

  !! arguments
  type(sknumerics_t), intent(out) :: selfn !< type which data will be assigned to
  integer, intent(in),optional :: kin   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='skyl_readcon' !< subroutine name
  logical :: skylight_lower !< use lower skylight
  integer(ki4), dimension(3) :: skylight_control !< skylight controls
  logical :: skylight_upper !< use upper skylight
  real(kr8), dimension(4) :: flux_extent !< extent in flux terms
  real(kr8), dimension(4) :: flux_bin_extent !< extent of bin in flux terms (INERT)
  integer(ki4), dimension(4) :: flux_extent_option !< choose how to calculate psiext (1 = user input, 2 = based on psi on geometry)
  integer(ki4), dimension(4) :: number_flux_boxes !< number of binning boxes in calculation of skylight(s)
  integer(ki4) :: number_of_divisions !< number of divisions in angle of skylight
  logical :: set_r_extent !< setting extent in \f$ R \f$
  logical :: set_z_extent !< setting extent in \f$ Z \f$
  real(kr8) :: flux_extent_margin !< margin for flux extent calculation of skylight(s)
  real(kr8) :: number_geoml_in_bin !< number of typical lengthscales in sample extent
  real(kr8) :: r_extent_min !< minimum value of \f$ R \f$ in extent calculation
  real(kr8) :: r_extent_max !< maximum value of \f$ R \f$ in extent calculation
  real(kr8) :: z_extent_min !< minimum value of \f$ Z \f$ in extent calculation
  real(kr8) :: z_extent_max !< maximum value of \f$ Z \f$ in extent calculation
  !> Distance that assumes intersection with skylight, in quantised units, in each direction
  !!  typically 1/1024 of geometry, default 1.
  real(kr8), dimension(3) :: intersect_tolerance !< .
  logical :: lpass !< input data is consistent
  logical :: lfail !< input data is not consistent

  !! skyl parameters
  namelist /skylparameters/ &
 &skylight_lower, skylight_upper, &
 &skylight_control, &
 &flux_extent_option, &
 &flux_bin_extent, &
 &flux_extent, &
 &set_r_extent, &
 &set_z_extent, &
 &flux_extent_margin, &
 &number_geoml_in_bin, &
 &r_extent_min, &
 &r_extent_max, &
 &z_extent_min, &
 &z_extent_max, &
 &intersect_tolerance, &
 &number_flux_boxes, &
 &number_of_divisions

  !! set default skyl parameters
  flux_extent=0._kr8
  flux_extent_option=2
  flux_bin_extent=0._kr8
  set_r_extent=.FALSE.
  set_z_extent=.FALSE.
  flux_extent_margin=0.1
  number_geoml_in_bin=10
  r_extent_min=0.1
  r_extent_max=0.1
  z_extent_min=0.
  z_extent_max=0.
  intersect_tolerance=1._kr8
  number_flux_boxes=20
  number_of_divisions=20

  skylight_lower=.FALSE.
  skylight_control=0
  skylight_upper=.FALSE.

  if(present(kin).AND.kin/=0) then
     !! assume unit already open and reading infile
     ninso=kin
  end if

  !!read skyl parameters
  read(ninso,nml=skylparameters,iostat=status)
  if(status/=0) then
     !!dump namelist contents to logfile to assist error location
     print '("Fatal error reading skyl parameters")'
     call log_getunit(ilog)
     write(ilog,nml=skylparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading skyl parameters')
  end if

  !! check for valid data
  do j=1,4
     if (flux_extent_option(j)<0.OR.flux_extent_option(j)>2) then
        call log_error(m_name,s_name,10,error_fatal,'invalid flux extent option')
     else
        if (flux_extent_option(j)==1) then
           if (skylight_lower) then
              lpass=(flux_extent(1)>0.AND.flux_extent(2)>0)
              lfail=(flux_bin_extent(1)<0.AND.number_flux_boxes(1)<1).OR. &
 &            (flux_bin_extent(2)<0.AND.number_flux_boxes(2)<1)
              if (.NOT.lpass.OR.lfail) call log_error(m_name,s_name,11,error_fatal,'invalid flux extent sampling lower')
           end if
           if (skylight_upper) then
              lpass=(flux_extent(3)>0.AND.flux_extent(4)>0)
              lfail=(flux_bin_extent(3)<0.AND.number_flux_boxes(3)<1).OR. &
 &            (flux_bin_extent(4)<0.AND.number_flux_boxes(4)<1)
              if (.NOT.lpass.OR.lfail) call log_error(m_name,s_name,12,error_fatal,'invalid flux extent sampling upper')
           end if
        end if
        if (flux_extent_option(j)==2) then
           if (skylight_lower) then
              lfail=(number_flux_boxes(1)<1.OR.number_flux_boxes(2)<1)
              if (lfail) call log_error(m_name,s_name,14,error_fatal,'invalid flux extent sampling lower')
           end if
           if (skylight_upper) then
              lfail=(number_flux_boxes(3)<1.OR.number_flux_boxes(4)<1)
              if (lfail) call log_error(m_name,s_name,15,error_fatal,'invalid flux extent sampling upper')
           end if
        end if
     end if
  end do

  if (any(skylight_control<0)) then
     call log_error(m_name,s_name,16,error_warning,'Only non-negative skylight controls expected')
  end if

  if (set_r_extent) then
     lfail=(r_extent_min-r_extent_max>=0.OR.r_extent_min<=0)
     if (lfail) call log_error(m_name,s_name,21,error_fatal,'invalid R extent sampling window')
  end if
  if (set_z_extent) then
     lfail=(z_extent_min-z_extent_max>=0)
     if (lfail) call log_error(m_name,s_name,22,error_fatal,'invalid Z extent sampling window')
  end if

  if (flux_extent_margin<0) then
     call log_error(m_name,s_name,23,error_fatal,'flux extent margin negative')
  end if
  if (number_geoml_in_bin<0) then
     call log_error(m_name,s_name,24,error_fatal,'number geoml in bin negative')
  end if

  !! store values
  selfn%skyltyp=0
  if (skylight_lower) then
     selfn%skyltyp=1
     if (skylight_upper) selfn%skyltyp=2
  else
     if (skylight_upper) selfn%skyltyp=3
  end if
  !! check consistency
  if (selfn%skyltyp==0) then
     call log_error(m_name,s_name,25,error_warning,'neither upper nor lower skylight specified')
  end if

  selfn%psiexts=flux_extent
  selfn%extsopt=flux_extent_option
  selfn%binsize=flux_bin_extent
  selfn%nexts=number_flux_boxes
  selfn%div=number_of_divisions
  selfn%lrext=set_r_extent
  selfn%extsdel=flux_extent_margin
  selfn%ngeoml=number_geoml_in_bin
  selfn%lzext=set_z_extent
  selfn%rextmin=r_extent_min
  selfn%rextmax=r_extent_max
  selfn%zextmin=z_extent_min
  selfn%zextmax=z_extent_max
  selfn%toli=intersect_tolerance
  selfn%control=skylight_control

end  subroutine skyl_readcon
!---------------------------------------------------------------------
!> fixup skyl numerics data structure
subroutine skyl_fixupn(selfn,klplot)

  !! arguments
  type(sknumerics_t), intent(inout) :: selfn !< object
  logical, dimension(2), intent(in) :: klplot   !< request for provisional skylight

  !! local
  character(*), parameter :: s_name='skyl_fixupn' !< subroutine name

  ! lower skylight
  if (klplot(1).AND..NOT.selfn%skyltyp<=2) then
     call log_error(m_name,s_name,1,error_warning,'Provisional lower skylight requested but lower_skylight not set')
  end if
  ! upper skylight
  if (klplot(2).AND..NOT.selfn%skyltyp>=2) then
     call log_error(m_name,s_name,2,error_warning,'Provisional upper skylight requested but upper_skylight not set')
  end if

  selfn%lprovis(1)=klplot(1).AND.selfn%skyltyp<=2
  selfn%lprovis(2)=klplot(2).AND.selfn%skyltyp>=2

end subroutine skyl_fixupn
!---------------------------------------------------------------------
!> fix up for skylight, replace missing data
subroutine skyl_fixup1(self,ktyps,kcall,kontrol)

  !! arguments
  type(skyl_t), intent(inout) :: self !< module object
  integer(ki4), intent(in) :: ktyps !<  lower (1) or upper (2) skylight type
  integer(ki4), intent(in) :: kcall !<  number of call
  integer(ki4), dimension(3), intent(in) :: kontrol !<  controls
  !! local
  character(*), parameter :: s_name='skyl_fixup' !< subroutine name
  integer(ki4) :: irid !< ! direction of travel in \f$ Z \f$ (+1) in lower, (-1) in upper
  integer(ki4) :: ibeg !< object array index

  ! direction of travel in Z (negative of idir used elsewhere)
  irid=-(2*ktyps-3)

  ! inner
  if (all(self%inboxr(:,kcall)>const_pushinf)) then
     call log_error(m_name,s_name,1,error_warning,'No inboard skylight')
     call log_value("call number ",kcall)
     self%inboxr(:,kcall)=self%n%rextmin
  else
     ibeg=1
     ! find first bin with OK contents
     do i=1,self%dimbox(1,kcall)
        if (self%inboxr(i,kcall)<const_pushinf) then
           ibeg=i
           exit
        end if
     end do
     ! fill all preceding bins
     do i=1,ibeg-1
        self%inboxr(i,kcall)=self%inboxr(ibeg,kcall)
        self%inboxz(i,kcall)=self%inboxz(ibeg,kcall)
     end do
     do i=ibeg+1,self%dimbox(1,kcall)
        if (self%inboxr(i,kcall)>const_pushinf) then
           self%inboxr(i,kcall)=self%inboxr(i-1,kcall)
           self%inboxz(i,kcall)=self%inboxz(i-1,kcall)
        end if
     end do
     ! make monotone in Z
     if (kontrol(2)>0) then
        do i=2,self%dimbox(1,kcall)
           if ((self%inboxz(i-1,kcall)-self%inboxz(i,kcall))*irid>0) then
              self%inboxz(i,kcall)=self%inboxz(i-1,kcall)
           end if
        end do
     end if
  end if

  ! outer
  if (all(self%ouboxr(:,kcall)>const_pushinf)) then
     call log_error(m_name,s_name,2,error_warning,'No outboard skylight')
     call log_value("call number ",kcall)
     self%ouboxr(:,kcall)=self%inboxr(1,kcall)
  else
     ibeg=1
     ! find first bin with OK contents
     do i=1,self%dimbox(2,kcall)
        if (self%ouboxr(i,kcall)<const_pushinf) then
           ibeg=i
           exit
        end if
     end do
     ! fill all preceding bins
     do i=1,ibeg-1
        self%ouboxr(i,kcall)=self%ouboxr(ibeg,kcall)
        self%ouboxz(i,kcall)=self%ouboxz(ibeg,kcall)
     end do
     do i=ibeg+1,self%dimbox(2,kcall)
        if (self%ouboxr(i,kcall)>const_pushinf) then
           self%ouboxr(i,kcall)=self%ouboxr(i-1,kcall)
           self%ouboxz(i,kcall)=self%ouboxz(i-1,kcall)
        end if
     end do
     ! make monotone in Z
     if (kontrol(3)>0) then
        do i=2,self%dimbox(2,kcall)
           if ((self%ouboxz(i-1,kcall)-self%ouboxz(i,kcall))*irid>0) then
              self%ouboxz(i,kcall)=self%ouboxz(i-1,kcall)
           end if
        end do
     end if
  end if


end subroutine skyl_fixup1
!---------------------------------------------------------------------
!> fix up for skylight, running 3-point extremum
subroutine skyl_fixup2(self,ktyps,kcall,kontrol)

  !! arguments
  type(skyl_t), intent(inout) :: self !< module object
  integer(ki4), intent(in) :: ktyps !<  lower (1) or upper (2) skylight type
  integer(ki4), intent(in) :: kcall !<  number of call
  integer(ki4), dimension(3), intent(in) :: kontrol !<  controls
  !! local
  character(*), parameter :: s_name='skyl_fixup2' !< subroutine name
  integer(ki4) :: irid !< ! direction of travel in \f$ Z \f$ (+1) in lower, (-1) in upper
  integer(ki4) :: ijm !< object array index
  integer(ki4) :: ijp !< object array index
  integer(ki4) :: iend !< object array index

  ! direction of travel in Z (negative of idir used elsewhere)
  irid=-(2*ktyps-3)

  ! allocate work1 array
  allocate(work1(self%dimbox(1,kcall)), stat=status)
  call log_alloc_check(m_name,s_name,1,status)

  ! fix up for upper (code for lower, then multiply by irid)
  ! thus find maximum z in lower, minimum in upper
  ! NB trust first (nearest centre) entry and last
  self%inboxz=irid*self%inboxz
  self%ouboxz=irid*self%ouboxz
  !! inner
  work1(1)=irid*self%inboxz(1,kcall)
  iend=self%dimbox(1,kcall)
  do i=2,iend-1
     ij=i
     ijm=ij-1
     ijp=min( ij+1,self%dimbox(1,kcall) )
     work1(ij)=irid*max( self%inboxz(ijm,kcall),self%inboxz(ij,kcall),self%inboxz(ijp,kcall) )
  end do
  work1(iend)=irid*self%inboxz(iend,kcall)
  self%inboxz(:,kcall)=irid*work1
  deallocate(work1)

  allocate(work1(self%dimbox(2,kcall)), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  !! outer
  work1(1)=irid*self%ouboxz(1,kcall)
  iend=self%dimbox(2,kcall)
  do i=2,iend-1
     ij=i
     ijm=ij-1
     ijp=min( ij+1,self%dimbox(2,kcall) )
     work1(ij)=irid*max( self%ouboxz(ijm,kcall),self%ouboxz(ij,kcall),self%ouboxz(ijp,kcall) )
  end do
  work1(iend)=irid*self%ouboxz(iend,kcall)
  self%ouboxz(:,kcall)=irid*work1
  deallocate(work1)

end subroutine skyl_fixup2
!---------------------------------------------------------------------
!> produce one array for each skylight by compressing in/oubox
subroutine skyl_provis(self,kcall)

  !! arguments
  type(skyl_t), intent(inout) :: self !< module object
  integer(ki4), intent(in) :: kcall !<  number of call
  !! local
  character(*), parameter :: s_name='skyl_provis' !< subroutine name
  integer(ki4) :: iw !< object array index
  integer(ki4) :: im !< object array previous index
  integer(ki4) :: iend !< object array extent
  real(kr8) :: delr !< difference in value of \f$ R \f$ of adjacent points
  real(kr8) :: delz !< difference in value of \f$ Z \f$ of adjacent points

  ! check needed
  if (.NOT.self%n%lprovis(kcall)) return

  ! allocate work2 array
  allocate(work2(self%dimbox(1,kcall)+self%dimbox(2,kcall),2), stat=status)
  call log_alloc_check(m_name,s_name,1,status)

  !! inner, in reverse order
  iw=1
  iend=self%dimbox(1,kcall)
  work2(iw,1)=self%inboxr(iend,kcall)
  work2(iw,2)=self%inboxz(iend,kcall)
  do i=iend,2,-1
     im=i-1
     delr=self%inboxr(i,kcall)-self%inboxr(im,kcall)
     delz=self%inboxz(i,kcall)-self%inboxz(im,kcall)
     if ( abs(delr)>self%eps .OR. abs(delz)>self%eps ) then
        iw=iw+1
        work2(iw,1)=self%inboxr(i,kcall)
        work2(iw,2)=self%inboxz(i,kcall)
     end if
  end do
  !! outer
  iw=iw+1
  work2(iw,1)=self%ouboxr(1,kcall)
  work2(iw,2)=self%ouboxz(1,kcall)
  iend=self%dimbox(2,kcall)
  do i=2,iend
     im=i-1
     delr=self%ouboxr(i,kcall)-self%ouboxr(im,kcall)
     delz=self%ouboxz(i,kcall)-self%ouboxz(im,kcall)
     if ( abs(delr)>self%eps .OR. abs(delz)>self%eps ) then
        iw=iw+1
        work2(iw,1)=self%ouboxr(i,kcall)
        work2(iw,2)=self%ouboxz(i,kcall)
     end if
  end do

  !! sort in increasing R order not a good idea
  ! call dsort(work2(1,1), work2(1,2), iw, 2)

  if (kcall==1) then
     ! allocate provisN array and assign
     allocate(self%provis1(2,iw), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     do l=1,2
        self%provis1(l,:)=work2(1:iw,l)
     end do
     self%nprovis1=iw
  else
     allocate(self%provis2(2,iw), stat=status)
     call log_alloc_check(m_name,s_name,11,status)
     do l=1,2
        self%provis2(l,:)=work2(1:iw,l)
     end do
     self%nprovis2=iw
  end if
  deallocate(work2)

end subroutine skyl_provis
!---------------------------------------------------------------------
!> read skyl data
subroutine skyl_read(self, infile)

  !! arguments
  type(skyl_t), intent(out) :: self   !< skyl data structure
  character(*), intent(in) :: infile !< name of input file

  !! local
  character(*), parameter :: s_name='skyl_read' !< subroutine name
  integer :: iread   !< output channel for skyl data structure
  integer(ki4), dimension(2) :: indbox   !<  extents of box arrays
  integer(ki4) :: ibdim !<  second dimension of boxr,z arrays

  !! sort out unit
  inquire(file=infile,number=iread,iostat=status)
  call log_read_check(m_name,s_name,90,status)
  if (iread<0) then
     call log_error(m_name,s_name,91,error_fatal,'Error reading skyl data structure file')
  end if

  read(iread,*,iostat=status) ibuff
  !write(*,*) 'dbg 1', ibuff
  call log_read_check(m_name,s_name,1,status)
  read(iread,*,iostat=status) self%n%skyltyp
  call log_read_check(m_name,s_name,2,status)
  if (self%n%skyltyp<=0) return
  ibdim=2-mod(self%n%skyltyp,2)

  read(iread,*,iostat=status) ibuff
  !write(*,*) ibdim,'dbg 2', ibuff
  call log_read_check(m_name,s_name,10,status)
  read(iread,*,iostat=status) ((self%dimbox(i,j),i=1,2),j=1,ibdim)
  call log_read_check(m_name,s_name,11,status)
  ! dimension and allocate boxrz arrays
  if (ibdim==2) then
     indbox(1)=max( self%dimbox(1,1),self%dimbox(1,2) )
     indbox(2)=max( self%dimbox(2,1),self%dimbox(2,2) )
  else
     indbox=self%dimbox(:,1)
  end if
  allocate(self%inboxr(indbox(1),ibdim),self%inboxz(indbox(1),ibdim), &
  self%ouboxr(indbox(2),ibdim),self%ouboxz(indbox(2),ibdim),&
  stat=status)

  call log_alloc_check(m_name,s_name,12,status)
  read(iread,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,13,status)
  read(iread,*,iostat=status) ((self%inboxr(i,j),i=1,indbox(1)),j=1,ibdim)
  call log_read_check(m_name,s_name,14,status)
  read(iread,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,15,status)
  read(iread,*,iostat=status) ((self%inboxz(i,j),i=1,indbox(1)),j=1,ibdim)
  call log_read_check(m_name,s_name,16,status)
  read(iread,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,17,status)
  read(iread,*,iostat=status) ((self%ouboxr(i,j),i=1,indbox(2)),j=1,ibdim)
  call log_read_check(m_name,s_name,18,status)
  read(iread,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,19,status)
  read(iread,*,iostat=status) ((self%ouboxz(i,j),i=1,indbox(2)),j=1,ibdim)
  call log_read_check(m_name,s_name,20,status)
  read(iread,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,30,status)
  read(iread,*,iostat=status) (((self%psilts(i,j,k),i=1,2),j=1,2),k=1,ibdim)
  call log_read_check(m_name,s_name,31,status)
  read(iread,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,32,status)
  read(iread,*,iostat=status) ((self%psidelta(i,j),i=1,2),j=1,ibdim)
  call log_read_check(m_name,s_name,33,status)
  read(iread,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,34,status)
  read(iread,*,iostat=status) (self%n%toli(j),j=1,3)
  call log_read_check(m_name,s_name,35,status)
  read(iread,*,iostat=status) ibuff
  call log_read_check(m_name,s_name,36,status)
  read(iread,*,iostat=status) self%geoml
  call log_read_check(m_name,s_name,37,status)

end subroutine skyl_read
!---------------------------------------------------------------------
!> write skyl data
subroutine skyl_write(self,kout)

  !! arguments
  type(skyl_t), intent(in) :: self   !< skyl data structure
  integer, intent(in), optional :: kout   !< output channel for skyl data structure

  !! local
  character(*), parameter :: s_name='skyl_write' !< subroutine name
  integer :: iout   !< output channel for skyl data structure
  integer(ki4), dimension(2) :: indbox   !<  extents of box arrays
  integer(ki4) :: ibdim !<  second dimension of boxr,z arrays

  !! sort out unit
  if(present(kout)) then
     iout=kout
  else
     iout=noutso
  end if

  write(iout,*,iostat=status) 'skyltyp'
  call log_write_check(m_name,s_name,1,status)
  write(iout,*,iostat=status) self%n%skyltyp
  call log_write_check(m_name,s_name,2,status)
  if (self%n%skyltyp<=0) return
  ibdim=2-mod(self%n%skyltyp,2)
  write(iout,*,iostat=status) 'dimbox'
  call log_write_check(m_name,s_name,11,status)
  write(iout,*,iostat=status) ((self%dimbox(i,j),i=1,2),j=1,ibdim)
  call log_write_check(m_name,s_name,12,status)
  if (ibdim==2) then
     indbox(1)=max( self%dimbox(1,1),self%dimbox(1,2) )
     indbox(2)=max( self%dimbox(2,1),self%dimbox(2,2) )
  else
     indbox=self%dimbox(:,1)
  end if
  write(iout,*,iostat=status) 'inboxr'
  call log_write_check(m_name,s_name,13,status)
  write(iout,*,iostat=status) ((self%inboxr(i,j),i=1,indbox(1)),j=1,ibdim)
  call log_write_check(m_name,s_name,14,status)
  write(iout,*,iostat=status) 'inboxz'
  call log_write_check(m_name,s_name,15,status)
  write(iout,*,iostat=status) ((self%inboxz(i,j),i=1,indbox(1)),j=1,ibdim)
  call log_write_check(m_name,s_name,16,status)
  write(iout,*,iostat=status) 'ouboxr'
  call log_write_check(m_name,s_name,17,status)
  write(iout,*,iostat=status) ((self%ouboxr(i,j),i=1,indbox(2)),j=1,ibdim)
  call log_write_check(m_name,s_name,18,status)
  write(iout,*,iostat=status) 'ouboxz'
  call log_write_check(m_name,s_name,19,status)
  write(iout,*,iostat=status) ((self%ouboxz(i,j),i=1,indbox(2)),j=1,ibdim)
  call log_write_check(m_name,s_name,20,status)
  write(iout,*,iostat=status) 'psilts'
  call log_write_check(m_name,s_name,30,status)
  write(iout,*,iostat=status) (((self%psilts(i,j,k),i=1,2),j=1,2),k=1,ibdim)
  call log_write_check(m_name,s_name,31,status)
  write(iout,*,iostat=status) 'psidelta'
  call log_write_check(m_name,s_name,32,status)
  write(iout,*,iostat=status) ((self%psidelta(i,j),i=1,2),j=1,ibdim)
  call log_write_check(m_name,s_name,33,status)
  write(iout,*,iostat=status) 'toli'
  call log_write_check(m_name,s_name,34,status)
  write(iout,*,iostat=status) (self%n%toli(j),j=1,3)
  call log_write_check(m_name,s_name,35,status)
  write(iout,*,iostat=status) 'geoml'
  call log_write_check(m_name,s_name,36,status)
  write(iout,*,iostat=status) self%geoml
  call log_write_check(m_name,s_name,37,status)

end subroutine skyl_write
!---------------------------------------------------------------------
!> write gnu skyl data structure
subroutine skyl_writeg(self,kchar,kout)

  !! arguments
  type(skyl_t), intent(inout) :: self !< object
  character(*), intent(in) :: kchar  !< case
  integer, intent(in) :: kout   !< output channel for gnuplot data

  !! local
  character(*), parameter :: s_name='skyl_writeg' !< subroutine name

  plot_type: select case (kchar)
  case('lower')
     ! positions in R-Z space
     do j=1,self%nprovis1
        write(kout,'(1x,'//cfmt2v,iostat=status) &
 &      1000*self%provis1(1,j),1000*self%provis1(2,j)
        call log_write_check(m_name,s_name,1,status)
     end do
  case('upper')
     ! positions in R-Z space
     do j=1,self%nprovis2
        write(kout,'(1x,'//cfmt2v,iostat=status) &
 &      1000*self%provis2(1,j),1000*self%provis2(2,j)
        call log_write_check(m_name,s_name,2,status)
     end do
  end select plot_type

end subroutine skyl_writeg
!---------------------------------------------------------------------
!> delete object
subroutine skyl_delete(self,ndebug)

  !! arguments
  type(skyl_t), intent(inout) :: self !< module object
  integer, intent(in), optional :: ndebug   !< debug channels to close

  !! local
  character(*), parameter :: s_name='skyl_delete' !< subroutine name
  integer :: idebug   !< debug channels to close

  idebug=0
  if (present(ndebug)) idebug=ndebug

  if (idebug==0) then
     deallocate(self%inboxr,self%inboxz)
     deallocate(self%ouboxr,self%ouboxz)
  else
     ! close units
     do j=1,self%ndskyln
        close(unit=self%ndskyl(j),iostat=status)
        if(status/=0)then
           !! error closing file
           print '("Warning message: Unable to close control file, ",a)',dbgfile(j)
           call log_error(m_name,s_name,1,error_warning,'Cannot close control data file')
        end if
     end do
  end if

end subroutine skyl_delete
!---------------------------------------------------------------------
!> close file
subroutine skyl_close

  !! local
  character(*), parameter :: s_name='skyl_close' !< subroutine name

  !! close file
  close(unit=ninso,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Warning message: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_warning,'Cannot close control data file')
  end if

end subroutine skyl_close

end module skyl_m
