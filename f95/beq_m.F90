module beq_m

  use const_kind_m
  use const_numphys_h
  use position_h
  use fmesh_h
  use beq_h
  use beqan_h
  use posang_h
  use position_m
  use log_m
  use misc_m
  use geobj_m
  use spl2d_m
  use spl3d_m
  use fmesh_m
  use posang_m
  use beqan_m
  use gfile_m

  implicit none
  private

! public subroutines
  public :: &
  beq_read, & !< read in data format
  beq_readequil, & !< read in non-strict EQDSK format
  beq_readequ, & !< read in EQU format
  beq_readana, & !< generate using ANA format
  beq_fixorigin,   & !< fix radial origin of  beq data structure
  beq_move, & !< move  beq data structure (controls)
  beq_init,   & !< create  beq data structure
  beq_sense,   & !< helical sense of field
  beq_delete,   & !< delete  beq data structure
  beq_readv, & !< read in visualisation format (TO DO)
  beq_readcheck, & !< check field as mapped or 3-cpt and for extras
  beq_readpart, & !< read field data needed by ITER
  beq_readplus, & !< read field data needed by 'msus','global'
  beq_writeg, & !< write in gnuplot format
  beq_writen, & !< write in nucode format
  beq_writev, & !< write in visualisation format
  beq_write, & !< write in data structure format
  beq_writepart, & !< write field data needed by ITER
  beq_deletepart, & !< delete field data used by ITER
  beq_writeplus, & !< write field data needed by 'msus','global'
  beq_deleteplus, & !< delete field data used by 'msus','global'
  beq_dia, &  !< calculate diagnostics and output
  beq_spl2d, & !< initialise spline data structure from beq
  beq_readcon, &  !< namelist input of controls
  beq_centre, &  !< calculate \f$ R_c, Z_c\f$
  beq_psilt, &  !< calculate \f$ \psi \f$ limits
  beq_thetalt, &  !< calculate \f$ \theta \f$ limits
  beq_xilt, &  !< calculate \f$ \xi \f$ limits
  beq_nset, &  !< calculate \f$ N_{\theta} \f$ and \f$ N_{\psi} \f$
  beq_bdryrb, &  !< calculate \f$ R_m \f$ and \f$ B_m \f$ from \f$ \psi_m \f$
  beq_psix, &  !< calculate  \f$ \psi \f$ at separatrix
  beq_rextrema, &  !< calculate \f$ r_{min} \f$ and \f$ r_{max} \f$ as functions of \f$ \theta_j \f$
  beq_rpsi, &  !< calculate \f$ r \f$ as a function of of \f$ \psi \f$
  beq_rispld, &  !< calculate spline coefficients for \f$ R/I \times \psi \f$ derivatives
  beq_rjpsi, &  !< calculate spline coefficients for \f$ R/J (\psi, \theta) \f$
  beq_ipsi, &   !< calculate spline coefficients for \f$ I(\psi) \f$ aka \f$ f(\psi) \f$
  beq_b, &   !< calculate B-field components at position given by posang (both cartesian)
  beq_ripple, &   !< calculate B-field ripple
  beq_ripple_calch1, & !> calculate B-field ripple as for MAST
  beq_rsig, &  !< return rsig to external module
  beq_rsigset, & !< set rsig, positive if psi increasing from centre
  beq_psicont, & !< calculate contour  \f$ \psi = \psi_{cont} \f$ as array of  \f$ R,Z \f$ values
  beq_findrzm, & !< calculate displacements (moves) based on argument and psibdry
  beq_ctrax, & !< calculate complete track of psi extremum through plasma centre
  beq_ctrack1r, & !< return track of psi extremum through plasma centre
  beq_ctrackc, & !< calculate track of psi extremum through plasma centre
  beq_ctrackcq, & !< query position relative to track of psi extremum through plasma centre
  beq_ctrackpt, & !< find \f$ (R, Z) \f$ position of psi on track through plasma centre
  beq_getunit,  & !< get unit number
  beq_ripple_h1 !< define ripple as for MAST

  private :: &
  beq_rextremum, &  !< calculate  extremal \f$ \psi \f$ in search direction
  beq_ctrack1 !< update of psi extremum through plasma centre


! public types and variables
!see beq_h file

! public variables

! private types

! private variables
  character(*), parameter :: m_name='beq_m' !< module name
  real(kr8), dimension(:), allocatable :: work1 !< 1D work array
  real(kr8), dimension(:), allocatable :: workr1 !< 1D work array
  real(kr8), dimension(:), allocatable :: workz1 !< 1D work array
  real(kr8), dimension(:), allocatable :: wvextn !< 1D work array extended size for nodes
  real(kr8), dimension(:), allocatable :: wvext !< 1D work array extended size for samples
  real(kr8), dimension(:), allocatable :: wvextd !< 1D work array extended size for derivative
  real(kr8), dimension(:,:), allocatable :: zwork2 !< scratch 2D work array
  real(kr8), dimension(:,:), allocatable :: work2 !< 2D work array
  real(kr8), dimension(:,:), allocatable :: workr2 !< 2D work array
  real(kr8), dimension(:,:), allocatable :: workz2 !< 2D work array
  real(kr8), dimension(:,:), allocatable :: worka !< 2D work array
  real(kr8), parameter :: epsg=1.e-6 !< tolerance for searches
  real(kr8), parameter :: epsr=0.01 !< relative tolerance for searches
  integer   :: status   !< error status
  integer, save :: nin=-1   !< input channel for beq data
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: idum !< dummy integer
  integer :: idunit !< debug unit
  real(kr8) :: zdum !< dummy real
  real(kr8), save :: rsig    !< sign of  \f$ d \psi/dr \f$, value  \f$ \pm 1 \f$
  character(len=80) :: ibuff !< buffer for input/output
  character(len=80) :: ibuf1 !< buffer for input/output
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> read in data format
subroutine beq_read(self,infile)

  !! arguments
  type(beq_t), intent(out) :: self   !< object data structure
  character(*),intent(in) :: infile !< name of input file

  !! local
  character(*), parameter :: s_name='beq_read' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

  !! open file
  call misc_getfileunit(nin)
  open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,1,error_fatal,'Error opening beq data structure file')
  else
     call log_error(m_name,s_name,2,log_info,'beq data structure file opened')
  end if

  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%f
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%rmin
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%rmax
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%mr
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%dr
  if(status/=0) then
     call log_error(m_name,s_name,5,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%zmin
  if(status/=0) then
     call log_error(m_name,s_name,6,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%zmax
  if(status/=0) then
     call log_error(m_name,s_name,7,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%mz
  if(status/=0) then
     call log_error(m_name,s_name,8,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%dz
  if(status/=0) then
     call log_error(m_name,s_name,9,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%i
  if(status/=0) then
     call log_error(m_name,s_name,10,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%psi,infile,nin )
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%dpsidr,infile,nin )
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%dpsidz,infile,nin )
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%r,infile,nin )
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%z,infile,nin )
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%rjac,infile,nin )
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%rcen
  if(status/=0) then
     call log_error(m_name,s_name,17,error_fatal,'Error reading spline data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%zcen
  if(status/=0) then
     call log_error(m_name,s_name,18,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%psimin
  if(status/=0) then
     call log_error(m_name,s_name,19,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%psimax
  if(status/=0) then
     call log_error(m_name,s_name,20,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%npsi
  if(status/=0) then
     call log_error(m_name,s_name,21,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%dpsi
  if(status/=0) then
     call log_error(m_name,s_name,22,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%thetamin
  if(status/=0) then
     call log_error(m_name,s_name,23,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%thetamax
  if(status/=0) then
     call log_error(m_name,s_name,24,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%ntheta
  if(status/=0) then
     call log_error(m_name,s_name,25,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%ntmin
  if(status/=0) then
     call log_error(m_name,s_name,26,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%ntmax
  if(status/=0) then
     call log_error(m_name,s_name,27,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%cenopt
  if(status/=0) then
     call log_error(m_name,s_name,28,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%psiopt
  if(status/=0) then
     call log_error(m_name,s_name,29,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%bdryopt
  if(status/=0) then
     call log_error(m_name,s_name,39,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%nopt
  if(status/=0) then
     call log_error(m_name,s_name,30,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%thetaopt
  if(status/=0) then
     call log_error(m_name,s_name,31,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%dtheta
  if(status/=0) then
     call log_error(m_name,s_name,32,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%psiref
  if(status/=0) then
     call log_error(m_name,s_name,33,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%thetaref
  if(status/=0) then
     call log_error(m_name,s_name,34,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%psibdry
  if(status/=0) then
     call log_error(m_name,s_name,35,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%psiqbdry
  if(status/=0) then
     call log_error(m_name,s_name,36,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%rbdry
  if(status/=0) then
     call log_error(m_name,s_name,37,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%bpbdry
  if(status/=0) then
     call log_error(m_name,s_name,38,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  if (adjustl(ibuff)=='btotbdry') then
     read(nin,*,iostat=status) self%btotbdry
     if(status/=0) then
        call log_error(m_name,s_name,39,error_fatal,'Error reading object data')
     end if
     read(nin,*,iostat=status) ibuff
  end if
  read(nin,*,iostat=status) self%srmin
  if(status/=0) then
     call log_error(m_name,s_name,40,error_fatal,'Error reading object data')
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%srmax
  if(status/=0) then
     call log_error(m_name,s_name,41,error_fatal,'Error reading object data')
  end if

end subroutine beq_read
!---------------------------------------------------------------------
!> read in non-strict EQDSK format
subroutine beq_readequil(self,infile,numerics)

  !! arguments
  type(beq_t), intent(out) :: self   !< object data structure
  character(*),intent(in) :: infile !< name of input file
  type(bnumerics_t), intent(in) :: numerics   !< numerics data structure

  !! local
  character(*), parameter :: s_name='beq_readequil' !< subroutine name
  character(len=8), parameter :: cfmt1='(5e16.9)'
  character(len=15) :: cfmtdh !< fixed format for header
  character(len=15) :: cfmtd !< fixed format for header
  logical, parameter :: debug=.TRUE. !< flag for e16.9 output
  !! logical :: unitused !< flag to test unit is available
  logical :: needfixup=.FALSE. !< flag special SMITER fix up
  integer :: istatus   !< inner status variable
  integer :: iin   !< input channel for object data structure
  integer(ki4) :: nw  !< EFIT number of grid points in R-direction
  integer(ki4) :: nh  !< EFIT number of grid points in Z-direction
  integer(ki4) :: ilen   !< length of string in buffer
  integer(ki4) :: irlen   !< reduced length of string in buffer
  integer(ki4) :: ifound   !< number found in line
  integer(ki4) :: inlin   !< number of valid lines found
  integer(ki4) :: ino  !< number of numbers on first line
  integer(ki4), dimension(2) :: ks=0 !< first index of number
  integer(ki4), dimension(2) :: js=0 !< second index of number
  integer(ki4) :: isfixedh  !< flag fixed-format input for header
  integer(ki4) :: isfixed  !< flag fixed-format input
  integer(ki4) :: iffiesta !< mistakes in fiesta eqdsk format
  real(kr8) :: bcentr !< EFIT vacuum toroidal field at r=rcentr ignored
  real(kr8) :: cpasma !< EFIT computed plasma current in A ignored
  real(kr8) :: rgrid !< EFIT Minimum R in metre of rectangular computational box
  real(kr8) :: rmaxis !< EFIT plasma centre (magnetic axis) in metre
  real(kr8) :: rzero !< EFIT variable ignored
  real(kr8) :: ssibry1 !< EFIT poloidal flux at the boundary (Wb/rad)
  real(kr8) :: ssibry2 !< EFIT poloidal flux at the boundary (Wb/rad) ignored
  real(kr8) :: ssimag1 !< EFIT poloidal flux at the magnetic axis (Wb/rad)
  real(kr8) :: ssimag2 !< EFIT poloidal flux at the magnetic axis (Wb/rad) ignored
  real(kr8) :: xdim !< EFIT Horizontal dimension in metre of computational box
  real(kr8) :: xdum !< EFIT variable ignored
  real(kr8) :: zdim !< EFIT Vertical dimension in metre of computational box
  real(kr8) :: zmaxis !< EFIT plasma centre (magnetic axis) in metre
  real(kr8) :: zmid !< EFIT centre of computational box in metre originally ignored
  real(kr8) :: zpsic !< diagnostic variable
  real(kr8) :: btcen !< diagnostic variable
  real(kr8) :: psifac !< scales big psi to usual psi in SMARDDA context
  integer(ki4) :: nbbbs  !< EFIT number of points describing boundary
  integer(ki4) :: limitr  !< EFIT number of points describing limitr
  integer(ki4) :: ncoil  !< EFIT number of points describing one of five coils
  integer :: igunit  !< unit for output of boundary data from eqdsk

  iffiesta=numerics%fiesta
  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then iin=i exit end if end do

  call misc_getfileunit(iin)
  open(unit=iin, file=infile,status='OLD',form='FORMATTED',iostat=status)
  call log_open_check(m_name,s_name,1,status)

  ! first line
  read(iin,'(a80)',iostat=status) ibuff
  if (status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error reading header data')
  end if
  ilen=len_trim(ibuff)

  irlen=ilen
  ! find last numeric
  ino=1
  do

     do j=irlen,1,-1
        js(ino)=j
        if ( lge(ibuff(j:j),'0').AND.lle(ibuff(j:j),'9') ) exit
     end do

     if (js(ino)==1) call log_error(m_name,s_name,3,error_fatal,'Not enough numbers in header data')

     ifound=0
     do k=js(ino),1,-1
        ks(ino)=k
        if ( lge(ibuff(k:k),'0').AND.lle(ibuff(k:k),'9') ) then
           cycle
        else if (ibuff(k:k)==' ') then
           ifound=1
           exit
        else
           exit
        end if
     end do
     if (ifound==0) call log_error(m_name,s_name,4,error_fatal,'Not enough numbers in header data')

     irlen=ks(ino)
     if (ino==2) exit

     ino=ino+1

  end do

  read( ibuff(ks(1):js(1)),* ) nh
  read( ibuff(ks(2):js(2)),* ) nw

  ! examine logical second and fifth lines

  inlin=0
  do
     read(iin,'(a80)',iostat=status)ibuff
     !!eof or error
     if (status/=0) then
        call log_error(m_name,s_name,5,error_fatal,'Error reading line of data after first')
     else if (ibuff(1:4)=='XXXX') then
        ! check for SMITER special fix-up
        needfixup=.TRUE.
        cycle
     else
        ! valid line
        inlin=inlin+1
        if (inlin==1) then
           call buffer_fixedfmt(ibuff,isfixedh,cfmtdh)
           call log_value("eqdsk header format ",cfmtdh)
        else if (inlin>=5) then
           call buffer_fixedfmt(ibuff,isfixed,cfmtd)
           call log_value("eqdsk body format ",cfmtd)
           exit
        end if
     end if
  end do
  ! now start again knowing the format
  rewind(iin)
  ! skip first line
  read(iin,'(a80)')ibuff
  if (needfixup) read(iin,'(a80)')ibuff

  ! EFIT G EQDSK defined at https://fusion.gat.com/theory/Efitgeqdsk
  if (isfixedh==0) then
     read(iin,*,iostat=istatus)xdim,zdim,rzero,rgrid,zmid
     if (debug) write(*,*) 'Header lines from eqdsk file---'
     if (debug) write(*,cfmt1)xdim,zdim,rzero,rgrid,zmid
     if(istatus/=0) then
        call log_error(m_name,s_name,10,error_fatal,'Error reading xdim etc.')
     end if
     read(iin,*,iostat=istatus)rmaxis,zmaxis,ssimag1,ssibry1,bcentr
     if (debug) write(*,cfmt1)rmaxis,zmaxis,ssimag1,ssibry1,bcentr
     if(istatus/=0) then
        call log_error(m_name,s_name,11,error_fatal,'Error reading rmaxis etc.')
     end if
     if (iffiesta>0) then
        ! ssimag1,ssibry1 wrong way round
        zdum=ssimag1
        ssimag1=ssibry1
        ssibry1=zdum
        call log_error(m_name,s_name,15,error_warning,'ssimag1,ssibry1 swapped')
     end if
     ! note use of multiple xdum's means some information is lost
     read(iin,*,iostat=istatus)cpasma,ssimag2,xdum,rmaxis,xdum
     if(istatus/=0) then
        call log_error(m_name,s_name,12,error_fatal,'Error reading cpasma etc.')
     end if
     if (debug) write(*,cfmt1)cpasma,ssimag2,xdum,rmaxis,xdum
     read(iin,*,iostat=istatus)zmaxis,xdum,ssibry2,xdum,xdum
     if(istatus/=0) then
        call log_error(m_name,s_name,13,error_fatal,'Error reading zmaxis etc.')
     end if
     if (debug) write(*,cfmt1)zmaxis,xdum,ssibry2,xdum,xdum
  else
     read(iin,cfmtdh,iostat=istatus)xdim,zdim,rzero,rgrid,zmid
     if(istatus/=0) then
        call log_error(m_name,s_name,50,error_fatal,'Error reading xdim etc.')
     end if
     if (debug) write(*,*) 'Header lines from eqdsk file---'
     if (debug) write(*,cfmt1)xdim,zdim,rzero,rgrid,zmid
     read(iin,cfmtdh,iostat=istatus)rmaxis,zmaxis,ssimag1,ssibry1,bcentr
     if(istatus/=0) then
        call log_error(m_name,s_name,51,error_fatal,'Error reading rmaxis etc.')
     end if
     if (debug) write(*,cfmt1)rmaxis,zmaxis,ssimag1,ssibry1,bcentr
     if (iffiesta>0) then
        ! ssimag1,ssibry1 wrong way round
        zdum=ssimag1
        ssimag1=ssibry1
        ssibry1=zdum
        call log_error(m_name,s_name,15,error_warning,'ssimag1,ssibry1 swapped')
     end if
     read(iin,cfmtdh,iostat=istatus)cpasma,ssimag2,xdum,rmaxis,xdum
     if(istatus/=0) then
        call log_error(m_name,s_name,52,error_fatal,'Error reading cpasma etc.')
     end if
     if (debug) write(*,cfmt1)cpasma,ssimag2,xdum,rmaxis,xdum
     read(iin,cfmtdh,iostat=istatus)zmaxis,xdum,ssibry2,xdum,xdum
     if(istatus/=0) then
        call log_error(m_name,s_name,53,error_fatal,'Error reading zmaxis etc.')
     end if
     if (debug) write(*,cfmt1)zmaxis,xdum,ssibry2,xdum,xdum
  end if

  call log_error(m_name,s_name,90,log_info,'Geometry parameters from EFIT')
  call log_value('Domain size in R xdim',xdim)
  call log_value('Domain size in Z zdim',zdim)
  call log_value('Domain start in R rgrid',rgrid)
  call log_value('Domain centre in Z (originally ignored) zmid',zmid)
  !call log_error(m_name,s_name,91,log_info,'Centre of domain in Z is assumed to be Z=0')
  call log_error(m_name,s_name,92,log_info,'Plasma parameters from EFIT')
  call log_value('Flux at centre ssimag1',ssimag1)
  call log_value('Flux at boundary ssibry1',ssibry1)
  call log_value('Plasma centre in R rmaxis',rmaxis)
  call log_value('Plasma centre in Z zmaxis',zmaxis)
  call log_value('Plasma current ',cpasma)
  call log_error(m_name,s_name,93,log_info,'Check consistency with other plasma data')
  if (debug) write(*,*) 'Geometry parameters from EFIT'
  if (debug) write(*,*) 'Domain size xdim=',xdim ,'zdim=',zdim
  if (debug) write(*,*) 'Domain start rgrid=',rgrid, 'Domain middle zmid=',zmid
  !if (debug) write(*,*) 'Centre of domain in Z is assumed to be Z=0'
  if (debug) write(*,*) 'Plasma parameters from EFIT'
  if (debug) write(*,*) 'Fluxes at centre ssimag1=',ssimag1, 'boundary ssibry1=',ssibry1
  if (debug) write(*,*) 'Plasma centre rmaxis=',rmaxis, 'zmaxis=',zmaxis
  if (debug) write(*,*) 'Check consistency with other plasma data'

  ! fpol array
  !! allocate fpol storage and read
  if(nw>0) then
     allocate(self%f(nw), stat=status)
     call log_alloc_check(m_name,s_name,15,status)
  else
     call log_error(m_name,s_name,16,error_fatal,'No 1D data')
  end if

  if (cfmtd=='*') then
     read(iin,*,iostat=status)(self%f(i),i=1,nw)
  else
     read(iin,cfmtd,iostat=status)(self%f(i),i=1,nw)
  end if
  if(status/=0) then
     call log_error(m_name,s_name,17,error_fatal,'Error reading fpol')
  end if
  print '("number of fpol values read = ",i10)',nw
  call log_value("number of fpol values read ",nw)
  btcen=abs(self%f(1))/max(rmaxis,1.e-9)
  call log_value("Estimated central B_T ",btcen)

  ! 1D work array
  !! allocate 1D work storage and read
  if(nw>0) then
     allocate(work1(nw), stat=status)
     call log_alloc_check(m_name,s_name,20,status)
  else
     call log_error(m_name,s_name,21,error_fatal,'No 1D data')
  end if

  do j=1,3
     if (cfmtd=='*') then
        read(iin,*,iostat=status)(work1(i),i=1,nw)
     else
        read(iin,cfmtd,iostat=status)(work1(i),i=1,nw)
     end if
     if(status/=0) then
        call log_error(m_name,s_name,21+j,error_fatal,'Error reading work1')
     end if
  end do

  self%mr=nw-1
  self%mz=nh-1
  ! work2 array
  !! allocate work2 storage
  if(nw*nh>0) then
     allocate(work2(nw,nh), stat=status)
     call log_alloc_check(m_name,s_name,40,status)
  else
     call log_error(m_name,s_name,41,error_fatal,'No psi data')
  end if

  if (cfmtd=='*') then
     read(iin,*,iostat=status) ((work2(i,j),i=1,nw),j=1,nh)
  else
     if (iffiesta>0) then
        do j=1,nh
           read(iin,cfmtd,iostat=status) (work2(i,j),i=1,nw)
        end do
     else
        read(iin,cfmtd,iostat=status) ((work2(i,j),i=1,nw),j=1,nh)
     end if
  end if
  if(status/=0) then
     call log_error(m_name,s_name,42,error_fatal,'Error reading psi')
  end if
  print '("number of psi values read = ",i10)',nw*nh
  call log_value("number of psi values read ",nw*nh)
  zpsic=work2(nw/2,nh/2)
  call log_value("array central raw psi",zpsic)
  if (zpsic*ssimag1<0) then
     call log_error(m_name,s_name,44,error_warning,'file may have inconsistent sign of central psi in header')
  end if
  !
  !      zpsimax=maxval(work2)
  !      zpsimin=minval(work2)
  !      write(*,*) 'zpsimin,zpsimax=',zpsimin,zpsimax
  !
  if (numerics%psibig>0) then
     psifac=1/(2*const_pid)
  else
     psifac=1
  end if
  self%psiaxis=ssimag1*psifac
  self%psiqbdry=ssibry1*psifac
  self%psibdry=self%psiqbdry
  self%rmin=rgrid
  self%rmax=rgrid+xdim
  self%zmin=zmid-zdim/2
  self%zmax=zmid+zdim/2
  self%rqcen=rmaxis
  self%zqcen=zmaxis

  if (.NOT.numerics%leqok) then
     if (BEQ_OVERRIDE_ITER) then
        call log_error(m_name,s_name,49,log_info,'override for ITER')
        ! special for ITER to align current and toroidal field
        if (ssimag1>ssibry1) then
           self%psiaxis=-ssimag1*psifac
           self%psiqbdry=-ssibry1*psifac
           self%psibdry=self%psiqbdry
           work2=-work2
        else
           ! possibly also in case where psi increases from axis
           self%f=-self%f
        end if
     else if (BEQ_OVERRIDE_FTU) then
        call log_error(m_name,s_name,48,log_info,'override for FTU')
        ! special for FTU modelling to align current and toroidal field
        ! possibly also in case where psi increases from axis
        self%f=-self%f
     end if
  end if

  !! now additional standard eqdsk stuff
  if (cfmtd=='*') then
     read(iin,*,iostat=status)(work1(i),i=1,nw)
  else
     read(iin,cfmtd,iostat=status)(work1(i),i=1,nw)
  end if
  if(status/=0) then
     call log_error(m_name,s_name,50,error_fatal,'Error reading qpsi')
  end if
  !DPR      write(*,*) 'qpsi', (i,work1(i),i=1,nw) !DPR
  !DPR      write(*,*) 'fpol', (i,self%f(i),i=1,nw) !DPR
  deallocate(work1)
  read(iin,*,iostat=istatus,end=1) nbbbs,limitr
  if(istatus/=0) then
     call log_error(m_name,s_name,51,error_fatal,'Error reading nbbbs,limitr')
  end if
  if(nbbbs>0) then
     allocate(workr1(nbbbs), workz1(nbbbs), stat=status)
     call log_alloc_check(m_name,s_name,52,status)
  else
     call log_error(m_name,s_name,53,error_warning,'No boundary data')
  end if
  if (cfmtd=='*') then
     read(iin,*,iostat=status)(workr1(i),workz1(i),i=1,nbbbs)
  else
     read(iin,cfmtd,iostat=status)(workr1(i),workz1(i),i=1,nbbbs)
  end if
  if(status/=0) then
     call log_error(m_name,s_name,54,error_fatal,'Error reading boundary data')
  end if
  if (numerics%eqbdry) then
     call gfile_init(numerics%eqbdryfile,'eqdsk (R,Z) boundary data',igunit)
     write(igunit,'(A,1PG12.5)') '#  ', nbbbs
     do i=1,nbbbs
        write(igunit,*) workr1(i),workz1(i)
     end do
     call gfile_close
  end if
  deallocate(workr1)
  deallocate(workz1)
  if(limitr>0) then
     allocate(workr1(limitr), workz1(limitr), stat=status)
     call log_alloc_check(m_name,s_name,55,status)
  else
     call log_error(m_name,s_name,56,error_warning,'No limiter data')
  end if
  if (cfmtd=='*') then
     read(iin,*,iostat=status)(workr1(i),workz1(i),i=1,limitr)
  else
     read(iin,cfmtd,iostat=status)(workr1(i),workz1(i),i=1,limitr)
  end if
  if(status/=0) then
     call log_error(m_name,s_name,57,error_fatal,'Error reading limiter data')
  end if
  if (numerics%eqbdry) then
     call gfile_init(numerics%eqltrfile,'eqdsk (R,Z) limiter data',igunit)
     write(igunit,'(A,1PG12.5)') '#  ', limitr
     do i=1,limitr
        write(igunit,*) workr1(i),workz1(i)
     end do
     call gfile_close
  end if
  deallocate(workr1)
  deallocate(workz1)
  !! now the changes for the kludged B output format
  ! error indicates not present
  !coil data
  read(iin,*,iostat=status) ncoil
  ! optionally skip B in file
  beq_nobinq=(status/=0).OR.numerics%skipb
  if(beq_nobinq) then
     call log_error(m_name,s_name,60,error_warning,'Unable to read ncoil')
     call log_value("Giving up on EQDSK for B values, status",status)
     !        deallocate(workr1)
     !        deallocate(workz1)
     return
  end if
  if(ncoil>0) then
     allocate(work1(5*ncoil), stat=status)
     call log_alloc_check(m_name,s_name,61,status)
  else
     call log_error(m_name,s_name,62,error_warning,'No coil data')
     allocate(work1(1), stat=status)
  end if
  if (cfmtd=='*') then
     read(iin,*,iostat=status)(work1(i),i=1,5*ncoil)
  else
     read(iin,cfmtd,iostat=status)(work1(i),i=1,5*ncoil)
  end if
  beq_nobinq=(status/=0)
  if(beq_nobinq) then
     call log_error(m_name,s_name,63,error_warning,'Error reading coils')
     call log_value("Giving up on EQDSK for B values, status",status)
     deallocate(work1)
     return
  end if
  deallocate(work1)
  ! finally the B fields
  ! workr2 array for Br
  !! allocate workr2 storage
  allocate(workr2(nw,nh), stat=status)
  call log_alloc_check(m_name,s_name,70,status)
  if (cfmtd=='*') then
     read(iin,*,iostat=status) ((workr2(i,j),i=1,nw),j=1,nh)
  else
     read(iin,cfmtd,iostat=status) ((workr2(i,j),i=1,nw),j=1,nh)
  end if
  beq_nobinq=(status/=0)
  if(beq_nobinq) then
     call log_error(m_name,s_name,71,error_warning,'Error reading Br')
     call log_value("Giving up on EQDSK for B values, status",status)
     deallocate(workr2)
     return
  end if
  print '("number of Br values read = ",i10)',nw*nh
  call log_value("number of Br values read ",nw*nh)
  ! workz2 array for Bz
  !! allocate workz2 storage
  allocate(workz2(nw,nh), stat=status)
  call log_alloc_check(m_name,s_name,72,status)
  if (cfmtd=='*') then
     read(iin,*,iostat=status) ((workz2(i,j),i=1,nw),j=1,nh)
  else
     read(iin,cfmtd,iostat=status,end=1) ((workz2(i,j),i=1,nw),j=1,nh)
  end if
  beq_nobinq=(status/=0)
  if(beq_nobinq) then
     call log_error(m_name,s_name,73,error_warning,'Error reading Bz')
     call log_value("Giving up on EQDSK for B values, status",status)
     deallocate(workr2)
     deallocate(workz2)
     return
  end if
  print '("number of Bz values read = ",i10)',nw*nh
  call log_value("number of Bz values read ",nw*nh)
  if (.NOT.numerics%leqok) then
     if (BEQ_OVERRIDE_ITER) then
        ! fix sign consistent with psi above
        if (ssimag1>ssibry1) then
           workr2=-workr2
           workz2=-workz2
        end if
     end if
  end if
  close(iin)

  return

1     continue
  print '("Eqdsk file ended status ",i10," but may have got enough data.")',status
  call log_error(m_name,s_name,80,error_warning,'Unexpected end of eqdsk file')
  call log_value("May have got enough data. Termination status ",status)

end subroutine beq_readequil
!---------------------------------------------------------------------
!> read in EQU format
subroutine beq_readequ(self,infile,numerics)

  !! arguments
  type(beq_t), intent(out) :: self   !< object data structure
  character(*),intent(in) :: infile !< name of input file
  type(bnumerics_t), intent(in) :: numerics   !< object control data structure

  !! local
  character(*), parameter :: s_name='beq_readequ' !< subroutine name
  character(8), parameter :: cfmt1='(5f17.8)'
  !! logical :: unitused !< flag to test unit is available
  integer :: iin   !< input channel for object data structure
  integer(ki4) :: jm  !< FIESTA number of grid points in R-direction
  integer(ki4) :: km  !< FIESTA number of grid points in Z-direction
  real(kr8) :: btf !< FIESTA vacuum toroidal field at r=rtf
  real(kr8) :: rtf !< FIESTA vacuum toroidal field position
  real(kr8) :: psib !< FIESTA poloidal flux at the boundary (Wb/rad)
  real(kr8) :: psic !< madeup poloidal flux at the plasma centre (Wb/rad), to force on-axis minimum/maximum

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then iin=i exit end if end do

  call misc_getfileunit(iin)
  open(unit=iin, file=infile,status='OLD',form='FORMATTED',iostat=status)
  !!eof or error
  if (status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error opening file')
  end if
  ! skip header data
  do
     read(iin,fmt='(a)',iostat=status) ibuff
     call log_read_check(m_name,s_name,2,status)
     if (adjustl(ibuff)=='jm') exit
  end do

  if (adjustl(ibuff)/='jm') then
     ! jm string not found
     call log_error(m_name,s_name,3,error_fatal,'Error reading object data')
  end if

  read(iin,*,iostat=status) jm
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error reading object data')
  end if
  read(iin,*,iostat=status) ibuff
  read(iin,*,iostat=status) km
  if(status/=0) then
     call log_error(m_name,s_name,5,error_fatal,'Error reading object data')
  end if
  read(iin,*,iostat=status) ibuff
  read(iin,*,iostat=status) psib
  if(status/=0) then
     call log_error(m_name,s_name,6,error_fatal,'Error reading object data')
  end if
  read(iin,*,iostat=status) ibuff
  read(iin,*,iostat=status) psic
  if(status/=0) then
     call log_error(m_name,s_name,7,error_fatal,'Error reading object data')
  end if
  read(iin,*,iostat=status) ibuff
  read(iin,*,iostat=status) btf
  if(status/=0) then
     call log_error(m_name,s_name,8,error_fatal,'Error reading object data')
  end if
  read(iin,*,iostat=status) ibuff
  read(iin,*,iostat=status) rtf
  if(status/=0) then
     call log_error(m_name,s_name,9,error_fatal,'Error reading object data')
  end if
  read(iin,*,iostat=status) ibuff
  !-----------------------------------------------------------------------
  !              Allocate 1D storage and read
  ! position 1 array
  !! allocate position 1 storage
  if(jm>0) then
     allocate(workr1(jm), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
  else
     call log_error(m_name,s_name,11,error_fatal,'No 1D data')
  end if

  read(iin,cfmt1,iostat=status)(workr1(i),i=1,jm)
  if(status/=0) then
     call log_error(m_name,s_name,12,error_fatal,'Error reading R values')
  end if
  print '("number of R values read = ",i10)',jm
  call log_value("number of R values read ",jm)
  read(iin,*,iostat=status) ibuff

  ! position 2 array
  !! allocate position 2 storage
  if(km>0) then
     allocate(workz1(km), stat=status)
     call log_alloc_check(m_name,s_name,15,status)
  else
     call log_error(m_name,s_name,16,error_fatal,'No 1D data')
  end if

  read(iin,cfmt1,iostat=status)(workz1(i),i=1,km)
  if(status/=0) then
     call log_error(m_name,s_name,17,error_fatal,'Error reading Z values')
  end if
  print '("number of Z values read = ",i10)',km
  call log_value("number of Z values read ",km)
  read(iin,*,iostat=status) ibuff

  !-----------------------------------------------------------------------
  !             Read 2D
  self%mr=jm-1
  self%mz=km-1
  ! work2 array
  !! allocate work2 storage
  if(jm*km>0) then
     allocate(work2(jm,km), stat=status)
     call log_alloc_check(m_name,s_name,20,status)
  else
     call log_error(m_name,s_name,21,error_fatal,'No psi data')
  end if

  read(iin,cfmt1,iostat=status) ((work2(i,j),i=1,jm),j=1,km)
  if(status/=0) then
     call log_error(m_name,s_name,22,error_fatal,'Error reading psi')
  end if
  print '("number of psi values read = ",i10)',jm*km
  call log_value("number of psi values read ",jm*km)

  self%psiaxis=psic
  self%psiqbdry=psib
  self%psibdry=self%psiqbdry
  self%rmin=workr1(1)
  self%rmax=workr1(jm)
  self%zmin=workz1(1)
  self%zmax=workz1(km)
  ! make up a value for plasma centre
  self%rqcen=(workr1(1)+workr1(jm))/2
  self%zqcen=(workz1(1)+workz1(km))/2 ! km not jm??

  beq_nobinq=.TRUE.

  ! set up fpol
  !! allocate fpol storage
  allocate(self%f(jm), stat=status)
  call log_alloc_check(m_name,s_name,30,status)
  !! vacuum field assumption
  self%f=rtf*btf
  self%ivac=rtf*btf

  fld_specn: select case (numerics%fldspec)
  case(1,3)
     if (.NOT.numerics%leqok) then
        if (BEQ_OVERRIDE_ITER) then
           call log_error(m_name,s_name,49,log_info,'override for ITER')
           ! special for ITER to align current and toroidal field
           if (psic>psib) then
              self%psiaxis=-psic
              self%psiqbdry=-psib
              self%psibdry=self%psiqbdry
              work2=-work2
           else
              ! possibly also in case where psi increases from axis
              self%f=-self%f
           end if
        else if (BEQ_OVERRIDE_FTU) then
           call log_error(m_name,s_name,50,log_info,'override for FTU')
           ! special for FTU modelling to align current and toroidal field
           ! possibly also in case where psi increases from axis
           self%f=-self%f
        end if
     end if
  case(2)
     ! Separating this case enables MAST test deck case to work,
     ! without setting BEQ_OVERRIDE_ITER=.FALSE., which is what
     ! really should be done, and no special fldspec test
     ! This works because .equ files only ever used for MAST cases
  end select fld_specn

  close(iin)

end subroutine beq_readequ
!---------------------------------------------------------------------
!> generate using ANA format
subroutine beq_readana(self,beqan,kfldspec)

  !! arguments
  type(beq_t), intent(out) :: self   !< object data structure
  type(beqan_t), intent(inout) :: beqan   !< object data structure
  integer(ki4),intent(in) :: kfldspec !< field as mapped or 3-cpt

  !! local
  character(*), parameter :: s_name='beq_readana' !< subroutine name
  integer(ki4) :: nw  !< number of grid points in R-direction
  integer(ki4) :: nh  !< number of grid points in Z-direction
  real(kr8) :: itf !< vacuum R * toroidal field at r=rtf
  real(kr8), dimension(3) :: zx !< position vector
  real(kr8), dimension(3) :: zB !< field vector
  real(kr8) :: zpsi    !<  \f$ \psi \f$
  real(kr8) :: psic    !<  \f$ \psi \f$ at centre
  real(kr8) :: psib    !<  \f$ \psi \f$ at boundary

  self%mr=beqan%mr
  self%mz=beqan%mz
  self%rmin=beqan%rmin
  self%rmax=beqan%rmax
  self%zmin=beqan%zmin
  self%zmax=beqan%zmax
  self%dr=(self%rmax-self%rmin)/self%mr
  self%dz=(self%zmax-self%zmin)/self%mz

  !-----------------------------------------------------------------------
  !              Allocate 1D storage
  nw=beqan%mr+1
  ! position 1 array
  !! allocate position 1 storage
  if(nw>0) then
     allocate(workr1(nw), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
  else
     call log_error(m_name,s_name,11,error_fatal,'No 1D data')
  end if

  nh=beqan%mz+1
  ! position 2 array
  !! allocate position 2 storage
  if(nh>0) then
     allocate(workz1(nh), stat=status)
     call log_alloc_check(m_name,s_name,15,status)
  else
     call log_error(m_name,s_name,16,error_fatal,'No 1D data')
  end if

  psic=beqan%psiaxis
  psib=beqan%psibdry

  !-----------------------------------------------------------------------
  !             2D array initialisation
  ! work2 array
  !! allocate work2 storage for psi
  if(nw*nh>0) then
     allocate(work2(nw,nh), stat=status)
     call log_alloc_check(m_name,s_name,20,status)
  else
     call log_error(m_name,s_name,21,error_fatal,'No psi data')
  end if

  beq_nobinq=.FALSE.
  ! finally the B fields
  ! workr2 array for Br
  !! allocate workr2 storage
  allocate(workr2(nw,nh), stat=status)
  call log_alloc_check(m_name,s_name,70,status)
  ! workz2 array for Bz
  !! allocate workz2 storage
  allocate(workz2(nw,nh), stat=status)
  call log_alloc_check(m_name,s_name,72,status)

  do j=1,nh
     workz1(j)=self%zmin+(j-1)*self%dz
  end do
  zx(2)=0

  do i=1,nw
     workr1(i)=self%rmin+(i-1)*self%dr
     zx(1)=workr1(i)
     do j=1,nh
        zx(3)=workz1(j)
        call beqan_solovev(beqan,zx,zB,zpsi)
        work2(i,j)=zpsi
        workr2(i,j)=zB(1)
        workz2(i,j)=zB(2)
     end do
  end do

  ! set up fpol
  !! allocate fpol storage
  allocate(self%f(nw), stat=status)
  call log_alloc_check(m_name,s_name,30,status)
  itf=zB(3)
  !! effectively vacuum field assumption
  self%f=itf
  self%ivac=itf

  ! make up a value for plasma centre
  self%rqcen=(workr1(1)+workr1(nw))/2
  self%zqcen=(workz1(1)+workz1(nh))/2

  self%psiaxis=psic
  self%psibdry=psib
  self%psiqbdry=self%psibdry

  if (BEQ_OVERRIDE_ITER) then
     call log_error(m_name,s_name,49,log_info,'override for ITER')
     ! special for ITER to align current and toroidal field
     if (psic>psib) then
        self%psiaxis=-psic
        self%psibdry=-psib
        self%psiqbdry=self%psibdry
        work2=-work2
     else
        ! possibly also in case where psi increases from axis
        self%f=-self%f
     end if
  else if (BEQ_OVERRIDE_FTU) then
     call log_error(m_name,s_name,50,log_info,'override for FTU')
     ! special for FTU modelling to align current and toroidal field
     ! possibly also in case where psi increases from axis
     self%f=-self%f
  end if
  self%psicen=self%psiaxis

end subroutine beq_readana
!---------------------------------------------------------------------
!> fix radial origin of  beq data structure
subroutine beq_fixorigin(self,numerics)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  type(bnumerics_t), intent(inout) :: numerics   !< object control data structure

  !! local
  character(*), parameter :: s_name='beq_fixorigin' !< subroutine name
  real(kr8) :: zdr !< mesh radial displacement
  integer(ki4) :: ipos !< new radial start
  integer(ki4) :: iw !< new radial total of samples
  integer(ki4) :: ih !< vertical dimension of array

  ! check for a problem, none if rmin>0
  if (self%rmin>epsg*(self%rmax-self%rmin)) return

  ! problem, strip out rows with zero and negative R
  call log_value("radial origin of beq too small, fixing up ",self%rmin)
  zdr=(self%rmax-self%rmin)/self%mr
  ipos=max(1.05_kr8,-self%rmin/zdr)
  iw=self%mr+1-ipos
  ih=self%mz+1

  ! first fix flux function array
  allocate(zwork2(iw,ih), stat=status)
  call log_alloc_check(m_name,s_name,10,status)
  do i=1,iw
     zwork2(i,:)=work2(i+ipos,:)
  end do
  deallocate(work2)
  allocate(work2(iw,ih))
  call log_alloc_check(m_name,s_name,11,status)
  work2=zwork2
  deallocate(zwork2)

  if (.NOT.beq_nobinq) then
     ! repeat fix for B component arrays
     allocate(zwork2(iw,ih), stat=status)
     call log_alloc_check(m_name,s_name,20,status)
     do i=1,iw
        zwork2(i,:)=workr2(i+ipos,:)
     end do
     deallocate(workr2)
     allocate(workr2(iw,ih))
     call log_alloc_check(m_name,s_name,21,status)
     workr2=zwork2
     deallocate(zwork2)
     !
     allocate(zwork2(iw,ih), stat=status)
     call log_alloc_check(m_name,s_name,22,status)
     do i=1,iw
        zwork2(i,:)=workz2(i+ipos,:)
     end do
     deallocate(workz2)
     allocate(workz2(iw,ih))
     call log_alloc_check(m_name,s_name,23,status)
     workz2=zwork2
     deallocate(zwork2)

  end if

  self%mr=iw-1
  self%rmin=self%rmin+ipos*zdr
  call log_value("radial origin now ",self%rmin)
  call log_value("number of columns of data removed ",ipos)

end subroutine beq_fixorigin
!---------------------------------------------------------------------
!> move  beq data structure (controls)
subroutine beq_move(self,numerics)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  type(bnumerics_t), intent(in) :: numerics   !< object control data structure

  !! local
  character(*), parameter :: s_name='beq_move' !< subroutine name
  real(kr8) :: zgfac    !< geometrical factor

  ! move (R,Z) values
  self%rmin=self%rmin+numerics%rmove
  self%rmax=self%rmax+numerics%rmove
  self%zmin=self%zmin+numerics%zmove
  self%zmax=self%zmax+numerics%zmove
  self%rqcen=self%rqcen+numerics%rmove
  self%zqcen=self%zqcen+numerics%zmove
  zgfac=self%rqcen/(self%rqcen-numerics%rmove)
  ! adjust \f$ f \f$ approximately (should use R instead of R_c,
  ! but then f becomes a function of \f$ \theta \f$)
  self%f=self%f*zgfac
  ! scale f
  self%f=self%f*numerics%fscale

end subroutine beq_move
!---------------------------------------------------------------------
!> create  beq data structure
subroutine beq_init(self,numerics,fmesh)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  type(bnumerics_t), intent(in) :: numerics   !< object control data structure
  type(fmesh_t), intent(in) :: fmesh   !< field mesh data structure


  !! local
  character(*), parameter :: s_name='beq_init' !< subroutine name

  ! initialise control structure
  self%n=numerics
  self%fmesh=fmesh
  self%replasi=.FALSE.

  ! set flag for whether psi increases or decreases with minor radius
  !     rsig=sign(1._kr8,self%psiqbdry-self%psiaxis)
  call beq_rsigset(self)
  ! set psi normalisation
  self%psinorm=abs(self%psiqbdry-self%psiaxis)/2

  ! resolve psi boundary definition as far as possible
  if(self%n%bdryopt==1.OR.self%n%bdryopt==5.OR.self%n%bdryopt==9) then
     self%psibdry=self%n%psiref
  else
     ! default is boundary-point value set in readq (EQDSK)
  end if
  ! initialise 2-D spline
  self%dr=(self%rmax-self%rmin)/self%mr
  self%dz=(self%zmax-self%zmin)/self%mz
  if (self%n%psibig>0) then
     ! scale psi
     work2=work2/(2*const_pid)
     call log_error(m_name,s_name,1,log_info,'Scaling psi by 2pi')
     ! and stop it happening again
     self%n%psibig=-self%n%psibig
  end if
  call spl2d_init(self%psi,work2,self%mr,self%mz,&
 &self%rmin,self%zmin,self%dr,self%dz,beq_spline_order)
  deallocate(work2)

  if (beq_nobinq) then
     call log_value("Calculating B from psi, rsig",rsig)
     ! and derivatives, for which need coeffs
     call spl2d_coeff(self%psi)
     call spl2d_deriv(self%psi,self%dpsidr,1)
     call spl2d_deriv(self%psi,self%dpsidz,2)
  else
     call log_value("Read B from disk, rsig",rsig)
     ! convert to dpsidz,dpsidr from Br,Bz respectively
     do j=1,self%mz+1
        do i=1,self%mr+1
           workr2(i,j)=-workr2(i,j)*(self%psi%pos1(i)-self%n%rmove)
           workz2(i,j)=workz2(i,j)*(self%psi%pos1(i)-self%n%rmove)
        end do
     end do
     call spl2d_init(self%dpsidz,workr2,self%mr,self%mz,&
 &   self%rmin,self%zmin,self%dr,self%dz,beq_spline_order)
     deallocate(workr2)
     call spl2d_init(self%dpsidr,workz2,self%mr,self%mz,&
 &   self%rmin,self%zmin,self%dr,self%dz,beq_spline_order)
     deallocate(workz2)
  end if

  ! initialise ripple field
  zdum=beq_ripple_h1(0._kr8,0._kr8,0._kr8,-1,self%n%mrip,self%ivac,self%n%arip)

end subroutine beq_init
!---------------------------------------------------------------------
!> helical sense of field
subroutine beq_sense(self,k3d)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  integer(ki4), intent(in) :: k3d   !< selector

  !! local
  character(*), parameter :: s_name='beq_sense' !< subroutine name
  type(posang_t) :: posang !< position and vector involving angles
  integer(ki4) :: isleft !< unity if left-handed helix

  !write(*,*) self%n%vacfile
  if ( (k3d==0.AND.self%n%vacfile/='null') .OR. &
 &(k3d==1.AND.self%n%vacfile=='null') ) return
  ! sense (R,Z) values
  posang%pos(1)=self%rmin+0.7*(self%rmax-self%rmin)
  posang%pos(2)=0
  posang%pos(3)=self%zmin+0.5*(self%zmax-self%zmin)
  posang%opt=0 ; posang%units=0
  call log_value("Sensing B at X ",posang%pos(1))
  call log_value("Sensing B at Z ",posang%pos(3))
  call beq_b(self,posang,0)
  call log_value("Sensed By, on Y=0, ",posang%vec(2))
  call log_value("Sensed Bz, at nominal central Z, ",posang%vec(3))
  isleft=sign(1.05_kr8,posang%vec(3)/posang%vec(2))
  if (isleft==1) then
     call log_error(m_name,s_name,1,log_info,'Left-handed helical field UNlike ITER')
  else
     call log_error(m_name,s_name,2,log_info,'Right-handed helical field like ITER')
  end if

end subroutine beq_sense
!---------------------------------------------------------------------
!> delete  beq data structure
subroutine beq_delete(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_delete' !< subroutine name

  deallocate(self%f)
  deallocate(self%i)

  call spl2d_delete(self%psi)
  call spl2d_delete(self%dpsidr)
  call spl2d_delete(self%dpsidz)
  fld_specn: select case (self%n%fldspec)
  case(1)
     call spl2d_delete(self%r)
     call spl2d_delete(self%z)
     call spl2d_delete(self%rjac)
  case(2)
  case(3)
     call spl2d_delete(self%rispldr)
     call spl2d_delete(self%rispldz)
  end select fld_specn

end subroutine beq_delete
!---------------------------------------------------------------------
!> read in visualisation format
subroutine beq_readv(self)

  !! arguments
  type(beq_t), intent(out) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_readv' !< subroutine name

  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading vec')
  end if

end subroutine beq_readv
!---------------------------------------------------------------------
!> check field as mapped or 3-cpt and for extras
subroutine beq_readcheck(self,infile,kextra)

  !! arguments
  type(beq_t), intent(out) :: self   !< object data structure
  character(*),intent(in) :: infile !< name of input file
  integer(ki4), intent(out) :: kextra !< extra info extracted

  !! local
  character(*), parameter :: s_name='beq_readcheck' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer(ki4) :: ifldspec !< field spec

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

  !! open file
  call misc_getfileunit(nin)
  open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,1,error_fatal,'Error opening beq data structure file')
  else
     call log_error(m_name,s_name,2,log_info,'beq data structure file opened')
  end if

  ! skip header data
  do i=1,20
     read(nin,fmt='(a)',iostat=status) ibuff
     call log_read_check(m_name,s_name,3,status)
     if (adjustl(ibuff)=='fldspec') exit
  end do

  if (adjustl(ibuff)/='fldspec') then
     ! fldspec string not found
     call log_error(m_name,s_name,3,error_warning,'Object data is mapped')
     ifldspec=1
  else

     read(nin,*,iostat=status) ifldspec
     call log_read_check(m_name,s_name,4,status)
  end if

  kextra=ifldspec/10
  self%n%fldspec=ifldspec-10*kextra

  close(nin)

end subroutine beq_readcheck
!---------------------------------------------------------------------
!> read field data needed by ITER
subroutine beq_readpart(self,infile)

  !! arguments
  type(beq_t), intent(out) :: self   !< object data structure
  character(*),intent(in) :: infile !< name of input file

  !! local
  character(*), parameter :: s_name='beq_readpart' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

  !! open file
  call misc_getfileunit(nin)
  open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,1,error_fatal,'Error opening beq data structure file')
  else
     call log_error(m_name,s_name,2,log_info,'beq data structure file opened')
  end if

  ! skip header data
  do
     read(nin,fmt='(a)',iostat=status) ibuff
     call log_read_check(m_name,s_name,3,status)
     if (adjustl(ibuff)=='mr') exit
  end do

  if (adjustl(ibuff)/='mr') then
     ! mr string not found
     call log_error(m_name,s_name,3,error_fatal,'Error reading object data')
  end if

  read(nin,*,iostat=status) self%mr
  call log_read_check(m_name,s_name,4,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%rcen
  call log_read_check(m_name,s_name,17,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%zcen
  call log_read_check(m_name,s_name,18,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%psiaxis
  call log_read_check(m_name,s_name,34,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%psibdry
  call log_read_check(m_name,s_name,35,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%psiqbdry
  call log_read_check(m_name,s_name,36,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%rbdry
  call log_read_check(m_name,s_name,37,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%bpbdry
  call log_read_check(m_name,s_name,38,status)
  read(nin,*,iostat=status) ibuff
  if (adjustl(ibuff)=='btotbdry') then
     read(nin,*,iostat=status) self%btotbdry
     call log_read_check(m_name,s_name,39,status)
     read(nin,*,iostat=status) ibuff
  end if
  !-----------------------------------------------------------------------
  !              Allocate 1D storage and read
  ! position 1 array
  !! allocate position 1 storage
  allocate(self%f(self%mr+1), stat=status)
  call log_alloc_check(m_name,s_name,50,status)

  read(nin,*,iostat=status) self%f
  call log_read_check(m_name,s_name,51,status)

  !-----------------------------------------------------------------------
  !             Read 2D splines
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%r,infile,nin )
  call spl2d_initpart( self%r )
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%z,infile,nin )
  call spl2d_initpart( self%z )
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%rjac,infile,nin )
  call spl2d_initpart( self%rjac )
  read(nin,*,iostat=status) ibuff

  if(status/=0)then
     call log_error(m_name,s_name,60,log_info,'beq read in from data file')
  end if

  ! added 'plus parts'
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%fldspec
  call log_read_check(m_name,s_name,4,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%vacfile
  call log_read_check(m_name,s_name,63,status)
  ! ripple field data
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%mrip
  call log_read_check(m_name,s_name,64,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%ivac
  call log_read_check(m_name,s_name,65,status)
  if (self%n%mrip/=0) then
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%n%arip
     call log_read_check(m_name,s_name,66,status)
  end if

  call log_error(m_name,s_name,70,log_info,'beq plus 2 read in from data file')

end subroutine beq_readpart
!---------------------------------------------------------------------
!> read field data needed by 'msus','global'
subroutine beq_readplus(self,infile)

  !! arguments
  type(beq_t), intent(out) :: self   !< object data structure
  character(*), intent(in) :: infile !< name of input file

  !! local
  character(*), parameter :: s_name='beq_readplus' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  integer(ki4) :: ifldspec !< field as mapped or 3-cpt + extra info
  integer(ki4) :: iextra !< extra info extracted

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

  !! open file
  call misc_getfileunit(nin)
  open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,1,error_fatal,'Error opening beq data structure file')
  else
     call log_error(m_name,s_name,2,log_info,'beq data structure file opened')
  end if

  ! skip header data
  do
     read(nin,fmt='(a)',iostat=status) ibuff
     call log_read_check(m_name,s_name,3,status)
     if (adjustl(ibuff)=='fldspec') exit
  end do

  if (adjustl(ibuff)/='fldspec') then
     ! fldspec string not found
     call log_error(m_name,s_name,3,error_fatal,'Error reading object data')
  end if

  read(nin,*,iostat=status) self%n%fldspec
  call log_read_check(m_name,s_name,4,status)
  ifldspec=self%n%fldspec
  iextra=ifldspec/10
  self%n%fldspec=ifldspec-10*iextra
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%mr
  call log_read_check(m_name,s_name,5,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%rcen
  call log_read_check(m_name,s_name,17,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%zcen
  call log_read_check(m_name,s_name,18,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%psiaxis
  call log_read_check(m_name,s_name,34,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%psibdry
  call log_read_check(m_name,s_name,35,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%psiqbdry
  call log_read_check(m_name,s_name,36,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%rbdry
  call log_read_check(m_name,s_name,37,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%bpbdry
  call log_read_check(m_name,s_name,38,status)
  read(nin,*,iostat=status) ibuff
  if (adjustl(ibuff)=='btotbdry') then
     read(nin,*,iostat=status) self%btotbdry
     call log_read_check(m_name,s_name,39,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%n%zetamin
     call log_read_check(m_name,s_name,40,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%n%zetamax
     call log_read_check(m_name,s_name,41,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%nzets
     call log_read_check(m_name,s_name,42,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%n%ximin
     call log_read_check(m_name,s_name,43,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%n%ximax
     call log_read_check(m_name,s_name,44,status)
     read(nin,*,iostat=status) ibuff
  end if

  !-----------------------------------------------------------------------
  !              Allocate 1D storage and read
  ! position 1 array
  !! allocate position 1 storage
  allocate(self%f(self%mr+1), stat=status)
  call log_alloc_check(m_name,s_name,50,status)

  read(nin,*,iostat=status) self%f
  call log_read_check(m_name,s_name,51,status)

  !-----------------------------------------------------------------------
  !             Read 2D splines
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%psi,infile,nin )
  call spl2d_initpart( self%psi )
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%dpsidr,infile,nin )
  call spl2d_initpart( self%dpsidr )
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) ibuff
  call spl2d_read( self%dpsidz,infile,nin )
  call spl2d_initpart( self%dpsidz )
  read(nin,*,iostat=status) ibuff

  fld_specn: select case (self%n%fldspec)
  case(3)
     read(nin,*,iostat=status) ibuff
     call spl2d_read( self%rispldr,infile,nin )
     call spl2d_initpart( self%rispldr )
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) ibuff
     call spl2d_read( self%rispldz,infile,nin )
     call spl2d_initpart( self%rispldz )
     read(nin,*,iostat=status) ibuff
  end select fld_specn


  ! ripple field data
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%mrip
  call log_read_check(m_name,s_name,60,status)
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%ivac
  call log_read_check(m_name,s_name,61,status)
  if (self%n%mrip/=0) then
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%n%arip
     call log_read_check(m_name,s_name,62,status)
  end if
  read(nin,*,iostat=status) ibuff
  read(nin,*,iostat=status) self%n%vacfile
  call log_read_check(m_name,s_name,63,status)

  if (iextra==2.OR.iextra==4) then
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%rxpt
     call log_read_check(m_name,s_name,64,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%zxpt
     call log_read_check(m_name,s_name,65,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%psixpt
     call log_read_check(m_name,s_name,66,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%zmin
     call log_read_check(m_name,s_name,67,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%zmax
     call log_read_check(m_name,s_name,68,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%mz
     call log_read_check(m_name,s_name,69,status)
     read(nin,*,iostat=status) ibuff
     read(nin,*,iostat=status) self%dz
     call log_read_check(m_name,s_name,70,status)
     read(nin,*,iostat=status) ibuff
     call log_read_check(m_name,s_name,71,status)
     read(nin,*,iostat=status) self%nctrack
     call log_read_check(m_name,s_name,72,status)
     read(nin,*,iostat=status) ibuff
     call log_read_check(m_name,s_name,73,status)
     allocate(self%ctrackrz(self%nctrack,2), stat=status)
     call log_alloc_check(m_name,s_name,74,status)
     read(nin,*,iostat=status) ((self%ctrackrz(i,j),i=1,self%nctrack),j=1,2)
     call log_read_check(m_name,s_name,75,status)
     self%n%skylpsi=.TRUE.
  else
     self%n%skylpsi=.FALSE.
     ! dummy definitions
     self%rxpt=0
     self%zxpt=0
     self%psixpt=0
     self%zmin=0
     self%zmax=1
     self%dz=1
     self%mz=1
  end if

  if (iextra==1) then
     self%n%duct=.TRUE.
     call fmesh_read(self%fmesh,infile,nin)
  else
     self%n%duct=.FALSE.
  end if

  call log_error(m_name,s_name,90,log_info,'beq read in from data file')

end subroutine beq_readplus
!---------------------------------------------------------------------
!> write in gnuplot format
subroutine beq_writeg(self,kchar,kout)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  character(*), intent(in) :: kchar  !< case
  integer, intent(in) :: kout   !< output channel for object data structure

  !! local
  character(*), parameter :: s_name='beq_writeg' !< subroutine name
  real(kr8) :: dzeta !< increment
  real(kr8) :: zr    !<   \f$ R(\psi_i,\theta_j) \f$
  real(kr8) :: zz    !<   \f$ Z(\psi_i,\theta_j) \f$
  real(kr8) :: zeta    !<  \f$ \zeta_k \f$
  real(kr8) :: zcos    !<  \f$ \cos(\zeta_k) \f$
  real(kr8) :: zsin    !<  \f$ \sin(\zeta_k) \f$
  type(posang_t) :: posang !< position and vector involving angles

  dzeta=2*const_pid/max(1,self%n%mzetag-1)

  plot_type: select case (kchar)
  case('R-Z')

     do i=1,self%r%n1p
        do j=1,self%r%n2p
           !     write(kout,'(2(1x,i4),4(2x,g12.5))',iostat=status) &
           write(kout,'(2(1x,i4),'//cfmt2v,iostat=status) &
 &         i-1,j-1,self%r%pos1(i),self%r%pos2(j),self%r%sampl(i,j),self%z%sampl(i,j)
           call log_write_check(m_name,s_name,1,status)
        end do
        write(kout,*,iostat=status) ' '
        call log_write_check(m_name,s_name,2,status)
     end do

  case('allcartesian')
     write(kout,'(a9,'//cfmt2v,iostat=status) &
 &   '# limits ',self%rmin,self%rmax,self%zmin,self%zmax
     call log_write_check(m_name,s_name,3,status)
     write(kout,'(a)',iostat=status) &
 &   '#X                Y               Z               BX               BY               BZ'
     call log_write_check(m_name,s_name,4,status)
     do k=1,self%n%mzetag
        zeta=(k-1)*dzeta
        zsin=sin(zeta)
        zcos=cos(zeta)
        do j=1,self%mz+1
           zz=self%zmin+(j-1)*self%dz
           do i=1,self%mr+1
              zr=self%rmin +(i-1)*self%dr
              posang%pos(1)=zr ; posang%pos(2)=zz ; posang%pos(3)=zeta
              posang%opt=1 ; posang%units=0
              ! posn to cartesians and output posn and vector
              call posang_tfm(posang,0)
              call beq_b(self,posang,0)
              call posang_units(posang,0)
              call posang_writev(posang,kout,3)
           end do
           write(kout,*,iostat=status) ' '
           call log_write_check(m_name,s_name,5,status)
        end do
        write(kout,*,iostat=status) ' '
        call log_write_check(m_name,s_name,6,status)
     end do

  case default

  end select plot_type

end subroutine beq_writeg
!---------------------------------------------------------------------
!> write in nucode format
subroutine beq_writen(self,kchar,kfield)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  character(*), intent(in) :: kchar  !< case
  integer :: kfield   !< output channel for nucode data

  !! local
  character(*), parameter :: s_name='beq_writen' !< subroutine name
  real(kr8) :: zx    !<  \f$ \x_i \f$
  real(kr8) :: zy    !<  \f$ \y_j \f$
  real(kr8) :: zz    !<  \f$ \z_k \f$
  type(posang_t) :: posang !< position and vector involving angles
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zpostfm !< local variable
  type(tfmdata_t) :: ztfmdata !< local variable

  !use duct not geometry (CAD)  coordinates
  ztfmdata=self%fmesh%tfmdata

  plot_type: select case (kchar)
  case('regular')
     ! regular for duct, which is actually irregular in that x and y loops transposed
     call log_error(m_name,s_name,1,log_info,'Output on Cartesian grid')
     call log_value("Number of sample points",self%fmesh%nxf*self%fmesh%nyf*self%fmesh%nzf)

     rewind(kfield)

     write(kfield,'(I8)') self%fmesh%nxf
     write(kfield,'(I8)') self%fmesh%nyf
     write(kfield,'(I8)') self%fmesh%nzf

     do k=1,self%fmesh%nzf
        zz=self%fmesh%zf(k)
        do j=1,self%fmesh%nxf
           zx=self%fmesh%xf(j)
           do i=1,self%fmesh%nyf
              zy=self%fmesh%yf(i)
              posang%pos(1)=zx ; posang%pos(2)=zy ; posang%pos(3)=zz
              posang%opt=0 ; posang%units=0
              ! transform point from duct to geometry (CAD) coordinates
              call posang_cartfm(posang,ztfmdata,0)
              ! calculate B in geometry (CAD)  coordinates/components
              call beq_b(self,posang,0)
              ! transform B to duct coordinates
              call posang_invcartfm(posang,ztfmdata,0)
              call posang_writev(posang,kfield,2)
           end do
        end do
     end do

  case default

     do k=1,self%fmesh%nzf
        zz=self%fmesh%zf(k)
        do j=1,self%fmesh%nyf
           zy=self%fmesh%yf(j)
           do i=1,self%fmesh%nxf
              zx=self%fmesh%xf(i)
              posang%pos(1)=zx ; posang%pos(2)=zy ; posang%pos(3)=zz
              posang%opt=0 ; posang%units=0
              ! transform point from duct to geometry (CAD) coordinates
              call posang_cartfm(posang,ztfmdata,0)
              ! calculate B in geometry (CAD)  coordinates/components
              call beq_b(self,posang,0)
              ! transform B to duct coordinates
              call posang_invcartfm(posang,ztfmdata,0)
              call posang_writev(posang,kfield,2)
           end do
        end do
     end do

  end select plot_type


end subroutine beq_writen
!---------------------------------------------------------------------
!> write in visualisation format
subroutine beq_writev(self,kchar,kplot)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  character(*), intent(in) :: kchar  !< case
  integer :: kplot   !< output channel for vis. data

  !! local
  character(*), parameter :: s_name='beq_writev' !< subroutine name
  integer(ki4) :: intheta    !< actual number of angles defining  \f$ R,Z(\psi,\theta) \f$
  real(kr8) :: zpsi    !<  \f$ \psi_i \f$
  real(kr8) :: ztheta    !<  \f$ \theta_j \f$
  real(kr8) :: zeta    !<  \f$ \zeta_k \f$
  real(kr8) :: zcos    !<  \f$ \cos(\zeta_k) \f$
  real(kr8) :: zsin    !<  \f$ \sin(\zeta_k) \f$
  real(kr8) :: ztmin    !<  minimum \f$ \theta \f$ for \f$ R/J(\psi,\theta) \f$
  real(kr8) :: zdpdr    !< \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !< \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zr    !<   \f$ R(\psi_i,\theta_j) \f$
  real(kr8) :: zz    !<   \f$ Z(\psi_i,\theta_j) \f$
  real(kr8) :: zc1    !<  1-component of output vector
  real(kr8) :: zc2    !<  2-component of output vector
  real(kr8) :: zc3    !<  3-component of output vector
  real(kr8) :: zf    !<   \f$ f(\psi) = RB_T \f$
  real(kr8) :: plotfac    !<   scale to mm
  type(posang_t) :: posang !< position and vector involving angles
  real(kr8) :: dzeta !< increment

  dzeta=2*const_pid/max(1,self%n%mzetav-1)
  plot_type: select case (kchar)
  case('cartesian')
     call log_error(m_name,s_name,1,log_info,'Plot parameter')
     call log_value("Number of toroidal points",self%n%mzetav)
     intheta=self%ntmax-self%ntmin+1
     write(kplot,'(''DATASET STRUCTURED_GRID'')')
     write(kplot,'(''DIMENSIONS '',I8,I8,I8)') self%n%npsi+1,intheta,self%n%mzetav
     write(kplot,'(''POINTS '',I8,'' float'')') (self%n%npsi+1)*intheta*self%n%mzetav

     do k=1,self%n%mzetav
        zeta=(k-1)*dzeta ; zsin=sin(zeta) ; zcos=cos(zeta)
        do j=self%ntmin,self%ntmax
           ztheta=self%n%thetamin+(j-1)*self%dtheta
           do i=1,self%n%npsi+1
              zpsi=self%n%psimin +(i-1)*self%dpsi
              call spl2d_evaln(self%r,zpsi,ztheta,1,zr)
              call spl2d_evaln(self%z,zpsi,ztheta,2,zz)
              posang%pos(1)=zr ; posang%pos(2)=zz ; posang%pos(3)=zeta
              posang%opt=1 ; posang%units=0
              ! posn to cartesians and output
              call posang_tfm(posang,-3)
              call posang_writev(posang,kplot)
           end do
        end do
     end do

     write(kplot,'(''POINT_DATA '',I8)') (self%n%npsi+1)*intheta*self%n%mzetav
     write(kplot,'(''VECTORS Bcart float'')')
     do k=1,self%n%mzetav
        zeta=(k-1)*dzeta
        zsin=sin(zeta)
        zcos=cos(zeta)
        do j=self%ntmin,self%ntmax
           ztheta=self%n%thetamin+(j-1)*self%dtheta
           do i=1,self%n%npsi+1
              zpsi=self%n%psimin +(i-1)*self%dpsi
              call spl2d_evaln(self%r,zpsi,ztheta,1,zr)
              call spl2d_evaln(self%z,zpsi,ztheta,2,zz)
              call spl2d_evaln(self%dpsidr,zr,zz,1,zdpdr)
              call spl2d_evaln(self%dpsidz,zr,zz,2,zdpdz)
              posang%pos(1)=zr ; posang%pos(2)=zz ; posang%pos(3)=zeta
              posang%units=0
              ! evaluate I aka f at psi
              call spleval(self%f,self%mr,self%psiaxis,self%psiqbdry,zpsi,zf,1)
              ! B in toroidal-cyl polars
              posang%vec(1)=-zdpdz/zr ; posang%vec(2)=zdpdr/zr; posang%vec(3)=zf/zr
              posang%opt=17
              ! B to cartesians and output
              call posang_tfm(posang,-3)
              call posang_writev(posang,kplot,2)
           end do
        end do
     end do

  case('allcartesian')
     call log_error(m_name,s_name,1,log_info,'Plot parameter')
     call log_value("Number of toroidal points",self%n%mzetav)

     write(kplot,'(''DATASET STRUCTURED_GRID'')')
     write(kplot,'(''DIMENSIONS '',I8,I8,I8)') self%mr+1,self%mz+1,self%n%mzetav
     write(kplot,'(''POINTS '',I8,'' float'')') (self%mr+1)*(self%mz+1)*self%n%mzetav

     do k=1,self%n%mzetav
        zeta=(k-1)*dzeta
        zsin=sin(zeta)
        zcos=cos(zeta)
        do j=1,self%mz+1
           zz=self%zmin+(j-1)*self%dz
           do i=1,self%mr+1
              zr=self%rmin +(i-1)*self%dr
              posang%pos(1)=zr ; posang%pos(2)=zz ; posang%pos(3)=zeta
              posang%opt=1 ; posang%units=0
              ! posn to cartesians and output
              call posang_tfm(posang,-3)
              call posang_writev(posang,kplot)
           end do
        end do
     end do

     write(kplot,'(''POINT_DATA '',I8)') (self%mr+1)*(self%mz+1)*self%n%mzetav
     write(kplot,'(''SCALARS Psi float'')')
     write(kplot,'(''LOOKUP_TABLE default'')')
     do k=1,self%n%mzetav
        call spl2d_writev(self%psi,'sampl',kplot)
     end do

     write(kplot,'(''VECTORS Bcart float'')')
     do k=1,self%n%mzetav
        zeta=(k-1)*dzeta
        do j=1,self%mz+1
           zz=self%zmin+(j-1)*self%dz
           do i=1,self%mr+1
              zr=self%rmin +(i-1)*self%dr
              posang%pos(1)=zr ; posang%pos(2)=zz ; posang%pos(3)=zeta
              posang%opt=1 ; posang%units=0
              ! posn to cartesians
              call posang_tfm(posang,-3)
              call beq_b(self,posang,0)
              call posang_writev(posang,kplot,2)
           end do
        end do
     end do

  case('cubeallcartesian')
     call log_error(m_name,s_name,1,log_info,'Visualisation output in duct coordinates')

     write(kplot,'(''DATASET RECTILINEAR_GRID'')')
     write(kplot,'(''DIMENSIONS '',I8,I8,I8)') &
 &   self%fmesh%nxf, self%fmesh%nyf, self%fmesh%nzf
     ! need to convert to mm
     if (self%fmesh%punit(1:2)/='me') then
        plotfac=1000.
     else
        plotfac=1.
     end if
     write(kplot,'(''X_COORDINATES '',I8,'' float'')') self%fmesh%nxf
     write(kplot,cfmtbs1) (plotfac*self%fmesh%xf(i), i=1,self%fmesh%nxf)
     write(kplot,'(''Y_COORDINATES '',I8,'' float'')') self%fmesh%nyf
     write(kplot,cfmtbs1) (plotfac*self%fmesh%yf(i), i=1,self%fmesh%nyf)
     write(kplot,'(''Z_COORDINATES '',I8,'' float'')') self%fmesh%nzf
     write(kplot,cfmtbs1) (plotfac*self%fmesh%zf(i), i=1,self%fmesh%nzf)

     write(kplot,'(''POINT_DATA '',I8)') &
 &   self%fmesh%nxf*self%fmesh%nyf*self%fmesh%nzf
     write(kplot,'(''VECTORS Bcart float'')')

     call beq_writen(self,'onlyvalues',kplot)

  case('psi-theta')
     call log_error(m_name,s_name,1,log_info,'Plot parameter')
     call log_value("Number of toroidal points",self%n%mzetav)

     intheta=self%ntmax-self%ntmin+1
     ztmin=self%n%thetamin+(self%ntmin-1)*self%dtheta
     write(kplot,'(''DATASET STRUCTURED_POINTS'')')
     write(kplot,'(''DIMENSIONS '',I8,I8,I8)') self%n%npsi+1,intheta,self%n%mzetav
     write(kplot,'(''ORIGIN '','//cfmt1v) self%n%psimin,ztmin,-const_pid
     write(kplot,'(''SPACING '','//cfmt1v) self%dpsi,self%dtheta,dzeta
     write(kplot,'(''POINT_DATA '',I8)') (self%n%npsi+1)*intheta*self%n%mzetav
     write(kplot,'(''VECTORS Bpt float'')')

     zc1=0
     do k=1,self%n%mzetav
        zeta=(k-1)*dzeta-const_pid ! not used
        do j=self%ntmin,self%ntmax
           ztheta=self%n%thetamin+(j-1)*self%dtheta
           do i=1,self%n%npsi+1
              zpsi=self%n%psimin +(i-1)*self%dpsi
              zc1=0._kr8
              call spl2d_eval(self%rjac,zpsi,ztheta,zc2)
              ! evaluate I aka f at psi
              call spleval(self%f,self%mr,self%psiaxis,self%psiqbdry,zpsi,zc3,1)
              write(kplot,cfmtb1v) zc1,zc2,zc3
           end do
        end do
     end do


  case('R-Z')

  case default

  end select plot_type

end subroutine beq_writev
!---------------------------------------------------------------------
!> calculate diagnostics and output
subroutine beq_dia(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_dia' !< subroutine name
  integer(ki4) :: intheta    !< actual number of angles defining  \f$ R,Z(\psi,\theta) \f$
  real(kr8) :: zpsi    !<  \f$ \psi_i \f$
  real(kr8) :: ztheta    !<  \f$ \theta_j \f$
  real(kr8) :: zr    !<   \f$ R(\psi_i,\theta_j) \f$
  real(kr8) :: zz    !<   \f$ Z(\psi_i,\theta_j) \f$
  real(kr8) :: zpsin    !<  recalculated \f$ \psi_i \f$
  real(kr8) :: zthetan    !<  recalculated \f$ \theta_j \f$
  real(kr8) :: zpsidiff    !< zpsin-zpsi
  real(kr8) :: zthetadiff    !< zthetan-ztheta
  real(kr8) :: zperr    !< maximum value of abs(zpsin-zpsi)
  real(kr8) :: zterr    !< maximum value of abs(zthetan-ztheta)

  !! allocate work2 storage
  intheta=self%ntmax-self%ntmin+1
  allocate(work2(self%n%npsi+1,intheta), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  zperr=0
  zterr=0
  ij=0
  do j=self%ntmin,self%ntmax
     ij=ij+1
     ztheta=self%n%thetamin+(j-1)*self%dtheta
     do i=1,self%n%npsi+1
        zpsi=self%n%psimin +(i-1)*self%dpsi
        call spl2d_evaln(self%r,zpsi,ztheta,1,zr)
        call spl2d_evaln(self%z,zpsi,ztheta,2,zz)
        call spl2d_eval(self%psi,zr,zz,zpsin)
        zpsidiff=zpsin-zpsi
        zthetan=atan2(zz-self%n%zcen,zr-self%n%rcen)
        if (zthetan<-const_pid/2) zthetan=2*const_pid+zthetan
        zthetadiff=zthetan-ztheta
        zperr=max( zperr,abs(zpsidiff) )
        zterr=max( zterr,abs(zthetadiff) )
        work2(i,ij)=zpsidiff
     end do
  end do

  call log_error(m_name,s_name,1,log_info,'Map data')
  call log_value("MAP LIMIT psimin ",self%n%psimin)
  call log_value("MAP LIMIT psimax ",self%n%psimax)
  ! call log_value("MAP LIMIT thetamin ",self%n%thetamin)
  ! call log_value("MAP LIMIT thetamax ",self%n%thetamax)
  call log_value("MAP LIMIT thetamin ",self%n%thetamin+(self%ntmin-1)*self%dtheta)
  call log_value("MAP LIMIT thetamax ",self%n%thetamax+(self%ntmax-1)*self%dtheta)
  call log_value("maximum error in psi ",zperr)
  call log_value("maximum error in theta ",zterr)
  !dbg      write(*,*) 'psi error array = ', work2

  deallocate(work2)

end subroutine beq_dia
!---------------------------------------------------------------------
!> write in data format
subroutine beq_write(self,kout)

  !! arguments
  type(beq_t), intent(in) :: self   !< object data structure
  integer, intent(in) :: kout   !< output channel for object data structure


  !! local
  character(*), parameter :: s_name='beq_write' !< subroutine name

  write(kout,*,iostat=status) 'rmin'
  write(kout,*,iostat=status) self%rmin
  if(status/=0) then
     call log_error(m_name,s_name,2,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'rmax'
  write(kout,*,iostat=status) self%rmax
  if(status/=0) then
     call log_error(m_name,s_name,3,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'mr'
  write(kout,*,iostat=status) self%mr
  if(status/=0) then
     call log_error(m_name,s_name,4,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'dr'
  write(kout,*,iostat=status) self%dr
  if(status/=0) then
     call log_error(m_name,s_name,5,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'zmin'
  write(kout,*,iostat=status) self%zmin
  if(status/=0) then
     call log_error(m_name,s_name,6,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'zmax'
  write(kout,*,iostat=status) self%zmax
  if(status/=0) then
     call log_error(m_name,s_name,7,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'mz'
  write(kout,*,iostat=status) self%mz
  if(status/=0) then
     call log_error(m_name,s_name,8,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'dz'
  write(kout,*,iostat=status) self%dz
  if(status/=0) then
     call log_error(m_name,s_name,9,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'rcen'
  write(kout,*,iostat=status) self%n%rcen
  if(status/=0) then
     call log_error(m_name,s_name,17,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'zcen'
  write(kout,*,iostat=status) self%n%zcen
  if(status/=0) then
     call log_error(m_name,s_name,18,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'psimin'
  write(kout,*,iostat=status) self%n%psimin
  if(status/=0) then
     call log_error(m_name,s_name,19,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'psimax'
  write(kout,*,iostat=status) self%n%psimax
  if(status/=0) then
     call log_error(m_name,s_name,20,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'npsi'
  write(kout,*,iostat=status) self%n%npsi
  if(status/=0) then
     call log_error(m_name,s_name,21,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'dpsi'
  write(kout,*,iostat=status) self%dpsi
  if(status/=0) then
     call log_error(m_name,s_name,22,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'thetamin'
  write(kout,*,iostat=status) self%n%thetamin
  if(status/=0) then
     call log_error(m_name,s_name,23,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'thetamax'
  write(kout,*,iostat=status) self%n%thetamax
  if(status/=0) then
     call log_error(m_name,s_name,24,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'ntheta'
  write(kout,*,iostat=status) self%n%ntheta
  if(status/=0) then
     call log_error(m_name,s_name,25,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'ntmin'
  write(kout,*,iostat=status) self%ntmin
  if(status/=0) then
     call log_error(m_name,s_name,26,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'ntmax'
  write(kout,*,iostat=status) self%ntmax
  if(status/=0) then
     call log_error(m_name,s_name,27,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'cenopt'
  write(kout,*,iostat=status) self%n%cenopt
  if(status/=0) then
     call log_error(m_name,s_name,28,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'psiopt'
  write(kout,*,iostat=status) self%n%psiopt
  if(status/=0) then
     call log_error(m_name,s_name,29,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'bdryopt'
  write(kout,*,iostat=status) self%n%bdryopt
  if(status/=0) then
     call log_error(m_name,s_name,39,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'nopt'
  write(kout,*,iostat=status) self%n%nopt
  if(status/=0) then
     call log_error(m_name,s_name,30,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'thetaopt'
  write(kout,*,iostat=status) self%n%thetaopt
  if(status/=0) then
     call log_error(m_name,s_name,31,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'dtheta'
  write(kout,*,iostat=status) self%dtheta
  if(status/=0) then
     call log_error(m_name,s_name,32,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'psiref'
  write(kout,*,iostat=status) self%n%psiref
  if(status/=0) then
     call log_error(m_name,s_name,33,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'thetaref'
  write(kout,*,iostat=status) self%n%thetaref
  if(status/=0) then
     call log_error(m_name,s_name,34,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'psibdry'
  write(kout,*,iostat=status) self%psibdry
  if(status/=0) then
     call log_error(m_name,s_name,35,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'psiqbdry'
  write(kout,*,iostat=status) self%psiqbdry
  if(status/=0) then
     call log_error(m_name,s_name,36,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'rbdry'
  write(kout,*,iostat=status) self%rbdry
  if(status/=0) then
     call log_error(m_name,s_name,37,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'bpbdry'
  write(kout,*,iostat=status) self%bpbdry
  if(status/=0) then
     call log_error(m_name,s_name,38,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'btotbdry'
  write(kout,*,iostat=status) self%btotbdry
  if(status/=0) then
     call log_error(m_name,s_name,39,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'f'
  write(kout,*,iostat=status) self%f
  if(status/=0) then
     call log_error(m_name,s_name,10,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'i'
  write(kout,*,iostat=status) self%i
  if(status/=0) then
     call log_error(m_name,s_name,11,error_fatal,'Error writing object data')
  end if
  write(kout,*,iostat=status) 'psi --- 2d spline'
  call spl2d_write( self%psi,kout )
  write(kout,*,iostat=status) 'end of 2d spline'
  write(kout,*,iostat=status) 'dpsidr --- 2d spline'
  call spl2d_write( self%dpsidr,kout )
  write(kout,*,iostat=status) 'end of 2d spline'
  write(kout,*,iostat=status) 'dpsidz --- 2d spline'
  call spl2d_write( self%dpsidz,kout )
  write(kout,*,iostat=status) 'end of 2d spline'
  write(kout,*,iostat=status) 'r --- 2d spline'
  call spl2d_write( self%r,kout )
  write(kout,*,iostat=status) 'end of 2d spline'
  write(kout,*,iostat=status) 'z --- 2d spline'
  call spl2d_write( self%z,kout )
  write(kout,*,iostat=status) 'end of 2d spline'
  write(kout,*,iostat=status) 'rjac --- 2d spline'
  call spl2d_write( self%rjac,kout )
  write(kout,*,iostat=status) 'end of 2d spline'

  write(kout,*,iostat=status) 'srmin'
  write(kout,*,iostat=status) self%srmin
  call log_write_check(m_name,s_name,40,status)
  write(kout,*,iostat=status) 'srmax'
  write(kout,*,iostat=status) self%srmax
  call log_write_check(m_name,s_name,41,status)

end subroutine beq_write
!---------------------------------------------------------------------
!> write field data needed by ITER
subroutine beq_writepart(self,kout)

  !! arguments
  type(beq_t), intent(in) :: self   !< object data structure
  integer, intent(in) :: kout   !< output channel for object data structure


  !! local
  character(*), parameter :: s_name='beq_writepart' !< subroutine name

  write(kout,*,iostat=status) 'mr'
  write(kout,*,iostat=status) self%mr
  call log_write_check(m_name,s_name,4,status)
  write(kout,*,iostat=status) 'rcen'
  write(kout,*,iostat=status) self%n%rcen
  call log_write_check(m_name,s_name,17,status)
  write(kout,*,iostat=status) 'zcen'
  write(kout,*,iostat=status) self%n%zcen
  call log_write_check(m_name,s_name,18,status)
  write(kout,*,iostat=status) 'psiaxis'
  write(kout,*,iostat=status) self%psiaxis
  call log_write_check(m_name,s_name,34,status)
  write(kout,*,iostat=status) 'psibdry'
  write(kout,*,iostat=status) self%psibdry
  call log_write_check(m_name,s_name,35,status)
  write(kout,*,iostat=status) 'psiqbdry'
  write(kout,*,iostat=status) self%psiqbdry
  call log_write_check(m_name,s_name,36,status)
  write(kout,*,iostat=status) 'rbdry'
  write(kout,*,iostat=status) self%rbdry
  call log_write_check(m_name,s_name,37,status)
  write(kout,*,iostat=status) 'bpbdry'
  write(kout,*,iostat=status) self%bpbdry
  call log_write_check(m_name,s_name,38,status)
  write(kout,*,iostat=status) 'btotbdry'
  write(kout,*,iostat=status) self%btotbdry
  call log_write_check(m_name,s_name,39,status)
  write(kout,*,iostat=status) 'f'
  write(kout,*,iostat=status) self%f
  call log_write_check(m_name,s_name,40,status)

  write(kout,*,iostat=status) 'r --- 2d spline'
  call spl2d_write( self%r,kout )
  write(kout,*,iostat=status) 'end of 2d spline'
  write(kout,*,iostat=status) 'z --- 2d spline'
  call spl2d_write( self%z,kout )
  write(kout,*,iostat=status) 'end of 2d spline'
  write(kout,*,iostat=status) 'rjac --- 2d spline'
  call spl2d_write( self%rjac,kout )
  write(kout,*,iostat=status) 'end of 2d spline'

  ! added 'plus parts'
  write(kout,*,iostat=status) 'fldspec'
  write(kout,*,iostat=status) self%n%fldspec
  call log_write_check(m_name,s_name,4,status)
  write(kout,*,iostat=status) 'vacfile'
  write(kout,*,iostat=status) self%n%vacfile
  call log_write_check(m_name,s_name,63,status)
  ! ripple field data
  write(kout,*,iostat=status) 'mrip'
  write(kout,*,iostat=status) self%n%mrip
  call log_write_check(m_name,s_name,64,status)
  write(kout,*,iostat=status) 'ivac'
  write(kout,*,iostat=status) self%ivac
  call log_write_check(m_name,s_name,65,status)
  if (self%n%mrip/=0) then
     write(kout,*,iostat=status) 'arip'
     write(kout,*,iostat=status) self%n%arip
     call log_write_check(m_name,s_name,66,status)
  end if


end subroutine beq_writepart
!---------------------------------------------------------------------
!> delete field data used by ITER
subroutine beq_deletepart(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='beq_deletepart' !< subroutine name
  !-----------------------------------------------------------------------
  !             Delete 2D splines
  call spl2d_delete( self%r )
  call spl2d_delete( self%z )
  call spl2d_delete( self%rjac )

end subroutine beq_deletepart
!---------------------------------------------------------------------
!> write field data needed by 'msus','global'
subroutine beq_writeplus(self,kout)

  !! arguments
  type(beq_t), intent(in) :: self   !< object data structure
  integer, intent(in) :: kout   !< output channel for object data structure

  !! local
  character(*), parameter :: s_name='beq_writeplus' !< subroutine name
  integer(ki4) :: ifldspec !< field as mapped or 3-cpt + extra info
  integer(ki4) :: iextra !< extra info extracted

  write(kout,*,iostat=status) 'fldspec'
  write(kout,*,iostat=status) self%n%fldspec
  call log_write_check(m_name,s_name,4,status)

  write(kout,*,iostat=status) 'mr'
  write(kout,*,iostat=status) self%mr
  call log_write_check(m_name,s_name,5,status)

  write(kout,*,iostat=status) 'rcen'
  write(kout,*,iostat=status) self%n%rcen
  call log_write_check(m_name,s_name,17,status)
  write(kout,*,iostat=status) 'zcen'
  write(kout,*,iostat=status) self%n%zcen
  call log_write_check(m_name,s_name,18,status)
  write(kout,*,iostat=status) 'psiaxis'
  write(kout,*,iostat=status) self%psiaxis
  call log_write_check(m_name,s_name,34,status)
  write(kout,*,iostat=status) 'psibdry'
  write(kout,*,iostat=status) self%psibdry
  call log_write_check(m_name,s_name,35,status)
  write(kout,*,iostat=status) 'psiqbdry'
  write(kout,*,iostat=status) self%psiqbdry
  call log_write_check(m_name,s_name,36,status)
  write(kout,*,iostat=status) 'rbdry'
  write(kout,*,iostat=status) self%rbdry
  call log_write_check(m_name,s_name,37,status)
  write(kout,*,iostat=status) 'bpbdry'
  write(kout,*,iostat=status) self%bpbdry
  call log_write_check(m_name,s_name,38,status)
  write(kout,*,iostat=status) 'btotbdry'
  write(kout,*,iostat=status) self%btotbdry
  call log_write_check(m_name,s_name,39,status)
  write(kout,*,iostat=status) 'zetamin'
  write(kout,*,iostat=status) self%n%zetamin
  call log_write_check(m_name,s_name,40,status)
  write(kout,*,iostat=status) 'zetamax'
  write(kout,*,iostat=status) self%n%zetamax
  call log_write_check(m_name,s_name,41,status)
  write(kout,*,iostat=status) 'nzets'
  write(kout,*,iostat=status) self%nzets
  call log_write_check(m_name,s_name,42,status)
  write(kout,*,iostat=status) 'ximin'
  write(kout,*,iostat=status) self%n%ximin
  call log_write_check(m_name,s_name,43,status)
  write(kout,*,iostat=status) 'ximax'
  write(kout,*,iostat=status) self%n%ximax
  call log_write_check(m_name,s_name,44,status)

  write(kout,*,iostat=status) 'f'
  write(kout,*,iostat=status) self%f
  call log_write_check(m_name,s_name,50,status)

  write(kout,*,iostat=status) 'psi --- 2d spline'
  call spl2d_write( self%psi,kout )
  write(kout,*,iostat=status) 'end of 2d spline'
  write(kout,*,iostat=status) 'dpsidr --- 2d spline'
  call spl2d_write( self%dpsidr,kout )
  write(kout,*,iostat=status) 'end of 2d spline'
  write(kout,*,iostat=status) 'dpsidz --- 2d spline'
  call spl2d_write( self%dpsidz,kout )
  write(kout,*,iostat=status) 'end of 2d spline'

  ifldspec=self%n%fldspec
  iextra=ifldspec/10
  ifldspec=self%n%fldspec-10*iextra
  fld_specn: select case (ifldspec)
  case(3)
     write(kout,*,iostat=status) 'rispldr --- 2d spline'
     call spl2d_write( self%rispldr,kout )
     write(kout,*,iostat=status) 'end of 2d spline'
     write(kout,*,iostat=status) 'rispldz --- 2d spline'
     call spl2d_write( self%rispldz,kout )
     write(kout,*,iostat=status) 'end of 2d spline'
  end select fld_specn

  ! ripple field data
  write(kout,*,iostat=status) 'mrip'
  write(kout,*,iostat=status) self%n%mrip
  call log_write_check(m_name,s_name,60,status)
  write(kout,*,iostat=status) 'ivac'
  write(kout,*,iostat=status) self%ivac
  call log_write_check(m_name,s_name,61,status)
  if (self%n%mrip/=0) then
     write(kout,*,iostat=status) 'arip'
     write(kout,*,iostat=status) self%n%arip
     call log_write_check(m_name,s_name,62,status)
  end if
  write(kout,*,iostat=status) 'vacfile'
  write(kout,*,iostat=status) self%n%vacfile
  call log_write_check(m_name,s_name,63,status)

  if (self%n%skylpsi) then
     write(kout,*,iostat=status) 'rxpt'
     write(kout,*,iostat=status) self%rxpt
     call log_write_check(m_name,s_name,64,status)
     write(kout,*,iostat=status) 'zxpt'
     write(kout,*,iostat=status) self%zxpt
     call log_write_check(m_name,s_name,65,status)
     write(kout,*,iostat=status) 'psixpt'
     write(kout,*,iostat=status) self%psixpt
     call log_write_check(m_name,s_name,66,status)
     write(kout,*,iostat=status) 'zmin'
     write(kout,*,iostat=status) self%zmin
     call log_write_check(m_name,s_name,67,status)
     write(kout,*,iostat=status) 'zmax'
     write(kout,*,iostat=status) self%zmax
     call log_write_check(m_name,s_name,68,status)
     write(kout,*,iostat=status) 'mz'
     write(kout,*,iostat=status) self%mz
     call log_write_check(m_name,s_name,69,status)
     write(kout,*,iostat=status) 'dz'
     write(kout,*,iostat=status) self%dz
     call log_write_check(m_name,s_name,70,status)
     write(kout,*,iostat=status) 'nctrack'
     call log_write_check(m_name,s_name,71,status)
     write(kout,*,iostat=status) self%nctrack
     call log_write_check(m_name,s_name,72,status)
     write(kout,*,iostat=status) 'ctrackrz'
     call log_write_check(m_name,s_name,73,status)
     write(kout,*,iostat=status) ((self%ctrackrz(i,j),i=1,self%nctrack),j=1,2)
     call log_write_check(m_name,s_name,74,status)
  end if

  if (self%n%duct) then
     call fmesh_write(self%fmesh,kout)
  end if

end subroutine beq_writeplus
!---------------------------------------------------------------------
!> delete field data used by 'msus','global'
subroutine beq_deleteplus(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='beq_deleteplus' !< subroutine name
  !-----------------------------------------------------------------------
  !             Delete 2D splines
  call spl2d_delete( self%psi )
  call spl2d_delete( self%dpsidr )
  call spl2d_delete( self%dpsidz )

  if (self%n%skylpsi) then
     deallocate( self%ctrackrz )
  end if

end subroutine beq_deleteplus
!---------------------------------------------------------------------
!> initialise spline data structure from beq
subroutine beq_spl2d(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_spl2d' !< subroutine name

  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading vec')
  end if

end subroutine beq_spl2d
!---------------------------------------------------------------------
!> calculate \f$ R_c, Z_c\f$
subroutine beq_centre(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_centre' !< subroutine name
  integer(ki4) :: igr !< \f$ R \f$ position index of global extremum
  integer(ki4) :: igz !< \f$ Z \f$ position index of global extremum
  integer(ki4) :: soutr=10 !< maximum number of outer searches
  integer(ki4) :: sinr=10 !< maximum number of inner searches
  integer(ki4) :: isr !< \f$ R \f$ search (increment) direction \f$ \pm 1 \f$
  integer(ki4) :: isz !< \f$ Z \f$ search (increment) direction \f$ \pm 1 \f$
  integer(ki4) :: ir !< current \f$ R \f$ position as index
  integer(ki4) :: iz !< current \f$ Z \f$ position as index
  real(kr8) :: zpp !< current value of \f$ \psi \f$
  real(kr8) :: zpg !< value of \f$ \psi \f$ at global extremum

  ! starting estimate of array index (igr,igz) where psi least
  start_centre : select case (self%n%cenopt)
  case(1) ! set by user, do nothing
     return
  case(2) ! use central index as guess
     igr=(self%rmax-self%rmin)/(2*self%dr)
     igz=(self%zmax-self%zmin)/(2*self%dz)
  case(3) ! use user hint
     igr=(self%n%rcen-self%rmin)/self%dr
     igz=(self%n%zcen-self%zmin)/self%dz
     ! fix up in case
     igr=min( max(igr,1) , self%mr )
     igz=min( max(igz,1) , self%mz )
  case(4) ! use value from EQDSK file
     self%n%rcen=self%rqcen
     self%n%zcen=self%zqcen
     call spl2d_eval(self%psi,self%rqcen,self%zqcen,zpp)
     goto 1
  end select start_centre

  ! set initially positive search  directions for increasing psi
  isr=1
  isz=1

  do j=1,soutr
     ir=igr
     iz=igz
     call spl2d_eval(self%psi,self%rmin+(ir-1)*self%dr,self%zmin+(iz-1)*self%dz,zpg)

     !! search in r
     ! change direction if necessary
     ir=igr+isr
     call spl2d_eval(self%psi,self%rmin+(ir-1)*self%dr,self%zmin+(iz-1)*self%dz,zpp)
     if ( (zpp-zpg)*rsig < 0 ) then
        igr=ir
        zpg=zpp
     else
        isr=-1
        ir=igr
     end if

     do i=1,sinr
        ir=ir+isr
        call spl2d_eval(self%psi,self%rmin+(ir-1)*self%dr,self%zmin+(iz-1)*self%dz,zpp)
        if ( (zpp-zpg)*rsig < 0 ) then
           igr=ir
           zpg=zpp
        else
           exit
        end if
     end do

     !! search in \f$ Z \f$
     ! change direction if necessary
     iz=igz+isz
     call spl2d_eval(self%psi,self%rmin+(ir-1)*self%dr,self%zmin+(iz-1)*self%dz,zpp)
     if ( (zpp-zpg)*rsig < 0 ) then
        igz=iz
        zpg=zpp
     else
        isz=-1
        iz=igz
     end if

     do i=1,sinr
        iz=iz+isz
        call spl2d_eval(self%psi,self%rmin+(ir-1)*self%dr,self%zmin+(iz-1)*self%dz,zpp)
        if ( (zpp-zpg)*rsig < 0 ) then
           igz=iz
           zpg=zpp
        else
           exit
        end if
     end do

     ! check global minimum
     call spl2d_eval(self%psi,self%rmin+(igr-2)*self%dr,self%zmin+(igz-1)*self%dz,zpp)
     if ( (zpp-zpg)*rsig < 0 ) cycle
     call spl2d_eval(self%psi,self%rmin+igr*self%dr,self%zmin+(igz-1)*self%dz,zpp)
     if ( (zpp-zpg)*rsig < 0 ) cycle
     call spl2d_eval(self%psi,self%rmin+(igr-1)*self%dr,self%zmin+(igz-2)*self%dz,zpp)
     if ( (zpp-zpg)*rsig < 0 ) cycle
     call spl2d_eval(self%psi,self%rmin+(igr-1)*self%dr,self%zmin+igz*self%dz,zpp)
     if ( (zpp-zpg)*rsig < 0 ) cycle
     exit

  end do

  self%n%rcen=self%rmin+(igr-1)*self%dr
  self%n%zcen=self%zmin+(igz-1)*self%dz

1     continue
  call log_error(m_name,s_name,1,log_info,'Equil. centre')
  call log_value("SMITER-GEOQ rcen ",self%n%rcen)
  call log_value("SMITER-GEOQ zcen ",self%n%zcen)
  self%psicen=zpp
  call log_value("SMITER-GEOQ psicen ",self%psicen)
  call log_value("SMITER-GEOQ psiaxis ",self%psiaxis)
  if (abs(self%psicen-self%psiaxis)> 0.1*self%psinorm) then
     call log_error(m_name,s_name,2,error_warning,'Eqdsk value of psi on axis very different from computed')
     call log_error(m_name,s_name,2,error_warning,'Using geoq computed value psicen')
     self%replasi=.TRUE.
     self%psiaxis=self%psicen
     call log_error(m_name,s_name,2,error_warning,'plasma boundary psi may change')
     boundary_type: select case (self%n%bdryopt)
     case(2,6,7,10,11,12)
        call log_error(m_name,s_name,3,error_warning,'beq_bdryopt may not work as expected')
     end select boundary_type
  end if

end subroutine beq_centre
!---------------------------------------------------------------------
!> namelist input of controls
subroutine beq_readcon(selfn,kin)

  !! arguments
  type(bnumerics_t), intent(out) :: selfn   !< object control data structure
  integer, intent(in) :: kin   !< input channel for object data structure


  !! local
  character(*), parameter :: s_name='beq_readcon' !< subroutine name
  integer(ki4) :: mzeta_vtk !< number of sample points in  \f$ \zeta \f$ for vtk
  integer(ki4) :: mzeta_gnup !< number of sample points in  \f$ \zeta \f$ for gnuplot
  integer(ki4) :: beq_cenopt !< namelist option for \f$ R_{cen}, Z_{cen} \f$
  integer(ki4) :: beq_psiopt !< namelist option for \f$ \psi{\min}, \psi_{\max}\f$
  integer(ki4) :: beq_bdryopt !< namelist option for setting reference boundary \f$ \psi \f$
 &! (reference boundary may be user given, X-point, limiter etc.)
  integer(ki4) :: beq_nopt !< namelist option for \f$ N_{\psi}, N_{\theta} \f$
  integer(ki4) :: beq_thetaopt !< namelist option for \f$ \theta_{\min}, \theta_{\max} \f$
  integer(ki4) :: beq_zetaopt !< namelist option for \f$ \zeta_{\min}, \zeta_{\max} \f$
  integer(ki4) :: beq_xiopt !< namelist option for \f$ \xi_{\min}, \xi_{\max} \f$

  real(kr8) :: beq_rcen !< namelist \f$ R_{cen} \f$
  real(kr8) :: beq_zcen !< namelist \f$ Z_{cen} \f$

  real(kr8) :: beq_psimin !< namelist \f$ \psi_{\min} \f$ is numerically smaller than \f$ \psi_{\max} \f$
  real(kr8) :: beq_psimax !< namelist \f$ \psi_{\max} > \psi_{\min}  \f$
  integer(ki4) :: beq_npsi !< namelist \f$ N_{\psi} \f$
  real(kr8) :: beq_delpsi !< namelist \f$ \Delta\psi \f$

  real(kr8) :: beq_thetamin !< namelist \f$ \theta_{\min} \f$
  real(kr8) :: beq_thetamax !< namelist \f$ \theta_{\max} \f$
  real(kr8) :: beq_deltheta !< namelist \f$ \Delta\theta \f$
  real(kr8) :: beq_zetamin !< namelist \f$ \zeta_{\min} \f$
  real(kr8) :: beq_zetamax !< namelist \f$ \zeta_{\max} \f$
  integer(ki4) :: beq_nzetap !< namelist \f$ N_{\zeta P} \f$
  logical :: skip_eqdsk_b !< read B cpts from eqdsk file unless .TRUE.
  real(kr8) :: beq_rmove !< namelist \f$ R_{mov} \f$ equilibrium displaced
  real(kr8) :: beq_zmove !< namelist \f$ Z_{mov} \f$ equilibrium displaced
  real(kr8) :: beq_fscale !< namelist \f$ f \f$ equilibrium scaled
  integer(ki4) :: beq_mrip !< \f$ N \f$ number of ripple coils
  real(kr8) :: beq_irip !< unused parameter for ripple coils
  logical :: beq_duct !< flag whether working in duct coordinates
  logical :: skylight_flux_limits !< flag whether any skylight(s) defined by flux limits
  logical :: skylight_centre_line !< flag whether any skylight(s) defined by centre_line
  logical :: centre_line_cutout !< flag whether centre line is to be a cutout
  integer(ki4) :: absorber_objects !< number of absorbers defined geometrically
  integer(ki4) :: invisible_objects !< number of invisible objects defined geometrically
  integer(ki4) :: skylight_objects !< number of skylights defined geometrically
  integer(ki4) :: beancan_objects !< number of beancan surfaces defined geometrically
  integer(ki4) :: cutout_objects !< number of cutouts defined geometrically
  integer(ki4) :: skylight_debug !< level of output to help understand skylight production
  real(kr8) :: beq_arip !< \f$ a \f$ for ripple coils
  real(kr8) :: beq_psiref !< namelist \f$ \psi_X \f$
  real(kr8) :: beq_thetaref !< namelist \f$ \theta_X \f$
  integer(ki4) :: beq_ntheta !< namelist \f$ N_{\theta} \f$
  integer(ki4) :: beq_fldspec !< field specification
  integer(ki4) :: beq_xsearch !< how to search for X-point (1 = within box)
  integer(ki4) :: limiter_search !< how to search for contact point (1 = within box)
  real(kr8) :: beq_psibig !< unity implies \f$ 2 \pi \f$ too big
  logical :: equil_helicity_ok !<  equilibrium helicity is ok
  real(kr8) :: beq_xrsta !< X-point box min R (m)
  real(kr8) :: beq_xrend !< X-point box max R (m)
  real(kr8) :: beq_xzsta !< X-point box min Z (m)
  real(kr8) :: beq_xzend !< X-point box max Z (m)
  real(kr8) :: search_r_start !< search box min R (m)
  real(kr8) :: search_r_end !< search box max R (m)
  real(kr8) :: search_z_start !< search box min Z (m)
  real(kr8) :: search_z_end !< search box max Z (m)
  real(kr8) :: z1 !< scratch
  real(kr8) :: z2 !< scratch
  character(len=80) :: beq_vacuum_field_file !< local variable

  !! beq parameters
  namelist /beqparameters/ &
 &beq_cenopt, beq_psiopt, &
 &mzeta_vtk, mzeta_gnup, &
 &beq_bdryopt, beq_nopt, &
 &beq_thetaopt, &
 &beq_zetaopt, beq_xiopt, &
 &beq_rcen, beq_zcen, &
 &beq_psimin, beq_psimax, beq_npsi, beq_delpsi,&
 &beq_thetamin, beq_thetamax, beq_deltheta,&
 &beq_rmove, beq_zmove,&
 &beq_fscale,&
 &beq_mrip, beq_irip, beq_arip,&
 &beq_duct,&
 &skylight_objects,&
 &absorber_objects,&
 &invisible_objects,&
 &beancan_objects,&
 &cutout_objects,&
 &skylight_flux_limits,&
 &skylight_centre_line,&
 &centre_line_cutout,&
 &skylight_debug,&
 &beq_psiref,&
 &beq_zetamin, beq_zetamax, beq_nzetap,&
 &skip_eqdsk_b,&
 &beq_thetaref, beq_ntheta, beq_fldspec,&
 &beq_psibig,&
 &equil_helicity_ok,&
 &beq_xsearch,beq_xrsta,beq_xrend,beq_xzsta,beq_xzend,&
 &limiter_search,search_r_start,search_r_end,search_z_start,search_z_end,&
 &beq_vacuum_field_file

  !! set default beq parameters
  mzeta_vtk=193
  mzeta_gnup=1
  ! default is user input
  beq_cenopt=1
  beq_psiopt=1
  beq_bdryopt=1
  beq_nopt=1
  beq_thetaopt=1
  beq_zetaopt=3
  beq_xiopt=2

  beq_rcen=1
  beq_zcen=0
  beq_psimin=0.9
  beq_psimax=1.1
  beq_delpsi=0.1
  beq_npsi=32
  beq_thetamin=const_pid/4
  beq_thetamax=3*const_pid/4
  beq_deltheta=0.1
  beq_zetamin=0
  beq_zetamax=2*const_pid
  beq_nzetap=0
  skip_eqdsk_b=.FALSE.
  beq_rmove=0.
  beq_zmove=0.
  beq_fscale=1
  beq_mrip=0
  beq_irip=0.
  beq_duct=.FALSE.
  absorber_objects=0
  invisible_objects=0
  skylight_objects=0
  beancan_objects=0
  cutout_objects=0
  skylight_flux_limits=.FALSE.
  skylight_centre_line=.FALSE.
  centre_line_cutout=.FALSE.
  skylight_debug=0
  beq_arip=0.
  beq_psiref=1
  beq_thetaref=0
  beq_ntheta=32
  beq_fldspec=1
  beq_psibig=0
  equil_helicity_ok=.FALSE.
  beq_xsearch=0
  beq_xrsta=0
  beq_xrend=0
  beq_xzsta=0
  beq_xzend=0
  limiter_search=0
  search_r_start=0
  search_r_end=0
  search_z_start=0
  search_z_end=0
  beq_vacuum_field_file='null'

  !!read beq parameters
  read(kin,nml=beqparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading beq parameters")'
     call log_getunit(ilog)
     write(ilog,nml=beqparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading beq parameters')
  end if


  !! check for valid data
  if(mzeta_vtk<=0) &
 &call log_error(m_name,s_name,2,error_fatal,'mzeta_vtk must be positive integer')
  if(mzeta_gnup<=0) &
 &call log_error(m_name,s_name,2,error_fatal,'mzeta_gnup must be positive integer')
  if(beq_cenopt<=0.OR.beq_cenopt>=5) &
 &call log_error(m_name,s_name,2,error_fatal,'beq_cenopt must be small positive integer')
  if(beq_cenopt==1.AND.beq_rcen<=0) &
 &call log_error(m_name,s_name,2,error_fatal,'beq_rcen must be > 0')
  if(beq_cenopt==1.AND.beq_zcen/=0) &
 &call log_error(m_name,s_name,2,error_warning,'beq_zcen non-zero')

  if(beq_psiopt<=0.OR.beq_psiopt>=4) &
 &call log_error(m_name,s_name,3,error_fatal,'beq_psiopt must be small positive integer')
  if(beq_psiopt==1.AND.beq_psimax-beq_psimin<=0) &
 &call log_error(m_name,s_name,3,error_fatal,'range of psi must be > 0')
  if(beq_psiopt==1.AND.beq_psimin<=0) &
 &call log_error(m_name,s_name,3,error_fatal,'beq_psimin must be > 0')

  if(beq_bdryopt<=0.OR.beq_bdryopt>=16) &
 &call log_error(m_name,s_name,6,error_fatal,'beq_bdryopt must be small positive integer')
  if(beq_nopt<=0.OR.beq_nopt>=4) &
 &call log_error(m_name,s_name,4,error_fatal,'beq_nopt must be small positive integer')
  if(beq_nopt==1.AND.beq_npsi*beq_ntheta<=1) &
 &call log_error(m_name,s_name,4,error_fatal,'dimensions of psi,theta space must be > 1')
  if(beq_nopt==1.AND.beq_npsi<=1) &
 &call log_error(m_name,s_name,4,error_fatal,'beq_npsi must be > 1')

  if(beq_thetaopt<=0.OR.beq_thetaopt>=4) &
 &call log_error(m_name,s_name,5,error_fatal,'beq_thetaopt must be small positive integer')
  if(beq_zetaopt<=0.OR.beq_zetaopt>=6) &
 &call log_error(m_name,s_name,5,error_fatal,'beq_zetaopt must be small positive integer')
  if(beq_xiopt<=0.OR.beq_xiopt>=6) &
 &call log_error(m_name,s_name,5,error_fatal,'beq_xiopt must be small positive integer')

  if(beq_thetaopt==1.AND.beq_thetamax-beq_thetamin<=0) &
 &call log_error(m_name,s_name,5,error_fatal,'range of theta must be > 0')
  if(beq_thetaopt==1.AND.beq_thetamin<=0) &
 &call log_error(m_name,s_name,5,error_fatal,'beq_thetamin must be > 0')
  if(beq_deltheta<0) &
 &call log_error(m_name,s_name,5,error_fatal,'theta margin must be >= 0')
  if(abs(beq_rmove)>100.) &
 &call log_error(m_name,s_name,5,error_warning,'rmove seems large')
  if(beq_rmove/=0.AND.abs(beq_rmove)<1.) &
 &call log_error(m_name,s_name,5,error_warning,'rmove seems small')
  if(abs(beq_zmove)>100.) &
 &call log_error(m_name,s_name,5,error_warning,'zmove seems large')
  if(beq_zmove/=0.AND.abs(beq_zmove)<1.) &
 &call log_error(m_name,s_name,5,error_warning,'zmove seems small')
  if(beq_zetaopt==1.AND.beq_zetamax-beq_zetamin<=0) &
 &call log_error(m_name,s_name,6,error_fatal,'range of zeta must be > 0')
  if(beq_zetaopt==1.AND.beq_zetamin<=0) &
 &call log_error(m_name,s_name,6,error_fatal,'beq_zetamin must be > 0')
  if(beq_nzetap<0) &
 &call log_error(m_name,s_name,6,error_fatal,'zeta scaling must be >= 0')
  if (beq_vacuum_field_file=='null') then
     if(beq_mrip<=0) &
 &   call log_error(m_name,s_name,7,error_warning,'negative or zero beq_mrip')
  end if
  if(beq_fldspec<=0.OR.beq_fldspec>=4) &
 &call log_error(m_name,s_name,7,error_fatal,'beq_fldspec must be small positive integer')
  if(beq_psibig<0.OR.beq_psibig>=2) &
 &call log_error(m_name,s_name,8,error_fatal,'beq_psibig must be small non-negative integer')
  if(.NOT.equil_helicity_ok) &
 &call log_error(m_name,s_name,8,error_warning,'Sense of helicity of field input may be changed - check log file')
  if(beq_xsearch<0.OR.beq_xsearch>=2) &
 &call log_error(m_name,s_name,9,error_fatal,'beq_xsearch must be small non-negative integer')
  if(limiter_search<0.OR.limiter_search>=2) &
 &call log_error(m_name,s_name,9,error_fatal,'limiter_search must be small non-negative integer')
  if(beq_xsearch==1) then
     z1=min(beq_xrsta,beq_xrend)
     z2=max(beq_xrsta,beq_xrend)
     if (abs(z1-z2)<1.e-3) &
 &   call log_error(m_name,s_name,10,error_fatal,'X-point search box too small in R')
     beq_xrsta=z1
     beq_xrend=z2
     z1=min(beq_xzsta,beq_xzend)
     z2=max(beq_xzsta,beq_xzend)
     if (abs(z1-z2)<1.e-3) &
 &   call log_error(m_name,s_name,11,error_fatal,'X-point search box too small in Z')
     beq_xzsta=z1
     beq_xzend=z2
  end if
  if(absorber_objects<0) &
 &call log_error(m_name,s_name,12,error_fatal,'absorber_objects must be non-negative integer')
  if(invisible_objects<0) &
 &call log_error(m_name,s_name,13,error_fatal,'invisible_objects must be non-negative integer')
  if(skylight_objects<0) &
 &call log_error(m_name,s_name,14,error_fatal,'skylight_objects must be non-negative integer')
  if(beancan_objects<0) &
 &call log_error(m_name,s_name,15,error_fatal,'beancan_objects must be non-negative integer')
  if(cutout_objects<0) &
 &call log_error(m_name,s_name,16,error_fatal,'cutout_objects must be non-negative integer')
  if(limiter_search==1) then
     z1=min(search_r_start,search_r_end)
     z2=max(search_r_start,search_r_end)
     if (abs(z1-z2)<1.e-3) &
 &   call log_error(m_name,s_name,20,error_fatal,'search box too small in R')
     search_r_start=z1
     search_r_end=z2
     z1=min(search_z_start,search_z_end)
     z2=max(search_z_start,search_z_end)
     if (abs(z1-z2)<1.e-3) &
 &   call log_error(m_name,s_name,20,error_fatal,'search box too small in Z')
     search_z_start=z1
     search_z_end=z2
  end if

  !! store values
  selfn%cenopt=beq_cenopt
  selfn%mzetav=mzeta_vtk
  selfn%mzetag=mzeta_gnup
  selfn%psiopt=beq_psiopt
  selfn%bdryopt=beq_bdryopt
  selfn%nopt=beq_nopt
  selfn%thetaopt=beq_thetaopt
  ! convert moves to metres
  selfn%rmove=.001_kr8*beq_rmove
  selfn%zmove=.001_kr8*beq_zmove
  ! scale f
  selfn%fscale=beq_fscale

  selfn%zetaopt=beq_zetaopt
  selfn%xiopt=beq_xiopt
  selfn%mrip=beq_mrip
  selfn%irip=beq_irip
  selfn%duct=beq_duct
  selfn%objadd=0
  selfn%objadd(GEOBJ_ABSORB)=absorber_objects
  selfn%objadd(GEOBJ_INVISI)=invisible_objects
  selfn%objadd(GEOBJ_SKYLIT)=skylight_objects
  selfn%objadd(GEOBJ_ERRLOS)=beancan_objects
  selfn%objadd(GEOBJ_CUTOUT)=cutout_objects
  selfn%skylpsi=skylight_flux_limits
  selfn%skyldbg=skylight_debug
  selfn%arip=beq_arip

  selfn%delpsi=beq_delpsi
  selfn%deltheta=beq_deltheta

  if(selfn%cenopt==1) then
     selfn%rcen=beq_rcen
     selfn%zcen=beq_zcen
  end if
  if(selfn%psiopt==1.OR.selfn%psiopt==3) then
     selfn%psimin=beq_psimin
     selfn%psimax=beq_psimax
  end if
  if(selfn%bdryopt==1.OR.selfn%bdryopt==5.OR.selfn%bdryopt==9) then
     selfn%psiref=beq_psiref
  end if
  if(selfn%nopt==1) then
     selfn%npsi=beq_npsi
     selfn%ntheta=beq_ntheta
  end if
  if(selfn%thetaopt==1) then
     selfn%thetamin=beq_thetamin
     selfn%thetamax=beq_thetamax
     selfn%thetaref=beq_thetaref
     ! move to usual polars, since thetamin/max set relative to X-point/vertical
     selfn%thetamin=selfn%thetamin-const_pid/2-selfn%thetaref
     selfn%thetamax=selfn%thetamax-const_pid/2-selfn%thetaref
  end if
  selfn%nzetp=beq_nzetap
  selfn%skipb=skip_eqdsk_b
  if(selfn%zetaopt==1) then
     selfn%zetamin=beq_zetamin
     selfn%zetamax=beq_zetamax
  else if(selfn%zetaopt==2) then
     if (beq_nzetap==0) then
        if (beq_mrip>0) then
           selfn%zetamin=-const_pid/beq_mrip
           selfn%zetamax=const_pid/beq_mrip
           selfn%nzetp=beq_mrip
        else
           selfn%zetamin=-const_pid
           selfn%zetamax=const_pid
        end if
     else
        selfn%zetamin=-const_pid/beq_nzetap
        selfn%zetamax=const_pid/beq_nzetap
     end if
  else if(selfn%zetaopt==3) then
     selfn%zetamin=-const_pid
     selfn%zetamax=const_pid
  else if(selfn%zetaopt==4) then
     selfn%zetamin=-3*const_pid/2
     selfn%zetamax=const_pid/2
  else if(selfn%zetaopt==5) then
     selfn%zetamin=-const_pid/2
     selfn%zetamax=3*const_pid/2
  end if
  selfn%fldspec=beq_fldspec
  selfn%psibig=beq_psibig
  selfn%leqok=equil_helicity_ok
  selfn%xsearch=beq_xsearch
  selfn%xrsta=beq_xrsta
  selfn%xrend=beq_xrend
  selfn%xzsta=beq_xzsta
  selfn%xzend=beq_xzend
  selfn%search=limiter_search
  selfn%lkrsta=search_r_start
  selfn%lkrend=search_r_end
  selfn%lkzsta=search_z_start
  selfn%lkzend=search_z_end
  selfn%vacfile=beq_vacuum_field_file

  selfn%skylcen=skylight_centre_line.OR.centre_line_cutout
  selfn%ctrackcode=GEOBJ_SKYLIT
  if (centre_line_cutout) then
     selfn%ctrackcode=GEOBJ_CUTOUT
  end if

  selfn%skyl=(selfn%objadd(GEOBJ_SKYLIT)>0)&
 &.OR.selfn%skylpsi.OR.(selfn%skyldbg>0)

end subroutine beq_readcon
!---------------------------------------------------------------------
!> calculate \f$ \psi \f$ limits
subroutine beq_psilt(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_psilt' !< subroutine name
  real(kr8) :: zdiff !< flux difference
  real(kr8) :: zdelta !< flux fuzzing factor
  real(kr8) :: zpsimin    !< user input minimum \f$ \psi \f$
  real(kr8) :: zpsimax    !< user input maximum \f$ \psi \f$

  if (self%n%psiopt==1) then
     ! scale user input from to range of EQDSK data
     ! NB psi may decrease outward, but always require  psimin<psimax
     ! NB psibdry may be different from psiqbdry from EQDSK (e.g. if limiter)
     zdiff=self%psibdry-self%psiaxis
     zpsimin=self%n%psimin
     zpsimax=self%n%psimax
     if (rsig>0) then
        !  psi increases outward
        self%n%psimin=self%psiaxis+zdiff*zpsimin
        self%n%psimax=self%psiaxis+zdiff*zpsimax
     else
        !  psi decreases outward
        self%n%psimin=self%psiaxis+zdiff*zpsimax
        self%n%psimax=self%psiaxis+zdiff*zpsimin
     end if
  else if (self%n%psiopt==2) then
     ! set limits based on geometry
     zdelta=self%n%delpsi*(self%psiotr-self%psiltr)
     if (rsig>0) then
        !  psi increases outward
        self%n%psimin=self%psiltr-zdelta
        self%n%psimax=self%psiotr+zdelta
     else
        !  psi decreases outward
        self%n%psimin=self%psiotr+zdelta
        self%n%psimax=self%psiltr-zdelta
     end if
  else if (self%n%psiopt==3) then
     zpsimin=min(self%n%psimin, self%n%psimax)
     zpsimax=max(self%n%psimin, self%n%psimax)
     !S     if (rsig>0) then
     !  psi increases outward
     self%n%psimin=zpsimin
     self%n%psimax=zpsimax
     !S     else
     !S        !  psi decreases outward
     !S        self%n%psimin=zpsimax
     !S        self%n%psimax=zpsimin
     !S     end if
  end if

end subroutine beq_psilt
!---------------------------------------------------------------------
!> calculate \f$ \theta \f$ limits
subroutine beq_thetalt(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_thetalt' !< subroutine name
  real(kr8) :: zdelta !< flux fuzzing factor

  if (self%n%thetaopt==1) then
     ! do nothing, user input
  else if (self%n%thetaopt==2) then
     ! set limits based on geometry
     zdelta=self%n%deltheta*(self%thetagmax-self%thetagmin)
     self%n%thetamin=self%thetagmin-zdelta
     self%n%thetamax=self%thetagmax+zdelta
  end if

end subroutine beq_thetalt
!---------------------------------------------------------------------
!> calculate \f$ \xi \f$ limits
subroutine beq_xilt(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='beq_xilt' !< subroutine name

  if(self%n%xiopt==1) then
     ! basically do nothing
     self%n%ximin=self%n%zetamin
     self%n%ximax=self%n%zetamax
     self%nzets=1
  else if(self%n%xiopt==2) then
     self%n%ximin=-const_pid
     self%n%ximax=const_pid
     if (self%n%nzetp/=0) then
        self%nzets=self%n%nzetp
     else
        self%nzets=self%n%mrip
     end if
  else if(self%n%xiopt==3) then
     call log_error(m_name,s_name,1,error_fatal,'xiopt not recognised')
  else if(self%n%xiopt==4) then
     self%n%ximin=-3*const_pid/2
     self%n%ximax=const_pid/2
     if (self%n%nzetp/=0) then
        self%nzets=self%n%nzetp
     else
        self%nzets=self%n%mrip
     end if
  else if(self%n%xiopt==5) then
     self%n%ximin=-const_pid/2
     self%n%ximax=3*const_pid/2
     if (self%n%nzetp/=0) then
        self%nzets=self%n%nzetp
     else
        self%nzets=self%n%mrip
     end if
  end if

end subroutine beq_xilt
!---------------------------------------------------------------------
!> calculate \f$ N_{\theta} \f$ and \f$ N_{\psi} \f$
subroutine beq_nset(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_nset' !< subroutine name

  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading vec')
  end if

end subroutine beq_nset
!---------------------------------------------------------------------
!> calculate  \f$ \psi \f$ at separatrix
subroutine beq_psix(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='beq_psix' !< subroutine name
  real(kr8), parameter :: zg1=1/(1+const_golden) !< factor for const_golden ratio search
  real(kr8), parameter :: zg2=1/const_golden !< factor for const_golden ratio search
  integer(ki4), parameter :: mcount=10000 !< limit on number of iterations in \f$ \theta \f$
  integer(ki4) :: mhemi=2 !< search below AND above midplane if two (xsearch=0)
  real(kr8) :: zepsg    !< normalised value of \f$ \epsilon_g \f$ tolerance
  real(kr8) :: zepsr    !< normalised value of \f$ \epsilon_r \f$ tolerance
  real(kr8) :: zeps    !< dynamically determined search tolerance
  real(kr8) :: zthet    !< \f$ \theta \f$ angle in mesh scan
  real(kr8) :: zdpdr    !< \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !< \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zsrr   !< estimate for maximum \f$ |R-R_c| \f$ in domain
  real(kr8) :: zszz   !< estimate for maximum \f$ |Z-Z_c| \f$ in domain
  real(kr8) :: zsrlt    !< minimum \f$ r \f$ or \f$ r^2 \f$ for search
  real(kr8) :: zsrsq    !< \f$ r^2 \f$ in search
  real(kr8) :: zgpsi    !< \f$ |\nabla\psi|^2 \f$
  real(kr8) :: zgpsimin    !< smallest \f$ |\nabla\psi|^2 \f$ in search
  integer(ki4) :: isrmin !< 1st array index of zgpsimin
  integer(ki4) :: jsrmin !< 2nd array index of zgpsimin
  real(kr8) :: zsrmin    !< minimum \f$ r^2 \f$ for search
  real(kr8) :: zsr1    !< minimum \f$ r \f$ bound for extremum
  real(kr8) :: zsr2    !< maximum \f$ r \f$ bound for extremum
  real(kr8) :: zsrxpt    !< \f$ r \f$ for extremum
  real(kr8) :: re    !<   \f$ R_i \f$
  real(kr8) :: ze    !<  \f$ Z_i \f$
  integer(ki4) :: i1 !< start index for mesh scan in \f$ R \f$
  integer(ki4) :: i2 !< stop index for mesh scan in \f$ R \f$
  integer(ki4) :: j1 !< start index for mesh scan in \f$ Z \f$
  integer(ki4) :: j2 !< stop index for mesh scan in \f$ Z \f$
  integer(ki4) :: jhemi !< loop over one or more hemispheres
  integer(ki4) :: jcount !< loop counter for search in \f$ \theta \f$
  real(kr8) :: zt1    !< \f$ \theta_1 \f$ search angle
  real(kr8) :: zt2    !< \f$ \theta_2 \f$ search angle
  real(kr8) :: zt3    !< \f$ \theta_3 \f$ search angle
  real(kr8) :: zt4    !< \f$ \theta_4 \f$ search angle
  real(kr8) :: zt5    !< \f$ \theta_5 \f$ search angle
  real(kr8) :: zpsi1    !< \f$ \psi_1 \f$ search result
  real(kr8) :: zpsi2    !< \f$ \psi_2 \f$ search result
  real(kr8) :: zpsi3    !< \f$ \psi_3 \f$ search result
  real(kr8) :: zpsi4    !< \f$ \psi_4 \f$ search result
  real(kr8) :: zpsi5    !< \f$ \psi_5 \f$ search result
  integer(ki4) :: is3chnged    !< \f$ \psi_3 \f$ needs changing to meet new tolerance
  integer(ki4) :: ierr    !< error code returned by r-extremum
  integer(ki4) :: ixf    !< code describing X-point find
  real(kr8) :: zsrm   !< extremal \f$ R \f$

  zsrmin=1.e+16
  zepsg=epsg*self%psinorm
  zepsr=epsr*self%psinorm
  ixf=0
  self%psixpt=rsig*100*self%psinorm
  zsrr=max( abs(self%rmax-self%n%rcen), abs(self%rmin-self%n%rcen) )
  zszz=max( abs(self%zmax-self%n%zcen), abs(self%zmin-self%n%zcen) )
  zsrlt=(zsrr**2+zszz**2)/100

  if (self%n%xsearch==1) mhemi=1
  do_hemi: do jhemi=0,mhemi-1

     if (self%n%xsearch==1) then
        zt2=-const_pid ; zt1=const_pid
        i1=1+max(int((self%n%xrsta-self%rmin)/self%dr),0)
        i2=2+min(int((self%n%xrend-self%rmin)/self%dr),self%mr-2)
        j1=1+max(int((self%n%xzsta-self%zmin)/self%dz),0)
        j2=2+min(int((self%n%xzend-self%zmin)/self%dz),self%mz-2)
     else
        zt2=(jhemi-1)*const_pid ; zt1=jhemi*const_pid
        j1=2+(self%mz/2)*jhemi
        j2=(self%mz/2)*(jhemi+1)-1
        i1=2
        i2=self%mr-1
     end if
     zsr1=zsrmin ; zsr2=-zsrmin
     isrmin=0 ; jsrmin=0
     zgpsimin=1.e+16

     ! step one, limiting r and theta for minimum of \f$ |\nabla\psi|^2 \f$
     ! one-a find mesh-point with smallest \f$ |\nabla\psi|^2 \f$
     do_meshj: do j=j1,j2
        ! note, skip edges of computational rectangle and avoid centre
        ze=self%zmin+(j-1)*self%dz
        do_meshi: do i=i1,i2
           re=self%rmin+(i-1)*self%dr
           zsrsq=(re-self%n%rcen)**2+(ze-self%n%zcen)**2
           if (zsrsq>zsrlt) then
              call spl2d_evaln(self%dpsidr,re,ze,1,zdpdr)
              call spl2d_evaln(self%dpsidz,re,ze,2,zdpdz)
              if (zdpdr*zdpdz<0.) then
                 zgpsi=zdpdr**2+zdpdz**2
                 if (zgpsi<zgpsimin) then
                    zgpsimin=zgpsi
                    isrmin=i
                    jsrmin=j
                 end if
              end if
           end if
        end do do_meshi
     end do do_meshj

     if (isrmin==0) then
        call log_error(m_name,s_name,1+jhemi,error_warning,'hemisphere: no X-point detected')
        exit do_hemi
     end if

     ! step one-b set r and theta limits for search using adjacent mesh-points
     do j=1,3
        ze=self%zmin+(jsrmin+j-3)*self%dz
        do i=1,3
           re=self%rmin+(isrmin+i-3)*self%dr
           zthet=atan2( ze-self%n%zcen, re-self%n%rcen )
           ! no shift if (zthet<-const_pid/2) zthet=2*const_pid+zthet
           zt1=min(zthet,zt1)
           zt2=max(zthet,zt2)
           zsrsq=(re-self%n%rcen)**2+(ze-self%n%zcen)**2
           zsr1=min(zsrsq,zsr1)
           zsr2=max(zsrsq,zsr2)
        end do
     end do
     zsr1=sqrt(zsr1)
     zsr2=sqrt(zsr2)

     ! step two, initialise for search
     zt3=zg1*zt1+zg2*zt2
     zeps=zepsr
     call beq_rextremum(self,zsr1,zsr2,zt1,zpsi1,zsrm,zeps,ierr)
     !DEXT write(*,*) 'zsr1,zsr2,zt1,zpsi1,zeps,ierr',zsr1,zsr2,zt1,zpsi1,zeps,ierr !DEXT
     if (ierr/=0) call log_error(m_name,s_name,10+ierr,error_fatal,'no extremum found where expected')
     call beq_rextremum(self,zsr1,zsr2,zt2,zpsi2,zsrm,zeps,ierr)
     !DEXT write(*,*) 'zsr1,zsr2,zt2,zpsi2,zeps,ierr',zsr1,zsr2,zt2,zpsi2,zeps,ierr !DEXT
     if (ierr/=0) call log_error(m_name,s_name,20+ierr,error_fatal,'no extremum found where expected')
     call beq_rextremum(self,zsr1,zsr2,zt3,zpsi3,zsrm,zeps,ierr)
     !DEXT write(*,*) 'zsr1,zsr2,zt3,zpsi3,zeps,ierr',zsr1,zsr2,zt3,zpsi3,zeps,ierr !DEXT
     if (ierr/=0) call log_error(m_name,s_name,30+ierr,error_fatal,'no extremum found where expected')
     ! step three, converge to X-point in theta, searching for a MAXimum in theta if rsig<0
     do_converge: do jcount=1,mcount
        is3chnged=0
        zeps=max( zepsg , epsr*(abs(zpsi3-zpsi1)+abs(zpsi3-zpsi2)) )
        zt4=zg1*zt1+zg2*zt3
        call beq_rextremum(self,zsr1,zsr2,zt4,zpsi4,zsrm,zeps,ierr)
        !DEXT write(*,*) 'zsr1,zsr2,zt4,zpsi4,zeps,ierr',zsr1,zsr2,zt4,zpsi4,zeps,ierr !DEXT
        if (ierr/=0) call log_error(m_name,s_name,40+ierr,error_fatal,'no extremum found where expected')
        if ((zpsi4-zpsi3)*rsig>0) then
           zpsi1=zpsi4
           zt1=zt4
        else
           zpsi2=zpsi3
           zpsi3=zpsi4
           zt2=zt3
           zt3=zt4
           is3chnged=1
        end if

        zt5=zg1*zt2+zg2*zt3
        call beq_rextremum(self,zsr1,zsr2,zt5,zpsi5,zsrm,zeps,ierr)
        !DEXT write(*,*) 'zsr1,zsr2,zt5,zpsi5,zeps,ierr',zsr1,zsr2,zt5,zpsi5,zeps,ierr !DEXT
        if (ierr/=0) call log_error(m_name,s_name,50+ierr,error_fatal,'no extremum found where expected')
        if ((zpsi5-zpsi3)*rsig>0) then
           zpsi2=zpsi5
           zt2=zt5
        else
           zpsi1=zpsi3
           zpsi3=zpsi5
           zt1=zt3
           zt3=zt5
           is3chnged=1
        end if
        ! reevaluate psi3 with tighter tolerance
        if (is3chnged==0) then
           call beq_rextremum(self,zsr1,zsr2,zt3,zpsi3,zsrm,zeps,ierr)
           !DEXT write(*,*) 'zsr1,zsr2,zt3,zpsi3,zeps,ierr',zsr1,zsr2,zt3,zpsi3,zeps,ierr !DEXT
           if (ierr/=0) call log_error(m_name,s_name,60+ierr,error_fatal,'no extremum found where expected')
        end if
        ! convergence test
        if (abs(zpsi3-zpsi1)+abs(zpsi3-zpsi2)<zepsg) exit

     end do do_converge

     if (jcount>=mcount) then
        call log_error(m_name,s_name,70+jhemi,error_warning,'no X-point found')
     else
        ixf=ixf+(jhemi+1)*3
        if ((self%psixpt-zpsi3)*rsig>0) then
           self%psixpt=zpsi3
           self%thetaxpt=zt3
           zsrxpt=(zsr1+zsr2)/2
        end if
     end if

  end do do_hemi

  if (ixf==0) then
     call log_error(m_name,s_name,80,error_fatal,'no X-point found')
  else
     call log_error(m_name,s_name,90,log_info,'X-point value')
     call log_value("SMITER-GEOQ psi-xpt ",self%psixpt)
     call log_value("SMITER-GEOQ theta-xpt ",self%thetaxpt)
     call log_value("psixpt find code",ixf)
     self%rxpt=self%n%rcen+zsrxpt*cos(self%thetaxpt)
     self%zxpt=self%n%zcen+zsrxpt*sin(self%thetaxpt)
     call log_value("SMITER-GEOQ R-xpt ",self%rxpt)
     call log_value("SMITER-GEOQ Z-xpt ",self%zxpt)
  end if

end subroutine beq_psix
!---------------------------------------------------------------------
!> calculate \f$ R_m \f$ and \f$ B_m \f$ from \f$ \psi_m \f$
subroutine beq_bdryrb(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='beq_bdryrb' !< subroutine name
  integer(ki4) :: nsrsamp    !< number of samples in \f$ r \f$
  real(kr8) :: ztheta    !< angle in  \f$ \theta \f$ of boundary point
  real(kr8) :: zdsrmin    !< floor to \f$ \Delta r_i \f$
  real(kr8) :: zp    !< \f$ \psi \f$
  real(kr8) :: zdpdr    !< \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !< \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zpsiinr    !<  estimate for inner limit of \f$ \psi \f$
  real(kr8) :: zpsioutr    !<  estimate for outer limit of \f$ \psi \f$
  real(kr8) :: zsrr   !< estimate for maximum \f$ |R-R_c| \f$ in domain
  real(kr8) :: zszz   !< estimate for maximum \f$ |Z-Z_c| \f$ in domain
  real(kr8) :: zsrsta   !< estimate for starting \f$ r \f$
  real(kr8) :: zsrmin   !< smaller \f$ r \f$ corresponding to zpsiinr
  real(kr8) :: zsrmax   !< larger \f$ r \f$ corresponding to zpsioutr
  real(kr8) :: zcos    !<  \f$ \cos(\theta=0) \f$
  real(kr8) :: zsin    !<  \f$ \sin(\theta=0) \f$
  real(kr8) :: re    !<   \f$ R_i \f$
  real(kr8) :: ze    !<  \f$ Z_i \f$
  real(kr8) :: zdpdsr    !<  \f$ \frac{\partial\psi}{\partial r}_i \f$
  real(kr8) :: zsr    !< \f$  r_i \f$
  real(kr8) :: zdsr    !< \f$ \Delta r_i \f$
  real(kr8) :: zdpdsrl    !< \f$ \frac{\partial\psi}{\partial r}_{i-1} \f$
  real(kr8) :: zsrinr    !<  smallest \f$ r \f$ in range
  real(kr8) :: zsroutr    !<  largest \f$ r \f$ in range
  integer(ki4) :: isr !< flag that value \f$ \psi < \psi_{\min} \f$ found
  integer(ki4) :: idplset !< flag that \f$ \frac{\partial\psi}{\partial r}_{i-1}  \f$ set
  integer(ki4) :: iext    !< maximum allowed number of knots for spline in \f$ r \f$
  integer(ki4) :: iknot    !< actual number of knots for spline in \f$ r \f$
  integer(ki4) :: intv    !< interval in which spline inverse found
  real(kr8) :: cpsi    !< constant for estimating \f$ \Delta r_i \f$
  real(kr8) :: zpsi    !<  \f$ \psi \f$
  real(kr8) :: zf    !<   \f$ f(\psi) = RB_T \f$
  real(kr8) :: cylj    !<  current estimated in cylindrical approx.
  real(kr8) :: zbr    !<  radial field component
  real(kr8) :: zbz    !<  vertical field component
  real(kr8) :: zbt    !<  toroidal field component

  pick_angle : select case (self%n%bdryopt)
  case(4,5,7,11,14) ! inboard point selected
     ztheta=const_pid
  case(8,9,10,12,15) ! outboard point selected
     ztheta=0.0_kr8
  case(1,3,13)
     ! check for replacement
     if (self%replasi) then
        call log_value('psi plasma boundary replaced with',self%psibdry)
        self%psiqbdry=self%psibdry
     end if
     return
  case default ! do nothing (assuming psiqbdry OK in eqdsk)
     return
  end select pick_angle

  zpsiinr=(self%psibdry+self%psiaxis)/2
  zsrr=max( abs(self%rmax-self%n%rcen), abs(self%rmin-self%n%rcen) )
  zszz=max( abs(self%zmax-self%n%zcen), abs(self%zmin-self%n%zcen) )
  zsrmax=sqrt(zsrr**2+zszz**2)
  nsrsamp=self%mr+self%mz
  zdsr=zsrmax/nsrsamp
  zsrsta=zsrmax/10
  zpsioutr=zpsiinr

  zcos=cos(ztheta)
  zsin=sin(ztheta)
  ! loop over distance from centre
  ! start a little way from origin
  zsr=zsrsta
  zsrinr=zsr
  isr=0
  idplset=0
  radial:do i=1,nsrsamp
     re=self%n%rcen+zsr*zcos
     if (re>=self%rmax.OR.re<=self%rmin) then
        !dbg         write(*,*) 're,ze=',re,ze
        !extremal R reached
        exit
     end if
     ze=self%n%zcen+zsr*zsin
     if (ze>=self%zmax.OR.ze<=self%zmin) then
        !dbg         write(*,*) 're,ze=',re,ze
        !extremal Z reached
        exit
     end if
     call spl2d_eval(self%psi,re,ze,zp)

     !dbg      write(*,*) j,i,re,ze,zp

     if (rsig>0) then
        ! if (self%psiaxis<self%psiqbdry)
        ! psi increasing outward
        if ((zp-zpsiinr)*rsig<0) then
           !record while psi < psi-inner
           isr=i
           zsrinr=zsr
        else
           if(isr==0) then
              ! sr sta is too large, reduce
              zsr=4*zsr/5
              cycle
           end if
        end if
     else
        ! psi decreasing outward
        if ((zp-zpsiinr)*rsig<0) then
           !if (zp>self%n%psimax)
           !record while psi > psi-inner
           isr=i
           zsrinr=zsr
        else
           if(isr==0) then
              ! sr sta is too large, reduce
              zsr=4*zsr/5
              cycle
           end if
        end if
     end if

     ! check monotone
     call spl2d_evaln(self%dpsidr,re,ze,1,zdpdr)
     call spl2d_evaln(self%dpsidz,re,ze,2,zdpdz)
     zdpdsr=zdpdr*zcos+zdpdz*zsin
     if(idplset==1) then
        if(zdpdsrl*zdpdsr<=0) then
           ! extremum in psi reached
           ! try to use the monotone range
           exit
        end if
     end if
     zsr=zsr+zdsr
     zdpdsrl=zdpdsr
     idplset=1
     zsroutr=zsr
     zpsioutr=zp
  end do radial

  if( (zpsiinr-zpsioutr)*rsig >= 0 ) then
     call log_error(m_name,s_name,1,error_fatal,'No suitable psi range exists')
  end if

  zsrmin=zsrinr
  zsrmax=zsroutr

  ! end of beq_rextrema equivalent code

  ! allocate generous amount of space for work (roughly half this should be needed, see  cpsi)
  iext=4*self%n%npsi
  allocate(wvextn(iext), stat=status)
  call log_alloc_check(m_name,s_name,11,status)
  allocate(wvext(iext), stat=status)
  call log_alloc_check(m_name,s_name,12,status)
  allocate(wvextd(iext), stat=status)
  call log_alloc_check(m_name,s_name,13,status)
  cpsi=abs(zpsioutr-zpsiinr)/(2*self%n%npsi)

  zdsrmin=(zsrmax-zsrmin)/(4*self%n%npsi-0.1)

  ! loop over r to define 1-D splines as a function of r
  zsr=zsrmin
  do i=1,iext
     re=self%n%rcen+zsr*zcos
     ze=self%n%zcen+zsr*zsin
     call spl2d_evaln(self%psi,re,ze,1,zp)
     call spl2d_evaln(self%dpsidr,re,ze,2,zdpdr)
     call spl2d_evaln(self%dpsidz,re,ze,2,zdpdz)
     zdpdsr=zdpdr*zcos+zdpdz*zsin
     if(i>1.AND.zdpdsrl*zdpdsr<=0) then
        ! should only be small non-monotone region, try to ignore
        call log_error(m_name,s_name,20,error_warning,'Lack of monotonicity')
        ! fix up
        zdpdsr=zdpdsrl
     end if
     wvextn(i)=zsr
     wvext(i)=zp
     wvextd(i)=zdpdsr
     iknot=i
     zdsr=abs(cpsi/zdpdsr)
     zsr=zsr+max(min(zdsr,3*zdsrmin),zdsrmin)
     if ( i>1 .AND. (zsr-zsrmax)>0 ) goto 2
     zdpdsrl=zdpdsr
  end do

  call log_error(m_name,s_name,25,error_fatal,'Too many nodes')

2     continue
  !DWVS write(*,*) s_name, ' debugging arrays' !DWVS
  !DWVS write(*,*) 'psi', self%psibdry !DWVS
  !DWVS do i=1,iknot !DWVS
  !DWVS write(*,*) i,wvextn(i),wvext(i),wvextd(i) !DWVS
  !DWVS end do !DWVS

  ! estimate interval where psi approx psim
  do j=1,iknot-1
     intv=j
     if ( (wvext(j)-self%psibdry)*(wvext(j+1)-self%psibdry) <= 0 ) exit
  end do
  zsr=wvextn(intv)

  zpsi=self%psibdry
  ! tc01a changes zsr to more exact value
  call tc01a(wvextn,wvext,wvextd,iknot,zsr,zpsi,intv)
  if (intv<=0) then
     call log_error(m_name,s_name,30,error_warning,'failure to improve r value')
  end if

  ! check for replacement
  if (self%replasi) then
     call log_value('psi plasma boundary replaced with',self%psibdry)
     self%psiqbdry=self%psibdry
  end if
  re=self%n%rcen+zsr*zcos
  ze=self%n%zcen+zsr*zsin
  call spl2d_evaln(self%dpsidr,re,ze,1,zdpdr)
  call spl2d_evaln(self%dpsidz,re,ze,2,zdpdz)
  self%rbdry=re
  self%bpbdry=(1/re)*sqrt( max(0.,(zdpdr**2+zdpdz**2)) )

  ! evaluate I aka f at psi
  call spleval(self%f,self%mr,self%psiaxis,self%psiqbdry,zpsi,zf,1)
  self%btotbdry=sqrt( max(0.,(self%bpbdry**2+(zf/re)**2)) )

  call log_error(m_name,s_name,2,log_info,'Reference boundary values')
  call log_value("SMITER-GEOQ psibdry ",self%psibdry)
  call log_value("SMITER-GEOQ rbdry ",self%rbdry)
  call log_value("SMITER-GEOQ zbdry ",ze)
  call log_value("SMITER-GEOQ bpbdry ",self%bpbdry)
  call log_value("SMITER-GEOQ btotbdry ",self%btotbdry)
  cylj=abs(zsr*self%bpbdry/2.e-7)
  call log_value("Estimated cylindrical current ",cylj)

  zbr=-(1/re)*zdpdz
  zbz=(1/re)*zdpdr
  zbt=zf/re
  call log_value("radial field component",zbr)
  call log_value("vertical field component",zbz)
  call log_value("toroidal field component",zbt)

  deallocate(wvextn)
  deallocate(wvext)
  deallocate(wvextd)

end subroutine beq_bdryrb
!---------------------------------------------------------------------
!> calculate \f$ r_{min} \f$ and \f$ r_{min} \f$ as functions of \f$ \theta_j \f$
subroutine beq_rextrema(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_rextrema' !< subroutine name
  integer(ki4) :: nsrsamp    !< number of samples in \f$ r \f$
  real(kr8) :: zp    !< \f$ \psi \f$
  real(kr8) :: zdpdr    !< \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !< \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: ztheta    !<  \f$ \theta_j \f$
  real(kr8) :: zsrr   !< estimate for maximum \f$ |R-R_c| \f$ in domain
  real(kr8) :: zszz   !< estimate for maximum \f$ |Z-Z_c| \f$ in domain
  real(kr8) :: zsrmin   !< estimate for starting \f$ r \f$
  real(kr8) :: zsrmax   !< estimate for maximum \f$ r \f$ in domain
  real(kr8) :: zcos    !<  \f$ \cos(\theta_j) \f$
  real(kr8) :: zsin    !<  \f$ \sin(\theta_j) \f$
  real(kr8) :: re    !<   \f$ R_i \f$
  real(kr8) :: ze    !<  \f$ Z_i \f$
  real(kr8) :: zdpdsr    !<  \f$ \frac{\partial\psi}{\partial r}_i \f$
  real(kr8) :: zsr    !< \f$  r_i \f$
  real(kr8) :: zdsr    !< \f$ \Delta r_i \f$
  real(kr8) :: zdpdsrl    !< \f$ \frac{\partial\psi}{\partial r}_{i-1} \f$
  real(kr8) :: zrpmin    !<  largest \f$ r \f$ giving  \f$  \psi < \psi_{\min} \f$
  real(kr8) :: zrpmax    !<  smallest \f$ r \f$ giving  \f$  \psi > \psi_{\max} \f$
  integer(ki4) :: isr !< flag that value \f$ \psi < \psi_{\min} \f$ found
  integer(ki4) :: idplset !< flag that \f$ \frac{\partial\psi}{\partial r}_{i-1}  \f$ set
  integer(ki4) :: im !< number of angles in largest interval of acceptable \f$ \theta \f$
  integer(ki4) :: i1 !< marks lower bound of current interval of acceptable \f$ \theta \f$
  integer(ki4) :: i1m !< marks lower bound of largest interval of acceptable \f$ \theta \f$

  ! work array, non-zero entry indicates non-monotone direction
  allocate(work1(self%n%ntheta+1), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  work1=0
  allocate(self%srmin(self%n%ntheta+1), stat=status)
  allocate(self%srmax(self%n%ntheta+1), stat=status)
  call log_alloc_check(m_name,s_name,2,status)

  self%dpsi=(self%n%psimax-self%n%psimin)/self%n%npsi
  zsrr=max( abs(self%rmax-self%n%rcen), abs(self%rmin-self%n%rcen) )
  zszz=max( abs(self%zmax-self%n%zcen), abs(self%zmin-self%n%zcen) )
  zsrmax=sqrt(zsrr**2+zszz**2)
  nsrsamp=self%mr+self%mz
  zdsr=zsrmax/nsrsamp
  zsrmin=zsrmax/10
  self%dtheta=(self%n%thetamax-self%n%thetamin)/self%n%ntheta

  ! loop over angle
  angle:do j=1,self%n%ntheta+1
     ztheta=self%n%thetamin+(j-1)*self%dtheta
     zcos=cos(ztheta)
     zsin=sin(ztheta)

     ! loop over distance from centre
     ! start a little way from origin
     zsr=zsrmin
     isr=0
     idplset=0
     radial:do i=1,nsrsamp
        re=self%n%rcen+zsr*zcos
        if (re>=self%rmax.OR.re<=self%rmin) then
           !dbg         write(*,*) 're,ze=',re,ze
           call log_error(m_name,s_name,10,error_warning,'rejecting direction-extremum not found')
           work1(j)=2
           exit
        end if
        ze=self%n%zcen+zsr*zsin
        if (ze>=self%zmax.OR.ze<=self%zmin) then
           !dbg         write(*,*) 're,ze=',re,ze
           call log_error(m_name,s_name,11,error_warning,'rejecting direction-extremum not found')
           work1(j)=2
           exit
        end if
        call spl2d_eval(self%psi,re,ze,zp)

        !dbg      write(*,*) j,i,re,ze,zp

        if (rsig>0) then
           ! if (self%psiaxis<self%psiqbdry)
           ! psi increasing outward
           if (zp<self%n%psimin) then
              !record
              isr=i
              zrpmin=zsr
           else
              if(isr==0) then
                 ! sr min is too large, reduce
                 zsr=4*zsr/5
                 cycle
              end if
           end if
           if (zp>self%n%psimax) then
              zrpmax=zsr
              exit
           end if
        else
           ! psi decreasing outward
           if (zp>self%n%psimax) then
              !record
              isr=i
              zrpmax=zsr
           else
              if(isr==0) then
                 ! sr min is too large, reduce
                 zsr=4*zsr/5
                 cycle
              end if
           end if
           if (zp<self%n%psimin) then
              zrpmin=zsr
              exit
           end if
        end if

        ! check monotone
        call spl2d_evaln(self%dpsidr,re,ze,1,zdpdr)
        call spl2d_evaln(self%dpsidz,re,ze,2,zdpdz)
        zdpdsr=zdpdr*zcos+zdpdz*zsin
        if(idplset==1) then
           if(zdpdsrl*zdpdsr<=0) then
              work1(j)=2
              exit
           end if
        end if
        zsr=zsr+zdsr
        zdpdsrl=zdpdsr
        idplset=1
     end do radial

     if (work1(j)<1) then
        self%srmin(j)=zrpmin
        self%srmax(j)=zrpmax
        ! reset sr min in terms of current value
        zsrmin=4*min(zrpmin,zrpmax)/5
     end if
  end do angle
  !dbg      write(*,*) work1

  ! find useable theta limits and save
  ! find largest interval with zero work1's
  i=0
  im=i
  i1=0
  i1m=i1
  do j=1,self%n%ntheta+1
     if (work1(j)<1) then
        i=i+1
        if (i>im) then
           ! interval is largest to date, so record
           i1m=i1
           im=i
        end if
     else
        ! reset counter and start marker
        i=0
        i1=j
     end if
  end do
  self%ntmin=i1m+1
  self%ntmax=i1m+im
  if (self%ntmin>=self%ntmax) then
     call log_error(m_name,s_name,20,error_fatal,'No directions contain psi interval')
  end if

  deallocate(work1)

end subroutine beq_rextrema
!---------------------------------------------------------------------
!> calculate \f$ r \f$ as a function of of \f$ \psi \f$
subroutine beq_rpsi(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='beq_rpsi' !< subroutine name
  integer(ki4) :: iext    !< maximum allowed number of knots for spline in \f$ r \f$
  integer(ki4) :: iknot    !< actual number of knots for spline in \f$ r \f$
  integer(ki4) :: intheta    !< actual number of angles defining  \f$ R,Z(\psi,\theta) \f$
  integer(ki4) :: intv    !< interval in which spline inverse found
  real(kr8) :: zsig    !< sign of  \f$ d \psi/dr \f$, value  \f$ \pm 1 \f$ (same as rsig)
  integer(ki4) :: isig    !< zero if  \f$  \psi \f$ decreases outward, else unity
  real(kr8) :: zdsrmin    !< floor to \f$ \Delta r_i \f$
  real(kr8) :: cpsi    !< constant for estimating \f$ \Delta r_i \f$
  real(kr8) :: ztheta    !<  \f$ \theta_j \f$
  real(kr8) :: ztmin    !<  minimum \f$ \theta \f$ for \f$ R,Z(\psi,\theta) \f$
  real(kr8) :: zdpdr    !< \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !< \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zsr    !< \f$ r \f$
  real(kr8) :: zcos    !<  \f$ \cos(\theta_j) \f$
  real(kr8) :: zsin    !<  \f$ \sin(\theta_j) \f$
  real(kr8) :: re    !<   \f$ R_i \f$
  real(kr8) :: ze    !<  \f$ Z_i \f$
  real(kr8) :: zdpdsr    !<  \f$ \frac{\partial\psi}{\partial r}_i \f$
  real(kr8) :: zdsr    !< \f$ \Delta r_i \f$
  real(kr8) :: zdpdsrl    !< \f$ \frac{\partial\psi}{\partial r}_{i-1} \f$
  real(kr8) :: zp    !<  \f$ \psi \f$
  real(kr8) :: zpsi    !<  \f$ \psi_i \f$
  real(kr8), dimension(:), allocatable :: swap !< workspace for reversing

  ! write(*,'(2G12.5)') (self%srmax(j),self%srmin(j),j=self%ntmin,self%ntmax)
  ! allocate generous amount of space for work (roughly half this should be needed, see  cpsi)
  iext=4*self%n%npsi
  allocate(wvextn(iext), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  allocate(wvext(iext), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  allocate(wvextd(iext), stat=status)
  call log_alloc_check(m_name,s_name,3,status)
  allocate(swap(iext), stat=status)
  call log_alloc_check(m_name,s_name,4,status)
  cpsi=(self%n%psimax-self%n%psimin)/(2*self%n%npsi)

  !! allocate work2 storage
  intheta=self%ntmax-self%ntmin+1
  if(self%n%npsi*intheta>0) then
     allocate(work2(self%n%npsi+1,intheta), stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     allocate(worka(self%n%npsi+1,intheta), stat=status)
     call log_alloc_check(m_name,s_name,11,status)
  else
     call log_error(m_name,s_name,12,error_fatal,'No theta data')
  end if
  ! loop over angle
  ij=0
  do j=self%ntmin,self%ntmax
     ij=ij+1
     ztheta=self%n%thetamin+(j-1)*self%dtheta
     zcos=cos(ztheta)
     zsin=sin(ztheta)
     ! srmin and srmax correspond to psimin and psimax, so possible that srmax <  srmin
     zdsrmin=abs((self%srmax(j)-self%srmin(j))/(4*self%n%npsi-1.1))
     zsig=sign(1._kr8,self%srmax(j)-self%srmin(j))
     isig=1
     if (zsig<0) isig=0
     !rev isig=0 ! reverse!!
     !rev if (zsig<0) isig=1 !!

     ! loop over r to define 1-D splines
     zsr=self%srmin(j)
     do i=1,iext
        re=self%n%rcen+zsr*zcos
        ze=self%n%zcen+zsr*zsin
        call spl2d_evaln(self%psi,re,ze,1,zp)
        call spl2d_evaln(self%dpsidr,re,ze,2,zdpdr)
        call spl2d_evaln(self%dpsidz,re,ze,2,zdpdz)
        zdpdsr=zdpdr*zcos+zdpdz*zsin
        if(i>1.AND.zdpdsrl*zdpdsr<=0) then
           ! should only be small non-monotone region, try to ignore
           call log_error(m_name,s_name,20,error_warning,'Lack of monotonicity')
           ! fix up
           zdpdsr=zdpdsrl
        end if
        wvextn(i)=zsr
        wvext(i)=zp
        wvextd(i)=zdpdsr
        iknot=i
        if ( i>1 .AND. (zsr-self%srmax(j))*(zsr-self%srmin(j)) > 0 ) goto 2
        zdsr=abs(cpsi/zdpdsr)
        zsr=zsr+zsig*max(zdsr,zdsrmin)
        !        write(*,*) 'zsr=',zsr
        zdpdsrl=zdpdsr
     end do
     !     write(*,*) 'i,j',i,j
     call log_error(m_name,s_name,21,error_fatal,'Too many nodes')

2        continue

     ! spline's independent variable must increase
     if (zsig<0) then
        call reverse(wvextn,swap,iknot)
        call reverse(wvext,swap,iknot)
        call reverse(wvextd,swap,iknot)
     end if

     ! loop over r using 1-D splines to calculate R and Z as functions (psi, theta)
     zsr=self%srmin(j)*isig+self%srmax(j)*(1-isig)
     do i=1,self%n%npsi+1
        zpsi=self%n%psimin +(i-1)*self%dpsi
        ! tc01a changes zsr to more exact value
        call tc01a(wvextn,wvext,wvextd,iknot,zsr,zpsi,intv)
        if (intv<=0) then
           call log_error(m_name,s_name,30,error_warning,'failure to improve r value')
        end if
        !dbg      if (ij==1) write(*,*) i,zsr
        work2(i,ij)=self%n%rcen+zsr*zcos
        worka(i,ij)=self%n%zcen+zsr*zsin
     end do

  end do

  ztmin=self%n%thetamin+(self%ntmin-1)*self%dtheta

  call spl2d_init(self%r,work2,self%n%npsi,intheta-1,&
 &self%n%psimin,ztmin,self%dpsi,self%dtheta,beq_spline_order)

  !dbg      write(*,*) ' check '

  call spl2d_init(self%z,worka,self%n%npsi,intheta-1,&
 &self%n%psimin,ztmin,self%dpsi,self%dtheta,beq_spline_order)

  deallocate(wvextn)
  deallocate(wvext)
  deallocate(wvextd)
  deallocate(swap)
  deallocate(work2)
  deallocate(worka)

end subroutine beq_rpsi
!---------------------------------------------------------------------
!> calculate spline coefficients for \f$ R/J (\psi, \theta) \f$
subroutine beq_rjpsi(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_rjpsi' !< subroutine name
  integer(ki4) :: intheta    !< actual number of angles defining  \f$ R,Z(\psi,\theta) \f$
  real(kr8) :: zpsi    !<  \f$ \psi_i \f$
  real(kr8) :: ztheta    !<  \f$ \theta_j \f$
  real(kr8) :: ztmin    !<  minimum \f$ \theta \f$ for \f$ R/J(\psi,\theta) \f$
  real(kr8) :: zr    !<   \f$ R(\psi_i,\theta_j) \f$
  real(kr8) :: zz    !<   \f$ Z(\psi_i,\theta_j) \f$
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R}_{i,j} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z}_{i,j} \f$
  real(kr8) :: zsrsq    !< \f$ r^2_{i,j} \f$
  real(kr8) :: zrgfac    !< reverse geometrical factor

  !! allocate work2 storage
  intheta=self%ntmax-self%ntmin+1
  allocate(work2(self%n%npsi+1,intheta), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  ! have adjusted \f$ f \f$ approximately so need to adjust back
  zrgfac=(self%rqcen-self%n%rmove)/self%rqcen

  ij=0
  do j=self%ntmin,self%ntmax
     ij=ij+1
     ztheta=self%n%thetamin+(j-1)*self%dtheta
     do i=1,self%n%npsi+1
        zpsi=self%n%psimin +(i-1)*self%dpsi
        call spl2d_evaln(self%r,zpsi,ztheta,1,zr)
        call spl2d_evaln(self%z,zpsi,ztheta,2,zz)
        call spl2d_evaln(self%dpsidr,zr,zz,1,zdpdr)
        call spl2d_evaln(self%dpsidz,zr,zz,2,zdpdz)
        zsrsq=(zr-self%n%rcen)**2+(zz-self%n%zcen)**2
        work2(i,ij)= zrgfac*(zr - self%n%rmove) &
 &      * (zdpdr*(zr-self%n%rcen)+zdpdz*(zz-self%n%zcen))/zsrsq
     end do
  end do

  ztmin=self%n%thetamin+(self%ntmin-1)*self%dtheta

  call spl2d_init(self%rjac,work2,self%n%npsi,intheta-1,&
 &self%n%psimin,ztmin,self%dpsi,self%dtheta,beq_spline_order)

  deallocate(work2)

end subroutine beq_rjpsi
!---------------------------------------------------------------------
!> calculate spline coefficients for \f$ R/I \times \psi \f$ derivatives
subroutine beq_rispld(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='beq_rispld' !< subroutine name
  real(kr8) :: zpsi    !<  \f$ \psi \f$
  real(kr8) :: zf    !<   \f$ f(\psi) = RB_T \f$
  real(kr8) :: zr    !<   \f$ R_i \f$
  real(kr8) :: zz    !<   \f$ Z_j \f$
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R}_{i,j} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z}_{i,j} \f$

  ! workr2 array for normalised Br
  !! allocate workr2 storage
  allocate(workr2(self%mr+1,self%mz+1), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  ! workz2 array for normalised Bz
  !! allocate workz2 storage
  allocate(workz2(self%mr+1,self%mz+1), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  ! evaluate R/I weighted derivatives
  do j=1,self%mz+1
     zz=self%zmin+(j-1)*self%dz
     do i=1,self%mr+1
        zr=self%rmin +(i-1)*self%dr
        call spl2d_evaln(self%dpsidr,zr,zz,1,zdpdr)
        call spl2d_evaln(self%dpsidz,zr,zz,2,zdpdz)
        call spl2d_evaln(self%psi,zr,zz,2,zpsi)
        ! evaluate I aka f at psi
        call spleval(self%f,self%mr,self%psiaxis,self%psiqbdry,zpsi,zf,1)
        workr2(i,j)=-zr*zdpdz/zf
        workz2(i,j)=zr*zdpdr/zf
     end do
  end do
  ! produce spline representation
  call spl2d_init(self%rispldz,workr2,self%mr,self%mz,&
 &self%rmin,self%zmin,self%dr,self%dz,beq_spline_order)
  deallocate(workr2)
  call spl2d_init(self%rispldr,workz2,self%mr,self%mz,&
 &self%rmin,self%zmin,self%dr,self%dz,beq_spline_order)
  deallocate(workz2)

end subroutine beq_rispld
!---------------------------------------------------------------------
!> reverse order of array (private)
subroutine reverse(parr,swap,kn)

  !! arguments
  real(kr8), dimension(kn), intent(inout) :: parr  !< array to be reversed
  real(kr8), dimension(kn), intent(inout) :: swap   !< swap space
  integer(ki4), intent(in) :: kn    !< number of entries in parr

  do k=1,kn
     l=kn+1-k
     swap(k)=parr(l)
  end do
  parr=swap

end subroutine reverse
!---------------------------------------------------------------------
!> calculate B-field components at position given by posang
subroutine beq_b(self,posang,kopt)

  !! arguments
  type(beq_t), intent(inout) :: self !< object data structure
  type(posang_t), intent(inout)  :: posang !< position and vector involving angles
  integer(ki4), intent(in) :: kopt   !< selector

  !! local
  character(*), parameter :: s_name='beq_b' !< subroutine name
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  real(kr8) :: zpsi    !<  \f$ \psi(R,Z) \f$
  real(kr8) :: zf    !<   \f$ f(\psi) = RB_T \f$
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zrfac    !<  \f$ R(x,y,z) \f$ factor
  !     real(kr8), dimension(3) :: zbrip    !<  Ripple field in Cartesians
  type(posang_t)  :: zparip !< position and vector for ripple
  real(kr8), dimension(3)  :: zvecf    !<  Vacuum field

  if(kopt==0) then

     ! calculate cartesian vector at cartesian posang
     ! first position to cyl polars to get psi derivs
     call posang_invtfm(posang,0)
     zr=posang%pos(1)
     zz=posang%pos(2)
     call spl2d_evaln(self%dpsidr,zr,zz,1,zdpdr)
     call spl2d_evaln(self%dpsidz,zr,zz,2,zdpdz)
     call spl2d_evaln(self%psi,zr,zz,2,zpsi)
     ! evaluate I aka f at psi
     call spleval(self%f,self%mr,self%psiaxis,self%psiqbdry,zpsi,zf,1)
     ! B in toroidal-cyl polars
     posang%vec(1)=-zdpdz/zr ; posang%vec(2)=zdpdr/zr ; posang%vec(3)=zf/zr
     posang%opt=17
     zparip=posang
     zparip%vec=0
     ! get ripple
     if (self%n%vacfile(1:4)=='null') then
        if (self%n%mrip/=0) then
           ! 'ripple' includes toroidal field
           posang%vec(3)=0
           zparip=posang
           call beq_ripple(self,zparip,1)
        end if
     else
        ! 'ripple' includes toroidal field
        posang%vec(3)=0
        call spl3d_evalm(self%vacfld,zparip%pos,zparip%vec)
        ! arranged so that spl3d stores F=R^(rpower)x(BR,BZ,BT/R=I/R^2) (rpower=1)
        !                                = (R.BR,R.BZ,I/R)
        ! following converts F to (BR,BZ,BT)
        zrfac=1/zr**spl3d_rpower
        zparip%vec=zrfac*zparip%vec ; zparip%vec(3)=zr*zparip%vec(3)
     end if
     ! (BR,BZ,BT) vector to cartesians
     call posang_tfm(posang,-3)
     call posang_tfm(zparip,-3)
     posang%vec=posang%vec+zparip%vec
     !        posang%vec=zparip%vec

  else if (kopt==1) then

     if (self%n%vacfile=='null') then
        call log_error(m_name,s_name,1,error_fatal,'No coding for this case')
     else
        ! calculate vector (RBR,RBZ,BT) at quantised coordinate
        zr=posang%pos(1)
        zz=posang%pos(2)
        call spl2d_evaln(self%dpsidr,zr,zz,1,zdpdr)
        call spl2d_evaln(self%dpsidz,zr,zz,2,zdpdz)
        call spl3d_evalm(self%vacfld,posang%pos,zvecf)
        posang%vec(1)=-zdpdz+zvecf(1)
        posang%vec(2)=zdpdr+zvecf(2)
        posang%vec(3)=zvecf(3)
     end if

  end if

end subroutine beq_b
!---------------------------------------------------------------------
!> calculate B-field ripple components at position given by posang
subroutine beq_ripple(self,posang,kresult,kopt)

  !! arguments
  type(beq_t), intent(inout) :: self !< object data structure
  type(posang_t), intent(inout)  :: posang !< position and vector
  !     type(posang_t), intent(out)  :: pbrip !< ripple vector
  integer(ki4), intent(in):: kresult   !< type of result
  integer(ki4), intent(in), optional :: kopt   !< optional selector

  !! local
  character(*), parameter :: s_name='beq_ripple' !< subroutine name

  if(.NOT.present(kopt)) then

     ! default is MAST ripple

     if (kresult==1) then
        ! need all three components
        if (posang%opt==17) then
           ! work in R,Z,zeta coordinates
           call beq_ripple_calch1(self,posang,1)
        end if
     end if
  end if

end subroutine beq_ripple
!---------------------------------------------------------------------
!> calculate B-field ripple as for MAST
subroutine beq_ripple_calch1(self,posang,ktype)

  !! arguments
  type(beq_t), intent(inout) :: self !< object data structure
  type(posang_t), intent(inout)  :: posang !< position and vector
  integer(ki4), intent(in):: ktype   !< type of result

  !! local
  character(*), parameter :: s_name='beq_ripple_calch1' !< subroutine name
  real(kr8) :: zr    !<   \f$ R \f$
  real(kr8) :: zz    !<   \f$ Z \f$
  real(kr8) :: zzeta    !<  \f$ \zeta \f$
  real(kr8) :: zh    !<   \f$ H \f$
  real(kr8) :: zdhdr    !<  \f$ \frac{\partial H}{\partial R} \f$
  real(kr8) :: zdhdp    !<  \f$ \frac{\partial H}{\partial \zeta} \f$
  real(kr8) :: zfac     !< normalisation factor

  zr=posang%pos(1)
  zz=posang%pos(2)
  zzeta=posang%pos(3)
  if (ktype==0) then
     ! return H and dhdR only
     zh=beq_ripple_h1(zr,zz,zzeta,0,self%n%mrip,self%ivac,self%n%arip)
     zdhdr=beq_ripple_h1(zr,zz,zzeta,1,self%n%mrip,self%ivac,self%n%arip)
     posang%vec(1)=zh ; posang%vec(2)=zdhdr ; posang%vec(3)=0
     posang%opt=17
  else if (ktype==1) then
     ! all three components in R-Z-zeta
     zdhdr=beq_ripple_h1(zr,zz,zzeta,1,self%n%mrip,self%ivac,self%n%arip)
     zdhdp=beq_ripple_h1(zr,zz,zzeta,2,self%n%mrip,self%ivac,self%n%arip)
     ! special normalisation for this ripple function
     ! minus which might be expected cancels with minus in definition of B as curl
     zfac=1._kr8/(real(self%n%mrip)*zr)
     ! B in toroidal-cyl polars
     posang%vec(1)=-zdhdp*zfac ; posang%vec(2)=0 ; posang%vec(3)=zr*zdhdr*zfac
     posang%opt=17
  end if

end subroutine beq_ripple_calch1
!---------------------------------------------------------------------
!> calculate complete track of psi extremum through plasma centre
subroutine beq_ctrax(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure

  !! local
  character(*), parameter :: s_name='beq_ctrack1' !< subroutine name
  real(kr8), dimension(:,:), allocatable :: ctrackrz !< central track
  integer(ki4) :: nctrack !<  size of ctrackrz array

  call beq_ctrackc(self,ctrackrz,nctrack,1)
  call beq_ctrack1(self,ctrackrz,nctrack,1)
  deallocate(ctrackrz)
  call beq_ctrackc(self,ctrackrz,nctrack,2)
  call beq_ctrack1(self,ctrackrz,nctrack,2)
  deallocate(ctrackrz)

end subroutine beq_ctrax
!---------------------------------------------------------------------
!> calculate track of psi extremum through plasma centre
subroutine beq_ctrackc(self,ctrackrz,nctrack,ktyp)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  real(kr8), dimension(:,:), allocatable, intent(out) :: ctrackrz !< central track
  integer(ki4), intent(out) :: nctrack !<  bound for ctrackrz array
  integer(ki4), intent(in) :: ktyp !<  lower or upper skylight type

  !! local
  character(*), parameter :: s_name='beq_ctrackc' !< subroutine name
  real(kr8) :: zepsg    !< normalised value of \f$ \epsilon_g \f$ tolerance
  integer(ki4) :: iwaydir !< ! direction of travel in \f$ Z \f$ (away from centre)
  integer(ki4) :: im !< ! size of array in direction of travel in \f$ Z \f$
  integer(ki4) :: igr !< \f$ R \f$ position index of global extremum
  integer(ki4) :: igz !< \f$ Z \f$ position index of global extremum
  integer(ki4) :: soutr=10 !< maximum number of outer searches
  integer(ki4) :: sinr=10 !< maximum number of inner searches
  integer(ki4) :: isr !< \f$ R \f$ search (increment) direction \f$ \pm 1 \f$
  integer(ki4) :: isz !< \f$ Z \f$ search (increment) direction \f$ \pm 1 \f$
  integer(ki4) :: ir !< current \f$ R \f$ position as index
  integer(ki4) :: iz !< current \f$ Z \f$ position as index
  real(kr8) :: zpp !< current value of \f$ \psi \f$
  real(kr8) :: zpg !< value of \f$ \psi \f$ at global extremum
  real(kr8) :: zt1    !< \f$ \theta_1 \f$ search angle
  real(kr8) :: zsr1    !< minimum \f$ r \f$ bound for extremum
  real(kr8) :: zsr2    !< maximum \f$ r \f$ bound for extremum
  real(kr8) :: zrm    !< \f$ r \f$ at extremum
  real(kr8) :: zpm    !< \f$ \psi \f$ at extremum
  real(kr8) :: zr    !<   \f$ R(\psi_i,\theta_j) \f$
  real(kr8) :: zz    !<   \f$ Z(\psi_i,\theta_j) \f$
  integer(ki4) :: ierr    !< error code returned by r-extremum

  ! direction of travel in \f$ Z \f$
  iwaydir=2*ktyp-3
  ! start  at centre
  zr=self%n%rcen
  zz=self%n%zcen
  igr=1+(zr-self%rmin)/self%dr
  igz=1+(zz-self%zmin)/self%dz

  !! allocate track storage
  im=(ktyp-1)*self%mz-iwaydir*igz
  allocate(ctrackrz(im,2), stat=status)
  call log_alloc_check(m_name,s_name,1,status)

  zepsg=epsg*self%psinorm
  ctrackrz(im,1)=zr
  ctrackrz(im,2)=zz

  iz=im
  loop_over_z: do k=2,im
     igz=igz+iwaydir
     iz=iz-1
     ! set initially positive search  directions for increasing psi
     isr=1
     isz=1

     dosoutr: do j=1,soutr
        ir=igr
        call spl2d_eval(self%psi,self%rmin+(ir-1)*self%dr,self%zmin+(igz-1)*self%dz,zpg)

        !! search in r
        ! change direction if necessary
        ir=igr+isr
        call spl2d_eval(self%psi,self%rmin+(ir-1)*self%dr,self%zmin+(igz-1)*self%dz,zpp)
        if ( (zpp-zpg)*rsig < 0 ) then
           igr=ir
           zpg=zpp
        else
           isr=-1
           ir=igr
        end if

        do i=1,sinr
           ir=ir+isr
           call spl2d_eval(self%psi,self%rmin+(ir-1)*self%dr,self%zmin+(igz-1)*self%dz,zpp)
           if ( (zpp-zpg)*rsig < 0 ) then
              igr=ir
              zpg=zpp
           else
              exit dosoutr
           end if
        end do

     end do dosoutr

     zt1=0
     zsr1=self%rmin+(igr-2)*self%dr
     zsr2=self%rmin+igr*self%dr
     zz=self%zmin+(igz-1)*self%dz
     call beq_rextremum(self,zsr1,zsr2,zt1,zpm,zrm,zepsg,ierr)
     if (ierr/=0) then
        call log_error(m_name,s_name,ierr,error_warning,'no extremum found where expected')
        call log_value('line Z = ',zz)
     else
        zr=zrm
     end if
     ctrackrz(iz,1)=zr
     ctrackrz(iz,2)=zz

  end do loop_over_z

  nctrack=im

  !CDBG  idunit=799 !CDBG
  !CDBG  !open(newunit=idunit)
  !CDBG  !open(unit=idunit, file='extr.gnu') !CDBG
  !CDBG  write(idunit,'(A)') '# track of extremum through plasma centre' !CDBG
  !CDBG  write(idunit,'(1P, 2(1X, G13.5))') (ctrackrz(k,1),ctrackrz(k,2),k=1,nctrack) !CDBG

end subroutine beq_ctrackc
!---------------------------------------------------------------------
!> update of psi extremum through plasma centre
subroutine beq_ctrack1(self,ptrack,kn,kcall)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  real(kr8), dimension(:,:), allocatable, intent(in) :: ptrack !< central track
  integer(ki4), intent(in) :: kn !<  bound for ptrack array
  integer(ki4), intent(in) :: kcall !<  number of call

  !! local
  character(*), parameter :: s_name='beq_ctrack1' !< subroutine name
  integer(ki4) :: jt !< ! loop variable

  if (kcall==1) then
     allocate(self%ctrackrz(self%mz,2), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     self%ctrackrz(1:kn,1:2)=ptrack
     self%nctrack=kn
  else
     do jt=1,kn
        self%ctrackrz(self%mz+1-jt,:)=ptrack(jt,:)
     end do
     self%nctrack=self%mz
  end if

end subroutine beq_ctrack1
!---------------------------------------------------------------------
!> return track of psi extremum through plasma centre
subroutine beq_ctrack1r(self,ptrack,kn,ktyp)

  !! arguments
  type(beq_t), intent(in) :: self   !< object data structure
  real(kr8), dimension(:,:), allocatable, intent(out) :: ptrack !< central track
  integer(ki4), intent(out) :: kn !<  bound for ptrack array
  integer(ki4), intent(in) :: ktyp !<  lower or upper skylight type

  !! local
  character(*), parameter :: s_name='beq_ctrack1r' !< subroutine name
  integer(ki4) :: iwaydir !< ! direction of travel in \f$ Z \f$ (away from centre)
  integer(ki4) :: igz !< \f$ Z \f$ position index of global extremum
  integer(ki4) :: im !< ! size of array in direction of travel in \f$ Z \f$
  integer(ki4) :: jt !< ! loop variable

  ! direction of travel in \f$ Z \f$
  iwaydir=2*ktyp-3
  ! start  at centre
  igz=1+(self%n%zcen-self%zmin)/self%dz

  ! direction of travel in \f$ Z \f$
  im=(ktyp-1)*self%mz-iwaydir*igz
  if (ktyp==1) then
     allocate(ptrack(im,2), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     ptrack=self%ctrackrz(1:im,:)
  else
     do jt=1,im
        ptrack(jt,:)=self%ctrackrz(self%mz+1-jt,:)
     end do
  end if
  kn=im

end subroutine beq_ctrack1r
!---------------------------------------------------------------------
!> query position relative to track of psi extremum through plasma centre
subroutine beq_ctrackcq(self,ctrackrz,nctrack,prz,pout)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  real(kr8), dimension(:,:), allocatable, intent(in) :: ctrackrz !< central track
  integer(ki4), intent(in) :: nctrack !<  bound for ctrackrz array
  real(kr8), dimension(2), intent(in) :: prz    !<  Query point position \f$ (R,Z) \f$
  real(kr8), intent(out) :: pout !<  point inside (negative, HFS) or outside (positive, LFS) track

  !! local
  character(*), parameter :: s_name='beq_ctrackcq' !< subroutine name
  real(kr8) :: zr    !<  Query point \f$ R \f$
  real(kr8) :: zz    !<  Query point \f$ Z \f$
  integer(ki4) :: iz    !<  array entry in \f$ Z \f$ direction
  real(kr8) :: remz    !<  remainder in \f$ Z \f$ direction
  real(kr8) :: zrc    !<  Central \f$ R \f$ at specified \f$ Z \f$

  zr=prz(1)
  zz=prz(2)
  iz=1+max(int((zz-self%zmin)/self%dz),0)
  remz=zz-self%zmin-(iz-1)*self%dz
  iz=min(iz,nctrack-1)
  zrc=(1-remz)*ctrackrz(iz,1)+remz*ctrackrz(iz+1,1)

  pout=zr-zrc

end subroutine beq_ctrackcq
!---------------------------------------------------------------------
!> find \f$ (R, Z) \f$ position of \f$ \psi \f$ value on track through plasma centre
subroutine beq_ctrackpt(self,ctrackrz,nctrack,prz,psi,ktyps)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  real(kr8), dimension(:,:), allocatable, intent(in) :: ctrackrz !< central track
  integer(ki4), intent(in) :: nctrack !<  bound for ctrackrz array
  real(kr8), dimension(2), intent(out) :: prz    !<  Query point position \f$ (R,Z) \f$
  real(kr8), intent(in) :: psi !<   \f$ \psi \f$ value
  integer(ki4), intent(in) :: ktyps !<  lower (1) or upper (2) half of cross-section

  !! local
  character(*), parameter :: s_name='beq_ctrackpt' !< subroutine name
  integer(ki4) :: idir !< ! direction of travel in \f$ Z \f$, (-1) in lower, (+1) in upper
  real(kr8) :: zra    !<  Old value of \f$ R \f$
  real(kr8) :: zza    !<  Old value of \f$ Z \f$
  real(kr8) :: zr    !<  Query point \f$ R \f$
  real(kr8) :: zz    !<  Query point \f$ Z \f$
  integer(ki4) :: ista    !<  loop start at plasma centre
  integer(ki4) :: iend    !<  loop end at flux array boundary
  integer(ki4) :: iz    !<  array entry in \f$ Z \f$ direction
  real(kr8) :: zpsia    !<  first value of \f$ \psi \f$
  real(kr8) :: zpsib    !<  second value of \f$ \psi \f$
  real(kr8) :: dpmin    !<  distance of point from target in \f$ \psi \f$ metric
  real(kr8) :: dpminew    !<  distance of point from target in \f$ \psi \f$ metric

  ! direction of travel in \f$ Z \f$
  idir=2*ktyps-3
  ! loop limits
  ista=1+max(int((self%n%zcen-self%zmin)/self%dz),0)
  iend=((1-idir)+(idir+1)*self%mz)/2

  ! main search loop
  zr=ctrackrz(ista,1)
  zz=ctrackrz(ista,2)
  call spl2d_eval(self%psi,zr,zz,zpsia)
  dpmin=abs(zpsia-psi)
  do i=ista+idir,iend,idir
     zra=zr
     zza=zz
     zr=ctrackrz(i,1)
     zz=ctrackrz(i,2)
     call spl2d_eval(self%psi,zr,zz,zpsib)
     !if ((zpsia-psi)*(zpsib-psi)<=0) exit
     dpminew=abs(zpsib-psi)
     if (dpminew>dpmin) exit
     dpmin=dpminew
     zpsia=zpsib
  end do

  if (i==iend+idir) then
     call log_error(m_name,s_name,1,error_warning,'no position found for flux')
  end if

  prz(1)=zra
  prz(2)=zza

end subroutine beq_ctrackpt
!---------------------------------------------------------------------
!> get unit number for beq writing
subroutine beq_getunit(kunit)

  !! arguments
  integer, intent(out) :: kunit    !< log unit number

  kunit=nin

end subroutine beq_getunit
!---------------------------------------------------------------------
!> define ripple as for MAST
function beq_ripple_h1(pr,pz,pzeta,kopt,km,pin,pa)

  !! arguments
  real(kr8), intent(in) :: pr    !<   \f$ R \f$
  real(kr8), intent(in) :: pz    !<   \f$ Z \f$
  real(kr8), intent(in) :: pzeta    !<  \f$ \zeta \f$
  integer(ki4), intent(in) :: kopt   !< option selector
  integer(ki4), intent(in), optional  :: km   !< number of coils
  real(kr8), intent(in), optional  :: pin    !< current/\f$ \zeta \f$ normalisation
  real(kr8), intent(in), optional  :: pa    !< ripple parameter/ \f$ R \f$ normalisation

  !! local
  character(*), parameter :: s_name='beq_ripple_h1' !< subroutine name
  real(kr8) :: beq_ripple_h1 !< local variable
  real(kr8) :: zfun    !<  function value
  integer(ki4), save :: im=1   !< number of coils
  real(kr8), save :: zi=0._kr8    !< plasma current
  real(kr8), save :: za=0._kr8    !< ripple parameter
  real(kr8), save :: zrnorm=1._kr8    !< \f$ R \f$ normalisation
  real(kr8), save :: zzetanorm=1._kr8    !< \f$ \zeta \f$ normalisation
  real(kr8), save :: zroff=0._kr8    !< \f$ R \f$ offset
  real(kr8), save :: zzetaoff=0._kr8    !< \f$ \zeta \f$ offset
  real(kr8), save :: zrnns=1._kr8    !< \f$ \zeta \f$ scaling
  real(kr8), save :: zarnns=0._kr8    !< \f$ \zeta \f$ component scaling
  real(kr8) :: zr    !< \f$ R \f$ normalised
  real(kr8) :: zzeta    !< \f$ \zeta \f$ normalised

  zfun=0
  if (kopt<0) then

     if(present(km)) then
        if(km==-3) then
           zrnns=pin
           zarnns=pa
        else if(km==-2) then
           zroff=pin
           zzetaoff=pa
        else if(km==-1) then
           zrnorm=pin
           zzetanorm=pa
        else
           !  save coil data
           im=km
           zi=pin
           za=pa
           zrnns=real(im)
           zarnns=za*real(im)
        end if
     end if

  else

     zr=(pr-zroff)*zrnorm
     zzeta=(pzeta-zzetaoff)*zzetanorm
     if(kopt==0) then
        !  function value
        zfun=real(im)*zi*log(zr)+za*(zr**im)*cos(zrnns*zzeta)
     else if (kopt==1) then
        ! derivative wrt \f$ R \f$
        zfun=real(im)/zr*(zi+za*(zr**im)*cos(zrnns*zzeta))
     else if (kopt==2) then
        ! derivative wrt \f$ \zeta \f$ strictly \f$ \xi = N_s \zeta \f$
        zfun=-zarnns*(zr**im)*sin(zrnns*zzeta)
     else if (kopt==3) then
        ! derivative wrt normalised \f$ R \f$
        zfun=zrnorm*real(im)/zr*(zi+za*(zr**im)*cos(zrnns*zzeta))
     else if (kopt==4) then
        ! derivative wrt normalised \f$ \xi \f$
        zfun=-zzetanorm*zarnns*(zr**im)*sin(zrnns*zzeta)
     end if
     !        write(*,*) kopt,zr,zzeta,zfun

  end if

  beq_ripple_h1=zfun

end function beq_ripple_h1
!---------------------------------------------------------------------
!> calculate spline coefficients for \f$ I (\psi) \f$
subroutine beq_ipsi(self)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure


  !! local
  character(*), parameter :: s_name='beq_ipsi' !< subroutine name
  real(kr8) :: zpsi    !<  \f$ \psi_i \f$

  !! allocate work1 storage
  allocate(work1(self%n%npsi+1), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  !! allocate self%i storage
  allocate(self%i(self%n%npsi+1), stat=status)
  call log_alloc_check(m_name,s_name,2,status)

  do i=1,self%n%npsi+1
     zpsi=self%n%psimin +(i-1)*self%dpsi
     work1(i)=zpsi
     self%i(i)=self%f(i)
  end do

  !     call spl2d_init(self%iac,work2,self%n%npsi,intheta-1,&
  !    &self%n%psimin,ztmin,self%dpsi,self%dtheta,beq_spline_order)

  deallocate(work1)

end subroutine beq_ipsi
!---------------------------------------------------------------------
!> calculate  extremal \f$ \psi \f$ in search direction
subroutine beq_rextremum(self,psr1,psr2,ptheta,psim,psrm,peps,kerr)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  real(kr8), intent(in) :: psr1   !< lower limiting value of \f$ r \f$
  real(kr8), intent(in) :: psr2   !< upper limiting value of \f$ r \f$
  real(kr8), intent(in) :: ptheta   !< search direction \f$ \theta \f$
  real(kr8), intent(out) :: psim   !< extremal \f$ \psi_m \f$
  real(kr8), intent(out) :: psrm   !< extremal \f$ R \f$
  real(kr8), intent(in) :: peps   !< convergence tolerance
  integer(ki4), intent(out) :: kerr    !< return code (non-zero means trouble)

  !! local
  character(*), parameter :: s_name='beq_rextremum' !< subroutine name
  real(kr8), parameter :: zg1=1/(1+const_golden) !< factor for const_golden ratio search
  real(kr8), parameter :: zg2=1/const_golden !< factor for const_golden ratio search
  integer(ki4), parameter :: micount=10000 !< limit on number of iterations in \f$ r \f$
  real(kr8) :: zcos    !<  \f$ \cos(\theta) \f$
  real(kr8) :: zsin    !<  \f$ \sin(\theta) \f$
  real(kr8) :: zsr1    !< \f$ r_1 \f$ search radius
  real(kr8) :: zsr2    !< \f$ r_2 \f$ search radius
  real(kr8) :: zsr3    !< \f$ r_3 \f$ search radius
  real(kr8) :: zsr4    !< \f$ r_4 \f$ search radius
  real(kr8) :: zsr5    !< \f$ r_5 \f$ search radius
  real(kr8) :: zpsi1    !< \f$ \psi_1 \f$ search result
  real(kr8) :: zpsi2    !< \f$ \psi_2 \f$ search result
  real(kr8) :: zpsi3    !< \f$ \psi_3 \f$ search result
  real(kr8) :: zpsi4    !< \f$ \psi_4 \f$ search result
  real(kr8) :: zpsi5    !< \f$ \psi_5 \f$ search result
  integer(ki4) :: jicount !< loop counter

  kerr=0
  ! initialise for search
  zcos=cos(ptheta)
  zsin=sin(ptheta)
  ! initial search range
  zsr1=psr1
  zsr2=psr2
  ! need 3 points to begin
  zsr3=zg1*zsr1+zg2*zsr2
  call spl2d_eval(self%psi,self%n%rcen+zsr1*zcos,self%n%zcen+zsr1*zsin,zpsi1)
  call spl2d_eval(self%psi,self%n%rcen+zsr2*zcos,self%n%zcen+zsr2*zsin,zpsi2)
  call spl2d_eval(self%psi,self%n%rcen+zsr3*zcos,self%n%zcen+zsr3*zsin,zpsi3)

  ! converge to extremum in radius, searching for a MINimum if rsig<0
  ! use a fourth point to reduce range where SINGLE extremum lies
  do_count: do jicount=1,micount
     ! insert fourth at left
     zsr4=zg1*zsr1+zg2*zsr3
     call spl2d_eval(self%psi,self%n%rcen+zsr4*zcos,self%n%zcen+zsr4*zsin,zpsi4)
     if ((zpsi4-zpsi3)*rsig<0) then
        zpsi1=zpsi4
        zsr1=zsr4
     else
        zpsi2=zpsi3
        zpsi3=zpsi4
        zsr2=zsr3
        zsr3=zsr4
     end if

     ! insert fourth at right
     zsr5=zg1*zsr2+zg2*zsr3
     call spl2d_eval(self%psi,self%n%rcen+zsr5*zcos,self%n%zcen+zsr5*zsin,zpsi5)
     if ((zpsi5-zpsi3)*rsig<0) then
        zpsi2=zpsi5
        zsr2=zsr5
     else
        zpsi1=zpsi3
        zpsi3=zpsi5
        zsr1=zsr3
        zsr3=zsr5
     end if
     ! convergence test
     if (abs(zpsi3-zpsi1)+abs(zpsi3-zpsi2)<peps) exit

  end do do_count

  if (jicount>=micount) then
     kerr=1
  end if

  psim=zpsi3
  psrm=zsr3

end subroutine beq_rextremum
!---------------------------------------------------------------------
!> return rsig to external module
  pure real(kr8) function beq_rsig()

  beq_rsig=rsig

end function beq_rsig
!---------------------------------------------------------------------
!> set up rsig, positive if psi increasing from centre
subroutine beq_rsigset(self)

  type(beq_t), intent(inout) :: self   !< object data structure

  rsig=sign(1._kr8,self%psiqbdry-self%psiaxis)
  !     write(*,*) 'rsig=',rsig

end subroutine beq_rsigset
!---------------------------------------------------------------------
!> calculate contour  \f$ \psi = \psi_{cont} \f$ as array of  \f$ R,Z \f$ values
subroutine beq_psicont(self,pcont,pr,pz)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  real(kr8), intent(in) :: pcont   !< value of \f$ \psi \f$ contour
  real(kr8), dimension(:), allocatable, intent(inout) :: pr !< \f$ R \f$ values of contour
  real(kr8), dimension(:), allocatable, intent(inout) :: pz !< \f$ Z \f$ values of contour

  !! local
  character(*), parameter :: s_name='beq_psicont' !< subroutine name
  real(kr8) :: zpsi    !<  \f$ \psi_{cont} \f$
  real(kr8) :: ztheta    !<  \f$ \theta_j \f$
  real(kr8) :: zrn   !<   \f$ R(\psi_b,\theta_j) \f$
  real(kr8) :: zzn   !<   \f$ Z(\psi_b,\theta_j) \f$

  !! allocate position storage
  if(self%r%n2p>0) then
     allocate(pr(self%r%n2p),pz(self%r%n2p),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No position data')
  end if

  !! loop over angles, evaluating positions

  zpsi=pcont
  do j=1,self%r%n2p
     ztheta=self%r%pos2(j)
     call spl2d_evaln(self%r,zpsi,ztheta,1,zrn)
     call spl2d_evaln(self%z,zpsi,ztheta,2,zzn)
     pr(j)=zrn
     pz(j)=zzn
  end do

  call log_error(m_name,s_name,1,log_info,'Contour values calculated')

end subroutine beq_psicont
!---------------------------------------------------------------------
!> calculate displacements (moves) based on argument and psibdry
subroutine beq_findrzm(self,pxl,pang)

  !! arguments
  type(beq_t), intent(inout) :: self   !< object data structure
  real(kr8), dimension(2) :: pxl !< nominal position vector in \f$ R,Z \f$
  real(kr8), intent(in) :: pang   !< angle \f$ \alpha \f$ made by normal at point

  !! local
  character(*), parameter :: s_name='beq_findrzm' !< subroutine name
  real(kr8) :: zpsi    !<  \f$ \psi_b \f$
  real(kr8) :: ztheta    !<  \f$ \theta_j \f$
  real(kr8) :: zr    !<   \f$ R(\psi_b,\theta_{j-1}) \f$
  real(kr8) :: zz    !<   \f$ Z(\psi_b,\theta_{j-1}) \f$
  real(kr8) :: zrn   !<   \f$ R(\psi_b,\theta_j) \f$
  real(kr8) :: zzn   !<   \f$ Z(\psi_b,\theta_j) \f$
  real(kr8) :: zdotprod    !<  \f$ \Delta_{\bf R} \cdot {\bf n} \f$
  real(kr8) :: zdotprodn    !<  \f$ \Delta_{\bf R} \cdot {\bf n} \f$
  real(kr8) :: zrdiff    !< \f$ \Delta R \f$
  real(kr8) :: zzdiff    !< \f$ \Delta Z \f$
  real(kr8) :: zcos    !<  \f$ \cos(\alpha) \f$
  real(kr8) :: zsin    !<  \f$ \sin(\alpha) \f$

  call log_error(m_name,s_name,1,log_info,'Position data')

  zcos=cos(pang)
  zsin=sin(pang)
  ij=1
  zpsi=self%psibdry

  ztheta=self%r%pos2(ij)
  call spl2d_evaln(self%r,zpsi,ztheta,1,zr)
  call spl2d_evaln(self%z,zpsi,ztheta,2,zz)
  zdotprod=0.
  do j=2,self%r%n2p
     ij=j
     ztheta=self%r%pos2(ij)
     call spl2d_evaln(self%r,zpsi,ztheta,1,zrn)
     call spl2d_evaln(self%z,zpsi,ztheta,2,zzn)
     zrdiff=zrn-zr
     zzdiff=zzn-zz
     zdotprodn=zrdiff*zcos+zzdiff*zsin
     if (zdotprodn*zdotprod<0) then
        self%n%rmove=pxl(1)-zr
        self%n%zmove=pxl(2)-zz
        call log_value("SMITER-GEOQ rmove ",1000*self%n%rmove)
        call log_value("SMITER-GEOQ zmove ",1000*self%n%zmove)
        call log_value("polar theta relative to qcen ",ztheta)
        call log_value("user polar theta ",ztheta+const_pid/2)
     end if
     zr=zrn
     zz=zzn
     zdotprod=zdotprodn
  end do

end subroutine beq_findrzm

subroutine buffer_fixedfmt(kbuff,kfixed,kcfmt)
  !! arguments
  character(*),intent(in) :: kbuff !< name of input file
  integer(ki4),intent(out) :: kfixed  !< flag fixed-format input
  character(*),intent(out) :: kcfmt !< format of buffer

  !! local
  integer(ki4) :: idec1!< positions of first and last decimal points in line
  integer(ki4) :: idec2  !< positions of first and last decimal points in line
  integer(ki4) :: ie1!< positions of first and last e/E  in line
  integer(ki4) :: ie2  !< positions of first and last e/E  in line
  logical, parameter :: lback=.TRUE. !< tests on character strings from far end
  logical :: lmatchdp !< do decimal point positions match

  kfixed=0
  !! find where first and last decimal points are
  idec2=scan(kbuff,'.',lback)
  if (idec2==64.OR.idec2==67.OR.idec2==68) then
     ! quick and dirty check OK, now do detailed examination
     idec1=scan(kbuff,'.')
     lmatchdp= (idec1==3.AND.idec2==67) .OR. (idec1==4.AND.idec2==68) &
 &   .OR. (idec1==4.AND.idec2==64) .OR. (idec1==7.AND.idec2==67)
     if ( lmatchdp ) then
        ! decimals in right place, check exponents
        ie1=scan(kbuff,'e')
        if (ie1==0) ie1=scan(kbuff,'E')
        ie2=scan(kbuff,'e',lback)
        if (ie2==0) ie2=scan(kbuff,'E',lback)
        if  (ie1==13.AND.ie2==77)  then
           ! found correct or almost fixed format
           if (idec1==3) then
              kcfmt='(5e16.9)'
              kfixed=1
           else if (idec1==4) then
              kcfmt='(5(1x,e15.8))'
              kfixed=2
           end if
        else if  (ie1==13.AND.ie2==73)  then
           ! found incorrect fixed format 1
           kcfmt='(5(1x,e14.8e1))'
           kfixed=3
        else if  (ie1==12.AND.ie2==72)  then
           ! found incorrect fixed format 2
           kcfmt='( 5(4x,e11.4) )'
           kfixed=4
        end if
     end if
  end if
  if (kfixed==0) kcfmt='*'

end subroutine buffer_fixedfmt

end module beq_m
