module vcontrol_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use misc_m
  use position_h
  use vcontrol_h
  use position_m
  use scontrol_m
!< access keys defined for sorting

  implicit none
  private

! public subroutines
  public :: &
 &vcontrol_init, &
 &vcontrol_read


! private variables
  character(*), parameter :: m_name='vcontrol_m' !< module name
  integer  :: status   !< error status
  integer  :: nin      !< control file unit number
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: islen !< string length
  integer(ki4), dimension(:), allocatable :: iwork1 !< work array
  integer(ki4), dimension(:), allocatable :: iwork2 !< work array
  character(len=80) :: root !< file root

  contains
!---------------------------------------------------------------------
!> open input control data file
subroutine vcontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='vcontrol_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  character(len=80) :: controlfile !< control file name


  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

  !! open file
  controlfile=trim(fileroot)//".ctl"
  root=fileroot
  call log_value("Control data file",trim(controlfile))
  call misc_getfileunit(nin)
  open(unit=nin,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open control data file')
     stop
  end if


end  subroutine vcontrol_init
!---------------------------------------------------------------------
!> read data for this run
subroutine vcontrol_read(file,numerics)

  !! arguments
  type(vfiles_t), intent(out) :: file !< file names
  type(vnumerics_t), intent(out) :: numerics !< input numerical parameters

  !!local
  character(*), parameter :: s_name='vcontrol_read' !< subroutine name
  integer(ki4), parameter  :: maximum_number_of_files=100 !< maximum number of files allowed
  integer(ki4), parameter  :: maximum_number_of_panels=100 !< maximum number of panels allowed
  character(len=80), dimension(maximum_number_of_files) :: vtk_input_file !< vtk format geometry data input file
  integer(ki4), dimension(maximum_number_of_files) :: number_of_copies !< number of copies to make of input file
  character(len=80), dimension(maximum_number_of_files) :: vtk_label  !< labels vtk format geometry data input file
  integer(ki4), dimension(maximum_number_of_files) :: vtk_label_number !< numeric label for input file (not used)
  character(len=20) :: angle_units  !< units, either radian(s) or degree(s)
  character(len=20) :: option !< for defining panels and their transforms
  logical :: new_controls !< true if vtktfmparameters namelist is present
  character(len=80) :: vtk_output_file  !< vtk format geometry data output file
  logical :: split_file !<  split file by attribute if set
  logical :: make_same !< make attribute take same value
  integer(ki4) :: same_value !<  value for attribute to take
  character(len=80) :: process_by_name !<  name of attribute to be split / homogenised
  integer(ki4) :: number_of_panels !< number of panels for which bodies defined
  integer(ki4) :: number_of_transforms !< number of panels for which transform defined
  integer(ki4) :: max_number_of_files !< number >= no. of files for which transform defined
  integer(ki4) :: max_number_of_panels !< number >= no. of panels for which transform defined
  integer(ki4) :: max_number_of_transforms !< number >= no. of panels for which transform defined
  integer(ki4), dimension(2*maximum_number_of_panels) :: panel_transform !< number of transform to apply
  character(len=20), dimension(maximum_number_of_panels) :: transform_id !< INACTIVE id of transform to apply
  integer(ki4), dimension(2*maximum_number_of_panels) :: panel_bodies !< bodies defining the geometry
  logical :: filefound !< true if file exists
  logical :: lflag !< flags if no panel transforms defined (two places)
  type(tfmdata_t) :: ztfmdata   !< position transform numeric controls
  integer(ki4) :: infil=0 !< number of files detected
  integer(ki4) :: inpan=0 !< number of panels detected
  integer(ki4) :: inbod=0 !< number of bodies detected
  integer(ki4) :: ipan=0 !< panel counter
  integer(ki4) :: intfm=0 !< number of transforms detected
  integer(ki4) :: iflag=0 !< allow position_readcon to return on error
  integer(ki4) :: icode !< body code when copies present
  logical :: ilcopy !< test whether there are copies in the input
  integer(ki4), dimension(2) :: iswap !< swap array
  integer(ki4) :: ierr !< error return code

  logical :: paneltfm !< apply transform if .true.
  integer(ki4) :: max_bods_index !< dimension of bods index array
  integer(ki4) :: max_bods_in_file !< used to generate unique bods numbers over many files
  logical :: preserve_internal_structure !< bods remain distinct
  logical :: extract !< extract objects according to key and limits
  character(len=80) :: extract_key !< key for extraction
  real(kr8), dimension(2) :: plasma_centre !< centre of discharge in \f$ (R,Z) \f$
  real(kr8) :: minimum_angle !< minimum angle for extraction
  real(kr8) :: maximum_angle !< maximum angle for extraction


  !> misc parameters, unusually comes first
  namelist /miscparameters/ &
 &option, new_controls, &
 &max_number_of_files, angle_units, &
 &max_number_of_panels,max_number_of_transforms,&
 &max_bods_index, max_bods_in_file, preserve_internal_structure,&
 &number_of_panels,number_of_transforms

  !> vtktfm parameters
  namelist /vtktfmparameters/ &
 &option, &
 &split_file, make_same, same_value, &
 &process_by_name, &
 &max_number_of_files, angle_units, &
 &max_number_of_panels,max_number_of_transforms,&
 &number_of_panels,number_of_transforms,&
 &paneltfm, extract, extract_key,&
 &max_bods_index, max_bods_in_file, preserve_internal_structure,&
 &plasma_centre, minimum_angle, maximum_angle

  !! file names
  namelist /vtkfiles/ &
 &vtk_input_file, number_of_copies, &
 &vtk_label, vtk_label_number, &
 &vtk_output_file

  !! panelarray parameters
  namelist /panelarrayparameters/ &
 &panel_bodies,panel_transform, &
 &transform_id

  !---------------------------------------------------------------------
  !! set default misc parameters
  new_controls = .false.
  max_number_of_panels = 1
  max_number_of_files = 1
  max_number_of_transforms = 1
  max_bods_index = 1000
  max_bods_in_file = 100
  preserve_internal_structure = .false.
  number_of_panels = 0
  number_of_transforms = 0
  angle_units = 'radians'
  option = 'split_panels'

  !!read misc parameters
  read(nin,nml=miscparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading misc parameters")'
     call log_getunit(ilog)
     write(ilog,nml=miscparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading misc parameters')
  end if

  split_file=.false.
  make_same=.true.
  same_value=1
  process_by_name='Body'
  paneltfm = .true.
  extract = .false.
  extract_key = 'null  '
  plasma_centre =  (/1,0/)
  minimum_angle = -22.5
  maximum_angle = 22.5

  if (new_controls) then
     !!read vtktfm parameters
     read(nin,nml=vtktfmparameters,iostat=status)
     if(status/=0) then
        print '("Fatal error reading vtktfm parameters")'
        call log_getunit(ilog)
        write(ilog,nml=vtktfmparameters)
        call log_error(m_name,s_name,2,error_fatal,'Error reading vtktfm parameters')
     end if
  end if

  !! check for valid data
  if(max_number_of_files<=0) &
 &call log_error(m_name,s_name,2,error_fatal,'max number of files must be >0')
  if(max_number_of_files>maximum_number_of_files) then
     call log_value("maximum_number_of_files parameter",maximum_number_of_files)
     call log_error(m_name,s_name,3,error_fatal,'too many files: increase parameter')
  end if
  !!added for upwards compatibility
  if (number_of_panels>0) max_number_of_panels=number_of_panels
  if(max_number_of_panels<=0) &
 &call log_error(m_name,s_name,11,error_fatal,'max number of panels must be >0')
  if(max_number_of_panels>maximum_number_of_panels) then
     call log_value("maximum_number_of_panels parameter",maximum_number_of_panels)
     call log_error(m_name,s_name,12,error_fatal,'too many panels: increase parameter')
  end if
  if (number_of_transforms>0) max_number_of_transforms=number_of_transforms
  if(max_number_of_transforms<=0) &
 &call log_error(m_name,s_name,13,error_fatal,'max number of transforms must be >0')
  if(angle_units(1:6)/='degree'.AND.angle_units(1:6)/='radian') then
     call log_value("angle_units not recognised",angle_units)
     call log_error(m_name,s_name,14,error_warning,'radians assumed')
  end if
  if(option(1:5)/='split'.AND.option(1:5)/='panel'.AND.option(1:5)/='tagge') then
     call log_value("option not recognised",option)
     call log_error(m_name,s_name,15,error_warning,'split panels assumed')
     option='split'
  end if
  if (extract) then
     if (maximum_angle-minimum_angle<=0) then
        call log_error(m_name,s_name,16,error_fatal,'invalid range for extraction')
     end if
     ! Check key
     numerics%key=extract_key
     call lowor(numerics%key,1,len_trim(extract_key))
     ierr=1
     do j=1,RECOGNISED_KEYS
        if (numerics%key(1:6)==reckey(j)) then
           ierr=0
           exit
        end if
     end do
     if (ierr==1) then
        call log_value("key ", numerics%key)
        call log_error(m_name,s_name,17,error_fatal,'Unrecognised key')
     end if
     numerics%centre  =  plasma_centre
     if(angle_units(1:6)=='degree') then
        numerics%angmin=minimum_angle*const_pid/180
        numerics%angmax=maximum_angle*const_pid/180
     else
        numerics%angmin=minimum_angle
        numerics%angmax=maximum_angle
     end if
  end if
  if(max_bods_index<0) &
 &call log_error(m_name,s_name,18,error_warning,'max size of index must be >=0')
  if(max_bods_in_file<=0) &
 &call log_error(m_name,s_name,19,error_fatal,'max number of bodies in file setting must be >0')

  !---------------------------------------------------------------------
  !! read input file names and associated data
  vtk_input_file='null'
  vtk_output_file='null'
  number_of_copies=1
  vtk_label='null'
  vtk_label_number=0
  read(nin,nml=vtkfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_getunit(ilog)
     write(ilog,nml=vtkfiles)
     call log_error(m_name,s_name,21,error_fatal,'Error reading input filenames')
  end if

  ! check number of files against actuals
  do j=1,maximum_number_of_files
     if (trim(vtk_input_file(j))=='null') exit
     infil=j
  end do
  file%nvtkdata=infil
  call log_value("actual number of files",infil)

  !!check files exist
  do j=1,infil
     call log_value("geometry data file number",j)
     call log_value("geometry data file, vtk_input_file",trim(vtk_input_file(j)))
     if(vtk_input_file(j)/='null') then
        inquire(file=vtk_input_file(j),exist=filefound)
        if(.not.filefound) then
           !! error opening file
           print '("Fatal error: Unable to find field data file, ",a)',vtk_input_file(j)
           call log_error(m_name,s_name,22,error_fatal,'field data file not found')
        end if
        if (number_of_copies(j)<0) then
           call log_error(m_name,s_name,23,error_fatal,'Number of copies is negative')
        end if
     end if
  end do

  !! check parameters for consistency
  if (file%nvtkdata==1) then
     if (make_same.EQV.preserve_internal_structure) then
        call log_error(m_name,s_name,24,error_warning,'make_same and preserve_internal_structure inputs inconsistent')
        preserve_internal_structure=.NOT.make_same
     else if (make_same.AND.split_file) then
        call log_error(m_name,s_name,25,error_warning,'make_same and split_file inputs inconsistent')
        split_file=.false.
     end if
  else
     !! more than one file in input, no sense splitting or extracting
     if (split_file) then
        call log_error(m_name,s_name,26,error_warning,'split_file must be set false')
        split_file=.false.
     end if
     if (extract) then
        call log_error(m_name,s_name,27,error_warning,'extract must be set false')
        extract=.false.
     end if
  end if


  ! allocate and assign files
  allocate(file%vtkdata(infil), file%vtkcopies(infil),&
 &file%vtklabel(infil), stat=status)
  call log_alloc_check(m_name,s_name,30,status)
  do j=1,infil
     file%vtkdata(j)   = vtk_input_file(j)
     file%vtklabel(j)   = vtk_label(j)
     file%vtkcopies(j)   = number_of_copies(j)
  end do

  if(vtk_output_file/='null') then
     !! strip any final '.vtk' string
     islen=len_trim(vtk_output_file)
     if (islen>4) then
        if (vtk_output_file(islen-3:islen)=='.vtk') then
           file%vtkout = vtk_output_file(:islen-4)
        end if
     end if
  else
     file%vtkout = trim(root)
  end if

  !> create output file names from root
  !! output file
  file%vout = trim(root)//"_v.out"
  filefound=.false.
  call log_value("vtktfm output file",trim(file%vout))
  inquire(file=file%vout,exist=filefound)
  if(filefound) then
     !! error opening file
     print '("Warning: output file ",a, " already exists")',trim(file%vout)
     !    print '("Remove it and restart run")'
     call log_error(m_name,s_name,31,error_warning,'Output file already exists')
  end if

  !---------------------------------------------------------------------
  panel_transform=-1
  panel_bodies=0
  transform_id='null'

  ! no need for another namelist
  if (.NOT.paneltfm) return

  !!read panelarray parameters
  read(nin,nml=panelarrayparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading array parameters")'
     call log_getunit(ilog)
     write(ilog,nml=panelarrayparameters)
     call log_error(m_name,s_name,32,error_fatal,'Error reading array parameters')
  end if

  !! check for valid panel data
  ! check number of panels against actuals
  do j=2*maximum_number_of_panels,1,-1
     inbod=j
     if (panel_bodies(j)>0) exit
     inbod=0
  end do
  if (inbod==0) &
 &call log_error(m_name,s_name,33,error_fatal,'No panels present')

  defn_option: select case (option(1:5))
  case ('split')
     inpan=(inbod+1)/2
     call log_value("actual number of panels",inpan)
     ipan=0
     numerics%npans=inpan
     do j=1,inbod,2
        ipan=ipan+1
        if (panel_transform(j)<0.AND.transform_id(ipan)=='null') then
           call log_value("panel number",ipan)
           call log_error(m_name,s_name,34,error_fatal,'Transform not defined')
        end if
        if (panel_bodies(j+1)<=0) panel_bodies(j+1)=panel_bodies(j)
     end do
     inbod=2*inpan ! make sure even

  case ('panel')
     inpan=inbod
     call log_value("actual number of panels",inpan)
     if (inpan>maximum_number_of_panels) &
 &   call log_error(m_name,s_name,36,error_fatal,'maximum_number_of_panels too small')
     numerics%npans=inpan
     do j=1,inpan
        if (panel_transform(j)<0.AND.transform_id(j)=='null') then
           call log_value("panel number",j)
           call log_error(m_name,s_name,37,error_fatal,'Transform not defined')
        end if
     end do
     ! double up panel bodies
     do j=inpan,1,-1
        panel_bodies(2*j)=panel_bodies(j)
        panel_bodies(2*j-1)=panel_bodies(j)
     end do
     inbod=2*inpan

  case ('tagge')
     inpan=(inbod+1)/2
     call log_value("actual number of panels",inpan)
     numerics%npans=inpan
     ipan=0
     lflag=(maxval(panel_transform)>0)
     do j=1,inbod,2
        ipan=ipan+1
        if (panel_bodies(j)<0) then
           call log_value("panel number",ipan)
           call log_error(m_name,s_name,41,error_fatal,'bodies not defined')
        end if
        if (panel_bodies(j+1)<0) then
           call log_value("panel number",ipan)
           call log_error(m_name,s_name,42,error_fatal,'Transform tag must be non-negative')
        end if
        if (lflag) then
           if (panel_transform(j)<0) then
              call log_value("panel number",ipan)
              call log_error(m_name,s_name,43,error_fatal,'Transform tag must be non-negative')
           end if
           if (panel_transform(j+1)<0.AND.transform_id(ipan)=='null') then
              call log_value("panel number",ipan)
              call log_error(m_name,s_name,44,error_fatal,'Transform not defined')
           end if
        else
           ! define panel_transform using body tags
           panel_transform(j)=panel_bodies(j+1)
           panel_transform(j+1)=panel_bodies(j+1)
        end if
     end do
     ! resolve tags
     allocate(iwork1(2*maximum_number_of_panels),iwork2(2*maximum_number_of_panels),stat=status)
     iwork1=0
     iwork2=-1
     call log_alloc_check(m_name,s_name,45,status)
     ipan=0
     do j=1,inbod,2
        ipan=ipan+1
        iwork1(ipan)=panel_bodies(j)
        do i=1,inbod,2
           if (panel_transform(i)==panel_bodies(j+1)) then
              iwork2(ipan)=panel_transform(i+1)
              exit
           end if
        end do
     end do
     ! test for error
     if (minval(iwork2(:inpan))<0) then
        call log_error(m_name,s_name,46,error_fatal,'Transforms not fully defined')
     end if
     ! double up panel bodies
     do i=inpan,1,-1
        panel_bodies(2*i)=iwork1(i)
        panel_bodies(2*i-1)=iwork1(i)
     end do
     panel_transform(:inpan)=iwork2(:inpan)
     inbod=2*inpan ! make sure even
     deallocate(iwork1,iwork2)
  end select defn_option

  !! check number of copies matches, i.e. if there are 10 copies of body 5, then
  !! bodies with codes 501...510 must appear in panel_bodies

  !write(*,*) inbod, (panel_bodies(l), l=1,inbod)
  ilcopy=.false.
  if (maxval(panel_bodies(:inbod))>100) then
     ! This implies copies exist, so bodies are not numbered contiguously
     ! so trick below to generate virtual bodies will not work
     ilcopy=.true.
     do j=1,infil
        do k=1,number_of_copies(j)
           icode=j*100+k

           do l=1,inbod
              lflag=(panel_bodies(l)==icode)
              if (lflag) exit
           end do
           if (.NOT.lflag) then
              call log_value("panel code number",icode)
              call log_error(m_name,s_name,49,error_fatal,'No matching body')
           end if

        end do
     end do
  end if
  !! allocate and assign
  !! panel transform parameters
  allocate(numerics%pantfm(maximum_number_of_panels), stat=status)
  call log_alloc_check(m_name,s_name,50,status)
  !! panel bodies
  allocate(numerics%panbod(2,maximum_number_of_panels), stat=status)
  call log_alloc_check(m_name,s_name,51,status)
  ! default
  numerics%pantfm=0
  ! redefine only those in input
  numerics%pantfm(:inpan)=panel_transform(:inpan)
  if (ilcopy) then
     do i=1,inpan
        numerics%panbod(1,i)=panel_bodies(2*i-1)
        numerics%panbod(2,i)=panel_bodies(2*i)
     end do
  else
     do i=1,maximum_number_of_panels
        numerics%panbod(:,i)=i
     end do
     do j=1,inpan
        do i=1,maximum_number_of_panels
           if (numerics%panbod(1,i)==panel_bodies(2*j-1)) then
              ! swap to make sure all panbod still defined
              iswap=numerics%panbod(:,j)
              numerics%panbod(1,j)=panel_bodies(2*j-1)
              numerics%panbod(2,j)=panel_bodies(2*j)
              numerics%panbod(:,i)=iswap
              exit
           end if
        end do
     end do
  end if

  !---------------------------------------------------------------------
  ! define transformations
  allocate(numerics%vptfm%scale(3,max_number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,55,status)
  allocate(numerics%vptfm%offset(3,max_number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,56,status)
  allocate(numerics%vptfm%matrix(3,3,max_number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,57,status)
  allocate(numerics%vptfm%ntfm(max_number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,58,status)
  !ID  allocate(numerics%vptfm%id(max_number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,59,status)
  !! transform array parameters
  do j=1,max_number_of_transforms
     call log_value('read namelist',j)
     call position_readcon(ztfmdata,nin,iflag)
     if (iflag/=0) exit
     intfm=j
     numerics%vptfm%scale(:,j)=ztfmdata%scale
     numerics%vptfm%offset(:,j)=ztfmdata%offset
     !F11 write(*,*) 'ztfmdata%offset',ztfmdata%offset
     numerics%vptfm%matrix(:,:,j)=ztfmdata%matrix
     numerics%vptfm%ntfm(j)=ztfmdata%ntfm
     !ID     numerics%vptfm%id(j)=ztfmdata%id
  end do
  call log_value('Return flag',iflag)
  call log_value('Number of transforms present',intfm)

  ! assign transform numbers
  do j=1,inpan
     if (transform_id(j)(1:4)=='unit') then
        numerics%pantfm(j)=0
        cycle
     end if
     !ID     do k=1,intfm
     !ID        if (transform_id(j)==numerics%vptfm%id(k)) then
     !ID           numerics%pantfm(j)=numerics%vptfm%ntfm(k)
     !ID           exit
     !ID        end if
     !ID     end do
  end do
  !! final check
  do j=1,inpan
     if (numerics%pantfm(j)<0.OR.&
 &   numerics%pantfm(j)>intfm) then
        call log_value("panel number",j)
        call log_error(m_name,s_name,60,error_fatal,'Transform not defined')
     end if
  end do
  numerics%vptfm%ntfmtyp=intfm
  !F11 write(*,*) 'numerics%pantfm(j)',(numerics%pantfm(j),j=1,inpan)

  numerics%paneltfm =  paneltfm
  numerics%maxindx =  max_bods_index
  numerics%maxbodsf =  max_bods_in_file
  numerics%preserve =  preserve_internal_structure
  numerics%extract  =  extract
  numerics%angles=angle_units(1:6)
  numerics%split=split_file
  numerics%same=make_same
  numerics%nvalue=same_value
  numerics%name=process_by_name

end subroutine vcontrol_read

end module vcontrol_m
