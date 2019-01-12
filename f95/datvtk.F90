program datvtk_p

  use const_kind_m
  use const_numphys_h
  use date_time_m
  use control_h
  use dcontrol_h
  use log_m
  use misc_m
  use clock_m
  use spl2d_m
  use spl3d_m
  use position_h
  use fmesh_h
  use beq_h
  use posang_h
  use posang_m

  use pcle_h
  use pcle_m
  use geobjlist_h
  use position_m
  use li_m
  use ld_m
  use ls_m
  use dbtree_h
  use dbtree_m
  use btree_m
  use geobj_m
  use query_m
  use geobjlist_m
  use indict_m

  use dcontrol_m
  use vcontrol_h
  use vfile_m
  use gfile_m
  use datline_h
  use datline_m
  use dfile_m
  use bods_h
  use bods_m
  use stlfile_m
  use stack_m

  implicit none


! Local variables
  character(*), parameter :: m_name='datvtk' !< module name
  type(dfiles_t)     :: file      !< names of files
  type(dnumerics_t)  :: numerics  !< numerical control parameters
  type(geobjlist_t)  :: geobjl      !< geometrical objects
  type(geobjlist_t)  :: geobj2      !< geometrical objects
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=20) :: buf1 !< buffer input
  character(len=20) :: buf2 !< buffer input
  character(len=3) :: optarg='nxx' !< optional argument
  character(len=256) :: vtkdesc !< descriptor line for vtk files
  character(len=80),save :: iched !< vtk field file descriptor

  type(bods_t) :: bods !< type of bods for geometrical objects
  integer:: nplot !< unit for vtk files
  integer:: nread !< unit for dat files
  integer:: nin !< unit for other data
  integer(ki4):: i !< loop variable
  integer(ki4):: j !< loop variable
  integer(ki4) :: islen   !< length of input field filename
  integer(ki4), dimension(2):: idum !< dummy array
  integer(ki4) :: inarg   !< number of command line arguments
  integer(ki4) :: inda   !< index of dash in command line argument
  integer(ki4) :: iopt   !< option for reading bods
!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(30)
  call clock_start(1,'datvtk run time')
!--------------------------------------------------------------------------
!! print header

  print *, '----------------------------------------------------'
  print *, 'datvtk: convert dat to vtk geometry file'
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long
!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run datvtk type at the command line:'
     print *, '   datvtk [-opt] fileroot '
     stop
  else
!!get fileroot
     inarg=command_argument_count()
     if (inarg==2)  then
        call get_command_argument(1,value=buf1)
! strip leading dash
        inda=index(buf1,'-')
        buf2=adjustl(buf1(inda+1:))
        i=min(len_trim(buf2),3)
        if (i>0) then
           optarg(1:i)=buf2
        else
!! no option specified after "-"
           print *, 'Fatal error: no option specified after "-".'
           print *, 'To run datvtk type at the command line:'
           print *, '   datvtk [-opt] fileroot '
           stop
        end if
     end if
     call get_command_argument(inarg,value=fileroot)
!! strip any final '.dat' or '.vtk' string
     islen=len_trim(fileroot)
     if (islen>4) then
        if (fileroot(islen-3:islen)=='.dat') then
           fileroot(islen-3:islen)='    '
        else if (fileroot(islen-3:islen)=='.vtk') then
           fileroot(islen-3:islen)='    '
        end if
     end if
  end if

!! start log
  call log_init(fileroot,timestamp)
  call log_value("option string for datvtk",optarg)
!  write(*,*) 'this output helps gfortran sometimes', fileroot
!--------------------------------------------------------------------------
!! read control file in case of silhouette data

  if (optarg(1:1)=='c') then
     call clock_start(2,'dcontrol_init time')
     call dcontrol_init(fileroot)
     call dcontrol_read(file,numerics)
     call clock_stop(2)
  end if
!--------------------------------------------------------------------------
!! do the main work
!! read  geobjl data

  call clock_start(4,'geobjlist_(d)read time')
  input_type: select case(optarg(1:1))
  case('v')
     call geobjlist_read(geobjl,trim(fileroot)//'.vtk',iched)
     inquire(file=trim(fileroot)//'.vtk',number=nread)
     iopt=1
     call vfile_iscalarread(bods%list,geobjl%ng,trim(fileroot)//'.vtk','Body',nread,iopt)
  case('c')
     iched='converted dat files in '//trim(fileroot)//'.ctl'
     call geobjlist_create3d(geobjl,numerics,0_ki2par)
     call bods_initlist(bods,geobjl,1)
  case default
     iched='converted dat file '//trim(fileroot)//'.dat'
     call dfile_init(fileroot,nread)
     call geobjlist_dread(geobjl,nread)
     call bods_initlist(bods,geobjl,1)
  end select input_type
  call clock_stop(4)
!--------------------------------------------------------------------------
!! shell/reorient  geobjl data as needed

  call clock_start(5,'geobjlist_orientri time')
  transform_type: select case(optarg(2:2))
  case('o')
     call geobjlist_orientri(geobjl)
  case('f')
     call geobjlist_fliptri(geobjl)
  case('s')
     call geobjlist_shelltets(geobjl,geobj2)
     call geobjlist_delete(geobjl)
     call geobjlist_copy(geobj2,geobjl,1)
     call geobjlist_delete(geobj2)
     call geobjlist_orientri(geobjl)
  end select transform_type
  call clock_stop(5)
!--------------------------------------------------------------------------
!! output file

  call clock_start(30,'outfile_init time')
!     write(*,*) 'iched=',iched
  divide_type: select case(optarg(3:3))
  case('d','s')
!! write out as separate files
     call bods_write(bods,geobjl,fileroot,'none','Body',1)
     call vfile_close
  case('t')
! stl, one set of triangles
     call stlfile_init(trim(fileroot),'triangles',nplot)
     call geobjlist_writestl(geobjl,'triangles',nplot)
     call stlfile_close
  case('b')
! stl, split by body number (TO DO)
     call stlfile_init(trim(fileroot),'bodies',nplot)
     call geobjlist_writestl(geobjl,'bodies',nplot)
     call stlfile_close
  case default
     call geobjlist_makehedline(geobjl,iched,vtkdesc)
     call vfile_init(trim(fileroot)//'_out',vtkdesc,nplot)
     call geobjlist_writev(geobjl,'geometry',nplot)
     if (optarg(3:3)=='m') bods%list=1 ! suppress body information
     call vfile_iscalarwrite(bods%list,bods%nbod,'Body','CELL',nplot,1)
     call vfile_close
  end select divide_type
  call clock_stop(30)
!--------------------------------------------------------------------------
!! cleanup and closedown
  call geobjlist_delete(geobjl)
  call bods_delete(bods)

  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete
!--------------------------------------------------------------------------

end program datvtk_p
