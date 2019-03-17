module smanal_h

  use const_kind_m
  use geobjlist_h
  use position_h
  use geobj_m
  use scontrol_h

! public types
!> analysis object
  type, public :: smanal_t

     type(geobjlist_t) :: geobjl !< geometrical objects
     type(snumerics_t)  :: n  !< numerical parameters
     integer(ki4), dimension(:), allocatable :: key !< integer key denoting e.g. array of bodies
     real(kr8), dimension(:), allocatable :: rkey !< real sort key denoting e.g heights of bodies
     integer(ki4), dimension(:), allocatable :: mkey !< array of more bodies for geometrical objects
     integer(ki4), dimension(:), allocatable :: dict !< array of with all body entries different
     integer(ki4), dimension(:), allocatable :: oldict !< keep array of with all body entries different
     real(kr8), dimension(:), allocatable :: rdict !< array setting interval boundaries for real sort keys
     integer(ki4) :: nkey !< number of body entries
     integer(ki4) :: ndict !< number of different bodies
     integer(ki4) :: nrdict !< number of different entries in real dictionary
     real(kr8), dimension(:), allocatable :: scal !< scalar to be analysed
     integer(ki4) :: nscal !< number of scalar entries
     real(kr8), dimension(:,:), allocatable :: stat !< statistics
     integer(ki4), dimension(:), allocatable :: nstat !< number of statistics in each entry of stat
     character(len=80) :: formula !< smanal formula
     real(kr8) :: f !< power split (math variable name allowed)
     integer(ki4) :: nrpams !< number of real parameters
     integer(ki4) :: nipams !< number of integer parameters
     real(kr8), dimension(:,:), allocatable   :: centr !< centroids of bodies
     real(kr4), dimension(:), allocatable   :: areab !< areas of bodies
     real(kr8), dimension(:), allocatable   :: bin !< rekey bins
     integer(ki4) :: inbin !< number of entries in bin
     logical  :: larea !< whether areas of geometrical objects calculated
     real(kr4), dimension(:), allocatable   :: area !< areas of geometrical objects
     integer(ki4), dimension(:), allocatable   :: newkey !< new set of keys, e.g. based on position
     real(kr8), dimension(:), allocatable   :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable   :: npar !< general integer parameters

  end type smanal_t

end module smanal_h
