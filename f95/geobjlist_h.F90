module geobjlist_h

  use const_kind_m
  use position_h
  use geobj_m


! public types

!> type storing geobj coord data and useful bits
  type, public :: geobjlist_t
     integer(ki4)  :: ngtype !< type of geobj (if all the same)
     type(posveclis_t) :: posl !< list of positions
     integer(ki4), dimension(:), allocatable :: nodl !< list of nodes
     integer(ki4)  :: np !< number of position vectors
     type(geobj1_t), dimension(:), allocatable :: obj !< type variable
     type(geobj2_t), dimension(:), allocatable :: obj2 !< type variable
     integer(ki4)  :: nnod !< no of nodes
     integer(ki4)  :: ng !< no of geobjs
     integer(ki4)  :: ngunassigned !< no of geobjs not assigned
     integer(ki4)  :: nsampl !< no of geobjs in sample
     integer(ki2)  :: nquant !< estimated then actual quantising number (log)
     integer(ki4)  :: minobj !< minimum number of geobj in a bin
     integer(ki2par)  :: numnparam=2 !< number of entries in nparam
   !> Integer parameter(s) of vtk file header
   !! Array index entries give
   !! 1. infilelevel (level of refinement in file), default 1
   !! 2. indicates triangles geometry coded as to type
   !! (so labelled scalar Code in file, default 0, no type code)
     integer(ki2par), dimension(2)  :: nparam =(/1,0/) !< .
     integer(ki2par)  :: numrparam=1 !< number of entries in rparam
     real(kr4), dimension(1) :: rparam=(/0._kr4/) !< real parameter(s) from file header (NOT USED)
     real(kr4), dimension(3,2) :: coordbb !< bb of geobj coord
     real(kr4), dimension(3,2) :: binbb !< bb for geobj binning
     real(kr4) :: tolerance !< max distance from face
     real(kr4) :: minmaxtolerance !< min max distance from face
     integer(ki4) :: nbdcub=0 !< number of extra bounding cubes
     integer(ki4) :: cornflag !< flags (=1) if corner vectors to be used
     real(kr4), dimension(3) :: lowcorner !< vector defining corner with lesser component values
     real(kr4), dimension(3) :: upcorner !< vector defining corner with greater component values
     real(kr4) :: dilen   !< inner cube separation from geometry
     real(kr4) :: dolen   !< outer cube separation from inner
     integer(ki4)  :: nwset !< nonzero if obj%weight set
     type(quantfm_t)    :: quantfm !< geobj \f$ x \f$ to mesh units  scaling
     type(tfmdata_t)    :: tfmdata !< position \f$ x \f$ to \f$ x \f$ scaling
  end type geobjlist_t

end module geobjlist_h
