module dbtree_h

  use const_kind_m
  use position_h
  use ls_m
  use li_m
  use ld_m

! public types

!> parameters describing how to construct object
  type, public :: dbnumerics_t
     integer(ki4) :: nsize !< size of binary tree
     integer(ki4) :: nsizep !< size of binary tree pter
     integer(ki4) :: nsizee !< size of exten array
     integer(ki4) :: nsizeh !< size of linked list hoc
     integer(ki4) :: nsizel !< size of linked list
     integer(ki2) :: ndepth !< max depth of binary tree
     integer(ki2) :: nttype !< type of binary tree
     integer(ki4) :: nttalg !< type of binary tree algorithm for top node
     integer(ki4) :: ntlist !< type of list structure used 0-ls,1-li,2-ld
     integer(ki2), dimension(3) :: nxyz !< dimensions of top node (nttype=2)
     real(kr4), dimension(3) :: hxyz !< specified mesh sizes
     integer(ki4) :: ngeobj !< number of geobj records to be read
     integer(ki4) :: mtype !< type of margin for quantisation
     integer(ki4) :: geomtype !< geobj type
     integer(ki4) :: maxinbin !< maximum number of geobj records in a bin
     integer(ki2) :: nquante !< estimated log number of mesh cells
     real(kr4) :: maxtolerance !< max max distance of geobj from its side
     real(kr4) :: mintolerance !< min max distance of geobj from its side
     integer(ki4) :: nbdcub !< number of extra bounding cubes
     integer(ki4) :: cornflag !< flags (=1) if corner vectors to be used
     real(kr4), dimension(3) :: lowcorner !< vector defining corner with lesser component values
     real(kr4), dimension(3) :: upcorner !< vector defining corner with greater component values
     real(kr4) :: dilen   !< inner cube separation from geometry
     real(kr4) :: dolen   !< outer cube separation from inner
     type(tfmdata_t) :: xtfm !< position \f$ x \f$ to \f$ x \f$ scaling
     integer(ki4) :: nqtfm !< type of quantum vector transform
     integer(ki4) :: splitalg !< decide BSP splitting  (1=most asymmetric, 3=least)
  end type dbnumerics_t

!> type which defines/instantiates the object
  type, public :: dbtree_t
     integer(ki4)  :: nt !< number of entries in tree
     integer(ki4), dimension(:,:), allocatable :: pter !< tree links
     integer(ki2), dimension(:,:), allocatable :: desc !< description
     integer(ki2), dimension(:,:), allocatable :: corner !< scaled position
     integer(ki2), dimension(:,:), allocatable :: exten !< array of exten descriptions
     integer(ki4) :: npter !< number of entries in pter
     integer(ki4) :: nexten !< number of entries in exten
     integer(ki4) :: nroot !< number of root record
     integer(ki4) :: nleaf !< number of leaves
     type(dbnumerics_t) :: n !< control  parameters
     type(ls_t) :: objectls !< numbered list object for children
     type(li_t) :: objectli !< simple linked list object for children
     type(ld_t) :: objectld !< doubly linked list object for children
     type(quantfm_t) :: quantfmq !< geobj \f$ x \f$ between quantised mesh units scaling
     type(quantfm_t) :: quantfm !< geobj \f$ x \f$ to mesh units scaling
     real(kr4), dimension(3,2) :: binbb !< bb for geobj binning
     integer(ki4) :: maxchildn !< maximum number of children of a node
     integer(ki4) :: minbdim !< log_2 (smallest box dimension in tree)
     integer(ki4) :: maxallb !< actual maximum number of objects in any node
     real(kr8), dimension(:), allocatable :: scal !< array with entry for each leaf
     real(kr8), dimension(:), allocatable :: scal2 !< 2nd array with entry for each leaf
     real(kr8), dimension(:), allocatable :: geom !< array with entry for each leaf
  end type dbtree_t

end module dbtree_h
