module namelist

  use const_kind_m

! public types
  integer(ki4), parameter :: MAX_NUMBER_OF_POSITIONS=100 !< maximum number of positions allowed
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  integer(ki4), parameter :: maxinp=100 !< Maximum number of inputs allowed for each variable
  integer(ki4), parameter  :: maximum_number_of_tracks=100 !< maximum number of tracks allowed
  integer(ki4), parameter  :: maximum_number_of_files=100 !< maximum number of files allowed
  integer(ki4), parameter  :: maximum_number_of_panels=100 !< maximum number of panels allowed
!---------------------------------------------------------------------
!> Entries in namelist analysisfiles
  type, public :: analysisfiles
     character(len=80), dimension(maximum_number_of_files) :: vtk_data_file !< vtk format geometry data input file, default 'null'
     character(len=80) :: vtk_input_file  !< smanal data input file _powx format, default 'null'
     character(len=80) :: gnu_output_file  !< smanal gnu/csv data output file, default 'null'
     character(len=80) :: gnure_output_file  !< smanal gnu/csv data output file, default 'null'
     character(len=80) :: vtk_small_output_file  !< smanal data output at centroids file, default 'null'
     character(len=80) :: vtk_output_file  !< smanal data output file, default 'null'
  end type analysisfiles
!---------------------------------------------------------------------
!> Entries in namelist beqparameters.
!! - Three different magnetic field configurations are recognised,
!! selected using parameter (beq_)fldspec,
!! corresponding to axisymmetry in either the flux-mapped case (1) or
!! cylindrical case (3), or nonaxisymmetry in cylindrical coordinates (2).
!! There is also the option to declare a duct calculation and specify
!! the field on a Cartesian grid.
!! - In the non-axisymmetric case, the geometry may be assumed to be
!! repeated in toroidal angle, e.g. \f$  N_{\zeta p} =18 \f$ repeats for ITER,
!! or M repeats due to the presence of a field ripple
!! with frequency M, \f$ B_R = - a_r R^{M-1} \sin M \zeta \f$
!! and corresponding \f$ B_\zeta \f$. The geometry may also be assumed
!! repeated in the axisymmetric cylindrical case.
!! - If an option is input which is not given below, then the code will stop.
!! - Coordinate directions \f$ 1,2,3 \f$ correspond to respectively
!! \f$ (\psi,\theta,\zeta) \f$ in case of flux coordinates.
!! - Coordinate directions \f$ 1,2,3 \f$ correspond to respectively
!! \f$ (R,Z,\xi) \f$ (N.B. \f$ \xi \f$ not \f$ \zeta \f$) in cylindrical case.
!! - \f$ \xi \f$ is \f$ \zeta \f$ scaled to  a region \f$ 2\pi \f$ in extent
!! \f$ \xi = N_{\zeta s} \zeta \f$, where \f$ N_{\zeta s} =N_{\zeta p} \f$
!! or \f$ M \f$ if there is a field ripple present.
  type, public :: beqparameters
   !> \f$ C_{opt}=1\f$ (default). Set position of plasma centre \f$ R_{cen}, Z_{cen} \f$,
   !! which becomes origin of flux-mapped coordinate system, and also sets \f$ Z_{cen} \f$
   !! as origin for determining flux expansion and hence deposited power levels.
   !! 1. User input of \f$ R_{cen}, Z_{cen} \f$
   !! 2. Find mesh-point with extremal psi value, start search at
   !!  \f$ (R_{\min}+R_{\max})/2, (Z_{\min}+Z_{\max})/2 \f$
   !! 3. Find mesh-point with extremal psi value, start search at
   !!  user-specified position \f$ R_{cen}, Z_{cen} \f$
   !! 4. Use values \f$ R_{qcen}, Z_{qcen} \f$ from EQDSK file
     integer(ki4) :: beq_cenopt !< .
   !> \f$ \psi_{opt} =1\f$ (default). Set range of flux to be used in flux mapping (if flux-mapped case)
   !! 1. User input of normalised \f$ \psi_{\min}\f$ and \f$ \psi_{\max}\f$ for mapping,
   !! normalised by X-point or limiter value, no margin \f$ \Delta \psi\f$.
   !! 2. \f$ \psi_{\min}\f$ and \f$ \psi_{\max}\f$ are calculated from input geometry,
   !! with margin \f$ \delta \psi = \Delta \psi  \times (\psi{\max}-\psi_{\min}) \f$
   !! 3. User input of \f$ \psi_{\min}\f$ and \f$ \psi_{\max}\f$ for mapping.
     integer(ki4) :: beq_psiopt !< .
   !> number of points for gnuplot file _cart.gnu
   !! display in direction of invariant coordinate \f$ \zeta \f$,
   !! default 193, if equal to unity, output plane at \f$ \zeta=0 \f$
     integer(ki4) :: mzeta_gnup !< .
   !> number of points for vtk file _geoqvolx.vtk
   !! display in direction of invariant coordinate \f$ \zeta \f$,
   !! default equal to unity, output plane at \f$ \zeta=0 \f$
     integer(ki4) :: mzeta_vtk !< .
   !> \f$ m_{opt} =1\f$ (default). Determine how to calculate value of flux at plasma boundary \f$ \psi_m \f$,
   !! and where to evaluate values of \f$ R_m, B_{pm} \f$ (and \f$ B_m \f$)
   !! used in power deposition formulae.
   !! If a limiter discharge, options 1-3 imply \f$ R_m, B_{pm} \f$ is evaluated
   !! at the contact point, otherwise they are evaluated on the midplane (through the
   !! plasma centre) at the point where \f$ \psi=\psi_m \f$, either
   !! inboard/HFS (options 4,5) or outboard/LFS (options 8,9).
   !! rsig is important since rsig>0(<0)  implies the plasma centre is at a minimum(maximum) of \f$ \psi \f$
   !! and the resulting handedness of the field-line "helices" matters to SMARDDA-PFC, although
   !! field direction does not. Below, \f$ \psi_X \f$ is the value of flux at the X-point which
   !! geoq has found to be nearest in physical space to the plasma centre.
   !! 1. User input of \f$ \psi_m=\psi_{ref} \f$
   !! 2. Use value \f$  \psi_{qbdry} \f$  from EQDSK file unless \f$ \psi_{ltr} \f$ is tighter
   !! Tighter implies \f$ \psi_m  = \max(\psi_{ltr}, \psi_{ref}) \f$ if rsig<0, vice versa if rsig>0
   !! 3. \f$ \psi_m=\psi_{ltr} \f$ is calculated using a hardwired spatial rate of at least deltal=10mm
   !! (to allow for use of silhouette).
   !! 4. \f$ \psi_m=\psi_X \f$ is used.
   !! Inboard values used for \f$ R_m, B_{pm}, B_m \f$
   !! 5. \f$ \psi_m=\psi_{ref} \f$ with inboard values used for \f$ R_m, B_{pm}, B_m \f$
   !! 6. \f$ \psi_m=\psi_{ltr} \f$ as 3 unless \f$ \psi_{qbdry} \f$ is tighter.
   !! Tighter implies \f$ \psi_m  = \max(\psi_{ltr}, \psi_{qbdry}) \f$ if rsig<0, vice versa if rsig>0
   !! 7. \f$ \psi_m=\psi_X \f$ is used unless \f$ \psi_{qbdry} \f$  is tighter.
   !! Tighter implies \f$ \psi_m  = \max(\psi_X, \psi_{qbdry}) \f$ if rsig<0, vice versa if rsig>0.
   !! Inboard values used for \f$ R_m, B_{pm}, B_m \f$
   !! 8. \f$ \psi_m=\psi_X \f$ is used.
   !! Outboard values used for \f$ R_m, B_{pm}, B_m \f$
   !! 9. \f$ \psi_m=\psi_{ref} \f$ with outboard values used for \f$ R_m, B_{pm}, B_m \f$
   !! 10. as 7, but with outboard values used for \f$ R_m, B_{pm}, B_m \f$
   !! 11. \f$ \psi_m=\psi_{qbdry} \f$  from EQDSK file regardless, inboard values
   !! used for \f$ R_m, B_{pm}, B_m \f$
   !! 12. \f$  \psi_m=\psi_{qbdry} \f$  from EQDSK file regardless, outboard values
   !! used for \f$ R_m, B_{pm}, B_m \f$
   !! 13. \f$  \psi_m=\psi_{ltr} \f$  boundary based on nodes of geometry only
   !! contact point values used for \f$ R_m, B_{pm}, B_m \f$
   !! 14. \f$  \psi_m=\psi_{ltr} \f$  boundary based on nodes of geometry only, inboard values
   !! used for \f$ R_m, B_{pm}, B_m \f$
   !! 15. \f$  \psi_m=\psi_{ltr} \f$  boundary based on nodes of geometry only, outboard values
   !! used for \f$ R_m, B_{pm}, B_m \f$
     integer(ki4) :: beq_bdryopt !< .
   !> \f$ N_{opt}=1 \f$ (default). Set numbers of points  \f$ N_{\psi}, N_{\theta} \f$ to be used in flux mapping
   !! 1. User input of \f$ N_{\psi}\f$ and \f$ N_{\theta}\f$  or defaults
   !! 2. Estimate sensible values of \f$ N_{\psi}, N_{\theta} \f$ proportional
   !! to \f$ N_R, N_Z \f$ from equilibrium file
     integer(ki4) :: beq_nopt !< .
   !> \f$ \theta_{opt}=1\f$ (default). Set range of \f$ \theta \f$ to be used in flux mapping.
   !! Origin of theta may be set using \f$ \theta_X \f$
   !! N.B. Full range may not be realised if there are local extrema of flux
   !! in the mapped region.
   !! 1. User input of \f$ \theta_{\min}\f$ and \f$ \theta_{\max}\f$ in coordinate
   !! system with zero at X-point, or relative to vertical, default unit radians,
   !! outboard is positive, no margin \f$ \delta \theta\f$.
   !! 2. \f$ \theta_{\min}\f$ and \f$ \theta_{\max}\f$ are calculated from input geometry,
   !! with margin \f$ \delta \theta =  \Delta \theta  \times (\theta_{\max}-\theta_{\min}) \f$
     integer(ki4) :: beq_thetaopt !< .
   !> \f$ \zeta_{opt} =3\f$ (default). Set range of \f$ \zeta \f$ to be used in mapping to cylindricals,
   !! note that these are only used if \f$ \xi_{opt}=1 \f$, otherwise
   !! \f$ \zeta_{opt} \f$ is effectively inert.
   !! 1. User input of \f$ \zeta_{\min}\f$ and \f$ \zeta_{\max}\f$, must be
   !! positive angles, default unit radians. Note \f$ \zeta = -\phi\f$.
   !! 2. \f$ \zeta_{\min}\f$ and \f$ \zeta_{\max}\f$ are chosen respectively as
   !! \f$ -\pi/m \f$ and \f$  +\pi/m \f$
   !! using m=\f$ N_{\zeta p} \f$
   !!  (if \f$ M_{\zeta p} \f$>0), or m=M (if M>0) or m=1 (if \f$ M_{\zeta p} \f$=M=0).
   !! 3. \f$ \zeta_{\min}\f$ and \f$ \zeta_{\max}\f$ are chosen respectively as
   !! \f$ -\pi \f$ and \f$  +\pi \f$ regardless
   !! 4. \f$ \zeta_{\min}\f$ and \f$ \zeta_{\max}\f$ are chosen respectively as
   !! \f$ -3\pi/2 \f$ and \f$  +\pi/2 \f$ regardless
   !! 5. \f$ \zeta_{\min}\f$ and \f$ \zeta_{\max}\f$ are chosen respectively as
   !! \f$ -\pi/2 \f$ and \f$  +3\pi/2 \f$ regardless
     integer(ki4) :: beq_zetaopt !< .
   !> \f$ \xi_{opt} =2\f$ (default). Set range of \f$ \xi \f$ to be used in mapping to cylindricals.
   !! The point is that toroidal phase matters if there is a ripple field.
   !! Otherwise, \f$ \xi_{opt}=2 \f$ is the best option, and if a segment only of
   !! geometry is to be treated, the segment should be centred on the \f$ X \f$ axis.
   !! 1. \f$ \xi_{\min}\f$ and \f$ \xi_{\max} =  \zeta_{\min}\f$ and \f$ \zeta_{\max}\f$
   !! respectively.
   !! 2. \f$ \xi_{\min}\f$ and \f$ \xi_{\max}\f$ are chosen respectively as
   !! \f$ -\pi \f$ and \f$  +\pi \f$, also  \f$ N_{\zeta s} = N_{\zeta p} \f$ (if \f$ N_{\zeta p}>0 \f$) or \f$ N_{\zeta s} = M \f$
   !! 3. Invalid
   !! 4. \f$ \xi_{\min}\f$ and \f$ \xi_{\max}\f$ are chosen respectively as
   !! \f$ -3\pi/2 \f$ and \f$  +\pi/2 \f$ regardless, also  \f$ N_{\zeta s} = N_{\zeta p} \f$
   !!  (if \f$ N_{\zeta p}>0 \f$) or \f$ N_{\zeta s} = M \f$
   !! 5. \f$ \xi_{\min}\f$ and \f$ \xi_{\max}\f$ are chosen respectively as
   !! \f$ -\pi/2 \f$ and \f$  +3\pi/2 \f$ regardless, also  \f$ N_{\zeta s} = N_{\zeta p} \f$
   !!  (if \f$ N_{\zeta p} > 0 \f$) or \f$ N_{\zeta s} = M \f$
     integer(ki4) :: beq_xiopt !< .
   !> Plasma centre coordinate, units metres \f$ R_{cen} >0 \f$,
   !! default 1
     real(kr8) :: beq_rcen !< .
   !> Plasma centre coordinate, units metres \f$ Z_{cen} \f$,
   !! default 0
     real(kr8) :: beq_zcen !< .
   !> User-specified limit for flux mapping
   !!  \f$ \psi_{\min} > 0 \f$ if \f$ \psi_{opt} =1 \f$,
   !! default 0.9
     real(kr8) :: beq_psimin !< .
   !> User-specified limit for flux mapping
   !! \f$ \psi_{\max} > \psi_{\min}  \f$ if \f$ \psi_{opt} =1 \f$,
   !! default 1.1
     real(kr8) :: beq_psimax !< .
   !> Number of sample points in flux coordinate \f$ N_{\psi} > 1 \f$,
   !! default 32
     integer(ki4) :: beq_npsi !< .
   !> Fractional margin for flux mapping \f$ \Delta\psi \geq 0\f$,
   !! default 0.1
     real(kr8) :: beq_delpsi !< .
   !> User-specified limit in poloidal angle for flux mapping
   !! \f$ \theta_{\min} > 0 \f$ if \f$ \theta_{opt} =1 \f$,
   !! default \f$ \pi/4 \f$
     real(kr8) :: beq_thetamin !< .
   !> User-specified limit in poloidal angle for flux mapping
   !! \f$ \theta_{\max} > \theta_{\min} \f$ if \f$ \theta_{opt} =1 \f$,
   !! default \f$ 3 \pi/4 \f$
     real(kr8) :: beq_thetamax !< .
   !> Fractional margin for flux mapping \f$ \Delta\theta \geq 0 \f$,
   !! default 0.1
     real(kr8) :: beq_deltheta !< .
   !> User-specified limit in toroidal angle for mapping to cylindricals
   !! \f$ \zeta_{\min}>0 \f$ if \f$ \zeta_{opt} =1 \f$,
   !! default 0
     real(kr8) :: beq_zetamin !< .
   !> User-specified limit in toroidal angle for mapping to cylindricals
   !! \f$ \zeta_{\max}>\zeta_{\min} \f$ if \f$ \zeta_{opt} =1 \f$,
   !! default \f$ 2 \pi \f$
     real(kr8) :: beq_zetamax !< .
   !> Number of repeats in toroidal direction for mapping to cylindricals
   !! \f$ N_{\zeta p} \f$, default 0
     integer(ki4) :: beq_nzetap !< .
   !> Equilibrium displacement \f$ R_{mov} \f$, units mm, default 0
     real(kr8) :: beq_rmove !< .
   !> Equilibrium displacement \f$ Z_{mov} \f$, units mm, default 0
     real(kr8) :: beq_zmove !< .
   !> Equilibrium \f$ I \f$ is scaled by this quantity, default 1
     real(kr8) :: beq_fscale !< .
   !> set .TRUE. to work in duct coordinates
     logical :: beq_duct !< .
     !> flag whether any skylight(s) defined by flux limits
     !! default .FALSE.
     logical :: skylight_flux_limits !< .
     !> flag whether any skylight(s) defined by plasma centre line
     !! default .FALSE.
     logical :: skylight_centre_line !< .
   !> number of absorbers defined geometrically, default 0.
   !! Following beqparameters in the file, there must be a
   !! absorber_objects x (datvtkparameters + optional progfiles)
  integer(ki4) :: absorber_objects !< .
  !> number of invisible objects defined geometrically, default 0.
   !! Following beqparameters in the file, there must be a
   !! invisible_objects x (datvtkparameters + optional progfiles)
  integer(ki4) :: invisible_objects !< .
   !> number of skylights defined geometrically, default 0.
   !! Following beqparameters in the file, there must be a
   !! skylparameters namelist and then skylight_objects x (datvtkparameters + optional progfiles)
     integer(ki4) :: skylight_objects !< .
  !> number of beancan surfaces defined geometrically, default 0.
   !! Following beqparameters in the file, there must be a
   !! beancan_objects x (datvtkparameters + optional progfiles)
  integer(ki4) :: beancan_objects !< .
  !> number of cutouts defined geometrically, default 0.
   !! A cutout is surface designed to absorb energy without trace,
   !! e.g. to stop fieldlines going "round the back" of surfaces.
   !! Following beqparameters in the file, there must be
   !! cutout_objects x (namelists datvtkparameters + optional progfiles)
     integer(ki4) :: cutout_objects !< .
   !> Level of output to help understand skylight production
   !! Default 0, no output.
     integer(ki4) :: skylight_debug !< .
   !> \f$ M \f$ number of ripple coils, default 0
   !! for field formula see header
     integer(ki4) :: beq_mrip !< .
   !> Unused parameter for ripple coils
     real(kr8) :: beq_irip !< .
   !> \f$ a_r \f$ for ripple coils, default 0
   !! for field formula see header
     real(kr8) :: beq_arip !< .
   !> \f$ \psi_{ref}=1 \f$ (default 1). User specified value of reference flux,
   !! either flux at X-point \f$ \psi_X \f$ or at limiter \f$ \psi_{ltr} \f$
     real(kr8) :: beq_psiref !< .
   !> \f$ \theta_X \f$ User specified offset value of poloidal angle of X-point
   !! For flux mapping, defaults to zero, meaning that angles are measured relative
   !! to the vertical at bottom, i.e. usual nominal location of X-point.
   !! \f$ \theta \f$ values are specified with respect to
   !! \f$ -\pi/2 - \theta_X \f$ as origin in the poloidal system normally used
   !! where zero corresponds to the horizontal through the plasma centre.
     real(kr8) :: beq_thetaref !< .
   !> Number of sample points in poloidal angle \f$ N_{\theta} >1 \f$
   !! for flux mapping, default 32
     integer(ki4) :: beq_ntheta !< .
   !> Specifies way equilibrium field described, default 1.
   !! 1. As \f$ R/J\f$  where \f$ J \f$ is Jacobian of flux map,
   !! implies mapped case (axisymmetric), equation from
   !! Arter et al, corrected (p)reprint at http://arxiv.org/abs/1403.7142
   !! \f{eqnarray*}{
   !! d\theta/d\phi &=&  (1/I)(R/J)\\
   !! I &=& R B_T,\;\;\;\psi=\mbox{const.}
   !! \f}
   !! 2. As \f$ \psi\f$ gradients plus vacuum field, implies mapping to cylindricals (general),
   !! 3-vector equation as in Arter at al. corrected (p)reprint.
   !! 3. As \f$ R/I \times \psi\f$ gradients implies mapping to cylindricals (axisymmetric),
   !! two fieldline equations
   !! \f{eqnarray*}{
   !! dR/d\phi &=&  -(R/I) \partial \psi/\partial Z \\
   !! dZ/d\phi &=&  (R/I) \partial \psi/\partial R \\
   !! I &=& R B_T
   !! \f}
     integer(ki4) :: beq_fldspec !< .
     logical :: skip_eqdsk_b !< read B cpts from eqdsk file unless .TRUE., default .FALSE.
   !> Specifies how to search for the X-point, hence \f$ \psi_X \f$, default 0 \n
   !! 0. Find the mesh-points with the smallest value of flux gradient
   !! above and below plasma
   !! centre and look near these two points. \n
   !! 1. Search within box specified by beq_xrsta, beq_xrend, beq_xzsta
   !! and beq_xzend. The box must be at least one mm large in each dimension, but
   !! the search will probably be unsuccessful unless the box is more than one
   !! mesh cell in size in each direction.
     integer(ki4) :: beq_xsearch !< .
     real(kr8) :: beq_xrsta !< X-point search box mininimum \f$ R \f$ (m), default 0
     real(kr8) :: beq_xrend !< X-point search box maximum \f$ R \f$ (m), default 0
     real(kr8) :: beq_xzsta !< X-point search box mininimum \f$ Z \f$ (m), default 0
     real(kr8) :: beq_xzend !< X-point search box maximum \f$ Z \f$ (m), default 0
   !> Specifies how to search for the contact point, hence \f$ \psi_{ltr} \f$, default 0 \n
   !! 0. Find the mesh-point with the value of flux  closest to that
   !! of the centre of the plasma. \n
   !! 1. Look only points within box specified by search_r_start, search_r_end, search_z_start
   !! and search_z_end. The box must be at least one mm large in each dimension, but
   !! the search will probably be unsuccessful unless the box is more than one
   !! mesh cell in size in each direction.
     integer(ki4) :: limiter_search !< .
     real(kr8) :: search_r_start !< search box mininimum \f$ R \f$ (m), default 0
     real(kr8) :: search_r_end !< search box maximum \f$ R \f$ (m), default 0
     real(kr8) :: search_z_start !< search box mininimum \f$ Z \f$ (m), default 0
     real(kr8) :: search_z_end !< search box maximum \f$ Z \f$ (m), default 0
   !> default 0, value of 1 implies \f$ \psi \f$ is \f$ 2 \pi \f$ too big
   !! relative to usual EQDSK convention.
     real(kr8) :: beq_psibig !< .
   !> equilibrium helicity is ok, no need for overrides which are ignored if .TRUE.
     logical :: equil_helicity_ok !< .
   logical :: equil_scale !<  default false, if true, scale flux by psic-psib (equ case)
   logical :: equil_mastequ !<  default true, MAST-U equ file has right helicity even if equil_helicity_ok=.false.
   !> Specifies file containing spl3d format vacuum field
     character(len=80) :: beq_vacuum_field_file !<  .
  end type beqparameters

!---------------------------------------------------------------------
!> Entries in namelist btreeparameters.
!! These define the fundamentals of the
!! hierarchical data structure (HDS) used in the
!! fieldline tracing, discussed in the ray-tracing context by the
!! preprint at http://arxiv.org/abs/1403.6750
!! The variables are of two types:
!! 1. Allocates an amount of static storage suitable for a geometry
!! with up to 1 million triangles approximately.
!! 2. Determine the type of HDS to be used, normally multi-octree.
!! - Vector components \f$ (x,y,z) \f$ correspond to those of
!! mapped coordinate system, respectively either
!! \f$ (\psi,\theta,\zeta) \f$ or \f$ (R,Z,\xi) \f$.
  type, public :: btreeparameters
   !> EXPERT  size of binary tree, must be positive,
   !! default 2000000
     integer(ki4):: btree_size !< .
   !> EXPERT  size of binary tree pter, must be positive,
   !! default 3
     integer(ki4):: btree_sizep !< .
   !> EXPERT  size of exten array, must be positive,
   !! default 1000
     integer(ki4):: btree_sizee  !< local variable
   !> EXPERT  size of list array hoc, must be positive,
   !! default 2
     integer(ki4):: btree_sizeh !< .
   !> EXPERT  size of list array, must be positive,
   !! default 2000000
     integer(ki4):: btree_sizel !< .
   !> EXPERT  max depth of tree, must be positive,
   !! default 30
     integer(ki4):: btree_depth !< .
   !> \f$ (n_x,n_y,n_z) \f$. Vector of no. of children at
   !! top of multi-octree, components must be positive,
   !! octree at all levels is default, vector (2,2,2),
   !! change only if ttalg=2
     integer(ki4), dimension(3) :: tree_nxyz !< .
   !> EXPERT  \f$ {\bf h}=(h_x,h_y,h_z) \f$ Use this vector of spacings to set
   !! octree quantisation, components must be positive
   !! default spacings (0.1,0.1,0.1)
     real(kr4), dimension(3) :: tree_hxyz !< .
   !> define quantisation algorithm if multi-octree
   !! 0. automatic
   !! 1. default, use hxyz (spacings) to set top
   !! 2. use nxyz, set split at top directly
     integer(ki4):: tree_ttalg !< .
   !> type of tree, must be positive
   !! 1. BSP  (EXPERT)
   !! 2. octree
   !! 3. multi-octree (octree with special top)
     integer(ki4):: tree_type !< .
  end type btreeparameters

!!---------------------------------------------------------------------
!> Entries in namelist equilparameters.
!! Normally equilibria are specified via an EQDSKG  file output by EFIT,
!! but there is also the possibility for the user to define equilibria
!! if the equilibrium file extension is given as 'ana'.
!! - For a description of EQDSKG format
!! see https://fusion.gat.com/theory/Efitgeqdsk
!! - File xxx.ana must contain namelists equilparameters and
!! meshparameters in that order
!! 1.  An analytic equilibrium may be specified following Solovev.
!! The Solovev solution of the Grad-Shafranov equation is given by
!! \f[
!! \psi  =  \psi_0  \chi  ,  0 \le \chi \le 1  ,  \chi =  \left [ (R^2
!! - R^2_0)^2  +  \frac{Z^2} {E^2}  (R^2 - R^2_x) \right ] \frac{4}
!! {(R^2_2 - R^2_1)^2}
!! \f]
!! ultimately determined by a set of four separate
!! positions \f$(R_1,R_2,R_m,Z_m)\f$ which define a D-shaped reference
!! surface going through the three points \f$(R_1,0), (R_2,0) \f$ and \f$ (R_m,Z_m)\f$.
!! (The solution is up-down symmetric: \f$\psi(R,Z) = \psi(R,-Z)\f$.)
!! These four positions are used to define the parameters via
!! \f{eqnarray*}{
!! R_0^2 &=& \frac {R_1^2 + R_2^2} {2},
!! R_x^2 &=& \frac{R_1^2 R_2^2 - R_{\rm m}^4} {R_1^2 + R_2^2 - 2R_ m^2,},
!!  E^2 &=& \frac {Z_{\rm m}^2} {R_1^2 + R_2^2 - 2 R_m^2}
!! \f}
!!      - plasma minor radius \f$a = \frac{R_2 - R_1} {2} \f$
!!      - it is up to the user to choose parameters sensibly.
!! Some indicative values are set in the source code beqan_m.f90
!! 2. User-defined  as special f95 subroutine coded by user, called
!! beqan_userdefined in file beqan_m.f90. The parameter arrays below are
!! provided for convenience.
!! 3. Remember that equilibria may have left or right-handed helicity, and
!! this may be changed by the code unless equil_helicity_ok=.TRUE.
  type, public :: equilparameters
   !> equilibrium formula, currently recognised values are
   !! - 'solovev'
   !! - 'userdefined'
     character(len=80) :: equilibrium_formula !< .
     real(kr8) :: BPhi !< \f$ B_T \f$ Value of toroidal field at plasma centre
     real(kr8) :: Psi0 !< \f$ \psi_0 \f$ Flux through plasma
     real(kr8) :: R1 !< \f$ R_1>0 \f$ inner cylindrical radius
     real(kr8) :: R2 !< \f$ R_2>R_1 \f$ outer cylindrical radius
     real(kr8) :: Rm !< \f$ R \f$ at maximum height, \f$ R_1<R_m<R_2 \f$
     real(kr8) :: Zm !< \f$ Z \f$ at maximum height
     real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters !< general real parameters
     integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters !< general integer_parameters
     integer(ki4) :: number_of_real_parameters !< actual  number of real parameters
     integer(ki4) :: number_of_integer_parameters !< actual  number of integer parameters
  end type equilparameters

!!---------------------------------------------------------------------
!> Entries in namelist meshparameters.
!! These define the mesh as a uniform lattice in R-Z
  type, public :: meshparameters
     real(kr8) :: rmin !< \f$ R_{\min} \f$ Minimum value of \f$ R \f$, \f$ R_{\max} \neq R_{\min} \f$
     real(kr8) :: rmax !< \f$ R_{\max} \f$ Maximum value of \f$ R \f$, \f$ R_{\max} \neq R_{\min} \f$
     integer(ki4) :: nw !< \f$ N_R>0 \f$ number of grid points in R-direction
     real(kr8) :: zmin !< \f$ Z_{\min} \f$ Minimum value of \f$ Z \f$, \f$ Z_{\max} \neq Z_{\min} \f$
     real(kr8) :: zmax !< \f$ Z_{\max} \f$ Maximum value of \f$ Z \f$, \f$ Z_{\max} \neq Z_{\min} \f$
     integer(ki4) :: nh !< \f$ N_Z>0 \f$ number of grid points in Z-direction
  end type meshparameters

!!---------------------------------------------------------------------
!> Entries in namelist hdsgenparameters.
!! These control subdivision that generates HDS.
!! A box corresponding to a leaf of the binary tree is referred
!! to as a 'bin'.
  type, public :: hdsgenparameters
   !> Stop dividing if less than this number of objects
   !! in 'bin', default 10.
     integer(ki4) :: limit_geobj_in_bin !< .
   !> INERT numerical parameter
     integer(ki4) :: nbins_child  ns_child!< .
   !> Type of geometry
   !! 1. Points (default)
   !! 2. Triangles
     integer(ki4) :: geometrical_type !< .
   !> INERT Positive tolerance <= max_tolerance
     real(kr4) :: min_tolerance !< .
   !> Positive tolerance, default 0.1.
   !! Expand bounding box by this margin in every direction, units of input geometry
     real(kr4) :: max_tolerance !< .
   !> Number of extra bounding cubes, no more than two, default 0
   !! If delta_inner_length>=0 then add this margin to volume,
   !! and if delta_outer_length>=0 add inner+outer_length margin.
   !! If delta_inner_length<0, then use corner values to define and only
   !! produce one cube. N.B. this will be expanded to include all objects
     integer(ki4) :: no_boundary_cubes !< .
   !> vector defining corner with lesser component values
   !! default origin (absolute coordinates)
   !! Only active if any of the components is less than the smallest
   !! corresponding coordinate in the geometry, and then only that coordinates are affected.
     real(kr4), dimension(3) :: lower_corner !< .
   !> vector defining corner with greater component values (absolute coordinates)
   !! Only active if any of the components is greater than the largest
   !! corresponding coordinate in the geometry, and then only those coordinates are affected.
     real(kr4), dimension(3) :: upper_corner !< .
   !> Inner bounding cube separation from geometry, mm,
   !! default -1
     real(kr4) :: delta_inner_length   !< .
   !> Outer bounding cube separation from inner, mm,
   !! default 0
     real(kr4) :: delta_outer_length   !< .
   !> \f$ N_Q \f$, default 1024
     integer(ki4) :: quantising_number !< .
   !> EXPERT Number of objects in tree, default 0, normally determined automatically
     integer(ki4) :: no_geobj_records !< .
   !> If positive, expand global bounding box bin by
   !! margin determined by other settings, default 0
     integer(ki4) :: margin_type !< .
   !> INERT type of scaling of geobj position values
     integer(ki4) :: type_geobj_coord_scaling !< .
  end type hdsgenparameters

!---------------------------------------------------------------------
!> Entries in namelist edgprofparameters.
!! These define the profile of edge power deposition.
!! -The standard and Eich formulae are defined in Arter et al, 2015,
!! preprint at http://arxiv.org/abs/1403.7142
!! using the same notation.
!! - The expdouble formula is taken from Kocan et al Eq(1),
!! paper at http://iopscience.iop.org/0029-5515/55/3/033019 as
!! \f$ q_{\parallel 0}\left(\exp(-\Delta r/\lambda_q) +
!! R_{near} \exp(-\Delta r/\lambda_{q,near}) \right)\f$
!! where \f$ \Delta r \f$ is the midplane separation between the flux
!! contour on the PFC and the LCFS (positive outside LCFS).
!! - One of \f$ P_{loss} \f$  and \f$ q_{\parallel 0} \f$ must be positive,
!! one which is negative is ignored
!! - Maximum number of positions at which to sample profile,
!! MAX_NUMBER_OF_POSITIONS=100
!! - Maximum number of parameters which may used to define fuction,
!! MAX_NUMBER_OF_PARAMETERS=100
  type, public :: edgprofparameters
   !> power deposition formula to be used.
   !! - 'unset' standard exponential profile
   !! - 'exp' standard exponential profile
   !! - 'expdouble' as two exponentials
   !! - 'eich' as defined by T. Eich
   !! - 'samples' as a set set of samples at specified positions.
   !! N.B. it is up to the user to arrange normalisation.
   !! - 'userdefined' as special f95 subroutine coded by user, called
   !! edgprof_userdefined in file edgprof_m.f90.
   !! N.B. it is up to the user to arrange normalisation.
     character(len=80) :: profile_formula !<.
   !> \f$ F \f$ power split ion to electron direction,
   !! must lie between 0 and 1, default 0.5
     real(kr8) :: power_split !<.
   !> \f$ \lambda_m \f$ power decay length at midplane (m),
   !! \f$ \lambda_m \geq 0 \f$, also doubles as \f$ \lambda_q \f$
   !! if Eich formula used, default 0.012 (12mm)
     real(kr8) :: decay_length !<.
   !> \f$ P_{loss} \f$ power crossing the separatrix (W),
   !! default 7.5e+06 (7.5MW)
     real(kr8) :: power_loss !<.
   !> \f$ q_{\parallel 0}  \f$ parallel reference flux (W/sqm).
   !! Default -1, if negative, ignored
     real(kr8) :: q_parallel0 !<.
   !>  \f$ \sigma \f$ Eich formula length (m)
     real(kr8) :: diffusion_length !<.
   !> \f$ R_{near} \f$, near-SOL parallel relative reference flux (dimensionless)
   !! for expdouble formula
     real(kr8) :: ratio_of_q_parallel0_near !<.
   !> \f$ \lambda_{q,near} \f$ near-SOL power decay length at midplane (m)
   !! for expdouble formula
     real(kr8) :: decay_length_near !<.
   !> Define profile as samples at these positions
     real(kr8), dimension(MAX_NUMBER_OF_POSITIONS) :: positions  !<.
   !> Profile determined as samples
     real(kr8), dimension(MAX_NUMBER_OF_POSITIONS) :: deposition_profile  !<.
   !> whether positions in distance or flux from SOL
   !! - 'radius'
   !! - 'flux'
     character(len=80) :: coord_of_positions !<.
   !> Number of positions actually used. If sampling,
   !! must be positive, not larger than MAX_NUMBER_OF_POSITIONS
     integer(ki4) :: number_of_positions  !<.
   !> Parameters which may be used in user-defined function
     real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters  !<.
   !> Parameters which may be used in user-defined function
     integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters  !<.
   !> Number of real parameters actually used. If defining function,
   !! must be positive, not larger than MAX_NUMBER_OF_PARAMETERS
     integer(ki4) :: number_of_real_parameters  !<.
   !> Number of integer parameters actually used. If defining function,
   !! must be positive, not larger than MAX_NUMBER_OF_PARAMETERS
     integer(ki4) :: number_of_integer_parameters  !<.
  end type edgprofparameters

!---------------------------------------------------------------------
!> Entries in namelist odesparameters.
!! The system integrated depends on the parameter fldspec, in fact has order
!! dependent on fldspec. (In axisymmetric flux coordinates, order is unity,
!! two in case of cylindrical axisymmetric, three in general cylindricals.)
!! - Minimum values for \f$ \epsilon_r \f$, \f$ \epsilon_a \f$ and
!! \f$ \Delta t_0 \f$ should in principle be much smaller in the
!! cylindrical mapped cases
!! than in flux-mapped, but in practice it makes little difference to the cost
!! of a flux-mapped run if smaller values than strictly necessary are used.
!! Values below single precision round-off are not however sensible.
!! - Coordinate directions \f$ 1,2,3 \f$ correspond to respectively
!! \f$ (\psi,\theta,\zeta) \f$ in case of flux coordinates.
!! - Coordinate directions \f$ 1,2,3 \f$ correspond to respectively
!! \f$ (R,Z,\xi) \f$ (N.B. \f$ \xi \f$ not \f$ \zeta \f$) in cylindrical case.
  type, public :: odesparameters
   !> EXPERT system order,
   !! 1. for  flux-mapped case
   !! 2. for axisymmetric field
   !! 3. for field with ripple
     integer(ki4) :: system_order !<.
   !> INERT start control
     integer(ki4) :: start_control !<.
   !> INERT termination control
   !! - Positive value corresponds to fieldline rising above plane,
   !! negative value corresponds to fieldline falling below plane.
   !! - Position of plane determined by termination_parameters(3).
     integer(ki4) :: termination_control !<.
   !> termination parameters
   !! - (1) Distance to avoid self-intersection, in quantised units,
   !!  typically 1/1024 of geometry, default 1.
   !! - (2) Maximum allowed length of one fieldline segment,
   !! units are domain size, default effectively \f$ \infty \f$.
   !! - (3) INERT Terminating plane position in Z, units metres,
   !! superseded by more flexible definitions in namelist termplaneparameters.
   !! (EXPERT Only active if termination_planes=F and parameter lpmidplane=T
   !! with fldspec=2, default effectively \f$ \infty \f$.)
     real(kr8), dimension(3)  :: termination_parameters !<.
   !> relative error \f$ \epsilon_r \f$,
   !! must be >=0 and <=1, default \f$ 10^{-4} \f$
     real(kr8) :: rel_error !<.
   !> absolute error \f$ \epsilon_a \f$,
   !! must be >=0 and <=1, default \f$ 10^{-4} \f$
     real(kr8) :: abs_error !<.
   !> Maximum number of steps for one fieldline integration,
   !! must be positive, default 100.
     integer(ki4) :: max_numsteps   !<.
   !> Initial step-size \f$ \Delta t_0 =0.1\f$ (default) for fieldline integration,
   !! i.e. increment in independent variable, must be positive.
   !! (Code will determine which is outgoing direction.)
     real(kr8) :: initial_dt !<.
   !> Maximum allowed change in \f$ \zeta \f$ (radians) over
   !! one fieldline integration, default 1.
     real(kr8) :: max_zeta !<.
  end type odesparameters

!---------------------------------------------------------------------
!> Entries in namelist powcalparameters.
!! This is the main parameter control for powcal, however there are also
!! subsidiary controls, which follow in order:
!! - termplaneparameters (must only  be present if termination_planes=.TRUE.)
!! - edgprofparameters (must only  be present if more_profiles=.TRUE.,
!! now controls power deposition formula)
!! - odesparameters
!!
!! Maximum number of tracks allowed for explicit selection, maximum_number_of_tracks=100.
  type, public :: powcalparameters
   !> analytic launch type, for EXPERT use only
     integer(ki4) :: analytic_launch_type !< local variable
   !> calculation type
   !! - 'local' shadowing only by  toroidally localised geometry
   !! - 'global' shadowing by fieldlines which may wrap around torus
   !! - 'middle' EXPERT start fieldlines at "midplane", Z=const.
     character(len=80) :: calculation_type !< local variable
   !> INERT if more_profiles, \f$ \lambda_m \f$ power decay length at midplane (m),
   !! default 0.012 (12mm)
     real(kr8) :: decay_length !< local variable
   !> INERT if more_profiles, \f$ \sigma \f$ Eich formula length (m),
   !! default 0
     real(kr8) :: diffusion_length !< local variable
   !> tracks selected by element number
     integer(ki4), dimension(maximum_number_of_tracks) :: element_numbers !< local variable
   !> If nonzero, compute only selected tracks, default 0.
   !! Tracks selected using element number. Remember to set plot_flinx/m .TRUE. to get
   !! track output , and setting plot_powx .FALSE. will save disk space
     integer(ki4) :: number_of_tracks !< local variable
   !> INERT if more_profiles, \f$ P_{loss} \f$ power crossing the separatrix (W) ,
   !! default 7.5e+06 (7.5MW)
     real(kr8) :: power_loss !< local variable
   !> power on launch geometry, for EXPERT use only
     integer(ki4) :: power_on_launch_geo !< local variable
   !> power on shadow geometry, for EXPERT use only
     integer(ki4) :: power_on_shadow_geo !< local variable
   !> INERT if more_profiles, \f$ F \f$ power split ion to electron direction,
   !! must lie between 0 and 1, default 0.5.
     real(kr8) :: power_split !< local variable
   !> INERT if more_profiles, \f$ q_{\parallel 0} \f$ parallel reference flux (W/sqm).
   !! Default is -1 and if negative is ignored.
     real(kr8) :: q_parallel0 !< local variable
   !> level of refinement to be used, always unity,
   !! EXPERT increasing to 2 or 3 increases number of samples by 4X or 16X
   !! if mesh has been appropriately prepared using CFMCNP
     integer(ki4) :: refine_level !< local variable
   !> control shadowing
   !! - (-2). "Paint power" onto surfaces
   !! - (-1). For testing purposes only
   !! - 0. No shadowing (default)
   !! - 1. Shadowing by shadow geometry
     integer(ki4) :: shadow_control !< local variable
   !> Default .FALSE. implies use of value psista (sic) for \f$ \psi \f$ on termination plane
   !! to evaluate \f$ Q \f$ in global calculation.  If .TRUE., use instead value of \f$ \psi \f$ on
   !! PFC to evaluate \f$ Q \f$. (Integration errors may lead to discrepancies in psista even in
   !! the case of exisymmetric field.)
     logical :: use_pfc_flux   !< .
   !> termination planes present, always set .TRUE. to use more
   !! sophisticated termination criteria in termplaneparameters namelist,
   !! default .FALSE.
     logical :: termination_planes !< local variable
   !> more profiles can be used, always set .TRUE. to use more
   !! power deposition profiles in edgprofparameters namelist,
   !! default .FALSE.
     logical :: more_profiles !< local variable
     logical :: skylight !< activate skylights if present, default .TRUE.
     logical :: power_in_pfr !<  default .TRUE., power in private flux region
     logical :: employ_skyl_pfr !<  INERT default .TRUE., employ skylights in private flux region
   !> Distance that assumes intersection with skylight, in quantised units, in each direction
   !!  typically 1/1024 of geometry, default 1.
     real(kr8), dimension(3) :: skylight_tolerance !< .
  end type powcalparameters

!---------------------------------------------------------------------
!> Entries in namelist plotselections.
!! There is a namelist plotselections for each ot the major component
!! programs, datvtk, geoq, hdsgen, smanal and powcal (although required for datvtk
!! if a .ctl file is used, the namelist controls have no effect). The plots that
!! can be produced vary with each program, but the following conventions
!! are observed:
!! 1. All variables are logicals, and must be set .TRUE. to produce a
!! data file for plotting, or be omitted, or set .FALSE. (otherwise an error will
!! terminate execution).
!! 2. All variable names begin with 'plot_', then followed by an identifier
!! for the particular program if a vtk file is wanted, thus
!!     - geoq or geofld for geoq
!!     - hds or geo for hdsgen
!!     - pow for powcal
!!     - if gnu appears, a file for gnuplot is produced
!! 3. All vtk files (for visualisation using ParaView) contain only surface
!! information, unless the string
!! 'vol' appears following the program identifier, or they show the HDS
!! which is inherently a 3-D wireframe.
!! 4. All gnu files (for visualisation using gnuplot) contain only
!! information projected in cross-section on a plane \f$ \zeta=const. \f$.
!! The main gnu files of interest are produced by geoq and may be
!! processed using the command geoqgnu root, where root.ctl is the
!! control file used to generate them
!! 5. The coordinate system used is indicated by the last letter before
!! the .vtk filetype, viz:
!!     - x for Cartesians (mm)
!!     - m for mapped coordinates (either flux coordinates or cylindricals
!! depending on fldspec, the way the magnetic field is specified). Angles
!! are in radians, and in the
!! case of cylindricals, units of metres are used for R and Z.
!!     - q for quantised coordinates, the system in which lines are actually
!! tracked, typically ranging over 0 to 1024.
!! 6. There are also names which are a legacy of previous developments,
!! if they are set .TRUE., a warning will be issued.
!! 7. Corresponding output files to variable plot_xxx have names root_xxx.vtk
!! or root_xxx.gnu or root_gnu(_opt).gnu, where root.ctl is the
!! control file used to generate them
!! 8. If it is not sensible to produce a plot file, its output will be
!! suppressed
  type, public :: plotselections
   !> geoq only, above conventions apply
     logical :: plot_geoqx !< .
   !> geoq only, above conventions apply
     logical :: plot_geoqvolx !< .
   !> geoq only, above conventions apply, except 'mm' implies field is mapped
   !! as well as coordinates
     logical :: plot_geoqmm !< .
   !> geoq only, above conventions apply
     logical :: plot_geoqm !< .
   !> geoq only, above conventions apply
     logical :: plot_geoqfldx !< .
   !> geoq only, above conventions apply. Files with _opt component will
   !! be produced, depending on fldspec, with content determined by opt
   !! - (omitted) \f$ \psi \f$ values in R-Z coordinates
   !! - dR \f$ \partial \psi/\partial R \f$ values in R-Z coordinates
   !! - dZ \f$ \partial \psi/\partial Z \f$ values in R-Z coordinates
     logical :: plot_gnu !< .
   !> geoq only, above conventions apply. Files with _opt component will
   !! be produced, depending on fldspec, with content determined by opt
   !! - R \f$ R \f$ values in \f$  \psi-\theta \f$  space
   !! - Z \f$ Z \f$ values in \f$  \psi-\theta \f$  space
   !! - RJ \f$ R/J \f$ values in \f$  \psi-\theta \f$  space
   !! - RZ \f$ R \f$ and \f$ Z \f$ values in \f$  \psi-\theta \f$  space
     logical :: plot_gnum !< .
   !> geoq only, above conventions apply. Points of the geometry
   !! projected in R-Z coordinates.
     logical :: plot_gnusil !< .
   !> geoq only, above conventions apply. Points of the geometry
   !! projected in \f$  \psi-\theta \f$  space
     logical :: plot_gnusilm !< .
   !> hdsgen only, above conventions apply.
     logical :: plot_hdsm !< .
   !> hdsgen only, above conventions apply.
     logical :: plot_hdsq !< .
   !> hdsgen only, plot geometry in quantised coordinates
     logical :: plot_geobjq !< .
   !> hdsgen only, plot points of geometry in quantised coordinates
     logical :: plot_geoptq !< .
     logical :: plot_cartv !<  magtfm only, vtk plot selector, default .FALSE.
     logical :: plot_allcartv !<  magtfm only, vtk plot selector, default .FALSE.
     logical :: plot_testtfm !<  magtfm only, do test of FFT analysis, default .FALSE.
     logical :: plot_bcartx !<  magtfm only, output scaled B in fmesh format in txt file, default .FALSE.
     logical :: plot_gnuv !<  magtfm only, gnu 2-D plot selector, default .FALSE.
     logical :: plot_modes !<  magtfm only, gnu 2-D plot selector, default .FALSE.
     logical :: plot_gnu !<  magtfm only, gnu 1-D plot selector, default .FALSE.
     logical :: plot_smanalout !< smanal output data selector, default  .FALSE.
     logical :: plot_anx !< smanal only, vtk plot selector, default  .FALSE.
     logical :: plot_ansmallx !< smanal only, vtk plot at centroids only, default  .FALSE.
     logical :: plot_angnu !< smanal only, gnuplot plot selector, default  .FALSE.
   !> powcal only, EXPERT use only, above conventions apply
     logical :: plot_powstatx !< .
   !> powcal only, above conventions apply.
     logical :: plot_powx !< .
   !> powcal only. Plot power on walls if 'middle' option set.
     logical :: plot_wall !< .
   !> powcal only, above conventions apply.  Files for fieldlines
   !! are produced, names beginning 'trackx', numbered by start element,
   !! with '01.vtk' terminating
     logical :: plot_flinx !< .
   !> powcal only. First and last points of
   !! each fieldline are output in gnuplot format in Cartesian coordinates
     logical :: plot_flinends !< .
   !> powcal only, above conventions apply.  Files for fieldlines
   !! are produced, names beginning 'trackm', numbered by start element,
   !! with '01.vtk' terminating
     logical :: plot_flinm !< .
   !> geoq only, produce files eqdsk_eqbdry.gnu and eqdsk_eqltr.gnu
   !! describing boundary and limiter from eqdsk file with root 'eqdsk'
     logical :: plot_eqdsk_boundary !< .
     logical :: plot_skylprovis_lower !< geoq only, provisional lower skylight plot selector
     logical :: plot_skylprovis_upper !< geoq only, provisional upper skylight plot selector
  end type plotselections

!---------------------------------------------------------------------
!> Entries in namelist termplaneparameters.
!! Terminate if fieldline intersects or hits a plane going in a certain direction,
!! possibly subject to an auxilliary condition on position.
!! Only activated if termination_planes = .TRUE. in namelist powcalparameters.
!! - Coordinate directions \f$ 1,2,3 \f$ correspond to respectively
!! \f$ (R,Z,\xi) \f$, N.B. \f$ \xi \f$ not \f$ \zeta \f$.
!! - Note arrangement that array entries align with termplane_position,
!! if setting termplane_direction, ensure corresponding termplane_intersection
!! is zero and vice versa.
!! - Maximum number of inputs allowed for each variable, maxinp=100.
  type, public :: termplaneparameters
   !> Terminate if fieldline hits a plane going in a certain direction,
   !! possibly subject to an auxilliary condition on position,
   !! terminating directions \f$ \pm 1,\pm 2, \pm 3 \f$ allowed,
   !! or dummy entry of zero.
   !! - Positive direction corresponds to fieldline rising above plane,
   !! negative direction corresponds to fieldline falling below plane.
   !! - Position of plane determined by termplane_position.
     integer(ki4), dimension(maxinp):: termplane_direction  !< .
   !> Terminate if fieldline intersects a plane,
   !! plane directions \f$ 1-3 \f$ allowed, or dummy entry of zero.
   !! - Position of plane determined by termplane_position.
     integer(ki4), dimension(maxinp):: termplane_intersection  !< .
   !> Terminate if fieldline reaches an extremum,
   !! terminating directions \f$ \pm 1,\pm 2, \pm 3 \f$ allowed.
   !! - Positive sign implies only maximum recognised,
   !! negative sign implies only minimum recognised,
     integer(ki4), dimension(maxinp):: termplane_extremum  !< .
   !> Position of plane in coordinate \f$ 1-3 \f$ as specified,
   !! units metres or radians
     real(kr8), dimension(maxinp):: termplane_position  !< .
   !> Terminate when fieldline intersects a plane provided
   !! an auxilliary condition on position in this direction has been met,
   !! terminating directions \f$ \pm 1,\pm 2, \pm 3 \f$ allowed,
   !! - Positive value corresponds to greater than criterion
   !! negative value corresponds to less than
   !! - Position of plane determined by termplane_condition.
     integer(ki4), dimension(maxinp):: termplane_condition_dir  !< .
   !> Position of auxilliary condition plane in coordinate \f$ 1-3 \f$ as specified,
   !! units metres or radians
     real(kr8), dimension(maxinp):: termplane_condition  !< .
   !> INERT number of active termination planes. This is determined
   !! automatically thanks to the arrangement with dummy zeroes described above.
     integer(ki4):: termplane_max  !< .
  end type termplaneparameters

!---------------------------------------------------------------------
!> Entries in namelist panelarrayparameters.
!! This assigns transforms to bodies, strictly objects defined by
!! the numbers defined by process_by_name in namelist
!! vtktfmparameters (if used). It is easiest to understand
!! if every body appears in the
!! panel_bodies list and has a corresponding transform,
!! unless the 'split' option is set when each transform corresponds
!! to two adjacent bodies in the list. However, the software will by
!! default enter all bodies into a list of length up to
!! maximum_number_of_panels and assign the
!! identity transformation to each.
!!
  type, public :: panelarrayparameters
   !> bodies defining the geometry,
   !! if no bodies are listed, then all are considered
     integer(ki4), dimension(2*maximum_number_of_panels) :: panel_bodies !< .
   !> number of transform to apply,
   !! default -1 gives error so must be changed
     integer(ki4), dimension(2*maximum_number_of_panels) :: panel_transform !< .
  end type panelarrayparameters

!---------------------------------------------------------------------
!> Entries in namelist vtkfiles.
!! This sets the names of the vtk format files to be used by vtktfm,
!! and how many copies in total are needed of each file.
!! Note that if a file is to be copied, it is assumed that it contains
!! only one body (use datvtk -xxm to remove body identifiers from dat files
!! when converting to vtk).
!! The copies of the first file in the vtk_input_file list are numbered
!! 101, 102, 103, ... etc, the second are numbered 201, 202, 203, ...etc.
!! (so that if 4 occurrences of the geometry in file 2 are needed,
!! set the number of copies to 4, and reference them as 201-204).
!!
  type, public :: vtkfiles
   !> vtk format geometry data input files,
   !! default 'null' gives error so must be changed
     character(len=80), dimension(maximum_number_of_files) :: vtk_input_file !< .
   !> number of copies needed of each corresponding input file,
   !! default 1 for each entry in vtk_input_file list
     integer(ki4), dimension(maximum_number_of_files) :: number_of_copies !< .
   !> labels vtk format geometry data input file (not used),
   !! default 'null'
     character(len=80), dimension(maximum_number_of_files) :: vtk_label  !< .
   !> numeric label for input file (not used),
   !! default 0
     integer(ki4), dimension(maximum_number_of_files) :: vtk_label_number !< .
   !> vtk format geometry data output file,
   !! defaults to root of ctl file (plus vtk suffix)
     character(len=80) :: vtk_output_file  !< .
  end type vtkfiles

!---------------------------------------------------------------------
!> Entries in namelist positionparameters.
!! These describe a transform of position in a 3-D coordinate system.
!!
  type, public :: positionparameters
   !> transformation 3x3 matrix,
   !! default is identity matrix
     real(kr4), dimension(3,3):: position_matrix  !< .
   !> transformation scale 3-vector,
   !! default \f$ (1.,  1.,  1. ) \f$
     real(kr4), dimension(3):: position_scale  !< .
   !> transformation offset 3-vector,
   !! default \f$ (0.,  0.,  0. ) \f$
     real(kr4), dimension(3):: position_offset  !< .
   !> transformation type,
   !! default 0 (identity),
   !! - 1. scale - scale in Cartesian coordinates
   !! - 2. offsca - scale and offset  in Cartesians
   !! - 6. poltil - poloidal tilt, rotate about an axis parallel to the
   !! toroidal direction through the centre of the body
   !! - 7. tortil - toroidal tilt, rotate about an axis in the
   !! vertical direction through the centre of the body
   !! - 12. torrot - rotate in toroidal direction
   !! - 22. polrot - rotate in poloidal direction
   !! - 42. radial - displace towards a point \f$ (R,Z) \f$, which will
   !! normally be the centre of the torus (i.e. displace in minor radius),
   !! negative displacements imply outwards translation. The first entry is
   !! the displacement, and the next two give the point coordinates
     integer(ki4):: position_transform  !< .
   !> describes transformation type,
   !! default 'cartesian_scale' (inert)
     character(len=20) :: transform_desc !< .
   !> identifies transformation,
   !! default 'unit' (inert)
     character(len=20) :: transform_id !< .
  end type positionparameters

!---------------------------------------------------------------------
!> Entries in namelist magfiles.
!! These define the input and output files for the vacuum magnetic
!! field data transformed by magtfm.
!!
  type, public :: magfiles
   !> txt format magnetic field data input file,
   !! default 'null' implies that there are no 3-D (non-axisymmetric)
   !! field effects to be considered.
   !! The precise format is defined in \ref magtfm.
     character(len=80) :: mag_input_file  !< .
   !> sp3 format magnetic field data output file,
   !! for input to powcal (internal to SMARDDA-PFC)
   !! setting to 'null' suppresses calculation as well as output
   !! of sp3 file, in which case it would be normal to have
   !! set plot_bcartx=T to produce txt format file
   !! for input to powcal
     character(len=80) :: mag_output_file  !< .
  end type magfiles

!---------------------------------------------------------------------
!> Entries in namelist progfiles.
!! These set the file(s) used by datvtk -c.
!! Either r_input_file and z_input_file  should contain
!! equal numbers of \f$ R \f$ and \f$ Z \f$ coordinates respectively, one per line, with
!! a null value for rz_input_file, or rz_input_file  should
!! contain the \f$ (R,Z) \f$ positions, one pair of numbers per line, and the other two
!! filenames should be 'null'.
!!
!! The file format used here is referred to as the 'gnu dat' format, to
!! distinguish it from the NASTRAN .dat format, as it is compatible with the
!! gnuplot plotting package. It allows for a single header line beginning '#'
!! with a closely defined layout to describe the content.
!!
  type, public :: progfiles
   !> position coordinate 1 \f$ (R) \f$ data input file,
   !! default 'null'
     character(len=80) :: r_input_file  !< .
   !> position coordinate 2 \f$ (Z) \f$ data input file,
   !! default 'null'
     character(len=80) :: z_input_file !< .
   !> position coordinate 1 and 2 \f$ (R,Z) \f$ data input file,
   !! default 'null'
     character(len=80) :: rz_input_file  !< .
  end type progfiles

!---------------------------------------------------------------------
!> Entries in namelist datvtkparameters.
!! This is the main parameter control for datvtk, when
!! a ctl file (datvtk -c) is used to generate a surface from file(s)
!! containing \f$ R \f$ and/or \f$ Z \f$ point values.
!! (Note that otherwise datvtk
!! acts as a filter, converting file.dat to file_out.vtk.)
!! This namelist also provides data to geoq to construct surfaces,
!! see Notes.
!!
!! Remember that angle \f$ \zeta = -\theta \f$ where \f$ \theta \f$ is
!! the more usual cylindrical polar angle. Despite the additional
!! restriction that \f$ |\zeta| < 360 \f$ (in degrees), it is possible
!! to cover any range of angles by suitable choice of limit angles.
!! As an example, to produce a cylinder with a gap of \f$ 15^o \f$ near
!! the negative \f$ y \f$ axis in \f$ x>0 \f$, set
!! start_angle=-280., finish_angle=65..
!! Order of start and finish matters in that the surface normal is
!! outward from the origin if the finish angle is larger than
!! the start angle, inward otherwise.
!! In the 'translate' case, finish_angle is ignored, translation is
!! the direction given by finish_position-start_position and by a
!! distance equal to the length of finish_position-start_position.
!!
!! Note 1) that if variables start_position and end_position are
!! used to define a surface, then
!! input variable line_divisions must be explicitly set positive. \n
!! Note for geoq 2) the ordering of namelists corresponding to each
!! description of geometry follows that of the list below, ordered
!! after skylparameters if the presence of skylight has been 
!! flagged in beqparameters \n
!! Note for geoq 3) .dat-style files may be used to define a surface if
!! filedata=.TRUE., these are defined using namelist progfiles
!! which comes AFTER namelist datvtkparameters
  type, public :: datvtkparameters
   !> surface interaction with particles
   !! This may be set as a check when geoq reads the .ctl file
   !! that the namelists are in the right order
   !! - 'absorb' first wall (absorbing boundary)
   !! - 'invisi' transparent geometry, totally ignored
   !! - 'skylit' skylight
   !! - 'escape' escape boundary (loss)
   !! - 'errlos' beancan (error loss)
   !! - 'cutout' cutout (expected loss)
   character(len=80) :: description  !< .
  !> If .TRUE., then data file specifies geometry
  !! The default depends on which  code is being used.
  !! if datvtk, then the default is .TRUE.
  !! if geoq, then the default is .FALSE.
  !! If filedata becomes .TRUE. when called by geoq,
  !! there needs to be a matching namelist progfiles following
  !! datvtkparameters in the file,
  !! which refers to .dat-style file(s) to define a surface
  logical :: filedata !< -
   !> angle units, either radian(s or degree(s),
   !! default 'degree'
     character(len=20) :: angle_units  !< .
   !> length units, either mm or me(tres)
   !! default 'mm'
     character(len=20) :: length_units  !< .
   !> transform type,
   !! default  'rotate', rotate about central axis of torus,
   !! alternative is 'translate'
     character(len=80) :: transform_type  !< .
   !> starting  toroidal angle \f$ \zeta \f$,
   !! will be overridden when constructing a skylight object
   !! default  0
     real(kr8) :: start_angle !< .
   !> finishing  toroidal angle \f$ \zeta \f$,
   !! will be overridden when constructing a skylight object
   !! default  90
     real(kr8) :: finish_angle !< .
   !> start position 
   !! default is polar \f$ (1,0,0) \f$
     real(kr8), dimension(3) :: start_position !< .
   !> finish position 
   !! default is polar \f$ (1,1,0) \f$
     real(kr8), dimension(3) :: finish_position !< .
   !> number of divisions around axis of rotation
   !! or along line of translation,
   !! default  10
     integer(ki4) :: number_of_divisions !< .
   !> coordinate system positive as posang (q.v.),
   !!default and only valid entry is unity for cylindrical polars
     integer(ki4) :: coordinate_system !< .
   !> number of divisions, default  0
   !! Number of points inserted between each point in .dat file, EXCEPT
   !! If line_divisions > number of positions in the .dat file, then
   !! replace those points with positions uniformly spaced on the line between
   !! the first and last \f$(R,Z)$\f$ points in the file
     integer(ki4) :: line_divisions !< .
  end type datvtkparameters

!---------------------------------------------------------------------
!> Entries in namelist inputfiles.
!! There is a namelist inputfiles for each of the major component
!! programs, geoq, hdsgen and powcal. The namelist
!! specifies names of key input files in each case.
!!
  type, public :: inputfiles
   !> hdsgen only, geobj data input file.
   !! geoq and powcal only,  legacy vtk format geometry data input file,
   !! default 'null',
   !! must be ASCII text and note that only a restricted set of objects
   !! such as triangles and tetrahedra is accepted.
     character(len=80) :: vtk_input_file  !< .
   !> geoq only,  EQDSK G format data input file,
   !! default 'null'
     character(len=80) :: eqdsk_input_file  !< .
   !> geoq and powcal only,  launch data input file,
   !! default 'null', EXPERT only
     character(len=80) :: lau_input_file  !< .
   !> geoq only,  generic equilibrium input file,
   !! default 'null'
     character(len=80) :: equil_input_file  !< .
   !> geoq only,  generic mesh input file used by nucode, on which field defined
   !! default 'null'
     character(len=80) :: mesh_input_file  !< .
   !> powcal only,  vtk format geometry data input file,
   !! describing target or results geometry, default 'null'
     character(len=80) :: vtkres_input_file  !< .
   !> powcal only,  special format output file from geoq
   !! containing equilibrium field as splines, other field
   !! information and key parameters, default 'null'
     character(len=80) :: geoq_input_file  !< .
   !> powcal only,  special format output file from hdsgen
   !! describing HDS and containing mapping parameters, default 'null'
     character(len=80) :: hds_input_file  !< .
  end type inputfiles

!---------------------------------------------------------------------
!> Entries in namelist miscparameters.
!! There is a namelist miscparameters for each of the major component
!! programs, geoq, magtfm, powcal, smanal and vtktfm.
!! The namelist was intended as a catch-all
!! for parameters that did not easily fit into any other category, and
!! primarily intended for expert use. However, for  magtfm and vtktfm
!! miscparameters is (mis)used to set critical user quantities.
!!
!! For magtfm indeed, the namelist determines the transformation of
!! the field to the special sp3 format, which implies a direct product
!! cubic spline in \f$ R \f$ and \f$ Z \f$ together with a Fourier transform in angle.
!! The components of \f$ {\bf B} \f$ are ordered \f$ (B_R, B_Z, B_\zeta) \f$.
!! For computational efficiency the Fourier transform needs to have a cut-off
!! as explained in the published paper, and most of the controls are
!! truncation criteria.
!!
!! For vtktfm, it helps to know a little history, in that for early
!! ITER work, each panel was regarded as composed of two separate bodies
!! either side of a central attachment (to avoid self-shadowing). For
!! SMARDDA-PFC use, the option 'panel' should be set so that there is a one-to-one
!! correspondence between isolated geometry bodies and panels.
!!
  type, public :: miscparameters
   !> magtfm only,  \f$ R \f$ value at which \f$ B_T R \f$ defined,
   !! default  0.8,
   !! the two parameters together determine the amplitude of the
   !! vacuum field.
     real(kr8) :: radius_of_axis !< .
   !> magtfm only,  how vacuum field is determined,
   !! default  1 (inert)
     integer(ki4) :: vacuum_field_axis_calculation    !< .
   !> magtfm only,  order of spline,
   !! default  4 (implies cubic spline, other values untested)
     integer(ki4) :: spline_order    !< .
   !> magtfm only,  order of transform splitting,
   !! default  2 (must be two, implying sine and cosine components)
     integer(ki4) :: split_transform    !< .
   !> Number of repeats of magnetic field in toroidal direction for mapping to cylindricals
   !! \f$ N_{\zeta p} \f$, default 0
     integer(ki4) :: nzetap !< .
   !> Toroidal angle \f$ zeta \f$ of first point in samples (degrees) \f$ \zeta_{0} \f$,
   !! defaukt 370 (so ignored)
     real(kr8) :: zeta_start  !< .
   !> magtfm only,  minimum wavenumber vector for evaluation,
   !! default  \f$ (1,1,0) \f$, filters out constant contributions
   !! to vacuum field \f$ R \f$ and \f$ Z \f$ components
     integer(ki4), dimension(3) :: min_wavenumber !< .
   !> magtfm only,  maximum wavenumber vector for evaluation,
   !! default  \f$ (1,1,1) \f$. It is not sensible to choose
   !! a value of more than half the number of samples in angle.
   !! Truncation is controlled by the mode_cutoff parameter if
   !! the wavevector is set to \f$ (-1,-1,-1) \f$
     integer(ki4), dimension(3) :: max_wavenumber !< .
   !> magtfm only,  parity (1=even,2=odd,3=both) for each field component,
   !! default  \f$ (2,2,1) \f$ is the usual symmetry.
     integer(ki4), dimension(3) :: mode_parity !< .
   !> magtfm only,  relative mode cut-off amplitude,
   !! default  .0001
     real(kr8) :: mode_cutoff !< .
   !> magtfm only,  \f$ B_T R \f$  to help determine with specified \f$ R \f$,
   !! the amplitude of the vacuum field.
   !! default  1.
     real(kr8) :: i_requested !< .
   !> magtfm only, layout of data within .txt input file
   !!  - 0. default, MAST-U periodic spec
   !!  - 1. MAST-U periodic spec after Geoff
   !!  - 2. ITER periodic spec
   !!  - 3. Muddled periodic spec
   !!  - 10. MAST-U full 360deg spec (TO DO)
   !!  - 12. ITER full 360deg spec
   !!  - 14. ITER full 360deg spec staggered in toroidal angle
   !!  - 23. muddled and duplicated 360deg ITER case
   !!  - 33. Incorrectly stacked datasets in zeta ordering
   !!  - 34. The ITER 360deg case (stacked datasets in phi ordering)
   !!  - 42. Incorrectly assumed data has angular data segments patched together
   integer(ki4) :: data_layout !< .
   !> vtktfm only,  for defining panels and their transforms,
   !! default  'split_panels', i.e. each panel consists of two disjoint
   !! bodies of geometry, unless set to 'panel' when it is assumed that
   !! there is one-to-one correspondence between bodies and panels.
   !! 'tagged' is for EXPERT use only
     character(len=20) :: option !< .
   !> vtktfm only,  set .TRUE. to read namelist vtktfmparameters
     logical :: new_controls !< .
   !> vtktfm only,  number \f$ \geq \f$ number of files for which transform defined,
   !! default  1, needs to be \f$ \leq \f$ maximum_number_of_files
     integer(ki4) :: max_number_of_files !< .
   !> vtktfm only,  units, either radian(s) or degree(s),
   !! default  'radians'
     character(len=20) :: angle_units  !< .
   !> vtktfm only,  number \f$ \geq \f$ number of panels in geometry,
   !! default  1, needs to be \f$ \leq \f$ maximum_number_of_panels
     integer(ki4) :: max_number_of_panels !< .
   !> vtktfm only,  number \f$ \geq \f$ number of transforms defined,
   !! default  1, no direct restriction on number
     integer(ki4) :: max_number_of_transforms !< .
   !> vtktfm only, dimension of bods index array, default 1000
     integer(ki4) :: max_bods_index !< .
   !> vtktfm only, used to generate unique bods numbers over many files, default 100
     integer(ki4) :: max_bods_in_file !< .
   !> vtktfm only, bods remain distinct, default .FALSE.
     logical :: preserve_input_structure !< .
   !> vtktfm only,  number of panels for which bodies defined,
   !! default  0, needs to be \f$ \leq \f$ maximum_number_of_panels
     integer(ki4) :: number_of_panels !< .
   !> vtktfm only,  number of transforms defined,
   !! default  0, no direct restriction on number
     integer(ki4) :: number_of_transforms !< .
  end type miscparameters

!---------------------------------------------------------------------
!> Entries in namelist smanalparameters.
!! Cluster analysis is performed using the sort_key, and the results will 
!! typically be rekeyed using new_key. The demonstration case classifies
!! bodies (sort_key) into clusters according to their centroids and then rekeys
!! using the poloidal angle of the centroid (new_key).
!! Analysis proceeds over each cluster.
!!
!! The centre about which the angle is taken is defined as the centroid of
!! all the geometry defined as a set of points. The centre
!! may be redefined in \f$ (R,Z) \f$ coordinates, using set_r_central, nominal_r_central, etc.
  type, public :: smanalparameters
     character(len=80) :: analysis_mode !< only one mode available, default 'regular'
   !> key to be used to sort scalar field, default 'Body', 'Surf' is probably another possibility,
   !! note that the leading capital letter is significant, these are in principle arbitrary labels.
     character(len=80) :: sort_key !< .
     character(len=80) :: analyse_scalar !< identifier of scalar field to analyse, default 'Q'
     integer(ki4) :: number_of_clusters !< sets cluster size by dividing range of key, default  18
   !> New key for sorting, always an angle,
   !! inert unless rekey is .TRUE.
   !! Only values 'angle', 'poloidal', 'toroidal', 'userdefined'  are allowed, default  'null',
   !! and 'angle' implies poloidal angle (only the first 6 characters signify)
     character(len=80) :: new_key !<.
     logical :: rekey !<  Only use new key if .TRUE., default  .FALSE.
   !> user defined dictionary
   !! Zero entry terminates. For example to group bodies 1,2,3 and 4, set first four entries
   !! to 1, then say bodies 5 and 6 are grouped by setting entries 5 and 6 each to 2 etc.
   !! If there are more bodies than entries in the dictionary, then they are all grouped
   !! together with key equal to number of entries plus one.
     integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS)  :: dictionary !< .
   !> Statistics to be calculated
   !! Allowed values are 'intg  ', 'intgsq', 'max   ', 'min   ', 'maxabs', and 'minabs'.
   !! intg performs the integral as sum of scalar x (triangle area) and
   !! intgsq is the integral of scalar squared, the others should have obvious meanings
     character(len=80), dimension(MAX_NUMBER_OF_PARAMETERS) :: required_statistics  !< .
     logical :: set_r_central !< user sets nominal central \f$ R \f$ for poloidal angle analysis
     logical :: set_z_central !< user sets nominal central \f$ Z \f$ for poloidal angle analysis
     !> user defined value of central \f$ R \f$,
     !! note ignored unless set_r_central=.TRUE.
     real(kr8) :: nominal_r_central !< .
     !> user defined value of central \f$ Z \f$,
     !! note ignored unless set_z_central=.TRUE.
     real(kr8) :: nominal_z_central !< .
  end type smanalparameters

!---------------------------------------------------------------------
!> Entries in namelist skylparameters
!!
!! Skylight generation is triggered by setting either 
!! the logical switch skylight_flux_limits or skylight_objects>0
!! in namelist beqparameters (q.v.).
!! A skylight in this context is a surface
!! guaranteed to let in the 'light', strictly particle
!! power flowing from the midplane.
!! The skylight is an axisymmetric surface defined in
!! \f$ (R,Z) \f$ coordinates.
!!
!! If based on flux limits, it will be calculated on the basis
!! of the shadow geometry, optionally confined to the rectangle of which the
!! opposite corners are \f$ (R_{extmin},Z_{extmin}) \f$ and
!! \f$ (R_{extmax},Z_{extmax}) \f$. The algorithm is straightforward, namely
!! bin the geometry according to flux \f$ \psi \f$ and (in the lower divertor case)
!! select the point with highest \f$ Z \f$ in each bin to form part of the skylight.
!!
!! The difficulty is to specify the extent of the flux range over which to
!! perform the binning.
!! There may be two skylights, although
!! since most machines only have a lower divertor, it is likely there will
!! be only one. Conceptually each skylight has an inner and an outer
!! section defined by which side of the line of extremal flux
!! running through the plasma center lies. The method of skylight
!! calculation may vary depending which section is considered, with
!! the convention regarding the array index
!! 1. inner, lower, (1,1)
!! 2. outer, lower, (2,1)
!! 3. inner, upper, (1,2)
!! 4. outer, upper, (2,2)
!!
!! A set of skylight_objects X flat skylight surfaces may also be defined
!! using namelist datvtkparameters, specifically variables start_position
!! and finish_position.
!!
  type, public :: skylparameters
     logical :: skylight_lower !< use lower skylight, default .FALSE.
     logical :: skylight_upper !< use upper skylight, default .FALSE.
   !> special controls for skylight
   !! first entry=1, use launch points instead of geometry points
   !! second entry=1, force inner \f$ Z(\psi) \f$ monotone
   !! third entry=1, force outer \f$ Z(\psi) \f$ monotone
     integer(ki4), dimension(3) :: skylight_control !< .
   !> choose how to calculate the flux extent, default 2
   !! (1 = user input, 2 = based on geometry)
     integer(ki4), dimension(4) :: flux_extent_option !< .
     real(kr8), dimension(4) :: flux_bin_extent !< extent of bin in flux terms, default 0._kr8 (INERT)
   !> extent of skylight in flux terms, default 0._kr8
   !! the range of flux is centred about the value of flux at the X-point.
     real(kr8), dimension(4) :: flux_extent !< .
     logical :: set_r_extent !< setting extent in \f$ R \f$, default .FALSE.
     logical :: set_z_extent !< setting extent in \f$ Z \f$, default .FALSE.
     real(kr8) :: flux_extent_margin !< margin for flux extent calculation of skylight(s), default 0.1
     real(kr8) :: r_extent_min !< \f$ R_{extmin} \f$  minimum extent in \f$ R \f$, default 0.1
     real(kr8) :: r_extent_max !< \f$ R_{extmax} \f$  maximum extent in \f$ R \f$, default 0.1
     real(kr8) :: z_extent_min !< \f$ Z_{extmin} \f$  minimum extent in \f$ Z \f$, default 0.
     real(kr8) :: z_extent_max !< \f$ Z_{extmax} \f$  maximum extent in \f$ Z \f$, default 0.
     integer(ki4), dimension(4) :: number_flux_boxes !< number of binning boxes in calculation of skylight(s), default 20
     integer(ki4) :: number_of_divisions !< number of divisions in angle of skylight, default 20
     real(kr8) :: number_geoml_in_bin !< number of typical lengthscales in sample extent
   !> Distance that assumes intersection with skylight, in quantised units, in each direction
   !!  typically 1/1024 of geometry, default 1.
   !! overwritten by input skylight_tolerance in namelist powcalparameters
     real(kr8), dimension(3) :: intersect_tolerance !< .
  end type skylparameters
!---------------------------------------------------------------------
!> Entries in namelist vtktfmparameters.
!!
!! It helps to know a little history, in that for early
!! ITER work, each panel was regarded as composed of two separate bodies
!! either side of a central attachment (to avoid self-shadowing). For
!! SMARDDA-PFC use, the option 'panel' should be set so that there is a one-to-one
!! correspondence between isolated geometry bodies and panels. Further,
!! it became possible for the code to distinguish any integer-labelled attribute
!! in the input, e.g. using set process_by_name=Surf for surfaces, hence the
!! term 'bods' is now preferred.
!!
!! There are important differences in behaviour depending on whether
!! the main geometry input is a single file or not. In the former case
!! it is possible to split the file by a named attribute (typically 'Body')
!! and/or to displace each body by a numbered transform. Provided
!! there are fewer than max_bods_in_file=100 bodies, then it is not necessary to specify all
!! the bodies in panel_bodies, those omitted will be assigned the
!! identity transform. Note the variable make_same must be explicitly set
!! .FALSE., also then geometry which is not assigned to a body
!! or surface (Body "0" or zero) cannot be transformed.
!!
!! If there is more than one input file,
!! then by default each file is treated as a single object with a body with
!! number greater than 100. The numbering convention should be obvious
!! from inspection of the examples in the smiter-aux test deck,
!! the first file is numbered 101, the second 201, the third 301, etc.
!! The number of copies may be specified for each file, default 1
!! (i.e. a single occurrence of the geometry). Other copies of
!! file 1 are numbered 102, 103, 104, etc, of file 2 are numbered 202,
!! 203,204, etc. max_bods_in_file may be freely changed on input to
!! accommodate more complicated geometries. If preserve_input_structure
!! is set .TRUE., then the bods in the input files are (re)numbered
!! consecutively starting with unity, but there must not be more than
!! max_bods_index bods. \n
!! WARNING (1) unless the bodies are numbered consecutively setting
!! preserve_input_structure=.TRUE. will usually cause a crash error.\n
!! (2) Most
!! combinations of preserve_input_structure, extract, split_file and
!! make_same are untested and may well also cause vtktfm to crash.
!!
  type, public :: vtktfmparameters
   !> Option for defining panels and their transforms
   !! 'split' default
   !! 'panel' bodies in pairs to define panel
   !! 'tagged' bodies tagged with their transform number
     character(len=20) :: option !< .
   !> number \f$ \geq \f$ number of files for which transform defined,
   !! default  1, needs to be \f$ \leq \f$ maximum_number_of_files
     integer(ki4) :: max_number_of_files !< .
   !> units, either radian(s) or degree(s),
   !! default  'radians'
     character(len=20) :: angle_units  !< .
   !> number \f$ \geq \f$ number of panels in geometry,
   !! default  1, needs to be \f$ \leq \f$ maximum_number_of_panels
     integer(ki4) :: max_number_of_panels !< .
   !> number \f$ \geq \f$ number of transforms defined,
   !! default  1, no direct restriction on number
     integer(ki4) :: max_number_of_transforms !< .
   !> number of panels for which bodies defined,
   !! default  0, needs to be \f$ \leq \f$ maximum_number_of_panels
     integer(ki4) :: number_of_panels !< .
   !> number of transforms defined,
   !! default  0, no direct restriction on number
     integer(ki4) :: number_of_transforms !< .
   !> set .TRUE. to split file by attribute name, default .FALSE.
     logical :: split_file !< .
   !> set .TRUE. to make all attribute names the same, default .TRUE.
     logical :: make_same !< .
   !>  value for attribute to take (e.g. body number), default 1
     integer(ki4) :: same_value !< .
   !>  name of attribute to be split / homogenised, default Body
     character(len=80) :: process_by_name !< .
     logical :: paneltfm !< apply transform only if .TRUE. (default)
     logical :: extract !< extract objects according to key and limits, default .FALSE.
     character(len=80) :: extract_key !< key for extraction, default 'toroidal'
   !> dimension of bods index array, default 1000
     integer(ki4) :: max_bods_index !< .
   !> used to generate unique bods numbers over many files, default 100
     integer(ki4) :: max_bods_in_file !< .
   !> bods remain distinct, default .FALSE.
     logical :: preserve_internal_structure !< .
     real(kr8), dimension(2) :: plasma_centre !< centre of discharge in \f$ (R,Z) \f$
     real(kr8) :: minimum_angle !< minimum angle for extraction, default -22.5
     real(kr8) :: maximum_angle !< maximum angle for extraction, default 22.5


  end type vtktfmparameters

end module namelist
