!>\mainpage SMARDDA-PFC User Manual
!! 
!! 
!! \section sec1 What are SMITER and SMARDDA-PFC?
!! 
!! Software for predicting power deposition on
!! accurate representations of first wall geometry. \n
!! (It is NOT for plotting magnetic islands or other core magnetic field.)
!! 
!! Technically: \n
!! \arg SMITER is a user-friendly GUI wrapper for SMARDDA-PFC,
!! using salome to assist in geometry preparation.
!! The kernel of SMITER is a selection of UKAEA's
!! SMARDDA object-oriented Fortran-95 modules
!! for ray, particle and fieldline tracing.
!! SMARDDA-PFC is primarily a script-based interface to the same modules
!! plus a bash dialog interface \Bold{smpfcint} primarily for teaching purposes.
!! Both SMITER and SMARDDA-PFC use ParaView for visualisation.
!! SMARDDA-PFC depends on the open-source SMARDDA-LIB ray-tracing core,
!! as indicated in the figure by the f95 and f77 lib box.
!! \image html layout.png "Logical arrangement of the SMARDDA software"
!! \image latex layout.eps "Logical arrangement of the SMARDDA software" width=5cm
!!
!! Pronunciation  for SMITER is "Sm" plus however you
!! pronounce ITER, N.B. SMARDDA is pronounced smarter.
!! 
!! SMARDDA-PFC technical features
!! 
!! \arg Accurate field-line integration \n
!!   -User controlled tolerance to check convergence  \n
!!   -Cubic spline interpolation of magnetic field  \n
!!   -Option to use flux coordinates for speed and accuracy in simple flux geometries (more complex coming)  \n
!!   -Vacuum field may be 3-D  \n
!! \arg Local or global calculations \n
!!   -Limiter or divertor \n
!! \arg Toroidally periodic feature \n
!!   -Only necessary to specify 1/18 th geometry \n
!! \arg On-going development programme \n
!!   -More complex flux geometries  \n
!!   -ANSYS integration  \n
!!   -GPU implementation  \n
!!   -Light-weight provenance capture \n
!!   -Developer documentation
!! 
!>
!! \image html tube.png
!! \image latex tube.eps "Flux tube from midplane to surface of PFC" width=5cm
!! SMall physics, In Tokamaks Exact Representation
!! \arg Flux tube
!! \arg Power in at midplane=power out at surface of PFC
!!
!! So Minimal InpuTs to Enable Runs:
!! \arg Accurate representations of first wall geometry to be studied
!! \arg Magnetic field equilibrium
!! \arg Optionally, ripple magnetic field representation
!! 
!! \subsection sec1coord Coordinate Systems
!! So Many InTERpretations of Space
!! \arg Accurate first wall geometry originates in CATIA-TM
!! Computer-Aided-Design (CAD) geometry database,
!! coordinates (X,Y,Z) in millimetres (mm).
!! \arg Fieldline will be approximated by a set of straight line
!! segments, so best coordinates for following fieldlines
!! are ones in which they are straightest
!! 
!! 
!! 
!>
!! \image html xy.png
!! \image latex xy.eps "Cartesian XYZ coordinates" width=5cm
!! Cartesian \f$ XYZ \f$  coordinates
!! 
!>
!! \image html rz.png
!! \image latex rz.eps "Cylindrical polars" width=5cm
!! Cylindrical polars \f$ RZ\zeta \f$, N.B. not \f$ R\phi Z \f$. Distances in
!! metres (m), angles in radians. Also \f$ RZ\xi \f$.
!>
!! \image html circ.png
!! \image latex circ.eps "Toroidal coordinates" width=5cm
!! Toroidal coordinates \f$ r\theta\zeta \f$, angles in radians. Angle
!! measured from vertical through bottom X-point. Also \f$ r\theta\xi \f$.
!!
!>
!! \image html fcirc.png
!! \image latex fcirc.eps "Flux coordinates" width=5cm
!! Flux coordinates \f$ \psi\theta\zeta \f$, angles in radians. Angle 
!! measured from vertical through bottom X-point. Also \f$ \psi\theta\xi \f$.
!! 
!! Angle \f$ \xi \f$ is \f$ \zeta \f$ scaled so that \f$ \xi \f$ runs from 
!! \f$ -\pi \f$  to \f$ +\pi \f$ in a toroidal segment, which is a
!! fraction \f$ 1/N_{\zeta p} \f$ of
!! the device geometry (often \f$ N_{\zeta p} = 18 \f$).
!! 
!! 
!! \subsection sec1surf Defining Surfaces Accurately
!! So Many InTERpretations of Surfaces
!! \arg CAD databases contain curves and surfaces defined
!! as NURBS, a special sort of rational polynomial that
!! enables exact representation of for example the conic sections.
!! \arg Intersections of NURBS surfaces are only
!! approximately represented by NURBS curves,
!! implying small gaps.
!! \arg Small gaps avoided by use of triangular faceting of
!! every surface. \n
!! Lots of triangles to test for intersection with fieldline
!! segments, unless you get SMARDDA.
!! 
!! \subsection sec1dnull Single and Double Null Equilibria
!! Many magnetic equilibria with divertors are single null (SN),
!! meaning the plasma region is restricted by an X-point
!! assumed to lie below at the bottom of the plasma, as sketched.
!!>
!! \image html BtopoSN.png
!! \image latex BtopoSN.eps "Single Null Equilibrium" width=5cm
!! Sketch of single null (SN) equilibrium. 
!! The appearance of the X-point in the contours
!! of flux in the symmetry plane implies that the restrictive
!! flux surface is singular.
!!
!! A double-null (DN) equilibrium is one where the plasma region is not
!! only restricted by a bottom X-point, but also by a second X-point above.
!! These will typically correspond to separatrices that are close in space,
!! and it will be assumed that no PFCs lie in the region between them.
!!>
!! \image html Btopo1.png 
!! \image latex Btopo1.eps "Double Null Equilibrium" width=5cm
!! Sketch of double null (DN) equilibrium. X-point "1" always
!! lies at the bottom, even if the flux surface through the upper X-point
!! contains a smaller region of plasma.
!!
!! It is expected that plasma is always confined to the
!! smaller region restricted by an X-point, so that to maximise
!! region size, equilibria are sought where both bottom
!! and top X-points lie on the same (singular) flux surface. Such equilibria
!! are described as symmetric double null (SDN).
!!>
!! \image html BtopoSDN.png
!! \image latex BtopoSDN.eps "Symmetric Double Null Equilibrium" width=5cm
!! Sketch of symmetric double null (SDN) equilibrium. X-point "1" always
!! lies at the bottom.
!!
!! The different flux surface topologies indicated by the sketches
!! make possible different physical behaviour, eg. in the SDN, the
!! scrape-off layer (SOL) is split into two disconnected parts
!! inboard and outboard of the X-points' common major radius.
!! Since it is generally agreed that plasma escapes more easily
!! on the low field side (LFS) of the tokamak, it is likely that not only
!! that more of the lost power will go outboard, but that the
!! loss profile across the SOL will be different in the two parts.
!! In fact, referring to the sketch of DN, there may be up to four distinct
!! profiles characterising the loss at the midplane, viz.
!! 1. Inboard of both singular flux surfaces
!! 2. Inboard, between the two singular flux surfaces
!! 3. Outboard, between the two singular flux surfaces
!! 4. Outboard of both singular flux surfaces
!!
!! The behaviour in each may be separately specified as exponential,
!! double-exponential, Eich etc.
!! 
!! To determine which profile is appropriate for a given fieldline
!! (supposed to correspond to flux surface value \f$ \psi \f$),
!! values of flux \f$ \psi_{x1} \f$ and \f$ \psi_{x2} \f$ corresponding
!! respectively to the two X-point values are first computed in the
!! preprocessing stage. To determine the power deposited,
!! the fieldline \f$ \psi \f$ is then compared with
!! \f$ \psi_{x1} \f$ and \f$ \psi_{x2} \f$, enabling its midpoint intersection
!! to be determined as lying either in SOL parts 1 and 4, or in parts 2 and 3.
!! To determine whether the intersection is inboard or outboard, the
!! coordinates of the position
!! of the fieldline intersection \f$ (R,Z) \f$  with the first wall
!! are then examined. If the point is below the plasma centre,
!! ie. \f$ Z<Z_{cen} \f$, then \f$ R \f$ is compared with \f$ R_{x1} \f$,
!! else \f$ R \f$ is compared with \f$ R_{x2} \f$.
!! In each case, should \f$ R \f$ be smaller, it is supposed to lie inboard.
!!
!! In part \f$ i \f$ of the SOL, the power density is assumed to have power
!! dependence
!! \f{equation}{
!! F_i Q_i f_i(\Delta R_i)
!! \f}
!! Where \f$ Q_i  \f$ is a normalisation constant, \f$ f_i \f$ is the profile
!! function, \f$ F_i=1/2 \f$ is the fraction of power directed up or down
!! (except that in the DN case \f$ F_2=F_3= 1 \f$), and generalising
!! Appendix A of Arter et al (corrected (p)reprint at http://arxiv.org/abs/1403.7142)
!! \f{equation}{
!! \Delta R_i= R_i-R, \;\;i=1,2 \\
!! \Delta R_i= R-R_i, \;\;i=3,4 \\
!! \f}
!! For coding purposes it is convenient to write \f$ \Delta R_i \f$
!! in terms of the value of the flux \f$ \psi \f$ on the 
!! surface and the value \f$ \psi_i \f$ at \f$ R=R_i \f$
!! on the midplane.  Following Arter et al (with further typo correction)
!! \f{equation}{
!! \Delta R_i= \frac{\Delta \psi_i}{|\partial \psi/\partial R|_{R=R_i}}
!! = \frac{\Delta \psi_i}{R_{i}B_{pmi}}
!! \f}
!! where \f$ \Delta \psi_i= \psi-\psi_i \f$, and \f$ B_{pmi} \f$
!! is the absolute value of the poloidal magnetic field at
!! \f$ R=R_i \f$ on the midplane.
!! 
!! Following Arter et al, the power density must be normalised by the
!! total loss, let this be \f$ P_{li} \f$ in part \f$ i \f$, then
!! approximating an element of midplane area by 
!! \f{equation}{
!! dA = d\phi R_i d (\Delta R_i)
!! \f}
!! it follows from axisymmetry that for \f$ i=3,4 \f$ (with obvious
!! index changes for \f$ i=1,2 \f$)
!! \f{equation}{
!! 2\pi R_i Q_i \int^{R=R_{i+1}}_{R=R_i} d (\Delta R_i)
!! f_i(\Delta R_i) = P_{li}
!! \f}
!! where \f$ R_5=\infty \f$, or in practice the major radius of the PFC
!! at outer midplane,
!! \f{equation}{
!! Q_i = \frac{P_{li}}{ 2\pi R_i \int^{|R_{i+1}-R_i|}_0 dx f_i(x)}
!! \f}
!! Similarly for \f$ i=1,2 \f$
!! \f{equation}{
!! Q_i = \frac{P_{li}}{ 2\pi R_i \int^{|R_{i-1}-R_i|}_0 dx f_i(x)}
!! \f}
!! where \f$ R_0=0 \f$ or in practice the major radius of the PFC
!! at inner midplane. For Eich the situation is slightly different as there is a
!! contribution from inside the PFR which may be accounted for by 
!! changing the lower limit of integration from \f$ 0 \f$ to \f$ -\infty \f$.
!!
!! For the exponential profile, \f$ f_i(x) = \exp{(-x/\lambda_q)}\f$,
!! and its integral may be recovered as a special case \f$ \rho=0 \f$ of the double
!! exponential \f$ f_i(x) = \exp{(-x/\lambda_q)} + \rho \exp{(-x/\lambda_n)}\f$,
!! which has integral
!! \f{equation}{
!! \lambda_q (1 - \exp{(-R'/\lambda_q)} + \rho \lambda_n (1-\exp{(-R'/\lambda_n)}
!! \f}
!! where \f$ R'=|R_{i+1}-R_i| \f$ or \f$ R'=|R_{i-1}-R_i| \f$.
!! The Eich profile is
!! \f{equation}
!! f_i(x)=\exp \left(s^2 -\frac{x}{\lambda_q} \right)
!! \text{erfc} \left( s-\frac{x}{2 \lambda_q s} \right)
!! \f}
!! where \f$ s=\frac{\sigma }{2\lambda_q } \f$,
!! and \f$ x=\frac{\Delta \psi }{R_m B_{pm}} \f$,
!! and the erfc function is defined as
!! \f{equation}
!! \text{erfc}(z) = 1- \frac{2}{\sqrt{\pi}} \int^z_0 \exp{-x^2}dx
!! \f}
!! The indefinite integral of the Eich profile is
!! \f{equation}
!! -\lambda_q \exp \left(s^2 -\frac{x}{\lambda_q} \right)
!! \text{erfc} \left( s-\frac{x}{2 \lambda_q s} \right)
!! +\lambda_q \text{erfc} \left( -\frac{x}{2 \lambda_q s} \right)
!! \f}
!! from which the normalisation follows by substituting \f$ x=R' \f$.
!!
!! One issue when there is more than one region is lack of continuity of 
!! the power deposition profile. In the case of four regions, this
!! implies that 
!! \f{equation}
!! Q_2 = \frac{F_1 Q_1 f_1(0)}{F_2f_2(R_2-R_1)}\\
!! Q_3 = \frac{F_4 Q_4 f_4(0)}{F_3f_3(R_4-R_3)}
!! \f}
!! 
!! \subsection sec1alg SMARDDA Algorithm
!! Reduces number of triangles to test by:
!! \arg Recursive division of 3-D space into cuboidal boxes, aka octree.
!! Recursion stops when number of triangles in a box aka bin
!! falls below a preset value. (Always know which bin your fieldline
!! segment is in.)
!! \arg Track segment from box to box by binary tricks, thanks to use of
!! quantised coordinates. The quantum lengths \f$ h_i \f$ are chosen so
!! that \f$ L_i=1024h_i \f$ or similar power of two. (The biggest box has
!! size \f$ L_1L_2L_3 \f$. There is a single precision restriction to exponent \f$ <32 \f$.)
!! \arg Jargon -  an octree is an example of an Hierarchical Data
!! Structure (HDS), designed to Help Do Surfaces quickly.
!! 
!! \subsection sec1morecoord Coordinates Summary
!! SuMmary of InTeRpretations of Space
!! 1. Cartesian \f$ XYZ \f$ mm from CAD
!! 2. Cylindrical polars \f$ RZ\zeta \f$ metres and radians (\f$ RZ\xi \f$)
!! 3. Toroidal polars \f$ r\theta\zeta \f$ metres and radians (\f$ r\theta\xi \f$)
!! 4. Flux coordinates \f$ \psi\theta\zeta \f$ (\f$ \psi\theta\xi \f$),
!! for cases of simple flux topology using involving limiters.
!! 5. Quantised coordinates \f$ xyz \f$ or \f$ 123 \f$, versions of 2 and 3 above
!! 
!! \subsection sec1work Workflow
!>
!! \image html flow.png
!! \image latex flow.eps "SMARDDA-PFC flow chart" width=10cm
!!
!! 
!! Steps MakIng TrianglEs Ready
!! 1. CAD database input -> CAD defeaturing of part body -> CAD
!! "healing" -> Meshing -> Mesh output of part body. \n
!! (.CATpart file
!! -> .dat file, produced by specially trained ITER Design
!! Office staff)
!! 2. Conversion of .dat format (proprietary MSC Nastran or
!! a project specific 'gnu dat' format \ref sec2formats) to
!! legacy .vtk format (OpenSource) by SMARDDA-PFC datvtk.
!! 3. Combination of part bodies to form PFC geometry under study, by
!! SMARDDA-LIB vtktfm.
!! 4. Conversion of STL and obj files to legacy .vtk format (OpenSource) by
!! ParaView, scripts stl2vtk and obj2vtk provided.
!! 5. Use of ParaView "Sources" menu  to produce simple objects in legacy .vtk,
!! format, other objects by meshing constant surfaces of 3-D scalar data.
!! 6. Check using ParaView that surface normals point outward, and if
!! necessary correct (a body at a time) using datvtk.
!!
!! Aside: the UKAEA approach: .CATpart .x_t or .step file imported by
!! CADfix-TM on Windows for user-controlled defeaturing, healing,
!! orient normals and meshing. cfmcn(p) converts CADfix proprietary format to .vtk
!! 
!! Second, Magnetic InpuT as Eqdsk Representation
!! SMARDDA-PFC allows
!! 1. For magnetic equilibrium after Solovev, but mostly
!! 2. Input as .eqdsk file.
!! \arg A. EQDSK G format is a "non-standard" standard for
!! solutions \f$ \psi(R,Z) \f$ of the Grad-Shafranov equation
!! \arg B. Owing to EFIT's use of a Fortran output
!! descriptor that need not have spaces between numbers.
!! \arg C. but SMARDDA-PFC expects reasonable sampling rate for flux.
!! 
!! 
!! SMARDDA-PFC and .eqdsk
!!
!! The basic standard for EQDSK G may be found at:
!! https://fusion.gat.com/theory/Efitgeqdsk (with password)
!! or try https://w3.pppl.gov/ntcc/TORAY/G_EQDSK.pdf
!! SMARDDA-PFC has been programmed to read many (more than seven)
!! different eqdsk formats, and negate flux \f$ \psi \f$ where necessary
!! to try to ensure that the field has the same handedness as ITER (the
!! direction of \f$ B \f$ is not important for SMARDDA-PFC).
!! The .eqdsk header defines the mesh-size, also the flux on axis
!! and at boundary. SMARDDA-PFC will try to use this information,
!! but it does no consistency checking on the read
!! contents (at UKAEA Glasser's dcon code is
!! (mis)used for this purpose, NB dcon assumes Zmid=0).
!! There are .eqdsk files with muddled signs,
!! factors of \f$ 2\pi \f$ in the flux, wrong helicity of field etc.
!! The geoq code can help resolve at least some of these issues, see \ref sec2geoq
!! producing a lengthy set of diagnostics which need inspecting carefully,
!! and may imply the need to change inputs to achieve the desired configuration.
!! One feature which currently confuses the software is the presence
!! of blank lines in the .eqdsk file, so the user will have to arrange
!! separately to detect and remove such lines.
!! 
!! \subsection sec1cli Command line interface
!!
!! (If SMARDDA-PFC has been installed in directory ~/smardda/smardda-pfc, the SMARDDA-PFC
!! codes may be added to the user's path with the shell commands: \n
!! export SMITER_DIR=~/smardda/smardda-pfc;source $SMITER_DIR/Extras/setup
!! )
!!
!! \arg The standard way to execute codes is \n
!! \Emph{code .ctl file}, \n
!! e.g. \Emph{geoq wall}, where file wall.ctl contains
!! inputs controlling geoq run, except
!! \arg datvtk was designed as a filter, so instead \n
!! datvtk -[opt] ".dat file", or datvtk -c ".ctl file" \n
!! In the latter case, the .ctl file
!! contains instructions to generate surfaces by rotating
!! straight lines defined by their end-points, and the points
!! may be defined by files in the gnu dat format (\ref sec2formats).
!! 
!! Aside: \Emph{opt=ABC}, where \Emph{A=v} -> .vtk input; \Emph{B=o} -> orient
!! triangles relative to first in file, \Emph{B=s}->shell tetrahedra (BROKEN, use ParaView
!! instead);\Emph{B=f} -> flip triangle orientation; \Emph{C=d} -> divide into
!! separate block body files, \Emph{C=m} ->suppress block
!! information, e.g. datvtk -xom FW34.dat (x as spacer)
!! 
!! The .ctl files
!! contain Fortran namelists which must appear in a
!! strict order although the variables within a namelist may be
!! (re)defined in any order.
!! Variables are documented using doxygen 1.8, which produces
!! .html files for input to firefox or more generally xdg-open,
!! which are linked to the \Bold{smpfcint} script \ref sec2smpfcint.
!! 
!! Example .ctl file input to geoq \n
!! &inputfiles \n
!! vtk_input_file='wall.vtk', \n
!! eqdsk_input_file='16_97s.eqdsk' \n
!! / \n
!! &miscparameters \n
!! / \n
!! &plotselections \n
!! plot_geoqx = .true., \n
!! plot_geoqvolx = .false., \n
!! / \n
!! &beqparameters \n
!! beq_rmove=-6., \n
!! beq_cenopt=4, \n
!! beq_psiopt=1, \n
!! beq_bdryopt=3, \n
!! beq_thetaopt=1, \n
!! /
!! 
!! Namelist &inputfiles  obvious \n
!! Namelist &miscparameters 
!! for future/expert use etc. \n
!! Namelist &plotselections 
!! almost obvious, see index of outputs in \ref sec1output \n
!! Namelist & beqparameters 
!! go to doxygen documentation headed "Data Types List", 
!! also accessible directly using
!! firefox or xdg-open ./smardda/smardda-pfc/doc/srcdoc/html/classes.html
!! (index under "Data Fields" sub-heading)
!! 
!! \subsection sec1output SMARDDA-PFC output files
!!
!! All the codes produce a log file, which contains typically time-and-date stamps,
!! information about the files used, a selection of key parameters, execution times and
!! at end a summary of errors and warnings, terminating with "END OF LOG FILE".
!! Output to the terminal (Fortran print command) consists of the code
!! header, errors associated with log-file writing and a small sub-set of the log-file output,
!! notably serious errors, and the execution times if the code
!! completes normally.
!!
!! The suite is set up so as to terminate codes when a serious error
!! occurs, other errors are classified as warnings, and execution continues.
!! Most also produce a .out file, which is primarily used for inter-code
!! communication, and also acts as a "lock-file", in that the same .ctl
!! file cannot be re-run until the corresponding .out file is deleted.
!! (Good practice is to change the .ctl file-name if the contents change.)
!!
!! Index of suffix strings of output files, either with
!! terminating .vtk or graphics (.png, .ps) suffixes.
!!
!! The .vtk files begin with the geometry definition, then the quantities in angled brackets,
!! which are cell-centred rather than point quantities (unless stated otherwise) in case this is
!! an option. J is the Jacobian of the mapping to flux coordinates (so R/J is the right-hand
!! side of the ODE representing the fieldline), \f$ I = R B_T \f$,
!! Bcart is B in Cartesian coordinates, Msign is mask (0=shadowed) with the sign of B.n
!! (where n is the surface normal), \f$ \psi_{start} \f$ is \f$ \psi \f$
!! at the nominal start of the fieldline (i.e. at the end which is NOT the PFC intersection).
!! L=lenpath and objhit are taken from trackx file header line described below.
!! Mapped geometry refers either to flux-mapped coordinates or cylindricals
!! depending on the type of calculation.
!! In the \Bold{smpfcint} single directory,
!! hdsgen and powcal outputs acquire initial distinguishing "_hds" and "_pow" suffices respectively.
!! \arg _angnu < statistics of field per geometry unit output in gnuplot format
!! \arg _angnure < statistics of rekeying of field per geometry unit output in gnuplot format
!! \arg _anx < statistics of field per geometry unit output on Cartesian geometry
!! \arg _ansmx < statistics of field per geometry unit output at centroids only in Cartesian geometry
!! \arg _beqvolx < \f$ \psi (point) \f$, Bcart (point) > geoq output on volume Cartesian geometry
!! \arg _geobjq < surfaces > hdsgen output in quantised geometry
!! \arg _geofldx < XYZ \f$ \psi (point)  \f$ Bcart (point)  Normal Bcell BRZeta > geoq output on Cartesian geometry (N.B. \f$ \zeta \f$ component corresponds to ParaView \f$ Z \f$.)
!! \arg _geoptq < points > hdsgen output in quantised geometry
!! \arg _geoqm < XYZ Bcart > geoq output on mapped geometry
!! \arg _geoqvolm <\f$ (0,J,I) (point) \f$> geoq output on volume flux-mapped geometry
!! \arg _hdsm < HDS > hdsgen output in mapped geometry
!! \arg _hdsq < HDS > hdsgen output in quantised geometry
!! \arg _powlenx < Q Msign \f$ \psi_{start} \f$ L objhit > _powx file augmented using \Bold{addlenvtk}
!! \arg _powm < Q \f$ \psi_{start} \f$ L Msign > powcal output on mapped geometry
!! \arg _powstatm < Q-avg Q-dev > powcal output on mapped geometry (EXPERT)
!! \arg _powstatx < Q-avg Q-dev > powcal output on Cartesian geometry (EXPERT)
!! \arg _powx < Q \f$ \psi_{start} \f$ L Msign > powcal output on Cartesian geometry
!! \arg _RZsil < points \f$ \psi \f$ > \Bold{geoqgnu} output in RZ geometry
!! \arg _silm < points \f$ \psi \f$ > \Bold{geoqgnu} output in flux-mapped geometry
!! 
!! The following files have no prefix, because unless stated, they are 2-D plots
!! of flux-related quantities, and of course the fieldline tracks only depend
!! on the geometry to the extent of their end-points, files :
!! \arg beq+wall <\f$ \psi \f$ silhouette> \Bold{geoqgnu} output in RZ geometry
!! \arg dR <\f$ \partial  \psi /\partial R \f$> \Bold{geoqgnu} output in RZ geometry
!! \arg dZ <\f$ \partial  \psi /\partial Z \f$> \Bold{geoqgnu} output in RZ geometry
!! \arg R <\f$ R \f$> \Bold{geoqgnu} output in flux-mapped geometry
!! \arg RJ <\f$ R/J \f$> \Bold{geoqgnu} output in flux-mapped geometry
!! \arg RZ <\f$ \psi \f$> \Bold{geoqgnu} output in RZ geometry showing flux-mapped region
!! \arg RZoom <\f$ \psi \f$> \Bold{geoqgnu} output in RZ geometry confined to flux-mapped region
!! \arg trackm.\Emph{n}01.vtk <element \f$ (n-1) \f$  fieldline> powcal output in mapped geometry
!! \arg trackx.\Emph{n}01.vtk <element \f$ (n-1) \f$  fieldline> powcal output in Cartesian geometry
!! \arg Z <\f$ Z \f$> \Bold{geoqgnu} output in flux-mapped geometry
!!
!! In track files, the header line contains \n
!! \arg \f$ elt=n \f$, i.e. 1+(identifier of triangular element line starts at)  
!! \arg sub=1 (EXPERT) 
!! \arg lenpath length of fieldline track in file in mm (accurate in trackx files only) \n
!! \arg objhit code, \n
!! >0 1+(element hit) \n
!! 0 no collision, ODE integration limits reached \n
!! -1 track left computational domain \n
!! -2 track hit termination plane (termplane) \n
!!
!! \subsection sec1resource SMARDDA-PFC Resources
!! 
!! \arg ITER git repository
!! \arg ./smardda/smardda-pfc/doc/srcdoc/html as above, includes links to SMARDDA papers on arXiv
!! \arg ./smardda/smardda-pfc/Extras \n
!!   .sgo files  scripts covering entire workflow \n
!!   .ctl files  skeleton control files used by .sgo files \n
!! \arg Test-deck files - smiter-aux  (tests.bash)
!! https://git.iter.org/projects/BND/repos/smiter-aux/browse
!! are data files (geometry and equilibria) used by the .sgo scripts.
!! The test-deck should be regarded as a valuable resource, providing
!! examples that can quickly be modified to address new problems.
!! \arg ParaView Guide - a printed version for ParaView 4.3
!! may be purchased via the web-site http://www.paraview.org, or an electronic version
!! of the latest release (4.4/5.0 as of March 2016) is a free download.
!! For hints for using ParaView with SMARDDA-PFC see \ref sec2paraview
!! 
!! \subsection sec1install SMARDDA-PFC Installation
!! 
!! First install SMARDDA-LIB which is available from CCFE gitlab or
!! IMAS, by
!! \arg Clone or pull from IMAS git repository
!! http://git.iter.org/projects/BND/repos/smardda/smardda-lib/browse
!! and follow instructions in README, noting that
!! \arg produce executables and optionally documentation (script will enquire)
!! by typing in the smardda/smardda-lib top directory \n
!! source ./exec/build intel \n
!! gfortran and the debug options  intel_dbg and gfortran_dbg are also available
!! \arg Ensure that for subsequent sessions environment variable
!! SMALIB_DIR is set correctly as returned by build, as this is where
!! SMARDDA-PFC will search for library code.
!!
!! Next install SMARDDA-PFC similarly
!! \arg Clone or pull from IMAS git repository
!! http://git.iter.org/projects/BND/repos/smardda/smardda-pfc/browse
!! and follow instructions in README, noting that
!! \arg produce executables and optionally documentation (script will enquire)
!! by typing in the smardda/smardda-pfc top directory \n
!! source ./exec/build intel \n
!! gfortran and the debug options  intel_dbg and gfortran_dbg are also available
!! \arg SMARDDA was developed in 2007 onwards, whereas the ITER-hpc cluster
!! runs Centos vintage 2005, hence special
!! modules are needed.
!! \arg ParaView needs a one-off re-setting of defaults (for
!! better \Bold{smpfcint} integration)
!! \arg SMARDDA-PFC has a test-deck (data files for which are in the smiter-aux repo),
!! see \ref sec3 for more details, which should be used to check the installation.
!!
!! It may be helpful to note also that
!! \arg the commands \n
!! ./exec/decompile and ./exec/compile \n
!! are provided for use after modifying source and/or when changing compilers
!! \arg after the initial build, the compiler may be changed by going to
!! the ./config directory and creating a new config.inc file (link), thus \n
!! cd ./config;ln -sf config_new-compiler.inc config.inc
!! 
!!
!! \section sec2 Description of Constituent Codes and Scripts
!!
!! \ref sec1work  indicates how the codes are combined to perform a
!! shadowing calculation, and the \Bold{smpfcint} script \ref sec2smpfcint
!! leads the user implicitly through this process. Each of the codes
!! is briefly described below. As described in the introduction, they are
!! driven ultimately by namelists in .ctl files. The namelists expected in the
!! respective .ctl files are listed below. Precise descriptions of 
!! the usage of the variables in each
!! namelist may be referenced in the online documentation, as data types
!! at ./smardda/smardda-pfc/doc/srcdoc/html/classes.html.
!!
!! \subsection sec2ctlgen ctlgen code
!! This code is not expected to be run directly by the user.
!!	
!! \subsection sec2datvtk datvtk code
!! datvtk was designed as a filter, so unlike the other main executables,
!! it accepts options immediately after the name, e.g. \n
!!  datvtk -[opt]  .dat file \n
!!  datvtk -c .ctl file \n
!! In the latter case, the .ctl file
!! contains instructions to generate surfaces by rotating
!! straight lines defined by their end-points, optionally these
!! points may be defined in separate gnu dat format files (\ref sec2formats).
!! Options: \Emph{opt=ABC}, where \n
!! \Emph{A=v} -> .vtk input\n
!! \Emph{B=o} -> orient triangles\n
!! \Emph{B=s} -> shell tetrahedra (BROKEN, use ParaView  instead)\n
!! \Emph{B=f} -> flip triangle orientation\n
!! \Emph{C=d} -> divide into separate block body files\n
!! \Emph{C=m} -> suppress block information \n
!! e.g. datvtk -xom FW34.dat (x as spacer).
!!
!! datvtk can also be as a filter to convert to STL format: \n
!! \Emph{C=b} INERT, output in STL format (divide into separate block body files) \n
!! \Emph{C=t} output as triangles in STL format (suppress block information) \n
!! e.g. datvtk -vxt FW34.vtk (x as spacer).
!!
!! Note that datvtk accepts only the .dat Small and Large Field Formats,
!! and not the Free Field Format, as described
!! in the MSC Nastran 2013 Quick Reference Guide. (Unfortunately this
!! document may be obtained only after
!! registering with the MSC Nastran web-site.) There is a further restriction to
!! the small set of "bulk data entries" output by CATIA mainly to define
!! simple meshes. Beware that CATIA is also capable of producing a more
!! compact, non-conforming .dat file output, which SMARDDA-PFC will not read.
!!
!! Namelist order in .ctl is as follows
!! \arg progfiles 
!! \arg datvtkparameters 
!! \arg plotselections 
!! 
!! \subsection sec2extrascripts  Extra scripts
!! The following, which reside in subdirectory Extras, are mentioned at
!! various points in the code description:
!! \arg addlenvtk take header lines from trackx, combine and add to end
!! of _powx.vtk -> _powlenx.vtk file
!! \arg brt_pvcalc.txt polar coordinates of B for input to ParaView calculator 
!! \arg chmodx ensures commands are executable after installation
!! \arg geoqgnu/ngeoqgnu plot with gnuplot .gnu files produced by geoq
!! for ITER/MAST
!! \arg obj2stl, obj2vtk, stl2vtk format conversions using pvpython
!! \arg setup load linux modules appropriate to system (if identified)
!! \arg trtognu convert track files to .gnu files
!! \arg vtkpng/cl generate png file from vtk/in close-up (only works for JET)
!! \arg jwall.txt, mwallr.txt, mwall.txt, wall.txt silhouette files for
!! JET, MAST-U upper Z, MAST-U lower Z, ITER
!! \arg EXTR.ctl -> vtktfm,  extract part of .vtk format geometry
!! \arg SENS.ctl -> vtktfm,  transform parts of .vtk format geometry
!! \arg SMAN.ctl -> smanal, analysis of results defined in .vtk file by row or otherwise
!! \arg VTKF.ctl -> vtktfm,  transform .vtk format geometry from mm to m \n
!! The use of other .ctl and .sgo files by the test-deck is described in detail
!! in the README for smiter-aux.
!!
!! \subsection sec2geoq geoq code
!! The geoq code may be used to analyse the flux function \f$ \psi \f$,
!! usually defined in a .eqdsk file. (A good of effort has been exerted to
!! ensure that SMARDDA-PFC can read the many variants of the EQDSK G "standard".)
!! To aid understanding, geoq produces as output .gnu files which may be plotted using the
!! \Bold{geoqgnu} script (from the smardda/smardda-pfc/Extras directory). Running geoq makes it possible to
!! \arg Produce contour plots of flux gradients (\f$ \propto B \f$ components)
!! \arg Overlay \f$ \psi \f$ contours with silhouette wall.vtk, and also 
!! CAD geometry as \f$ (R,Z) \f$ points
!! \arg Obtain the extrema of flux on geometry \f$ \psi_{ltr} \f$, flux values from .eqdsk file \f$ \psi_{q} \f$ and the flux
!! value at nearest X-point \f$ \psi_{X} \f$ as labelled output in the .log file.
!! This output the user
!! should inspect, and if necessary set up new input file with appropriate \f$ \psi_{m} \f$.
!! Sometimes it may be necessary to define a special search box to
!! locate the correct X-point, by setting beq_xsearch=1 in namelist beqparameters.
!! The .ctl file also controls whether inboard (HFS)  or outboard (LFS) values of
!! \f$ R \f$ and \f$ B \f$ are used in the power
!! calculation (to account for the effects of flux expansion) via input
!! beq_bdryopt in namelist beqparameters.
!! 
!! It is also possible to define the flux function analytically, if the
!! the variable equil_input_file has a suffix .ana, this file must
!! contain namelists as set out below.
!! The suffix .equ is also recognised and corresponds to an output
!! format related to that produced by the FIESTA code, see example
!! the MAST-U test case \ref sec3T2.
!!
!! A further addition is the ability to add a 'skylight' for divertor
!! discharges, ie. a set consisting of surfaces which if fieldlines reach 
!! any of them, the start triangles on the PFCs are guaranteed 
!! not to be shadowed, usually sitting well
!! below the midplane in diverted discharges. It will often
!! be better to ask additionally for the production of a skylight
!! based on flux limits of the launch points.
!!
!! Namelist order in .ctl as follows
!! \arg inputfiles 
!! \arg miscparameters 
!! \arg plotselections 
!! \arg beqparameters 
!! \arg skylparameters  (if skylight_flux_limits TRUE  or skylight_objects>0 in beqparameters)
!! \arg datvtkparameters x skylight_objects set in beqparameters
!!
!! Namelist order in .ana as follows
!! \arg equilparameters 
!! \arg meshparameters 
!!
!! \subsection sec2powcal powcal code
!! This is the code that actually performs the power surface deposition
!! calculation, following fieldlines using transformed geometry from geoq and the HDS from
!! hdsgen. There is the capability to select from a range of power
!! deposition profiles, and even define new ones, provided by namelist
!! edgprofparameters (see the documentation for namelist edgprofparameters).
!! For global calculations, a wide range of fieldline
!! termination criteria may be set using namelist termplaneparameters,
!! including multiple simultaneous criteria, useful for private flux region
!! calculations.
!!
!! Partly to indicate that execution is continuing successfully, powcal
!! outputs to the terminal, the numbers of elements for which the starting
!! direction is uncertain, together with the values of B.n in the
!! Cartesian co-ordinate system and in the mapped system which because of
!! their different signs gave rise to the condition (these elements are
!! regarded as shadowed). powcal local or "escape" calculations give rise
!! to warnings in the .log file that "Object not found in binary tree",one
!! for each element that is thus regarded as illuminated or "wetted".
!!
!! powcal may also produce the track files described in \ref sec1output if
!! the right plotselections are made. The utility \Bold{addlenvtk} is provided
!! to process \f$ L \f$, connection length information, and append to the corresponding
!! _powx.vtk file (the root of which must be specified), producing
!! a file root_powlenx.vtk. Alternatively the _powx.vtk file contains \f$ L \f$ as a
!! scalar field defined for all triangles.
!!
!! Namelist order in .ctl as follows
!! \arg inputfiles 
!! \arg miscparameters 
!! \arg plotselections 
!! \arg powcalparameters 
!! \arg (termplaneparameters if termination_planes=T)
!! \arg (edgprofparameters if more_profiles=T)
!! \arg odesparameters 
!! 
!! \subsection sec2hdsgen hdsgen code
!! This computes the multi-octree HDS introduced in \ref sec1alg,
!! which is designed to accelerate the computation of track/ray-triangle intersection.
!! The user should not normally need to be concerned with details of
!! hdsgen operation, except for the need to increase the parameter limit_geobj_in_bin for
!! geometries with a large number of triangles, i.e. greater than
!! approximately 100000 triangles, depending on their degree of clustering.
!!
!! Namelist order in .ctl as follows
!! \arg inputfiles 
!! \arg hdsgenparameters 
!! \arg btreeparameters 
!! \arg positionparameters 
!! \arg plotselections 
!! 
!! \subsection sec2smpfcint smpfcint script
!! \Bold{smpfcint} is the main GUI for SMARDDA-PFC, designed to help the beginning user:
!! \arg Assemble geometry from meshed parts or bodies
!! \arg Run and analyse simple problems interactively
!! and start:
!! \arg Constructing scripts for investigations of problem sequences
!! \arg Editing .ctl files to access more advanced features
!! 
!! Primarily \Bold{smpfcint} is a bash 3.0 script using Linux dialog, with features
!! \arg Self-contained as an introduction to SMARDDA-PFC
!! \arg Point-and-click is optional
!! \arg Help on dialog usage from top-level menu
!! \arg Firefox start-up with links to web-based documentation is also available from within script
!!
!! (bash 3.0 is described in the book Learning the bash shell 3rd Edition, by
!! C.Newham and B.Rosenblatt, O'Reilly 2005)
!! 
!! \Bold{smpfcint} start
!! \arg Create directory and preferably copy geometry
!! and equilibrium files into it (although you can
!! copy later): Linux  \Emph{mkdir},  \Emph{cp} and  \Emph{cd} commands.
!! \arg Attach (\Emph{cd}) to directory, type \Bold{smpfcint}
!! \arg First time, click on Help button \n
!! Note: Ctrl/C (undocumented) can get you out of
!! problems on menus below top menu
!! 
!! \Bold{smpfcint} finish
!! \arg \Bold{smpfcint} as it finishes tells you
!! 1. Where it logged the codes it ran, this file is potential starting-point
!! for a control script (contains !* where errors occurred)
!! 2. Which scratch directory it used (usually ./ismtemp),
!! contains files logging \Bold{smpfcint} script input data for
!! reuse and potentially provenancing.
!! \arg Directory now full of input and output files.
!! vtk and png files may be visualised at any time
!! using ParaView, or use eog for png.
!! \arg To visualise fieldlines, the recommended approach is to use ParaView to
!! select the geometry on which they terminate, see \ref sec2paraview for details.
!!
!! In the example below, the shadow geometry was assembled from files 
!! partshad1_R_modifier.dat \n
!! partshad1_Z_modifier.dat \n
!! partshad2_R.dat \n
!! partshad2_Z.dat \n
!! beancan_RZ.dat \n
!! results.dat, plus transformations
!! 
!! and the target or results geometry from results.dat.
!! 
!! The beancan is a simple geometrical device to trap any anomalous fieldline
!! integrations within a "beancan", which should lie outside and preferably
!! enclose the volume of interest.
!! 
!! The resulting files in the \Bold{smpfcint} directory were: \n
!! beancan.ctl control file for datvtk execution \n
!! beancan.log file logging output of datvtk \n
!! beancan_RZ.dat point data used to define a surface by datvtk \n
!! beancan.vtk ASCII file defining geometry as set of oriented triangles \n
!! beq+wall.png 2-D plot of equilibrium contours and wall silhouette \n
!! beq+wall.ps 2-D plot of equilibrium contours and wall silhouette \n
!! cmdwrap\Emph{nnn}.log start-point for script construction \n
!! dR.png 2-D plot of contours of R derivative of flux from \Bold{geoqgnu} \n
!! dR.ps 2-D plot of contours of R derivative of flux from \Bold{geoqgnu} \n
!! dZ.png 2-D plot of contours of Z derivative of flux from \Bold{geoqgnu} \n
!! dZ.ps 2-D plot of contours of Z derivative of flux from \Bold{geoqgnu} \n
!! HR.eqdsk equilibrium in EQDSK G format \n
!! ismfiles.sav records principal geometry and equilibrium file(s) used by \Bold{smpfcint} \n
!! ismform.sav records principal parameters used by \Bold{smpfcint} \n
!! ismtemp scratch directory \n
!! partshad1_modifier.ctl control file for vtktfm execution \n
!! partshad1_modifier.log file logging output of vtktfm \n
!! partshad1_modifier.vtk ASCII file defining geometry as set of oriented triangles \n
!! partshad1_R_modifier.dat point data used to define a surface by datvtk \n
!! partshad1_Z_modifier.dat point data used to define a surface by datvtk \n
!! partshad2.ctl control file for datvtk execution \n
!! partshad2.log file logging output of datvtk \n
!! partshad2_R.dat point data used to define a surface by datvtk \n
!! partshad2.vtk ASCII file defining geometry as set of oriented triangles \n
!! partshad2_Z.dat point data used to define a surface by datvtk \n
!! results.ctl control file for geoq execution \n
!! results.dat mesh data converted to .vtk surface format by datvtk \n
!! results_geofldx.vtk geoq output of geometry and magnetic field in mapped coordinates \n
!! results_geoqm.vtk geoq output of geometry in mapped coordinates \n
!! results_geoq.out output file from geoq \n
!! results_gnu_dR.gnu file for input to \Bold{geoqgnu} to produce dR graphics \n
!! results_gnu_dZ.gnu file for input to \Bold{geoqgnu} to produce dZ graphics \n
!! results_gnu.gnu file for input to \Bold{geoqgnu} to produce beq+wall graphics \n
!! results_gnu_cart.gnu file for input to gnuplot with all 3 Cartesian components of B-field \n
!! results_gnusil.gnu file for input to \Bold{geoqgnu} to produce RZsil graphics \n
!! results.log file logging output of geoq \n
!! results_pow.ctl control file for powcal execution \n
!! results_pow.log file logging output of powcal \n
!! results_pow_powcal.out output file from powcal \n
!! results_pow_powx.vtk output of power deposition on Cartesian geometry \n
!! results_tfm1.ctl control file for vtktfm execution \n
!! results_tfm1.ctlin control file for ctlgen execution \n
!! results_tfm1.log file logging output of vtktfm \n
!! results_tfm1.vtk vtktfm output of geometry after transformation in Cartesians \n
!! results_tfm2.ctl control file for vtktfm execution \n
!! results_tfm2.ctlin control file for ctlgen execution \n
!! results_tfm2.log file logging output of vtktfm \n
!! results_tfm2.vtk vtktfm output of geometry after transformation in Cartesians \n
!! results.vtk ASCII file defining geometry as set of oriented triangles \n
!! results_RZsil.png 2-D plot of nodes of geometry with flux contour background from \Bold{geoqgnu} \n
!! results_RZsil.ps 2-D plot of nodes of geometry with flux contour background from \Bold{geoqgnu} \n
!! shadow.ctl control file for geoq execution \n
!! shadow.ctlin control file for ctlgen execution \n
!! shadow_geofldx.vtk geoq output of geometry and magnetic field in mapped coordinates \n
!! shadow_geoqm.vtk geoq output of geometry in mapped coordinates \n
!! shadow_geoq.out output file from geoq \n
!! shadow_gnu_dR.gnu file for input to \Bold{geoqgnu} to produce dR graphics \n
!! shadow_gnu_dZ.gnu file for input to \Bold{geoqgnu} to produce dZ graphics \n
!! shadow_gnu.gnu file for input to \Bold{geoqgnu} to produce beq+wall graphics \n
!! shadow_gnu_cart.gnu file for input to gnuplot with all 3 Cartesian components of B-field \n
!! shadow_gnusil.gnu file for input to \Bold{geoqgnu} to produce RZsil graphics \n
!! HR_eqbdry.gnu file for input to \Bold{geoqgnu} to produce _bdry graphics \n
!! HR_eqltr.gnu file for input to \Bold{geoqgnu} to produce _bdry graphics \n
!! shadow_hds.ctl control file for hdsgen execution \n
!! shadow_hds_geobjq.vtk hdsgen output of geometry in quantised space \n
!! shadow_hds_geoptq.vtk hdsgen output of geometry as points only  in quantised space \n
!! shadow_hds_hdsgen.out output file from hdsgen \n
!! shadow_hds_hds.hds file defining HDS produced by hdsgen \n
!! shadow_hds_hdsm.vtk hdsgen geometric output of HDS in mapped coordinates \n
!! shadow_hds_hdsq.vtk hdsgen geometric output of HDS in quantised space \n
!! shadow_hds.log file logging output of hdsgen \n
!! shadow.log file logging output of geoq \n
!! shadow.vtk ASCII file defining geometry as set of oriented triangles \n
!! shadow_RZsil.png 2-D plot of nodes of geometry with flux contour background from \Bold{geoqgnu} \n
!! shadow_RZsil.ps 2-D plot of nodes of geometry with flux contour background from \Bold{geoqgnu} \n
!! equil_bdry.ps 2-D plot of silhouette of geometry and plasma boundary from \Bold{geoqgnu} \n
!! table.dat scratch files used by \Bold{geoqgnu} \n
!! tabledR.dat scratch files used by \Bold{geoqgnu} \n
!! tabledZ.dat scratch files used by \Bold{geoqgnu} \n
!! trackm.\Emph{n}01.vtk track in mapped coordinates corresponding to element \Emph{n-1} \n
!! trackx.\Emph{n}01.vtk track in Cartesian coordinates corresponding to element \Emph{n-1} \n
!! wall.txt (R,Z) coordinates of ITER silhouette, cold wall
!! 
!! \subsection sec2magtfm magtfm code
!! The vacuum magnetic field is defined in a .txt file, in what will be referred to as mag format. Details of the 
!! required contents may be deduced by inspection of the example file 
!! ./smiter-aux/Data/Equilibrium/testov5.txt, given that:
!!
!! mag format data has to overlap in toroidal angle  \f$ \zeta \f$, i.e. field data for the 
!! first point in  \f$ \zeta \f$  is repeated on the last point.  Note that  \f$ \zeta \f$  is 
!! the negative of the usual cylindrical angle coordinate, so that  \f$ \zeta \f$  
!! increasing from zero implies the toroidal sector in negative \f$ Y \f$
!! adjoining \f$ X=0\f$. The mesh used is a direct product mesh, i.e. specified 
!! by one set of sample points in each coordinate direction  \f$ \zeta \f$, \f$ R \f$ and \f$ Z \f$.
!! The sampling of  \f$ \zeta \f$  must be uniform, units radians, the \f$ R \f$ and \f$ Z \f$ sample 
!! points are specified in metres, one sample value per line. Field \f${\bf  B}\f$ is 
!! given in Cartesian components \f$ (X,Y,Z) \f$, one vector per line. The \f${\bf  B}\f$ 
!! vectors are ordered with  \f$ \zeta \f$  varying the fastest, \f$ R \f$ next and then \f$ Z \f$ 
!! varying the slowest. Each of the 4 items ( \f$ \zeta \f$  grid, \f$ R \f$ grid, \f$ Z \f$ grid, \f${\bf  B}\f$) in the 
!! file is separated from the next by a blank line, each item starts with 
!! a description line, followed by the array size, except in the case of \f${\bf  B}\f$.
!!
!! Namelist order in .ctl as follows
!! \arg magfiles 
!! \arg plotselections 
!! \arg miscparameters 
!!
!! \subsection sec2vtktfm vtktfm code
!! The code transforms and/or combines .vtk format surface geometry files (ASCII format).
!! In the \Bold{smpfcint} script, vtktfm is used either to apply a solid body
!! transformation or to combine files, but in fact as shown by the validation
!! cases, the code can simultaneously transform and combine, although no more
!! than one transformation on any body in one code usage.
!! It can also or instead  extract geometry according to angle constraints,
!! useful for cutting down geometry to  a periodic cell.
!!
!! Namelist order in .ctl as follows
!! \arg miscparameters 
!! \arg vtktfmparameters (only if new_controls=.TRUE. in miscparameters)
!! \arg vtkfiles 
!! \arg panelarrayparameters 
!! \arg positionparameters 
!! 
!! \subsection sec2smanal smanal code
!! The code analyses the _powx file output by powcal, using the information on the
!! structure of the geometry in the results geometry input file.
!! Analysis, e.g. max or integral, is conducted on the basis of the geometrical unit specified.
!! Usually 'Body' information is defined in the input file, and the scalar 'Q'
!! from the _powx file is analysed, for e.g. local maxima or integrals.
!! The results may be displayed in .vtk format at cell centroids or as a table/.csv file
!! suitable for input to gnuplot. A refinement for reactor work is the capability to group
!! the geometry further in poloidal angle, so that by suitable choice of
!! cluster parameter, e.g. the max Q on a particular row of tiles may be
!! plotted/tabulated. The _gnure.gnu file provides details about the re-keying process.
!!
!! To obtain useful plots with _ansm.vtk files, the point size
!! should be increased to 10 or so. To get labelling to appear, it is
!! necessary to make a selection of cells using the icons immediately
!! above the plot window, then check the appropriate label under the "Cell
!! labels" drop-down menu. The default double precision labelling is not
!! very readable, to change click on the gearwheel adjacent to "Selection
!! Color", and alter the "Cell label format" to %.3e for example.
!! (Note for Fortran programmers, these are C printf formats, which are
!! perverse in that the type specifier (e, d or g) comes at the end and a % at beginning.)
!!
!! Namelist order in .ctl as follows
!! \arg analysisfiles
!! \arg miscparameters 
!! \arg plotselections 
!! \arg smanalparameters 
!! 
!! \subsection sec2paraview ParaView and SMARDDA-PFC
!! The initial chapters of the official guide explain basic concepts such as the 
!! graphics pipeline and should be at least skim-read before reading these
!! SMARDDA-PFC-specific hints. (As of 5.2, there is an initial pop-up window
!! with information for the beginner.) This preliminary should help the user  acquire
!! the necessary familiarity with basic ParaView GUI 
!! layout, knowledge of when to use the "Apply" button, and thereby avoid
!! the commonest mistakes and sources of confusion.
!!
!! General points worth noting are that
!! \arg ParaView is continually being developed, and the hints here assume
!! version 4.4 or later.
!! \arg SMARDDA-PFC in normal use only generates surface information, so many ParaView
!! functions which require volume data are not available or needed.
!! \arg ParaView starts counting from zero,
!! whereas SMARDDA-PFC is more conventional, so that for example, the first track
!! starts from triangle zero.
!! \arg Toggling the gear-wheel symbol gives access to a wider and often useful
!! range of extra functionality.
!! 
!! The functionality most often needed for SMARDDA-PFC concerns the making and saving
!! of selections of the data in the .vtk files.
!! A specific example is the **visualisation of fieldlines emanating from part
!! of the geometry.** 
!! The selection process is activated by toggling the "Select Cells on" icon
!! immediately above the main "render" window ("Select Cells With Polygon"
!! offers finer control). A selection box is defined by the mouse and the
!! selected cells (triangles) are highlighted. The selection may be inspected
!! in more detail, by use of the Selection Display Inspector
!! (View -> Selection Display Inspector). If satisfactory, the "Extract
!! Selection" icon which is on the menu bar just above the Select Cells icons,
!! should be clicked upon. File->Save data can then be used
!! to save the selection to a Legacy vtk (.vtk) format file, using the 
!! options to save all timesteps to the Ascii file type. The new .vtk file
!! thereby produced can be used as an input to a new SMARDDA-PFC run with
!! field-line plotting activated in the .ctl file.
!! datvtk preserves block information from each .dat file when it
!! is converted to .vtk (unfortunately it is discarded in subsequent
!! processing). If block numbers are known, or identified from
!! ParaView visulation, their geometry can be selected by first
!! using Edit -> "Find data", where "is one of" is the most flexible
!! selection option. Having created the selection use the "Extract 
!! Selection" icon as above, clicking on "Copy active selection".
!! 
!! Most other functionality is accessed by the Filters menu. Of the myriad
!! options, the ones most relevant to SMARDDA-PFC are
!! \arg Integrate variables, for total power.
!! \arg Calculator, operates on fields such as Bcart and Q, to give B.n etc.
!! \arg Point data -> Cell data, and vice versa, since Calculator only operates
!! on variables of the same data type.
!! \arg Cell centres, when producing spread-sheet .csv output.
!! \arg Generate Ids, to help locate elements and points on the surface.
!! \arg Append attributes, enables data from different vtk files with same
!! geometry to be combined.
!! \arg Plot on intersection curves, for line plots of data as functions
!! of arclength on surface, see "Line plots" below.
!! 
!! Views in the render window  may be saved to .png file for plotting
!! using File->Save Screenshot.
!! The same view may be applied to a separate session by using File->Save State,
!! then File->Load State in the new ParaView session (same release required).
!!
!! Spread-sheets may be produced using File->Save Data, for which .csv is
!! the default option. For point data this is straighforward, as the
!! position is saved with the data (consistent with ParaView usage, the 
!! first vector component is zero). Cell centre points need to be
!! explicitly defined by the user with the Cell centre filter before using
!! File->Save Data with cell data such as Q.
!!
!! Line plots using ParaView
!!
!! In ParaView notation, SMARDDA-PFC calculates Q values at cell centres, but ParaView
!! only plots point values, so it interpolates in a manner than does not seem to
!! documented officially. Checks indicate that the point value is an average over
!! cell-centred values for all triangles meeting at the point (probably
!! area-weighted), hence shadow edges get blurred since typically only 2 cells
!! out of 6 or 8 have non-zero Q.  The intersection curve is likely to
!! be generated by the intersection of a plane and the surface, and
!! when generically the intersection between the plotting plane does not
!! align with the surface triangulation points, there is a further interpolation
!! to points on the intersection. The intersection curve is not sampled uniformly
!! but at points which seem to correspond to the projections of the vertices
!! of intersected triangles. \n
!! Hints for line plots\n
!! It is often easier to reduce the size of the geometry in the window
!! by using  "select cells through" icon just above followed by  \n
!! Extract selection, (on bar just above) then \n
!! Cell data-> point data \n
!! The filter is "Plot on intersection curves", when
!! the menu brought up is comprehensive if the gear-wheel is clicked.
!! However it will often be found helpful to save the data to a .csv file
!! for graphing using gnuplot (comment out first line with leading "#", replace
!! commas with spaces) or spreadsheet software.
!!
!! \subsection sec2sensit Sensitivity studies using SMARDDA-PFC
!!
!! A file SMLAB.ctl is provided in directory $SMARDDA-PFC_DIR/Extras to enable smanal
!! to produce an _ansmx.vtk file which displays labels 
!! of bodies or surfaces useful for sensitivity work. 
!! The string RES needs editing so the name of the geometry file appears instead.
!! To visualise the cell data in the _ansmx.vtk file which provide the block/body numbers, it will be
!! necessary to view as "Points" in ParaView and to select "Point Size" of approx. 5.
!! It may be necessary to launch the ParaView Selection Display Inspector from "View"
!! in ParaView at this point.
!! To see the values associated with the cell centres, make the 
!! selection "Select Cells Through" (N.B. **Cells**) of all the points, then
!! use the menus in the Selection Display Inspector to view the 
!! appropriate cell label, "max", a good cell label format is "%4.3g".
!! (May be set using the gear-wheel in the Selection Display Inspector
!! which gives boxes for setting "Selection Label Properties".)
!!
!! The file SENS.ctl provided in directory $SMARDDA-PFC_DIR/Extras indicates the
!! transformations that may be applied as part of a sensitivity analysis.
!! The string RES needs editing so the name of the geometry file appears instead.
!! The surface/body numbers displayed in the SMLAB/_ansmx.vtk file 
!! may be input to variable panel_bodies so as to apply the transform given by 
!! the corresponding entry to variable panel_transform, see the documentation
!! for the namelists in SENS.ctl.
!! Note that the example requires that there be at least 18 bodies. A further
!! point is that surfaces/bodies to which different transformations are applied
!! must be disjoint, or strange effects will result at their overlap.
!!
!! \subsection sec2formats Format of gnu dat and .vtk file headers
!!
!! Special header lines have been introduced for gnu dat and .vtk files 
!! to encode information about the geometrical data contained, as to
!! which coordinate system and units are used therein. In a few cases, it
!! may help the user to know the format of these headers, hence a
!! description follows. There are numerous examples in the test deck
!! 
!! For gnu dat files, the header should begin with '#', then \f$29\f$ characters
!! of free format, non-blank description terminated by space. After the
!! space is the string 'Number_Parameters=', followed by the number of
!! integer parameters in the rest of the line. Currently, this number is
!! six in Fortran I3 integer format, followed by numbers in 1X,I4 format, meanings as follows
!! 1. default 1
!! 2. descriptive code for geometry, see geobj for full list, 
!! eg. 0. geometry, 4. skylight, 5. terminating geometry, 6. beancan, 7. cutout.
!! 3. units (consistent with posang, -3 for mm, 0 for m)
!! 4. coordinate system, after posang, namely \n
!!   0 - position in Cartesians \f$ (x,y,z) \f$ \n
!!   1 - position in cylindrical polars \f$ (R,Z,\zeta), |\zeta| \leq \pi \f$ \n
!!   2 - position in flux coordinates  \f$ (\psi,\theta,\zeta) \f$ \n
!! 5. default 0
!! 6. default 0
!!
!! For .vtk files, the header is actually the second line, starting with
!! \f$30\f$ characters of free format, non-blank description terminated by space. After the
!! space is the string 'Number_Parameters=', followed by the number of
!! integer parameters in the rest of the line. Currently, this number is
!! six in Fortran I3 integer format, followed by numbers in 1X,I4 format, meanings as follows
!! 1. infilelevel (level of refinement in file), default 1
!! 2. indicates triangles geometry coded as to type
!! 3. units (consistent with posang, -3 for mm, 0 for m)
!! 4. coordinate system, after posang, namely \n
!!   0 - position in Cartesians \f$ (x,y,z) \f$ \n
!!   1 - position in cylindrical polars \f$ (R,Z,\zeta), |\zeta| \leq \pi \f$ \n
!!   2 - position in flux coordinates  \f$ (\psi,\theta,\zeta) \f$ \n
!! 5. quantised (if 1, else 0)
!! 6. (number of transforms applied) - (number of inverse transforms applied)
!!
!!
!! \section sec3 Test-Deck
!! 
!! \subsection sec3intro Introduction
!! 
!! The test-deck, including the validation test cases,  is designed to
!! demonstrate correct functioning of all the commonly used parts of
!! SMARDDA-PFC. Data for the test-deck is in the smiter-aux directory,
!! together with files "_notes.txt" on how to produce geometry for the validation
!! test cases.
!!
!! The .sgo scripts have extensive internal documentation (q.v.).
!! It is here noted that the .sgo scripts in ./smardda/smardda-pfc/Extras arrange
!! files in a subdirectory structure within the main run directory
!! for each case, unlike \Bold{smpfcint} where files are created directly in one
!! run directory. These .sgo subdirectories are named M, G, S, H and P
!! and correspond to execution of the code with the corresponding
!! initial letter, thus "G" for geoq, "P" for powcal, etc., except that
!! since geoq has to be run separately for results and shadowing geometry,
!! "S" is used for the latter.
!!
!! \subsection sec3cases Test and Validation Cases
!! A brief description of each case follows so that the
!! user may determine which is closest to the problem of current
!! interest and therefore might be most easily modified to solve it. 
!! Test cases from the SMARDDA development are performed
!! in directories starting "Test-", whereas SMARDDA-PFC validation cases
!! begin "Val-". The remaining components of the names in order are
!! an identifier of the equilibrium, the shadowing geometry and
!! the results geometry respectively. An optional 4th component
!! indicates whether a local or global problem was solved.
!! Thus the geometries may be easily identified among the .vtk files
!! in ./smiter-aux/Data/VTK and visualised using ParaView, although note that
!! the .vtk files for the validation cases have to be generated from the
!! .dat files in ./smiter-aux/Data/Geometry.
!!
!! For the two ITER test cases, the first wall silhouette is used to determine
!! the position of the contact point, hence the value of flux \f$ \psi \f$
!! of the last closed flux-surface (LCFS) and thence the power normalisation.
!! In the MAST-U test case \ref sec3T2 , the LCFS is known to correspond
!! to \f$ \psi=0 \f$. In the validation work, the \f$ \psi \f$ of the LCFS
!! is taken from the eqdsk file, except for \ref sec3V1g, the LCFS is 
!! determined using the displaced geometry, hence is closer to the centre
!! of the plasma increasing the peak power deposition by a factor of two.
!! In the original VDE validation case, there was an issue regarding
!! mesh resolution, hence the examination of a much finer grid.
!!
!! \subsubsection sec3T1 Test-EQ3-inrshad1-inres1
!! \image html 1.png
!! \image latex 1.eps "Power deposition on inres1 geometry" width=10cm
!! This test-case is taken from the AFWS project,
!! https://idm.f4e.europa.eu/?uid=24NDFU&version=v1.2 
!! which saw the development
!! of SMARDDA primarily for ITER limiter problems.
!! It consists of a typical ITER equilibrium with inner (high-field side HFS) limiting,
!! with contact on the left-hand body of panel 3, and shadowing by the 8
!! adjacent semi-panels.
!! Flux coordinates are used and this is a local or escape calculation,
!! i.e. if the fieldlines do not hit the shadowing geometry before they leave
!! the volume covered by the HDS, they are assumed to connect with the midplane.
!! Note that since the shadowing geometry is disjoint from the target geometry,
!! the calculation relies on an accurate silhouette of the vessel
!! wall to obtain the limiter flux value.
!!
!!
!! \subsubsection sec3T2 Test-MNEW8-sha12-res12
!! \image html 2.png
!! \image latex 2.eps "Power deposition on res12 geometry" width=10cm
!! This test-case is taken from the MAST-U project, which saw the development
!! of SMARDDA for divertor modelling in the presence of strong ripple field,
!! here defined in file testov5.txt.
!! The axisymmetric contribution consists of a conventional MAST DN
!! (double-null) discharge.
!! The effect on the planned MAST-U divertor geometry, represented by
!! a 30 degree toroidal sector is examined. The results res12 geometry is tiles T1 and T2
!! adjacent to the MAST centre-column and the power deposition profile is a 
!! simple exponential without Eich's modification.
!! Cylindrical coordinates allowing for non-axisymmetry are used and this is
!! a global calculation. However to reduce computational cost, fieldlines are
!! only followed up to Z=-1.29, whence it is assumed that they will act as
!! conduits of power from the midplane.
!!
!! \subsubsection sec3T3g Test-EQX40-shado14r-resf14r-3g
!! \image html 3.png
!! \image latex 3.eps "Power deposition on resf14r geometry" width=10cm
!! This test-case is taken from the AFWS project, which saw the development
!! of SMARDDA primarily for ITER limiter problems,
!! but this example demonstrates just how small the deposited power is if
!! the discharge is detached, as in this case. The plot shows the (semi-)panel with maximum
!! power deposition, the right-hand of panel 14, subject to shadowing by 4 adjacent rows 
!! each represented by up to 5 panels (there are gaps corresponding to the
!! equatorial ports). This is a global calculation using cylindrical
!! coordinates, but the geometry is not replicated to cover 360 degrees,
!! performed to validate the axisymmetric field solver introduced in the
!! original smiter development project.
!!
!! \subsubsection sec3T3l Test-EQX40-shado14r-resf14r
!! \image html 4.png
!! \image latex 4.eps "Power deposition on resf14r geometry" width=10cm
!! This test-case is taken from the AFWS project, which saw the development
!! of SMARDDA primarily for ITER limiter problems,
!! but this example demonstrates just how small the deposited power is if
!! the discharge is detached.
!! The difference between this case and that of \ref sec3T3g is that this is a local
!! calculation, using flux coordinates.
!!
!! \subsubsection sec3V1l Val-F11-fw264t-fw4t-1l
!! \image html 5.png
!! \image latex 5.eps "Power deposition on fw4t geometry" width=10cm
!! This represents the power deposited on slightly misaligned panel 4 geometry
!! by an equilibrium as described in Kocan et al., Nucl. Fusion 55(2015)033019(16pp).
!! http://dx.doi.org/10.1088/0029-5515/55/3/033019
!! Shadowing is by 11 panels on each of rows 2-6, centred in toroidal angle
!! on the target.
!! This is a local calculation using flux coordinates.
!!
!! \subsubsection sec3V1g Val-F11-fw264t-fw4t
!! \image html 6.png
!! \image latex 6.eps "Power deposition on fw4t geometry" width=10cm
!! The difference between this case and that of \ref sec3V1l is that this
!! calculation uses cylindrical coordinates and is a
!! global calculation requiring fieldlines to rise to a z=+3 "midplane"
!! (strictly termination plane or termplane).
!! The difference in peak power is explained by the different assumption
!! used to determine the LCFS (by the displaced geometry here).
!!
!! \subsubsection sec3V3 Val-vde-bot-baf
!! \image html 7.png
!! \image latex 7.eps "Power deposition on baf geometry" width=10cm
!! This involves a VDE equilibrium on grid of \f$ 65 \times 129 \f$ which is limited
!! low down, on the outer baffle, see Carpentier-Chouchana et al.
!! https://user.iter.org/?uid=D43N8F
!! This is a global calculation, using cylindrical coordinates,
!! in a 20 degree toroidal sector, i.e. 18-fold symmetry.
!!
!! \subsubsection sec3V3sy Val-vde-botsy-baf
!! \image html 8.png
!! \image latex 8.eps "Power deposition on baf geometry" width=10cm
!! The difference between this case and that of \ref sec3V3 is that this
!! calculation uses a 60 degree toroidal sector, i.e. sixfold symmetry.
!!
!! \subsubsection sec3V4 Val-HR-bott-baff
!! \image html 9.png
!! \image latex 9.eps "Power shadow on baff/align geometry" width=10cm
!! This uses a VDE equilibrium on a fine grid of \f$ 257 \times 513 \f$
!! originally from Y.Gribov https://user.iter.org/?uid=2V2XYR,
!! supplied by Carpentier-Chouchana in January 2016.
!! It is used on a geometry identical to that of \ref sec3V3 except without
!! misalignment, in a global calculation, using cylindrical coordinates,
!! in a 20 degree toroidal sector, i.e. 18-fold symmetry.
!!
!! \subsubsection sec3sm Val-smanal-mjg
!! \image html 10.png
!! \image latex 10.eps "Max power per row on reactor geometry" width=10cm
!! This processes data from a reactor model, further details of which
!! may be deduced from the references https://idm.euro-fusion.org/?uid=2N8RMN (SOF equilibrium) and 
!! https://idm.euro-fusion.org/?uid=2MJGX9 for the PFC geometry, which is a
!! full 20 degree toroidal sector albeit with a simplified divertor.
!! To visualise the cell data in the _ansmx.vtk file to produce the above plot, it will be
!! necessary to view as "Points" in ParaView and to select "Point Size" of approx. 5.
!! To see the values associated with the cell centres, make the
!! selection "Select Cells Through" (N.B. **Cells**) of all the points and view the
!! appropriate cell label, e.g. "max", a good cell label format is "%4.3e".
!! (May be set using the gear-wheel in the Selection Display Inspector
!! which gives boxes for setting "Selection Label Properties".)
