#! /bin/bash --norc
## script for SMITER data processing to be called by ismiter
#
## Steps 1 & 2 are assumed completed using the geometry menu item
## (In other workflows Step 1 would consist of simplifying and
## meshing the geometry using a CAD package. In this workflow, this step
## is regarded as performed by a Drawing Office to give a set of CATIA dat files.)
## 1. Convert elemental file meshes from dat to vtk format
## 2. Combine the elemental vtk files to give
## the shadowing ($gshad) and results ($gres) meshes.
## 2.5 Set up ctl files
## 3. Run geoq code to generate geometry and equilibrium data
## 3a. the option to calculate boundary psi from silhouette is not used
## 3b. process shadow geometry
## 3c. process launch aka results geometry
## 4. Run hdsgen code to generate HDS to accelerate line following
## 5. Run powcal code to produce power deposition profiles.
#
## directories
## HS set to smiter top directory
if [[ -z "$HS" ]] ; then
 if [[ -n "$SMITER_DIR" ]] ; then HS=$SMITER_DIR; else HS=$HOME/smiter; fi
fi
if [ ! -d $HS ] ; then
  echo "Neither directory \$SMITER_DIR=$SMITER_DIR nor $HOME/smiter exists - installation corrupt ?"
  echo "Script quitting"; exit
fi
export HS
## directory for interactive scripts
INT=$HS/int
#INT=~/dial # fix for testing
## scratch directory for script
tempdir=./ismtemp
## scratch file for editing .ctl files
temp=./ismtemp/tempctl
## Top level controls
## 1+0 = set up ctl files for first time only
## 0+1 = continue execution
## 2+2 = force overwrite of ctl and all execution outputs
sw=2+2
## global or local calculation
cal=local
## Inputs, defaults
## A. delivered EQDSK file and identifier/name for processed version
eqfile=16_97s.eqdsk
eqid=eqm
## B. geometry files
gres=fw264t
gshad=fw264t
gwall=wall  # not an input keyword, file containing ITER wall silhouette
## HDS number of objects in a bin
nbin=20
## C. control parameters
## set beq_fldspec (see namelist beqparameters documentation)
flds=3
## set beq_bdryopt (see namelist beqparameters documentation)
bopt=2
## boundary options: in, out or x
bdryl=in # now a label since $bopt is set
## power deposition parameters
decl=.05
powe=5.e+06 # note $powe is in W, whereas $pow in MW
## ripple
vfld=
arip=0
mrip=0
## toroidal geometry is periodic with this mode number
nzetp=1
#
## set variables from list key=arg
for i in $*; do 
key=${i%%=*}
arg=${i##*=}
case "$key" in
  tempdir) tempdir=${arg};;
  sw) sw=${arg//+/ };;
  flds) flds=${arg//+/ };;
  cal) cal=${arg//+/ };;
  bopt) bopt=${arg//+/ };;
  eqfile) eqfile=${arg//+/ };;
  eqid) eqid=${arg//+/ };;
  gres) gres=${arg//+/ };;
  gshad) gshad=${arg//+/ };;
  gwall) gwall=${arg//+/ };;
  bdryl) bdryl=${arg//+/ };;
  powe) powe=${arg};;
  decl) decl=${arg};;
  vfld) vfld=${arg//+/ };;
  arip) arip=${arg//+/ };;
  mrip) mrip=${arg//+/ };;
  nzetp) nzetp=${arg//+/ };;
  nbin) nbin=${arg//+/ };;
  ?) echo "unknown key";exit 2;;
esac 
done
## derived parameters
runid="$eqid"-"$gshad"-$gres
rundir=$PWD
runlog=$tempdir/$runid.log
## Logging
echo "Start script" $(date) > $runlog
echo "field_mapping = " $flds >> $runlog
echo "calculation_type = " $cal >> $runlog
echo "boundary_option = " $bopt >> $runlog
echo "equilibrium = " $eqfile >> $runlog
echo "equil_id = " $eqid >> $runlog
echo "results_geometry = " $gres >> $runlog
echo "shadow_geometry = " $gshad >> $runlog
echo "(wall_geometry) = " $gwall >> $runlog
echo "boundary_type = " $bdryl >> $runlog
echo "power_loss = " $powe >> $runlog
echo "decay_length = " $decl >> $runlog
echo "vacuum_field = " $vfld >> $runlog
echo "ripple_amplitude = " $arip >> $runlog
echo "ripple_modenumber = " $mrip >> $runlog
echo "geometry_toroidal_copies = " $nzetp >> $runlog
echo "objects__grouped_in_HDS = " $nbin >> $runlog
echo "run identifier = " $runid >> $runlog
## Input OK and logged, start by creating $tempdir
if [ ! -d $tempdir ] ; then rm -f $tempdir; mkdir $tempdir;fi
#
## Step 2.5 set up data files
swa=(${sw/+/ })
#echo "swa" ${swa[0]} ${swa[1]} ${swa[2]}
if [ ${swa[0]} != 0 ] ; then
## Hidden file protects ctl files unless first switch is 2
  if [ ${swa[0]} == 2 ] || [ ! -e ".STOP-"$runid"-ctl" ] ; then
  ## Step 2.5-3b  Shadow geometry
    rm -f $gshad.ctl
    sed -e "s/EQDSK/$eqid/" -e "s/GEOM/$gshad/" -e "s/BOPT/$bopt/" \
    -e "s/FLDS/$flds/" -e "s/NZETP/$nzetp/" \
    -e "s/VFLD/$vfld/" -e "s/ARIP/$arip/" -e "s/MRIP/$mrip/" < $INT/skel_SHAD.ctl > $gshad.ctl
  ## Step 2.5-3c Results geometry
    rm -f $gres.ctl
    sed -e "s/EQDSK/$eqid/" -e "s/GEOM/$gres/" -e "s/BOPT/$bopt/" \
    -e "s/FLDS/$flds/" -e "s/NZETP/$nzetp/" \
    -e "s/VFLD/$vfld/" -e "s/ARIP/$arip/" -e "s/MRIP/$mrip/" < $INT/skel_SHAD.ctl > $gres.ctl
  ## Step 2.5-4 HDS generation
    rm -f "$gshad"_hds.ctl
    sed -e "s/GEOM/$gshad/" -e "s/NBIN/$nbin/" < $INT/skel_HDS.ctl > "$gshad"_hds.ctl
  ## Step 2.5-5 Power deposition
    rm -f $temp "$gres"_pow.ctl
    cat $INT/skel_POWstart.ctl $INT/skel_POWtermplane.ctl $INT/skel_POWedgprof.ctl $INT/skel_POWodes.ctl > $temp
    sed -e "s/GEOMH/GEOM_hds/" -e "s/GEOM/$gshad/" -e "s/RES/$gres/" -e "s/CAL/$cal/" \
    -e "s/POWE/$powe/" -e "s/DECL/$decl/" < $temp > "$gres"_pow.ctl
    touch ./".STOP-"$runid"-ctl"
  fi
fi
if [ ${swa[1]} == 0 ] ; then
  echo "Now inspect and possibly edit .ctl files before proceeding"
  echo "***rsmiter completed successfully***" 
  exit 2
fi
#
## Step 3. geoq processing
## Step 3b.  Shadow geometry
if [ ${swa[1]} == 2 ] || [ ! -e "$gshad"_geoq.out  ] ; then
  if [ "$eqfile" !=  "$eqid.eqdsk" ] ; then ln -sf $eqfile $eqid.eqdsk;fi
  rm -f "$gshad"_geoq.out
  $INT/cmdwrap $tempdir geoq $gshad
  if [ $? -ne 0 ] ; then exit 1 ;fi
  #
  ## field diagnostics
  ln -sf $HS/Extras/$gwall.txt
  $HS/Extras/geoqgnu $gshad &> /dev/null
fi
## Step 3c. Results geometry
if [ ${swa[1]} == 2 ] || [ ! -e "$gres"_geoq.out ] ; then
  rm -f "$gres"_geoq.out
  $INT/cmdwrap $tempdir geoq $gres
  if [ $? -ne 0 ] ; then exit 1 ;fi
fi
#
## Step 4. HDS generation
if [ ${swa[1]} == 2 ] || [ ! -e "$gshad"_hds_hdsgen.out ] ; then
  rm -f "$gshad"_hds_hdsgen.out
  $INT/cmdwrap $tempdir hdsgen "$gshad"_hds
  if [ $? -ne 0 ] ; then exit 1 ;fi
fi
#
## Step 5. Power deposition results
if [ ${swa[1]} == 2 ] || [ ! -e "$gres"_pow_powcal.out ] ; then
  rm -f "$gres"_pow_powcal.out
  $INT/cmdwrap $tempdir powcal "$gres"_pow
  if [ $? -ne 0 ] ; then exit 1 ;fi
fi
echo "End script" $(date) >> $runlog
echo ""
echo "***rsmiter completed successfully***" 
