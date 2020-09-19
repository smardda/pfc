#! /bin/bash --norc
#script to guide user to, through and beyond SMARDDA-PFC execution
##check sanity of installation
if [[ -n "$SMITER_DIR" ]] ; then HS=$SMITER_DIR; else HS=$HOME/smardda-pfc; fi
if [ ! -d $HS ] ; then
  echo "Neither directory \$SMITER_DIR=$SMITER_DIR nor $HOME/smardda-pfc exists - installation corrupt ?"
  echo "Script quitting"; exit
fi
#set defaults
export HS
INT=$HS/int
#INT=~/dial #fix for testing
# --version and --help handling then quit, or set --tmp
$INT/ismiterhelp.bash $*
retval=$?
##directories
case $retval in
  0) tempdir=./ismtemp
  if [ ! -d $tempdir ] ; then rm -f $tempdir;mkdir $tempdir;fi #preserve existing directory
  ;;
  1) exit ;;
  2) tempdir=/tmp/ism"$LOGNAME"$$.d # uses scratch
    rm -rf $tempdir;mkdir $tempdir ;;
  3) tempdir=/tmp/ism"$LOGNAME"$$.d # uses scratch and
    ## pare back to shell built-ins
    POSIXLY_CORRECT=1
    \unset -f help read unset
    \unset POSIXLY_CORRECT
    while \read cmd; do 
        [[ "$cmd" =~ ^([a-z]+): ]] && \unset ${BASH_REMATCH[1]}; 
    done < <( \help -s "*" )
    echo "Reverting to shell built-ins"
    rm -rf $tempdir;mkdir $tempdir ;;
  *) echo "Leaving subshell with return code $retval and quitting"; exit;;
esac
if [ ! -d $tempdir ] ; then
  echo "Unable to create scratch directory $tempdir - quitting"; exit
fi
##file handling
dfile=$tempdir/dfile; dfiled=$tempdir/dfiled; dfilev=$tempdir/dfilev
ifile=$tempdir/ifile; ifiles=$tempdir/ifiles; iform=$tempdir/iform
doutop=$tempdir/doutop; doutg=$tempdir/doutg; doutt=$tempdir/doutt
dout=$tempdir/dout; dout0=$tempdir/dout0; dout1=$tempdir/dout1
dout2=$tempdir/dout2; dout3=$tempdir/dout3; dout4=$tempdir/dout4
doutsw=$tempdir/doutsw; doutrip=$tempdir/doutrip; doutpow=$tempdir/doutpow
filerm=($dfile $dfiled $dfilev $ifile $doutop $doutg $doutt $dout $dout0 $dout1 $dout2 $dout3 $dout4 $doutsw $doutrip $doutpow) #keep $ifiles and $iform
for i in ${filerm[@]};do rm -f $i;done
#try to speed loading
hash paraview
##Useful web bits
SMDOCWEB="http://arxiv.org/abs/1403.6750 http://arxiv.org/abs/1403.7142" #SMARDDA preprints
PVDOCWEB=http://www.paraview.org/paraview-guide #for download
PVDOC=$SMITER_DOC/TheParaViewGuide-v4.3-CC-Edition.pdf #once downloaded
BROWSER=xdg-open # does not always work
BROWSER=firefox
#check terminal window big enough
rm -f $dout;dialog --print-maxsize 2> $dout # have to do it this way - dialog o/ps non-printing nasties
vara=($(cat $dout));rm -f $dout
winheight=${vara[1]//,/}
winwidth=${vara[2]//[^[0-9]]/}
#echo window size $winheight $winwidth
if [ $winheight -lt 30 ] || [ $winwidth -lt 80 ]  ; then
echo "terminal window is too small for smpfcint"
echo "window must be at least" ; echo "80 characters wide" ; echo "30 characters high" ;
echo "please resize and/or" ; echo "adjust font size"
echo "(Edit->Pref..->App..) and" ; echo "start again"
exit
fi
##smardda-pfc parameters
sp3f=in_eqdsk
gwall=wall
bdry=out
cal=local
pow=1.0
decl=0.010
nbin=20
fldspec=cyl
nzetp=1
arip=0;mrip=0
ulr=_R; ulz=_Z; ulrz=_RZ
angle_units_deg=on; length_units_mm=on
powe=1.e+06; decl=0.01
##inline functions
function isnumber {
parse number $1 2> /dev/null
}
function isposnumber {
parse numbergt0 $1 2> /dev/null
}
function isposvector {
parse vectorgt0 $* 2> /dev/null
}
function isvector {
parse vector $* 2> /dev/null
}
#
#main (top) loop, will normally terminate with exit after Quit
declare -a afiles
chkinput=0
while true ; do #main
  ##Checklist at top level
  dialog --backtitle "SMARDDA-PFC in $PWD" \
  --ok-label "Proceed" --extra-button --extra-label "Help" --cancel-label "Quit" \
  --file $INT/topselmenu.txt 2> $doutop
  retval=$?
  case $retval in
    0) ;;
    1) echo "Interactive script ends ; key commands and outputs logged in $tempdir"
    cp $tempdir/cmdwrap.log ./cmdwrap$$.log
    cp ismfiles.sav ismfiles$$ ; cp ismform.sav ismform$$
    echo "Commands also saved in ./cmdwrap$$.log "
    echo "Other settings in files ./ismfiles$$ and ./ismform$$";exit;;
    3) dialog --begin 1 2 --tailboxbg $INT/helptop.txt 16 80 --and-widget \
    --begin 17 2 --yesno "Do you want still more information?" 0 0
    retval1=$?
    case $retval1 in
      0) dialog --ok-label "Web-based help" --extra-button --extra-label "Quit help" --textbox $INT/morehelptop.txt 12 80
      retval2=$?
      case $retval2 in
        0) $BROWSER $SMITER_DOC/index.html $SMITER_DOC/classes.html  $PVDOCWEB &> /dev/null & ;;
        *) ;;
      esac ;;
      1) dialog --begin 1 2 --infobox "More info not wanted" 6 30; sleep 1;;
      *) dialog --begin 1 2 --infobox "Unexpected input gives return code $retval1 - ignored" 6 50; sleep 1 ;;
    esac ;;
    *) echo "Unexpected return code $retval - script abandoned";exit;;
  esac #retval
  alltopsel=$(cat $doutop)
  if [ "$alltopsel" == "" ] ; then
  dialog --yesno "No options specified - do you want to quit?" 0 0
  retval=$?
  case $retval in
    0) echo "Script quitting ";exit;;
    *) ;;
  esac
  fi
  ##defaults for define surfaces
  angs=-10;angf=10;ndiv=12
  ##defaults for assemble geometry
  tfmtype=0; vis=1; vtkf= ; datf=
  declare -a vecsta ; vecsta=(0 0 0) 
  declare -a vecfin ; vecfin=(1000 1000 1000)
  declare -a vecrip ; vecrip=(0 0 0)
  declare -a vecpow ; vecpow=(0 0 0)
  declare -a tfmsca ; tfmsca=(1 1 1) 
  declare -a tfmoff ; tfmoff=(0 0 0)
  declare -a tfmat1 ; tfmat1=(1 0 0)
  declare -a tfmat2 ; tfmat2=(0 1 0)
  declare -a tfmat3 ; tfmat3=(0 0 1)
  declare -a tfmout 
  # loop over top menu entries
  for itop in $alltopsel;do #top
    case $itop in #top
      ##Assemble geometry 
      1t) while true; do #ag
      ##Checklist for assemble geometry
      dialog --backtitle "SMARDDA-PFC in directory $PWD" \
      --ok-label "Proceed" --extra-button --extra-label "Help" --cancel-label "Back to top" \
      --file $INT/geoselmenu.txt 2> $doutg
      retval=$?
      case $retval in
        0) ;;
        3) dialog --begin 1 2 --tailboxbg $INT/helpgeo.txt 12 80 --and-widget \
        --begin 14 2 --yesno "Do you want still more information?" 0 0
        retval1=$?
        case $retval1 in
          0) dialog --textbox $INT/morehelpgeo.txt 21 80;;
          1) dialog --begin 1 2 --infobox "More info not wanted" 6 30; sleep 1;;
          *) dialog --begin 1 2 --infobox "Unexpected input gives return code $retval1 - ignored" 6 50; sleep 1 ;;
        esac ;;
        *) dialog --begin 1 2 --infobox "Returning to top menu" 6 30; sleep 1 ; break 2 ;;
      esac #retval
      allgeosel=$(cat $doutg)
      if [ "$allgeosel" == "" ] ; then
      dialog --yesno "No options specified - finished with assemble geometry?" 0 0
      retval3=$?
      case $retval3 in
        0) dialog --begin 1 2 --infobox "Returning to top menu" 6 30; sleep 1 ; break ;;
        *) ;;
      esac #retval3
      fi
      #loop over secondary menu items
      for ii in $allgeosel;do #sm
      case $ii in
        #Change visualisation setting 
        vg) stop_visualise=on;vis=1
        dialog --msgbox --trim "3-D visualisation by ParaView.  \
        \nThe online manual for ParaView may be downloaded from \n $PVDOCWEB \
        \nor a local version may be available as file \n  $PVDOC" 10 0
        dialog --title "Change visualisation setting:" \
        --checklist "\nUncheck to keep visualisation" 0 0 12 \
        0 "Stop visualisation on select/assemble :" $stop_visualise \
        2> $dout0
        retval=$?
        case $retval in
          0) defselv=$(cat $dout0)
            for iii in $defselv;do
            case $iii in
            V) stop_visualise=off;vis=0;;
            esac
            done ;;
          *) dialog --infobox "Invalid response - ignored" 6 30; sleep 1;;
        esac ;; #retval
        #Change default geometry settings 
        dg) angle_units_deg=on; length_units_mm=on
        dialog --title "Change default assemble geometry settings:" \
        --checklist "\nUncheck to keep default" 0 0 12 \
        R "Angle in radians (default degrees):" $angle_units_deg \
        M "Length in metres (default millimetres):" $length_units_mm \
        2> $dout0
        retval=$?
        case $retval in
          0) defsel=$(cat $dout0)
            for iii in $defsel;do
            case $iii in
            R) angle_units_deg=off;;
            M) length_units_mm=off;;
            esac
            done ;;
          *) dialog --infobox "Invalid response - ignored" 6 30; sleep 1;;
        esac ;; #retval
        #Define surfaces
        0g) afiles=($(ls -L 2> /dev/null  *$ulr*.dat));ndat=${#afiles[@]} 
        if [ $ndat == 0 ] ; then
          dialog --begin 1 2 --msgbox "Copy $ulr.dat file(s) into the current directory $PWD." 6 50; break 
        else
          okinput=0
          until [ $okinput == 1 ] ; do
            $INT/buil.bash rdat $dfiled
            $tempdir/filsel.bash # tempdir=${dfiled%/*}
            retvalo=$?
            case $retvalo in
              0) datr=$(cat $dfiled); datra=($datr);;
              1) break 2 ;; #2
              3) dialog --textbox $INT/help$ulr.txt 12 80 ; break 2 ;; #help
              *) datr=;dialog --infobox "Invalid response - try again" 6 30; sleep 1;break;;
            esac #retvalo
            ndatr=${#datra[@]}
            if [ $ndatr -eq 1 ] ; then
              if [[ $datr == *"$ulr"* ]] ; then
                if [[ $datr == *"$ulrz"* ]] ; then
                  datf=${datr//_RZ/} #file root
                  datrz=$datr.dat;datr=null;datz=null;okinput=1
                else
                  datf=${datr//_R/} #file root
                  datz=${datr//_R/_Z}.dat; datr=$datr.dat ;datrz=null
                  if [ ! -e $datz ] ; then
                     dialog --msgbox "Matching file $datz does not exist " 6 50
                  else
                     okinput=1
                  fi
                fi 
              fi 
            elif [ $ndatr -gt 1 ] ; then 
              dialog --msgbox "Only select one $ulr.dat file" 6 50
            else
              dialog --infobox "No $ulr.dat file selected" 6 50; sleep 1
            fi
          done #okinput
          #build menu
          if [ $angle_units_deg == on ] ; then angu=degree; else angu=radian;fi
          if [ $length_units_mm == on ] ; then lenu=mm; else lenu=metres;fi
#----------------------------------------------------------
          dialog --yes-label "Angle" --no-label "Vector" \
          --yesno "Sweep in angle or along vector?" 0 0
          retvali=$?
          case $retvali in
            0) tfmtyp=rotate;;
            *) tfmtyp=translate;;
          esac
          okinput=0  
          until [ $okinput == 1 ] ; do
            rm -f $doutsw
            case $tfmtyp in
              rotate) dialog --form "Sweep parameters"  0 0 12 \
              "start angle:" 1 1 "$angs" 1 25 25 30 \
              "stop angle:" 2 1 "$angf" 2 25 25 30 \
              "divisions in angle:" 3 1 "$ndiv" 3 25 25 30 \
              2> $doutsw
              retvali=$?
              ## if any entry is empty substitute string 'FILL-IN'
              rm -f $dout4;sed -e 's/^ *$/FILL-IN/' < $doutsw > $dout4
              diff $doutsw $dout4 >/dev/null;difft=$?
              if [ $difft == 1 ] ; then 
                dialog --begin 1 2 --infobox "Missing entry " 6 30; sleep 1
              else
                case $retvali in
                  0) sweepsdirt=($(cat $dout4))
                  sweeps=(${sweepsdirt[@]//[^[:alnum:]+.-]/})
                  #Decho ${sweepsdirt[@]} ${sweeps[@]} > doutsw2
                  if  $(isposnumber ${sweeps[2]}) ; then
                    if  $(isvector ${sweeps[*]}) ; then okinput=1 
                      else  dialog --infobox "Numeric input missing or corrupt" 6 35; sleep 1;fi
                    else  dialog --infobox "Positive number of divisions wanted" 6 40; sleep 1;fi
                  angs=${sweeps[0]}; angf=${sweeps[1]}; ndiv=${sweeps[2]};;
                  1) break ;;
                  *) dialog --infobox "Unexpected return code $retvali" 6 40; sleep 1 ;;
                esac #retvali
              fi
              ;;
              translate) dialog --form "Sweep parameters"  0 0 12 \
              "start in x direction (mm):" 1 1 "${vecsta[0]}" 1 25 25 30 \
              "start in y direction (mm):" 2 1 "${vecsta[1]}" 2 25 25 30 \
              "start in z direction (mm):" 3 1 "${vecsta[2]}" 3 25 25 30 \
              "finish in x direction (mm):" 4 1 "${vecfin[0]}" 4 25 25 30 \
              "finish in y direction (mm):" 5 1 "${vecfin[1]}" 5 25 25 30 \
              "finish in z direction (mm):" 6 1 "${vecfin[2]}" 6 25 25 30 \
              "divisions in sweep:" 7 1 "$ndiv" 7 25 25 70 \
              2> $doutsw
              retvali=$?
              ## if any entry is empty substitute string 'FILL-IN'
              rm -f $dout4;sed -e 's/^ *$/FILL-IN/' < $doutsw > $dout4
              diff $doutsw $dout4 >/dev/null;difft=$?
              if [ $difft == 1 ] ; then 
                dialog --begin 1 2 --infobox "Missing entry " 6 30; sleep 1
              else
                case $retvali in
                  0) tfmdirt=($(cat $dout4))
                    tfmclean=(${tfmdirt[@]//[^[:alnum:]+.-]/})
                    vecsta=(${tfmclean[0]} ${tfmclean[1]} ${tfmclean[2]})
                    vecfin=(${tfmclean[3]} ${tfmclean[4]} ${tfmclean[5]})
                    ndiv=${tfmclean[6]}
                    if  $(isposnumber $ndiv) ; then
                      if  $(isvector ${vecsta[*]}) && $(isvector ${vecfin[*]}); then okinput=1 
                      else  dialog --infobox "Numeric input missing or corrupt" 6 35; sleep 1;fi
                    else  dialog --infobox "Positive number of divisions wanted" 6 40; sleep 1;fi;;
                  1) break ;;
                  *) dialog --infobox "Unexpected return code $retvali" 6 50; sleep 1 ;;
                esac #retvali
              fi
              ;;
              *) break;;
            esac
          done #okinput
          rm -f $dout2
          tfmst1=${vecsta[@]};tfmfi1=${vecfin[@]}
          sed -e "s/ANGU/$angu/" -e "s/LENU/$lenu/" -e "s/TFMTYP/$tfmtyp/" \
          -e "s/DATRZ/$datrz/" -e "s/DATR/$datr/" -e "s/DATZ/$datz/"  \
          -e "s/TFMSTA/$tfmst1/" -e "s/TFMFIN/$tfmfi1/" \
          -e "s/ANGS/$angs/" -e "s/ANGF/$angf/" -e "s/NDIV/$ndiv/"  < $INT/skel$ulr.ctl > $dout2
          retvals=1
          until [ $retvals == 0 ] ; do
            dialog --ok-label "Proceed" --extra-button --extra-label "Web-based help" --cancel-label "Ignore edits" \
            --title "Edit ctl file at own risk" --editbox $dout2 24 100 2> $dout3
            retvals=$?
              case $retvals in
                0) rm -f $datf.ctl; mv -f $dout3 $datf.ctl;;
                1) break;;
                3) $BROWSER $SMITER_DOC/classes.html $SMITER_DOC/index.html &> /dev/null &;;
              esac
          done #retvals
          $INT/cmdwrap $tempdir datvtk -c $datf  2>&1 | dialog --programbox "Generating surface vtk file" 20 80
          if [ -e $tempdir/cmdwrap_datvtk.err ] ; then dialog --msgbox "Error - try again" 6 30;break;else
          if [ $vis == 1 ] ; then paraview --data="$datf".vtk ; fi; fi
        fi ;;
        #Use dat files
        1g) afiles=($(ls -L 2> /dev/null  *.dat))
        bfiles=
        for i in ${afiles[@]}; do
           if [[ $i == *"$ulr"* ]] || [[ $i == *"$ulz"* ]]; then continue;else bfiles=(${bfiles[@]} $i);fi
        done
        ndat=${#bfiles[@]} 
        if [ $ndat == 0 ] ; then
          dialog --begin 1 2 --msgbox "Copy dat file(s) into the current directory $PWD." 6 60; sleep 1;break 
        else
          okinput=0
          until [ $okinput == 1 ] ; do
            $INT/buil.bash Xdat $dfiled
            $tempdir/filsel.bash # tempdir=${dfiled%/*}
            retvalo=$?
            case $retvalo in
              0) datf=$(cat $dfiled);;
              1) break 2 ;; #2
              *) datf=;dialog --msgbox "Invalid response - try again" 6 30;break;;
            esac #retvalo
            if [ ${#datf} -gt 0 ] ; then okinput=1
              else dialog --infobox "No file selected " 6 30; sleep 1;fi
          done #okinput
          $INT/cmdwrap $tempdir datvtk -xxm $datf  2>&1 | dialog --programbox "Converting dat file(s)" 20 80
          if [ -e $tempdir/cmdwrap_datvtk.err ] ; then dialog --msgbox "Error - try again" 6 30;break;fi
        fi ;;
        #Select geometry file(s) 
        2g) afiles=($(ls -L 2> /dev/null  *.vtk)) ; nvtk=${#afiles[@]}
        if [ $nvtk == 0 ] ; then
          dialog --begin 1 2 --msgbox "Copy vtk file(s) into the current directory $PWD" 6 60; sleep 1;break
        else
          okinput=0
          until [ $okinput == 1 ] ; do
            $INT/buil.bash vtk $dfilev
            $tempdir/filsel.bash  # tempdir=${dfilev%/*}
            retvalo=$?
            case $retvalo in
              0) vtkf=$(cat $dfilev);vtkfa=($vtkf);;
              1) break 2;; #2
              *) vtkf="";dialog --infobox "Invalid response - try again" 6 30; sleep 1;break;;
            esac #retvalo
            cvtkf=${#vtkf}
            wvtkf=${#vtkfa[@]}
              #D svtkf=${vtkf// /};wvtkf=$((${#vtkf}-${#svtkf})) #(counts number of blanks)
              #D dialog --infobox "$vtkf,$cvtkf, $wvtkf" 6 30; sleep 2
            if [ $cvtkf -gt 0 ] && [ $cvtkf -lt 133 ] ; then 
              if [ $wvtkf -gt 0 ] && [ $wvtkf -lt 21 ] ; then okinput=1; fi
            elif [ $cvtkf -gt 132 ] || [ $wvtkf -gt 20 ] ; then
              dialog --msgbox "Too many files at once - no more than 20 files \
               \n or 132 characters" 6 50; break
            else 
              dialog --infobox "No files selected" 6 30;sleep 1;break
            fi
          done #okinput
        vtkfil=$(sed -e "s/ *$/.vtk;/" -e "s/  */.vtk;/g" < $dfilev)
        if [ $vis == 1 ] && [ $okinput == 1 ] ; then paraview --data="$vtkfil" ;fi
        fi ;; #nvtk
        #Define transform
        3g) if [ "$vtkf" == "" ] ; then
        dialog --msgbox "vtk file must be selected before transform defined - try again" 6 50; break;fi
        fileroot=${vtkf// /+}
        tfname=
        until [ "$tfname" != "" ] ; do
          dialog --inputbox "Name of transform? (Alphanumerics and underscore only)" 6 50 "tfm_1" 2> $dout 
          dirty=$(cat $dout)
          tfname=${dirty//[^[:alnum:]_-]/}
        done #tfname
        tfmtype=
        until [ "$tfmtype" != "" ] ; do
          dialog --begin 1 2 --title "Transform $tfname" \
          --file $INT/tfmenu.txt 2> $doutt
          retval=$?
          case $retval in
            0) tfmtype=$(cat $doutt) ;;
            *) break;;
          esac
        done #tfmtype
        tfmfil=$fileroot"_"$tfname
        rm -f $tfmfil.ctlin; echo 'STA vtktfm' > $tfmfil.ctlin;  echo "DEF $defsel" >> $tfmfil.ctlin;\
        echo "FIL $vtkf" >> $tfmfil.ctlin; echo "TFM $tfmtype" >> $tfmfil.ctlin
        okinput=0
        until [ $okinput == 1 ] ; do
          case $tfmtype in
            0) break;;
            1) dialog --form "Transform parameters"  0 0 12 \
            "scale in x direction:" 1 1 "${tfmsca[0]}" 1 25 25 30 \
            "scale in y direction:" 2 1 "${tfmsca[1]}" 2 25 25 30 \
            "scale in z direction:" 3 1 "${tfmsca[2]}" 3 25 25 30 \
            2> $dout
            retvali=$?
            case $retvali in
              0) tfmsca=($(cat $dout))
              if  $(isposvector ${tfmsca[*]}) ; then okinput=1 
              echo "SCA ${tfmsca[*]}" >> $tfmfil.ctlin
              else  dialog --infobox "Positive parameters wanted" 6 30; sleep 1;fi;;
              1) break ;;
              *) dialog --infobox "Unexpected return code $retvali" 6 50; sleep 1 ;;
            esac #retvali
            ;;
            2) dialog --form "Transform parameters"  0 0 12 \
            "offset in x direction:" 1 1 "${tfmoff[0]}" 1 25 25 30 \
            "offset in y direction:" 2 1 "${tfmoff[1]}" 2 25 25 30 \
            "offset in z direction:" 3 1 "${tfmoff[2]}" 3 25 25 30 \
            "scale in x direction:" 4 1 "${tfmsca[0]}" 4 25 25 30 \
            "scale in y direction:" 5 1 "${tfmsca[1]}" 5 25 25 30 \
            "scale in z direction:" 6 1 "${tfmsca[2]}" 6 25 25 30 \
            2> $dout
            retvali=$?
            case $retvali in
              0) tfmout=($(cat $dout))
                tfmoff=(${tfmout[0]} ${tfmout[1]} ${tfmout[2]})
                tfmsca=(${tfmout[3]} ${tfmout[4]} ${tfmout[5]})
                if  $(isposvector ${tfmsca[*]}) && $(isvector ${tfmoff[*]}); then okinput=1 
                echo "SCA ${tfmsca[*]}" >> $tfmfil.ctlin
                else  dialog --infobox "Positive scales and numeric offsets wanted" 6 50; sleep 1;fi;;
              1) break ;;
              *) dialog --infobox "Unexpected return code $retvali" 6 50; sleep 1 ;;
            esac #retvali
            ;;
            6) dialog --form "Transform parameter"  0 0 12 \
            "Tilt angle in poloidal direction:" 1 1 "${tfmoff[0]}" 1 40 40 45 \
            2> $dout
            retvali=$?
            case $retvali in
              0) tfmoff=($(cat $dout) 0 0)
              if  $(isvector ${tfmoff[*]}) ; then okinput=1 
              else  dialog --infobox "Numeric input missing or corrupt" 6 35; sleep 1;fi;;
              1) break ;;
              *) dialog --infobox "Unexpected return code $retvali" 6 50; sleep 1 ;;
            esac #retvali
            ;;
            7) dialog --form "Transform parameter"  0 0 12 \
            "Tilt angle in toroidal direction:" 1 1 "${tfmoff[0]}" 1 40 40 45 \
            2> $dout
            retvali=$?
            case $retvali in
              0) tfmoff=($(cat $dout) 0 0)
              if  $(isvector ${tfmoff[*]}) ; then okinput=1 
              else  dialog --infobox "Numeric input missing or corrupt" 6 35; sleep 1;fi;;
              1) break ;;
              *) dialog --infobox "Unexpected return code $retvali" 6 50; sleep 1 ;;
            esac #retvali
            ;;
            12) dialog --form "Transform parameter"  0 0 12 \
            "Rotation angle in toroidal direction:" 1 1 "${tfmoff[2]}" 1 40 40 45 \
            2> $dout
            retvali=$?
            case $retvali in
              0) tfmoff=(0 0 $(cat $dout))
              if  $(isvector ${tfmoff[*]}) ; then okinput=1 
              else  dialog --infobox "Numeric input missing or corrupt" 6 35; sleep 1;fi;;
              1) break ;;
              *) dialog --infobox "Unexpected return code $retvali" 6 50; sleep 1 ;;
            esac #retvali
            ;;
            22) dialog --form "Transform parameter"  0 0 12 \
            "Displacement in poloidal direction:" 1 1 "${tfmoff[2]}" 1 40 40 45 \
            2> $dout
            retvali=$?
            case $retvali in
              0) tfmoff=(0 0 $(cat $dout))
              if  $(isvector ${tfmoff[*]}) ; then okinput=1 
              else  dialog --infobox "Numeric input missing or corrupt" 6 35; sleep 1;fi;;
              1) break ;;
              *) dialog --infobox "Unexpected return code $retvali" 6 50; sleep 1 ;;
            esac #retvali
            ;;
            42) tfmoff[1]=6200;tfmoff[2]=0
            if [ $length_units_mm == off ] ; then tfmoff[1]="6.200";fi
            #dialog --form "Negative implies displacement towards centre"  0 0 12 \
            dialog --form "Transform parameters"  0 0 12 \
            "Displacement in minor radius :" 1 1 "${tfmoff[0]}" 1 40 40 45 \
            "Negative implies displacement towards centre" 2 1 "" 2 50 0 1 \
            "Central position R coordinate :" 3 1 "${tfmoff[1]}" 3 40 40 45 \
            "Central position Z coordinate :" 4 1 "${tfmoff[2]}" 4 40 40 45 \
            2> $dout
            retvali=$?
            case $retvali in
              0) tfmoff=($(cat $dout))
              if  $(isvector ${tfmoff[*]}) ; then okinput=1 
              else  dialog --infobox "Numeric input missing or corrupt" 6 35; sleep 1;fi;;
              1) break ;;
              *) dialog --infobox "Unexpected return code $retvali" 6 50; sleep 1 ;;
            esac #retvali
            ;;
            *) break;;
          esac
        #write to transform file
        done #okinput
        if [ $okinput == 1 ] ; then echo "OFS ${tfmoff[*]}" >> $tfmfil.ctlin;fi
        ;;
        #Apply transform
        4g)  if [ "$vtkf" == "" ] ; then
          dialog --infobox "Select file to transform" 6 30; sleep 1;break
        fi
        $INT/cmdwrap $tempdir ctlgen -v $tfmfil  2>&1 | dialog --programbox "Preparing .ctl file" 12 80
        if [ -e $tempdir/cmdwrap_ctlgen.err ] ; then dialog --msgbox "Error - try again" 6 30;break;fi
        tfmctl=$tfmfil.ctl
        #rm -f $dout[23]; sed -e 's/.*=/\L&/' -e 's/^ *\&.*/\L&/' < $tfmctl >$dout2
        rm -f $dout[23]; sed -e 's/.*=/\L&/' -e 's/^ *\&.*/\L&/' -e 's/  */ /g' < $tfmctl >$dout2
        retval=1
        until [ $retval == 0 ] ; do
          dialog --ok-label "Proceed" --extra-button --extra-label "Web-based help" --cancel-label "Ignore edits" \
          --title "Edit ctl file at own risk" --editbox $dout2 24 100 2> $dout3
          retval=$?
          case $retval in
            0) rm -f $tfmctl; mv -f $dout3 $tfmctl;;
            1) break;;
            3) $BROWSER $SMITER_DOC/classes.html $SMITER_DOC/index.html &> /dev/null &;;
          esac
        done #retval
        $INT/cmdwrap $tempdir vtktfm $tfmfil  2>&1 | dialog --programbox "Transforming vtk file(s)" 20 80
        if [ -e $tempdir/cmdwrap_vtktfm.err ] ; then dialog --infobox "Error - try again" 6 30;sleep 1;break;fi
        vtkfilt=$tfmfil.vtk
        if [ $vis == 1 ] ; then paraview --data="$vtkfil;$vtkfilt";fi
        #suppress combine geometry after transform
        break ;; 
        #Combine geometry
        5g) if [ "$vtkf" == "" ] ; then
          dialog --infobox "Select files to combine" 6 30; sleep 1;break;fi
        geoname=
        until [ "$geoname" != "" ] ; do
          dialog --inputbox "Name of combined geometry? (Alphanumerics and underscore only)" 6 50 "your_choice" 2> $dout 
          dirty=$(cat $dout)
          geoname=${dirty//[^[:alnum:]_-]/}
        done #geoname
        rm -f $geoname.ctlin; echo 'STA vtktfm' > $geoname.ctlin;  echo "DEF $defsel" >> $geoname.ctlin;\
        echo "FIL $vtkf" >> $geoname.ctlin; echo "TFM 0" >> $geoname.ctlin
        combctl=$geoname.ctl;rm -f $combctl
        $INT/cmdwrap $tempdir ctlgen -v $geoname  2>&1 | dialog --programbox "Preparing .ctl file" 12 80
        if [ -e $tempdir/cmdwrap_ctlgen.err ] ; then dialog --msgbox "Error - try again" 6 30;break;fi
        #rm -f $dout[23]; sed -e 's/.*=/\L&/' -e 's/^ *\&.*/\L&/' < $combctl >$dout2
        rm -f $dout[23]; sed -e 's/.*=/\L&/' -e 's/^ *\&.*/\L&/' -e 's/  */ /g' < $combctl >$dout2
        retval=1
        until [ $retval == 0 ] ; do
          dialog --ok-label "Proceed" --extra-button --extra-label "Web-based help" --cancel-label "Ignore edits" \
          --title "Edit ctl file at own risk" --editbox $dout2 24 100 2> $dout3
          retval=$?
          case $retval in
            0) rm -f $combctl; mv -f $dout3 $combctl;;
            1) break;;
            3) $BROWSER $SMITER_DOC/classes.html $SMITER_DOC/index.html &> /dev/null &;;
            *) dialog --infobox "Unexpected code $retval - try again (1t)" 6 50;sleep 1;break;;
          esac
        done #retval
        $INT/cmdwrap $tempdir vtktfm $geoname  2>&1 | dialog --programbox "Combining vtk file(s)" 20 80
        if [ -e $tempdir/cmdwrap_vtktfm.err ] ; then dialog --infobox "Error - try again" 6 30;sleep 1;break;fi
        vtkfil=$geoname.vtk
        if [ $vis == 1 ] ; then paraview --data="$vtkfil";fi
        ;;
        *) dialog --msgbox "Unrecognised selection from a lower menu" 6 50;;
      esac
      done #sm  
      done;;  #ag
      ##Select geometry and equil files
      2t) ##see whether have already selected files
      if [ ! -e $ifiles ] && [ -e ismfiles.sav ] ; then cp ismfiles.sav $ifiles; fi
      if [ -e $ifiles ] ; then
        dialog --yesno "Use saved input data file names?" 0 0 
        retval=$?
        case $retval in
          0) for iii in $(cat $ifiles); do
            key=${iii%%=*}
            arg=${iii##*=}
            case "$key" in
              vtk1f) vtk1f=$arg;;
              vtk1) vtk1=$arg;;
              vtk2f) vtk2f=$arg;;
              vtk2) vtk2=$arg;;
              eqdskf) eqdskf=$arg;;
              eqdsk) eqdsk=$arg;;
              sp3f) sp3f=$arg;;
              sp3) sp3=$arg;;
              *) dialog --infobox "key $key in file $ifiles ignored " 6 50;sleep 1;;
            esac 
          done;;
          1) dialog --infobox "Ignoring saved file names" 6 30; sleep 1 ;;
          *) dialog --infobox "Unexpected code $retval - try again (2t)" 6 50;sleep 1;break;;
        esac
      else
        retval=1
      fi
      if [ $retval == 1 ]; then #start of defining input data filess
        ##checklist for input data files
        dialog --file $INT/chkinp.txt 2> $dfiled
        afiles=($(ls -L 2> /dev/null  *.vtk)) ; nvtk=${#afiles[@]}
        afiles=($(ls -L 2> /dev/null  *.sp3)) ; nsp3=${#afiles[@]}
        afiles=($(ls -L 2> /dev/null  *.eqdsk)) ; neqdsk=${#afiles[@]}
        lis1=$(cat $dfiled)
        if [ "$lis1" == "" ] || [[ "$lis1" != "1"* ]]  || [ $nvtk == 0 ] ; then
          dialog --msgbox  "Geometry must defined either using dat (MSC Nastran ex CATIA) \
          \nor Legacy vtk format files (ANSYS APDL script) for surfaces \
          \nEFIT++ and many other codes generate eqdsk files \
          \nSMARDDA-PFC/magtfm can help generate vacuum field file with ripple \
          \nYou must copy at least one dat or vtk file to this directory" 8 80;break
        fi
        vtk1=on; vtk2=off; eqdsk=off; sp3=off
        for iii in $lis1; do
          case $iii in
            2) vtk2=on;;
            3) eqdsk=on;;
            4) sp3=on;;
          esac
        done
        ##immediate breaks
        if [ $nvtk == 1 ] && [ $vtk2 == on ] ; then
          dialog --msgbox "Only one vtk file found where two expected\
          \nCopy file into current directory" 6 50;break
        fi
        if [ $neqdsk == 0 ] && [ $eqdsk == on ]  ; then
          dialog --msgbox "No eqdsk file found where one expected\
          \nCopy file into current directory" 6 50;break
        fi
        if [ $nsp3 == 0 ] && [ $sp3 == on ]  ; then
          dialog --msgbox "No sp3 file found where one expected\
          \nCopy file into current directory" 6 50;break
        fi
        # vtk file(s)
        if [ $nvtk == 1 ]  ; then
          vtk1f=$(ls -L 2> /dev/null  *.vtk)
          vtk2f=$vtk1f
        else
          $INT/menu.bash vtk.1 $dfile
          $tempdir/filsel.bash # tempdir=${dfile%/*}
          vtk1f=$(cat $dfile)'.vtk'
          if [ $vtk2 == on ] ; then
            $INT/menu.bash vtk.2 $dfile
            $tempdir/filsel.bash # tempdir=${dfile%/*}
            vtk2f=$(cat $dfile)'.vtk'
          fi
        fi
        # eqdsk file
        if [ $neqdsk == 1 ]  ; then
          eqdskf=$(ls -L 2> /dev/null  *.eqdsk)
        elif [ $eqdsk == on ] ; then
          $INT/menu.bash eqdsk $dfile
          $tempdir/filsel.bash # tempdir=${dfile%/*}
          eqdskf=$(cat $dfile)'.eqdsk'
        fi
        ## need to do something here about analytic equil ??
        # sp3 file
        if [ $nsp3 == 1 ]  ; then
          sp3f=$(ls -L 2> /dev/null  *.sp3)
        elif [ $sp3 == on ] ; then
          $INT/menu.bash sp3 $dfile
          $tempdir/filsel.bash # tempdir=${dfile%/*}
          sp3f=$(cat $dfile)'.sp3'
        fi
        rm -f $ifiles ismfiles.sav
        printf "vtk1f=$vtk1f  \
        \nvtk2f=$vtk2f  \
        \neqdskf=$eqdskf  \
        \nsp3f=$sp3f" > $ifiles
        cp $ifiles ismfiles.sav
      fi ;;
      ##Set run type
      3t) #Checklist for ctl files
      ##anti-defaults - checked by default
      flds=1;fldspec=flux;limiter=0;cal=global
      eqdodgy=1;powexp=0;bdry=in
      ##defaults -  unchecked by default
      nzetp=1;analrip=0
      dialog --file $INT/chkctl.txt 2> $dfiled
      retvalr=$?
      case $retvalr in
        0) allctla=($(cat $dfiled))
        for ctl in ${allctla[@]}; do 
          case $ctl in
            0c) flds=3;fldspec=cyl;;
            1c) limiter=1;;
            2c) nzetp=18;;
            3c) cal=local;;
            4c) eqdodgy=0;;
            5c) powexp=1;;
            6c) bdry=out;;
            7c) analrip=1;;
          esac
        done
        if [ $eqdodgy == 1 ] ; then
          dialog --textbox $INT/helpeqdsk.txt 12 80;sleep 1
        fi
        if [ $powexp == 0 ] ; then
          dialog --textbox $INT/helpow.txt 12 80;sleep 1
        elif [ $powexp == 1 ] ; then
          #menu to specify powe and decl
          okin == 0
          until [ $okin == 1 ] ; do
            vecpow[0]=$powe;vecpow[1]=$decl;vecpow[2]=1
            dialog --form "Simple exponential profile"  0 0 12 \
            "Total power in W :" 1 1 "${vecpow[0]}" 1 40 40 45 \
            "Layer width in m :" 2 1 "${vecpow[1]}" 2 40 40 45 \
            2> $doutpow
            retvali=$?
            case $retvali in
              0) vecpow=($(cat $doutpow) 1)
              if  $(isposvector ${vecpow[*]}) ; then okin=1 
              powe=${vecpow[0]};decl=${vecpow[1]}
              else  dialog --infobox "Numeric input missing or corrupt" 6 35; sleep 1;fi;;
              1) break ;;
              *) dialog --infobox "Unexpected return code $retvali" 6 50; sleep 1 ;;
            esac #retvali
          done
        fi
        if [ $analrip == 1 ] ; then
          #menu to specify arip and mrip
          okin == 0
          until [ $okin == 1 ] ; do
            vecrip[0]=$arip;vecrip[1]=$mrip;vecrip[2]=0
            dialog --form "Amplitude is in T"  0 0 12 \
            "Amplitude of ripple :" 1 1 "${vecrip[0]}" 1 40 40 45 \
            "Toroidal mode number of ripple :" 2 1 "${vecrip[1]}" 2 40 40 45 \
            2> $doutrip
            retvali=$?
            case $retvali in
              0) vecrip=($(cat $doutrip) 0)
              if  $(isvector ${vecrip[*]}) ; then okin=1 
              arip=${vecrip[0]};mrip=${vecrip[1]}
              else  dialog --infobox "Numeric input missing or corrupt" 6 35; sleep 1;fi;;
              1) break ;;
              *) dialog --infobox "Unexpected return code $retvali" 6 50; sleep 1 ;;
            esac #retvali
          done
          if [ ! $mrip == 0 ] ; then
            sp3f="in_eqdsk"
            flds=2;fldspec=cyl_nonsymm
            dialog --msgbox "analytic ripple definition used instead of data file" 6 80
          fi
        fi;;
        1) okinput=1 ;;
        3) $BROWSER $SMITER_DOC/classes.html $SMITER_DOC/index.html &> /dev/null &;;
        *) dialog --infobox "Unexpected code $retvalr - try again (3)" 6 50;sleep 1 ;;
      esac #retvalr
      eqid=${eqdskf%%.eqdsk}
      gwall=wall
      ;;
      ##Change default parameters 
      4t) ##Form for rsmardda-pfc control parameters
      ##see whether have already selected form parameters
      if [ ! -e $iform ] && [ -e ismform.sav ] ; then cp ismform.sav $iform; fi
      if [ -e $iform ] ; then
        dialog --yesno "Use saved input form parameters?" 0 0 
        retval=$?
        case $retval in
          0) for iii in $(cat $iform); do
            key=${iii%%=*}
            arg=${iii##*=}
              case "$key" in
              eqid) eqid=$arg;;
              fldspec) fldspec=$arg;;
              bdry) bdry=$arg;;
              cal) cal=$arg;;
              powe) powe=$arg;;
              decl) decl=$arg;;
              nzetp) nzetp=$arg;;
              arip) arip=$arg;;
              mrip) mrip=$arg;;
              nbin) nbin=$arg;;
              *) dialog --infobox "key $key in file $iform ignored " 6 50;sleep 1;;
            esac 
          done;;
          1) dialog --infobox "Ignoring saved form parameters" 6 40; sleep 1 ;;
          *) dialog --infobox "Unexpected code $retval - try again (2t)" 6 50;sleep 1;break;;
        esac
      else
        retval=1
      fi
      if [ $retval == 0 ]  || [ $retval == 1 ] ; then #start of (re)defining form paraemters
        chkinput=0
        until [ $chkinput == 1 ] ; do
          dialog --ok-label "Proceed" --extra-button --extra-label "Web-based help" --cancel-label "Ignore changes" \
          --title "Calculation parameters:" \
          --form "\nChange parameter values at own risk" 0 0 15 \
          "vacuum field definition:" 1 1 "$sp3f" 1 25 35 70\
          "equilibrium field:" 2 1 "$eqdskf" 2 25 35 70\
          "short id for equil:" 3 1 "$eqid" 3 25 25 30 \
          "shadowed geometry:" 4 1 "$vtk1f" 4 25 35 70\
          "shadowing geometry:" 5 1 "$vtk2f" 5 25 35 70\
          "mapping type:" 6 1 "$fldspec" 6 25 25 30 \
          "Q from boundary:" 7 1 "$bdry" 7 25 25 30 \
          "calculation type:" 8 1 "$cal" 8 25 25 30 \
          "total power in W:" 9 1 "$powe" 9 25 25 30 \
          "layer width in m:" 10 1 "$decl" 10 25 25 30 \
          "toroidal periodicity n=:" 11 1 "$nzetp" 11 25 25 30 \
          "ripple amplitude in T:" 12 1 "$arip" 12 25 25 30 \
          "ripple mode number n=:" 13 1 "$arip" 13 25 25 30 \
          "HDS objects grouped in:" 14 1 "$nbin" 14 25 25 30 \
          "Try increasing x 10 if hdsgen fails" 15 1 "" 15 40 0 1 \
          2> $dout2
          retvalr=$?
          case $retvalr in
            0) ## if any entry is empty substitute string 'FILL_IN'
            rm -f $dout3;sed -e 's/^ *$/FILL_IN/' < $dout2 > $dout3
            diff $dout2 $dout3 >/dev/null;difft=$?
            if [ $difft == 1 ] ; then 
            dialog --begin 1 2 --infobox "Missing entry " 6 30; sleep 1;fi
            TFORM=($(cat $dout3))
            sp3f=${TFORM[0]//[^[:alnum:]_.-/]/}
            eqdskf=${TFORM[1]//[^[:alnum:]_.-/]/}
            eqid=${TFORM[2]//[^[:alnum:]_.-/]/}
            vtk1f=${TFORM[3]//[^[:alnum:]_.-/]/}
            vtk2f=${TFORM[4]//[^[:alnum:]_.-/]/}
            fldspec=${TFORM[5]//[^[:alnum:]_.-/]/}
            bdry=${TFORM[6]//[^[:alnum:]_.-/]/}
            cal=${TFORM[7]//[^[:alnum:]_.-/]/}
            powe=${TFORM[8]//[^[:alnum:].+-]/}
            decl=${TFORM[9]//[^[:alnum:].+-]/}
            nzetp=${TFORM[10]//[^[:alnum:].+]/} #strip minus
            arip=${TFORM[11]//[^[:alnum:].+-]/} 
            mrip=${TFORM[12]//[^[:alnum:].+]/} #strip minus
            nbin=${TFORM[13]//[^[:alnum:].+]/} #strip minus
            oknum=0;if $(isposnumber $powe) &&  $(isposnumber $decl) && $(isnumber $nzetp) && \
               $(isnumber $arip) && $(isnumber $mrip) && $(isposnumber $nbin)  ; then oknum=1 ;fi
            if [ $oknum == 0 ] ; then dialog --begin 1 2 --infobox "Numeric input missing or corrupt" 6 35; sleep 1;fi
            #if  $(isposnumber $powe) && $(isposnumber $decl) ; then oknum=1
            flds=0; if [ $fldspec == flux ] ; then flds=1
              elif [ $fldspec == cyl ] ; then flds=3
              elif [ $fldspec == cyl_nonsymm ] ; then flds=2
              else dialog --begin 1 2 --msgbox "$fldspec is invalid mapping" 6 30; sleep 1;fi
            if [ $difft == 0 ] && [ $oknum == 1 ] && [ $flds -gt 0 ] ; then chkinput=1 ;fi
            ;;
            1) chkinput=1 ;;
            3) $BROWSER $SMITER_DOC/classes.html $SMITER_DOC/index.html &> /dev/null &;;
            *) dialog --infobox "Unexpected code $retvalr - try again (4t)" 6 50;sleep 1 ;;
          esac #retvalr
        done #chkinput
        rm -f $iform ismform.sav
        printf "eqid=$eqid \
        \nfldspec=$fldspec \
        \nbdry=$bdry \
        \ncal=$cal \
        \npowe=$powe \
        \ndecl=$decl \
        \nnzetp=$nzetp \
        \narip=$arip \
        \nmrip=$mrip \
        \nnbin=$nbin" > $iform
        cp $iform ismform.sav
      fi ;;
      ##Run SMARDDA-PFC
      5t) ##Checklist for rsmardda-pfc
      if [ $chkinput == 0 ] ; then dialog --msgbox "Checklist was not successfully completed \
      \nSelect (at least) item t4 on Main Menu" 6 50;break ;fi
      #code for bopt setting consistent with parameters set above
      if [ $limiter == 0 ] ; then
        bdryl=x
        if [ $bdry == in ] ; then bopt=7;else bopt=10; fi
      else
        bdryl=$bdry;bopt=2
      fi
      if [ $sp3f == in_eqdsk ] ; then vfld=null;else vfld=$sp3f;fi
      eqfile=$eqdskf; gres=${vtk1f%%.vtk}; gshad=${vtk2f%%.vtk} 
      smdone=9 
      until [ $smdone == 0 ] ; do
        dialog --item-help --checklist "SMARDDA-PFC execution : The computation will " 0 0 6 \
        1 "Generate ctl files" on "option to edit ctl (outside SMARDDA-PFC) before running " \
        2 "Run SMARDDA-PFC" on "executes all SMARDDA-PFC codes to calculate power" \
        3 "Ignore previous" on "performs operation(s) 1/2 overwriting files" \
        2> $ifile
        retvalsm=$?
        case $retvalsm in
          0) reply=$(cat $ifile);;
          1) dialog --infobox "Abort" 6 50;sleep 1;break;;
          *) dialog --infobox "Unexpected code $retvalsm - try again (5t)" 6 50;sleep 1;break 2;;
        esac #retvalsm
        if [ "$reply" == "" ] || [ "${reply:0:1}" == "3" ] ; then
          dialog --infobox "No action requested" 6 50;sleep 1;break;fi
        case $reply in
          1) sw=1+0 ;;
          2) sw=0+1 ;;
          "1 2") sw=1+1 ;;
          "1 3") sw=2+0 ;;
          "2 3") sw=0+2 ;;
          "1 2 3") sw=2+2 ;;
          *) dialog --infobox "Option $reply not implemented yet" 6 50;sleep 1;break;;
        esac #reply
        #Advance to generation of ctl files and possibly run software
        $INT/rsmiter.bash tempdir=$tempdir sw=$sw \
        flds=$flds bopt=$bopt cal=$cal eqfile=$eqfile eqid=$eqid  \
        gres=$gres  gshad=$gshad bdryl=$bdryl powe=$powe decl=$decl \
        vfld=$vfld arip=$arip mrip=$mrip nzetp=$nzetp nbin=$nbin \
        2>&1 | dialog --programbox "rsmardda-pfc execution - do not close box until see \
        \n***rsmardda-pfc completed successfully***" 30 80
        smdone=$?
        if [ ${sw:2:3} == 0 ] ;then
          respow="$gres"_pow; shadhds="$gshad"_hds
          dialog --title "ctl file names" --msgbox \
          "$gres.ctl input to geoq for target geometry \
          \n$gshad.ctl input to geoq for shadow geometry \
          \n$shadhds.ctl input to hdsgen (shadow geometry) \
          \n$respow.ctl input to powcal (target geometry)" 10 70
        fi
      done #smdone
      ;;
      ##Analyse output
      6t) #convert ps to png so ParaView can read
      afiles=($(ls -L 2> /dev/null  *.ps)) ; nps=${#afiles[@]}
      pngfile=null
      if [ ! $nps == 0 ] ; then
        for i in ${afiles[@]}; do fi=${i%.ps}
          if [ ! -e $fi.png ] || [ $fi.ps -nt $fi.png ] ; then
            convert -rotate 90 $fi.ps $fi.png 2> /dev/null
            if [ $? == 0 ] ; then pngfile=$fi.png;fi
          fi
        done
      fi
      dataf="$gres"_pow_powx.vtk
      if [ -e $dataf ] ; then
        dialog --infobox "Launching ParaView" 6 50;sleep 1
        paraview --data=$dataf
      elif  [ ! $pngfile == null ] ; then 
        dialog --infobox "Launching ParaView to view 2-D graphics" 6 50;sleep 1
        paraview --data=$pngfile
      else
        dialog --msgbox "No output found - try again (6t)" 6 50
      fi
      ;;
      *) dialog --msgbox "Unrecognised selection from top menu" 6 50;;
    esac #top
  done #top
done #main
# end loop over top menu entries
