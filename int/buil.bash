#! /bin/bash --norc
# produce file selection boxes for first argument type files
if [ $# -eq 2 ] ; then
dout=$2
tempdir=${dout%/*}
fout=$tempdir/filsel.bash
else
dout=dout
fout=filsel.bash
fi
ulr=_R;ulz=_Z
rm -f $dout $fout
test=$1
case "${test%%.[0-9]}" in
       vtk) 
       echo '#! /bin/bash --norc' > $fout 
       echo 'dialog --buildlist "Select files by mouse-click or with space-bar"  0 80 10 \' >> $fout
       for i in $(ls -L *.vtk 2>/dev/null);do echo "${i%%.vtk} $i" 'off \';done  >> $fout
       ## remove unwanted files from list
       mv $fout $dout
       sed -e "/^track[xm][0-9][0-9]/d" -e "/_hds_/d" -e "/_pow_/d" -e "/_geo[a-z][a-z]*.vtk/d" < $dout >$fout
       echo "2> $dout" >> $fout
       chmod a+x $fout;;
       dat) 
       echo '#! /bin/bash --norc' > $fout 
       echo 'dialog --buildlist "Select files by mouse-click or with space-bar"  0 80 10 \' >> $fout
       for i in $(ls -L *.dat 2>/dev/null);do echo "${i%%.dat} $i" 'off \';done  >> $fout
       echo "2> $dout" >> $fout
       chmod a+x $fout;;
       rdat) 
       echo '#! /bin/bash --norc' > $fout 
       #echo 'dialog --extra-button --extra-label "Help" \
       echo 'dialog \
       --buildlist "Select one _R.dat file by mouse-click or with space-bar"  0 80 10 \' >> $fout
       for i in $(ls -L *_R*.dat 2>/dev/null);do echo "${i%%.dat} $i" 'off \';done >> $fout
       echo "2> $dout" >> $fout
       chmod a+x $fout;;
       Xdat)
       echo '#! /bin/bash --norc' > $fout 
       #echo 'dialog --extra-button --extra-label "Help" \
       echo 'dialog \
       --buildlist "Select files by mouse-click or with space-bar"  0 80 10 \' >> $fout
       afiles=($(ls -L 2> /dev/null  *.dat))
       bfiles=
       for i in ${afiles[@]}; do
          if [[ $i == *"$ulr"* ]] || [[ $i == *"$ulz"* ]]; then continue;else bfiles=(${bfiles[@]} $i);fi
       done
       for i in ${bfiles[@]};do echo "${i%%.dat} $i" 'off \';done >> $fout
       echo "2> $dout" >> $fout
       chmod a+x $fout;;
esac
