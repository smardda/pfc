#! /bin/bash --norc
# produce file selection boxes for one first argument type file
if [ $# -eq 2 ] ; then
dout=$2
tempdir=${dout%/*}
fout=$tempdir/filsel.bash
elif [ $# -eq 1 ] ; then
dout=dout
fout=filsel.bash
else
echo "Usage $0 filetype [dir/workfile]"
fi
ftype=$1
if [ $ftype == vtk.1 ] ; then dtitl="results geometry file"
elif [ $ftype == vtk.2 ] ; then dtitl="shadowing geometry file"
else dtitl="one file";fi
rm -f $dout $fout
case "${ftype%%.[0-9]}" in
       vtk) 
       echo '#! /bin/bash --norc' > $fout
       echo 'dialog --menu "Select' $dtitl 'using arrow keys"  0 80 10 \' >> $fout
       for i in $(ls -L *.vtk 2>/dev/null);do echo "${i%%.vtk} $i" ' \';done >> $fout 
       ## remove unwanted files from list
       mv $fout $dout
       sed -e "/^track[xm][0-9][0-9]/d" -e "/_hds_/d" -e "/_pow_/d" -e "/_geo[a-z][a-z]*.vtk/d" < $dout >$fout
       echo "2> $dout" >> $fout
       chmod a+x $fout;;
       eqdsk) 
       echo '#! /bin/bash --norc' > $fout
       echo 'dialog --menu "Select magnetic equilibrium file  using arrow keys"  0 80 10 \' >> $fout
       for i in $(ls -L *.eqdsk 2>/dev/null);do echo "${i%%.eqdsk} $i" ' \';done >> $fout 
       echo "2> $dout" >> $fout
       chmod a+x $fout;;
       sp3) 
       echo '#! /bin/bash --norc' > $fout
       echo 'dialog --menu "Select vacuum magnetic file  using arrow keys"  0 80 10 \' >> $fout
       for i in $(ls -L *.sp3 2>/dev/null);do echo "${i%%.sp3} $i" ' \';done >> $fout 
       echo "2> $dout" >> $fout
       chmod a+x $fout;;
esac
