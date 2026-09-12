#! /bin/bash
# first wall silhouette : convert list of (R,Z) points in m to
# lines in (X,0,Z) in mm in vtk format, and wrap to close
if [ $# -ne 1 ] ; then        # Eactly  1 argument must be present.
   echo "Usage: $0 <txt file root name>"
   exit
fi
rm -f numtxt2vtk.txt $1.vtk $1.vtk_bit
echo -e "# vtk DataFile Version 3.0 \nsilhouette \nASCII \nDATASET UNSTRUCTURED_GRID" > $1.vtk
cat $1.txt | awk 'END{ print NR }'> numtxt2vtk.txt
nline=$(cat numtxt2vtk.txt | awk 'END{ print $1 }')
nlinep=$(cat numtxt2vtk.txt | awk 'END{ print $1+1 }')
nline3=$(cat numtxt2vtk.txt | awk 'END{ print 3*$1 }')
echo "POINTS     $nlinep float" >> $1.vtk
cat $1.txt | awk '{ print $1*1000 "  0.  " $2*1000}'  >> $1.vtk_bit
cat $1.vtk_bit  >> $1.vtk
# repeat 1st point as last to wrap
head -1 $1.vtk_bit  >> $1.vtk
echo " " >> $1.vtk
echo "CELLS $nline $nline3" >> $1.vtk
cat $1.txt | awk '{ print "2  "NR-1" "NR }'  >> $1.vtk
echo " " >> $1.vtk
echo "CELL_TYPES $nline" >> $1.vtk
cat $1.txt | awk '{ print "3" }'  >> $1.vtk
