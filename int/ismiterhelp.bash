#! /bin/bash --norc
#help and version handling
ixit=0
for i in $*; do tst=$i #echo $tst
if [ "${tst#-}" != $tst ] ; then ixit=1
  if [ "${tst##--}" == version ] ;  then
    echo "Introduces/interacts with SMITER 1.1"
    echo "(c) UKAEA and IO 2016"
  elif [ "${tst##--}" == help ] ;  then
    echo "Introduces the occasional linux user to fieldline tracing"
    echo "using SMITER, SMARDDA for ITER"
    echo "with web-links to extensive more advanced material."
    echo "Report issues using the IMAS system."
  elif [ "${tst##--}" == tmp ] ;  then
    ixit=2
  elif [ "${tst##--}" == force ] ;  then
    ixit=3
  else
    echo "Valid options are --help and --version"
    echo "EXPERT options are --tmp and --force"
  fi
fi
done
if [ $ixit == 2 ] ; then exit 2;fi
if [ $ixit == 3 ] ; then exit 3;fi
if [ $ixit != 0 ] ; then exit 1;fi
