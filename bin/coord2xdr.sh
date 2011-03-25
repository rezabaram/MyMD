#!/bin/bash

if [  $# != 2 ] 
then
echo "USAGE: coord2sdr.sh <inputfile> <resolution>"
exit 1
fi 
input=$1
res=$2
bin/coord2xdr $input $res
#bin/VTKheader $2 > $input.vtk
#cat $input.xdr >> $input.vtk




