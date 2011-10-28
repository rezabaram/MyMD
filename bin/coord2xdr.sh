#!/bin/bash

#USAGE EXAMPE: bash coord2xdr.sh `tree -if mono/ | grep /out$`
# this looks for all the "out" files and make their corresponding .xdr file for given resolutions
#to clean the created file:  rm -f `tree -if mono/ | grep .xdr`

res="256"
input=$*
suffix="z"

current=`pwd`
bin=`dirname $0`
for i in $input
do
	cd $current
	inputdir=$current/`dirname $i`
	inputname=`basename $i`
	for r in $res
	do

	if [ -f $inputdir/$inputname$r$suffix.xdr ] ;
	then 
		continue
	fi

	$bin/coord2xdr $inputdir/$inputname $r
	cd $inputdir
	mv $inputname.xdr $inputname$r$suffix.xdr
	zip $inputname$r$suffix.xdr.zip $inputname$r$suffix.xdr
	rm $inputname$r$suffix.xdr
	done
done



#bin/VTKheader $2 > $input.vtk
#cat $input.xdr >> $input.vtk
