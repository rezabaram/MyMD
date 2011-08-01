#!/bin/bash

#USAGE EXAMPE: bash coord2xdr.sh `tree -if mono/ | grep /out$`
# this looks for all the "out" files and make their corresponding .xdr file for given resolutions
#to clean the created file:  rm -f `tree -if mono/ | grep .xdr`

bin=`dirname $0`
res="516"
input=$*

for i in $input
do
	inputdir=`dirname i`
	for r in $res
	do

	if [ -f $inputdir/$i$r.xdr ] ;
	then 
		continue
	fi

	cd $inputdir
	echo $i
	$bin/coord2xdr $i $r
	mv $i.xdr $inputdir/$i$r.xdr
	zip $inputdir/$i$r.xdr.zip $inputdir/$i$r.xdr
	rm $inputdir/$i$r.xdr
	done
done



#bin/VTKheader $2 > $input.vtk
#cat $input.xdr >> $input.vtk
