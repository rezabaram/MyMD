#!/bin/bash
#a function for simple calculations
function calc(){
	echo "scale=8;$*" | bc
	}


function cal_density(){
	packings=$*
	f=0
	f2=0
	for p in $packings;
	do
		echo $p
		d=`./packing_density $p` 
		f=`calc $d + $f `
		f2=`calc $f2 + $d*$d`
	done

	f=`calc $f/$#`
	f2=`calc $f2/$#`
	std=`calc "sqrt($f2-$f*$f)" `
}

asp="-0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8"
aspw="0.1"
rootdir=$1

#rm -f density$aspw.dat
for a in $asp
do
	files="`ls $rootdir/e_$a"_"$aspw/*/out00031` "
	cal_density $files 
	echo $a $f $std >> densityT30_$aspw.dat
done
