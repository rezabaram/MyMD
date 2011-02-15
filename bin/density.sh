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
		d=`bin/packing_density $p` 
		f=`calc $d + $f `
		f2=`calc $f2 + $d*$d`
	done

	f=`calc $f/$#`
	f2=`calc $f2/$#`
	std=`calc "sqrt($f2-$f*$f)" `
}

#asp="-0.9 -0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8 0.9"
asp="1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09"
aspw="0.0"
rootdir=$1

#rm -f density$aspw.dat
for a in $asp
do
	#files="`ls $rootdir/e_$a"_"$aspw/*/out00040` "
	files="`ls $rootdir/asph_$a"_"$aspw/*/out00035` "
	cal_density $files 
	echo $a $f $std >> density_gen.dat
done
