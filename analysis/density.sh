
#a function for simple calculations
function calc(){
	echo "scale=6;$*" | bc
	}

packings=$*

wc $packings
f=0
for p in $packings;
do
	echo $p
	d=`./packing_density $p`
	f=`calc $d + $f `
done

calc $f/$#
