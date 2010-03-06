
files=$*
w=500
h=500

echo "Generating jpg frames "
for file in $files
do
	#/Users/reza/bin/test.perl -x-0.13 -y-0.25 -s1.3  -rx20  $file | render -size $w"x"$h -jpeg > $file.jpg 
	coord2pr3d -x0 -y0 -s1 $file | render -size $w"x"$h -jpeg > $file.jpg 

done

./encodejpg.sh -w=$w -h=$h 
