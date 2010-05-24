
files=$*
w=800
h=800

echo "Generating jpg frames "
for file in $files
do
	#/Users/reza/bin/test.perl -x0 -y0 -s1  -rx75  $file | render -size $w"x"$h -jpeg > $file.jpg 
	coord2pr3d -x0 -y0 -s1 $file | render -size $w"x"$h -jpeg > $file.jpg 

done

./encodejpg.sh -w=$w -h=$h 
