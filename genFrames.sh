
files=$*
w=500
h=500

echo "Generating jpg frames "
for file in $files
do
	#/Users/reza/bin/test.perl -x0 -y0 -s1  -rx75  $file | render -size $w"x"$h -jpeg > $file.jpg 
	if [ -f _$file.jpg ]
	then
		continue
	fi
	coord2pr3d -x0 -y0 -s1 $file | render -size $w"x"$h -jpeg > _$file.jpg 

done

./encodejpg.sh -w=$w -h=$h 
