
files=$*
w=500
h=500

echo "Generating jpg frames "
for file in $files
do
	#/Users/reza/bin/test.perl -x0 -y0 -s1  -rx75  $file | render -size $w"x"$h -jpeg > $file.jpg 
	filename=`basename $file`
	filedir=`dirname $file`
	if [ -f _$filename.jpg ]
	then
		continue
	fi
	analysis/coord_convert -r $file > $filename".r"
	bin/coord2pr3d -x0 -y0 -z0  -s1 $filename".r" | render -size $w"x"$h -jpeg > _$filename.jpg 
	rm -f $filename".r"
	#coord2pr3d -x0 -y0 -s1 $file | render -size $w"x"$h -jpeg > _$file.jpg 
	

done

./encodejpg.sh -w=$w -h=$h 
