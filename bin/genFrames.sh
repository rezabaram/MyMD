
files=$*
w=500
h=500

rasterconverter=$HOME/workstation/MD/bin/coord_convert
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
	$rasterconverter -r $file > $file".r"
	coord2pr3d -x0 -y0 -z0  -s1 $file".r" | render -size $w"x"$h -jpeg > _$file.jpg 
	rm -f $file".r"
	#coord2pr3d -x0 -y0 -s1 $file | render -size $w"x"$h -jpeg > _$file.jpg 
	

done

bin/encodejpg.sh -w=$w -h=$h 
