rootdir=$1
base=$2
outputdir="densities"
mkdir -p $outputdir
asph="1 2 3 4 5 6 7 8 9"
#asph="0"
sample="1 2 3 4 5 6 7 8 9 10"
rm $outputdir/dens
for a in $asph
do
for s in $sample
do
continue
bash bin/multidensity.sh $rootdir/asph_1.0$a"_0.0"/$s/$base* > $outputdir/d$a"_"$s
done
	avg 1 $outputdir/d$a* > $outputdir/avg$a
	temp=` head -n42 $outputdir/avg$a | tail -n1`
	echo 1.0$a $temp >> $outputdir/dens
done

