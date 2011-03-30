rootdir=$1
base=$2
outputdir=./$3
mkdir -p $outputdir
rm -f $outputdir/dens
asph="01 02 03 04 05 06 07 08 09 10"
sample="1 2 3 4"
for a in $asph
do
for s in $sample
do
bash bin/multidensity.sh $rootdir/asph_1.$a"_0.0"/$s/$base* > $outputdir/d$a"_"$s
done
	avg 1 $outputdir/d$a* > $outputdir/avg$a
	temp=` head -n24 $outputdir/avg$a | tail -n1`
	echo 1.$a $temp >> $outputdir/dens
done

