outputdir="densities"
mkdir -p $outputdir
#asph="1 2 3 4 5 6 7 8 9"
asph="0"
sample="1 2 3 4 5 6 7 8 9 10"
for a in $asph
do
for s in $sample
do
bash bin/multidensity.sh results/mono_volume/asph_1.0$a"_0.0"/$s/relax* > $outputdir/d$a"_"$s
done
done

for a in $asph
do
continue
	avg 1 $outputdir/d$a* > $outputdir/avg$a
done

for a in $asph
do
	temp=` head -n22 avg$a | tail -n1`
	echo 1.0$a $temp >> $outputdir/dens
done
