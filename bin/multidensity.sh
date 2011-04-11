#!/bin/bash

dir=`dirname $0`
for file in $*
do
	$dir/packing_density $file
done
