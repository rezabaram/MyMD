#!/bin/bash

if [ -f "out00000" ]
then
	echo "Please remove old output files if you want to run a new simulation."
	exit 1;
fi
./a.out 
