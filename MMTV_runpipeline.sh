#!/bin/bash

# Run pipeline for all Marco data

for file in *.fna
	do
	base=${file/.fna/''}
	base=${base/*\./''}
	base=${base/' '/''}
	echo $file
	echo $base
	mkdir $base
	mv $file $base
	cd $base
	~/marco/pipeline/pipeline_MMTV.sh $file $base 1
	cd ..
done
