#! /usr/bin/bash

module load vital-it
module add UHTS/Analysis/deepTools/2.5.4;


for bamFile in "$@";
do	
	echo $bamFile
	fileBase=`basename $bamFile`
	bwFile=${fileBase%noOL.bam}cov.bw
	echo $bwFile
	bamCoverage -b $bamFile -o ./bigwig/${bwFile} --binSize 10  --numberOfProcessors 4
done
