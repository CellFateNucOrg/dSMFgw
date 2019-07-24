#! /usr/bin/bash

shortName=$1
longName=$2
shift 2

#make directory
mkdir -p  ./bigwig/${shortName}/ce11

#make genomes.txt file
echo -e "genome ce11\ntrackDb ce11/trackDb.txt\n" > ./bigwig/${shortName}/genomes.txt

#make hub.txt file
echo -e "hub "${shortName} > ./bigwig/${shortName}/hub.txt
echo -e "shortLabel "${shortName} >> ./bigwig/${shortName}/hub.txt
echo -e "longLabel "${longName} >> ./bigwig/${shortName}/hub.txt
echo -e "genomesFile genomes.txt" >> ./bigwig/${shortName}/hub.txt
echo -e "email peter.meister@izb.unibe.ch" >> ./bigwig/${shortName}/hub.txt
echo -e "descriptionUrl "${shortName}".html" >> ./bigwig/${shortName}/hub.txt

#make trackDb.txt file
if [[ -f ./bigwig/${shortName}/ce11/trackDb.txt ]]
then
	rm ./bigwig/${shortName}/ce11/trackDb.txt
fi
touch ./bigwig/${shortName}/ce11/trackDb.txt

for bwFile in "$@";
do
	fileBase=`basename ${bwFile}`
	shortFileLabel=${fileBase%.bw}
	echo -e "track "${shortFileLabel} >> ./bigwig/${shortName}/ce11/trackDb.txt
	echo -e "bigDataUrl "${fileBase} >> ./bigwig/${shortName}/ce11/trackDb.txt
	echo -e "shortLabel "${fileBase} >> ./bigwig/${shortName}/ce11/trackDb.txt
	echo -e "longLabel "${fileBase} >> ./bigwig/${shortName}/ce11/trackDb.txt
	echo -e "type bigWig" >> ./bigwig/${shortName}/ce11/trackDb.txt
	echo -e "maxHeightPixels 20" >> ./bigwig/${shortName}/ce11/trackDb.txt
	echo -e "color 75,0,130\n" >> ./bigwig/${shortName}/ce11/trackDb.txt
	mv $bwFile ./bigwig/${shortName}/ce11/
done

sed s%http://www.meister.izb.unibe.ch/ucsc/%http://www.meister.izb.unibe.ch/ucsc/${shortName}/ce11/%g < ./bigwig/urlsToUpload.txt > ./bigwig/urlsToUploadTmp.txt
mv ./bigwig/urlsToUploadTmp.txt ./bigwig/urlsToUpload.txt
