#! /usr/bin/bash
## mapping dSMF-gw sequences with bwa-meth
## required software: fastqc, cutadapt, trimmomatic, bwa-meth, samtools, picard, qualimap, bamutils, bamtools

###############################
########### VARIABLES #########
###############################
# variable are set in the varSettings.sh file and values read in from command line

source ./varSettings.sh

#the start of the basename of the fastq files
bname=$1

#the biological test group that you will later compare (used for preparing the QuasR input file)
testGroup=$2

# number of threads
numThreads=$3


toCombine=( /gpfs/homefs/izb/semple/myData/20190218_dSMFv016v020_N2gw_2x150pe/properpair /home/ubelix/izb/semple/myData/20190531_dSMFv016v020_N2gw_2x150pe/properpair )


#########################################################
#### merge datasets                                    ##
#########################################################
#
#mkdir -p aln
#samtools merge  ./aln/${bname}_${seqDate}.noOL.bam /gpfs/homefs/izb/semple/myData/20190218_dSMFv016v020_N2gw_2x150pe/properpair/aln/${bname}_20190218.noOL.bam /home/ubelix/izb/semple/myData/20190531_dSMFv016v020_N2gw_2x150pe/properpair/aln/${bname}_20190531.noOL.bam
#
#
#
## index bam files for methyldackel input
#samtools index aln/${bname}_${seqDate}.noOL.bam
#
#
#
#########################################################
#### Get stats on clipped reads                       ##
#########################################################
#
## get alignment stats again post-clipping
#samtools flagstat aln/${bname}_${seqDate}.noOL.bam  > qc/aln/report_flagstat_5_${bname}_${seqDate}_noOL.txt
#
### Get insert size statistics and plots with picard and qualimap post-filtering
#java -Xms1g -Xmx8g -jar ${picardDIR}/picard.jar CollectInsertSizeMetrics I=aln/${bname}_${seqDate}.noOL.bam O=qc/aln/${bname}_${seqDate}_noOL_picard_insert_size_metrics.txt H=qc/aln/${bname}_${seqDate}_noOL_picard_insert_size_histogram.pdf
#
#
#mkdir -p qc/aln/noOL_${bname}
#qualimap bamqc -bam aln/${bname}_${seqDate}.noOL.bam -c --java-mem-size=8G -outdir qc/aln/noOL_${bname} -outfile ${bname}_${seqDate}_noOL_report_qualimap.pdf -outformat PDF
#
#
#########################################################
#### get median coverage                               ##
#########################################################
#
#samtools depth -a aln/${bname}_${seqDate}.noOL.bam | cut -f3  > qc/aln/${bname}_${seqDate}_depthCol.txt
#
#echo "min\tmax\tmedian\tmean" > qc/aln/${bname}_${seqDate}_depthStats.txt
#./R/mmmm.r < qc/aln/${bname}_${seqDate}_depthCol.txt >> qc/aln/${bname}_${seqDate}_depthStats.txt
#        #depthStats=`./R/mmmm.r < $^`
#        #echo ${depthStats}
#
#rm qc/aln/${bname}_${seqDate}_depthCol.txt
#
#
#########################################################
#### make input file for R                             ##
#########################################################
#
#if [[ -e ./txt/bwameth_Aligned.txt ]]
#then
#	echo -e $PWD/aln/${bname}_${seqDate}.noOL.bam"\t"${bname}"\t"${testGroup} >> txt/bwameth_Aligned.txt
#else
#	mkdir -p txt
#	echo -e "FileName\tSampleName\tTestGroup" > txt/bwameth_Aligned.txt
#	echo -e $PWD/aln/${bname}_${seqDate}.noOL.bam"\t"${bname}"\t"${testGroup} >> txt/bwameth_Aligned.txt
#fi





#######################################################
## extract methylation with MethylDackel             ##
#######################################################

source ${HOME}/.bashrc
source ${CONDA_ACTIVATE} bwameth

mkdir -p methCalls
mkdir -p perRead
mkdir -p mbias

if [[ "$dataType" = "gw" ]]
then
	MethylDackel extract --CHH --CHG -o methCalls/${bname}_${seqDate} -d 1 -@ ${numThreads} ${genomefile} aln/${bname}_${seqDate}.noOL.bam
	MethylDackel perRead -@ ${numThreads} -o perRead/${bname}_${seqDate} ${genomefile} aln/${bname}_${seqDate}.noOL.bam
	MethylDackel mbias ${genomefile} aln/${bname}_${seqDate}.noOL.bam mbias/${bname}_${seqDate}
else
	# do not remove duplicates so change the default ignoreFlags (-F)
	MethylDackel extract --CHH --CHG -o methCalls/${bname}_${seqDate} -F 2816 -d 1 -@ ${numThreads} ${genomefile} aln/${bname}_${seqDate}.noOL.bam
	MethylDackel perRead -@ ${numThreads} -F 2816 -o perRead/${bname}_${seqDate} ${genomefile} aln/${bname}_${seqDate}.noOL.bam
	MethylDackel mbias ${genomefile} -F 2816 aln/${bname}_${seqDate}.noOL.bam mbias/${bname}_${seqDate}
fi


conda deactivate


#######################################################
## Prepare RDS of genome CG and GC motifs            ##
#######################################################

if [[ ! -f ${genomefile%.fa}.CGGC_motifs.RDS ]]
then
	Rscript getGenomeMotifs.R ${genomefile}
fi

