#! /usr/bin/bash
## mapping dSMF-gw sequences with bwa-meth
## required software: fastqc, cutadapt, trimmomatic, bwa-meth, samtools, picard, qualimap

###############################
########### VARIABLES #########
###############################
# variable are set in the varSettings.sh file and values read in from command line

source ./varSettings.sh

#the start of the basename of the fastq files
bname=$1

# number of threads
numThreads=$4

# get foward and reverse read files for this sample
fileList=( `ls ../rawData/${bname}*.fastq.gz` )

#######################################################
## get initial read stats                            ##
#######################################################

#run fastqc on sequences
mkdir -p ./fastQC/rawData

#for f in ${fileList[@]}:
#do
#	fastqc -t ${numThreads} $f -o ./fastQC/rawData
#done

fastqc ${fileList[@]} -o ./fastQC/rawData

#######################################################
## trim adaptors with cutadapt                       ##
#######################################################

# use cutadapt to trim
mkdir -p cutadapt
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
                -o cutadapt/${bname}_${seqDate}_R1.fastq.gz -p cutadapt/${bname}_${seqDate}_R2.fastq.gz \
                ${fileList[@]} 


#redo fastQC on trimmed reads
mkdir -p fastQC/cutadapt
fastqc cutadapt/${bname}_${seqDate}_R?.fastq.gz -o ./fastQC/cutadapt 


#######################################################
## quality trim reads with Trimmomatic               ##
#######################################################

#graphical parameter for bash shell
export DISPLAY=:0

#use trimmomatic to trim
mkdir -p trim
mkdir -p fastQC/trim
echo $numThreads
echo $trimmomaticDIR

java -Xms1g -Xmx5g -jar ${trimmomaticDIR}/trimmomatic-0.36.jar PE -threads ${numThreads} cutadapt/${bname}_${seqDate}_R1.fastq.gz cutadapt/${bname}_${seqDate}_R2.fastq.gz trim/${bname}_${seqDate}_forward_paired.fq.gz trim/${bname}_${seqDate}_forward_unpaired.fq.gz trim/${bname}_${seqDate}_reverse_paired.fq.gz trim/${bname}_${seqDate}_reverse_unpaired.fq.gz ILLUMINACLIP:${trimAdapterFile}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 2> fastQC/trim/report_${bname}_${seqDate}_trimmomatic.txt

# redo fastQC on trimmed reads	
fastqc trim/${bname}_${seqDate}_*.fq.gz -o fastQC/trim

#######################################################
## align to genome with BWA-meth and convert to bam  ##
#######################################################

source activate bwaMeth
	
# convert and index genome file for bwameth alignment
if [[ ! -f ${genomefile}.bwameth.ct ]];
then
	${BWAMETH} index ${genomefile}
fi

# align sequences to meth converted genome with bwameth
mkdir -p aln
${BWAMETH} --threads ${numThreads} --reference ${genomefile} trim/${bname}_${seqDate}_forward_paired.fq.gz trim/${bname}_${seqDate}_reverse_paired.fq.gz > aln/${bname}_${seqDate}.sam

source deactivate

# use samtools to convert to bam and sort
samtools view -u aln/${bname}_${seqDate}.sam | samtools sort -o aln/${bname}_${seqDate}.sorted.bam -@ ${numThreads}
rm aln/${bname}_${seqDate}.sam


#######################################################
## Get alignment stats                               ##
#######################################################

# 	get alignment stats
mkdir -p fastQC/aln/prefilt
samtools flagstat  aln/${bname}_${seqDate}.sorted.bam > fastQC/aln/prefilt/report_${bname}_${seqDate}_flagstats.txt
samtools stats aln/${bname}_${seqDate}.sorted.bam > fastQC/aln/prefilt/report_${bname}_${seqDate}_stats.txt
	
#samtools view -cF 0x100 accepted_hits.bam

# Get insert size statistics and plots with picard and qualimap (set path to $QUALIMAP in .bash_profile)
java -Xms1g -Xmx5g -jar picard.jar CollectInsertSizeMetrics I=aln/${bname}_${seqDate}.sorted.bam \
  O=fastQC/aln/prefilt/${bname}_${seqDate}_picard_insert_size_metrics.txt \
  H=fastQC/aln/prefilt/${bname}_${seqDate}_picard_insert_size_histogram.pdf

qualimap bamqc -bam aln/${bname}_${seqDate}.sorted.bam -c -outdir fastQC/aln/prefilt -outfile ${bname}_${seqDate}_report_qualimap.pdf -outformat PDF


#######################################################
## Mark duplicates and filter reads. Then redo stats ##
#######################################################

# mark duplicates with picard (path to picard should be set in $PICARD variable in .bash_profile or in session)
mkdir -p fastQC/aln/postfilt
java -Xmx5g -jar picard.jar MarkDuplicates I=aln/${bname}_${seqDate}.sorted.bam O=aln/${bname}_${seqDate}.dup.bam M=fastQC/aln/postfilt/report_${bname}_${seqDate}_picard.txt

# 	remove mitochondrial reads
samtools view -q 30 -F 1804 -b aln/${bname}_${seqDate}.dup.bam > aln/${bname}_${seqDate}.filt.bam
rm aln/${bname}_${seqDate}.dup.bam

# NOTE: sam flag 1804 means the following:
# read unmapped
# mate unmapped
# not primary alignment
# read fails platform/vendor quality checks
# read is PCR or optical duplicate

# 	get alignment stats again post-filtering
samtools flagstat aln/${bname}_${seqDate}.filt.bam  > fastQC/aln/postfilt/repot_aln/${bname}_${seqDate}_flagstats.txt

samtools stats aln/${bname}_${seqDate}.filt.bam  > fastQC/aln/postfilt/report_${bname}_${seqDate}_stats.txt 

# Get insert size statistics and plots with picard and qualimap post-filtering
java -Xms1g -Xmx5g -jar picard.jar CollectInsertSizeMetrics I=aln/${bname}_${seqDate}.filt.bam \
  O=fastQC/aln/postfilt/${bname}_${seqDate}_picard_insert_size_metrics.txt \
  H=fastQC/aln/postfilt/${bname}_${seqDate}_picard_insert_size_histogram.pdf

qualimap bamqc -bam aln/${bname}_${seqDate}.filt.bam -c -outdir fastQC/aln/postfilt -outfile ${bname}_${seqDate}_report_qualimap.pdf -outformat PDF



#######################################################
## get median coverage 				     ##
#######################################################

samtools depth -a aln/${bname}_${seqDate}.filt.bam | cut -f3  > fastQC/aln/postfilt/${bname}_${seqDate}_depthCol.txt

echo "min\tmax\tmedian\tmean" > fastQC/aln/postfilt/${bname}_${seqDate}_depthStats.txt
./R/mmmm.r < fastQC/aln/postfilt/${bname}_${seqDate}_depthCol.txt >> fastQC/aln/postfilt/${bname}_${seqDate}_depthStats.txt
	#depthStats=`./R/mmmm.r < $^`
	#echo ${depthStats}
	#echo "${depthStats}" >> $@


