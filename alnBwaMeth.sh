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
numThreads=$2

# get foward and reverse read files for this sample
fileList=( `ls ../rawData/${bname}*.fastq.gz` )

if [[ ! trimmed ]] 
then

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

fi # end trimmed brackets

#######################################################
## remove duplicate fastq reads                      ##
#######################################################
#bname=dS20N2gw
#seqDate=20190206
mkdir -p uniq
gunzip cutadapt/${bname}_${seqDate}_R?.fastq.gz
ls cutadapt/${bname}_${seqDate}_R?.fastq > cutadapt/${bname}_${seqDate}.fofn
$FASTUNIQ -i cutadapt/${bname}_${seqDate}.fofn -t q -o uniq/${bname}_${seqDate}_R1.fastq -p uniq/${bname}_${seqDate}_R2.fastq

# recompress the files
#gzip cutadapt/${bname}_${seqDate}_R?.fastq
gzip uniq/${bname}_${seqDate}_R?.fastq

#redo fastQC on uniq reads
mkdir -p fastQC/uniq
fastqc uniq/${bname}_${seqDate}_R?.fastq.gz -o ./fastQC/uniq

#######################################################
## merge overlapping reads                           ##
#######################################################

flash2 -O -r 150 -f 300 -s 30 -M 100 -o ${bname}_${seqDate} -d cutadapt -z uniq/${bname}_${seqDate}_R1.fastq.gz uniq/${bname}_${seqDate}_R2.fastq.gz 2>&1 | tee fastQC/uniq/flash_${bname}_${seqDate}.log

##graphical parameter for bash shell
##export DISPLAY=:0
#
##use trimmomatic to trim
#mkdir -p trim
#mkdir -p fastQC/trim
#echo $numThreads
#echo $trimmomaticDIR
#echo $trimAdapterFile
#
#java -Xms1g -Xmx5g -jar ${trimmomaticDIR}/trimmomatic-0.36.jar PE -threads ${numThreads} cutadapt/${bname}_${seqDate}_R1.fastq.gz cutadapt/${bname}_${seqDate}_R2.fastq.gz trim/${bname}_${seqDate}_forward_paired.fq.gz trim/${bname}_${seqDate}_forward_unpaired.fq.gz trim/${bname}_${seqDate}_reverse_paired.fq.gz trim/${bname}_${seqDate}_reverse_unpaired.fq.gz ILLUMINACLIP:${trimAdapterFile}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 2> fastQC/trim/report_${bname}_${seqDate}_trimmomatic.txt
#
## redo fastQC on trimmed reads	
#fastqc trim/${bname}_${seqDate}_*.fq.gz -o fastQC/trim


if [[ ! TRUE  ]] 
then

#######################################################
## align to genome with BWA-meth and convert to bam  ##
#######################################################
# consider fastx-collapser to get rid of duplicates (works on SE reads!) http://hannonlab.cshl.edu/fastx_toolkit/

source activate bwaMeth
	
# convert and index genome file for bwameth alignment
if [[ ! -f ${genomefile}.bwameth.c2t ]]
then
	${BWAMETH} index ${genomefile}
fi

# align sequences to meth converted genome with bwameth
mkdir -p aln
${BWAMETH} --threads ${numThreads} --reference ${genomefile} trim/${bname}_${seqDate}_forward_paired.fq.gz trim/${bname}_${seqDate}_reverse_paired.fq.gz > aln/${bname}_${seqDate}.sam

source deactivate


# use samtools to convert to bam and sort by name for picard tools duplicate marking
samtools sort -n -o aln/${bname}_${seqDate}.bam -@ ${numThreads} aln/${bname}_${seqDate}.sam
#rm aln/${bname}_${seqDate}.sam


#######################################################
## Get alignment stats                               ##
#######################################################

# get alignment stats
mkdir -p fastQC/aln/prefilt
samtools sort -@ ${numThreads} aln/${bname}_${seqDate}.bam | samtools flagstat  > fastQC/aln/prefilt/report_${bname}_${seqDate}_flagstat_aln.txt



#######################################################
## Remove duplicates with Picard                     ##
#######################################################

# mark duplicates with picard (path to picard should be set in $PICARD variable in .bash_profile or in session)
# Note that to mark unmapped mates of mapped records and supplementary/secondary alignments as duplicates the bam
# file must be querysorted (by name) not by coordinate. Consider using SortSam from Picard if these are not being flagged
java -Xmx5g -jar ${picardDIR}/picard.jar MarkDuplicates \
        I=aln/${bname}_${seqDate}.sorted.bam \
        O=aln/${bname}_${seqDate}.dup.bam \
        M=fastQC/aln/postfilt/report_${bname}_${seqDate}_picard.txt \
        REMOVE_DUPLICATES=true

# if you wish to see the duplicate tags using the following option instead of REMOVE_DUPLICATES
# The tags will go in the DT field
#--TAGGING_POLICY=All
# the you would have to use samtools view -F 1024


#######################################################
## Sort and get alignment stats                      ##
#######################################################

# sort by position
samtools sort -o aln/${bname}_${seqDate}.sorted.bam -@ ${numThreads} aln/${bname}_${seqDate}.bam

# get alignment stats
mkdir -p fastQC/aln/prefilt
samtools flagstat  aln/${bname}_${seqDate}.sorted.bam > fastQC/aln/prefilt/report_${bname}_${seqDate}_flagstat_noDup.txt
	


#######################################################
## Filter by read orientation                        ##
#######################################################

# take reads only in the right orientation
samtools view -q 30 -F 3852 -f 97 -b aln/${bname}_${seqDate}.dup.bam > aln/${bname}_${seqDate}.fr1.bam
samtools view -q 30 -F 3852 -f 145 -b aln/${bname}_${seqDate}.dup.bam > aln/${bname}_${seqDate}.fr2.bam
samtools view -q 30 -F 3852 -f 81 -b aln/${bname}_${seqDate}.dup.bam > aln/${bname}_${seqDate}.rf1.bam
samtools view -q 30 -F 3852 -f 161 -b aln/${bname}_${seqDate}.dup.bam > aln/${bname}_${seqDate}.rf2.bam

listBams=( aln/${bname}_${seqDate}.fr1.bam aln/${bname}_${seqDate}.fr2.bam aln/${bname}_${seqDate}.rf1.bam aln/${bname}_${seqDate}.rf2.bam )

samtools merge -f aln/${bname}_${seqDate}.filt.bam  ${listBams[@]}

# remove duplicate reads
#samtools view -q 30 -F 3852 -b aln/${bname}_${seqDate}.dup.bam > aln/${bname}_${seqDate}.filt.bam
#rm aln/${bname}_${seqDate}.dup.bam
#rm ${listBams[@]}

# NOTE: sam flag 3852 (if want supl alignments, use 1804) means excluding and of the following:
# 4    read unmapped
# 8    mate unmapped
# 256  not primary alignment
# 512  read fails platform/vendor quality checks
# 1024 read is PCR or optical duplicate
# 2048 Supplementary alignment



#######################################################
## Get stats on filtered reads                       ##
#######################################################


mkdir -p fastQC/aln/postfilt
# 	get alignment stats again post-filtering
samtools flagstat aln/${bname}_${seqDate}.filt.bam  > fastQC/aln/postfilt/report_${bname}_${seqDate}_flagstat_filt.txt

# Get insert size statistics and plots with picard and qualimap post-filtering
java -Xms1g -Xmx5g -jar ${picardDIR}/picard.jar CollectInsertSizeMetrics I=aln/${bname}_${seqDate}.filt.bam \
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



#######################################################
## index bam files for QuasR input                   ##
#######################################################

samtools index aln/${bname}_${seqDate}.filt.bam



#######################################################
## clip overlap between reads                        ##
#######################################################

$BAMUTILDIR clipOverlap --in aln/${bname}_${seqDate}.filt.bam --out aln/${bname}_${seqDate}.noOL.bam  --stats  2>&1 fastQC/aln/clipOl_${bname}_${seqDate}.txt

# get alignment stats again post-filtering
samtools flagstat aln/${bname}_${seqDate}.noOL.bam  > fastQC/aln/postfilt/report_${bname}_${seqDate}_flagstat_noOL.txt

samtools index aln/${bname}_${seqDate}.noOL.bam


fi
