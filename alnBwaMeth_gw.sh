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

# get foward and reverse read files for this sample
fileList=( `ls ../rawData/${bname}*.fastq.gz` )

# setup up a conditional statement to avoid repeating already executed steps
if [[ "$trimmed" = "FALSE" ]] 
then

#######################################################
## get initial read stats                            ##
#######################################################

#run fastqc on sequences
mkdir -p ./fastQC/rawData

fastqc -t ${numThreads} ${fileList[@]} -o ./fastQC/rawData

#######################################################
## trim adaptors with cutadapt                       ##
#######################################################

# use cutadapt to trim
mkdir -p cutadapt
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
                -o cutadapt/${bname}_${seqDate}_R1.fastq.gz -p cutadapt/${bname}_${seqDate}_R2.fastq.gz \
                ${fileList[@]} 


#redo fastQC on cut reads
mkdir -p fastQC/cutadapt
fastqc cutadapt/${bname}_${seqDate}_R?.fastq.gz -o ./fastQC/cutadapt 

fi # end trimmed brackets

fi # end trimmed brackets

#######################################################
## quality trim reads with Trimmomatic               ##
#######################################################

#graphical parameter for bash shell
export DISPLAY=:0

#use trimmomatic to trim
mkdir -p trim
mkdir -p fastQC/trim

java -Xms1g -Xmx8g -jar ${trimmomaticDIR}/trimmomatic-0.36.jar PE -threads ${numThreads} cutadapt/${bname}_${seqDate}_R1.fastq.gz cutadapt/${bname}_${seqDate}_R2.fastq.gz -baseout trim/${bname}_${seqDate}.fq.gz ILLUMINACLIP:${trimAdapterFile}:2:30:10:3:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 2> fastQC/trim/report_${bname}_${seqDate}_trimmomatic.txt

# redo fastQC on trimmed reads	
fastqc trim/${bname}_${seqDate}_*.fq.gz -o fastQC/trim


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
${BWAMETH} --threads ${numThreads} --reference ${genomefile} trim/${bname}_${seqDate}_1P.fq.gz trim/${bname}_${seqDate}_2P.fq.gz > aln/${bname}_${seqDate}.sam
#${BWAMETH} --threads ${numThreads} --reference ${genomefile} cutadapt/${bname}_${seqDate}_R1.fastq.gz cutadapt/${bname}_${seqDate}_R2.fastq.gz > aln/${bname}_${seqDate}.sam

source deactivate

# convert to bam file to save space
samtools view -b -o aln/${bname}_${seqDate}.bam -@ ${numThreads} aln/${bname}_${seqDate}.sam
rm aln/${bname}_${seqDate}.sam



#######################################################
## Get alignment stats                               ##
#######################################################

# get alignment stats (need to sort temporarily first)
mkdir -p fastQC/aln
samtools sort -o aln/${bname}_${seqDate}.sort.bam -@ ${numThreads} aln/${bname}_${seqDate}.bam 
samtools flagstat aln/${bname}_${seqDate}.sort.bam > fastQC/aln/report_flagstat_1_${bname}_${seqDate}_bam.txt
#rm aln/${bname}_${seqDate}.sort.bam



#######################################################
## Remove duplicates with Picard                     ##
#######################################################

#sort by query name for duplicate removal
java -Xms1g -Xmx8g -jar ${picardDIR}/picard.jar SortSam I=aln/${bname}_${seqDate}.bam O=aln/${bname}_${seqDate}.qsort.bam SORT_ORDER=queryname TMP_DIR=${TMPDIR}


# remove duplicates with picard (path to picard should be set in $PICARD variable in .bash_profile or in session)
# Note that to mark unmapped mates of mapped records and supplementary/secondary alignments as duplicates the bam
# file must be querysorted (by name) not by coordinate. 
java -Xms1g -Xmx8g -jar ${picardDIR}/picard.jar MarkDuplicates I=aln/${bname}_${seqDate}.qsort.bam O=aln/${bname}_${seqDate}.noDup.bam M=fastQC/aln/report_${bname}_${seqDate}_picard.txt REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true  ASSUME_SORT_ORDER=queryname TMP_DIR=${TMP}

#-XX:ParallelGCThreads=${numThreads}

# if you wish to see the duplicate tags using the following option instead of REMOVE_DUPLICATES
# The tags will go in the DT field
#--TAGGING_POLICY=All
# the you would have to use samtools view -F 1024

rm aln/${bname}_${seqDate}.qsort.bam

#######################################################
## Sort and get stats                                ##
#######################################################

# sort by position
samtools sort -o aln/${bname}_${seqDate}.sorted.bam -@ ${numThreads} aln/${bname}_${seqDate}.noDup.bam

# get alignment stats
samtools flagstat  aln/${bname}_${seqDate}.sorted.bam > fastQC/aln/report_flagstat_2_${bname}_${seqDate}_noDup.txt
	
rm aln/${bname}_${seqDate}.noDup.bam



#######################################################
## Filter by mapping score, orientation, same chr    ##
#######################################################

# keep only reads that have Q>=30, both are mapped and in a FR or RF orientation.
bamtools filter -in aln/${bname}_${seqDate}.sorted.bam -out aln/${bname}_${seqDate}.filt2.bam  -script myBamFilters.json

# keep only reads that map to the same chromosome
# write header to file temporarily
samtools view -H aln/${bname}_${seqDate}.filt2.bam >  ${bname}_${seqDate}.header.sam
# extract rows where the 7th column has "=" (same chromosome) and the 9th column has insert length!=0 (found in wrongly oriented pairs), 
# and combine with header into a new bam file.
samtools view aln/${bname}_${seqDate}.filt2.bam | awk '($7=="=" && $9!="0" )' | cat ${bname}_${seqDate}.header.sam - | samtools view -b - -o aln/${bname}_${seqDate}.filt3.bam
rm ${bname}_${seqDate}.header.sam

#rm aln/${bname}_${seqDate}.sorted.bam

# for simple samtools filtering:
# NOTE: sam flag 3852 (if want supl alignments, use 1804) means excluding and of the following:
# 4    read unmapped
# 8    mate unmapped
# 256  not primary alignment
# 512  read fails platform/vendor quality checks
# 1024 read is PCR or optical duplicate
# 2048 Supplementary alignment



########################################################
### Get stats on filtered reads                       ##
########################################################

# get alignment stats again post-filtering
samtools flagstat aln/${bname}_${seqDate}.filt2.bam  > fastQC/aln/report_flagstat_3_${bname}_${seqDate}_filt2.txt
samtools flagstat aln/${bname}_${seqDate}.filt3.bam  > fastQC/aln/report_flagstat_4_${bname}_${seqDate}_filt3.txt

## Get insert size statistics and plots with picard and qualimap post-filtering
java -Xms1g -Xmx8g -jar ${picardDIR}/picard.jar CollectInsertSizeMetrics I=aln/${bname}_${seqDate}.filt2.bam O=fastQC/aln/${bname}_${seqDate}_filt2_picard_insert_size_metrics.txt H=fastQC/aln/${bname}_${seqDate}_filt2_picard_insert_size_histogram.pdf
java -Xms1g -Xmx8g -jar ${picardDIR}/picard.jar CollectInsertSizeMetrics I=aln/${bname}_${seqDate}.filt3.bam O=fastQC/aln/${bname}_${seqDate}_filt3_picard_insert_size_metrics.txt H=fastQC/aln/${bname}_${seqDate}_filt3_picard_insert_size_histogram.pdf


mkdir -p fastQC/aln/file2_${bname}
qualimap bamqc -bam aln/${bname}_${seqDate}.filt2.bam -c --java-mem-size=8G -outdir fastQC/aln/file2_${bname} -outfile ${bname}_${seqDate}_filt2_report_qualimap.pdf -outformat PDF
mkdir -p fastQC/aln/file3_${bname}
qualimap bamqc -bam aln/${bname}_${seqDate}.filt3.bam -c --java-mem-size=8G -outdir fastQC/aln/filt3_${bname} -outfile ${bname}_${seqDate}_filt3_report_qualimap.pdf -outformat PDF

#rm aln/${bname}_${seqDate}.filt2.bam


########################################################
### index bam files for QuasR input                   ##
########################################################
samtools index aln/${bname}_${seqDate}.sorted.bam
samtools index aln/${bname}_${seqDate}.filt2.bam
samtools index aln/${bname}_${seqDate}.filt3.bam


########################################################
### clip overlap between reads                        ##
########################################################

${BAMUTIL} clipOverlap --in aln/${bname}_${seqDate}.filt3.bam --out aln/${bname}_${seqDate}.noOL.bam --stats &> fastQC/aln/clipOl_${bname}_${seqDate}.txt

# index bam files for QuasR input
samtools index aln/${bname}_${seqDate}.noOL.bam

#rm aln/${bname}_${seqDate}.filt3.bam

########################################################
### Get stats on filtered reads                       ##
########################################################

# get alignment stats again post-filtering
samtools flagstat aln/${bname}_${seqDate}.noOL.bam  > fastQC/aln/report_flagstat_5_${bname}_${seqDate}_noOL.txt

## Get insert size statistics and plots with picard and qualimap post-filtering
java -Xms1g -Xmx8g -jar ${picardDIR}/picard.jar CollectInsertSizeMetrics I=aln/${bname}_${seqDate}.noOL.bam O=fastQC/aln/${bname}_${seqDate}_noOL_picard_insert_size_metrics.txt H=fastQC/aln/${bname}_${seqDate}_noOL_picard_insert_size_histogram.pdf

mkdir -p fastQC/aln/noOL_${bname}
qualimap bamqc -bam aln/${bname}_${seqDate}.noOL.bam -c --java-mem-size=8G -outdir fastQC/aln/noOL_${bname} -outfile ${bname}_${seqDate}_noOL_report_qualimap.pdf -outformat PDF


########################################################
### get median coverage                               ##
########################################################

samtools depth -a aln/${bname}_${seqDate}.noOL.bam | cut -f3  > fastQC/aln/${bname}_${seqDate}_depthCol.txt

echo "min\tmax\tmedian\tmean" > fastQC/aln/${bname}_${seqDate}_depthStats.txt
./R/mmmm.r < fastQC/aln/${bname}_${seqDate}_depthCol.txt >> fastQC/aln/${bname}_${seqDate}_depthStats.txt
        #depthStats=`./R/mmmm.r < $^`
        #echo ${depthStats}
        #echo "${depthStats}" >> $@

rm fastQC/aln/${bname}_${seqDate}_depthCol.txt


########################################################
### get multiqc report                                ##
########################################################

#multiqc ./fastQC


########################################################
### make input file for quasR                         ##
########################################################

if [[ -e ./txt/QuasR_Aligned.txt ]]
then
	echo -e $PWD/aln/${bname}_${seqDate}.noOL.bam"\t"${testGroup}"_"${bname} >> txt/QuasR_Aligned.txt
else
	mkdir -p txt
	echo -e "FileName\tSampleName" > txt/QuasR_Aligned.txt
	echo -e $PWD/aln/${bname}_${seqDate}.noOL.bam"\t"${testGroup}"_"${bname} >> txt/QuasR_Aligned.txt
fi
