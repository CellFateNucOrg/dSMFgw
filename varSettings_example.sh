#! /bin/bash
# This file contains variable values that are important for running the scipts but that
# can change, and are therefore externalised from the script. Copy the varSettings_example.sh
# to varSettings.sh, then change verSettings.sh as necessary.

###################################################
# sample sepcific information (need to change for each run)
###################################################

#type of dSMF library
dataType=gw #either genome-wide ("gw") or amplicon ("amp")

#start of the fastq filename that contains the sample name and is common between R1 and R2
sampleNames=( dS16N2gw dS20N2gw ) #these are the beginning of the fastq filenames that specify the different samples
testGroups=( N2_dS16gw N2_dS20gw ) # this is the biological group that you wish to test for contrasts (it must have the same length as the number of samples)

# date the sequencing was performed
seqDate="20190206"

# if script crashes or runs out of time you can avoid re-running cutadapt and trimmomatic by setting this to TRUE
# normally this should be set to FALSE
trimmed=TRUE


###################################################
# paths to software or files required by software (need to change for each cluster environment)
###################################################

# genome file used for alignment
genomefile=${HOME}/genomeVer/ws260/sequence/c_elegans.PRJNA13758.WS260.genomic.fa  # on ubelix
#genomefile=/data/projects/p025/Jenny/genomeVer/ws265/c_elegans.PRJNA13758.WS265.genomic.fa  # on bioinformatics cluster
#genomefile= ${HOME}/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa # local

# path to trimmomatic program
trimmomaticDIR='/software/UHTS/Analysis/trimmomatic/0.36/bin'
# path to trimmomatic adaptor file
trimAdapterFile='./TruSeq_2-3_PE.fa'

# install bwa-meth and bamutil programmes and then set the following variables in the ~/.bashrc 
# using the path to the installation location. e.g.
#export BWAMETHDIR=/home/ubelix/izb/semple/mySoftware/bwa-meth-master
export BAMUTIL=/home/ubelix/izb/semple/mySoftware/bamUtil/bin/bam
picardDIR='/software/UHTS/Analysis/picard-tools/2.18.11/bin'



