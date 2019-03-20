#! /bin/bash
# This file contains variable values that are important for running the scipts but that
# can change, and are therefore externalised from the script. Copy the varSettings_example.sh
# to varSettings.sh, then change verSettings.sh as necessary.


genomefile=${HOME}/genomeVer/ws260/sequence/c_elegans.PRJNA13758.WS260.genomic.fa  # on ubelix
#genomefile=/data/projects/p025/Jenny/genomeVer/ws265/c_elegans.PRJNA13758.WS265.genomic.fa  # on bioinformatics cluster
#genomefile= ${HOME}/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa # local

# path to cutadapt program
#cutadaptPath='/software/UHTS/Quality_control/cutadapt/1.13/bin/cutadapt'
# path to trimmomatic program
trimmomaticDIR='/software/UHTS/Analysis/trimmomatic/0.36/bin'
# path to trimmomatic adaptor file
trimAdapterFile='./TruSeq_2-3_PE.fa'
#trimAdapterFile='/software/UHTS/Analysis/trimmomatic/0.36/bin/adapters/TruSeq2-PE.fa'
#trimmomaticDIR := ${HOME}/Trimmomatic-0.36
#trimAdapterFile := ${trimmomaticDIR}/adapters/TruSeq_2-3_PE.fa

# install bwa-meth and bamutil programmes and then set the following variables in the ~/.bashrc 
# using the path to the installation location. e.g.
#export BWAMETHDIR=/home/ubelix/izb/semple/mySoftware/bwa-meth-master
#export BAMUTIL=/home/ubelix/izb/semple/mySoftware/bamUtil/bin/bam


picardDIR='/software/UHTS/Analysis/picard-tools/2.18.11/bin'

trimmed=TRUE

#start of the fastq filename that contains the sample name and is common between R1 and R2
sampleNames=( dS16N2gw dS20N2gw )

# date the sequencing was performed
seqDate="20190206"


########## don't edit below this line ###############

#export NUMSAMPLES=${#sampleNames[@]}
