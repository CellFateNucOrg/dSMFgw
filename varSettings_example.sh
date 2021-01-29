##! /bin/bash
# This file contains variable values that are important for running the scipts but that
# can change, and are therefore externalised from the script. Copy the varSettings_example.sh
# to varSettings.sh, then change verSettings.sh as necessary.

###################################################
# sample sepcific information (need to change for each run)
###################################################

#type of dSMF library
dataType=amp #either genome-wide ("gw") or amplicon ("amp")

#start of the fastq filename that contains the sample name and is common between R1 and R2
sampleNames=( dSAlph1 dSAlph2 dSDRB1 dSDRB2 dSFlp1 dSFlp2 dSN2noT1 dSN2noT2 dSTrp1 dSTrp2 ) #these are the beginning of the fastq filenames that specify the different samples
testGroups=( Alph Alph DRB DRB Flp Flp N2 N2 Trp Trp  ) # this is the biological group that you wish to test for contrasts (it must have the same lenght as the number of sample, and must be completely included in respective sampleNames)


# date the sequencing was performed
seqDate=20191023
expName=dSampchem

# if script crashes or runs out of time you can avoid re-running cutadapt and trimmomatic by setting this to TRUE
# normally this should be set to FALSE
trimmed=FALSE
aligned=FALSE



###################################################
# paths to software or files required by software (need to change for each cluster environment)
###################################################

genomeVer=WS235
# genome file used for alignment
genomefile=/home/ubelix/izb/bi18k694/genomeversion/ws270/c_elegans.PRJNA13758.WS270.genomic.fa  # on ubelix
#genomefile=/home/ubelix/izb/semple/genomeVer/ws260/sequence/c_elegans.PRJNA13758.WS260.genomic.fa
#genomefile=/data/projects/p025/Jenny/genomeVer/ws265/c_elegans.PRJNA13758.WS265.genomic.fa  # on bioinformatics cluster
#genomefile= ${HOME}/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa # local

# path to trimmomatic program
trimmomaticDIR='/software/UHTS/Analysis/trimmomatic/0.36/bin'
# path to trimmomatic adaptor file
trimAdapterFile=/software/UHTS/Analysis/trimmomatic/0.36/bin/adapters/TruSeq3-PE-2.fa

# install bwa-meth and bamutil programmes and then set the following variables in the ~/.bashrc 
# using the path to the installation location. e.g.
#export BWAMETHDIR=/home/ubelix/izb/semple/mySoftware/bwa-meth-master
export BAMUTIL=/home/ubelix/izb/bi18k694/software/bamUtil/bin/bam
#export BAMUTIL=/home/ubelix/izb/semple/mySoftware/bamUtil/bin/bam
picardDIR='/software/UHTS/Analysis/picard-tools/2.18.11/bin'

