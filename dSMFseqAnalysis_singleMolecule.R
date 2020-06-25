## mapping with bwaMeth
####### Align bisulfite sequencing of genomeWide dSMF  to C elegans genome
# date: 2019-05-25
# author: Jennnifer Semple
#
# input required:
# bam files aligned with bwa-meth in aln/ directory
# bwameth_Aligned.txt with FileName SampleName headings in txt directory
# suggested sampleName structure:
# strainName_dSMFvXXXgw where XXX is the unique sample number, and strainName is the biological
# difference we wish to compare.
# Methyl-Dackel called methylation frequency (CHH CHG and CpG)
# Path to genome file and an RDS with unique non-overlapping CG and GC motifs in genome


###############################################################################################
#####################             SETUP                                 #######################
###############################################################################################


#############
# Libraries
#############


suppressMessages(library(rtracklayer))
suppressMessages(library("RColorBrewer"))
suppressMessages(library(GGally))
suppressMessages(library(gridExtra))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggpubr))
suppressMessages(library(methMatrix))



#####################
# Custom functions
#####################

source('./R/dSMF_auxiliaryFunctions.R')
source('./R/dSMF_bigwig.R')
source('./R/variableSettings.R')



##############
# Variables
##############
options("scipen"=8, "digits"=4)
options(device=pdf)

# see ./R/variableSettings.R file. Some of these variables need to be adjusted before running
# the script. the variableSettings_example.R file downloaded from the repo should be correctly
# filled out and saved without the _example extension.


# standard variables that should probably not be changed
minConversionRate=0.8
maxNAfraction=0.2


# assume motif file and bed files are located in the same directory as the genomeFile and their name is
# derived from it (as created by the getGenomeMotifs.R script
motifFile=gsub("\\.fa","\\.CGGC_motifs.RDS",genomeFile)
genomeMotifGR<-readRDS(motifFile)
bedFilePrefix=gsub("\\.fa","",genomeFile)


# assume sample table is stored in ./txt/bwameth_Aligned.txt
fileList<-read.table("txt/bwameth_Aligned.txt",stringsAsFactors=F,header=T)
fileList<-fileList[match(sampleNames, fileList$SampleName),]
write.table(fileList,file=paste0(path,"/txt/bwameth_Aligned.txt"), quote=F,
sep="\t",row.names=F)
samples<-fileList$SampleName


# read in task id from command line
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
	stop("Set number of SLURM tasks to match the number of samples")
}
taskID=as.numeric(args[1])
nThreads=as.numeric(args[2])

#################
# Directories
################

#setup directory structure this to desired location of your alignments
makeDirs(path,c("plots","methCalls","csv","bigwig","rds"))


if(dataType=="amp"){ # only execute of it is amplicon data
    ###############################################################################################
    #####################         Extract single read data for amplicons     ######################
    ###############################################################################################

    ###################################################
    # get single read matrices for all Cs
    ###################################################


    # to get individual methylation matrices for an amplicon
    # make a list (by sample) of lists of matrices (by amplicon)
    # e.g
    # [1] sample1
    #     [1] amplicon1 matrix of reads x Cpositions
    #     [2] amplicon2 matrix of reads x Cpositions
    # [2] sample2
    #     [1] amplicon1 matrix of reads x Cpositions
    #     [2] amplicon2 matrix of reads x Cpositions

    #######
    # Extract all C positions within each amplicon (getCmethMatrix function)
    #######
    regionType="rawAmp"

    allSampleMats<-getSingleMoleculeMatrices(sampleTable=fileList[taskID,], genomeFile=genomeFile, 
						regionGRs=amplicons, regionType=regionType, 
						genomeMotifGR=genomeMotifGR, minConversionRate=minConversionRate, 
					   	maxNAfraction=maxNAfraction, bedFilePrefix=NULL, path=path, 
						convRatePlots=TRUE, nThreads=nThreads)

    #saveRDS(allSampleMats,paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,"_", 
    #		samples[taskID],".rds"))

} # end of amplicon only loop





##############################################################################################
#####################    Extract single read data centered on amplicon TSS      ###############
###############################################################################################


# to get individual methylation matrices for an TSS
# make a list (by sample) of lists of matrices (by TSS)
# e.g
# [1] sample1
#     [1] TSS1 matrix of reads x Cpositions
#     [2] TSS2 matrix of reads x Cpositions
# [2] sample2
#     [1] TSS1 matrix of reads x Cpositions
#     [2] TSS2 matrix of reads x Cpositions

#######
# Extract all C positions within each amplicon
#######
if(dataType=="amp"){
    winSize=500
} else {
    winSize=300
}

regionType="ampTSS"
tssWin<-ampTSS
mcols(tssWin)$TSS<-start(tssWin)
tssWin<-resize(tssWin,width=winSize,fix="center")


allSampleMats<-getSingleMoleculeMatrices(sampleTable=fileList[taskID,], genomeFile=genomeFile, 
					regionGRs=tssWin,regionType=regionType, genomeMotifGR=genomeMotifGR, 
					 minConversionRate=minConversionRate, maxNAfraction=maxNAfraction, 
					 bedFilePrefix=NULL, path=path, convRatePlots=TRUE, nThreads=nThreads)


#saveRDS(allSampleMats,paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,"_",
		samples[taskID],".rds"))



if(dataType=="gw"){
    ###############################################################################################
    #####################    Extract single read data centered on high confidence TSS   ###########
    ###############################################################################################

    ###################################################
    # get single read matrices for all Cs
    ###################################################

    # to get individual methylation matrices for an TSS
    # make a list (by sample) of lists of matrices (by TSS)
    # e.g
    # [1] sample1
    #     [1] TSS1 matrix of reads x Cpositions
    #     [2] TSS2 matrix of reads x Cpositions
    # [2] sample2
    #     [1] TSS1 matrix of reads x Cpositions
    #     [2] TSS2 matrix of reads x Cpositions

    #######
    # Extract all C positions within each amplicon
    #######

    winSize=300
    regionType="hcTSS"
    tssWin<-highConfTSS
    mcols(tssWin)$TSS<-start(tssWin)
    tssWin<-resize(tssWin,width=winSize,fix="center")

    allSampleMats<-getSingleMoleculeMatrices(sampleTable=fileList[taskID,], genomeFile=genomeFile, 
					     regionGRs=tssWin, regionType=regionType, 
					     genomeMotifGR=genomeMotifGR, 
					     minConversionRate=minConversionRate, maxNAfraction=maxNAfraction,
					     bedFilePrefix=NULL, path=path, convRatePlots=FALSE)


    #saveRDS(allSampleMats,paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,"_",
    #						samples[taskID],".rds"))




    ###############################################################################################
    #####################    Extract single read data centered on less confidence TSS   ###########
    ###############################################################################################

    # to get individual methylation matrices for an TSS
    # make a list (by sample) of lists of matrices (by TSS)
    # e.g
    # [1] sample1
    #     [1] TSS1 matrix of reads x Cpositions
    #     [2] TSS2 matrix of reads x Cpositions
    # [2] sample2
    #     [1] TSS1 matrix of reads x Cpositions
    #     [2] TSS2 matrix of reads x Cpositions

    #######
    # Extract all C positions within each amplicon (getCmethMatrix function)
    #######

    winSize=300
    regionType="lcTSS"
    tssWin<-lessConfTSS
    mcols(tssWin)$TSS<-start(tssWin)
    tssWin<-resize(tssWin,width=winSize,fix="center")

    allSampleMats<-getSingleMoleculeMatrices(sampleTable=fileList[taskID,], genomeFile=genomeFile, 
					     regionGRs=tssWin, regionType=regionType, 
					     genomeMotifGR=genomeMotifGR, minConversionRate=minConversionRate,
					     maxNAfraction=maxNAfraction, bedFilePrefix=NULL, path=path, 
					     convRatePlots=FALSE)


    #saveRDS(allSampleMats,paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,
    #					"_",samples[taskID],".rds"))

} # end of gw data only processing

