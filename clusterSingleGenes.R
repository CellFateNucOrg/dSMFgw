## EM clustering of single gene profiles
####### Align bisulfite sequencing of genomeWide dSMF  to C elegans genome
# date: 2010-01-10
# author: Jennnifer Semple
#
# input required:
# single gene dSMF matrices saved in RDS files and a csv table with list of all matrices


###############################################################################################
#####################             SETUP                                 #######################
###############################################################################################


#############
# Libraries
#############


suppressMessages(library(methMatrix))
suppressMessages(library(EMclassifieR))


#####################
# Custom functions
#####################

#source('./R/dSMF_auxiliaryFunctions.R')
#source('./R/dSMF_bigwig.R')
source('./R/variableSettings.R')



##############
# Variables
##############
options("scipen"=8, "digits"=4)
options(device=pdf)


# standard variables that should probably not be changed
minConversionRate=0.8
maxNAfraction=0.2


# assume motif file and bed files are located in the same directory as the genomeFile and their name is
# derived from it (as created by the getGenomeMotifs.R script
#motifFile=gsub("\\.fa","\\.CGGC_motifs.RDS",genomeFile)
#genomeMotifGR<-readRDS(motifFile)
#bedFilePrefix=gsub("\\.fa","",genomeFile)


# assume sample table is stored in ./txt/bwameth_Aligned.txt
#fileList<-read.delim("./txt/bwameth_Aligned.txt",stringsAsFactors=F)
#fileList<-fileList[ordered(fileList$SampleName),]
#samples<-fileList$SampleName



###################################################
# load single read data centered on amplicon TSS
###################################################

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


matTable<-readRDS(paste0(path,"/rds/allSampleRelCoordMats_",regionType,"_",seqDate,"_",expName,".rds"))

#TSSrelCoord<-convertGRtoRelCoord(ampTSS,1,anchorPoint="middle")
#tssWinRelCoord<-convertGRtoRelCoord(ampTSS,winSize,anchorPoint="middle")

#samples<-unique(allSampleRelCoordMats$sample)

################
# learn classes for single gene
################

################
# parameters
################
k_range = 2:8      # Number of classes to be found
maxIterations = 100 # number of iterations of EM clustering to perform if it does not converge
convergenceError = 10e-6
numRepeats=10 # number of repeats of clustering each matrix (to account for fraction of methylation)
xRange=c(-250,250)
maxB=100 # Number of randomised matrices to generate
outPath="./resultsEM"

for (i in 1:nrow(matTable)) {
 
  ################
  # process matrix
  ################
  set.seed(1)
  # read in first matrix
  regionName=matTable$region[i]
  sampleName=matTable$sample[i]
  outFileBase=paste(sampleName, regionName, sep="_")
  print(paste("Clust", outFileBase))
  dataMatrix<-readRDS(matTable[i, "filename"])
  # remove NAs
  dim(dataMatrix)
  dataMatrix<-removeAllNArows(dataMatrix)
  dim(dataMatrix)
  
  tryCatch( 
    {
       runEMrangeClassNum(dataMatrix, k_range, convergenceError, maxIterations,
                     repeats=numRepeats, outPath=outPath, xRange=xRange, 
                     outFileBase=paste(sampleName, regionName, sep="_"),
                     doIndividualPlots=FALSE)
  
       plotClusteringMetrics(dataMatrix, k_range, maxB, convergenceError,
                       maxIterations, outPath, outFileBase)
    },
    error=function(e){"Matrix not valid"}
  ) 
}








