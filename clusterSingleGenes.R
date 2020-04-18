#########################
## EM clustering of single gene profiles
# date: 2020-04-18
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
suppressMessages(library(foreach))
suppressMessages(library(doParallel))

#####################
# Custom functions
#####################

#source('./R/dSMF_auxiliaryFunctions.R')
#source('./R/dSMF_bigwig.R')
source('./R/variableSettings.R')


##############
# Variables
##############
args = commandArgs(trailingOnly=TRUE)
taskId=as.numeric(args[1])
maxTasks=as.numeric(args[2])
nThreads=as.numeric(args[3])
print(paste0("taskId:",taskId," maxTasks:",maxTasks," nThreads:",nThreads))

options("scipen"=8, "digits"=4)
options(device=pdf)


# standard variables that should probably not be changed
minConversionRate=0.8
maxNAfraction=0.2



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
outPath=paste0(path,"/EMres")
setSeed=FALSE

if (!dir.exists(outPath)){
  dir.create(outPath)
}

#split table indecies into nTasks number of groups
taskSubList<-split(1:nrow(matTable),sort(1:nrow(matTable)%%maxTasks))

set.seed(200413)
#for (i in 1:nrow(matTable)) {
for (i in taskSubList[[taskId]]){

  ################
  # process matrix
  ################
  # read in first matrix
  print(paste("matrix number: ", i))
  regionName=matTable$region[i]
  sampleName=matTable$sample[i]
  outFileBase=paste(sampleName, regionName, sep="_")
  print(paste("Clust", outFileBase))
  dataMatrix<-readRDS(matTable[i, "filename"])
  # remove NAs
  print("dimensions of dataMatrix: ")
  print(dim(dataMatrix))
  dataMatrix<-removeAllNArows(dataMatrix)
  print("dimensions of dataMatrix after NA row removal:  ")
  print(dim(dataMatrix))
  
  allClassMeans<-tryCatch( {
  	print("running EM for a range of class sizes")
	allClasssMeans<-runEMrangeClassNum(dataMatrix, k_range, convergenceError, 
    			maxIterations, repeats=numRepeats, outPath=outPath, xRange=xRange, 
			outFileBase=paste(sampleName, regionName, sep="_"),
			doIndividualPlots=FALSE)
  },
   	error=function(e){"Matrix not valid"}
  )

  if(is.list(dim(allClassMeans))){
     	saveRDS(allClassMeans,paste0(outPath,"/allClassMeans_",outFileBase,".rds"))
  } else {
     	print(allClassMeans)
  }
	
  clustMetrics<-tryCatch( {
	print("plotting clustering metrics for a range of class sizes")
	plotClusteringMetrics(dataMatrix, k_range, maxB, convergenceError,
		maxIterations, outPath, outFileBase, nThreads, setSeed)
  },
  error=function(e){"Matrix not valid"}
  )
  print(clustMetrics)
}








