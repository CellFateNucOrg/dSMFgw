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
options(warn=-1)

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
head(matTable)
matTable<-matTable[!is.na(matTable$filename),]
sampleName=unique(matTable$sample)[taskId]
print(sampleName)

###################################################
# create multi gene matrix
###################################################
multiGeneMat<-NULL
genesIncluded<-0
minReads<-200
#make multigene matrix from only one sample at a time
for(i in 1:nrow(matTable[matTable$sample==sampleName,])){
  regionName=matTable$region[i]
  outFileBase=paste(sampleName, regionName, sep="_")
  dataMatrix<-readRDS(matTable$filename[i])
  # remove rows with too many NAs
  dataMatrix<-removeNArows(dataMatrix, maxNAfraction=maxNAfraction) 
 
  subMatrix<-selectReadsFromMatrix(dataMatrix,minReads=minReads,
                                 addToReadName=outFileBase,
                                 preferBest=T)
  if(!is.null(subMatrix)){
    fullMatrix<-getFullMatrix(subMatrix)
    winMatrix<-prepareWindows(fullMatrix)
    genesIncluded<-genesIncluded+1
    if(is.null(multiGeneMat)){
      multiGeneMat<-winMatrix
    } else {
      multiGeneMat<-rbind(multiGeneMat,winMatrix)
    }
  }
}
print(paste(genesIncluded,"genes included in the multi gene matrix"))

#multiGeneMat<-rescale_minus1To1(multiGeneMat)
#multiGeneMat<-rescale_0To1(multiGeneMat)



################
# learn classes for multiple genes
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
outPath=paste0(path,"/EMres_cosine_50reads_m1to1_rpt2")
setSeed=FALSE
distMetric=list(name="cosineDist",rescale=T)


if (!dir.exists(outPath)){
  dir.create(outPath)
}

#split table indecies into nTasks number of groups
#taskSubList<-split(1:nrow(matTable),sort(1:nrow(matTable)%%maxTasks))

regionName="multiGene"
outFileBase=paste(sampleName, regionName, sep="_")
print(paste("Clustering", outFileBase))
dataMatrix<-multiGeneMat
dim(dataMatrix)

#set.seed(200413) #rpt1
set.seed(210820) #rpt2
################
# process matrix
################

allClassMeans<-tryCatch( {
	print("running EM for a range of class sizes")
	runEMrangeClassNum(dataMatrix, k_range, convergenceError, 
  			maxIterations, EMrepeats=numRepeats, outPath=outPath, xRange=xRange, 
			outFileBase=paste(sampleName, regionName, sep="_"),
			doIndividualPlots=FALSE, distMetric=distMetric)
},
 	error=function(e){"Matrix not valid"}
)

if(is.list(allClassMeans)){
   	saveRDS(allClassMeans,paste0(outPath,"/allClassMeans_",outFileBase,".rds"))
} else {
   	print(allClassMeans)
}
	
clustMetrics<-tryCatch( {
	print("plotting clustering metrics for a range of class sizes")
	plotClusteringMetrics(dataMatrix, k_range, maxB, convergenceError,
		maxIterations, outPath, outFileBase, EMrep=NULL, nThreads=nThreads, 
		setSeed=setSeed, distMetric=distMetric)
  },
  error=function(e){"Matrix not valid"}
)
if(length(clustMetrics)==1){
	print(clustMetrics)
}

pcaPlots<-tryCatch( {
 	print("plotting PCA of clusters")
	plotPCAofMatrixClasses(k_range, outPath, outFileBase)
  },
  error=function(e){"Matrix not valid"}
)
if(length(pcaPlots)==1) {
	print(pcaPlots)
}

umapPlots<-tryCatch(
  {
    print("plotting UMAP of clusters")
    plotUMAPofMatrixClasses(k_range, outPath, outFileBase)
  },
  error=function(e){"Matrix not valid"}
)
if(length(umapPlots)==1) {
  print(umapPlots)
}

print("plotting classes per gene")
plotGenesPerClass(k_range, outPath, outFileBase)





