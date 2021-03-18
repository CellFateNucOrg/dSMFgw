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
#nThreads=1
#maxTasks=1
#taskId=1
print(paste0("taskId:",taskId," maxTasks:",maxTasks," nThreads:",nThreads))

options("scipen"=8, "digits"=4)
options(device=pdf)
options(warn=-1)

# standard variables that should probably not be changed
minConversionRate=0.8
maxNAfraction=0.2

# parameters for multigeneMatrix generation
minReads<-50 # number of reads sampled from each gene
binSize<-30 # size of window used in multigene matrix

# parameters for multigene clustering
#k_range = 4:10      # Number of classes to be found
#maxIterations = 100 # number of iterations of EM clustering to perform if it does not converge
#convergenceError = 10e-6
numRepeats=10 # number of repeats of clustering each matrix (to account for fraction of methylation)
xRange=c(-250,250)
#maxB=100 # Number of randomised matrices to generate
#setSeed=FALSE # refers to seed within the functions, only necessary when doing automatic testing

#distMetric=list(name="cosineDist", valNA=0, rescale=F)
#distMetric=list(name="cosineDist",rescale=T)
#distMetric=list(name="canberraDist",rescale=T)
distMetric=list(name="euclidean", valNA=0.5, rescale=T)
#distMetric=list(name="correlationDist", valNA=0.5, rescale=T)
#distMetric=list(name="mutualInformation",valNA=0.5,rescale=F)

rndSeed=267413
#rndSeed=181965
set.seed(rndSeed)

outPath=paste0(path,"/EMm_",substr(distMetric$name,1,3),"_",minReads,"r_w",binSize,"_NA",maxNAfraction*100,"_",rndSeed,"_", ifelse(distMetric$rescale,"m1To1","0To1"),"_vNA",distMetric$valNA)


print(paste("outPath:", outPath))
print(paste("winSize:", binSize))
print(paste("maxNAfraction:", maxNAfraction))
print(paste("minReads:", minReads))
print(paste("distMetric:",distMetric$name))
print(paste("randSeed:",rndSeed))



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
print(paste("sampleName:",sampleNames))

matTable<-matTable[matTable$sample==sampleName,]

###################################################
# create multi gene matrix
###################################################

multiGeneMat<-NULL
genesIncluded<-0
#make multigene matrix from only one sample at a time
for(i in 1:nrow(matTable)){
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
    winMatrix<-prepareWindows(fullMatrix, binSize=binSize)
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
# classify by known classes
################


classes<-as.matrix(readRDS(paste0(path,"/classMeans_w30_euc.rds")))

numClasses=nrow(classes)


if (!dir.exists(outPath)){
  dir.create(outPath)
}

regionName="multiGene"
outFileBase=paste(sampleName, regionName, sep="_")
print(paste("Clustering", outFileBase))
dataMatrix<-multiGeneMat
dim(dataMatrix)

#pdf(file=paste0("hist_w",binSize,"_NA",maxNAfraction,"_genes",genesIncluded,".pdf"), width=5.5, height=4)
#hist(rowSums(is.na(dataMatrix)),main=paste0("binSize: ",binSize,", maxNA: ",maxNAfraction,
#                                              ", genesIncluded:",genesIncluded),
#     xlab="number of NA windows per molecule",xlim=c(0,500), col="grey80")
#abline(v=mean(rowSums(is.na(dataMatrix))),lty=2,col="grey20")
#dev.off()


################
# process matrix
################

dataOrderedByClass<-tryCatch( {
	print("running EM for a range of class sizes")
	runClassLikelihoodRpts(dataMatrix, classes,  numRepeats=numRepeats, 
			outPath=outPath, xRange=xRange, 
			outFileBase=outFileBase, figFormat="pdf", 
			distMetric=distMetric)
},
 	error=function(e){"Matrix not valid"}
)

	
pcaPlots<-tryCatch( {
 	print("plotting PCA of clusters")
	plotPCAofMatrixClasses(c(numClasses), outPath, outFileBase)
  },
  error=function(e){"Matrix not valid"}
)
if(length(pcaPlots)==1) {
	print(pcaPlots)
}


umapPlots<-tryCatch(
  {
    print("plotting UMAP of clusters")
    plotUMAPofMatrixClasses(c(numClasses), outPath, outFileBase)
  },
  error=function(e){"Matrix not valid"}
)
if(length(umapPlots)==1) {
  print(umapPlots)
}

print("plotting classes per gene")
plotGenesPerClass(c(numClasses), outPath, outFileBase)





