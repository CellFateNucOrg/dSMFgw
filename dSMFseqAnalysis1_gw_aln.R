## mapping with QuasR
####### Align bisulfite sequencing of genomeWide dSMF  to C elegans genome
# date: 2018-08-22
# author: Jennnifer Semple
# Scripts adapted from Drosophila scripts kindly provided by Arnaud Krebs
#
# input required:
# fastq files from sequencing in ../rawData directory
# sampleList.txt with FileName1 FileName2 SampleName headings in rawData directory
# suggested sampleName structure:
# strainName_dSMFvXXXgw where XXX is the unique sample number, and strainName is the biological
# difference we wish to compare.
#
# output files include:
# alignment .bam files (contain a unique automatically generated ref num) in aln/ directory
# QC plot: ./plots/QC_QualityTrimmed_"date".pdf in plots/ directory
# list of aligned files to recall the project in other scripts: QuasR_Aligned.txt in current dir


###############################################################################################
#####################             SETUP                                 #######################
###############################################################################################


#############
# Libraries
#############


library(QuasR)
library(rtracklayer)
library("RColorBrewer")
library(GGally)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggpubr)


# collect citations for packages used
packageBib<-toBibtex(c(citation("QuasR"),
                       citation("rtracklayer"),
                       citation("RColorBrewer"),
                       citation("GGally"),
                       citation("gridExtra"),
                       citation("dplyr"),
                       citation("ggpubr")))


#####################
# Custom functions
#####################

source('./R/dSMF_auxiliaryFunctions.R')
source('./R/dSMF_bigwig.R')
source('./R/variableSettings.R')



##############
# Variables
##############

args = commandArgs(trailingOnly=TRUE)
genomeFile=args[1]
# see ./R/variableSettings.R file. Some of these variables need to be adjusted before running
# the script. the variableSettings_example.R file downloaded from the repo should be correctly
# filled out and saved without the _example extension.

#################
# Directories
################

#setup directory structure this to desired location of your alignments
if (!dir.exists(paste0(path,"/plots"))) {
  dir.create(paste0(path,"/plots"))
}
if (!dir.exists(paste0(path,"/txt"))) {
  dir.create(paste0(path,"/txt"))
}
if (!dir.exists(paste0(path,"/methylation_calls"))){
  dir.create(paste0(path,"/methylation_calls"))
}
if (!dir.exists(paste0(path,"/csv"))){  # for various summary tables of information
  dir.create(paste0(path,"/csv"))
}
if (!dir.exists(paste0(path,"/bigwig"))){
  dir.create(paste0(path,"/bigwig"))
}


###############################################################################################
#####################               CALL METHYLATION on all Cs           ######################
###############################################################################################

cluObj=makeCluster(threadNum)
# point the QuasR project at BAM files from now on
dSMFproj=qAlign(sampleFile=paste0(path,'/txt/QuasR_Aligned.txt'),
                genome=genomeFile,
                paired="fr",
                bisulfite="undir",
                projectName=projectName,
                clObj=cluObj)

# read sample names from dSMFproj
samples<-dSMFproj@alignments$SampleName

########################################################
# Call all methylted Cs in data (context independant)
########################################################

# call methylation on allCs within genome
meth_gr <- qMeth(dSMFproj,mode="allC",clObj=cluObj)
#16994

# qMeth command produces a genomic ranges object with all C positions in the genome in the query regions.
# The mcols contain two columns for each sample with total (_T) counts of reads at
# that position, or methylated counts at that positions (_M). Also includes positions
# without coverage.

## save as rds for future access
saveRDS(meth_gr,paste0(path,'/methylation_calls/dSMFproj_allCs_gw_counts.rds'))

########################################################
# plots of all C methylation coverage
########################################################

#meth_gr<-readRDS(paste0(path,'/methylation_calls/dSMFproj_allCs_gw_counts.rds'))

# Make some histograms of overall C methylation in the genome
pdf(paste0(path,"/plots/hist_allC_gw_coverage.pdf"),width=8,height=11,paper="a4")
par(mfrow=c(2,1))
for (s in samples) {
  columnTotal<-paste0(s,"_T")
  columnMeth<-paste0(s,"_M")
  hist(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnTotal],breaks=10000,xlim=c(1,200),
       main=paste(s, ": total coverage"),xlab="read counts")
  hist(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnMeth],breaks=10000,xlim=c(1,200),
       main=paste(s, ": counts of methylated Cs"),xlab="read counts")
}
dev.off()


########################################################
# tables of all C methylation coverage
########################################################

# make a table summarising total and methylated C coverage.
# Cs0freq is number of cytosines with no coverage. rest of data is based on Cs with coverage (meaningful for amplicon data)
Ccoverage<-data.frame(sampleNames=samples,Cs0freq=0,meanCoverage=0,medianCoverage=0,stdevCoverage=0,
                      meanMethCount=0,medianMethCount=0,stdevMethCount=0)
#excluding 0s
for (s in samples) {
  columnTotal<-paste0(s,"_T")
  columnMeth<-paste0(s,"_M")
  Ccoverage[Ccoverage$sampleNames==s,"Cs0freq"]<-sum(mcols(meth_gr)[,columnTotal]==0)/length(mcols(meth_gr)[,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"meanCoverage"]<-mean(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"medianCoverage"]<-median(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"stdevCoverage"]<-sd(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"meanMethCount"]<-mean(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnMeth])
  Ccoverage[Ccoverage$sampleNames==s,"medianMethCount"]<-median(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnMeth])
  Ccoverage[Ccoverage$sampleNames==s,"stdevMethCount"]<-sd(mcols(meth_gr)[mcols(meth_gr)[,columnTotal]!=0,columnMeth])
}

write.csv(Ccoverage, file=paste0(path,"/csv/CytosineCoverage_gw_excluding0.csv"))

# now do the same but include all 0s in the stats (meaningful for genome wide data)
# including 0s
Ccoverage<-data.frame(sampleNames=samples,Cs0freq=0,meanCoverage=0,medianCoverage=0,stdevCoverage=0,
                      meanMethCount=0,medianMethCount=0,stdevMethCount=0)
for (s in samples) {
  columnTotal<-paste0(s,"_T")
  columnMeth<-paste0(s,"_M")
  Ccoverage[Ccoverage$sampleNames==s,"Cs0freq"]<-sum(mcols(meth_gr)[,columnTotal]==0)/length(mcols(meth_gr)[,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"meanCoverage"]<-mean(mcols(meth_gr)[,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"medianCoverage"]<-median(mcols(meth_gr)[,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"stdevCoverage"]<-sd(mcols(meth_gr)[,columnTotal])
  Ccoverage[Ccoverage$sampleNames==s,"meanMethCount"]<-mean(mcols(meth_gr)[,columnMeth])
  Ccoverage[Ccoverage$sampleNames==s,"medianMethCount"]<-median(mcols(meth_gr)[,columnMeth])
  Ccoverage[Ccoverage$sampleNames==s,"stdevMethCount"]<-sd(mcols(meth_gr)[,columnMeth])
}

write.csv(Ccoverage, file=paste0(path,"/csv/CytosineCoverage_gw_including0.csv"))



###############################################################################################
#####################         CALL METHYLATION in CpG/GpC context        ######################
###############################################################################################

# find sequence context of Cs using call_context_methylation
c0=10 # minimal read coverage for a C (low coverage discarded)
methFreq_grl<-call_context_methylation(meth_gr,c0,genome=genome)

# call_context_methylation returns list of two matrices, "CG" and "GC" in which
# sample names (_M) columns contain fraction methylation and type column contains C context.
# There has been no filtering for rows with no coverage.
saveRDS(methFreq_grl,paste0(path,"/methylation_calls/dSMFproj_allCG-GC_gw_methFreq_grl.rds"))


##################################################################
# histograms of all C methylation frequency in different contexts
##################################################################

# plots histograms of frequency of methylation of CG and GC
pdf(paste0(path,"/plots/hist_CG-GC_gw_freq.pdf"),width=8,height=11,paper="a4")
par(mfrow=c(2,1))
for (s in samples){
  columnMeth<-paste0(s,"_M")
  hist(unlist(mcols(methFreq_grl[["CG"]])[,columnMeth]),breaks=100,
       main=paste(s,"CmG frequency"),xlab="fraction methylated")
  hist(unlist(mcols(methFreq_grl[["GC"]])[,columnMeth]),breaks=100,
       main=paste(s,"GCm frequency"),xlab="fraction methylated")
}
dev.off()

###################################################
# plot number of CG/GCs with meth data per amplicon
###################################################

#methFreq_grl<-readRDS(paste0(path,"/methylation_calls/dSMFproj_allCG-GC_gw_methFreq_grl.rds"))

# merge and sort the CG and GC matrices
methFreq_gr=sort(unlist(GRangesList(methFreq_grl)))

# save merged and sorted methylation gr
saveRDS(methFreq_gr,paste0(path,"/methylation_calls/dSMFproj_allCG-GC_gw_methFreq_gr.rds"))

#methFreq_gr<-readRDS(paste0(path,"/methylation_calls/dSMFproj_allCG-GC_gw_methFreq_gr.rds"))


pdf(paste0(path,"/plots/ampliconCoverageWithMethdata.pdf"),width=8,height=11,paper="a4")
par(mfrow=c(2,1))

plotMethCoveragePerAmplicon<-function(allCcount,methCcount,sites="",samples="") {
  barplot(t(cbind(allCcount-methCcount,methCcount)),beside=F,xlab="amplicons",col=c("aliceblue","azure4"),
          main=paste0("Number of ",sites," with meth data per amplicon ",samples))
  legend("topleft",fill=c("aliceblue","azure4"), bg="white",
         legend=c(paste0(sites,c(" without meth data", " with meth data"))))
  barplot(t(cbind((allCcount-methCcount)/allCcount,methCcount/allCcount)),beside=F,xlab="amplicons",col=c("aliceblue","azure4"),
          main=paste0("Freq of ",sites," with meth data per amplicon ",samples))
  legend("topleft",fill=c("aliceblue","azure4"), bg="white",
         legend=c(paste0(sites,c(" without meth data", " with meth data"))))
}

# filter out rows that are only NAs to see how many amplicons have data
naix=rowSums(is.na(mcols(methFreq_gr)))==dim(mcols(methFreq_gr))[2]-1
allCcount<-countOverlaps(amplicons, methFreq_gr)
methCcount<-countOverlaps(amplicons, methFreq_gr[!naix]) ##
plotMethCoveragePerAmplicon(allCcount,methCcount,sites="CG/GCs",samples="(at least 1 sample)")

# same but only looking at Cs which have data for all samples
methCcount<-countOverlaps(amplicons, methFreq_gr[complete.cases(mcols(methFreq_gr))])
plotMethCoveragePerAmplicon(allCcount,methCcount,sites="CG/GCs",samples="(in all samples)")

### do it separately for CG and GC - only looking at samples with data for all samples
# CG sites with data in all samples
allCcount<-countOverlaps(amplicons, methFreq_grl[["CG"]])
methCcount<-countOverlaps(amplicons, methFreq_grl[["CG"]][complete.cases(mcols(methFreq_grl[["CG"]]))])
plotMethCoveragePerAmplicon(allCcount,methCcount,sites="CGs",samples="(in all samples)")

# GC sites with data in all samples
allCcount<-countOverlaps(amplicons, methFreq_grl[["GC"]])
methCcount<-countOverlaps(amplicons, methFreq_grl[["GC"]][complete.cases(mcols(methFreq_grl[["GC"]]))])
plotMethCoveragePerAmplicon(allCcount,methCcount,sites="GCs",samples="(in all samples)")

dev.off()


###################################################
# plot correlation between duplicates
###################################################

pdf(paste0(path,"/plots/CorrelationBetweenSamples.pdf"),width=8,height=11,paper="a4")

#avoid plotting too many points - randomly choose 10,000 if more than that
set.seed(1)
if (length(methFreq_gr)>1e4) {
  sampleCs<-sample(1:length(methFreq_gr),1e4)
} else {
  sampleCs<-1:length(methFreq_gr)
}

lapply(testGroup,function(g) {
    bioReplicates<-c(grep(g,names(mcols(methFreq_gr))))
    p<-ggpairs(as.data.frame(1-as.matrix(mcols(methFreq_gr[sampleCs, bioReplicates]))),
             title=paste0("Correlation between biological replicates: ",g))
    print(p)
})


# correlation plot for all samples
typeColumn=which(names(mcols(methFreq_gr))=="type")
ggpairs(as.data.frame(1-as.matrix(mcols(methFreq_gr[sampleCs,-typeColumn]))),
        title="Correlation between all samples")

# simple colour box correlation plot
ggcorr(as.data.frame(1-as.matrix(mcols(methFreq_gr[,-typeColumn]))),
       low = "steelblue", mid = "white", high = "darkred", label = T,
       label_size = 8, label_color = "white",label_round=2, name="rho")
dev.off()


###################################################
# make bigwig tracks of data
###################################################

# reading in combined GC GC list
#methFreq_gr<-readRDS(paste0(path,"/methylation_calls/dSMFproj_allCG-GC_gw_methFreq_gr.rds"))

w=10 #winSize for smoothing

if (GenomeInfoDb::seqlevelsStyle(methFreq_gr)!="UCSC") {
  methFreq_gr_u<-wbToUcscGR(methFreq_gr)
} else {
  methFreq_gr_u<-methFreq_gr
}

dataCols<-names(mcols(methFreq_gr_u))[grep("_M$",names(mcols(methFreq_gr_u)))]

dsmf_all<-methToDSMFgr(methFreq_gr_u,dataCols)

# export raw data with no smoothing
grToBw(dsmf_all,dataCols,bwPath=paste0(path,"/bigwig"),
       filenamePrefix="rawDSMF_",
       urlPrefixForUCSC="http://www.meister.izb.unibe.ch/ucsc/")

#smoothe data
smDSMF<-smootheGRdata(dsmf_all,dataCols,winSize=w,winStep=1)
#export smoothed tracks
grToBw(smDSMF,dataCols,bwPath=paste0(path,"/bigwig"),
       filenamePrefix=paste0("w",w,"smDSMF_"),
       urlPrefixForUCSC="http://www.meister.izb.unibe.ch/ucsc/")




###############################################################################################
#####################    Extract single read data centered on amplicon TSS      ###############
###############################################################################################

###################################################
# get single read matrices for all Cs
###################################################

# point the QuasR project at BAM files from now on
dSMFproj=qAlign(sampleFile=paste0(path,'/txt/QuasR_Aligned.txt'),
                genome=genomeFile,
                paired="fr",
                bisulfite="dir",
                projectName=projectName,
                clObj=cluObj)

# read sample names from dSMFproj
samples<-dSMFproj@alignments$SampleName


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
winSize<-300
TSS<-ampTSS
tssWin<-TSS
mcols(tssWin)$TSS<-start(TSS)
tssWin<-resize(tssWin,width=winSize,fix="center")
names(mcols(tssWin))[1]<-"ID"
allCmats=list()
allSampleCmats=list()
for (currentSample in samples) {
  for (i in seq_along(tssWin)) {
    print(i)
    mat<-getCmethMatrix(dSMFproj,tssWin[i],currentSample)
    allCmats[[tssWin[i]$ID]]<-mat
  }
  allSampleCmats[[currentSample]]<-allCmats
}
saveRDS(allSampleCmats,paste0(path,"/methylation_calls/allSampleCmats_TSS_amp.rds"))



###################################################
# get separate GC and CG (and GCGC) matrices within each TSS region (getGCmatrix function)
###################################################

#allSampleCmats<-readRDS(paste0(path,"/methylation_calls/allSampleCmats_TSS_amp.rds"))

allSampleGCmats<-list()
for (i in seq_along(samples)) {
  allSampleGCmats[[samples[i]]]<-getGCmatrix1(matList=allSampleCmats[[i]], ampliconGR=tssWin, genome=genome,
                                              conv.rate=80, sampleName=names(allSampleCmats)[i],destrand=F,plotData=F) # for Amplicon data use destrand=F !!!
}
saveRDS(allSampleGCmats,paste0(path,"/methylation_calls/allSampleCGmatGCmat_TSS_amp.rds"))


###################################################
# merge GC and CG (and GCGC) matrices within each TSS region (mergeGC_CGmats) to get a single matrix
###################################################

#allSampleGCmats<-readRDS(paste0(path,"/methylation_calls/allSampleCGmatGCmat_TSS_amp.rds"))


allSampleMergedMats<-list()
for (i in seq_along(samples)) {
  allSampleMergedMats[[samples[i]]]<-mergeGC_CGmats(matList=allSampleGCmats[[samples[i]]])
}
saveRDS(allSampleMergedMats,paste0(path,"/methylation_calls/allSampleMergedMats_TSS_amp.rds"))


###################################################
# convert merged matrices to relative coordinates
###################################################

#allSampleMergedMats<-readRDS(paste0(path,"/methylation_calls/allSampleMergedMats_TSS_amp.rds"))

allSampleRelCoordMats<-list()
for (i in seq_along(samples)) {
  allSampleRelCoordMats[[samples[i]]]<-getRelativeCoordMats(matList=allSampleMergedMats[[samples[i]]],
                                                            grs=tssWin,anchorCoord=winSize/2)
}
saveRDS(allSampleRelCoordMats,paste0(path,"/methylation_calls/allSampleRelCoordMats_TSS_amp.rds"))


#######
# plot single read matrices for all regions
#######


#allSampleRelCoordMats<-readRDS(paste0(path,"/methylation_calls/allSampleRelCoordMats_TSS_amp.rds"))

TSSrelCoord<-TSS
strand(TSSrelCoord)<-"*"
names(mcols(TSSrelCoord))[1]<-"ID"
start(TSSrelCoord)<-0
end(TSSrelCoord)<-0

tssWinRelCoord<-tssWin
start(tssWinRelCoord)<- -winSize/2
end(tssWinRelCoord)<- winSize/2


allAmp2plot<-unique(unlist(lapply(allSampleRelCoordMats,function(x){names(x)})))

# plot single molecule matrices on their own
for (i in allAmp2plot) {
  if (!dir.exists(paste0(path,"/plots/singleMoleculePlots"))){  # for alignments
    dir.create(paste0(path,"/plots/singleMoleculePlots"))
  }
  plotList=list()
  print(paste("plotting", i))
  for (j in seq_along(samples)) {
    mat<-allSampleRelCoordMats[[samples[j]]][[i]]
    maxReads=10000
    if (!is.null(dim(mat))) {
      if (dim(mat)[1]>maxReads) { # if matrix contains more than 10000 reads, do a random subsample
        chooseRows<-sample(1:dim(mat)[1],maxReads)
        mat<-mat[chooseRows,]
      }
      p<-plotSingleMoleculesAmp(mat=mat, regionName=i,regionGRs=tssWinRelCoord,
                                featureGRs=TSSrelCoord,title=samples[j])
      plotList[[samples[j]]]<-p
    }
  }
  regionGR<-tssWinRelCoord[match(i,tssWinRelCoord$ID)]
  XorA<-ifelse(seqnames(regionGR)=="chrX","X","A")
  title<-paste0(i, ": ",seqnames(regionGR)," ",strand(regionGR),"ve strand")
  mp<-marrangeGrob(grobs=plotList,nrow=2,ncol=2,top=title)
  ggsave(paste0(path,"/plots/singleMoleculePlots/TSS_",XorA,"_",i,".png"),
         plot=mp, device="png", width=29, height=20, units="cm")
}



# plot single molecule matrices together with average dSMF
for (i in allAmp2plot) {
  if (!dir.exists(paste0(path,"/plots/singleMoleculePlotsWithAvr"))){  # for alignments
    dir.create(paste0(path,"/plots/singleMoleculePlotsWithAvr"))
  }
  plotList=list()
  print(paste("plotting",i))
  for (j in seq_along(samples)) {
    mat<-allSampleRelCoordMats[[samples[j]]][[i]]
    maxReads=10000
    if (!is.null(dim(mat))) {  # only do if there is a matrix
      if (dim(mat)[1]>maxReads) { # if matrix contains more than 10000 reads, do a random subsample
        chooseRows<-sample(1:dim(mat)[1],maxReads)
        mat<-mat[chooseRows,]
      }
      p<-plotSingleMoleculesWithAvr(mat=mat, regionName=i,regionGRs=tssWinRelCoord,
                                    featureGRs=TSSrelCoord,title=samples[j])
      plotList[[samples[j]]]<-p
    }
  }
  regionGR<-tssWinRelCoord[match(i,tssWinRelCoord$ID)]
  XorA<-ifelse(seqnames(regionGR)=="chrX","X","A")
  title<-paste0(i, ": ",seqnames(regionGR)," ",strand(regionGR),"ve strand")
  mp<-marrangeGrob(grobs=plotList,nrow=2,ncol=2,top=title)
  ggsave(paste0(path,"/plots/singleMoleculePlotsWithAvr/TSS_",XorA,"_",i,".png"),
         plot=mp, device="png", width=20, height=29, units="cm")
}



###################################################
# summarise data for Metagene
###################################################

# point the QuasR project at BAM files from now on
dSMFproj=qAlign(sampleFile=paste0(path,'/txt/QuasR_Aligned.txt'),
                genome=genomeFile,
                paired="fr",
                bisulfite="dir",
                projectName=projectName,
                clObj=cluObj)

# read sample names from dSMFproj
samples<-dSMFproj@alignments$SampleName

#allSampleRelCoordMats<-readRDS(paste0(path,"/methylation_calls/allSampleRelCoordMats_TSS_amp.rds"))

tssWin<-TSS
mcols(tssWin)$TSS<-start(TSS)
tssWin<-resize(tssWin,width=winSize,fix="center")
names(mcols(tssWin))[1]<-"ID"
tssWinRelCoord<-tssWin
start(tssWinRelCoord)<--winSize/2
end(tssWinRelCoord)<-winSize/2

rm("allSampleMetaMethFreqDF")
for (i in seq_along(samples)) {
  metaMethFreqDF<-getMetaMethFreq(matList=allSampleRelCoordMats[[samples[i]]],regionGRs=tssWinRelCoord,minReads=10)
  print(samples[i])
  metaMethFreqDF$sample<-samples[i]
  if(exists("allSampleMetaMethFreqDF")) {
    allSampleMetaMethFreqDF<-rbind(allSampleMetaMethFreqDF,metaMethFreqDF)
  } else {
   allSampleMetaMethFreqDF<-metaMethFreqDF
  }
}
# convert position from factor to numeric
allSampleMetaMethFreqDF$position<-as.numeric(as.character(allSampleMetaMethFreqDF$position))
saveRDS(allSampleMetaMethFreqDF,paste0(path,"/methylation_calls/allSampleMetaMethFreqDF_TSS_amp.rds"))

#allSampleMetaMethFreqDF<-readRDS(paste0(path,"/methylation_calls/allSampleMetaMethFreqDF_TSS_amp.rds"))

### plot metagene by sample
# subsample if too many points
if (nrow(allSampleMetaMethFreqDF)>10000) {
  idx<-sample(1:nrow(allSampleMetaMethFreqDF),10000)
 } else {
  idx<-1:nrow(allSampleMetaMethFreqDF)
 }

p1<-ggplot(allSampleMetaMethFreqDF[idx,],aes(x=position,y=1-methFreq)) +
  theme_light(base_size=16) + ylim(0,1) +
  xlab("Position relative to TSS") + ylab("dSMF (1-%methylation)") +
  geom_linerange(aes(x=0, y=NULL, ymin=0, ymax=1),color="steelblue",size=1) +
  geom_point(alpha=0.1) +
  geom_smooth(colour="red",fill="red") +
  facet_wrap(~sample)

p2<-ggplot(allSampleMetaMethFreqDF,aes(x=position,y=1-methFreq,colour=sample)) +
  theme_light(base_size=16) + ylim(0,1) +
  xlab("Position relative to TSS") + ylab("dSMF (1-%methylation)") +
  geom_linerange(aes(x=0, y=NULL, ymin=0, ymax=1),color="steelblue",size=1) +
  geom_smooth(se=FALSE)

ml <- marrangeGrob(list(p1,p2), nrow=1, ncol=1)
ggsave(paste0(path,"/plots/metaGenePlots_TSS_amp.pdf"),plot=ml,device="pdf",
       width=20,height=20,units="cm")




#} # end of brackets for already aligned...



###############################################################################################
#####################    Extract single read data centered on high confidence TSS   ###########
###############################################################################################

###################################################
# get single read matrices for all Cs
###################################################

# point the QuasR project at BAM files from now on
dSMFproj=qAlign(sampleFile=paste0(path,'/txt/QuasR_Aligned.txt'),
                genome=genomeFile,
                paired="fr",
                bisulfite="dir",
                projectName=projectName,
                clObj=cluObj)

# read sample names from dSMFproj
samples<-dSMFproj@alignments$SampleName


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
TSS<-highConfTSS
tssWin<-TSS
mcols(tssWin)$TSS<-start(TSS)
tssWin<-resize(tssWin,width=winSize,fix="center")
#names(mcols(tssWin))[1]<-"ID"
allCmats=list()
allSampleCmats=list()
for (currentSample in samples) {
  for (i in seq_along(tssWin)) {
    print(i)
    mat<-getCmethMatrix(dSMFproj,tssWin[i],currentSample)
    allCmats[[tssWin[i]$ID]]<-mat
  }
  allSampleCmats[[currentSample]]<-allCmats
}
saveRDS(allSampleCmats,paste0(path,"/methylation_calls/allSampleCmats_TSS_hc.rds"))



###################################################
# get separate GC and CG (and GCGC) matrices within each TSS region (getGCmatrix function)
###################################################

#allSampleCmats<-readRDS(paste0(path,"/methylation_calls/allSampleCmats_TSS_hc.rds"))

allSampleGCmats<-list()
for (i in seq_along(samples)) {
  allSampleGCmats[[samples[i]]]<-getGCmatrix1(matList=allSampleCmats[[i]], ampliconGR=tssWin, genome=genome,
                                              conv.rate=80, sampleName=names(allSampleCmats)[i],destrand=F,plotData=F) # for Amplicon data use destrand=F !!!
}
saveRDS(allSampleGCmats,paste0(path,"/methylation_calls/allSampleCGmatGCmat_TSS_hc.rds"))


###################################################
# merge GC and CG (and GCGC) matrices within each TSS region (mergeGC_CGmats) to get a single matrix
###################################################

#allSampleGCmats<-readRDS(paste0(path,"/methylation_calls/allSampleCGmatGCmat_TSS_hc.rds"))


allSampleMergedMats<-list()
for (i in seq_along(samples)) {
  allSampleMergedMats[[samples[i]]]<-mergeGC_CGmats(matList=allSampleGCmats[[samples[i]]])
}
saveRDS(allSampleMergedMats,paste0(path,"/methylation_calls/allSampleMergedMats_TSS_hc.rds"))


###################################################
# convert merged matrices to relative coordinates
###################################################

#allSampleMergedMats<-readRDS(paste0(path,"/methylation_calls/allSampleMergedMats_TSS_hc.rds"))

allSampleRelCoordMats<-list()
for (i in seq_along(samples)) {
  allSampleRelCoordMats[[samples[i]]]<-getRelativeCoordMats(matList=allSampleMergedMats[[samples[i]]],
                                                            grs=tssWin,anchorCoord=winSize/2)
}
saveRDS(allSampleRelCoordMats,paste0(path,"/methylation_calls/allSampleRelCoordMats_TSS_hc.rds"))



###################################################
# summarise data for Metagene
###################################################

# point the QuasR project at BAM files from now on
dSMFproj=qAlign(sampleFile=paste0(path,'/txt/QuasR_Aligned.txt'),
                genome=genomeFile,
                paired="fr",
                bisulfite="dir",
                projectName=projectName,
                clObj=cluObj)

# read sample names from dSMFproj
samples<-dSMFproj@alignments$SampleName

#allSampleRelCoordMats<-readRDS(paste0(path,"/methylation_calls/allSampleRelCoordMats_TSS_hc:.rds"))

tssWin<-TSS
mcols(tssWin)$TSS<-start(TSS)
tssWin<-resize(tssWin,width=winSize,fix="center")
#names(mcols(tssWin))[1]<-"ID"
tssWinRelCoord<-tssWin
start(tssWinRelCoord)<--winSize/2
end(tssWinRelCoord)<-winSize/2

rm("allSampleMetaMethFreqDF")
for (i in seq_along(samples)) {
  metaMethFreqDF<-getMetaMethFreq(matList=allSampleRelCoordMats[[samples[i]]],regionGRs=tssWinRelCoord,minReads=10)
  print(samples[i])
  metaMethFreqDF$sample<-samples[i]
  if(exists("allSampleMetaMethFreqDF")) {
    allSampleMetaMethFreqDF<-rbind(allSampleMetaMethFreqDF,metaMethFreqDF)
  } else {
   allSampleMetaMethFreqDF<-metaMethFreqDF
  }
}
# convert position from factor to numeric
allSampleMetaMethFreqDF$position<-as.numeric(as.character(allSampleMetaMethFreqDF$position))
saveRDS(allSampleMetaMethFreqDF,paste0(path,"/methylation_calls/allSampleMetaMethFreqDF_TSS_hc.rds"))

#allSampleMetaMethFreqDF<-readRDS(paste0(path,"/methylation_calls/allSampleMetaMethFreqDF_TSS_hc.rds"))

### plot metagene by sample
# subsample if too many points
if (nrow(allSampleMetaMethFreqDF)>10000) {
  idx<-sample(1:nrow(allSampleMetaMethFreqDF),10000)
 } else {
  idx<-1:nrow(allSampleMetaMethFreqDF)
 }

p1<-ggplot(allSampleMetaMethFreqDF[idx,],aes(x=position,y=1-methFreq)) +
  theme_light(base_size=16) + ylim(0,1) +
  xlab("Position relative to TSS") + ylab("dSMF (1-%methylation)") +
  geom_linerange(aes(x=0, y=NULL, ymin=0, ymax=1),color="steelblue",size=1) +
  geom_point(alpha=0.1) +
  geom_smooth(colour="red",fill="red") +
  facet_wrap(~sample)

p2<-ggplot(allSampleMetaMethFreqDF,aes(x=position,y=1-methFreq,colour=sample)) +
  theme_light(base_size=16) + ylim(0,1) +
  xlab("Position relative to TSS") + ylab("dSMF (1-%methylation)") +
  geom_linerange(aes(x=0, y=NULL, ymin=0, ymax=1),color="steelblue",size=1) +
  geom_smooth(se=FALSE)

ml <- marrangeGrob(list(p1,p2), nrow=1, ncol=1)
ggsave(paste0(path,"/plots/metaGenePlots_TSS_hc.pdf"),plot=ml,device="pdf",
       width=20,height=20,units="cm")




###############################################################################################
#####################    Extract single read data centered on less confidence TSS   ###########
###############################################################################################

###################################################
# get single read matrices for all Cs
###################################################

# point the QuasR project at BAM files from now on
dSMFproj=qAlign(sampleFile=paste0(path,'/txt/QuasR_Aligned.txt'),
                genome=genomeFile,
                paired="fr",
                bisulfite="dir",
                projectName=projectName,
                clObj=cluObj)

# read sample names from dSMFproj
samples<-dSMFproj@alignments$SampleName


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
TSS<-lessConfTSS
tssWin<-TSS
mcols(tssWin)$TSS<-start(TSS)
tssWin<-resize(tssWin,width=300,fix="center")
#names(mcols(tssWin))[1]<-"ID"
allCmats=list()
allSampleCmats=list()
for (currentSample in samples) {
  for (i in seq_along(tssWin)) {
    print(i)
    mat<-getCmethMatrix(dSMFproj,tssWin[i],currentSample)
    allCmats[[tssWin[i]$ID]]<-mat
  }
  allSampleCmats[[currentSample]]<-allCmats
}
saveRDS(allSampleCmats,paste0(path,"/methylation_calls/allSampleCmats_TSS_lc.rds"))



###################################################
# get separate GC and CG (and GCGC) matrices within each TSS region (getGCmatrix function)
###################################################

#allSampleCmats<-readRDS(paste0(path,"/methylation_calls/allSampleCmats_TSS_lc.rds"))

allSampleGCmats<-list()
for (i in seq_along(samples)) {
  allSampleGCmats[[samples[i]]]<-getGCmatrix1(matList=allSampleCmats[[i]], ampliconGR=tssWin, genome=genome,
                                              conv.rate=80, sampleName=names(allSampleCmats)[i],destrand=F,plotData=F) # for Amplicon data use destrand=F !!!
}
saveRDS(allSampleGCmats,paste0(path,"/methylation_calls/allSampleCGmatGCmat_TSS_lc.rds"))


###################################################
# merge GC and CG (and GCGC) matrices within each TSS region (mergeGC_CGmats) to get a single matrix
###################################################

#allSampleGCmats<-readRDS(paste0(path,"/methylation_calls/allSampleCGmatGCmat_TSS_lc.rds"))


allSampleMergedMats<-list()
for (i in seq_along(samples)) {
  allSampleMergedMats[[samples[i]]]<-mergeGC_CGmats(matList=allSampleGCmats[[samples[i]]])
}
saveRDS(allSampleMergedMats,paste0(path,"/methylation_calls/allSampleMergedMats_TSS_lc.rds"))


###################################################
# convert merged matrices to relative coordinates
###################################################

#allSampleMergedMats<-readRDS(paste0(path,"/methylation_calls/allSampleMergedMats_TSS_lc.rds"))

allSampleRelCoordMats<-list()
for (i in seq_along(samples)) {
  allSampleRelCoordMats[[samples[i]]]<-getRelativeCoordMats(matList=allSampleMergedMats[[samples[i]]],
                                                            grs=tssWin,anchorCoord=winSize/2)
}
saveRDS(allSampleRelCoordMats,paste0(path,"/methylation_calls/allSampleRelCoordMats_TSS_lc.rds"))



###################################################
# summarise data for Metagene
###################################################

# point the QuasR project at BAM files from now on
dSMFproj=qAlign(sampleFile=paste0(path,'/txt/QuasR_Aligned.txt'),
                genome=genomeFile,
                paired="fr",
                bisulfite="dir",
                projectName=projectName,
                clObj=cluObj)

# read sample names from dSMFproj
samples<-dSMFproj@alignments$SampleName

#allSampleRelCoordMats<-readRDS(paste0(path,"/methylation_calls/allSampleRelCoordMats_TSS_lc.rds"))
winSize
tssWin<-TSS
mcols(tssWin)$TSS<-start(TSS)
tssWin<-resize(tssWin,width=winSize,fix="center")
#names(mcols(tssWin))[1]<-"ID"
tssWinRelCoord<-tssWin
start(tssWinRelCoord)<--winSize/2
end(tssWinRelCoord)<-winSize/2

rm("allSampleMetaMethFreqDF")
for (i in seq_along(samples)) {
  metaMethFreqDF<-getMetaMethFreq(matList=allSampleRelCoordMats[[samples[i]]],regionGRs=tssWinRelCoord,minReads=10)
  print(samples[i])
  metaMethFreqDF$sample<-samples[i]
  if(exists("allSampleMetaMethFreqDF")) {
    allSampleMetaMethFreqDF<-rbind(allSampleMetaMethFreqDF,metaMethFreqDF)
  } else {
   allSampleMetaMethFreqDF<-metaMethFreqDF
  }
}
# convert position from factor to numeric
allSampleMetaMethFreqDF$position<-as.numeric(as.character(allSampleMetaMethFreqDF$position))
saveRDS(allSampleMetaMethFreqDF,paste0(path,"/methylation_calls/allSampleMetaMethFreqDF_TSS_lc.rds"))

#allSampleMetaMethFreqDF<-readRDS(paste0(path,"/methylation_calls/allSampleMetaMethFreqDF_TSS_lc.rds"))

### plot metagene by sample
# subsample if too many points
if (nrow(allSampleMetaMethFreqDF)>10000) {
  idx<-sample(1:nrow(allSampleMetaMethFreqDF),10000)
 } else {
  idx<-1:nrow(allSampleMetaMethFreqDF)
 }

p1<-ggplot(allSampleMetaMethFreqDF[idx,],aes(x=position,y=1-methFreq)) +
  theme_light(base_size=16) + ylim(0,1) +
  xlab("Position relative to TSS") + ylab("dSMF (1-%methylation)") +
  geom_linerange(aes(x=0, y=NULL, ymin=0, ymax=1),color="steelblue",size=1) +
  geom_point(alpha=0.1) +
  geom_smooth(colour="red",fill="red") +
  facet_wrap(~sample)

p2<-ggplot(allSampleMetaMethFreqDF,aes(x=position,y=1-methFreq,colour=sample)) +
  theme_light(base_size=16) + ylim(0,1) +
  xlab("Position relative to TSS") + ylab("dSMF (1-%methylation)") +
  geom_linerange(aes(x=0, y=NULL, ymin=0, ymax=1),color="steelblue",size=1) +
  geom_smooth(se=FALSE)

ml <- marrangeGrob(list(p1,p2), nrow=1, ncol=1)
ggsave(paste0(path,"/plots/metaGenePlots_TSS_lc.pdf"),plot=ml,device="pdf",
       width=20,height=20,units="cm")


