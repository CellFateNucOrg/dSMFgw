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
#write.table(fileList,file=paste0(path,"/txt/bwameth_Aligned.txt"), quote=F,
#sep="\t",row.names=F)
samples<-fileList$SampleName


#################
# Directories
################

#setup directory structure this to desired location of your alignments
makeDirs(path,c("plots","methCalls","csv","bigwig","rds"))





if(dataType=="amp"){ # only execute of it is amplicon data
    ###############################################################################################
    #####################         Plot single read data for amplicons     ######################
    ###############################################################################################

    #################
    # Merge matrix lists that were separated by sample
    ################
    regionType="rawAmp"

    allSampleMats<-mergeSampleMats(path, regionType, samples, deleteSplitFiles=F)
    
    saveRDS(allSampleMats, paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,".rds"))

    #######
    # plot single read matrices for all amplicons
    #######

    plotAllMatrices(allSampleMats, samples, regionGRs=amplicons, featureGRs=ampTSS, featureLabel="TSS",
                    regionType=regionType, maxNAfraction=maxNAfraction, withAvr=FALSE, includeInFileName=seqDate,
                    drawArrow=TRUE)

    plotAllMatrices(allSampleMats, samples, regionGRs=amplicons, featureGRs=ampTSS, featureLabel="TSS",
                    regionType=regionType, maxNAfraction=maxNAfraction, withAvr=TRUE, includeInFileName=seqDate,
                    drawArrow=TRUE)


} # end of amplicon only loop




###############################################################################################
#####################    Metagene centered on amplicon TSS                      ###############
###############################################################################################

if(dataType=="amp"){
    winSize=500
} else {
    winSize=300
}

regionType="ampTSS"
tssWin<-ampTSS
mcols(tssWin)$TSS<-start(tssWin)
tssWin<-resize(tssWin,width=winSize,fix="center")

#################
# Merge matrix lists that were separated by sample
################

allSampleMats<-mergeSampleMats(path, regionType, samples, deleteSplitFiles=F)

saveRDS(allSampleMats, paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,".rds"))


###################################################
# convert merged matrices to relative coordinates
###################################################

#allSampleMats<-readRDS(paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,".rds"))

allSampleRelCoordMats<-getRelativeCoordMats(matList=allSampleMats, regionGRs=tssWin, regionType=regionType, anchorCoord=winSize/2)

saveRDS(allSampleRelCoordMats,paste0(path,"/rds/allSampleRelCoordMats_",regionType,"_",seqDate,"_",expName,".rds"))


#######
# plot single read matrices for all regions
#######


#allSampleRelCoordMats<-readRDS(paste0(path,"/rds/allSampleRelCoordMats_",regionType,"_",seqDate,"_",expName,".rds"))

TSSrelCoord<-convertGRtoRelCoord(ampTSS,1,anchorPoint="middle")
tssWinRelCoord<-convertGRtoRelCoord(ampTSS,winSize,anchorPoint="middle")

plotAllMatrices(allSampleRelCoordMats, samples, regionGRs=tssWinRelCoord, featureGRs=TSSrelCoord,
                featureLabel="TSS", regionType=regionType, maxNAfraction=maxNAfraction, withAvr=FALSE,
                includeInFileName=seqDate, drawArrow=FALSE)

plotAllMatrices(allSampleRelCoordMats, samples, regionGRs=tssWinRelCoord, featureGRs=TSSrelCoord,
                featureLabel="TSS", regionType=regionType, maxNAfraction=maxNAfraction, withAvr=TRUE,
                includeInFileName=seqDate, drawArrow=FALSE)


###################################################
# summarise data for Metagene
###################################################

#allSampleRelCoordMats<-readRDS(paste0(path,"/rds/allSampleRelCoordMats_",regionType,"_",seqDate,"_",expName,".rds"))

tssWinRelCoord<-convertGRtoRelCoord(ampTSS,winSize,anchorPoint="middle")

allSampleMetaMethFreqDF<-getAllSampleMetaMethFreq(allSampleRelCoordMats,samples,regionGRs=tssWinRelCoord,
                                                   minReads=10)

saveRDS(allSampleMetaMethFreqDF,paste0(path,"/rds/allSampleMetaMethFreqDF_",regionType,"_",seqDate,"_",expName,".rds"))


### plot metagene by sample
mp<-plotDSMFmetageneDF(metageneDF=allSampleMetaMethFreqDF,maxPoints=10000)

ggsave(paste0(path,"/plots/metaGenePlots_",regionType,"_",seqDate,"_",expName,".pdf"),plot=mp,device="pdf",
       width=20,height=20,units="cm")



if(dataType=="gw"){
    ###############################################################################################
    #####################    Metagene centered on high confidence TSS                   ###########
    ###############################################################################################

    winSize=300    
    regionType="hcTSS"
    tssWin<-highConfTSS
    mcols(tssWin)$TSS<-start(tssWin)
    tssWin<-resize(tssWin,width=winSize,fix="center")

    #################
    # Merge matrix lists that were separated by sample
    ################

    allSampleMats<-mergeSampleMats(path, regionType, samples, deleteSplitFiles=F)

    saveRDS(allSampleMats, paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,".rds"))


    ###################################################
    # convert merged matrices to relative coordinates
    ###################################################

    #allSampleMats<-readRDS(paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,".rds"))


    allSampleRelCoordMats<-getRelativeCoordMats(matList=allSampleMats, regionGRs=tssWin, regionType=regionType,  anchorCoord=winSize/2)

    saveRDS(allSampleRelCoordMats,paste0(path,"/rds/allSampleRelCoordMats_",regionType,"_",seqDate,"_",expName,".rds"))


    ###################################################
    # summarise data for Metagene
    ###################################################

    #allSampleRelCoordMats<-readRDS(paste0(path,"/rds/allSampleRelCoordMats_",regionType,"_",seqDate,"_",expName,".rds"))

    tssWinRelCoord<-convertGRtoRelCoord(highConfTSS,winSize,anchorPoint="middle")

    allSampleMetaMethFreqDF<-getAllSampleMetaMethFreq(allSampleRelCoordMats,samples,regionGRs=tssWinRelCoord,
                                                     minReads=10)

    saveRDS(allSampleMetaMethFreqDF,paste0(path,"/rds/allSampleMetaMethFreqDF_",regionType,"_",seqDate,"_",expName,".rds"))


    ### plot metagene by sample
    mp<-plotDSMFmetageneDF(metageneDF=allSampleMetaMethFreqDF,maxPoints=10000)

    ggsave(paste0(path,"/plots/metaGenePlots_",regionType,"_",seqDate,"_",expName,".pdf"),plot=mp,device="pdf",
           width=20,height=20,units="cm")




    ###############################################################################################
    #####################    Metagene centered on less confidence TSS                   ###########
    ###############################################################################################

    winSize=300
    regionType="lcTSS"
    tssWin<-lessConfTSS
    mcols(tssWin)$TSS<-start(tssWin)
    tssWin<-resize(tssWin,width=winSize,fix="center")

    #################
    # Merge matrix lists that were separated by sample
    ################

    allSampleMats<-mergeSampleMats(path, regionType, samples, deleteSplitFiles=F)

    saveRDS(allSampleMats, paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,".rds"))

    ###################################################
    # convert merged matrices to relative coordinates
    ###################################################

    #allSampleMats<-readRDS(paste0(path,"/rds/allSampleMats_",regionType,"_",seqDate,"_",expName,".rds"))


    allSampleRelCoordMats<-getRelativeCoordMats(matList=allSampleMats, regionGRs=tssWin, regionType=regionType, anchorCoord=winSize/2)


    saveRDS(allSampleRelCoordMats,paste0(path,"/rds/allSampleRelCoordMats_",regionType,"_",seqDate,"_",expName,
						".rds"))


    ###################################################
    # summarise data for Metagene
    ###################################################
    allSampleRelCoordMats<-readRDS(paste0(path,"/rds/allSampleRelCoordMats_",regionType,"_",seqDate,"_",expName,".rds"))

    tssWinRelCoord<-convertGRtoRelCoord(lessConfTSS,winSize,anchorPoint="middle")

    allSampleMetaMethFreqDF<-getAllSampleMetaMethFreq(allSampleRelCoordMats,samples,regionGRs=tssWinRelCoord,
                                                      minReads=10)

    saveRDS(allSampleMetaMethFreqDF,paste0(path,"/rds/allSampleMetaMethFreqDF_",regionType,"_",seqDate,"_",expName,".rds"))


    ### plot metagene by sample
    mp<-plotDSMFmetageneDF(metageneDF=allSampleMetaMethFreqDF,maxPoints=10000)

    ggsave(paste0(path,"/plots/metaGenePlots_",regionType,"_",seqDate,"_",expName,".pdf"),plot=mp,device="pdf",
           width=20,height=20,units="cm")



} # end of gw data only processing

