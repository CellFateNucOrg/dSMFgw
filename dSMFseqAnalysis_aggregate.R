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


#################
# Directories
################

#setup directory structure this to desired location of your alignments
makeDirs(path,c("plots","methCalls","csv","bigwig","rds"))


###############################################################################################
#####################               PLOT C METHYLATION CALLS               ######################
###############################################################################################


########################################################
# read in methylation calls
########################################################

# read methylation calls
methFreqGR<-list()
for (sampleName in samples) {
  baseFileName<-paste0(sampleName,"_",seqDate)
  methFreqGR[[sampleName]]<-getMethFreqGR(baseFileName,paste0(path,"/methCalls"),motifFile,minDepth=10)

  # make bigwig coverage file for all Cs
  Ccovbw<-sort(c(methFreqGR[[sampleName]]$CG,methFreqGR[[sampleName]]$GC, methFreqGR[[sampleName]]$C))
  if("UCSC" %in% GenomeInfoDb::seqlevelsStyle(Ccovbw)){
	Ccovbw_u<-Ccovbw  	
  } else {
	Ccovbw_u<-wbToUcscGR(Ccovbw)
  }

  mcols(Ccovbw_u)<-NULL
  Ccovbw_u$score<-mcols(Ccovbw)$readDepth
  seqlengths(Ccovbw_u)<-seqlengths(Celegans)[1:length(seqlengths(Ccovbw_u))]
  # export coverage data with no smoothing
  rtracklayer::export.bw(Ccovbw_u,con=paste0(path,"/bigwig/Ccov_",sampleName,".bw"))
  line=paste0(urlPrefixForUCSC="http://www.meister.izb.unibe.ch/ucsc/Ccov_",sampleName,".bw")
  write(line,file=paste0(path,"/bigwig/urlsToUpload.txt"),append=TRUE)
  rm(list=c("Ccovbw","Ccovbw_u"))
}



## save as rds for future access
saveRDS(methFreqGR,paste0(path,"/rds/methFreqAllCsGRs_",seqDate,"_",dataType,".rds"))
#methFreqGR<-readRDS(paste0(path,"/rds/methFreqAllCsGRs_",seqDate,"_",dataType,".rds"))



########################################################
# plots of all C methylation coverage
########################################################

# plot C conversion rate for each sample
Cs<-grlToDf(methFreqGR,samples,"C")
Cs$testGroup<-factor(fileList$TestGroup[match(Cs$sample,fileList$SampleName)])
Cs<- Cs %>% group_by(sampleName) %>% mutate(medianCoverage=median(readDepth),
                                               medianMethFreq=median(methylated/readDepth))


Cstats<-Cs %>% group_by(sampleName) %>% summarize(meanCoverage=mean(readDepth),
                                                      medianCoverage=median(readDepth),
                                                      meanMeth=mean(methylated),
                                                      medianMeth=median(methylated),
                                                      allCcount=length(methPercent),
                                                      noMethylatedcount=sum(methPercent==0),
                                                      allMethylatedcount=sum(methPercent==100))
write.csv(Cstats,paste0(path,"/csv/allCs_stats_",seqDate,"_",expName,".csv"),row.names=F)


p1<-ggplot(Cs,aes(x=readDepth,fill=testGroup))+
  geom_histogram(bins=50)+
  facet_wrap(~sampleName) +
  ggtitle("Read coverage of non-CpG/GpC Cs")+
  geom_vline(aes(xintercept=medianCoverage,group=sampleName))

p2<-ggplot(Cs,aes(x=methylated/readDepth,fill=testGroup))+
  geom_histogram(bins=50,binwidth=0.02)+
  geom_histogram(data=subset(Cs,(methylated/readDepth)>0.2&(methylated/readDepth)<=1),fill="black",binwidth=0.02,bins=10)+
  facet_wrap(~sampleName)+
  ggtitle("Fraction non-converted Cs in non-CpG/GpC sites") +
  geom_vline(xintercept=0.2,color="black",linetype="dashed") +
  annotate("text", label = "Low bisulfite\nconversion", x = 0.6, y = dim(Cs)[1]/50, color = "black",size=3)
# None of these sites should be methylated!! High methylation rates mean poor bisulfite conversion

p3<-ggplot(Cs,aes(x=methylated/readDepth,y=readDepth))+
  geom_count(alpha=0.1,col="blue") +
  ggtitle("Fraction non-converted Cs in non-CpG/GpC sites vs coverage")

p<-marrangeGrob(list(p1,p2,p3),ncol=1,nrow=3)
ggsave(paste0(path,"/plots/allC_convRateStats_",seqDate,"_",expName,".pdf"),plot=p,device="pdf",width=20,height=29,units="cm")





#bad<-Cs[Cs$methPercent<80,]
#bad$score<-bad$methPercent
#good<-Cs[Cs$methPercent>=80,]
#good$score<-good$methPercent
#export(bad,paste0(path,"/bedgraph/badCs.bedgraph"),format="bedgraph")
#export(good,paste0(path,"/bedgraph/goodCs.bedgraph"),format="bedgraph")

rm(list=c("Cs","p","p1","p2","p3"))


###############################################################################################
#####################         PLOT METHYLATION in CpG context        ######################
###############################################################################################

# plot CpG methylation rate for each sample
CGs<-grlToDf(methFreqGR,samples,"CG")
CGs$testGroup<-factor(fileList$TestGroup[match(CGs$sample,fileList$SampleName)])
CGs<- CGs %>% group_by(sampleName) %>% mutate(medianCoverage=median(readDepth),
                                            medianMethFreq=median(methylated/readDepth))


CGstats<-CGs %>% group_by(sampleName) %>% summarize(meanCoverage=mean(readDepth),
                                                  medianCoverage=median(readDepth),
                                                  meanMeth=mean(methylated),
                                                  medianMeth=median(methylated),
                                                  allCcount=length(methPercent),
                                                  noMethylatedcount=sum(methPercent==0),
                                                  allMethylatedcount=sum(methPercent==100))
write.csv(CGstats,paste0(path,"/csv/CGs_stats_",seqDate,"_",expName,".csv"),row.names=F)



##################################################################
# histograms of all C methylation frequency in CpG context
##################################################################

p1<-ggplot(CGs,aes(x=readDepth,fill=testGroup))+
  geom_histogram(bins=50)+
  facet_wrap(~sampleName) +
  ggtitle("Read coverage of CGs")+
  geom_vline(aes(xintercept=medianCoverage,group=sampleName))

p2<-ggplot(CGs,aes(x=methylated/readDepth,fill=testGroup))+
  geom_histogram(bins=50,binwidth=0.02)+
  facet_wrap(~sampleName)+
  ggtitle("Fraction methylation CGs") +
  geom_vline(aes(xintercept=medianMethFreq,group=sampleName))

not100<-CGs[CGs$methPercent!=100,]
not100<- not100 %>% group_by(sampleName) %>% mutate(medianCoverage=median(readDepth),
                                              medianMethFreq=median(methylated/readDepth))
p3<-ggplot(not100,aes(x=methylated/readDepth,fill=testGroup))+
  geom_histogram(bins=50,binwidth=0.02)+
  facet_wrap(~sampleName)+
  ggtitle("Fraction methylation CGs excluding 100% methylated positions")+
  geom_vline(aes(xintercept=medianMethFreq,group=sampleName))
# None of these sites should be methylated!! High methylation rates mean poor bisulfite conversion


p<-marrangeGrob(list(p1,p2,p3),ncol=1,nrow=3)
ggsave(paste0(path,"/plots/CGs_methStats_",seqDate,"_",expName,".pdf"),plot=p,device="pdf",width=20,height=29,units="cm")

#bad<-CGs[CGs$methPercent==100,]
#bad$score<-bad$methPercent
#good<-CGs[CGs$methPercent<100,]
#good$score<-good$methPercent
#export(bad,paste0(path,"/bedgraph/badCGs.bedgraph"),format="bedgraph")
#export(good,paste0(path,"/bedgraph/goodCGs.bedgraph"),format="bedgraph")

rm(list=c("CGs","p","p1","p2","p3"))



###############################################################################################
#####################         PLOT METHYLATION in GpC context        ######################
###############################################################################################

# plot GpC methylation rate for each sample
GCs<-grlToDf(methFreqGR,samples,"GC")
GCs$testGroup<-factor(fileList$TestGroup[match(GCs$sample,fileList$SampleName)])
GCs<- GCs %>% group_by(sampleName) %>% mutate(medianCoverage=median(readDepth),
                                              medianMethFreq=median(methylated/readDepth))


GCstats<-GCs %>% group_by(sampleName) %>% summarize(meanCoverage=mean(readDepth),
                                                    medianCoverage=median(readDepth),
                                                    meanMeth=mean(methylated),
                                                    medianMeth=median(methylated),
                                                    allCcount=length(methPercent),
                                                    noMethylatedcount=sum(methPercent==0),
                                                    allMethylatedcount=sum(methPercent==100))
write.csv(GCstats,paste0(path,"/csv/GCs_stats_",seqDate,"_",expName,".csv"),row.names=F)



##################################################################
# histograms of all C methylation frequency in GpC context
##################################################################

p1<-ggplot(GCs,aes(x=readDepth,fill=testGroup))+
  geom_histogram(bins=50)+
  facet_wrap(~sampleName) +
  ggtitle("Read coverage of GCs")+
  geom_vline(aes(xintercept=medianCoverage,group=sampleName))

p2<-ggplot(GCs,aes(x=methylated/readDepth,fill=testGroup))+
  geom_histogram(bins=50,binwidth=0.02)+
  facet_wrap(~sampleName)+
  ggtitle("Fraction methylation GCs") +
  geom_vline(aes(xintercept=medianMethFreq,group=sampleName))

not100<-GCs[GCs$methPercent!=100,]
not100<- not100 %>% group_by(sampleName) %>% mutate(medianCoverage=median(readDepth),
                                                    medianMethFreq=median(methylated/readDepth))
p3<-ggplot(not100,aes(x=methylated/readDepth,fill=testGroup))+
  geom_histogram(bins=50,binwidth=0.02)+
  facet_wrap(~sampleName)+
  ggtitle("Fraction methylation CGs excluding 100% methylated positions")+
  geom_vline(aes(xintercept=medianMethFreq,group=sampleName))
# None of these sites should be methylated!! High methylation rates mean poor bisulfite conversion


p<-marrangeGrob(list(p1,p2,p3),ncol=1,nrow=3)
ggsave(paste0(path,"/plots/GCs_methStats_",seqDate,"_",expName,".pdf"),plot=p,device="pdf",width=20,height=29,units="cm")

#bad<-GCs[GCs$methPercent==100,]
#bad$score<-bad$methPercent
#good<-GCs[GCs$methPercent<100,]
#good$score<-good$methPercent
#export(bad,paste0(path,"/bedgraph/badGCs.bedgraph"),format="bedgraph")
#export(good,paste0(path,"/bedgraph/goodGCs.bedgraph"),format="bedgraph")

rm(list=c("GCs","p","p1","p2","p3"))



###############################################################################################
#####################         combine CG and GC metrices                 ######################
###############################################################################################

cggcFreqGR<-combineCGGCgr(methFreqGR,samples)
# fix change of - to . in sample names when merging gr
colnames(mcols(cggcFreqGR))<-gsub("\\.","-",colnames(mcols(cggcFreqGR)))

# save merged and sorted methylation gr
saveRDS(cggcFreqGR,paste0(path,"/rds/methFreqCombinedCGGC_",seqDate,"_",expName,"_GR.rds"))

rm(list=c("methFreqGR"))
cggcFreqGR<-readRDS(paste0(path,"/rds/methFreqCombinedCGGC_",seqDate,"_",expName,"_GR.rds"))


# keep only methyltion frequency data
onlyM<-grep("_M$",colnames(mcols(cggcFreqGR)))
mcols(cggcFreqGR)<-mcols(cggcFreqGR)[,onlyM]

###################################################
# plot number of CG/GCs with meth data per amplicon
###################################################

if(dataType=="amp"){
  pdf(paste0(path,"/plots/ampliconCoverageWithMethdat_",seqDate,"_",expName,"a.pdf"),width=8,height=11,paper="a4")
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
  allCcount<-countOverlaps(amplicons, cggcFreqGR)
  methCcount<-countOverlaps(amplicons, cggcFreqGR) ##
  plotMethCoveragePerAmplicon(allCcount,methCcount,sites="CG/GCs",samples="(at least 1 sample)")

  # same but only looking at Cs which have data for all samples
  methCcount<-countOverlaps(amplicons, cggcFreqGR[complete.cases(mcols(cggcFreqGR))])
  plotMethCoveragePerAmplicon(allCcount,methCcount,sites="CG/GCs",samples="(in all samples)")

  dev.off()
}

###################################################
# plot correlation between duplicates
###################################################

pdf(paste0(path,"/plots/CorrelationBetweenSamples_",seqDate,"_",expName,".pdf"),width=8,height=11,paper="a4")
cggcFreqGR<-cggcFreqGR[which(rowSums(as.data.frame(mcols(cggcFreqGR))!=1)==length(samples)),]
#avoid plotting too many points - randomly choose 10,000 if more than that
set.seed(1)
if(length(cggcFreqGR)>1e4){
  sampleCs<-sample(1:length(cggcFreqGR),1e4)
} else {
  sampleCs<-1:length(cggcFreqGR)
}



print(testGroup)
lapply(testGroup,function(g) {
    bioReplicates<-c(grep(g,names(mcols(cggcFreqGR))))
    print(bioReplicates)
    p<-ggpairs(as.data.frame(1-as.matrix(mcols(cggcFreqGR[sampleCs, bioReplicates]))),
             title=paste0("Correlation between biological replicates: ",g))
    print(p)
})


# correlation plot for all samples
ggpairs(as.data.frame(1-as.matrix(mcols(cggcFreqGR[sampleCs,]))),
        title="Correlation between all samples")

# simple colour box correlation plot
ggcorr(as.data.frame(1-as.matrix(mcols(cggcFreqGR))),
       low = "steelblue", mid = "white", high = "darkred", label = T,
       label_size = 8, label_color = "white",label_round=2, name="rho")
dev.off()


###################################################
# make bigwig tracks of data
###################################################


w=10 #winSize for smoothing

if(GenomeInfoDb::seqlevelsStyle(cggcFreqGR)!="UCSC"){
  cggcFreqGR_u<-wbToUcscGR(cggcFreqGR)
} else {
  cggcFreqGR_u<-cggcFreqGR
}

dataCols<-names(mcols(cggcFreqGR_u))[grep("_M$",names(mcols(cggcFreqGR_u)))]

dsmf_all<-methToDSMFgr(cggcFreqGR_u,dataCols)
seqlengths(dsmf_all)<-seqlengths(Celegans)

# export raw data with no smoothing
grToBw(dsmf_all,dataCols,bwPath=paste0(path,"/bigwig"),
       filenamePrefix=paste0("rawDSMF_"),
       urlPrefixForUCSC="http://www.meister.izb.unibe.ch/ucsc/")

#smoothe data
smDSMF<-smootheGRdata(dsmf_all,dataCols,winSize=w,winStep=1)
#export smoothed tracks
grToBw(smDSMF,dataCols,bwPath=paste0(path,"/bigwig"),
       filenamePrefix=paste0("w",w,"smDSMF_"),
       urlPrefixForUCSC="http://www.meister.izb.unibe.ch/ucsc/")

rm(list=c("smDSMF","dataCols","dsmf_all","cggcFreqGR"))

