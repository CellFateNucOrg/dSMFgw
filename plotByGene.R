library(grid)
library(png)
library(ggplot2)
library(gridExtra)

source("./R/variableSettings.R")
regionType="ampTSS"
matTable<-readRDS(paste0(path,"/rds/allSampleRelCoordMats_",regionType,"_",seqDate,"_",expName,".rds"))
head(matTable)
matTable<-matTable[!is.na(matTable$filename),]
empath=paste0(path,"/EMs_euc_50r_w30_NA20_267413_m1To1_vNA0.5")
if(!dir.exists(paste0(empath,"/combinedByRegion"))){
  dir.create(paste0(empath,"/combinedByRegion"), showWarnings=F)
}

for(region in unique(matTable$region)){
  pngList<-list.files(path=empath,pattern=paste0("finalClassReads_.*_",region,"_K10.png"),full.names=T) 
  if(length(pngList)>0){
    plots <- lapply(pngList,function(x){
      img <- as.raster(readPNG(x))
      rasterGrob(img, interpolate = FALSE)
    })
    chr<-as.character(seqnames(ampTSS)[grep(region,ampTSS$ID)]) 
    ggsave(paste0(empath,"/combinedByRegion/finalClassReads_",region,"_K10.png"),width=8.5, height=11, 
        dev="png",marrangeGrob(grobs=plots, nrow=length(sampleNames)/2, ncol=2,top=paste0("chr: ",chr)))
    idx<-NULL
  }
}
