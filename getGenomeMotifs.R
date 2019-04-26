
source("./R/dSMF.R")
source("./R/grangesUtils.R")

args<-commandArgs(trailingOnly=TRUE)

genomeFile<-args[1]

genomeFile<-"/home/ubelix/izb/semple/genomeVer/ws260/sequence/c_elegans.PRJNA13758.WS260.genomic.fa"

gnmgr<-findGenomeMotifs(genomeFile)

outdir<-dirname(genomeFile)
outfile<-gsub("\\.fa","\\.CGGC_motifs.RDS",basename(genomeFile))
saveRDS(gnmgr,paste0(outdir,"/",outfile))


