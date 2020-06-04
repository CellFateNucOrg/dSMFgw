
source("./R/dSMF.R")
source("./R/grangesUtils.R")

args<-commandArgs(trailingOnly=TRUE)

genomeFile<-args[1]

#genomeFile<-"/home/ubelix/izb/semple/genomeVer/ws260/sequence/c_elegans.PRJNA13758.WS260.genomic.fa"

gnmgr<-findGenomeMotifs(genomeFile)

outdir<-dirname(genomeFile)
outfile<-gsub("\\.fa","\\.CGGC_motifs.RDS",basename(genomeFile))
saveRDS(gnmgr,paste0(outdir,"/",outfile))


if (class(genomeFile)=="BSgenome"){
        DNAss<-BSgenomeToDNAStringSet(genomeFile)
        #fa<-genomeFile
} else if (is.character(genomeFile)){
        DNAss<-Biostrings::readDNAStringSet(genomeFile)
        #fa<-Rsamtools::FaFile(genomeFile)
}


# make bed files for pileup single molecule methylation calling

# Cs in non-methylated context
#gnmgr<-readRDS(paste0(outdir,"/",outfile))
Cs<-mIdxToGR(Biostrings::vmatchPattern("C",DNAss))
ol<-GenomicRanges::findOverlaps(Cs,gnmgr)
Cs<-Cs[-S4Vectors::queryHits(ol),]
#refbase<-Biostrings::getSeq(fa,Cs)
GenomicRanges::strand(Cs)<-"+"
outfile<-gsub("\\.fa","\\.C.bed",basename(genomeFile))
rtracklayer::export.bed(Cs,paste0(outdir,"/",outfile))

# on rev strand
Gs<-mIdxToGR(Biostrings::vmatchPattern("G",DNAss))
ol<-GenomicRanges::findOverlaps(Gs,gnmgr)
Gs<-Gs[-S4Vectors::queryHits(ol),]
#refbase<-Biostrings::getSeq(fa,Cs)
GenomicRanges::strand(Gs)<-"+"
outfile<-gsub("\\.fa","\\.G.bed",basename(genomeFile))
rtracklayer::export.bed(Gs,paste0(outdir,"/",outfile))



# Non unique CGs
CGs<-mIdxToGR(Biostrings::vmatchPattern("CG",DNAss))
fCGs<-CGs
GenomicRanges::strand(fCGs)<-"+"
rCGs<-CGs
GenomicRanges::strand(rCGs)<-"-"
CGs<-sort(c(fCGs,rCGs),ignore.strand=T)
CGs<-GenomicRanges::resize(CGs,width=1,fix="start")
outfile<-gsub("\\.fa","\\.CG.bed",basename(genomeFile))
rtracklayer::export.bed(CGs,paste0(outdir,"/",outfile))


# Non unique GCs
GCs<-mIdxToGR(Biostrings::vmatchPattern("GC",DNAss))
fGCs<-GCs
GenomicRanges::strand(fGCs)<-"+"
rGCs<-GCs
GenomicRanges::strand(rGCs)<-"-"
GCs<-sort(c(fGCs,rGCs),ignore.strand=T)
GCs<-GenomicRanges::resize(GCs,width=1,fix="start")
outfile<-gsub("\\.fa","\\.GC.bed",basename(genomeFile))
rtracklayer::export.bed(GCs,paste0(outdir,"/",outfile))
