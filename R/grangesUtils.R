# general useful functions for granges

#' Convert DNAStringSet to GenomicRanges
#' @param DNAss A DNAStringSet
#' @param defaultStrand One of "+", "-" or "*" for the strand from which the sequences come from (default="*")
#' @return A GenomicRanges object
#' @examples
#' DNAstringsetToGR(Biostrings::DNAStringSet(x=c(chr1="agctgtagct",chr2="agagagagttt"))
#' @export
DNAStringSetToGR<-function(DNAss,defaultStrand="*") {
  #extract only the first word in fasta '>' row to be the seq name
  chrNames<-sapply(strsplit(names(DNAss)," "),'[[',1)
  gr<-GenomicRanges::GRanges(seqnames=chrNames,
                             ranges=IRanges::IRanges(start=1,width=Biostrings::width(DNAss)),
                             strand=defaultStrand)
  return(gr)
}



#' Convert BSgenome to DNAStringSet
#' @param BSgnm A BSgenome object
#' @return A DNAStringSet object
#' @examples
#' BSgenomeToDNAStringSet(BSgenome.Celegans.UCSC.ce11::Celegans)
#' @export
BSgenomeToDNAStringSet<-function(BSgnm){
  listOfChr<-sapply(seqnames(BSgnm),function(x){BSgnm[[x]]})
  DNAss<-as(listOfChr,"DNAStringSet")
  return(DNAss)
}



#' Convert BSgenome to GenomicRanges
#' @param BSgnm A BSgenome object
#' @param defaultStrand One of "+", "-" or "*" for the strand from which the sequences come from (default="*")
#' @return A GenomicRanges object
#' @examples
#' BSgenomeToGR(BSgenome.Celegans.UCSC.ce11::Celegans)
#' @export
BSgenomeToGR<-function(BSgnm,defaultStrand="*") {
  gr<-GenomicRanges::GRanges(seqnames=BSgenome::seqnames(BSgnm),
                             ranges=IRanges::IRanges(start=1,end=GenomeInfoDb::seqlengths(BSgnm)),
                             strand="*")
  return(gr)
}



#' Convert UCSC to wormbase chromosome names
#' @param ucscGR A GRanges object with ucsc chomosome names ("chrI"..."chrM")
#' @return A GenomicRanges object
#' @examples
#' ucscToWbGR(BSgenomeToGR(BSgenome.Celegans.UCSC.ce11::Celegans))
#' @export
ucscToWbGR<-function(ucscGR) {
  wbGR<-ucscGR
  GenomeInfoDb::seqlevels(wbGR)<-gsub("chr","",GenomeInfoDb::seqlevels(wbGR))
  GenomeInfoDb::seqlevels(wbGR)<-gsub("M","MtDNA",GenomeInfoDb::seqlevels(wbGR))
  if(class(ucscGR)=="BSgenome") {
    wbGR@provider<-"Wormbase"
    wbGR@provider_version<-"WS235"
  }
  return(wbGR)
}



#' Convert wormbase to ucsc chromosome names
#' @param wbGR A GRanges object with ucsc chomosome names ("chrI"..."chrM")
#' @return A GenomicRanges object
#' @examples
#' wbToUcscGR(ucscToWbGR(BSgenomeToGR(BSgenome.Celegans.UCSC.ce11::Celegans)))
#' @export
wbToUcscGR<-function(wbGR) {
  ucscGR<-wbGR
  GenomeInfoDb::seqlevels(ucscGR)<-gsub("MtDNA","M",GenomeInfoDb::seqlevels(ucscGR))
  GenomeInfoDb::seqlevels(ucscGR)<-paste0("chr",GenomeInfoDb::seqlevels(ucscGR))
  return(ucscGR)
}



#' Convert MIndex object to Genomic Ranges object
#' @param mIdx A MIndex object
#' @return A GenomicRanges object
#' @examples
#' mIdxToGR(Biostrings::vmatchPattern("ATTTAGGGTTTTAGAATACTGCCATTAATTAAAAAT",ttTi5605dna))
#' @export
mIdxToGR<-function(mIdx) {
  #function to convert Mindex object (obtained matching patterns on DNAstringset) to genomic ranges
  allGR<-GenomicRanges::GRanges()
  GenomeInfoDb::seqlevels(allGR)<-names(mIdx)
  for (n in names(mIdx)) {
    if (length(mIdx[[n]])>0) {
      gr<-GenomicRanges::GRanges(seqnames=S4Vectors::Rle(c(n),length(mIdx[[n]])),
                     ranges=mIdx[[n]], strand="*")
      allGR<-append(allGR,gr)
    }
  }
  return(allGR)
}




#' Apply a function on a data from gr2 with windows provided by gr1
#' @param gr1 A GenomicRanges object with windows of interest
#' @param gr2 A GenomicRanges object with data of interest
#' @param applyTo Column name in gr2 on which to apply the function
#' @param fun Function to apply (e.g. sum, mean, paste0)
#' @return Ranges from gr1 with summed values of metadata column from gr2
#' @examples
#' gr1 <- GRanges("chr1", IRanges(c(1,3,7), c(5,6,10),names=paste0("win", letters[1:3])), score=4:6)
#' gr2 <- GRanges("chr1", IRanges(c(1, 3, 8), c(1, 3, 8),names=paste0("dataID:", letters[1:3])), score=c(10,20,30))
#' applyGRonGR(gr1,gr2,applyTo="score",fun=sum)
#' @export
applyGRonGR<-function(gr1,gr2,applyTo,fun,...){
  newGR<-IRanges::subsetByOverlaps(gr1,gr2)
  ol<-IRanges::findOverlaps(gr1,gr2)
  newData<-sapply(split(mcols(gr2)[S4Vectors::subjectHits(ol),applyTo],S4Vectors::queryHits(ol)),fun,...)
  mcols(newGR)[,get("applyTo")]<-newData
  return(newGR)
}




