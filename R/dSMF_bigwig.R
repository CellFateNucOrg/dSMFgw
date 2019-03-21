# functions for smoothing and exporting meth data as bigwig

#' Remove rows which have only NAs from QuasR methlyation genomic ranges
#' @param grObj A QuasR methyaltion genomic ranges object with metadata columns for frequency of methylation in 1 or more samples
#' @param datacols The column names that contain methylation frequency data
#' @return A QuasR methylation genomic in which rows which are all NAs have been removed
#' @export
removeAllNAs<-function(grObj,dataCols) {
  # function to remove lines in a GRanges object where the given data columns all
  # have NA values. Input should be a GRanges object, and a vector with the names of the metadata
  # columns that you want to check for NAs (default is first mcols column).
  # Note, that this will not remove rows where at least on sample has a numeric value, in order to
  # avoid a single sample with poor coverage blanking out the table
  numCols<-length(dataCols)
  data2check<-as.data.frame(mcols(grObj)[,dataCols])
  i<-rowSums(is.na(data2check))==numCols
  return(grObj[!i])
}



#' Remove rows which have only NAs from QuasR methlyation genomic ranges
#'
#' \code{padEnds} takes a vector of values that have been smoothed with a window of size winSize
#' and pads the end of the vector with a mirror image of the end values in order to maintain the total
#' vector length. Here it is assumed that the data was smoothed with rollapply aligned to the center
#' of the window (so padding occurs at both ends of the vector)
#'
#' @param vec A vector of values that have been smoothed
#' @param winSize The Size of the window the data was smoothed by.
#' @return Vector of smoothed values with end values repeated to give a vector of the same total length
padEnds<-function(vec,winSize){
  #for roll apply you always lose winSize-1 values.
  padSize<-winSize-1
  leftPad<-padSize%/%2
  rightPad<-padSize-leftPad
  paddedVec<-c(rep(vec[1],leftPad),vec,rep(vec[length(vec)],rightPad))
  return(paddedVec)
}


#' Remove rows which have only NAs from QuasR methlyation genomic ranges
#'
#' \code{smootheGRdata} takes a genomic ranges object with metadata columns and smoothes
#' the columns listed in sampleNames using a sliding window with size winSize and step of
#' winStep. The function NAs in two steps in order to avoid loosing too much data due to low coverage
#' samples. First all rows that have only NAs in them are removed. Next before smoothing individual samples,
#' the NAs in that sample are removed.
#'
#' @param grObj A QuasR methyaltion genomic ranges object with metadata columns for frequency of methylation in 1 or more samples
#' @param sampleNames The column names that contain methylation frequency data that you wish to smoothe.
#' @param winSize The Size of the window the data was smoothed by (default = 1)
#' @param winStep The Size of the step by which to move the sliding window (default = 1)
#' @return Vector of smoothed values with end values repeated to give a vector of the same total length
#' @examples
#' \dontrun{
#' smootheGRdata(methGR,sampleNames,winSize=10,winStep=1)
#' }
#' @export
smootheGRdata<-function(grObj,sampleNames,winSize=10,winStep=1) {
  # First remove all rows that have only NAs in all the sample columns of interest
  grObj_noNA<-removeAllNAs(grObj,sampleNames)
  #cycle through the samples and smooothe the rle and convert back to GRanges
  for (s in sampleNames){
    dataCol<-as.vector(mcols(grObj_noNA)[,s])
    smData<-zoo::rollapply(dataCol,width=winSize,by=winStep, FUN=mean, na.rm=T)
    #pad smData ends with first and last values to get same length vector as input
    mcols(grObj_noNA)[,s]<-padEnds(smData,winSize)
  }
  return(grObj_noNA)
}



#' Convert QuasR genomic ranges with methylation frequency data to DSMF values
#'
#' \code{methToDSMFgr} Converts methylation frequency data to DSMF values (1-methylationFrequency).
#' Note that metadata columns that are all NAs will be removed, so Granges object might change size.
#' NAs in individual samples are ignored during conversion.
#'
#' @param methGR A QuasR methyaltion with metadata columns for frequency of methylation in 1 or more samples
#' @param sampleNames The column names that you wish to convert.
#' @return A genomic ranges
#' @export
methToDSMFgr<-function(methGR,sampleNames) {
  # First remove all rows that have only NAs in all the sample columns of interest
  grObj_noNA<-removeAllNAs(methGR,sampleNames)
  #cycle through the samples and remove NAs so be able to subtract frequencies
  for (s in sampleNames){
    dataCol<-as.vector(mcols(grObj_noNA)[,s])
    idx<-is.na(dataCol)
    dsmf<-1-dataCol[!idx]
    #pad smData ends with first and last values to get same length vector as input
    mcols(grObj_noNA)[!idx,s]<-dsmf
  }
  return(grObj_noNA)
}



#' Export genomic ranges as bigwig file
#'
#' \code{grToBw} exports metadata columns listed in sampleNames as bigwig files. The files will be
#' saved with the an filenamePrefix given and the name of the sample from the metadata column.
#' This script can also create a file called urlsToUpload.txt with a list if filenames. The filenames
#' can also be joined to the url where the files will be stored if provided as a urlPrefixForUCSC.
#'
#' @param grObj genomicRanges object with metadata columns for 1 or more samples
#' @param sampleNames The column names that you wish to export.
#' @param bwPath Path to directory where files should be saved (default = "./bigwig")
#' @param filenamePrefix Prefix to add to the file name before the name of the sample (default = "")
#' @param urlPrefixForUCSC url where new bigwigs will be located for upload to ucsc as custom tracks
#' @return A bigwig file for each metadata column listed in sampleNames
#' @examples
#' \dontrun{
#'  grToBw(dsmf_all,dataCols,bwPath=paste0(path,"/bigwig"),filenamePrefix="rawDSMF_")
#' }
#' @export
grToBw<-function(grObj,sampleNames,bwPath="./bigwig",filenamePrefix="",urlPrefixForUCSC="") {
  if (!dir.exists(bwPath)) {
    dir.create(bwPath,recursive=TRUE)
  }
  if (!file.exists(paste0(bwPath,"/urlsToUpload.txt"))) {
    file.create(paste0(bwPath,"/urlsToUpload.txt"))
  }
  IRanges::width(grObj)<-1
  for (s in sampleNames){
    idx<-is.na(mcols(grObj)[,s])
    newGR<-grObj[!idx]
    mcols(newGR)<-NULL
    newGR$score<-mcols(grObj)[!idx,s]
    shortName<-gsub("_M$","",s)
    rtracklayer::export.bw(newGR,con=paste0(bwPath,"/",filenamePrefix,shortName,".bw"))
    line=paste0(urlPrefixForUCSC,filenamePrefix,shortName,".bw")
    write(line,file=paste0(bwPath,"/urlsToUpload.txt"),append=TRUE)
  }
  return(print("bigwig files written successfully"))
}








######################################### from different source



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

