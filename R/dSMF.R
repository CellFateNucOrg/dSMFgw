#' Do a single molecule plot of dSMF data
#' @param mat A methylation matrix
#' @param regionName The name of the region to be plotted (usually gene name). Must match Ids in reagionGRs
#' @param regionGRs Genomic Ranges which include region to be plotted (ID metadata column must match regionName)
#' @param featureGRs Genomic Ranges of any features of interest (e.g. TSS) that you wish to be plotted as vertical red line
#' @param myXlab Text to be used for X axis label (default="CpG/GpC position")
#' @param featureLabel  Text to be used to label special feature you are plotting (default="TSS")
#' @param title Text to be used as plot title (default=NULL will results in the plot named accorind to the regionGR)
#' @param baseFontSize The base font size to be used for the plot (default is 12)
#' @return A ggplot2 plotting object
#' @examples
#' plotSingleMoleculesAmp(methMat,"ttTi5605",ttTi5605gr,featureGRs=c(),myXlab="CpG position",featureLabel="TSS",title=NULL,baseFontSize=12)
#' @export
plotSingleMoleculesAmp<-function(mat,regionName,regionGRs,featureGRs=c(),myXlab="CpG/GpC position",featureLabel="TSS",title=NULL,
                                 baseFontSize=12) {
  ### single molecule plot. mat is matrix containing methylation values at different postions (columns) in
  # individual reads (rows). regionName is the ID of the amplicon or genomic region being plotted. regionGRs is a
  # genomicRanges object containing the region being plotted. one of its mcols must have a name "ID" in which the
  # same ID as in regionName appears. featureGRs is genomic ranges object for plotting location of some feature in
  # the region, such as the TSS. myXlab is the X axis label. featureLabel is the label for the type of feature that
  # will be plotted underneath the feature
  if(dim(mat)[1]>10) {
    regionGR<-regionGRs[match(regionName,regionGRs$ID)]
    if (length(featureGRs)>0) {
      featGR<-featureGRs[match(regionName,featureGRs$ID)]
    }
    na.matrix<-is.na(mat)
    mat[na.matrix]<- -1
    # try to perform heirarchical clustering
    hc <- try(
      hclust(dist(apply(mat,2,as.numeric))),
      silent = TRUE)
    mat[na.matrix]<-NA
    if (class(hc) == "try-error") {
      df<-as.data.frame(mat)
      print("hclust failed. Matrix dim: ")
      print(dim(mat))
    } else {
      df<-as.data.frame(mat[hc$order,])
    }

    reads<-row.names(df)
    d<-tidyr::gather(df,key=position,value=methylation)
    d$molecules<-seq_along(reads)
    d$methylation<-as.character(d$methylation)
    d$position<-as.numeric(d$position)
    if (is.null(title)) {
      title=paste0(regionName, ": ", GenomeInfoDb::seqnames(regionGR)," ", GenomicRanges::strand(regionGR),"ve strand")
    }
    p<-ggplot(d,aes(x=position,y=molecules,width=2)) +
      geom_tile(aes(width=6*GenomicRanges::width(regionGR)/500,fill=methylation),alpha=0.8) +
      scale_fill_manual(values=c("0"="black","1"="grey80"),na.translate=F,na.value="white",
                        labels=c("protected","accessible"),name="dSMF") + theme_light(base_size=baseFontSize) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold",hjust = 0.5),
            legend.position="bottom", legend.key.height = unit(0.1, "cm")) +
      ggtitle(title) +
      xlab(myXlab) + ylab("Single molecules") + xlim(GenomicRanges::start(regionGR),GenomicRanges::end(regionGR)+10)
    if(length(featureGRs)>0) {
      p<-p+geom_linerange(aes(x=GenomicRanges::start(featGR), y=NULL, ymin=0, ymax=length(reads)+max(3,0.04*length(reads))),col="red") +
        ggplot2::annotate("segment", x = GenomicRanges::start(featGR), xend = GenomicRanges::start(featGR)+20*ifelse(GenomicRanges::strand(featGR)=="-",-1,1),
                          y = length(reads)+max(3,0.04*length(reads)), yend =length(reads)+max(3,0.04*length(reads)), colour = "red",
                          arrow=arrow(length = unit(0.3, "cm")), size=0.7) +
        ggplot2::annotate(geom="text", x=GenomicRanges::start(featGR), y=-max(2,0.03*length(reads)), label=featureLabel,color="red")
    }
  } else {
    p<-NULL
  }
  return(p)
}



#' Adds overlapping sites from first GRanges to the second
#' @param gr1 A GenomicRanges object
#' @param gr2 A GenomicRanges object
#' @return A GenomicRanges object containing all of gr2, and any GRanges from gr1 that overlap with gr1
#' @examples
#' getOverlappingSites(gr1,gr2)
#' @export
getOverlappingSites<-function(gr1,gr2) {
  # adds any sites from gr1 that overlap with gr2 to gr2
  ol<-GenomicRanges::findOverlaps(gr1,gr2,ignore.strand=T)
  if (length(ol)>0) {
    overlapping<-c(gr2,gr1[S4Vectors::queryHits(ol)])
  } else {
    overlapping<-gr2
  }
  return(overlapping)
}


#' Adds unique sites from first GRanges that are not in the second
#' @param gr1 A GenomicRanges object
#' @param gr2 A GenomicRanges object
#' @return A GenomicRanges object containing all of gr1 ranges that do not overlap with gr2
#' @examples
#' getIsolatedSites(gr1,gr2)
#' @export
getIsolatedSites<-function(gr1,gr2) {
  # removes from gr1 and sites which overlap with gr2
  ol<-GenomicRanges::findOverlaps(gr1,gr2,ignore.strand=T)
  if (length(ol)>0) {
    isolated<-gr1[-S4Vectors::queryHits(ol)]
  } else {
    isolated<-gr1
  }
  return(isolated)
}

#' Splits GenomicRange into multiple ranges of width 2.
#' @param gr A GenomicRanges object
#' @return A GenomicRanges object with multiple non-overlapping ranges width 2.
#' @examples
#' splitRangeInPairs(gr1)
#' @export
splitRangeInPairs<-function(gr) {
  starts<-seq(GenomicRanges::start(gr),GenomicRanges::end(gr),by=2)
  newGR<-GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gr),ranges=IRanges::IRanges(start=starts,width=2),strand="*")
  return(newGR)
}


#' Splits GenomicRange into multiple ranges of width 2 or 3.
#' @param gr A GenomicRanges object
#' @return A GenomicRanges object with multiple non-overlapping ranges width 2 or 3.
#' @examples
#' splitGCGCruns(gr1)
#' @export
splitGCGCruns<-function(gr) {
  # function to (arbitrarily) decompose overlapping GCs and CGs to triplets and doublets
  # create new gr for split up ranges
  olGCGC<-GenomicRanges::GRanges()
  # find ranges that do not need splitting
  width3<-GenomicRanges::width(gr)==3
  olGCGC<-append(olGCGC,gr[width3])
  # find ranges that need splitting
  toSplit<-gr[!width3]
  if (length(toSplit)>0) {
    for (i in 1:length(toSplit)) {
      isOdd<-GenomicRanges::width(toSplit[i])%%2==1
      if (isOdd) {
        # if there are an odd number of bases, chop off the first three as a triplet, and the rest in pairs
        firstTriplet<-GenomicRanges::resize(toSplit[i],3,"start")
        evenRange<-GenomicRanges::resize(toSplit[i],GenomicRanges::width(toSplit[i])-3,"end")
        olGCGC<-append(olGCGC,firstTriplet)
        olGCGC<-append(olGCGC,splitRangeInPairs(evenRange))
      } else {
        # if there are an even number of bases, chop the gr into smaller gr length2
        olGCGC<-append(olGCGC,splitRangeInPairs(toSplit[i]))
      }
    }
  }
  return(sort(olGCGC))
}


#' Finds all CG, GC and long GCGC stretches in DNAStringSet
#' @param DNAss A DNAStringSet object
#' @return A list of three genomic ranges for isolated CGs, isolated GCs, and CG/GC/GCG/CGC from longer runs
#' @examples
#' findNonOverlappingMotifs(DNAss)
#' @export
findNonOverlappingMotifs<-function(DNAss) {
  #search for patterns in sequence
  CGs<-mIdxToGR(Biostrings::vmatchPattern("CG",DNAss))
  GCs<-mIdxToGR(Biostrings::vmatchPattern("GC",DNAss))
  GCGs<-mIdxToGR(Biostrings::vmatchPattern("GCG",DNAss))
  CGCs<-mIdxToGR(Biostrings::vmatchPattern("CGC",DNAss))
  GCGCs<-sort(c(GCGs,CGCs))

  # move any GC sites that overlap with GCG to the GCG list (removing them from GC list)
  GCGCs<-getOverlappingSites(GCs,GCGCs) # always get overlapping first, otherwise will lose them
  GCs<-GenomicRanges::reduce(getIsolatedSites(GCs,GCGCs),ignore.strand=T) # collapse GR on opposite strands
  GCs$context<-"GCH"

  # move any CG sites that overlap with GCG to the GCG list (removing them from CG list)
  GCGCs<-getOverlappingSites(CGs,GCGCs) # always get overlapping first, otherwise will lose them
  CGs<-GenomicRanges::reduce(getIsolatedSites(CGs,GCGCs),ignore.strand=T) # collapse GR on opposite strands
  CGs$context<-"HCG"

  # reduce all the overlapping set to a smallest merged set
  GCGCs<-GenomicRanges::reduce(GCGCs, ignore.strand=T, min.gapwidth=0L)
  # then arbitrarily split them into non-overlapping GR 2-3 bp long
  GCGCs<-splitGCGCruns(GCGCs)
  GCGCs$context<-"GCGorCGC"
  allGR<-list(CGs,GCs,GCGCs)
  names(allGR)<-c("CG","GC","GCG")
  return(allGR)
}


#' Find non overlapping CG, GC and GCG/CGC motifs in genome
#' @param genomeFile A DNAstringSet, path to FASTA file or BSgenome for the genome sequence.
#' @return A GenomicRanges containing non ovelapping CG GC and GCGorCGC motifs. The context is shown in the metadata
#' @examples
#' findGenomeMotifs(ttTi5605dna)
#' @export
findGenomeMotifs<-function(genome){
  # deal with different input genome formats
  if (class(genome)=="BSgenome"){
    genome<-BSgenomeToDNAStringSet(genome)
  } else if (is.character(genome)){
    genome<-Biostrings::readDNAStringSet(genome)

  }
  #strip off anything after a space to deal with additional fasta header info
  names(genome)<-gsub("[[:space:]].*$","",names(genome),perl=F)
  grl<-findNonOverlappingMotifs(genome)
  names(grl)<-NULL
  gr<-do.call("c",grl)
  gr<-sort(gr,ignore.strand=T)
  return(gr)
}



