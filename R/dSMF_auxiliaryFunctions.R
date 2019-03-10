# functions used in dSMF analysis scripts
library(QuasR)
library(data.table)
library(tidyr)

######################################################################################################
##################### dSMF specific functions ########################################################
######################################################################################################

################################################################
######## functions for extracting meth Genomic Ranges ##########
################################################################

miObj2gr<-function(miObj) {
  #function to convert Mindex object (obtained from applying GR on DNAstringset) to genomic ranges
  # this function is required by call_context_methylation
  allGR<-GRanges()
  seqlevels(allGR)<-names(miObj)
  for (n in names(miObj)) {
    if (length(miObj[[n]])>0) {
      grObj<-GRanges(seqnames=Rle(c(n),length(miObj[[n]])),
                     ranges=miObj[[n]], strand="*")
      allGR<-append(allGR,grObj)
    }
  }
  return(allGR)
}

#miObj2gr1<-function(miObj) {
#  #function to convert Mindex object (obtained from applying GR on DNAstringset) to genomic ranges
#  # this function is required by call_context_methylation
#  allGR<-GRanges()
#  seqlevels(allGR)<-seqlevels(miObj)
#  for (n in names(miObj)) {
#    if (length(miObj[[n]])>0) {
#      grObj<-GRanges(seqnames=Rle(c(n),length(miObj[[n]])),
#                     ranges=miObj[[n]], strand="*")
#      allGR<-append(allGR,grObj)
#    }
#  }
#  return(allGR)
#}
#
#
call_context_methylation<-function(meth_gr,c0,genome=genome){
  # call_context_methylation takes in a genomic ranges object will all C positions in the genome
  # and mcols containing the total coveage (_T) and methylation count (_M) for each sample.
  # it returns a list of three GRanges objects "CG" "GC" and "GCG"
  # in each of these, the V1 column contains fraction methylation and the 'type' column shows C context

  # this function can use either a BSgenome or a fasta file with genome sequence
  # fasta file needs to be read in. BSgenome should be preloaded in the environment.
  if (class(genome)=="BSgenome"){
    fastaFlag=FALSE
  } else if (is.character(genome)){
    genome<-readDNAStringSet(genome)
    fastaFlag=TRUE
  } else {
    print("genome must be BSgenome or path to fasta file")
  }

  # find CGs and GCs in genome
  genome_CGs <- vmatchPattern("CG",genome)
  genome_GCs <- vmatchPattern("GC",genome)

  # fix problems with counts depending on genome input method
  if (fastaFlag==FALSE) {
    #if the genome is a BSgenome object then you have to get rid
    #of duplicate matches on positive and negative strand
    genome_CGs <- genome_CGs[strand(genome_CGs)=="+",]
    strand(genome_CGs) <- "*"

    genome_GCs <- genome_GCs[strand(genome_GCs)=="+",]
    strand(genome_GCs) <- "*"
  } else {
    # if the genome is fasta file, we need to convert Mindex object
    #  to GRanges object
    genome_CGs<-miObj2gr(genome_CGs)
    genome_GCs<-miObj2gr(genome_GCs)
    # get rid of additional info in fasta file header and keep only sequence name
    seqlevels(genome_CGs)=strsplit(seqlevels(genome_CGs)," ")[[1]][1]
    seqlevels(genome_GCs)=strsplit(seqlevels(genome_GCs)," ")[[1]][1]
  }

  # get indeces of overlaps between allCs in meth_gr and Cs in CpG or GpC context
  sel_CGs <- meth_gr %over% genome_CGs
  sel_GCs <- meth_gr %over% genome_GCs

  # subset meth_gr (which contains read count data) into matrices for each context
  meth_CGs_gr <- meth_gr[sel_CGs]
  meth_GCs_gr <- meth_gr[sel_GCs]

  ##################
  # collapse strands
  ##################
  # due to symmetric nature of CG, we have counts of two Cs from opposite strands.
  # Previous approach of taking just every other row does not work if half CG or GC sites
  # are present (due to amplicon boundaries). Instead, CG/GC sites will be expanded in a strand
  # specific way and the overlaps combined

  # get an indicator variable for the strand of the gr
  plusIndex<-strand(meth_CGs_gr)=="+"
  # extend the GR by 1 in different directions, depending on the strand
  end(meth_CGs_gr[plusIndex])<-end(meth_CGs_gr[plusIndex])+1
  start(meth_CGs_gr[!plusIndex])<-start(meth_CGs_gr[!plusIndex])-1
  # find 2bp overlaps between gr on different strands.
  ol<-findOverlaps(meth_CGs_gr[plusIndex],meth_CGs_gr[!plusIndex], minoverlap=2,ignore.strand=TRUE)
  # sum the methylation counts from both strands for a given CG site.
  meth_CGsCol_gr<-meth_CGs_gr[plusIndex][queryHits(ol)]
  strand(meth_CGsCol_gr)<-"*"
  values(meth_CGsCol_gr)<-as.matrix(values(meth_CGs_gr[plusIndex][queryHits(ol)]))+as.matrix(values(meth_CGs_gr[!plusIndex][subjectHits(ol)]))
  # we could add back the half CG sites, but as they are at the edge, they probably are distorted by the primer
  # seq. (in two cases i checked, all columns had 0 methylation). The command to add them back would be:
  # meth_CGsCol_gr<-c(meth_CGsCol_gr,meth_CGs_gr[plusIndex][!(seq(length(meth_CGs_gr[plusIndex])) %in% queryHits(ol))],
  #                  meth_CGs_gr[!plusIndex][!(seq(length(meth_CGs_gr[!plusIndex])) %in% subjectHits(ol))])


  # do the same to the GC matrix
  # get an indicator variable for the strand of the gr
  plusIndex<-strand(meth_GCs_gr)=="+"
  # extend the GR by 1 in different directions, depending on the strand (extends opposite way from CG)
  start(meth_GCs_gr[plusIndex])<-start(meth_GCs_gr[plusIndex])-1
  end(meth_GCs_gr[!plusIndex])<-end(meth_GCs_gr[!plusIndex])+1
  # find 2bp overlaps between gr on different strands.
  ol<-findOverlaps(meth_GCs_gr[plusIndex],meth_GCs_gr[!plusIndex], minoverlap=2,ignore.strand=TRUE)
  # sum the methylation counts from both strands for a given CG site.
  meth_GCsCol_gr<-meth_GCs_gr[plusIndex][queryHits(ol)]
  strand(meth_GCsCol_gr)<-"*"
  values(meth_GCsCol_gr)<-as.matrix(values(meth_GCs_gr[plusIndex][queryHits(ol)]))+as.matrix(values(meth_GCs_gr[!plusIndex][subjectHits(ol)]))


  #####################
  #filter for coverage
  #####################

  # identify columns that have total and methylated C counts
  Tcounts=grep('_T\\>',colnames(elementMetadata(meth_CGsCol_gr)))
  Mcounts=grep('_M\\>',colnames(elementMetadata(meth_CGsCol_gr)))

  ######
  #CGs
  ######

  # calculate frequency of methyltion
  CG.met.mat=as.matrix(elementMetadata(meth_CGsCol_gr)[,Mcounts])/as.matrix(elementMetadata(meth_CGsCol_gr)[,Tcounts])
  # filter for coverage (only positions that are covered by >c0 reads)
  CovFilter=as.matrix(elementMetadata(meth_CGsCol_gr)[,Tcounts])>c0
  # put NA into positions where coverage is <c0
  CG.met.mat[!CovFilter]<-NA
  #bind the GRanges with the scores
  CG.met=meth_CGsCol_gr
  elementMetadata(CG.met)=CG.met.mat


  ######
  #GCs
  ######

  # calculate frequency of methyltion
  GC.met.mat=as.matrix(elementMetadata(meth_GCsCol_gr)[,Mcounts])/as.matrix(elementMetadata(meth_GCsCol_gr)[,Tcounts])
  # filter for coverage (only positions that are covered by >c0 reads)
  CovFilter=as.matrix(elementMetadata(meth_GCsCol_gr)[,Tcounts])>c0
  # put NA into positions where coverage is <c0
  GC.met.mat[!CovFilter]<-NA
  #bind the GRanges with the scores
  GC.met=meth_GCsCol_gr
  elementMetadata(GC.met)=GC.met.mat

  ###########################
  #create an unified object
  ###########################

  # need to label GCGs that will appear in both matrices as they form GCH and HCG
  # this is to avoid double counting of GCGs later on when CpGs and GpCs are combined
  oGCG=as.matrix(findOverlaps(resize(GC.met,1,fix='end'),resize(CG.met,1,fix='start'),type='equal'))
  GC.met$type='GCH'
  CG.met$type='CGH'
  GC.met$type[oGCG[,1]]='GCG'
  CG.met$type[oGCG[,2]]='GCG'

  # extract different types of sites
  GCG.met<-sort(c(CG.met[CG.met$type=="GCG"],GC.met[GC.met$type=="GCG"]))
  CG.met<-sort(CG.met[CG.met$type!="GCG"])
  GC.met<-sort(GC.met[GC.met$type!="GCG"])

  # combine all three  GR objects into a single list
  umet=list(CG.met,GC.met,GCG.met)
  names(umet)=c('CG','GC','GCG')
  return(umet)
}



call_context_methylation_counts<-function(meth_gr,c0,genome=genome){
  # call_context_methylation takes in a genomic ranges object will all C positions in the genome
  # and mcols containing the total coveage (_T) and methylation count (_M) for each sample.
  # it returns a list of two GRanges objects "CG" and "GC"
  # in each of these, the V1 column contains fraction methylation and the 'type' column shows C context

  # this function can use either a BSgenome or a fasta file with genome sequence
  # fasta file needs to be read in. BSgenome should be preloaded in the environment.
  if (class(genome)=="BSgenome"){
    fastaFlag=FALSE
  } else if (is.character(genome)){
    genome<-readDNAStringSet(genome)
    fastaFlag=TRUE
  } else {
    print("genome must be BSgenome or path to fasta file")
  }

  # find CGs and GCs in genome
  genome_CGs <- vmatchPattern("CG",genome)
  genome_GCs <- vmatchPattern("GC",genome)

  # fix problems with counts depending on genome input method
  if (fastaFlag==FALSE) {
    #if the genome is a BSgenome object then you have to get rid
    #of duplicate matches on positive and negative strand
    genome_CGs <- genome_CGs[strand(genome_CGs)=="+",]
    strand(genome_CGs) <- "*"

    genome_GCs <- genome_GCs[strand(genome_GCs)=="+",]
    strand(genome_GCs) <- "*"
  } else {
    # if the genome is fasta file, we need to convert Mindex object
    #  to GRanges object
    genome_CGs<-miObj2gr(genome_CGs)
    genome_GCs<-miObj2gr(genome_GCs)
  }

  # get indeces of overlaps between allCs in meth_gr and Cs in CpG or GpC context
  sel_CGs <- meth_gr %over% genome_CGs
  sel_GCs <- meth_gr %over% genome_GCs

  # subset meth_gr (which contains read count data) into matrices for each context
  meth_CGs_gr <- meth_gr[sel_CGs]
  meth_GCs_gr <- meth_gr[sel_GCs]

  ##################
  # collapse strands
  ##################
  # due to symmetric nature of CG, we have counts of two Cs from opposite strands.
  # Previous approach of taking just every other row does not work if half CG or GC sites
  # are present (due to amplicon boundaries). Instead, CG/GC sites will be expanded in a strand
  # specific way and the overlaps combined

  # get an indicator variable for the strand of the gr
  plusIndex<-strand(meth_CGs_gr)=="+"
  # extend the GR by 1 in different directions, depending on the strand
  end(meth_CGs_gr[plusIndex])<-end(meth_CGs_gr[plusIndex])+1
  start(meth_CGs_gr[!plusIndex])<-start(meth_CGs_gr[!plusIndex])-1
  # find 2bp overlaps between gr on different strands.
  ol<-findOverlaps(meth_CGs_gr[plusIndex],meth_CGs_gr[!plusIndex], minoverlap=2,ignore.strand=TRUE)
  # sum the methylation counts from both strands for a given CG site.
  meth_CGsCol_gr<-meth_CGs_gr[plusIndex][queryHits(ol)]
  strand(meth_CGsCol_gr)<-"*"
  values(meth_CGsCol_gr)<-as.matrix(values(meth_CGs_gr[plusIndex][queryHits(ol)]))+as.matrix(values(meth_CGs_gr[!plusIndex][subjectHits(ol)]))
  # we could add back the half CG sites, but as they are at the edge, they probably are distorted by the primer
  # seq. (in two cases i checked, all columns had 0 methylation). The command to add them back would be:
  # meth_CGsCol_gr<-c(meth_CGsCol_gr,meth_CGs_gr[plusIndex][!(seq(length(meth_CGs_gr[plusIndex])) %in% queryHits(ol))],
  #                  meth_CGs_gr[!plusIndex][!(seq(length(meth_CGs_gr[!plusIndex])) %in% subjectHits(ol))])


  # do the same to the GC matrix
  # get an indicator variable for the strand of the gr
  plusIndex<-strand(meth_GCs_gr)=="+"
  # extend the GR by 1 in different directions, depending on the strand (extends opposite way from CG)
  start(meth_GCs_gr[plusIndex])<-start(meth_GCs_gr[plusIndex])-1
  end(meth_GCs_gr[!plusIndex])<-end(meth_GCs_gr[!plusIndex])+1
  # find 2bp overlaps between gr on different strands.
  ol<-findOverlaps(meth_GCs_gr[plusIndex],meth_GCs_gr[!plusIndex], minoverlap=2,ignore.strand=TRUE)
  # sum the methylation counts from both strands for a given CG site.
  meth_GCsCol_gr<-meth_GCs_gr[plusIndex][queryHits(ol)]
  strand(meth_GCsCol_gr)<-"*"
  values(meth_GCsCol_gr)<-as.matrix(values(meth_GCs_gr[plusIndex][queryHits(ol)]))+as.matrix(values(meth_GCs_gr[!plusIndex][subjectHits(ol)]))


  #####################
  #filter for coverage
  #####################

  # identify columns that have total and methylated C counts
  Tcounts=grep('_T\\>',colnames(elementMetadata(meth_CGsCol_gr)))
  Mcounts=grep('_M\\>',colnames(elementMetadata(meth_CGsCol_gr)))

  ######
  #CGs
  ######

  # calculate frequency of methyltion
  CG.met.mat<-as.matrix(mcols(meth_CGsCol_gr))
  # filter for coverage (only positions that are covered by >c0 reads)
  CovFilter<-as.matrix(mcols(meth_CGsCol_gr)[,Tcounts])>c0
  # put NA into positions where coverage is <c0
  CG.met.mat[,Tcounts][!CovFilter]<-NA
  CG.met.mat[,Mcounts][!CovFilter]<-NA
  #bind the GRanges with the scores
  CG.met<-meth_CGsCol_gr
  mcols(CG.met)<-CG.met.mat


  ######
  #GCs
  ######

  # calculate frequency of methyltion
  GC.met.mat<-as.matrix(mcols(meth_GCsCol_gr))
  # filter for coverage (only positions that are covered by >c0 reads)
  CovFilter<-as.matrix(mcols(meth_GCsCol_gr)[,Tcounts])>c0
  # put NA into positions where coverage is <c0
  GC.met.mat[,Tcounts][!CovFilter]<-NA
  GC.met.mat[,Mcounts][!CovFilter]<-NA
  #bind the GRanges with the scores
  GC.met<-meth_GCsCol_gr
  mcols(GC.met)<-GC.met.mat

  ###########################
  #create an unified object
  ###########################

  # need to label GCGs that will appear in both matrices as they form GCH and HCG
  # this is to avoid double counting of GCGs later on when CpGs and GpCs are combined
  oGCG=as.matrix(findOverlaps(resize(GC.met,1,fix='end'),resize(CG.met,1,fix='start'),type='equal')) # do i miss CGCs?
  GC.met$type='GCH'
  CG.met$type='CGH'
  GC.met$type[oGCG[,1]]='GCG'
  CG.met$type[oGCG[,2]]='GCG'

  # extract different types of sites
  GCG.met<-sort(c(CG.met[CG.met$type=="GCG"],GC.met[GC.met$type=="GCG"]))
  CG.met<-sort(CG.met[CG.met$type!="GCG"])
  GC.met<-sort(GC.met[GC.met$type!="GCG"])


  # combine both GR objects into a single list
  umet=sort(c(CG.met,GC.met,GCG.met))
  return(umet)
}






################################################################
######## functions for extracting C meth matrices  #############
################################################################

getCmethMatrix<-function(proj,gr,sampleName){
  # function returns methylation matrix for a genomic range (gr) for a particular
  # sample. rows of the matrix are individual reads, columns are methylation positions
  # input is a QuasR project, a GRanges object and name of the sample for the qProject
  # Gets every single C, not just CpG/GpC
  Cs=qMeth(proj, query=gr,mode="allC",reportLevel="alignment")
  # use data.table to get a 1,0 matrix of methylation profiles
  out<-if (length(Cs[[sampleName]]$Cid)==0) {
    print("no methylation")
    NA
  } else {
    print("processing mat")
    all.cids=unique(Cs[[sampleName]]$Cid) # get all possible C locations
  # make the data.table object with methylation status (meth), read name (aid), and
  # C genome position (Cid)
    dt=data.table(meth=Cs[[sampleName]]$meth ,aid=Cs[[sampleName]]$aid ,cid=Cs[[sampleName]]$Cid)
  # spread positions into columns and make a matrix
    dtm<-spread(dt,key=cid,value=meth)
    CpGm<-as.matrix(dtm[,2:dim(dtm)[2]])
    row.names(CpGm)<-dtm$aid
    CpGm
  }
  out
}


fuseReadMat <- function(resp,resm){ #assumes that reads are no on different strands
  uReads=unique(c(rownames(resp),rownames(resm)))
  uCs=unique(c(colnames(resp),colnames(resm)))
  matCGo=matrix(nrow=length(uReads),ncol=length(uCs))
  colnames(matCGo)=uCs
  rownames(matCGo)=uReads
  matCGo[rownames(resp),colnames(resp)]=resp #fill up positive reads
  #add negative reads
  NAi=apply(resm,1,function(x){sum(is.na(x))==length(x)})
  matCGo[rownames(resm[!NAi,]),colnames(resm[!NAi,])]=resm[!NAi,]
  matCGo
}




getIsolatedSites<-function(gr1,gr2) {
  # removes from gr1 and sites which overlap with gr2
  ol<-findOverlaps(gr1,gr2,ignore.strand=T)
  if (length(ol)>0) {
    isolated<-gr1[-queryHits(ol)]
  } else {
    isolated<-gr1
  }
  return(isolated)
}


getOverlappingSites<-function(gr1,gr2) {
  # adds any sites from gr1 that overlap with gr2 to gr2
  ol<-findOverlaps(gr1,gr2,ignore.strand=T)
  if (length(ol)>0) {
    overlapping<-c(gr2,gr1[queryHits(ol)])
  } else {
    overlapping<-gr2
  }
  return(overlapping)
}


splitRangeInPairs<-function(gr) {
  starts<-seq(start(gr),end(gr),by=2)
  newGR<-GRanges(seqnames=seqnames(gr),ranges=IRanges(start=starts,width=2),strand="*")
  return(newGR)
}


splitGCGCruns<-function(gr) {
  # function to (arbitrarily) decompose overlapping GCs and CGs to triplets and doublets
  # create new gr for split up ranges
  olGCGC<-GRanges()
  # find ranges that do not need splitting
  width3<-width(gr)==3
  olGCGC<-append(olGCGC,gr[width3])
  # find ranges that need splitting
  toSplit<-gr[!width3]
  if (length(toSplit)>0) {
    for (i in 1:length(toSplit)) {
      isOdd<-width(toSplit[i])%%2==1
      if (isOdd) {
        # if there are an odd number of bases, chop off the first three as a triplet, and the rest in pairs
        firstTriplet<-resize(toSplit[i],3,"start")
        evenRange<-resize(toSplit[i],width(toSplit[i])-3,"end")
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


destrandReads<-function(allCGGCgr,CposGR,matC) {
  ##### not sure if destranding works. As i have no real data to test it on ######
  # function to deal with counts coming from adjacent sites on
  # complementary strands. both strands values within a 2bp motif
  # are shifted to the position of the first base in the motif.
  # with directional sequencing this is not a problem as the reads
  # from the other strand are empty.

  # we will ignore third position in triplet as it is hard to figure out how to
  #distribute its values between F and R reds
  triplets<-which(width(allCGGCgr)==3)
  #### NOTE: test correlation between positions in triplets when i have stranded data
  allCGGCgrPos1<-resize(allCGGCgr[triplets],1,"start")
  allCGGCgrPos2<-resize(allCGGCgr[triplets],1,"center")
  allCGGCgrPos3<-resize(allCGGCgr[triplets],1,"end")
  allCGGCgr[triplets]<-resize(allCGGCgr[triplets],2,"start")

  # create stranded version of CPosGR
  CposGRstranded<-CposGR
  strand(CposGRstranded)<-ifelse(as.character(getSeq(genome,CposGR))=="C","+","-")

  #####
  # get vectors for triplets to look at correlation
  ol1<-findOverlaps(CposGR,allCGGCgrPos1)
  ol2<-findOverlaps(CposGR,allCGGCgrPos2)
  ol3<-findOverlaps(CposGR,allCGGCgrPos3)
  p1<-as.character(start(CposGR[queryHits(ol1)]))
  p2<-as.character(start(CposGR[queryHits(ol2)]))
  p3<-as.character(start(CposGR[queryHits(ol3)]))
  matp1<-matC[,p1,drop=FALSE]
  matp2<-matC[,p2,drop=FALSE]
  matp3<-matC[,p3,drop=FALSE]
  tripletDF<-data.frame("p1"=as.vector(matp1),"p2"=as.vector(matp2),"p3"=as.vector(matp3))
  ####


  #get vector for which Cs are on the positive strand
  CposPlusi<-as.character(strand(CposGRstranded))=="+"

  # change gr in CposGR  to destrand. Always taking first position in motif as common location
  ol<-findOverlaps(CposGR[CposPlusi],allCGGCgr)
  start(CposGR[CposPlusi][queryHits(ol)])<-start(allCGGCgr[subjectHits(ol)])
  end(CposGR[CposPlusi][queryHits(ol)])<-start(allCGGCgr[subjectHits(ol)])

  #same for the reads on the negative strand
  ol<-findOverlaps(CposGR[!CposPlusi],allCGGCgr)
  start(CposGR[!CposPlusi][queryHits(ol)])<-start(allCGGCgr[subjectHits(ol)])
  end(CposGR[!CposPlusi][queryHits(ol)])<-start(allCGGCgr[subjectHits(ol)])

  #take all sites on +ive strand in the matrix and remove empty rows
  matCplus<-matC[,CposPlusi,drop=FALSE]
  matCplus<-matCplus[rowSums(!is.na(matCplus))>0,,drop=FALSE]
  colnames(matCplus)<-as.character(start(CposGR)[CposPlusi])

  #take all sites on -ive strand and remove empty rows
  matCminus<-matC[,!CposPlusi,drop=FALSE]
  matCminus<-matCminus[rowSums(!is.na(matCminus))>0,,drop=FALSE]
  colnames(matCminus)<-as.character(start(CposGR)[!CposPlusi])

  # combine the matrices again
  newColNames<-as.character(sort(as.numeric(unique(c(colnames(matCplus),colnames(matCminus))))))
  newMatC<-matrix(nrow=dim(matCplus)[1]+dim(matCminus)[1],ncol=length(newColNames))
  colnames(newMatC)<-newColNames

  # add in the matCplus
  idx<-match(colnames(matCplus),newColNames)
  if(length(idx)>0) {
    newMatC[1:dim(matCplus)[1],idx]<-matCplus
  }
  # add in the matCminus
  idx<-match(colnames(matCminus),newColNames)
  if(length(idx)>0) {
  newMatC[(dim(matCplus)[1]+1):dim(newMatC)[1],idx]<-matCminus
  }
  # add the read names
  rownames(newMatC)<-c(rownames(matCplus),rownames(matCminus))

  # remove redundent GR from CposGR
  CposGR<-disjoin(CposGR)
  # double check the names are the same
  stopifnot(sum(as.character(start(CposGR))!=colnames(newMatC))==0)
  return(list(CposGR,newMatC,tripletDF))
}




getGCmatrix1<-function(matList, ampliconGR, genome=Celegans, conv.rate=80, destrand=TRUE, sampleName="",plotData=F){
  # this funciton produces CG and GC matrices for individual amplicons as well as all triplet GCG/CGC and overlapping
  # CGs and GCs (in GCGC objects). It plots histograms of the conversion rates of different reads for each amplicon,
  # and produces a table with various data to get an overview of how the matrices are being processed

  # This function assumes the ampliconGR has the name of the amplicon in the first column of the metadata.

  if (plotData==T) {
    # create data frame to store some numbers
    GCCGdata<-data.frame(amplicons=ampliconGR$ID, chr=as.vector(seqnames(ampliconGR)), ampliconOri=strand(ampliconGR),
                       totalReads=0, gt80convRate=0, GCs=0, CGs=0, GCGCs=0)

    # open PDF for histograms of C conversion rate in reads
    pdf(paste0(path,"/plots/conversionRateCutoff_",sampleName,".pdf"),paper="a4",height=8,width=8)
    par(mfrow=c(2,2))
  }

  # scroll through matList an create CG GC and GCGC matrices
  matGClist<-list()
  for (i in seq_along(matList)) {
    # get matrix for individual amplicon
    matC<-matList[[i]]

    # extract chr info from amplicon genomic ranges object
    ampliconName=names(matList)[i]
    stopifnot(mcols(ampliconGR)[i,"ID"]==ampliconName) # ensure the ampliconGR and matrix list are in the same order
    chr<-as.character(seqnames(ampliconGR)[i])

    # stop prossing if matC is not a matrix (no reads for this amplicon)
    if (!is.matrix(matC)) {
      next()
    }
    print(paste(i,"processing matrices"))
    ########
    # make non-overlapping GRanges for GCs CGs and longer triplets or GCGC stretches
    ########

    # first get all C positions in the matrix
    Cpos=as.numeric(colnames(matC))
    CposGR=GRanges(rep(chr,length(Cpos)),IRanges(Cpos,Cpos))
    # also create stranded version
    CposGRstranded<-CposGR
    strand(CposGRstranded)<-ifelse(as.character(getSeq(genome,CposGR))=="C","+","-")

    # get squence context around C to find CGs and GCs and GCG/CGC triplets
    cPosSeq=getSeq(genome,resize(CposGR,3,fix='center'))

    # create genomic ranges for the different types of sites
    GCgr=resize(CposGRstranded[vcountPattern('GC',cPosSeq)==1],2,fix='end')
    CGgr=resize(CposGRstranded[vcountPattern('CG',cPosSeq)==1],2,fix='start')
    GCGgr=resize(CposGRstranded[vcountPattern('GCG',cPosSeq)==1],3,fix='center')
    CGCgr=resize(CposGRstranded[vcountPattern('CGC',cPosSeq)==1],3,fix='center')
    GCGCgr<-c(GCGgr,CGCgr) # both strands

    # move any GC sites that overlap with GCG to the GCG list (removing them from GC list)
    GCGCgr<-getOverlappingSites(GCgr,GCGCgr) # always get overlapping first, otherwise will lose them
    GCgr<-reduce(getIsolatedSites(GCgr,GCGCgr),ignore.strand=T) # collapse GR on opposite strands

    # move any CG sites that overlap with GCG to the GCG list (removing them from CG list)
    GCGCgr<-getOverlappingSites(CGgr,GCGCgr) # always get overlapping first, otherwise will lose them
    CGgr<-reduce(getIsolatedSites(CGgr,GCGCgr),ignore.strand=T) # collapse GR on opposite strands

    # reduce all the overlapping set to a smallest merged set
    GCGCgr<-reduce(GCGCgr,ignore.strand=T)
    # then arbitrarily split them into non-overlapping GR 2-3 bp long
    GCGCgr<-splitGCGCruns(GCGCgr)

    if (plotData==T) {
      #save data about number of sites
      GCCGdata[i,"totalReads"]<-dim(matC)[1]
      GCCGdata[i,"GCs"]<-length(GCgr)
      GCCGdata[i,"CGs"]<-length(CGgr)
      GCCGdata[i,"GCGCs"]<-length(GCGCgr)
    }
    ########
    # destrand if sequencing is non-directional (reads from both strands)
    ########

    if (destrand==T) {
      # make unified object with all independant GR
      allCGGCgr<-c(GCgr,CGgr,GCGCgr)
      # when reads come from both strands, it is necessary to shift meth score by one bp to align
      # the methylation per site on F and R reads (otherwise the matrix is full of NAs)
      res<-destrandReads(allCGGCgr,CposGR,matC)
      CposGR<-res[[1]]
      matC<-res[[2]]
      tripletDF<-res[[3]]
      tp<-tripletDF[complete.cases(tripletDF[,c(1,3)]),c(1,3)]
      print(paste("fraction pos 1&3 that differ:",sum(tp$p1!=tp$p3)/dim(tp)[1]))
    }

    ########
    # filter reads based on non-CG-GC C->T conversion rate
    ########

    # make indeces for the various subtypes
    GCi<-queryHits(findOverlaps(CposGR,GCgr))
    CGi<-queryHits(findOverlaps(CposGR,CGgr))
    if (destrand==T) {
      GCGCi<-queryHits(findOverlaps(CposGR,resize(GCGCgr,2,"start"))) # ignore third position in triplet
    } else {
      GCGCi<-queryHits(findOverlaps(CposGR,GCGCgr))
    }

    # get index of all non-CG/GC Cs to calculate conversion rate.
    mcols(CposGR)$type="C"
    mcols(CposGR)$type[GCi]<-"GC"
    mcols(CposGR)$type[CGi]<-"CG"
    mcols(CposGR)$type[GCGCi]<-"GCGC"
    Ci=which(mcols(CposGR)$type=="C")

    # filter based on conversion
    convRate=100-(rowMeans(matC[,Ci,drop=FALSE],na.rm=T)*100)
    Convi=convRate>conv.rate

    if(plotData==T){
      GCCGdata[i,"gt80convRate"]<-sum(Convi) #number of reads with convRate > 80%

      #plot histograms of conversion rate
      p<-hist(convRate, breaks=50, col="grey", xlim=c(0,100), main=paste0(mcols(ampliconGR)[i,1],
                                " (amplicon strand: ",ifelse(strand(ampliconGR[i])=="+","pos","neg"),")"))
      abline(v=80,lty=2,ylim=c(0,max(p$counts)),col="light grey")
    }
    ########
    # make CG GC matrices
    ########

    # get matrices
    matGC=matC[Convi,GCi,drop=FALSE]
    matCG=matC[Convi,CGi,drop=FALSE]
    matGCGC=matC[Convi,GCGCi,drop=FALSE]

    # Remove the empty columns
    matGC=matGC[,colSums(!is.na(matGC))>0,drop=FALSE] # remove no coverage columns
    matCG=matCG[,colSums(!is.na(matCG))>0,drop=FALSE]
    matGCGC=matGCGC[,colSums(!is.na(matGCGC))>0,drop=FALSE]
    # Remove the empty rows
    matGC=matGC[rowSums(!is.na(matGC))>0,,drop=FALSE] # remove no coverage rows
    matCG=matCG[rowSums(!is.na(matCG))>0,,drop=FALSE]
    matGCGC=matGCGC[rowSums(!is.na(matGCGC))>0,,drop=FALSE]

    # add matrices as list to list of amplicons
    matGClist[[ampliconName]]<-list(matGC=matGC,matCG=matCG,matGCGC=matGCGC)
  }
  if (plotData==T) {
    dev.off() # close histogram plot
    write.csv(GCCGdata,paste0(path,"/csv/GCCGdata_",sampleName,".csv"))
  }
  return(matGClist)
}






findCGGCinGenome<-function(genome){
  CGgr<-vmatchPattern("CG",genome)
  GCgr<-vmatchPattern("GC",genome)
  return(list(GCgenome=GCgr,CGgenome=CGgr))
}





mergeGC_CGmats<-function(matList){
  mergedMatList<-lapply(seq_along(matList),function(x) {
    mats<-matList[[x]]
    matCG=mats$matCG
    matGC=mats$matGC
    matGCGC=mats$matGCGC
    uReads=unique(c(rownames(matCG),rownames(matGC),rownames(matGCGC)))
    uCs=sort(unique(c(colnames(matCG),colnames(matGC),colnames(matGCGC))))
    matCGo=matrix(nrow=length(uReads),ncol=length(uCs))
    colnames(matCGo)=uCs
    rownames(matCGo)=uReads
    matCGo[rownames(matCG),colnames(matCG)]<-matCG#get CGs in
    matCGo[rownames(matGC),colnames(matGC)]<-matGC#get GCs in
    matCGo[rownames(matGCGC),colnames(matGCGC)]<-matGCGC#get GCGCs in
    matCGo
  })
  names(mergedMatList)<-names(matList)
  return(mergedMatList)
}


getRelativeCoord<-function(mat,regionGR,invert=F){
  # converts matrix from absolute genome coordinates to
  # relative coordinates within a genomic Range
  pos<-as.numeric(colnames(mat))
  regionStart<-start(regionGR)
  regionEnd<-end(regionGR)
  if (invert==F) {
    newPos<-pos-regionStart
    colnames(mat)<-as.character(newPos)
  } else {
    newPos<-regionEnd-pos
    colnames(mat)<-newPos
    mat<-mat[,order(colnames(mat)),drop=F]
    colnames(mat)<-as.character(colnames(mat))
  }
  return(mat)
}


changeAnchorCoord<-function(mat,anchorCoord=0) {
  # changes the 0 coordinate position (anchorCoord) of
  # a matrix. e.g. sets position 250 to 0 in a 500 region around TSS
  # to get +-250 bp around TSS
  pos<-as.numeric(colnames(mat))
  newPos<-pos-anchorCoord
  colnames(mat)<-as.character(newPos)
  return(mat)
}


getRelativeCoordMats<-function(matList,grs,anchorCoord=0) {
  newMatList<-lapply(seq_along(matList),function(x){
    print(x)
    mat<-matList[[x]]
    if(sum(dim(mat)==c(0,0))<1) {
      regionID<-names(matList)[x]
      regionGR<-grs[grs$ID==regionID]
      newMat<-getRelativeCoord(mat,regionGR,invert=ifelse(strand(regionGR)=="+",F,T))
      newMat<-changeAnchorCoord(mat=newMat,anchorCoord=anchorCoord)
      newMat
    } else {
      newMat<-mat
    }
  })
  names(newMatList)<-names(matList)
  return(newMatList)
}



getMetaMethFreq<-function(matList,regionGRs,minReads=50) {
  for (i in seq_along(matList)) {
    if (dim(matList[[i]])[1]>minReads) {
      vecSummary<-colMeans(matList[[i]],na.rm=T)
      df<-data.frame("position"=names(vecSummary),"methFreq"=vecSummary)
      df$ID<-names(matList)[i]
      df$chr<-as.character(seqnames(regionGRs)[match(names(matList)[i],regionGRs$ID)])
      if (exists("methFreqDF")){
        methFreqDF<-rbind(methFreqDF,df)
      } else {
        methFreqDF<-df
      }
    }
  }
  return(methFreqDF)
}



#################### plotting functions #########################

grToIgv<-function(gr,sampleName,fileName){
  # this function converts data from GR to track for IGV
  # The frequencies are converted to rounded counts from 100
  scores=100-round((mcols(gr)[,sampleName]*100))
  NAi=is.na(scores) #remove uncovered
  grf=gr[!NAi]
  df <- data.frame(Chromosome=seqnames(grf),
                   Start=start(grf)-1,
                   End=end(grf),
                   Feature=c(rep(sampleName, length(grf))),
                   R1=scores[!NAi])
  write.table(df, file=fileName, quote=F, sep="\t", row.names=F, col.names=F)
}


grToBw<-function(gr,sampleName,fileName){
  # this function converts data from GR to track for IGV
  # The frequencies are converted to rounded counts from 100
  newGR<-gr
  newGR$score=1-mcols(gr)[,sampleName]
  NAi=is.na(mcols(newGR)$score) #remove uncovered
  newGR=newGR[!NAi]
  rtracklayer::export.bw(newGR, con=fileName)
}


# plotBinaryMatrix<-function(mat) {
#   mat[is.na(mat)]<- -1
#   mm<-as.data.frame(mat)
#   mm$readName<-1:dim(mm)[1]
#   mmv<-gather(mm,key="Cposition",value="methylation",-"readName")
#   ggplot(mmv,aes(x=Cposition,y=readName)) +
#     geom_tile(aes(fill=methylation)) +
#     scale_fill_gradient(low="white",high="steelblue")
# }



plotAmpliconCGGC<-function(matClist,matGClist,sampleName="",ampliconGR) {
  # plot the position of informative Cs in the CG and GC matrices. Also plot position of
  # GCGs (excluded?)
  pdf(paste0(path,"/plots/GCCGloss_",sampleName,".pdf"),paper="a4r",height=8,width=11)
  for (ampliconName in names(matGClist)) {
    print(paste("plotting",ampliconName))

    #extract C matrix from list
    matC<-matClist[[ampliconName]]

    # extract GC CG and GCGC matrices
    GCCGmats<-matGClist[[ampliconName]]
    matGC<-GCCGmats[["matGC"]]
    matCG<-GCCGmats[["matCG"]]
    matGCGC<-GCCGmats[["matGCGC"]]

    # only look at reads that are in Convi filtered matrices
    validReads<-row.names(matGC)
    idx<-match(validReads,row.names(matC))
    matC<-matC[idx,,drop=FALSE]

    #extract Cs from matC that have some % methylation
    informativeCs<-matC[,colSums(!is.na(matC[,,drop=FALSE]))>0,drop=FALSE]
    if (sum(dim(informativeCs))>0) {
      plotTable<-data.frame(position=as.numeric(colnames(informativeCs)),
                            methylatedCs=as.vector(colSums(informativeCs,na.rm=T)),
                            methylatedGCs=NA, methylatedCGs=NA,methylatedGCGs=NA)

      #filter GC positions for plotting
      matGC<-matGC[,(as.numeric(as.vector(colnames(matGC))) %in% plotTable$position),drop=F] # remove positions not in informativeCs
      idx<-plotTable$position %in% as.numeric(as.vector(colnames(matGC))) # make sure at least 1 C is in the table
      if(sum(idx)>0) { plotTable[idx,"methylatedGCs"]<-as.vector(colSums(matGC,na.rm=T))}

      #filter CG positions for plotting
      matCG<-matCG[,(as.numeric(as.vector(colnames(matCG))) %in% plotTable$position),drop=F]
      idx<-plotTable$position %in% as.numeric(as.vector(colnames(matCG)))
      if(sum(idx)>0) {plotTable[idx,"methylatedCGs"]<-as.vector(colSums(matCG,na.rm=T))-0.1}

      #filter GCG positions for plotting
      matGCGC<-matGCGC[,(as.numeric(as.vector(colnames(matGCGC))) %in% plotTable$position),drop=F]
      idx<-plotTable$position %in% as.numeric(as.vector(colnames(matGCGC)))
      if(sum(idx)>0) {plotTable[idx,"methylatedGCGs"]<-as.vector(colSums(matGCGC,na.rm=T))-0.2}

      # rearrange data for plotting
      plotTableLong<-gather(plotTable,key="siteType",value="methCount",c("methylatedGCs","methylatedCGs","methylatedGCGs"))
      if(sum(!is.na(plotTableLong$methCount))==0) {plotTableLong$methCount<-0}

      p<-ggplot(plotTableLong,aes(x=position,y=methylatedCs/3))+geom_bar(stat="identity") +
        geom_point(aes(y=methCount,col=siteType),alpha=0.5,size=2) +
        scale_color_manual(values=c("methylatedGCs"="blue","methylatedCGs"="red","methylatedGCGs"="darkolivegreen4")) +
        ggtitle(paste0(ampliconName," (amplicon strand: ",
                       ifelse(strand(ampliconGR[match(ampliconName,ampliconGR$ID)])=="+","pos","neg"),")"))
      print(p)
    }
  }
  dev.off()
}





plotSingleMoleculesAmp<-function(mat,regionName,regionGRs,featureGRs,myXlab="CpG/GpC position",featureLabel="TSS",title=NULL,
                                 baseFontSize=12) {
  ### single molecule plot. mat is matrix containing methylation values at different postions (columns) in
  # individual reads (rows). regionName is the ID of the amplicon or genomic regoin being plotted. regionGRs is a
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
    d<-gather(df,key=position,value=methylation)
    d$molecules<-seq_along(reads)
    d$methylation<-as.character(d$methylation)
    d$position<-as.numeric(d$position)
    if (is.null(title)) {
      title=paste0(regionName, ": ",seqnames(regionGR)," ",strand(regionGR),"ve strand")
    }
    p<-ggplot(d,aes(x=position,y=molecules,width=2)) +
      geom_tile(aes(width=6,fill=methylation),alpha=0.8) +
      scale_fill_manual(values=c("0"="black","1"="grey80"),na.translate=F,na.value="white",
                        labels=c("protected","accessible"),name="dSMF") + theme_light(base_size=baseFontSize) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold",hjust = 0.5),
            legend.position="bottom", legend.key.height = unit(0.1, "cm")) +
      ggtitle(title) +
      xlab(myXlab) + ylab("Single molecules") + xlim(start(regionGR),end(regionGR)+10)
    if(length(featureGRs)>0) {
      p<-p+geom_linerange(aes(x=start(featGR), y=NULL, ymin=0, ymax=length(reads)+max(3,0.04*length(reads))),col="red") +
        ggplot2::annotate("segment", x = start(featGR), xend = start(featGR)+20*ifelse(strand(featGR)=="-",-1,1),
                    y = length(reads)+max(3,0.04*length(reads)), yend =length(reads)+max(3,0.04*length(reads)), colour = "red",
                    arrow=arrow(length = unit(0.3, "cm")), size=0.7) +
        ggplot2::annotate(geom="text", x=start(featGR), y=-max(2,0.03*length(reads)), label=featureLabel,color="red")

    }
  } else {
    p<-NULL
  }
  return(p)
}





plotSingleMoleculesWithAvr<-function(mat,regionName,regionGRs,featureGRs,myXlab="CpG/GpC position",featureLabel="TSS",
                                     title=NULL, baseFontSize=11) {
  if(dim(mat)[1]>10) {
    regionGR<-regionGRs[match(regionName,regionGRs$ID)]
    if (length(featureGRs)>0) {
    featGR<-featureGRs[match(regionName,featureGRs$ID)]
    }
    na.matrix<-is.na(mat)
    mat[na.matrix]<--1
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
    d<-gather(df,key=position,value=methylation)
    d$molecules<-seq_along(reads)
    d$methylation<-as.character(d$methylation)
    d$position<-as.numeric(d$position)
    if (is.null(title)) {
      title=paste0(regionName, ": ",seqnames(regionGR)," ",strand(regionGR),"ve strand")
    }
    dAvr<-data.frame(position=as.numeric(colnames(df)),dSMF=1-colSums(df,na.rm=T)/length(reads))

    p1<-ggplot(dAvr,aes(x=position,y=dSMF,group=1)) + geom_point()+
      geom_line(size=1,show.legend=F) + guides(fill=FALSE, color=FALSE) +
      theme_light(base_size=baseFontSize) + ylab("Mean dSMF") +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
      ylim(0,1) + xlim(start(regionGR),end(regionGR)+20)
    if (length(featureGRs)>0) { # plot feature if present
      p1<-p1 + geom_linerange(aes(x=start(featGR), y=NULL, ymin=0, ymax=1),col="red",size=0.7) +
        ggplot2::annotate("segment", x = start(featGR), xend = start(featGR)+20*ifelse(strand(featGR)=="-",-1,1),
               y = 1, yend = 1, colour = "red", size=0.7, arrow=arrow(length = unit(0.2, "cm")))
    }

    p2<-ggplot(d,aes(x=position,y=molecules,width=2)) +
      geom_tile(aes(width=6,fill=methylation),alpha=0.8) +
      scale_fill_manual(values=c("0"="black","1"="grey80"),na.translate=F,na.value="white",
                        labels=c("protected","accessible"),name="dSMF") + theme_light(base_size=baseFontSize) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_blank(),
            legend.position="bottom", legend.key.height = unit(0.1, "cm")) +
      xlab(myXlab) + ylab("Single molecules") + xlim(start(regionGR),end(regionGR)+20)
    if (length(featureGRs)>0) {
      p2<-p2+geom_linerange(aes(x=start(featGR), y=NULL, ymin=0, ymax=length(reads)+max(3,0.04*length(reads))),col="red") +
        ggplot2::annotate("segment", x = start(featGR), xend = start(featGR)+20*ifelse(strand(featGR)=="-",-1,1),
               y = length(reads)+max(3,0.04*length(reads)), yend =length(reads)+max(3,0.04*length(reads)),
               arrow=arrow(length = unit(0.2, "cm")), colour="red", size=0.7) +
        ggplot2::annotate(geom="text", x=start(featGR), y=-max(2,0.03*length(reads)), label=featureLabel, color="red")
    }
    figure<-ggarrange(p1, p2, heights = c(0.5, 2),
                      ncol = 1, nrow = 2, align = "v")
    figure<-annotate_figure(figure, top = text_grob(title, face = "bold"))
  } else {
    figure=NULL
  }
  return(figure)
}

