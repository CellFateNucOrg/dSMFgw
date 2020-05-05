# settings of some variables required by dSMFseqAnalysis1_amp_aln.R and dSMFseqAnalysis2_amp_XvA.R

# sourcing variables from varSettings.sh
# number of threads you wish to use
#threadNum=Sys.getenv("threadNum")
getUnixVar<-function(varName){
	varLine<-system(paste0("grep ^",varName," varSettings.sh"),intern=TRUE)
	varVal<-gsub(paste0("^",varName,"\\s*=\\s*"),"",varLine)
	varVal<-gsub("#.*$","",varVal)
	varVal<-gsub("\\s*$","",varVal)	
	return(varVal)
}


dataType=getUnixVar("dataType")
genomeFile=getUnixVar("genomefile")
seqDate=getUnixVar("seqDate")
expName=getUnixVar("expName")


# default path is the project working directory.
# Note that rawData should be one level up from this location (../rawData)
path<-getwd()

# labels for main biological comparison you wish to make (labels should be included at the start of the
# sample name in the settings file)
#testGroup<-c("N2")
#extract testgroups from varSettings.sh file
testGroup<-system("grep testGroups varSettings.sh",intern=TRUE)
testGroup<-gsub("^testGroups=\\s*\\(\\s*","",testGroup)
testGroup<-gsub("\\s*\\).*","",testGroup)
testGroup<-unique(unlist(strsplit(testGroup,"\\s+")))


#load main genome to which you wish to align sequences and store it under the geneirc "genome" variable
library("BSgenome.Celegans.UCSC.ce11")
genome<-Celegans

#genomeVer="WS235"
genomeVer<-getUnixVar("genomeVer")

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

if (genomeVer=="WS235") {
  genome<-ucscToWbGR(genome)
}

#genomeFile<-"/home/ubelix/izb/semple/genomeVer/ws260/sequence/c_elegans.PRJNA13758.WS260.genomic.fa"


#load secondary genome to which you wish to align sequences and store it under the geneirc "auxGenome" variable
#library("BSgenome.Ecoli.NCBI.20080805")
#auxGenomeName<-c("Escherichia_coli_ensembl", "PhiX174_NC_001422.1", "LambdaPhage_NC_001416.1")
#auxGenomeFile<-c("~/genomeVer/ecoli/Escherichia_coli.HUSEC2011CHR1.dna.fa", 
#        "~/genomeVer/phiX/phiX.fasta", "~/genomeVer/lambda/lambdaPhage.fasta") 


# file with genomic ranges for amplicons. Must contain a metadata column called "ID" with a unique name for
# each amplicon (e.g. gene name)
amplicons<-readRDS('/home/ubelix/izb/semple/myData/usefulObjects/ampliconGR.RDS')
if (genomeVer=="WS235") {
  amplicons<-ucscToWbGR(amplicons)
}
names(mcols(amplicons))[1]<-"ID"

# file with genomic ranges for TSS (or other genomic feature). Must contain a metadata column called "ID" with a unique name for
# each TSS (e.g. gene name). This ID should be the same as for the amplicons
ampTSS<-readRDS('/home/ubelix/izb/semple/myData/usefulObjects/ampliconMaxTSSgr.RDS')
if (genomeVer=="WS235") {
  ampTSS<-ucscToWbGR(ampTSS)
}
names(mcols(ampTSS))[1]<-"ID"

#files with TSS for gw alignments
#highConfTSS where all three datasets agree on the TSS with the maxTSS (872 genes)
highConfTSS<-readRDS('/home/ubelix/izb/semple/genomeVer/ws260/rds/ChenKreusSaitoTSS_highConf_872.RDS')
highConfTSS$ID<-names(highConfTSS)
if (genomeVer!="WS235") {
    seqlevels(highConfTSS)<-paste0("chr",seqlevels(highConfTSS))
}
#lessConfTSS which are maxTSS from combining the three TSS datasets but were not identical (1955 genes)
lessConfTSS<-readRDS('/home/ubelix/izb/semple/genomeVer/ws260/rds/ChenKreusSaitoTSS_lessConf_1955.RDS')
lessConfTSS$ID<-names(lessConfTSS)
if (genomeVer!="WS235") {
    seqlevels(lessConfTSS)<-paste0("chr",seqlevels(lessConfTSS))
}


