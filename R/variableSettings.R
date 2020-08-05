# settings of some variables required by dSMFseqAnalysis1_amp_aln.R and dSMFseqAnalysis2_amp_XvA.R


print("importing unix variables")
# sourcing variables from varSettings.sh
# number of threads you wish to use
#threadNum=Sys.getenv("threadNum")
getUnixVar<-function(varName){
	varLine<-system(paste0("grep ^",varName,"= varSettings.sh"),intern=TRUE)
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


getUnixArray<-function(varName){
        varLine=system(paste0("grep ^",varName,"= varSettings.sh"),intern=TRUE)
        varLine=gsub(paste0("^",varName,"=\\s*\\(\\s*"),"",varLine) # clean up beginning of line
        varLine<-gsub("\\s*\\).*","",varLine) # clean up end of line
        varLine<-unique(unlist(strsplit(varLine,"\\s+"))) # make into a list
        return(varLine)
}

# labels for main biological comparison you wish to make (labels should be included at the start of the
# sample name in the settings file)
#testGroup<-c("N2")
#extract testgroups from varSettings.sh file
#testGroup<-system("grep testGroups varSettings.sh",intern=TRUE)
#testGroup<-gsub("^testGroups=\\s*\\(\\s*","",testGroup)
#testGroup<-gsub("\\s*\\).*","",testGroup)
#testGroup<-unique(unlist(strsplit(testGroup,"\\s+")))
testGroup<-getUnixArray("testGroups")

#extract sampleNames from varSettings.sh file
sampleNames<-getUnixArray("sampleNames")


#load main genome to which you wish to align sequences and store it under the geneirc "genome" variable
suppressMessages(library("BSgenome.Celegans.UCSC.ce11"))
genome<-Celegans

#genomeVer="WS235"
genomeVer<-getUnixVar("genomeVer")

print(genomeVer)
if (grepl("^WS",genomeVer)) {
	GenomeInfoDb::seqlevelsStyle(genome)<-"Ensembl"
}

genomeFile="/home/ubelix/izb/bi18k694/genomeversion/ws270/c_elegans.PRJNA13758.WS270.genomic.fa"
#genomeFile<-"/home/ubelix/izb/semple/genomeVer/ws260/sequence/c_elegans.PRJNA13758.WS260.genomic.fa"


#load secondary genome to which you wish to align sequences and store it under the geneirc "auxGenome" variable
#library("BSgenome.Ecoli.NCBI.20080805")
#auxGenomeName<-c("Escherichia_coli_ensembl", "PhiX174_NC_001422.1", "LambdaPhage_NC_001416.1")
#auxGenomeFile<-c("~/genomeVer/ecoli/Escherichia_coli.HUSEC2011CHR1.dna.fa", 
#        "~/genomeVer/phiX/phiX.fasta", "~/genomeVer/lambda/lambdaPhage.fasta") 

print("importing genomic ranges for amplicons")
# file with genomic ranges for amplicons. Must contain a metadata column called "ID" with a unique name for
# each amplicon (e.g. gene name)
amplicons<-readRDS('/home/ubelix/izb/bi18k694/usefulfiles/ampliconGR.RDS')
if (grepl("^WS",genomeVer)) {
  GenomeInfoDb::seqlevelsStyle(amplicons)<-"Ensembl"
}
names(mcols(amplicons))[1]<-"ID"

# file with genomic ranges for TSS (or other genomic feature). Must contain a metadata column called "ID" with a unique name for
# each TSS (e.g. gene name). This ID should be the same as for the amplicons
ampTSS<-readRDS('/home/ubelix/izb/bi18k694/usefulfiles/ampliconMaxTSSgr.RDS')
if (grepl("^WS", genomeVer)) {
  GenomeInfoDb::seqlevelsStyle(ampTSS)<-"Ensembl"
}
names(mcols(ampTSS))[1]<-"ID"

print("importing genomic ranges for genome wide TSSs")
#files with TSS for gw alignments
#highConfTSS where all three datasets agree on the TSS with the maxTSS (872 genes)
highConfTSS<-readRDS('/home/ubelix/izb/bi18k694/usefulfiles/ChenKreusSaitoTSS_highConf_872.RDS')
highConfTSS$ID<-names(highConfTSS)
if (grepl("^WS",genomeVer)) {
  GenomeInfoDb::seqlevelsStyle(highConfTSS)<-"Ensembl"
} else {
  GenomeInfoDb::seqlevelsStyle(highConfTSS)<-"UCSC"
}

#lessConfTSS which are maxTSS from combining the three TSS datasets but were not identical (1955 genes)
lessConfTSS<-readRDS('/home/ubelix/izb/bi18k694/usefulfiles/ChenKreusSaitoTSS_lessConf_1955.RDS')
lessConfTSS$ID<-names(lessConfTSS)
if (grepl("^WB",genomeVer)) {
  GenomeInfoDb::seqlevelsStyle(lessConfTSS)<-"Ensembl"
} else {
  GenomeInfoDb::seqlevelsStyle(lessConfTSS)<-"UCSC"
}

print("finished importing variables")

