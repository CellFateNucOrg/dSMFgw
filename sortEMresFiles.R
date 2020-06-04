dirToSort<-"./EMres"
fileNameLst<-list.files(dirToSort,pattern="WBGene")
geneNameLst<-unique(unlist(regmatches(fileNameLst,gregexpr("WBGene[[:digit:]]{8}",fileNameLst))))

for(gene in geneNameLst) {
  dir.create(paste0(dirToSort,"/",gene))
  idx<-grep(gene,fileNameLst)
  file.copy(paste0(dirToSort,"/",fileNameLst[idx]),
  paste0(dirToSort,"/",gene,"/",fileNameLst[idx]))
  file.remove(paste0(dirToSort,"/",fileNameLst[idx]))
}

