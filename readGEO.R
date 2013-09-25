#Read GEO series_matrix.txt files
#Author: Tai-Hsien Ou Yang, Columbia University

library("GEOquery")
library("hgu133plus2.db")
library(impute)


#fileName="cytonormal.LAML.GSE12417-GPL570_series_matrix.txt"

readGEO<-function( fileName,saveRda=FALSE,transformLog=TRUE){

	#Read from series_matrix.txt file
	gse<-getGEO(filename= fileName)
	e<-exprs(gse)
	clinical<-phenoData(gse)
	clinical<-clinical@data
	#sample_id<-rownames(clinical@data)

	#imputation
	e<-impute.knn(e)$data
	#Annotation
	probeIDs<-ls(hgu133plus2SYMBOL)
	symbols<-mget(probeIDs,hgu133plus2SYMBOL)
	GPL570Symbols<-unlist(symbols)
	genes.shared=intersect(rownames(e),names(GPL570Symbols))
	e.reduced<-e[genes.shared,]
	rownames(e.reduced)<- GPL570Symbols[genes.shared]

	#log transformation
	if(transformLog==TRUE){
		e.reduced=log2(e.reduced+0.5)
	}

	geo.data<-list()
	geo.data$e<-e.reduced
	geo.data$clinical<-clinical

	if(saveRda==TRUE){
		save(geo.data,  file=paste(fileName,".rda",sep="") )
	}

	return(geo.data)
}