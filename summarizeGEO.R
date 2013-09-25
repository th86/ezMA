#Summarize subset files parsed by bigGEO
#Tai-Hsien Ou Yang
#Genomic Information Systems Laboratory
#Department of Electrical Engineering
#Columbia University
#2960 Broadway, NY, USA

data.sources = list.files(pattern="*.rda")

load("header.rda")
sample_id<-header$sample_id

#load( data.sources[1] ) #GPL platform data as the first
#geneSymbol<-make.unique(as.character(Table(gsmlist[[ 1 ]]@dataTable)[,11]))

load( data.sources[1] )
probeset_ID<-Table(gsmlist[[2]]@dataTable)[,"ID_REF"]
ge<-matrix(NA,length(probeset_ID), 1  )
rownames(ge) <- probeset_ID



for( file.itr in 1:(length(data.sources)-1) ){ #1 is GPL platform data
	message(sprintf("Loading %d split...", file.itr ))

	load( data.sources[file.itr] )
	for( sample.itr in 1:length(gsmlist) ){
		if( names(Table(gsmlist[[ sample.itr ]]@dataTable))[2]== "VALUE" ){
			ge.vec<-Table(gsmlist[[ sample.itr ]]@dataTable)[,"VALUE"]
			ge=cbind(ge, ge.vec)
			colnames(ge)[ncol(ge)]=names(gsmlist)[sample.itr]
		}
	}

}
ge=ge[,-1]


save(ge,file="ge.rda")