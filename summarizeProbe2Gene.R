#Summarize Top 2 Correlated Genes to a Gene 
#Author: Tai-Hsien Ou Yang, Columbia University

summarizeProbe2Gene<-function(ge){

ge_name_unique<-unique(rownames(ge))

ge_sum<-matrix(0, length(ge_name_unique), ncol(ge)) 
rownames(ge_sum)=ge_name_unique
colnames(ge_sum)=colnames(ge)

for( gene_symbol in 1:length(ge_name_unique)){
	
	gene_probes_dup<-which(rownames(ge)==ge_name_unique[gene_symbol]  ) 

	if(length(gene_probes_dup)>2){
		#cat("Found probes with duplicated gene symbol at row", gene_probes_dup,"\n")
	cor_genes<-0
	cor_genes.max<-0
	cor_genes.pairs<-NULL
	for(i in 1:length(gene_probes_dup)){
		for(j in i:length(gene_probes_dup)){
			if(j>i){
			cor_genes<-cor(ge[gene_probes_dup[i],] ,ge[gene_probes_dup[j],],) 
				if( cor_genes>cor_genes.max ){
					#cat(cor_genes,"\n") 
					cor_genes.max<-cor_genes
					cor_genes.pairs<-c(gene_probes_dup[i],gene_probes_dup[j])
				}
			}
		}
	}

			ge_sum[ge_name_unique[gene_symbol],]<-colSums(ge[cor_genes.pairs,])/2

	}else{

		if(length(gene_probes_dup)==2){
			ge_sum[ge_name_unique[gene_symbol],]<-colSums(ge[gene_probes_dup,])/2
		}else{
			ge_sum[ge_name_unique[gene_symbol],]<-ge[gene_probes_dup,]
		}

	}

}

return(ge_sum)
}