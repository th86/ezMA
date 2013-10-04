#Hypergeometric Test
#Tai-Hsien Ou Yang


fisher.test<-function(k,n,K,N){
	#N = population, K = signature set size, n = significant sample size, k = sample enriched in signature set 
	p = K/N		#uniform distribution
	# calculate p-values
	pvals = pbinom(k-1, prob=p, size=n, lower.tail=F)
	return(pvals)
}

#Affymetrix 54675
#illumina 20530
enrichment.pval<-function( gene.list, std.gene.list, probeSize=54675   ) {
		k=length(which( gene.list  %in%  std.gene.list ))
		n=length(gene.list)
		K=length(std.gene.list)
		N=probeSize
		return(fisher.test( k, n, K, N))
}
