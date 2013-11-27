#Time-split Sliding Window Survival Curve
#Tai-Hsien Ou Yang
#Nov. 26, 2013
#Columbia University

surv.rate<-function( pred ,surv, surv_thres=3650, windowRange=100 ){

surv.10yr.d<-rownames(surv)[which(surv[,1]<surv_thres & surv[,2]==1)]
surv.10yr.s<-rownames(surv)[c( which(surv[,1]<surv_thres & surv[,2]==1), 
								which(surv[,1]>=surv_thres)   ) ]

surv.gr<-list()
surv.rate<-rep(0,length(pred))

pred.all.rank<-rank(pred)
pred.all.rank<-sort(pred.all.rank)
surv.10yr.sample<-names(pred.all.rank)
			
	for( i in 1:length(pred.all.rank) ) {
		surv.gr[[i]]<-unique(c(names(pred.all.rank[max(1,i-(windowRange)):i]  ), 
			                   names(pred.all.rank[i:min( length(pred.all.rank) ,i+(windowRange) ) ] ) ))
		surv.rate[i]<-1- length(intersect( surv.10yr.d, surv.gr[[i]]))/length(intersect( surv.10yr.s, surv.gr[[i]]) )	     
	}

	survObj<-list(	rank=pred.all.rank/length(pred.all.rank), 
					value=sort(pred),
					group.list=surv.gr , 
					rate=surv.rate  )
	return(survObj)
}