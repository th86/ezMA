#Three-Way Scatter plot
#Tai-Hsien Ou Yang

scplot3<-function( ge, x, y, z, nbrk=50 , title=""){

tiff(file = paste( x,"x",y, "x", z ,".tiff",sep=""), width =7.3, height = 7.3, units = "in" ,res = 300,  compression = "lzw") 

par(mar=c(4,4,4,4)) #oma=c( 4,2,2,3)
par( mfrow = c( 3, 2 ) )

cc <- colorRampPalette(c("blue","red"))(nbrk)

plot(ge[x,] , ge[y,], col=cc[cut(ge[z,],nbrk)],
	main=title,
	xlab=x,ylab=y,
	pch=16 )
mtext(z,4,cex=0.6)
image.plot( legend.only=TRUE, zlim= c(min(ge[z,]),max(ge[z,])), smallplot=c(0.9,0.91, .2,.8)) 


dev.off()
}
