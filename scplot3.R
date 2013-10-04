

scplot3<-function( ge, x, y, z, nbrk=50 , title=""){

tiff(file = paste( x,"x",y, "x", z ,".tiff",sep=""), width =7.3, height = 7.3, units = "in" ,res = 300,  compression = "lzw") 

cc <- colorRampPalette(c("blue","red"))(nbrk)

plot(ge[x,] , ge[y,], col=cc[cut(ge[z,],nbrk)],
	main=title,
	xlab=x,ylab=y,
	pch=16 )
mtext(z,4,cex=0.6)


dev.off()
}
