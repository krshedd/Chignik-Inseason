#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### SK ####
setwd("V:/Documents/_DOC_JunkDrawer/Power point presentations/Sockeye/Chignik/2015 CRAA")
require(xlsx)

# Create function to calculate proprtions based on Birch's logistic function coefficients (kappa, gamma method vs. a, b)
# NOTE: kappa = b, gamma = -a/b
logistic=function(x,kappa,gamma){
  p=1/(1+exp(-kappa*(x-gamma)))
  return(p)
}


par(mar=c(3.1, 4.6, 1.1, 1.1))
par(family="serif")

Birch.WLS.RegCoeffs=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="Birch.WLS.RegCoeffs",stringsAsFactors=FALSE,colClasses=c(rep("numeric",3)))
Birch.WLS.RegCoeffs$col <- colorRampPalette(c("black", "skyblue"))(7)
Birch.WLS.RegCoeffs$lty <- c(1:6, 1)

plot(-5,cex=1.5,xlim=c(1,76),ylim=c(0,1),bty="n",axes=FALSE,cex.lab=2,xlab="",ylab="Proportion Black Lake",col.lab=1)
mtext(text = "Date", side = 1, line = 2, cex = 2)
rect(xleft=29,ybottom=0,xright=69,ytop=1,col="grey90")
axis(side=1,at=seq(1,71,by=10),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3"),cex.axis=1.5,pos=0,lwd=4,col=1,col.axis=1)
axis(side=2,cex.axis=1.5,lwd=4,col=1,col.axis=1)
apply(Birch.WLS.RegCoeffs, 1, function(yr) {
  lines(x = 1:76, y = 1-logistic(x = 1:76, kappa = as.numeric(yr[2]), gamma = as.numeric(yr[3])), type = "l", lwd = 3, col = yr[4], lty = as.numeric(yr[5]))
} )
abline(h=c(0),lwd=4,col=1)
legend(x=0,y=1,bty="n",legend=2010:2016,lwd=4,col=Birch.WLS.RegCoeffs$col, lty = Birch.WLS.RegCoeffs$lty,cex=1.5,text.col=1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
