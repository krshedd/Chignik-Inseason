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

Birch.WLS.RegCoeffs=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="Birch.WLS.RegCoeffs",stringsAsFactors=FALSE,colClasses=c(rep("numeric",3)))
Birch.WLS.RegCoeffs$col <- colorRampPalette(c("black", "skyblue"))(7)
Birch.WLS.RegCoeffs$lty <- c(1:6, 1)

png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2017/Logistic Regression Summary Plot.png", width = 7, height = 6.5, units = "in", res = 400)

par(mar=c(3.1, 4.6, 1.1, 1.1))
par(family="serif")
par(bg=rgb(red=31,green=73,blue=125,maxColorValue=255))

plot(-5,cex=1.5,xlim=c(1,76),ylim=c(0,1),bty="n",axes=FALSE,cex.lab=2,xlab="",ylab="",col.lab=1)
mtext(text = "Date", side = 1, line = 2, cex = 2, col = "white")
mtext(text = "Proportion Black Lake", side = 2, line = 3, cex = 2, col = "white")
rect(xleft=29,ybottom=0,xright=69,ytop=1,col="white")
axis(side=1,at=seq(1,71,by=10),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,cex.axis=1.5,lwd=4,col="white",col.axis="white")
apply(Birch.WLS.RegCoeffs, 1, function(yr) {
  lines(x = 1:76, y = 1-logistic(x = 1:76, kappa = as.numeric(yr[2]), gamma = as.numeric(yr[3])), type = "l", lwd = 3, col = yr[4], lty = as.numeric(yr[5]))
} )
yr2 <- as.vector(Birch.WLS.RegCoeffs[7, ])
lines(x = 1:76, y = 1-logistic(x = 1:76, kappa = as.numeric(yr2[2]), gamma = as.numeric(yr2[3])), type = "l", lwd = 10, col = as.character(yr2[4]), lty = as.numeric(yr2[5]))

abline(h=c(0),lwd=4,col="white")
legend(x=0,y=1,bty="n",legend=2010:2016,lwd=c(rep(4, 6), 10),col=Birch.WLS.RegCoeffs$col, lty = Birch.WLS.RegCoeffs$lty,cex=1.5,text.col="white")

dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
