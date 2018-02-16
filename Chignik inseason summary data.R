## 2014 Chignik Escapement mixture analysis - Inseason! - CRAA Presentation
## Kyle Shedd Fri Jan 23 14:35:04 2015

date()

ls()
rm(list=ls(all=TRUE))
search()
getwd()
# setwd("V:/DOC/Power point presentations/Sockeye/Chignik/2015 CRAA")
# setwd("V:/Documents/_DOC_JunkDrawer/Power point presentations/Sockeye/Chignik/2015 CRAA")
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
opar=par()

## save.image("V:/DOC/Power point presentations/Sockeye/Chignik/2015 CRAA/Chignik inseason summary data.RData")
##                                                                                                                                              

# Determine Julian date                                                                                                          
x=readClipboard()                   
as.Date(x,format="%m/%d/%Y")                                       
as.POSIXlt(as.Date(x,format="%m/%d/%Y"))$yday+1
writeClipboard(as.character(as.POSIXlt(as.Date(x,format="%m/%d/%Y"))$yday+1))

# Determine "Day" with 5/25 as day 1, NOTE, still need to fix leap year (2012)
as.POSIXlt(as.Date("5/25/2010",format="%m/%d/%Y"))$yday-143
as.POSIXlt(as.Date("7/31/2010",format="%m/%d/%Y"))$yday-143
as.numeric(as.POSIXlt(as.Date("8/31/2010",format="%m/%d/%Y"))-as.POSIXlt(as.Date("5/24/2010",format="%m/%d/%Y")))
x=readClipboard()
writeClipboard(as.character(as.POSIXlt(as.Date(x,format="%m/%d/%Y"))$yday-143)) # NEED TO CHANGE 2012 (SUBTRACT 1 DAY)


### Read in data
require(xlsx)

  Estimates.2010=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2010",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
  Estimates.2011=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2011",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
  Estimates.2012=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2012",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
  Estimates.2013=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2013",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
  Estimates.2014=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2014",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
  Estimates.2015=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2015",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
  Estimates.2016=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2016",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
  
# Add model assumption of 0% May 25, 100% July 31, 
# Weighting for these points was determined to have a SD of 0.025 via Birch Foster
# I found 0.0095 when taking the mean of estimates close to 1 or 0 ...
# Perhpas a better way is to take the SD from 100% proof tests! which would give
  AnchorSDs=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="AnchorSDs",stringsAsFactors=FALSE,colClasses=c(rep("numeric",2)))[,2:3]
  EarlyAssumption=c(NA,145,1,NA,1.0,1+qnorm(0.05)*AnchorSDs[1,1],1,AnchorSDs[1,1],0,0,0+qnorm(0.95)*AnchorSDs[1,1],AnchorSDs[1,1])
  LateAssumption=c(NA,212,68,NA,0,0,0+qnorm(0.95)*AnchorSDs[1,2],AnchorSDs[1,2],1.0,1+qnorm(0.05)*AnchorSDs[1,2],1,AnchorSDs[1,2])

  for(i in 0:6){
    assign(paste("Estimates.201",i,sep=""),rbind(get(paste("Estimates.201",i,sep="")),EarlyAssumption,LateAssumption))
  }

  Estimates.2010[(length(Estimates.2010$Julian)-1):length(Estimates.2010$Julian),1]=c("5/25","7/31")
  Estimates.2011[(length(Estimates.2011$Julian)-1):length(Estimates.2011$Julian),1]=c("5/25","7/31")
  Estimates.2012[(length(Estimates.2012$Julian)-1):length(Estimates.2012$Julian),1]=c("5/25","7/31")
  Estimates.2013[(length(Estimates.2013$Julian)-1):length(Estimates.2013$Julian),1]=c("5/25","7/31")
  Estimates.2014[(length(Estimates.2014$Julian)-1):length(Estimates.2014$Julian),1]=c("5/25","7/31")
  Estimates.2015[(length(Estimates.2015$Julian)-1):length(Estimates.2015$Julian),1]=c("5/25","7/31")

# Determine Inverse variance
  Estimates.2010$InvVar=round(1/(Estimates.2010$ChigSD)^2,2)
  Estimates.2011$InvVar=round(1/(Estimates.2011$ChigSD)^2,2)
  Estimates.2012$InvVar=round(1/(Estimates.2012$ChigSD)^2,2)
  Estimates.2013$InvVar=round(1/(Estimates.2013$ChigSD)^2,2)
  Estimates.2014$InvVar=round(1/(Estimates.2014$ChigSD)^2,2)
  Estimates.2015$InvVar=round(1/(Estimates.2015$ChigSD)^2,2)

# Create version of Birch input data
  Birch.Estimates.2010=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2010.Birch",stringsAsFactors=FALSE,colClasses=c("Date",rep("numeric",3)))
  Birch.Estimates.2011=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2011.Birch",stringsAsFactors=FALSE,colClasses=c("Date",rep("numeric",3)))
  Birch.Estimates.2012=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2012.Birch",stringsAsFactors=FALSE,colClasses=c("Date",rep("numeric",3)))
  Birch.Estimates.2013=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2013.Birch",stringsAsFactors=FALSE,colClasses=c("Date",rep("numeric",3)))
  Birch.Estimates.2014=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2014.Birch",stringsAsFactors=FALSE,colClasses=c("Date",rep("numeric",3)))
  Birch.Estimates.2015=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2015.Birch",stringsAsFactors=FALSE,colClasses=c("Date",rep("numeric",3)))

# Get coefficients from Birch's WLS method
  Birch.WLS.RegCoeffs=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="Birch.WLS.RegCoeffs",stringsAsFactors=FALSE,colClasses=c(rep("numeric",3)))
  Birch.WLS.RegCoeffs=as.matrix(Birch.WLS.RegCoeffs[,2:3])
  rownames(Birch.WLS.RegCoeffs)=2010:2016

# Create function to calculate proprtions based on Birch's logistic function coefficients (kappa, gamma method vs. a, b)
# NOTE: kappa = b, gamma = -a/b
logistic=function(x,kappa,gamma){
  p=1/(1+exp(-kappa*(x-gamma)))
  return(p)
}
# Test
logistic(30,Birch.WLS.RegCoeffs["2015","kappa"],Birch.WLS.RegCoeffs["2015","gamma"]) # Works

# Check out w/o anchors
  data.Estimates.2010=Estimates.2010[1:(length(Estimates.2010$Julian)-2),]
  data.Estimates.2011=Estimates.2011[1:(length(Estimates.2011$Julian)-2),]
  data.Estimates.2012=Estimates.2012[1:(length(Estimates.2012$Julian)-2),]
  data.Estimates.2013=Estimates.2013[1:(length(Estimates.2013$Julian)-2),]
  data.Estimates.2014=Estimates.2014[1:(length(Estimates.2014$Julian)-2),]
  data.Estimates.2015=Estimates.2015[1:(length(Estimates.2015$Julian)-2),]


## Compare 2014 real quick for "Birch vs. R" worksheet
Birch.Estimates.2014
Birch.fit.w.Var.2014=glm(ChigMean~Day,data=Birch.Estimates.2014,weights=InvVar,family=binomial(logit))
Birch.fit.w.Var.2014$coefficients
sumtab=summary(Birch.fit.w.Var.2014)$coefficients[1:2,1:2]
sumtab[1,1] = -sumtab[1,1]/sumtab[2,1]
sumtab
writeClipboard(as.character(rev(sumtab[,1])))
writeClipboard(as.character(rev(sumtab[,2])))

writeClipboard(as.character(Birch.fit.w.Var.2014$coefficients))
Birch.model.fit.2014=data.frame(Day=c(1:69))
Birch.model.fit.2014$Estimate.w.Var=round(predict(Birch.fit.w.Var.2014,Birch.model.fit.2014,type="response"),4)
writeClipboard(as.character(Birch.model.fit.2014$Estimate.w.Var))


# Paste into excel sheet
  for(yr in 2010:2014){  
    write.xlsx(get(paste("Estimates.",yr,sep="")),file="Chignik inseason summary data.xlsx",sheetName=paste(yr,"_Estimates",sep=""),append=TRUE,row.names=FALSE)
  }
write.xlsx(get(paste("Estimates.2015", sep = "")),file="Chignik inseason summary data.xlsx", sheetName = paste("2015_Estimates", sep = ""), append = TRUE, row.names = FALSE)


### Chignik Lake
## Fit Inverse Variance weighted Logistic Model with and without assumptions of 0% on 5/25 and 100% on 7/31
## NOTE THAT THESE GLM MODELS ARE FIT BY MINIMIZING THE ITERATIVELY WEIGHTED LEAST SQUARES (IWLS)
# 2010
    fit.w.Var.2010=glm(ChigMean~Day,data=Estimates.2010,weights=InvVar,family=binomial(logit))
    data.fit.w.Var.2010=glm(ChigMean~Day,data=data.Estimates.2010,weights=InvVar,family=binomial(logit))
    Birch.fit.w.Var.2010=glm(ChigMean~Day,data=Birch.Estimates.2010,weights=n,family=binomial(logit))
    #round(fit.w.Var.2010$fitted,3)
# 2011
    fit.w.Var.2011=glm(ChigMean~Day,data=Estimates.2011,weights=InvVar,family=binomial(logit))
    data.fit.w.Var.2011=glm(ChigMean~Day,data=data.Estimates.2011,weights=InvVar,family=binomial(logit))
    Birch.fit.w.Var.2011=glm(ChigMean~Day,data=Birch.Estimates.2011,weights=n,family=binomial(logit))
    #round(fit.w.Var.2011$fitted,3)
# 2012
    fit.w.Var.2012=glm(ChigMean~Day,data=Estimates.2012,weights=InvVar,family=binomial(logit))
    data.fit.w.Var.2012=glm(ChigMean~Day,data=data.Estimates.2012,weights=InvVar,family=binomial(logit))
    Birch.fit.w.Var.2012=glm(ChigMean~Day,data=Birch.Estimates.2012,weights=n,family=binomial(logit))
    #round(fit.w.Var.2012$fitted,3)
# 2013
    fit.w.Var.2013=glm(ChigMean~Day,data=Estimates.2013,weights=InvVar,family=binomial(logit))
    data.fit.w.Var.2013=glm(ChigMean~Day,data=data.Estimates.2013,weights=InvVar,family=binomial(logit))
    Birch.fit.w.Var.2013=glm(ChigMean~Day,data=Birch.Estimates.2013,weights=n,family=binomial(logit))
    #round(fit.w.Var.2013$fitted,3)
# 2014
    fit.w.Var.2014=glm(ChigMean~Day,data=Estimates.2014,weights=InvVar,family=binomial(logit))
    data.fit.w.Var.2014=glm(ChigMean~Day,data=data.Estimates.2014,weights=InvVar,family=binomial(logit))
    Birch.fit.w.Var.2014=glm(ChigMean~Day,data=Birch.Estimates.2014,weights=InvVar,family=binomial(logit))
    #round(fit.w.Var.2014$fitted,3)
# 2015
  fit.w.Var.2015=glm(ChigMean~Day,data=Estimates.2015,weights=InvVar,family=binomial(logit))
  data.fit.w.Var.2015=glm(ChigMean~Day,data=data.Estimates.2015,weights=InvVar,family=binomial(logit))
  Birch.fit.w.Var.2015=glm(ChigMean~Day,data=Birch.Estimates.2015,weights=InvVar,family=binomial(logit))
  #round(fit.w.Var.2015$fitted,3)

## Verify goodness of fit
# Summary shows regression coefficients w/ SE and z-test
# Chi-sq test of residual deviance with residual df tests whether model is significantly different that data, if low, then model is not a good fit
# $model gives the exact inputs used
  summary(fit.w.Var.2010)
    1-pchisq(fit.w.Var.2010$deviance,fit.w.Var.2010$df.residual)
  summary(fit.w.Var.2011)
    1-pchisq(fit.w.Var.2011$deviance,fit.w.Var.2011$df.residual)
  summary(fit.w.Var.2012)
    1-pchisq(fit.w.Var.2012$deviance,fit.w.Var.2012$df.residual)
  summary(fit.w.Var.2013)
    1-pchisq(fit.w.Var.2013$deviance,fit.w.Var.2013$df.residual)
  summary(fit.w.Var.2014)
    1-pchisq(fit.w.Var.2014$deviance,fit.w.Var.2014$df.residual)

# Extract Null Deviance
comp.models=matrix(data=c(fit.w.Var.2010$deviance,data.fit.w.Var.2010$deviance,Birch.fit.w.Var.2010$deviance,
                          fit.w.Var.2011$deviance,data.fit.w.Var.2011$deviance,Birch.fit.w.Var.2011$deviance,
                          fit.w.Var.2012$deviance,data.fit.w.Var.2012$deviance,Birch.fit.w.Var.2012$deviance,
                          fit.w.Var.2013$deviance,data.fit.w.Var.2013$deviance,Birch.fit.w.Var.2013$deviance,
                          fit.w.Var.2014$deviance,data.fit.w.Var.2014$deviance,Birch.fit.w.Var.2014$deviance),
                   nrow=5,ncol=3,byrow=TRUE,dimnames=list(c(2010:2014),c("Clean","NoAnchor","Birch")))
# Extract Residual DFs
comp.models.df=matrix(data=c(fit.w.Var.2010$df.residual,data.fit.w.Var.2010$df.residual,Birch.fit.w.Var.2010$df.residual,
                             fit.w.Var.2011$df.residual,data.fit.w.Var.2011$df.residual,Birch.fit.w.Var.2011$df.residual,
                             fit.w.Var.2012$df.residual,data.fit.w.Var.2012$df.residual,Birch.fit.w.Var.2012$df.residual,
                             fit.w.Var.2013$df.residual,data.fit.w.Var.2013$df.residual,Birch.fit.w.Var.2013$df.residual,
                             fit.w.Var.2014$df.residual,data.fit.w.Var.2014$df.residual,Birch.fit.w.Var.2014$df.residual),
                      nrow=5,ncol=3,byrow=TRUE,dimnames=list(c(2010:2014),c("Clean","NoAnchor","Birch")))
# Extract AIC
comp.models.aic=matrix(data=c(fit.w.Var.2010$aic,data.fit.w.Var.2010$aic,Birch.fit.w.Var.2010$aic,
                              fit.w.Var.2011$aic,data.fit.w.Var.2011$aic,Birch.fit.w.Var.2011$aic,
                              fit.w.Var.2012$aic,data.fit.w.Var.2012$aic,Birch.fit.w.Var.2012$aic,
                              fit.w.Var.2013$aic,data.fit.w.Var.2013$aic,Birch.fit.w.Var.2013$aic,
                              fit.w.Var.2014$aic,data.fit.w.Var.2014$aic,Birch.fit.w.Var.2014$aic),
                        nrow=5,ncol=3,byrow=TRUE,dimnames=list(c(2010:2014),c("Clean","NoAnchor","Birch")))
  comp.models
  comp.models.df
  comp.models.aic

# Determin if model is significantly different than data
  1-pchisq(comp.models,comp.models.df)

##  WOW, it looks like weighint the models is BAAAAD


# Anova shows the reduction in deviance for each model term added sequentially from first to last w/ Chi-sq test
  anova(fit.w.Var.2010, test = "Chisq")  
  anova(fit.w.Var.2011, test = "Chisq")  
  anova(fit.w.Var.2012, test = "Chisq")  
  anova(fit.w.Var.2013, test = "Chisq")  
  anova(fit.w.Var.2014, test = "Chisq")  


  par(mfrow=c(2,2))
  for(i in 2010:2014){
    plot(get(paste("fit.w.Var.",i,sep="")))
  }
  par(mfrow=c(1,1))

  par(mfrow=c(2,2))
  for(i in 2010:2014){
    plot(get(paste("data.fit.w.Var.",i,sep="")))
  }
  par(mfrow=c(1,1))

  par(mfrow=c(2,2))
  for(i in 2010:2014){
    plot(get(paste("Birch.fit.w.Var.",i,sep="")))
  }
  par(mfrow=c(1,1))

## Predict daily values - August 31 is 243, September 30 is 273
# 2010
    model.fit.2010=data.frame(Day=c(1:99))
    model.fit.2010$Date=as.Date(model.fit.2010$Day,origin=as.Date("2010-05-24"))
    model.fit.2010$Estimate.w.Var=round(predict(fit.w.Var.2010,model.fit.2010,type="response"),4)
    model.fit.2010$data.Estimate.w.Var=round(predict(data.fit.w.Var.2010,model.fit.2010,type="response"),4)
    model.fit.2010$BirchIWLS.Estimate=round(predict(Birch.fit.w.Var.2010,model.fit.2010,type="response"),4)
    model.fit.2010$BirchWLS.Estimate=round(logistic(x=model.fit.2010$Day,kappa=Birch.WLS.RegCoeffs["2010","kappa"],gamma=Birch.WLS.RegCoeffs["2010","gamma"]),4)
#    preds=predict(fit.w.Var.2010, newdata = data.frame(Day=model.fit.2010$Day), type = "response", se.fit = TRUE)
#   write.xlsx(x=model.fit.2010,file="Chignik inseason summary data.xlsx",sheetName="2010_Logistic",append=TRUE,row.names=FALSE)
# 2011
    model.fit.2011=data.frame(Day=c(1:99))
    model.fit.2011$Date=as.Date(model.fit.2011$Day,origin=as.Date("2011-05-24"))
    model.fit.2011$Estimate.w.Var=round(predict(fit.w.Var.2011,model.fit.2011,type="response"),4)
    model.fit.2011$data.Estimate.w.Var=round(predict(data.fit.w.Var.2011,model.fit.2011,type="response"),4)
    model.fit.2011$BirchIWLS.Estimate=round(predict(Birch.fit.w.Var.2011,model.fit.2011,type="response"),4)
    model.fit.2011$BirchWLS.Estimate=round(logistic(x=model.fit.2011$Day,kappa=Birch.WLS.RegCoeffs["2011","kappa"],gamma=Birch.WLS.RegCoeffs["2011","gamma"]),4)
#   write.xlsx(x=model.fit.2011,file="Chignik inseason summary data.xlsx",sheetName="2011_Logistic",append=TRUE,row.names=FALSE)
# 2012, modified Day date due to leap year so that July 4 is 185 across all years
    model.fit.2012=data.frame(Day=c(1:99))
    model.fit.2012$Date=as.Date(model.fit.2012$Day,origin=as.Date("2012-05-24"))
    model.fit.2012$Estimate.w.Var=round(predict(fit.w.Var.2012,model.fit.2012,type="response"),4)
    model.fit.2012$data.Estimate.w.Var=round(predict(data.fit.w.Var.2012,model.fit.2012,type="response"),4)
    model.fit.2012$BirchIWLS.Estimate=round(predict(Birch.fit.w.Var.2012,model.fit.2012,type="response"),4)
    model.fit.2012$BirchWLS.Estimate=round(logistic(x=model.fit.2012$Day,kappa=Birch.WLS.RegCoeffs["2012","kappa"],gamma=Birch.WLS.RegCoeffs["2012","gamma"]),4)
#   write.xlsx(x=model.fit.2012,file="Chignik inseason summary data.xlsx",sheetName="2012_Logistic",append=TRUE,row.names=FALSE)
# 2013
    model.fit.2013=data.frame(Day=c(1:99))
    model.fit.2013$Date=as.Date(model.fit.2013$Day,origin=as.Date("2013-05-24"))
    model.fit.2013$Estimate.w.Var=round(predict(fit.w.Var.2013,model.fit.2013,type="response"),4)
    model.fit.2013$data.Estimate.w.Var=round(predict(data.fit.w.Var.2013,model.fit.2013,type="response"),4)
    model.fit.2013$BirchIWLS.Estimate=round(predict(Birch.fit.w.Var.2013,model.fit.2013,type="response"),4)
    model.fit.2013$BirchWLS.Estimate=round(logistic(x=model.fit.2013$Day,kappa=Birch.WLS.RegCoeffs["2013","kappa"],gamma=Birch.WLS.RegCoeffs["2013","gamma"]),4)
#   write.xlsx(x=model.fit.2013,file="Chignik inseason summary data.xlsx",sheetName="2013_Logistic",append=TRUE,row.names=FALSE)
# 2014
    model.fit.2014=data.frame(Day=c(1:99))
    model.fit.2014$Date=as.Date(model.fit.2014$Day,origin=as.Date("2014-05-24"))
    model.fit.2014$Estimate.w.Var=round(predict(fit.w.Var.2014,model.fit.2014,type="response"),4)
    model.fit.2014$data.Estimate.w.Var=round(predict(data.fit.w.Var.2014,model.fit.2014,type="response"),4)
    model.fit.2014$BirchIWLS.Estimate=round(predict(Birch.fit.w.Var.2014,model.fit.2014,type="response"),4)
    model.fit.2014$BirchWLS.Estimate=round(logistic(x=model.fit.2014$Day,kappa=Birch.WLS.RegCoeffs["2014","kappa"],gamma=Birch.WLS.RegCoeffs["2014","gamma"]),4)
#   write.xlsx(x=model.fit.2014,file="Chignik inseason summary data.xlsx",sheetName="2014_Logistic",append=TRUE,row.names=FALSE)
# 2015
    model.fit.2015=data.frame(Day=c(1:99))
    model.fit.2015$Date=as.Date(model.fit.2015$Day,origin=as.Date("2015-05-24"))
    model.fit.2015$Estimate.w.Var=round(predict(fit.w.Var.2015,model.fit.2015,type="response"),4)
    model.fit.2015$data.Estimate.w.Var=round(predict(data.fit.w.Var.2015,model.fit.2015,type="response"),4)
    model.fit.2015$BirchIWLS.Estimate=round(predict(Birch.fit.w.Var.2015,model.fit.2015,type="response"),4)
    model.fit.2015$BirchWLS.Estimate=round(logistic(x=model.fit.2015$Day,kappa=Birch.WLS.RegCoeffs["2015","kappa"],gamma=Birch.WLS.RegCoeffs["2015","gamma"]),4)
write.xlsx(x=model.fit.2015,file="Chignik inseason summary data.xlsx",sheetName="2015_Logistic",append=TRUE,row.names=FALSE)

### Plot data!
dir.create(path="Plots")
setwd("Plots")

psh=c(15:20)
clr=c("white",2:6)
par(mar=c(5.1,5.1,2.1,2.1))
wid=4
par(bg=rgb(red=31,green=73,blue=125,maxColorValue=255))

# 2010
    plot(BlackMean~Day,data=Estimates.2010[1:11,],pch=psh[1],col="white",cex=1.5,xlim=c(1,76),ylim=c(0,1),bty="n",axes=FALSE,cex.lab=2,xlab="Date",ylab="Proportion Black Lake",col.lab="white")
    axis(side=1,at=seq(1,71,by=10),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
    axis(side=2,cex.axis=1.5,lwd=4,col="white",col.axis="white")
    abline(h=c(0),lwd=4,col="white")
    arrows(x0 = Estimates.2010$Day[1:11], y0 = Estimates.2010$BlackMean[1:11], x1 = Estimates.2010$Day[1:11], y1 = Estimates.2010$Black5CI[1:11], angle = 90, length = 0.08, lwd=wid, col="white")
    arrows(x0 = Estimates.2010$Day[1:11], y0 = Estimates.2010$BlackMean[1:11], x1 = Estimates.2010$Day[1:11], y1 = Estimates.2010$Black95CI[1:11], angle = 90, length = 0.08, lwd=wid, col="white")
    lines(model.fit.2010$Day,(1-model.fit.2010$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col="white")
    #arrows(x0=-1.5, y0=0.5, x1=48.6, y1=0.5, col="darkorange2", lwd=4, angle=45)
    #arrows(x0=48.6, y0=0.5, x1=48.6, y1=0.01, col="darkorange2", lwd=4, angle=45)
#  polygon(c(rev(model.fit.2010$Day), model.fit.2010$Day), c((rev(preds[[1]])+1.96*rev(preds[[2]])),(preds[[1]]-1.96*preds[[2]])),col=8,border=NA) # polygon of 95%CI for regression

# 2011
    points(BlackMean~Day,data=Estimates.2011,pch=psh[2],col=clr[2],cex=1.5)
    arrows(x0 = Estimates.2011$Day, y0 = Estimates.2011$BlackMean, x1 = Estimates.2011$Day, y1 = Estimates.2011$Black5CI, angle = 90, length = 0.08, lwd=wid, col=clr[2])
    arrows(x0 = Estimates.2011$Day, y0 = Estimates.2011$BlackMean, x1 = Estimates.2011$Day, y1 = Estimates.2011$Black95CI, angle = 90, length = 0.08, lwd=wid, col=clr[2])
    lines(model.fit.2011$Day,(1-model.fit.2011$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[2])
# 2012
    points(BlackMean~Day,data=Estimates.2012,pch=psh[3],col=clr[3],cex=1.5)
    arrows(x0 = Estimates.2012$Day, y0 = Estimates.2012$BlackMean, x1 = Estimates.2012$Day, y1 = Estimates.2012$Black5CI, angle = 90, length = 0.08, lwd=wid, col=clr[3])
    arrows(x0 = Estimates.2012$Day, y0 = Estimates.2012$BlackMean, x1 = Estimates.2012$Day, y1 = Estimates.2012$Black95CI, angle = 90, length = 0.08, lwd=wid, col=clr[3])
    lines(model.fit.2012$Day,(1-model.fit.2012$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[3])
# 2013
    points(BlackMean~Day,data=Estimates.2013,pch=psh[4],col=clr[4],cex=1.5)
    arrows(x0 = Estimates.2013$Day, y0 = Estimates.2013$BlackMean, x1 = Estimates.2013$Day, y1 = Estimates.2013$Black5CI, angle = 90, length = 0.08, lwd=wid, col=clr[4])
    arrows(x0 = Estimates.2013$Day, y0 = Estimates.2013$BlackMean, x1 = Estimates.2013$Day, y1 = Estimates.2013$Black95CI, angle = 90, length = 0.08, lwd=wid, col=clr[4])
    lines(model.fit.2013$Day,(1-model.fit.2013$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[4])
# 2014
    points(BlackMean~Day,data=Estimates.2014,pch=psh[5],col=clr[5],cex=1.5)
    arrows(x0 = Estimates.2014$Day, y0 = Estimates.2014$BlackMean, x1 = Estimates.2014$Day, y1 = Estimates.2014$Black5CI, angle = 90, length = 0.08, lwd=wid, col=clr[5])
    arrows(x0 = Estimates.2014$Day, y0 = Estimates.2014$BlackMean, x1 = Estimates.2014$Day, y1 = Estimates.2014$Black95CI, angle = 90, length = 0.08, lwd=wid, col=clr[5])
    lines(model.fit.2014$Day,(1-model.fit.2014$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[5])
# 2015
    points(BlackMean~Day,data=Estimates.2015,pch=psh[6],col=clr[6],cex=1.5)
    arrows(x0 = Estimates.2015$Day, y0 = Estimates.2015$BlackMean, x1 = Estimates.2015$Day, y1 = Estimates.2015$Black5CI, angle = 90, length = 0.08, lwd=wid, col=clr[6])
    arrows(x0 = Estimates.2015$Day, y0 = Estimates.2015$BlackMean, x1 = Estimates.2015$Day, y1 = Estimates.2015$Black95CI, angle = 90, length = 0.08, lwd=wid, col=clr[6])
    lines(model.fit.2015$Day,(1-model.fit.2015$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[6])
  abline(h=c(0),lwd=4,col="white")

legend(x=50,y=1,bty="n",legend=2010:2015,lwd=4,col=clr,cex=2)


### View individually
for(yr in 2010:2015){
#yr=2014
  plot(ChigMean~Day,data=get(paste("Estimates.",yr,sep="")),pch=psh[yr-2009],col=clr[yr-2009],cex=1.5,xlim=c(1,76),bty="n",axes=FALSE,cex.lab=2,xlab="Date",ylab="Proportion Chignik Lake",main=yr,cex.main=3)
  axis(side=1,at=seq(145,215,by=10),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3"),cex.axis=1.5,pos=0,lwd=4)
  axis(side=2,cex.axis=1.5,lwd=4)
  
  arrows(x0 = get(paste("Estimates.",yr,sep=""))$Day, y0 = get(paste("Estimates.",yr,sep=""))$ChigMean, x1 = get(paste("Estimates.",yr,sep=""))$Day, y1 = get(paste("Estimates.",yr,sep=""))$Chig5CI, angle = 90, length = 0.08, lwd=wid, col=clr[yr-2009])
  arrows(x0 = get(paste("Estimates.",yr,sep=""))$Day, y0 = get(paste("Estimates.",yr,sep=""))$ChigMean, x1 = get(paste("Estimates.",yr,sep=""))$Day, y1 = get(paste("Estimates.",yr,sep=""))$Chig95CI, angle = 90, length = 0.08, lwd=wid, col=clr[yr-2009])
  
  lines(get(paste("model.fit.",yr,sep=""))$Day,get(paste("model.fit.",yr,sep=""))$Estimate.w.Var,type="l",lwd=wid,lty=1,col=clr[yr-2009])
  lines(get(paste("model.fit.",yr,sep=""))$Day,get(paste("model.fit.",yr,sep=""))$data.Estimate.w.Var,type="l",lwd=wid,lty=2,col=clr[yr-2009])
  lines(get(paste("model.fit.",yr,sep=""))$Day,get(paste("model.fit.",yr,sep=""))$BirchIWLS.Estimate,type="l",lwd=wid,lty=3,col=clr[yr-2009])
  lines(get(paste("model.fit.",yr,sep=""))$Day,get(paste("model.fit.",yr,sep=""))$BirchWLS.Estimate,type="l",lwd=wid,lty=4,col=clr[yr-2009])
  
  abline(h=c(0,1),lwd=4)
  legend(x=-5,y=1,bty="n",legend=c("IWLS InvVar w/ Anch","IWLS InvVar w/o Anch","IWLS n","WLS n"),lwd=4,col=clr[yr-2009],lty=1:4,cex=2)
  
  segments(x0 = trans.anchor[yr-2009], y0 = 0, x1 = trans.anchor[yr-2009], y1 = get(paste("model.fit.",yr,sep=""))$Estimate.w.Var[match(trans.anchor[yr-2009],get(paste("model.fit.",yr,sep=""))$Day)], lwd=4)
  
  par(mfrow=c(2,2))
    plot(get(paste("fit.w.Var.",yr,sep="")))
    plot(get(paste("data.fit.w.Var.",yr,sep="")))
  par(mfrow=c(1,1))
}
par(mfrow=c(1,1))


######### Explore residuals for 2014

# IWLS, ordered by Day
resid(Birch.fit.w.Var.2014,type="resp")[order(Birch.fit.w.Var.2014$data$Day)]
(Birch.fit.w.Var.2014$data$ChigMean-Birch.fit.w.Var.2014$fitted.values)[order(Birch.fit.w.Var.2014$data$Day)]

# WLS, ordered by Day
(logistic(x=Birch.fit.w.Var.2014$data$Day,kappa=Birch.WLS.RegCoeffs["2014","kappa"],gamma=Birch.WLS.RegCoeffs["2014","gamma"])-Birch.fit.w.Var.2014$data$ChigMean)[order(Birch.fit.w.Var.2014$data$Day)]

# Plot residuals
plot(resid(Birch.fit.w.Var.2014,type="resp")[order(Birch.fit.w.Var.2014$data$Day)]~Birch.fit.w.Var.2014$data$Day[order(Birch.fit.w.Var.2014$data$Day)],
     pch=16,cex=2,bty="n",axes=TRUE,xlab="Date",ylab="Residuals",type="b",lwd=6)
#segments(x0=Birch.fit.w.Var.2014$data$Day,x1=Birch.fit.w.Var.2014$data$Day,y0=0,y1=Birch.fit.w.Var.2014$fitted.values-Birch.fit.w.Var.2014$data$ChigMean,lwd=4)
points((Birch.fit.w.Var.2014$data$ChigMean-logistic(x=Birch.fit.w.Var.2014$data$Day,kappa=Birch.WLS.RegCoeffs["2014","kappa"],gamma=Birch.WLS.RegCoeffs["2014","gamma"]))[order(Birch.fit.w.Var.2014$data$Day)]~
         Birch.fit.w.Var.2014$data$Day[order(Birch.fit.w.Var.2014$data$Day)],pch=15,col=4,cex=2,type="b",lwd=6)
#segments(x0=Birch.fit.w.Var.2014$data$Day,x1=Birch.fit.w.Var.2014$data$Day,
 #        y0=0,y1=logistic(x=Birch.fit.w.Var.2014$data$Day,kappa=Birch.WLS.RegCoeffs["2014","kappa"],gamma=Birch.WLS.RegCoeffs["2014","gamma"])-Birch.fit.w.Var.2014$data$ChigMean,
  #       lwd=4,col=4)
abline(h=0,lwd=6,lty=2,col=2)

#IWLS
sum(abs(resid(Birch.fit.w.Var.2014,type="resp")[order(Birch.fit.w.Var.2014$data$Day)]))
sum((resid(Birch.fit.w.Var.2014,type="resp")[order(Birch.fit.w.Var.2014$data$Day)])^2)
#WLS
sum(abs((Birch.fit.w.Var.2014$data$ChigMean-logistic(x=Birch.fit.w.Var.2014$data$Day,kappa=Birch.WLS.RegCoeffs["2014","kappa"],gamma=Birch.WLS.RegCoeffs["2014","gamma"]))[order(Birch.fit.w.Var.2014$data$Day)]))
sum(((Birch.fit.w.Var.2014$data$ChigMean-logistic(x=Birch.fit.w.Var.2014$data$Day,kappa=Birch.WLS.RegCoeffs["2014","kappa"],gamma=Birch.WLS.RegCoeffs["2014","gamma"]))[order(Birch.fit.w.Var.2014$data$Day)])^2)



######## Deeper look at residuals

for(yr in 2010:2014){
  plot(resid(get(paste("fit.w.Var.",yr,sep="")),type="resp")[order(get(paste("fit.w.Var.",yr,sep=""))$data$Day)]~get(paste("fit.w.Var.",yr,sep=""))$data$Day[order(get(paste("fit.w.Var.",yr,sep=""))$data$Day)],pch=16,col=1,cex=3,type="b",lwd=4,
       xlab="Day",ylab="Residuals",cex.axis=2,cex.lab=1.8,main=yr,cex.main=2,ylim=c(-.15,.2))
  points(resid(get(paste("data.fit.w.Var.",yr,sep="")),type="resp")[order(get(paste("data.fit.w.Var.",yr,sep=""))$data$Day)]~get(paste("data.fit.w.Var.",yr,sep=""))$data$Day[order(get(paste("data.fit.w.Var.",yr,sep=""))$data$Day)],pch=15,col=2,cex=3,type="b",lwd=4)
  points(resid(get(paste("Birch.fit.w.Var.",yr,sep="")),type="resp")[order(get(paste("Birch.fit.w.Var.",yr,sep=""))$data$Day)]~get(paste("Birch.fit.w.Var.",yr,sep=""))$data$Day[order(get(paste("Birch.fit.w.Var.",yr,sep=""))$data$Day)],pch=17,col=3,cex=3,type="b",lwd=4)
  points((get(paste("Birch.fit.w.Var.",yr,sep=""))$data$ChigMean-logistic(x=get(paste("Birch.fit.w.Var.",yr,sep=""))$data$Day,kappa=Birch.WLS.RegCoeffs[paste(yr),"kappa"],gamma=Birch.WLS.RegCoeffs[paste(yr),"gamma"]))[order(get(paste("Birch.fit.w.Var.",yr,sep=""))$data$Day)]~
           get(paste("Birch.fit.w.Var.",yr,sep=""))$data$Day[order(get(paste("Birch.fit.w.Var.",yr,sep=""))$data$Day)],pch=18,col=4,cex=3,type="b",lwd=4)
  abline(h=0,lwd=6)
  legend("topleft",bty="n",legend=c("IWLS InvVar w/ Anch","IWLS InvVar w/o Anch","IWLS n","WLS n"),pch=c(16,15,17,18),col=1:4,cex=1.5)
}









### View just lines
clr=c(4,2,3,1,5, 6)
plot(NA,xlim=c(1,76),ylim=c(0,1),bty="n",axes=FALSE,cex.lab=2,xlab="Date",ylab="Proportion Chignik Lake",main="",cex.main=3)
axis(side=1,at=seq(1,71,by=10),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3"),cex.axis=1.5,pos=0,lwd=4)
axis(side=2,cex.axis=1.5,lwd=4)
abline(h=c(0,1),lwd=4)
for(yr in 2010:2015){
  lines(get(paste("model.fit.",yr,sep=""))$Day,get(paste("model.fit.",yr,sep=""))$Estimate.w.Var,type="l",lwd=wid,lty=1,col=clr[yr-2009])
# lines(get(paste("model.fit.",yr,sep=""))$Day,get(paste("model.fit.",yr,sep=""))$data.Estimate.w.Var,type="l",lwd=wid,lty=2,col=clr[yr-2009])
}
legend(x=1,y=1,bty="n",legend=2010:2015,lwd=4,col=clr,cex=2)

### Determine 50% transition date
  trans=0.5
  # W/ anchors
    transition=matrix(data=NA,ncol=length(2010:2015),nrow=4,dimnames=list(c("IWLS Anchor InvVar","IWLS NoAnchor InvVar","IWLS Birch n","WLS Birch n"),2010:2015))
    for(yr in 2010:2015){
      for(method in 3:6){
        transition[(method-2),(yr-2009)]=get(paste("model.fit.",yr,sep=""))$Day[which.min(abs(get(paste("model.fit.",yr,sep=""))[,method]-trans))]
      }
    }
    
    writeClipboard(as.character(as.Date(transition[,"2010"],origin="2010-05-24")))
    writeClipboard(as.character(as.Date(transition[,"2011"],origin="2011-05-24")))
    writeClipboard(as.character(as.Date(transition[,"2012"],origin="2012-05-24")))
    writeClipboard(as.character(as.Date(transition[,"2013"],origin="2013-05-24")))
    writeClipboard(as.character(as.Date(transition[,"2014"],origin="2014-05-24")))
    writeClipboard(as.character(as.Date(transition[,"2015"],origin="2015-05-24")))


    as.Date(transition[1,],origin=c("2010-05-24","2011-05-24","2012-05-24","2013-05-24","2014-05-24", "2015-05-24"))
    as.Date(transition[2,],origin=c("2010-05-24","2011-05-24","2012-05-24","2013-05-24","2014-05-24", "2015-05-24"))
    as.Date(transition[3,],origin=c("2010-05-24","2011-05-24","2012-05-24","2013-05-24","2014-05-24", "2015-05-24"))
    as.Date(transition[4,],origin=c("2010-05-24","2011-05-24","2012-05-24","2013-05-24","2014-05-24", "2015-05-24"))

    apply(transition,1,function(row){as.Date(row,origin=c("2010-05-24","2011-05-24","2012-05-24","2013-05-24","2014-05-24", "2015-05-24"))})








## Some years the residuals look really bad (especially with artificial anchors), maybe logarithmic distribution isn't perfect?
## Try polynomial

poly.fit.w.Var.2014=lm(ChigMean~poly(Day,3,raw=TRUE),data=Estimates.2014,weights=InvVar)
summary(poly.fit.w.Var.2014)
plot(poly.fit.w.Var.2014)

poly.model.fit.2014=data.frame(Day=c(160:220))
poly.model.fit.2014$Estimate.w.Var=round(predict(poly.fit.w.Var.2014,poly.model.fit.2014,type="response"),3)
lines(get(paste("poly.model.fit.",yr,sep=""))$Day,get(paste("poly.model.fit.",yr,sep=""))$Estimate.w.Var,type="l",lwd=wid,lty=1,col=clr[1])

## Very bad





## View daily escapement
# Daily escapement numbers for 5/28/2010 to 8/31/2010
BlackDailyEscp.2010 = as.numeric(readClipboard())
ChignikDailyEscp.2010 = as.numeric(readClipboard())

plot(y=BlackDailyEscp.2010, x=4:99, type="l",col="skyblue",cex=1.5,bty="n",axes=FALSE,cex.lab=2,xlab="Date",ylab="Daily Escapement",col.lab="white",lwd=4)
axis(side=1,at=c(seq(1,91,by=10),99),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,at=seq(0,50000,by=10000),cex.axis=1.5,lwd=4,col="white",col.axis="white")
points(y=ChignikDailyEscp.2010, x=4:99, type="l",col="red",lwd=4,add=TRUE)
abline(h=c(0),lwd=4,col="white")
legend("topright", bty="n", legend=c("Black Lake", "Chignik Lake"), lwd=4, col=c("skyblue", "red"), text.col="white", cex=2)



## View all of Birch's regressions
dir.create("New for AFS")
clr <- c("white", 2:6)
plot(-5,pch=psh[1],col="white",cex=1.5,xlim=c(1,76),ylim=c(0,1),bty="n",axes=FALSE,cex.lab=2,xlab="Date",ylab="Proportion Black Lake",col.lab="white")
rect(xleft=29,ybottom=0,xright=69,ytop=1,col=1)
axis(side=1,at=seq(1,71,by=10),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,cex.axis=1.5,lwd=4,col="white",col.axis="white")
abline(h=c(0),lwd=4,col="white")
lines(model.fit.2010$Day,(1-model.fit.2010$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col="white")
lines(model.fit.2011$Day,(1-model.fit.2011$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[2])
lines(model.fit.2012$Day,(1-model.fit.2012$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[3])
lines(model.fit.2013$Day,(1-model.fit.2013$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[4])
lines(model.fit.2014$Day,(1-model.fit.2014$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[5])
lines(model.fit.2015$Day,(1-model.fit.2015$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[6])
abline(h=c(0),lwd=4,col="white")
legend(x=55,y=1,bty="n",legend=2010:2015,lwd=4,col=clr,cex=1.8,text.col="white")
arrows(x0=-1,y0=0.5,x1=55,y1=0.5,lwd=wid,col=rgb(red=238,green=118,blue=0,maxColorValue=255),angle=45)  # 55
arrows(x0=41,y0=0.5,x1=55,y1=0.5,lwd=wid,col=rgb(red=238,green=118,blue=0,maxColorValue=255),angle=45,code=3)
arrows(x0=55,y0=0.5,x1=55,y1=0.01,lwd=wid,col=rgb(red=238,green=118,blue=0,maxColorValue=255),angle=45)
arrows(x0=41,y0=0.5,x1=41,y1=0.01,lwd=wid,col=rgb(red=238,green=118,blue=0,maxColorValue=255),angle=45)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### SK ####
setwd("V:/Documents/_DOC_JunkDrawer/Power point presentations/Sockeye/Chignik/2015 CRAA")

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



lines(model.fit.2010$Day,(1-model.fit.2010$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=1)
lines(model.fit.2011$Day,(1-model.fit.2011$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[2])
lines(model.fit.2012$Day,(1-model.fit.2012$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[3])
lines(model.fit.2013$Day,(1-model.fit.2013$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[4])
lines(model.fit.2014$Day,(1-model.fit.2014$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[5])
lines(model.fit.2015$Day,(1-model.fit.2015$BirchWLS.Estimate),type="l",lwd=wid,lty=1,col=clr[6])
abline(h=c(0),lwd=4,col=1)
legend(x=55,y=1,bty="n",legend=2010:2015,lwd=4,col=clr,cex=1.8,text.col=1)
arrows(x0=-1,y0=0.5,x1=55,y1=0.5,lwd=wid,col=rgb(red=238,green=118,blue=0,maxColorValue=255),angle=45)  # 55
arrows(x0=41,y0=0.5,x1=55,y1=0.5,lwd=wid,col=rgb(red=238,green=118,blue=0,maxColorValue=255),angle=45,code=3)
arrows(x0=55,y0=0.5,x1=55,y1=0.01,lwd=wid,col=rgb(red=238,green=118,blue=0,maxColorValue=255),angle=45)
arrows(x0=41,y0=0.5,x1=41,y1=0.01,lwd=wid,col=rgb(red=238,green=118,blue=0,maxColorValue=255),angle=45)




## Show that escapement magnitude matters!!!
# Show how in 2012, July 4th works well!

# Daily Escapement w/ transition period
DailyEscp.2012=as.numeric(readClipboard())

plot(y=DailyEscp.2012, x=4:sum(length(DailyEscp.2012),3), type="h", col="white", cex=1.5, xlim=c(4,sum(length(DailyEscp.2012),3)), ylim=c(0,26000), bty="n", axes=FALSE, cex.lab=2, xlab="Date", ylab="Daily Escapement", col.lab="white", lwd=6)
rect(xleft=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.05))],
     ybottom=0,xright=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.95))],ytop=26000,col=1)
points(y=DailyEscp.2012, x=4:sum(length(DailyEscp.2012),3), type="h",col="white",cex=1.5, lwd=6)
axis(side=1,at=c(seq(1,91,by=10),99),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,at=seq(0,50000,by=5000),cex.axis=1.45,lwd=4,col="white",col.axis="white")
segments(x0=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.05))],
         x1=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.05))],y0=0,y1=25900, col=1, lwd=4)
segments(x0=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.95))],
         x1=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.95))],y0=0,y1=25900, col=1, lwd=4)
abline(h=c(0),lwd=4,col="white")


# July 4th Escapement
plot(y=DailyEscp.2012[1:37], x=4:40, type="h", col="skyblue", cex=1.5, xlim=c(4,sum(length(DailyEscp.2012),3)), ylim=c(0,26000), bty="n", axes=FALSE, cex.lab=2, xlab="Date", ylab="Daily Escapement", col.lab="white", lwd=6)
rect(xleft=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.05))],
     ybottom=0,xright=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.95))],ytop=26000,col=1)
points(y=DailyEscp.2012[1:37], x=4:40, type="h", col="skyblue", lwd=6)
points(y=DailyEscp.2012[38:length(DailyEscp.2012)], x=41:sum(length(DailyEscp.2012),3), type="h", col="red", lwd=6)
axis(side=1,at=c(seq(1,91,by=10),99),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,at=seq(0,50000,by=5000),cex.axis=1.45,lwd=4,col="white",col.axis="white")
segments(x0=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.05))],
         x1=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.05))],y0=0,y1=25900, col=1, lwd=4)
segments(x0=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.95))],
         x1=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.95))],y0=0,y1=25900, col=1, lwd=4)
segments(x0=41,x1=41,y0=0,y1=25900, lwd=6, lty=1, col=rgb(red=31,green=73,blue=125,maxColorValue=255))
abline(h=c(0),lwd=4,col="white")
legend("topright", bty="n", legend=c("Black Lake", "Chignik Lake"), lwd=4, col=c("skyblue", "red"), text.col="white", cex=1.7)


# Genetics Escapement
plot(y=DailyEscp.2012, x=4:sum(length(DailyEscp.2012),3), type="h", col="skyblue", cex=1.5, xlim=c(4,sum(length(DailyEscp.2012),3)), ylim=c(0,26000), bty="n", axes=FALSE, cex.lab=2, xlab="Date", ylab="Daily Escapement", col.lab="white", lwd=6)
rect(xleft=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.05))],
     ybottom=0,xright=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.95))],ytop=26000,col=1)
points(y=DailyEscp.2012, x=4:sum(length(DailyEscp.2012),3), type="h",col="skyblue",cex=1.5, lwd=6)
points(y=(DailyEscp.2012*model.fit.2012$BirchWLS.Estimate[4:sum(length(DailyEscp.2012),3)]), x=4:sum(length(DailyEscp.2012),3), type="h",col="red",cex=1.5, lwd=6)
axis(side=1,at=c(seq(1,91,by=10),99),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,at=seq(0,50000,by=5000),cex.axis=1.45,lwd=4,col="white",col.axis="white")
segments(x0=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.05))],
         x1=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.05))],y0=0,y1=25900, col=1, lwd=4)
segments(x0=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.95))],
         x1=model.fit.2012$Day[which.min(abs(model.fit.2012$BirchWLS.Estimate-0.95))],y0=0,y1=25900, col=1, lwd=4)
segments(x0=41,x1=41,y0=0,y1=25900, lwd=6, lty=1, col=rgb(red=31,green=73,blue=125,maxColorValue=255))
abline(h=c(0),lwd=4,col="white")
legend("topright", bty="n", legend=c("Black Lake", "Chignik Lake"), lwd=4, col=c("skyblue", "red"), text.col="white", cex=1.7)


## But not for 2010
# Genetics Escapement
DailyEscp.2010=as.numeric(readClipboard())

plot(y=DailyEscp.2010[1:94], x=4:sum(length(DailyEscp.2012),3), type="h", col="skyblue", cex=1.5, xlim=c(4,sum(length(DailyEscp.2012),3)), ylim=c(0,55000), bty="n", axes=FALSE, cex.lab=2, xlab="Date", ylab="Daily Escapement", col.lab="white", lwd=6)
rect(xleft=model.fit.2010$Day[which.min(abs(model.fit.2010$BirchWLS.Estimate-0.05))],
     ybottom=0,xright=model.fit.2010$Day[which.min(abs(model.fit.2010$BirchWLS.Estimate-0.95))],ytop=55000,col=1)
points(y=DailyEscp.2010[1:94], x=4:sum(length(DailyEscp.2012),3), type="h",col="skyblue",cex=1.5, lwd=6)
points(y=(DailyEscp.2010*model.fit.2010$BirchWLS.Estimate[4:sum(length(DailyEscp.2010),3)]), x=4:sum(length(DailyEscp.2010),3), type="h",col="red",cex=1.5, lwd=6)
axis(side=1,at=c(seq(1,91,by=10),99),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,at=seq(0,50000,by=5000),cex.axis=1.45,lwd=4,col="white",col.axis="white")
segments(x0=model.fit.2010$Day[which.min(abs(model.fit.2010$BirchWLS.Estimate-0.05))],
         x1=model.fit.2010$Day[which.min(abs(model.fit.2010$BirchWLS.Estimate-0.05))],y0=0,y1=54900, col=1, lwd=4)
segments(x0=model.fit.2010$Day[which.min(abs(model.fit.2010$BirchWLS.Estimate-0.95))], 
         x1=model.fit.2010$Day[which.min(abs(model.fit.2010$BirchWLS.Estimate-0.95))], y0=0, y1=54900, col=1, lwd=4)
segments(x0=41,x1=41,y0=0,y1=54900, lwd=6, lty=1, col=rgb(red=31,green=73,blue=125,maxColorValue=255))
abline(h=c(0),lwd=4,col="white")
legend("topright", bty="n", legend=c("Black Lake", "Chignik Lake"), lwd=4, col=c("skyblue", "red"), text.col="white", cex=1.7)



## Or for 2014
# Genetics Escapement
DailyEscp.2014=as.numeric(readClipboard())

plot(y=DailyEscp.2014[1:97], x=1:sum(length(DailyEscp.2012),3), type="h", col="skyblue", cex=1.5, xlim=c(4,sum(length(DailyEscp.2012),3)), ylim=c(0,26000), bty="n", axes=FALSE, cex.lab=2, xlab="Date", ylab="Daily Escapement", col.lab="white", lwd=6)
rect(xleft=model.fit.2014$Day[which.min(abs(model.fit.2014$BirchWLS.Estimate-0.05))],
     ybottom=0,xright=model.fit.2014$Day[which.min(abs(model.fit.2014$BirchWLS.Estimate-0.95))],ytop=26000,col=1)
points(y=DailyEscp.2014[1:97], x=1:sum(length(DailyEscp.2012),3), type="h",col="skyblue",cex=1.5, lwd=6)
points(y=(DailyEscp.2014*model.fit.2014$BirchWLS.Estimate[1:sum(length(DailyEscp.2014),3)]), x=1:sum(length(DailyEscp.2014),3), type="h",col="red",cex=1.5, lwd=6)
axis(side=1,at=c(seq(1,91,by=10),99),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,at=seq(0,50000,by=5000),cex.axis=1.45,lwd=4,col="white",col.axis="white")
segments(x0=model.fit.2014$Day[which.min(abs(model.fit.2014$BirchWLS.Estimate-0.05))],
         x1=model.fit.2014$Day[which.min(abs(model.fit.2014$BirchWLS.Estimate-0.05))],y0=0,y1=25900, col=1, lwd=4)
segments(x0=model.fit.2014$Day[which.min(abs(model.fit.2014$BirchWLS.Estimate-0.95))],
         x1=model.fit.2014$Day[which.min(abs(model.fit.2014$BirchWLS.Estimate-0.95))],y0=0,y1=25900, col=1, lwd=4)
segments(x0=41,x1=41,y0=0,y1=25900, lwd=6, lty=1, col=rgb(red=31,green=73,blue=125,maxColorValue=255))
abline(h=c(0),lwd=4,col="white")
legend("topright", bty="n", legend=c("Black Lake", "Chignik Lake"), lwd=4, col=c("skyblue", "red"), text.col="white", cex=1.7)



## Or for 2015
# Genetics Escapement
DailyEscp.2015=as.numeric(readClipboard())

plot(y=DailyEscp.2015[1:97], x=1:sum(length(DailyEscp.2012),3), type="h", col="skyblue", cex=1.5, xlim=c(4,sum(length(DailyEscp.2012),3)), ylim=c(0,60000), bty="n", axes=FALSE, cex.lab=2, xlab="Date", ylab="Daily Escapement", col.lab="white", lwd=6)
rect(xleft=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.05))],
     ybottom=0,xright=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.95))],ytop=60000,col=1)
points(y=DailyEscp.2015[1:97], x=1:sum(length(DailyEscp.2012),3), type="h",col="skyblue",cex=1.5, lwd=6)
points(y=(DailyEscp.2015*model.fit.2015$BirchWLS.Estimate[1:sum(length(DailyEscp.2015),3)]), x=1:sum(length(DailyEscp.2015),3), type="h",col="red",cex=1.5, lwd=6)
axis(side=1,at=c(seq(1,91,by=10),99),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,at=seq(0,50000,by=5000),cex.axis=1.45,lwd=4,col="white",col.axis="white")
segments(x0=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.05))],
         x1=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.05))],y0=0,y1=60000, col=1, lwd=4)
segments(x0=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.95))],
         x1=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.95))],y0=0,y1=60000, col=1, lwd=4)
segments(x0=41,x1=41,y0=0,y1=60000, lwd=6, lty=1, col=rgb(red=31,green=73,blue=125,maxColorValue=255))
abline(h=c(0),lwd=4,col="white")
legend("topright", bty="n", legend=c("Black Lake", "Chignik Lake"), lwd=4, col=c("skyblue", "red"), text.col="white", cex=1.7)



model.fit.2010$Day[which.min(abs(model.fit.2010$BirchWLS.Estimate-0.05))]
model.fit.2010$Day[which.min(abs(model.fit.2010$BirchWLS.Estimate-0.95))]



# Escapement goals
black.lower.escp.goal <- as.numeric(readClipboard())
black.lower.escp.goal <- rep(black.lower.escp.goal, each = 2)
writeClipboard(as.character(black.lower.escp.goal))

black.upper.escp.goal <- as.numeric(readClipboard())
black.upper.escp.goal <- rep(black.upper.escp.goal, each = 2)
writeClipboard(as.character(black.upper.escp.goal))


setwd("V:/Documents/_DOC_JunkDrawer/Power point presentations/Sockeye/Chignik/2015 CRAA")
detail2015 <- read.xlsx (file="Chignik inseason summary data.xlsx",sheetName="2015 Escapement Data",stringsAsFactors=FALSE,colClasses=c("Date",rep("numeric",16)))
str(detail2015)

plot(0, xlim = c(145, 266), ylim = c(0, 600000), bty="n", axes=FALSE, cex.lab=2, xlab="Date", ylab="Cumulative Escapement (1000's)", col.lab="white", lwd=6)
axis(side=1,at=c(seq(145,235,by=10),243, 255, 265),labels=rep("",13),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
text(x = c(seq(145,235,by=10),243, 255, 265), y = -50000, labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31", "9/12", "9/22"), cex = 1.5, col = "white", srt = 45, xpd = TRUE)
axis(side=2,at=seq(0,600000,by=100000), labels = seq(0, 600, by = 100), cex.axis=1.45,lwd=4,col="white",col.axis="white")

# Fishery Open Dates
points(x = detail2015$Day[detail2015$Fishery.Open], y = rep(600000, sum(detail2015$Fishery.Open)), type = "h", col = 1, lwd = 5)

legend("topleft", legend = c("Escapement Goal Range", "Black Lake", "Chignik Lake"), col = c("white", "skyblue", "red"), lwd = 5, bty = "n", text.col = "white", cex = 1.5)
# Early Escapement Goal Range
points(x = detail2015$Day[which(!is.na(detail2015$Black.Lower))], y = detail2015$Black.Lower[which(!is.na(detail2015$Black.Lower))], type = "l", col = "white", lwd = 5)
points(x = detail2015$Day[which(!is.na(detail2015$Black.Upper))], y = detail2015$Black.Upper[which(!is.na(detail2015$Black.Upper))], type = "l", col = "white", lwd = 5)
# July 4
points(x = detail2015$Day[which(!is.na(detail2015$CumJuly4Early))], y = detail2015$CumJuly4Early[which(!is.na(detail2015$CumJuly4Early))], type = "l", col = "skyblue", lwd = 5)
# Genetics
points(x = detail2015$Day, y = detail2015$CumGeneticsEarly, type = "l", col = "skyblue", lwd = 5)

# Late Escapement Goal Range
points(x = detail2015$Day[which(!is.na(detail2015$Chignik.Lower))], y = detail2015$Chignik.Lower[which(!is.na(detail2015$Chignik.Lower))], type = "l", col = "white", lwd = 5)
points(x = detail2015$Day[which(!is.na(detail2015$Chignik.Upper))], y = detail2015$Chignik.Upper[which(!is.na(detail2015$Chignik.Upper))], type = "l", col = "white", lwd = 5)
# July 4
points(x = detail2015$Day[which(!is.na(detail2015$CumJuly4Late))], y = detail2015$CumJuly4Late[which(!is.na(detail2015$CumJuly4Late))], type = "l", col = "red", lwd = 5)
# Genetics
points(x = detail2015$Day[which(!is.na(detail2015$CumGeneticsLate))], y = detail2015$CumGeneticsLate[which(!is.na(detail2015$CumGeneticsLate))], type = "l", col = "red", lwd = 5)
# X-axis
segments(x0 = 138, y0 = 0, x1 = 266, y1 = 0, col = "white", lwd = 5, xpd = TRUE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2017 CRAA Figures for Feb 2018 Meeting####
# Thu Feb 15 15:51:16 2018
date()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create function to calculate proprtions based on Birch's logistic function coefficients (kappa, gamma method vs. a, b)
# NOTE: kappa = b, gamma = -a/b
logistic=function(x,kappa,gamma){
  p=1/(1+exp(-kappa*(x-gamma)))
  return(p)
}

require(xlsx)
Birch.WLS.RegCoeffs=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="Birch.WLS.RegCoeffs",stringsAsFactors=FALSE,colClasses=c(rep("numeric",3)))
Birch.WLS.RegCoeffs$col <- c(colorRampPalette(c("black", "skyblue"))(7), "black")
Birch.WLS.RegCoeffs$lty <- c(1:7, 1)
# yr2017 <- as.numeric(readClipboard())  # 5/25-8/8 from Chignik escapement 2017.xlsx

png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2018/Figures/Logistic Regression Summary Plot.png", width = 7, height = 6.5, units = "in", res = 400)

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
yr2 <- Birch.WLS.RegCoeffs[8, ]
lines(x = 1:76, y = 1-logistic(x = 1:76, kappa = as.numeric(yr2[2]), gamma = as.numeric(yr2[3])), type = "l", lwd = 10, col = 1, lty = 1)  # y = 1-logistic(x = 1:76, kappa = as.numeric(yr2[2]), gamma = as.numeric(yr2[3]))

abline(h=c(0),lwd=4,col="white")
legend(x=0,y=1,bty="n",legend=2010:2017,lwd=c(rep(4, 7), 10),col=Birch.WLS.RegCoeffs$col, lty = Birch.WLS.RegCoeffs$lty,cex=1.5,text.col="white")

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2017
# Genetics Escapement
require(xlsx)
detail2017 <- read.xlsx (file="Chignik inseason summary data.xlsx",sheetName="2017 Escapement Data",stringsAsFactors=FALSE,colClasses=c("Date",rep("numeric",16)))


png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2018/Figures/2017 Detail Daily Escapement Plot.png", width = 8, height = 6.5, units = "in", res = 400)

par(mar=c(4.6, 4.6, 1.1, 1.1))
par(family="serif")
par(bg=rgb(red=31,green=73,blue=125,maxColorValue=255))

plot(y = detail2017$Escapement[3:101], x = 1:99, type = "h", col = "skyblue", cex = 1.5, xlim = c(0, 102), ylim = c(0, 30000), bty = "n", axes = FALSE, cex.lab = 2, xlab = "Date", ylab = "Daily Escapement", col.lab = "white", lwd = 6)
rect(xleft = which.min(abs(detail2017$GeneticsPropLate[3:101] - 0.05)), ybottom = 0, xright = which.min(abs(detail2017$GeneticsPropLate[3:101] - 0.95)), ytop = 30000, col = 1)
# points(x = c(1:99)[detail2017$Fishery.Open[3:101]], y = rep(600000, sum(detail2017$Fishery.Open[3:101])), type = "h", col = 1, lwd = 7)
segments(x0 = as.Date("2017-07-04") - as.Date("2017-05-24"), x1 = as.Date("2017-07-04") - as.Date("2017-05-24"), 
         lwd = 6, col = rgb(red=31,green=73,blue=125,maxColorValue=255), y0 = 0, y1 = 30000)
points(y = detail2017$Escapement[3:101], x = 1:99, type = "h", col = "skyblue", lwd = 6)
points(y = detail2017$DailyGeneticsLate[3:101], x = 1:99, type = "h", col = "red", lwd = 6)
axis(side=1,at=c(seq(1,91,by=10),99),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31"),cex.axis=1.45,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,at=seq(0,30000,by=5000),cex.axis=1.45,lwd=4,col="white",col.axis="white")
abline(h=c(0),lwd=4,col="white")
legend("topright", bty="n", legend=c("Black Lake", "Chignik Lake"), lwd=4, col=c("skyblue", "red"), text.col="white", cex=1.7)

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2017
# Cumulative Genetics Escapement
png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2018/Figures/2017 Detail Escapement Range Plot.png", width = 8, height = 6.5, units = "in", res = 400)
png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2018/Figures/2017 Detail July 4 Plot.png", width = 8, height = 6.5, units = "in", res = 400)
png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2018/Figures/2017 Detail Genetics Plot.png", width = 8, height = 6.5, units = "in", res = 400)
png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2018/Figures/2017 Detail Genetics Fishing Periods Plot.png", width = 8, height = 6.5, units = "in", res = 400)

par(mar=c(4.6, 4.6, 1.1, 1.1))
par(family="serif")
par(bg=rgb(red=31,green=73,blue=125,maxColorValue=255))

plot(0, xlim = c(144, 267), ylim = c(0, 450000), bty="n", axes=FALSE, cex.lab=2, xlab="Date", ylab="Cumulative Escapement (1000's)", col.lab="white", lwd=6)
axis(side=1,at=c(seq(146,236,by=10),244, 256, 267),labels=rep("",13),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
text(x = c(seq(146,236,by=10),244, 256, 267), y = -35000, labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31", "9/12", "9/23"), cex = 1.5, col = "white", srt = 45, xpd = TRUE)
axis(side=2,at=seq(0,450000,by=50000), labels = seq(0, 450, by = 50), cex.axis=1.3,lwd=4,col="white",col.axis="white")

# Fishery Open Dates
points(x = detail2017$Day[detail2017$Fishery.Open], y = rep(600000, sum(detail2017$Fishery.Open)), type = "h", col = 1, lwd = 6)

legend("topleft", legend = c("Escapement Goal Range", "Black Lake", "Chignik Lake"), col = c("white", "skyblue", "red"), lwd = 5, bty = "n", text.col = "white", cex = 1.3)
# Early Escapement Goal Range
points(x = detail2017$Day[which(!is.na(detail2017$Black.Lower))], y = detail2017$Black.Lower[which(!is.na(detail2017$Black.Lower))], type = "l", col = "white", lwd = 5)
points(x = detail2017$Day[which(!is.na(detail2017$Black.Upper))], y = detail2017$Black.Upper[which(!is.na(detail2017$Black.Upper))], type = "l", col = "white", lwd = 5)
# July 4
points(x = detail2017$Day[which(!is.na(detail2017$CumJuly4Early))], y = detail2017$CumJuly4Early[which(!is.na(detail2017$CumJuly4Early))], type = "l", col = "skyblue", lwd = 5)
# Genetics
points(x = detail2017$Day, y = detail2017$CumGeneticsEarly, type = "l", col = "skyblue", lwd = 5)

# Late Escapement Goal Range
points(x = detail2017$Day[which(!is.na(detail2017$Chignik.Lower))], y = detail2017$Chignik.Lower[which(!is.na(detail2017$Chignik.Lower))], type = "l", col = "white", lwd = 5)
points(x = detail2017$Day[which(!is.na(detail2017$Chignik.Upper))], y = detail2017$Chignik.Upper[which(!is.na(detail2017$Chignik.Upper))], type = "l", col = "white", lwd = 5)
# July 4
points(x = detail2017$Day[which(!is.na(detail2017$CumJuly4Late))], y = detail2017$CumJuly4Late[which(!is.na(detail2017$CumJuly4Late))], type = "l", col = "red", lwd = 5)
# Genetics
points(x = detail2017$Day[which(!is.na(detail2017$CumGeneticsLate))], y = detail2017$CumGeneticsLate[which(!is.na(detail2017$CumGeneticsLate))], type = "l", col = "red", lwd = 5)
# X-axis
segments(x0 = 138, y0 = 0, x1 = 267, y1 = 0, col = "white", lwd = 5, xpd = TRUE)

dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2016 CRAA Figures ####
detail2016 <- read.xlsx (file="Chignik inseason summary data.xlsx",sheetName="2016 Escapement Data",stringsAsFactors=FALSE,colClasses=c("Date",rep("numeric",16)))
str(detail2016)


png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2017/2016 Detail Escapement Range Plot.png", width = 8, height = 6.5, units = "in", res = 400)
png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2017/2016 Detail July 4 Plot.png", width = 8, height = 6.5, units = "in", res = 400)
png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2017/2016 Detail Genetics Plot.png", width = 8, height = 6.5, units = "in", res = 400)
png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2017/2016 Detail Genetics Fishing Periods Plot.png", width = 8, height = 6.5, units = "in", res = 400)

par(mar=c(4.6, 4.6, 1.1, 1.1))
par(family="serif")
par(bg=rgb(red=31,green=73,blue=125,maxColorValue=255))

plot(0, xlim = c(144, 267), ylim = c(0, 450000), bty="n", axes=FALSE, cex.lab=2, xlab="Date", ylab="Cumulative Escapement (1000's)", col.lab="white", lwd=6)
axis(side=1,at=c(seq(146,236,by=10),244, 256, 267),labels=rep("",13),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
text(x = c(seq(146,236,by=10),244, 256, 267), y = -35000, labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31", "9/12", "9/23"), cex = 1.5, col = "white", srt = 45, xpd = TRUE)
axis(side=2,at=seq(0,450000,by=50000), labels = seq(0, 450, by = 50), cex.axis=1.3,lwd=4,col="white",col.axis="white")

# Fishery Open Dates
points(x = detail2016$Day[detail2016$Fishery.Open], y = rep(600000, sum(detail2016$Fishery.Open)), type = "h", col = 1, lwd = 6)

legend("topleft", legend = c("Escapement Goal Range", "Black Lake", "Chignik Lake"), col = c("white", "skyblue", "red"), lwd = 5, bty = "n", text.col = "white", cex = 1.3)
# Early Escapement Goal Range
points(x = detail2016$Day[which(!is.na(detail2016$Black.Lower))], y = detail2016$Black.Lower[which(!is.na(detail2016$Black.Lower))], type = "l", col = "white", lwd = 5)
points(x = detail2016$Day[which(!is.na(detail2016$Black.Upper))], y = detail2016$Black.Upper[which(!is.na(detail2016$Black.Upper))], type = "l", col = "white", lwd = 5)
# July 4
points(x = detail2016$Day[which(!is.na(detail2016$CumJuly4Early))], y = detail2016$CumJuly4Early[which(!is.na(detail2016$CumJuly4Early))], type = "l", col = "skyblue", lwd = 5)
# Genetics
points(x = detail2016$Day, y = detail2016$CumGeneticsEarly, type = "l", col = "skyblue", lwd = 5)

# Late Escapement Goal Range
points(x = detail2016$Day[which(!is.na(detail2016$Chignik.Lower))], y = detail2016$Chignik.Lower[which(!is.na(detail2016$Chignik.Lower))], type = "l", col = "white", lwd = 5)
points(x = detail2016$Day[which(!is.na(detail2016$Chignik.Upper))], y = detail2016$Chignik.Upper[which(!is.na(detail2016$Chignik.Upper))], type = "l", col = "white", lwd = 5)
# July 4
points(x = detail2016$Day[which(!is.na(detail2016$CumJuly4Late))], y = detail2016$CumJuly4Late[which(!is.na(detail2016$CumJuly4Late))], type = "l", col = "red", lwd = 5)
# Genetics
points(x = detail2016$Day[which(!is.na(detail2016$CumGeneticsLate))], y = detail2016$CumGeneticsLate[which(!is.na(detail2016$CumGeneticsLate))], type = "l", col = "red", lwd = 5)
# X-axis
segments(x0 = 138, y0 = 0, x1 = 267, y1 = 0, col = "white", lwd = 5, xpd = TRUE)

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 2016
# Genetics Escapement
png(filename = "V:/Presentations/Regional/4_Westward/Sockeye/CRAA/2017/2016 Detail Daily Escapement Plot.png", width = 8, height = 6.5, units = "in", res = 400)

par(mar=c(4.6, 4.6, 1.1, 1.1))
par(family="serif")
par(bg=rgb(red=31,green=73,blue=125,maxColorValue=255))

plot(y = detail2016$Escapement[3:101], x = 1:99, type = "h", col = "skyblue", cex = 1.5, xlim = c(0, 102), ylim = c(0, 30000), bty = "n", axes = FALSE, cex.lab = 2, xlab = "Date", ylab = "Daily Escapement", col.lab = "white", lwd = 6)
rect(xleft = which.min(abs(detail2016$GeneticsPropLate[3:101] - 0.05)), ybottom = 0, xright = which.min(abs(detail2016$GeneticsPropLate[3:101] - 0.95)), ytop = 30000, col = 1)
# points(x = c(1:99)[detail2016$Fishery.Open[3:101]], y = rep(600000, sum(detail2016$Fishery.Open[3:101])), type = "h", col = 1, lwd = 7)
segments(x0 = as.Date("2016-07-04") - as.Date("2016-05-24"), x1 = as.Date("2016-07-04") - as.Date("2016-05-24"), 
         lwd = 6, col = rgb(red=31,green=73,blue=125,maxColorValue=255), y0 = 0, y1 = 30000)
points(y = detail2016$Escapement[3:101], x = 1:99, type = "h", col = "skyblue", lwd = 6)
points(y = detail2016$DailyGeneticsLate[3:101], x = 1:99, type = "h", col = "red", lwd = 6)
axis(side=1,at=c(seq(1,91,by=10),99),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31"),cex.axis=1.45,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,at=seq(0,30000,by=5000),cex.axis=1.45,lwd=4,col="white",col.axis="white")
abline(h=c(0),lwd=4,col="white")
legend("topright", bty="n", legend=c("Black Lake", "Chignik Lake"), lwd=4, col=c("skyblue", "red"), text.col="white", cex=1.7)

dev.off()





DailyEscp.2015=as.numeric(readClipboard())

plot(y=DailyEscp.2015[1:97], x=1:sum(length(DailyEscp.2012),3), type="h", col="skyblue", cex=1.5, xlim=c(4,sum(length(DailyEscp.2012),3)), ylim=c(0,60000), bty="n", axes=FALSE, cex.lab=2, xlab="Date", ylab="Daily Escapement", col.lab="white", lwd=6)
rect(xleft=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.05))],
     ybottom=0,xright=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.95))],ytop=60000,col=1)
points(y=DailyEscp.2015[1:97], x=1:sum(length(DailyEscp.2012),3), type="h",col="skyblue",cex=1.5, lwd=6)
points(y=(DailyEscp.2015*model.fit.2015$BirchWLS.Estimate[1:sum(length(DailyEscp.2015),3)]), x=1:sum(length(DailyEscp.2015),3), type="h",col="red",cex=1.5, lwd=6)
axis(side=1,at=c(seq(1,91,by=10),99),labels=c("5/25","6/4","6/14","6/24","7/4","7/14","7/24","8/3","8/13","8/23","8/31"),cex.axis=1.5,pos=0,lwd=4,col="white",col.axis="white")
axis(side=2,at=seq(0,50000,by=5000),cex.axis=1.45,lwd=4,col="white",col.axis="white")
segments(x0=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.05))],
         x1=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.05))],y0=0,y1=60000, col=1, lwd=4)
segments(x0=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.95))],
         x1=model.fit.2015$Day[which.min(abs(model.fit.2015$BirchWLS.Estimate-0.95))],y0=0,y1=60000, col=1, lwd=4)
segments(x0=41,x1=41,y0=0,y1=60000, lwd=6, lty=1, col=rgb(red=31,green=73,blue=125,maxColorValue=255))
abline(h=c(0),lwd=4,col="white")
legend("topright", bty="n", legend=c("Black Lake", "Chignik Lake"), lwd=4, col=c("skyblue", "red"), text.col="white", cex=1.7)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2015 Start Prior for inseason ####
mean(c(model.fit.2010[35, "BirchWLS.Estimate"],
       model.fit.2011[35, "BirchWLS.Estimate"],
       model.fit.2012[35, "BirchWLS.Estimate"],
       model.fit.2013[35, "BirchWLS.Estimate"],
       model.fit.2014[35, "BirchWLS.Estimate"]))











## Chignik Lake
# Inverse SD Weighted Logistic Model
Estimates.2010$InvSD=1/(Estimates.2010$ChigSD)
fit.w.SD=glm(ChigMean~Julian,data=Estimates.2010,weights=InvSD,family=binomial(logit))
summary(fit.w.SD)
round(predict(fit.w.SD,data.frame(Julian=seq(190,195,by=.1)),type="response"),3)
model.fit$Estimate.w.SD=round(predict(fit.w.SD,model.fit,type="response"),3)
lines(model.fit$Julian,model.fit$Estimate.w.SD,type="l",lwd=wid,lty=2,col=2)



## Black Lake
plot(BlackMean~Julian,data=Estimates.2010,pch=16,cex=1.5,xlim=c(160,220))

# Inverse Variance Weighted Logistic Model
Estimates.2010$InvVar=1/(Estimates.2010$BlackSD)^2
fit.w.Var=glm(BlackMean~Julian,data=Estimates.2010,weights=InvVar,family=binomial(logit))
summary(fit.w.Var)

# Predicted values for measured Julian Date
round(fit.w.Var$fitted,3)
lines(Estimates.2010$Julian,fit.w.Var$fitted,type="l",col=2,lwd=4)

# Determine daily predicted values for every Julian Date
model.fit=data.frame(Julian=c(160:220))
model.fit$Estimate.w.Var=round(predict(fit.w.Var,model.fit,type="response"),3)
wid=4
lines(model.fit$Julian,model.fit$Estimate.w.Var,type="l",lwd=wid,lty=2,col=2)
abline(h=c(0,1),lwd=4)


# Inverse SD Weighted Logistic Model
Estimates.2010$InvSD=1/(Estimates.2010$BlackSD)
fit.w.SD=glm(BlackMean~Julian,data=Estimates.2010,weights=InvSD,family=binomial(logit))
summary(fit.w.SD)
round(predict(fit.w.SD,data.frame(Julian=seq(190,195,by=.1)),type="response"),3)
model.fit$Estimate.w.SD=round(predict(fit.w.SD,model.fit,type="response"),3)
lines(model.fit$Julian,model.fit$Estimate.w.SD,type="l",lwd=wid,lty=2,col=2)


### Modify inseason image ####
Genetics.Estimates <- rbind(Estimates.2010, Estimates.2011, Estimates.2012, Estimates.2013, Estimates.2014)
Genetics.Estimates$Year <- c(rep("2010", nrow(Estimates.2010)), rep("2011", nrow(Estimates.2011)), rep("2012", nrow(Estimates.2012)), rep("2013", nrow(Estimates.2013)), rep("2014", nrow(Estimates.2014)))

plot(Genetics.Estimates$BlackMean ~ Genetics.Estimates$Day, pch = 16)
arrows(x0 = Genetics.Estimates$Day, y0 = Genetics.Estimates$Black5CI, x1 = Genetics.Estimates$Day, y1 = Genetics.Estimates$Black95CI, angle = 90, code = 3)
lines(x = dd$day, y = dd$pred)
lines(x = dd$day, y = dd$lower)
lines(x = dd$day, y = dd$upper)

overall.mod <- glm(BlackMean ~ Day, data = Genetics.Estimates, weights = InvVar, family = binomial(logit))
summary(overall.mod)

overall.year.mod <- glm(BlackMean ~ Day + Year, data = Genetics.Estimates, weights = InvVar, family = binomial(logit))
summary(overall.year.mod)

logodds.day <- coef(overall.mod)['Day']
prop.increase <- exp(logodds.day) - 1

days <- seq(from = 1, to = 68, by = 1)
overall.logistic <- data.frame("Day" = days)

preds <- predict(overall.mod, newdata = overall.logistic, se.fit = TRUE)

invl <- function(x){ 1/(1+exp(-x)) }

dd <- with(preds, data.frame(day = overall.logistic$Day,
                             lower = invl(fit - 1.64 * se.fit),
                             pred = invl(fit), 
                             upper = invl(fit + 1.64 * se.fit))
           )

dd

require(ggplot2)
ggplot(dd[34:56, ], aes(day, pred)) + 
  geom_ribbon(aes(x=day, ymin=lower, ymax=upper, fill="grey40", linetype=NA), 
              alpha=1) +
  geom_line(size=1) + 
  theme_bw() + 
  theme(legend.position="none") +
  labs(x = "Day", y = "Black (Early Run) percentage of sample")


# Jim's idea; determining overall parameters from current output ####
est.Day <- setNames(object = sapply(2010:2014, function (year) {coef(get(paste("Birch.fit.w.Var.", year, sep = "")))["Day"]}), nm = 2010:2014)

se.Day <- setNames(object = sapply(2010:2014, function (year) {coef(summary(get(paste("Birch.fit.w.Var.", year, sep = ""))))["Day", "Std. Error"]}), nm = 2010:2014)

CV.Day <- se.Day / est.Day

# Looks too tight, lack of faith in the SE around beta (slope)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Jim's better idea, sample from the MCMC output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))
require(beepr)
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")

ChignikGroups <- dget(file = "V:/WORK/Sockeye/Chignik/2013 Chignik Inseason/Objects/ChignikGroups.txt")

setwd("V:/DOC/Power point presentations/Sockeye/Chignik/2015 CRAA")
require(xlsx)
Estimates.2010=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2010",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
Estimates.2011=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2011",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
Estimates.2012=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2012",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
Estimates.2013=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2013",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))
Estimates.2014=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="2014",stringsAsFactors=FALSE,colClasses=c("character",rep("numeric",10)))

# Add model assumption of 0% May 25, 100% July 31, 
# Weighting for these points was determined to have a SD of 0.025 via Birch Foster
# I found 0.0095 when taking the mean of estimates close to 1 or 0 ...
# Perhpas a better way is to take the SD from 100% proof tests! which would give
AnchorSDs=read.xlsx(file="Chignik inseason summary data.xlsx",sheetName="AnchorSDs",stringsAsFactors=FALSE,colClasses=c(rep("numeric",2)))[,2:3]
EarlyAssumption=c(NA,145,1,NA,1.0,1+qnorm(0.05)*AnchorSDs[1,1],1,AnchorSDs[1,1],0,0,0+qnorm(0.95)*AnchorSDs[1,1],AnchorSDs[1,1])
LateAssumption=c(NA,212,68,NA,0,0,0+qnorm(0.95)*AnchorSDs[1,2],AnchorSDs[1,2],1.0,1+qnorm(0.05)*AnchorSDs[1,2],1,AnchorSDs[1,2])

for(i in 0:4){
  assign(paste("Estimates.201",i,sep=""),rbind(get(paste("Estimates.201",i,sep="")),EarlyAssumption,LateAssumption))
}

Estimates.2010[(length(Estimates.2010$Julian)-1):length(Estimates.2010$Julian),1]=c("5/25","7/31")
Estimates.2011[(length(Estimates.2011$Julian)-1):length(Estimates.2011$Julian),1]=c("5/25","7/31")
Estimates.2012[(length(Estimates.2012$Julian)-1):length(Estimates.2012$Julian),1]=c("5/25","7/31")
Estimates.2013[(length(Estimates.2013$Julian)-1):length(Estimates.2013$Julian),1]=c("5/25","7/31")
Estimates.2014[(length(Estimates.2014$Julian)-1):length(Estimates.2014$Julian),1]=c("5/25","7/31")

# Determine Inverse variance
Estimates.2010$InvVar=round(1/(Estimates.2010$ChigSD)^2,2)
Estimates.2011$InvVar=round(1/(Estimates.2011$ChigSD)^2,2)
Estimates.2012$InvVar=round(1/(Estimates.2012$ChigSD)^2,2)
Estimates.2013$InvVar=round(1/(Estimates.2013$ChigSD)^2,2)
Estimates.2014$InvVar=round(1/(Estimates.2014$ChigSD)^2,2)

## Read in estimates objects ####
# 2010
setwd("V:/WORK/Sockeye/Chignik/2010 Chignik Escapement")
mixvecs <- list.dirs(path = paste(getwd(),"/BAYES/Output", sep = ""), recursive = FALSE, full.names = FALSE)

SCHIG10_1_Jun14_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[1], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG10_2_Jun21_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[2], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG10_3_Jun27_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[3], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG10_4_Jul01_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[4], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG10_5_Jul05_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[5], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG10_6_Jul0809_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[6], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG10_7_Jul11_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[7], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG10_8_Jul14_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[8], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG10_9_Jul1819_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[9], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG10_10_Jul23_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[10], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG10_11_Jul30_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[11], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# 2011
setwd("V:/WORK/Sockeye/Chignik/2011 Chignik Escapement")
mixvecs <- list.dirs(path = paste(getwd(),"/BAYES/Output", sep = ""), recursive = FALSE, full.names = FALSE)

SCHIG11_1_Jun10_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[1], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG11_2_Jun17_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[2], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG11_3_Jun24_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[3], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG11_4_Jun28_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[4], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG11_5_Jul02_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[5], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG11_6_Jul05_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[6], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG11_7_Jul0910_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[7], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG11_8_Jul12_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[8], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG11_9_Jul14_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[9], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG11_10_Jul21_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[10], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
SCHIG11_11_Jul28_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = mixvecs[11], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# 2012
setwd("V:/WORK/Sockeye/Chignik/2012 Chignik Escapement/Estimate objects")
for(file in list.files()) {
  assign(x = unlist(strsplit(x = file, split = ".txt")), value = dget(file = file))
}; rm(file)

SCHIG12_1_Jun11_Estimates <- get(objects(pattern = "12_Estimates")[1])
SCHIG12_2_Jun18_Estimates <- get(objects(pattern = "12_Estimates")[2])
SCHIG12_3_Jun25_Estimates <- get(objects(pattern = "12_Estimates")[3])
SCHIG12_4_Jul01_Estimates <- get(objects(pattern = "12_Estimates")[4])
SCHIG12_5_Jul05_Estimates <- get(objects(pattern = "12_Estimates")[5])
SCHIG12_6_Jul0809_Estimates <- get(objects(pattern = "12_Estimates")[6])
SCHIG12_7_Jul11_Estimates <- get(objects(pattern = "12_Estimates")[7])
SCHIG12_8_Jul14_Estimates <- get(objects(pattern = "12_Estimates")[8])
SCHIG12_9_Jul17_Estimates <- get(objects(pattern = "12_Estimates")[9])
SCHIG12_10_Jul20_Estimates <- get(objects(pattern = "12_Estimates")[10])
SCHIG12_11_Jul28_Estimates <- get(objects(pattern = "12_Estimates")[11])

rm(list = objects(pattern = "12_Estimates")[-12])

# 2013
setwd("V:/WORK/Sockeye/Chignik/2013 Chignik Inseason/Estimate objects")
for(file in list.files()) {
  assign(x = unlist(strsplit(x = file, split = ".txt")), value = dget(file = file))
}; rm(file)

rm(SCHIG12_1_Jun27_Estimates)

# 2014
setwd("V:/WORK/Sockeye/Chignik/2014 Chignik Inseason/Kyle/Estimates objects")
for(file in list.files()) {
  assign(x = unlist(strsplit(x = file, split = ".txt")), value = dget(file = file))
}; rm(file)

# Put all estimates objects names into a big list, ordered by strata
estimates <- list("2010" = objects(pattern = "SCHIG10")[order(as.numeric((sapply(objects(pattern = "SCHIG10"), function(est) {unlist(strsplit(x = unlist(strsplit(x = est, split = "SCHIG10_"))[2], split = "_"))[1]}))))],
                  "2011" = objects(pattern = "SCHIG11")[order(as.numeric((sapply(objects(pattern = "SCHIG11"), function(est) {unlist(strsplit(x = unlist(strsplit(x = est, split = "SCHIG11_"))[2], split = "_"))[1]}))))],
                  "2012" = objects(pattern = "SCHIG12")[order(as.numeric((sapply(objects(pattern = "SCHIG12"), function(est) {unlist(strsplit(x = unlist(strsplit(x = est, split = "SCHIG12_"))[2], split = "_"))[1]}))))],
                  "2013" = objects(pattern = "SCHIG13")[order(as.numeric((sapply(objects(pattern = "SCHIG13"), function(est) {unlist(strsplit(x = unlist(strsplit(x = est, split = "SCHIG13_"))[2], split = "_"))[1]}))))],
                  "2014" = objects(pattern = "SCHIG14")[order(as.numeric((sapply(objects(pattern = "SCHIG14"), function(est) {unlist(strsplit(x = unlist(strsplit(x = est, split = "SCHIG14_"))[2], split = "_"))[1]}))))]
                  )

n.estimates <- lapply(estimates, length)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Extract coefficients at each MCMC rep within a year ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## With hard anchors
# 2010
begin <- proc.time()
coef.2010 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2010$Day),"BlackEst" = c(1, sapply(estimates[["2010"]], function(est) {get(est)[[2]][[1]][i, 1]}), 0))
  coef.2010[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)

# 2011
begin <- proc.time()
coef.2011 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2011$Day),"BlackEst" = c(1, sapply(estimates[["2011"]], function(est) {get(est)[[2]][[1]][i, 1]}), 0))
  coef.2011[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)

# 2012
begin <- proc.time()
coef.2012 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2012$Day),"BlackEst" = c(1, sapply(estimates[["2012"]], function(est) {get(est)[[2]][[1]][i, 1]}), 0))
  coef.2012[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)

# 2013
begin <- proc.time()
coef.2013 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2013$Day),"BlackEst" = c(1, sapply(estimates[["2013"]], function(est) {get(est)[[2]][[1]][i, 1]}), 0))
  coef.2013[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)

# 2014
begin <- proc.time()
coef.2014 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2014$Day),"BlackEst" = c(1, sapply(estimates[["2014"]], function(est) {get(est)[[2]][[1]][i, 1]}), 0))
  coef.2014[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)


coef.all <- list("2010" = coef.2010, "2011" = coef.2011, "2012" = coef.2012, "2013" = coef.2013, "2014" = coef.2014)
dput(x = coef.all, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/logistic.coef.all.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Without hard anchors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2010
begin <- proc.time()
coef.na.2010 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2010$Day)[-c(1, dim(Estimates.2010)[1])],"BlackEst" = c(sapply(estimates[["2010"]], function(est) {get(est)[[2]][[1]][i, 1]})))
  coef.na.2010[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)

rm(list = objects(pattern = "SCHIG10"))

# 2011
begin <- proc.time()
coef.na.2011 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2011$Day)[-c(1, dim(Estimates.2011)[1])],"BlackEst" = c(sapply(estimates[["2011"]], function(est) {get(est)[[2]][[1]][i, 1]})))
  coef.na.2011[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)

rm(list = objects(pattern = "SCHIG11"))
save.image("C:/Users/krshedd/Desktop/temp.RData")

# 2012
begin <- proc.time()
coef.na.2012 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2012$Day)[-c(1, dim(Estimates.2012)[1])],"BlackEst" = c(sapply(estimates[["2012"]], function(est) {get(est)[[2]][[1]][i, 1]})))
  coef.na.2012[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)

rm(list = objects(pattern = "SCHIG12"))
save.image("C:/Users/krshedd/Desktop/temp.RData")

# 2013
begin <- proc.time()
coef.na.2013 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2013$Day)[-c(1, dim(Estimates.2013)[1])],"BlackEst" = c(sapply(estimates[["2013"]], function(est) {get(est)[[2]][[1]][i, 1]})))
  coef.na.2013[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)

rm(list = objects(pattern = "SCHIG13"))
save.image("C:/Users/krshedd/Desktop/temp.RData")

# 2014
begin <- proc.time()
coef.na.2014 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2014$Day)[-c(1, dim(Estimates.2014)[1])],"BlackEst" = c(sapply(estimates[["2014"]], function(est) {get(est)[[2]][[1]][i, 1]})))
  coef.na.2014[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)

rm(list = objects(pattern = "SCHIG14"))
save.image("C:/Users/krshedd/Desktop/temp.RData")

coef.na.all <- list("2010" = coef.na.2010, "2011" = coef.na.2011, "2012" = coef.na.2012, "2013" = coef.na.2013, "2014" = coef.na.2014)
dput(x = coef.na.all, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/logistic.coef.na.all.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Determine global parameters ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## With 5/25 and 7/31 anchors

coef.all <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/logistic.coef.all.txt")
coef.all <- dget(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2015/Objects/logistic.coef.all.txt")
coef.na.all <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/logistic.coef.na.all.txt")

coef <- coef.all

alpha.global <- mean(sapply(coef, function(year) {mean(year[, "Intercept"])}))
beta.global <- mean(sapply(coef, function(year) {mean(year[, "Day"])}))

beta.var <- sum(sapply(coef, function(year) {(mean(year[, "Day"]) - beta.global)^2} ) / (length(coef) - 1)) +
  sum(sapply(coef, function(year) {var(year[, "Day"])} ) / length(coef)^2)

abs(sqrt(beta.var) / beta.global)

alpha.var <- sum(sapply(coef, function(year) {(mean(year[, "Intercept"]) - alpha.global)^2} ) / (length(coef) - 1)) +
  sum(sapply(coef, function(year) {var(year[, "Intercept"])} )) / length(coef)^2

abs(sqrt(alpha.var) / alpha.global)

# Determine 90% CI from global parameters ####
#logistic.f  <- function(x, kappa, gamma){
  #p = 1 / (1 + exp(-kappa * (x - gamma)))
  #return(p)
#}

#logistic.f(x = 33:55, kappa = beta.global, gamma = alpha.global)

logistic.f  <- function(x, alpha, beta){
  p = exp(alpha + beta * x) / (1 + exp(alpha + beta * x))
  return(p)
}

# Plot results
Genetics.Estimates <- rbind(Estimates.2010, Estimates.2011, Estimates.2012, Estimates.2013, Estimates.2014)
Genetics.Estimates$Year <- c(rep("2010", nrow(Estimates.2010)), rep("2011", nrow(Estimates.2011)), rep("2012", nrow(Estimates.2012)), rep("2013", nrow(Estimates.2013)), rep("2014", nrow(Estimates.2014)))


# Using standard error where n = 5 (so it is sqrt(var) / sqrt(5))
# alpha
alpha.sd <- sqrt(alpha.var)
alpha.se <- alpha.sd / sqrt(5)

alpha.5.sd <- alpha.global + qt(p = .05, df = 100000) * alpha.sd
alpha.95.sd <- alpha.global + qt(p = .95, df = 100000) * alpha.sd

alpha.5.se <- alpha.global + qt(p = .05, df = 100000) * alpha.se
alpha.95.se <- alpha.global + qt(p = .95, df = 100000) * alpha.se

# beta
beta.sd <- sqrt(beta.var)
beta.se <- beta.sd / sqrt(5)

beta.5.sd <- beta.global + qt(p = .05, df = 100000) * beta.sd
beta.95.sd <- beta.global + qt(p = .95, df = 100000) * beta.sd

beta.5.se <- beta.global + qt(p = .05, df = 100000) * beta.se
beta.95.se <- beta.global + qt(p = .95, df = 100000) * beta.se


plot(Genetics.Estimates$BlackMean ~ Genetics.Estimates$Day, pch = 16, bty = "n", xlab = "Day", ylab = "Proportion Black Lake", col = as.numeric(Genetics.Estimates$Year) - 2009, cex = 2, xlim = c(22, 68))
arrows(x0 = Genetics.Estimates$Day, y0 = Genetics.Estimates$Black5CI, x1 = Genetics.Estimates$Day, y1 = Genetics.Estimates$Black95CI, angle = 90, code = 3)
points(x = 22:68, y = logistic.f(x = 22:68, alpha = alpha.global, beta = beta.global), type = "l", lwd = 5, col = "red", add = TRUE)
points(x = 22:68, y = logistic.f(x = 22:68, alpha = alpha.5.se, beta = beta.5.se), type = "l", lwd = 5, lty = 2, add = TRUE)
points(x = 22:68, y = logistic.f(x = 22:68, alpha = alpha.95.se, beta = beta.95.se), type = "l", lwd = 5, lty = 2, add = TRUE)

legend("bottomleft", legend = 2010:2014, fill = 1:5, bty = "n")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Instead of calculating variance, just take the 5% and 95% of sorted alpha/beta posterior values ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coef.mat <- do.call(what = rbind, coef)
hist(coef.mat[, "Day"], col = 8)
hist(coef.mat[, "Intercept"], col = 8)

sort(coef.mat[, "Day"])[dim(coef.mat)[1] * 0.05]
sort(coef.mat[, "Day"])[dim(coef.mat)[1] * 0.95]

sort(coef.mat[, "Intercept"])[dim(coef.mat)[1] * 0.05]
sort(coef.mat[, "Intercept"])[dim(coef.mat)[1] * 0.95]


plot(Genetics.Estimates$BlackMean ~ Genetics.Estimates$Day, pch = 16, bty = "n", xlab = "Day", ylab = "Proportion Black Lake", col = as.numeric(Genetics.Estimates$Year) - 2009, cex = 2, xlim = c(22, 68))
arrows(x0 = Genetics.Estimates$Day, y0 = Genetics.Estimates$Black5CI, x1 = Genetics.Estimates$Day, y1 = Genetics.Estimates$Black95CI, angle = 90, code = 3)
points(x = 22:68, y = logistic.f(x = 22:68, alpha = alpha.global, beta = beta.global), type = "l", lwd = 5, xlim = c(21 ,67), ylim = c(0, 1), add = TRUE)
points(x = 22:68, y = logistic.f(x = 22:68, alpha = sort(coef.mat[, "Intercept"])[dim(coef.mat)[1] * 0.05], beta = sort(coef.mat[, "Day"])[dim(coef.mat)[1] * 0.05]), type = "l", lwd = 5, lty = 2, add = TRUE)
points(x = 22:68, y = logistic.f(x = 22:68, alpha = sort(coef.mat[, "Intercept"])[dim(coef.mat)[1] * 0.95], beta = sort(coef.mat[, "Day"])[dim(coef.mat)[1] * 0.95]), type = "l", lwd = 5, lty = 2, add = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Try modeling all years together, with equal weights across years ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Genetics.Estimates
Genetics.Estimates$Year <- as.factor(Genetics.Estimates$Year)

# Weight should be flat across years, and scaled to InvSD within years
yr.weight <- 1 / length(unique(Genetics.Estimates$Year))
fit.weights <- unlist(by(data = Genetics.Estimates, INDICES = Genetics.Estimates[, "Year"], function(yr) yr.weight * (1 / yr$BlackSD) / sum(yr$BlackSD), simplify = FALSE))

mod <- glm(BlackMean ~ Day, data = Genetics.Estimates, family = binomial(logit), weights = fit.weights)
mod2 <- glm(BlackMean ~ Day + Year, data = Genetics.Estimates, family = binomial(logit), weights = fit.weights)
mod3 <- glm(BlackMean ~ Day * Year, data = Genetics.Estimates, family = binomial(logit), weights = fit.weights)

plot(Genetics.Estimates$BlackMean ~ Genetics.Estimates$Day, pch = 16, bty = "n", xlab = "Day", ylab = "Proportion Black Lake", col = as.numeric(Genetics.Estimates$Year), cex = 2, xlim = c(22, 68))
arrows(x0 = Genetics.Estimates$Day, y0 = Genetics.Estimates$Black5CI, x1 = Genetics.Estimates$Day, y1 = Genetics.Estimates$Black95CI, angle = 90, code = 3)
legend("bottomleft", legend = unique(Genetics.Estimates$Year), fill = seq(unique(Genetics.Estimates$Year)), bty = "n")

# 90% Confidence interval
newdata <- data.frame("Day" = 22:68)
pred <- predict(object = mod, newdata = newdata, type = "link", se.fit = TRUE)
critval <- qt(p = .95, df = Inf)
upr <- pred$fit + (critval * pred$se.fit)
lwr <- pred$fit - (critval * pred$se.fit)
fit <- pred$fit

avg.mod <- cbind(22:68, mod$family$linkinv(fit), mod$family$linkinv(lwr), mod$family$linkinv(upr))
colnames(avg.mod) <- c("Day", "Fit", "Lower 5%", "Upper 95%")
lines(x = avg.mod[, "Day"], y = avg.mod[, "Fit"], lwd = 5)
lines(x = avg.mod[, "Day"], y = avg.mod[, "Lower 5%"], lwd = 5, lty = 2)
lines(x = avg.mod[, "Day"], y = avg.mod[, "Upper 95%"], lwd = 5, lty = 2)

# 90% Prediction interval
upr.pred <- pred$fit + (critval * sqrt(pred$se.fit^2 + 1)) # prediction interval adds 1 to the SE
lwr.pred <- pred$fit - (critval * sqrt(pred$se.fit^2 + 1)) # prediction interval adds 1 to the SE
avg.mod.pred <- cbind(22:68, mod$family$linkinv(fit), mod$family$linkinv(lwr.pred), mod$family$linkinv(upr.pred))
colnames(avg.mod.pred) <- c("Day", "Fit", "Lower 5%", "Upper 95%")
lines(x = avg.mod[, "Day"], y = avg.mod.pred[, "Lower 5%"], lwd = 5, lty = 3)
lines(x = avg.mod[, "Day"], y = avg.mod.pred[, "Upper 95%"], lwd = 5, lty = 3)

# Birch's method 
# I struggled a bit to come up with a good method of showing the model central tendency over the years we have.  
# If we used all years as a gauge the spread would be quite wide.  
# I Looked at the daily model percentages and averaged them over 2010-2014 (note the mean and median are virtually identical in this case).  
# Then as a conservative measure of the uncertainty, not so much model uncertainty but point estimate uncertainties based on our data, I created a lower and upper bound for the curve by +/- 2 SDs 
birch.lwr <- as.numeric(readClipboard())
birch.upr <- as.numeric(readClipboard())

lines(x = 22:68, y = birch.lwr, lwd = 5, lty = 4)
lines(x = 22:68, y = birch.upr, lwd = 5, lty = 4)

birch.lwr <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/birch.lwr.txt")
birch.upr <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/birch.upr.txt")

birch.lwr <- dget(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/Objects/birch.lwr.txt")
birch.upr <- dget(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/Objects/birch.upr.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot everything!!! ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(Genetics.Estimates$BlackMean ~ Genetics.Estimates$Day, pch = 16, bty = "n", xlab = "Day", ylab = "Proportion Black Lake", col = as.numeric(Genetics.Estimates$Year), cex = 2, xlim = c(22, 68))

legend("bottomleft", legend = unique(Genetics.Estimates$Year), fill = seq(unique(Genetics.Estimates$Year)), bty = "n")
# 90% Prediction Interval
polygon(x = c(22:68, 68:22), y = c(avg.mod.pred[, "Upper 95%"], rev(avg.mod.pred[, "Lower 5%"])), col = "grey80", border = FALSE)
# Birch's 95% something
polygon(x = c(22:68, 68:22), y = c(birch.upr, rev(birch.lwr)), col = "grey50", border = FALSE)
# 90% Confidence Interval
polygon(x = c(22:68, 68:22), y = c(avg.mod[, "Upper 95%"], rev(avg.mod[, "Lower 5%"])), col = "grey30", border = FALSE)
# Best fit from global model
lines(x = avg.mod[, "Day"], y = avg.mod[, "Fit"], lwd = 5)
# Best fit line from Jims mean parameter values (mean of posterior)
alpha.global <- mean(sapply(coef, function(year) {mean(year[, "Intercept"])}))
beta.global <- mean(sapply(coef, function(year) {mean(year[, "Day"])}))
lines(x = 22:68, y = logistic.f(x = 22:68, alpha = alpha.global, beta = beta.global), lwd = 5, lty = 2)
lines(x = 22:68, y = logistic.f(x = 22:68, alpha = median(coef.mat[, "Intercept"]), beta = median(coef.mat[, "Day"])), lwd = 5, lty = 3)

points(Genetics.Estimates$BlackMean ~ Genetics.Estimates$Day, pch = 16, col = as.numeric(Genetics.Estimates$Year), cex = 2)
arrows(x0 = Genetics.Estimates$Day, y0 = Genetics.Estimates$Black5CI, x1 = Genetics.Estimates$Day, y1 = Genetics.Estimates$Black95CI, angle = 90, code = 3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spoke with Hamachan and need to predict values from each regression (thin from 500K to 1K) and then calculate prediction interval around that ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

coef <- coef.all
coef.mat <- do.call(what = rbind, coef)
str(coef.mat)

logistic.f  <- function(x, alpha, beta){
  p = exp(alpha + beta * x) / (1 + exp(alpha + beta * x))
  return(p)
}

# Equal weight for all years
x <- 22:68
reps <- 10000
out <- matrix(data = NA, nrow = reps, ncol = length(x), dimnames = list(seq(reps), x))
samp <- dim(coef.mat)[1]/reps

for(i in seq(reps)) {
  out[i, ] <- logistic.f(x = x, alpha = coef.mat[i*samp, "Intercept"], beta = coef.mat[i*samp, "Day"])
}; beep(2)

# Weighting years by number of samples taken/year
x <- 22:68
reps <- 10000
out <- matrix(data = NA, nrow = reps, ncol = length(x), dimnames = list(seq(reps), x))

samp.yr <- table(Genetics.Estimates$Year) - 2 # subtract 2 to remove anchors

iter <- c(sample(x = 1:100000, size = reps * samp.yr[1] / sum(samp.yr), replace = FALSE),
          sample(x = 100001:200000, size = reps * samp.yr[2] / sum(samp.yr), replace = FALSE),
          sample(x = 200001:300000, size = reps * samp.yr[3] / sum(samp.yr), replace = FALSE),
          sample(x = 300001:400000, size = reps * samp.yr[4] / sum(samp.yr), replace = FALSE),
          sample(x = 400001:500000, size = reps * samp.yr[5] / sum(samp.yr), replace = FALSE))

for(i in seq(reps)) {
  out[i, ] <- logistic.f(x = x, alpha = coef.mat[iter[i], "Intercept"], beta = coef.mat[iter[i], "Day"])
}; beep(2)

# Check
str(out)
tail(out)

# 90% prediction interval
bayes.lwr <- apply(out, 2, function(day) {sort(day)[length(day)*0.05]})
bayes.upr <- apply(out, 2, function(day) {sort(day)[length(day)*0.5]})

# Bayesian 90% Prediction Interval
polygon(x = c(22:68, 68:22), y = c(bayes.upr, rev(bayes.lwr)), col = "orange", border = FALSE)
lines(x = 22:68, y = birch.upr, lwd = 3)
lines(x = 22:68, y = birch.lwr, lwd = 3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Jim suggests to plot log-odds so as to avoid "over-precision" in tails
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
year <- "2013"

P <- Genetics.Estimates[which(Genetics.Estimates$Year == year)[seq(samp.yr[year])], "BlackMean"]
t <- Genetics.Estimates[which(Genetics.Estimates$Year == year)[seq(samp.yr[year])], "Day"]

plot(log(P / (1-P)) ~ t, pch = 16, cex = 3, main = year)

lm(log(P / (1-P)) ~ t)


log((1-5e-7) / 5e-7)


# 2010
begin <- proc.time()
coef.logodds.2010 <- matrix(NA, nrow = 100000, ncol = 2, dimnames = list(1:100000, c("Intercept", "Day")))
for(i in seq(100000)) {
  rep.df <- data.frame("Day" = sort(Estimates.2010$Day),"BlackEst" = c(1, sapply(estimates[["2010"]], function(est) {get(est)[[2]][[1]][i, 1]}), 0))
  rep.df$BlackEst[which(rep.df$BlackEst == 1)] = 1 - 5e-7
  rep.df$BlackEst[which(rep.df$BlackEst == 0)] = 5e-7
  rep.df$logodd <- log(rep.df$BlackEst / (1 - rep.df$BlackEst))
  plot(BlackEst ~ Day, data = rep.df, pch = 16, cex = 3)
  fit <- glm(BlackEst ~ Day, data = rep.df, family = binomial(logit))
  newdata <- data.frame("Day" = 22:68)
  pred <- predict(object = fit, newdata = newdata, type = "response", se.fit = TRUE)
  lines(y = pred$fit, x = 22:68)
  coef.logodds.2010[i, ] <- coef(glm(BlackEst ~ Day, data = rep.df, family = binomial(logit)))
}; proc.time() - begin; beep(2)


# Show that log of odds is still bad for tails
p <- seq(from = 0.99, to = 0.01, by = 0.01)
plot(log(p / (1-p)) ~ rev(p), type = "l")
