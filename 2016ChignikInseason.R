## 2016 Chignik Escapement mixture analysis - Inseason!
## Kyle Shedd Mon Jun 13 09:43:05 2016

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 June 27 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016")
## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/2016ChignikInseason_1.RData")
## load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/2016ChignikInseason_1.RData")

#This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")
#source("C:/Users/thdann/Documents/R/Functions.GCL.R")

username = "krshedd"
#username = "thdann"

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = username, password = password)

## Save original LocusControl
loci24 <- LocusControl$locusnames
mito.loci24 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl_1_Jun27.txt")
dput(x = loci24, file = "Objects/loci24.txt")
dput(x = mito.loci24, file = "Objects/mito.loci24.txt")

## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "SCHIG16", username = username, password = password)
rm(username, password)

## Save unaltered .gcl's as back-up:
dput(x = SCHIG16.gcl, file = paste("Raw genotypes/SCHIG16_1_Jun27.gcl.txt", sep = ''))
SCHIG16.gcl$n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2016.1 June 27
unique(SCHIG16.gcl$attributes$CAPTURE_DATE)
SCHIG16_1_Jun27IDs <- AttributesToIDs.GCL(silly = "SCHIG16", attribute = "CAPTURE_DATE", matching = unique(SCHIG16.gcl$attributes$CAPTURE_DATE)[1])

SCHIG16_1_Jun27IDs <- list(na.omit(SCHIG16_1_Jun27IDs))
names(SCHIG16_1_Jun27IDs) <- "SCHIG16"

PoolCollections.GCL("SCHIG16", loci = loci24, IDs = SCHIG16_1_Jun27IDs, newname = "SCHIG16_1_Jun27")
SCHIG16_1_Jun27.gcl$n ## 190


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG16_1_Jun27_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG16_1_Jun27", loci24)
min(Original_SCHIG16_1_Jun27_SampleSizebyLocus) ## Good? Yes, 185.
apply(Original_SCHIG16_1_Jun27_SampleSizebyLocus, 1, min) / apply(Original_SCHIG16_1_Jun27_SampleSizebyLocus, 1, max)

## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG16_1_Jun27_ColSize <- SCHIG16_1_Jun27.gcl$n

## Remove individuals with >20% missing data
SCHIG16_1_Jun27_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG16_1_Jun27", proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG16_1_Jun27_PostMissLoci <- SCHIG16_1_Jun27.gcl$n

SCHIG16_1_Jun27_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, dimnames = list("SCHIG16_1_Jun27", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG16_1_Jun27_SampleSizes[, "Genotyped"] <- Original_SCHIG16_1_Jun27_ColSize
SCHIG16_1_Jun27_SampleSizes[, "Missing"] <- Original_SCHIG16_1_Jun27_ColSize - ColSize_SCHIG16_1_Jun27_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG16_1_Jun27_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SCHIG16_1_Jun27", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG16_1_Jun27_RemovedDups <- RemoveDups.GCL(SCHIG16_1_Jun27_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG16_1_Jun27_PostDuplicate <- SCHIG16_1_Jun27.gcl$n

SCHIG16_1_Jun27_SampleSizes[, "Duplicate"] <- ColSize_SCHIG16_1_Jun27_PostMissLoci-ColSize_SCHIG16_1_Jun27_PostDuplicate
SCHIG16_1_Jun27_SampleSizes[, "Final"] <- ColSize_SCHIG16_1_Jun27_PostDuplicate

write.xlsx(SCHIG16_1_Jun27_SampleSizes, file = "Output/SCHIG16_1_Jun27_SampleSizes.xlsx")
dput(x = SCHIG16_1_Jun27.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG16_1_Jun27_IDs.txt")

## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG16_1_Jun27", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG16_1_Jun27", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG16_1_Jun27", loci = loci24[-mito.loci24], path = "Genepop/SCHIG16_1_Jun27_23nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/SCHIG16_1_Jun27_23nuclearloci.txt.P")

# Plot Fis values
plot(sort(HWE[, "WC Fis"]), type = "h", lwd = 10, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40"); abline(h = 0, lwd = 5)

# Look at data for any markers out of HWE
HWE[HWE$PValue < 0.05, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Baseline 2012")

## Get baseline objects needed for MSA
Chignik7Populations <- dget(file = "Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "Objects/Groupvec7.txt")
Chignik24BaselineFormat <- dget(file = "Objects/Chignik24BaselineFormat.txt")
Inits <- dget(file = "Objects/Inits.txt")
loci24MSA <- dget(file = "Objects/loci24.txt")
ChignikGroups <- dget(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2013/Objects/ChignikGroups.txt")

## Defining the random seeds as the same as WASSIP mixtures for repeatability.
WASSIPSockeyeSeeds <- dget(file = "V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt")

setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG16_1_Jun27", loci = loci24MSA, IDs = NULL, mixname = "SCHIG16_1_Jun27",
                                          dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "Objects/ChignikMixtureFormat.txt")


## Start Priors are based upon best available information, Birch's logistic fit of the past five years for the 28th: 
birch.lwr <- dget(file = "Objects/birch.lwr.txt")
birch.upr <- dget(file = "Objects/birch.upr.txt")
start.prior.upr <- mean(x = c(birch.lwr["06/27"], birch.upr["06/27"]))
Chignik2016StartPriorWeights <- c(start.prior.upr, 1-start.prior.upr)
Chignik2016StartPrior <- Prior.GCL(groupvec = Groupvec7, groupweights = Chignik2016StartPriorWeights, minval = 0.01)
dput(x = Chignik2016StartPrior, file = "Objects/Chignik2016StartPrior.txt")

## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG16_1_Jun27", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2016StartPrior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

## Create output directory
dir.create("BAYES/Output/SCHIG16_1_Jun27")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG16_1_Jun27_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG16_1_Jun27",
                                                          prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
dput(x = SCHIG16_1_Jun27_Estimates, file="Estimates objects/SCHIG16_1_Jun27_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG16_1_Jun27_Estimates[[1]][[1]][, "GR"]

# View traceplot
par(mfrow = c(2, 1), mar = c(1.1, 4.1, 4.1, 2.1))
plot(SCHIG16_1_Jun27_Estimates$Output$SCHIG16_1_Jun27[seq(from = 1, to = 100000, by = 10), 1], type = "l", ylim = c(0, 1), ylab = "", main = "Black Lake", xlab = "")
abline(v = seq(from = 0, to = 10000, by = 2000))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(SCHIG16_1_Jun27_Estimates$Output$SCHIG16_1_Jun27[seq(from = 1, to = 100000, by = 10), 2], type = "l", ylim = c(0, 1), ylab = "", main = "Chignik Lake", xlab = "Repetitions", cex.lab = 1.2)
abline(v = seq(from = 0, to = 10000, by = 2000))
mtext(text = "Posterior", side = 2, outer = TRUE, line = -1, cex = 1.2)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG16_1_Jun27_Estimates$Stats[[1]],
           file="Estimates tables/SCHIG16_1_Jun27_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="Objects/ChignikInseasonReport.f.txt")

SCHIG16_1_Jun27_SampleSizes

ChignikInseasonReport.f(NewData = SCHIG16_1_Jun27_Estimates, Period = 1, NumSampled = 190, NumAnalyzed = SCHIG16_1_Jun27_SampleSizes[1, "Genotyped"],
                        Included = SCHIG16_1_Jun27_SampleSizes[1, "Final"], Month = "June", Day = 27)


## Prior for next round
Chignik2016Period2Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG16_1_Jun27_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2016Period2Prior, file = "Objects/Chignik2016Period2Prior.txt")

## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/2016ChignikInseason_1.RData")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 July 02 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016")
## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/2016ChignikInseason_2.RData")
## load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/2016ChignikInseason_2.RData")

#This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")
#source("C:/Users/thdann/Documents/R/Functions.GCL.R")

username = "krshedd"
#username = "thdann"

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = username, password = password)

## Save original LocusControl
loci24 <- LocusControl$locusnames
mito.loci24 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl_2_Jul02.txt")
dput(x = loci24, file = "Objects/loci24.txt")
dput(x = mito.loci24, file = "Objects/mito.loci24.txt")

## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "SCHIG16", username = username, password = password)
rm(username, password)

## Save unaltered .gcl's as back-up:
dput(x = SCHIG16.gcl, file = paste("Raw genotypes/SCHIG16_2_Jul02.gcl.txt", sep = ''))
SCHIG16.gcl$n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2016.1 June 27
unique(SCHIG16.gcl$attributes$CAPTURE_DATE)
SCHIG16_2_Jul02IDs <- AttributesToIDs.GCL(silly = "SCHIG16", attribute = "CAPTURE_DATE", 
                                          matching = unique(SCHIG16.gcl$attributes$CAPTURE_DATE)[2])

SCHIG16_2_Jul02IDs <- list(na.omit(SCHIG16_2_Jul02IDs))
names(SCHIG16_2_Jul02IDs) <- "SCHIG16"

PoolCollections.GCL("SCHIG16", loci = loci24, IDs = SCHIG16_2_Jul02IDs, newname = "SCHIG16_2_Jul02")
SCHIG16_2_Jul02.gcl$n ## 190
table(SCHIG16_2_Jul02.gcl$attributes$CAPTURE_DATE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG16_2_Jul02_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG16_2_Jul02", loci24)
min(Original_SCHIG16_2_Jul02_SampleSizebyLocus) ## Good? Not really, 142.
apply(Original_SCHIG16_2_Jul02_SampleSizebyLocus, 1, min) / SCHIG16_2_Jul02.gcl$n

Original_SCHIG16_2_Jul02_PercentbyLocus <- apply(Original_SCHIG16_2_Jul02_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(Original_SCHIG16_2_Jul02_PercentbyLocus < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_SCHIG16_2_Jul02_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG16_2_Jul02_ColSize <- SCHIG16_2_Jul02.gcl$n

## Remove individuals with >20% missing data
SCHIG16_2_Jul02_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG16_2_Jul02", proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG16_2_Jul02_PostMissLoci <- SCHIG16_2_Jul02.gcl$n

SCHIG16_2_Jul02_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, 
                                      dimnames = list("SCHIG16_2_Jul02", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG16_2_Jul02_SampleSizes[, "Genotyped"] <- Original_SCHIG16_2_Jul02_ColSize
SCHIG16_2_Jul02_SampleSizes[, "Missing"] <- Original_SCHIG16_2_Jul02_ColSize - ColSize_SCHIG16_2_Jul02_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG16_2_Jul02_DuplicateCheck95MinProportion <- 
  CheckDupWithinSilly.GCL(sillyvec = "SCHIG16_2_Jul02", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG16_2_Jul02_RemovedDups <- RemoveDups.GCL(SCHIG16_2_Jul02_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG16_2_Jul02_PostDuplicate <- SCHIG16_2_Jul02.gcl$n

SCHIG16_2_Jul02_SampleSizes[, "Duplicate"] <- ColSize_SCHIG16_2_Jul02_PostMissLoci-ColSize_SCHIG16_2_Jul02_PostDuplicate
SCHIG16_2_Jul02_SampleSizes[, "Final"] <- ColSize_SCHIG16_2_Jul02_PostDuplicate

SCHIG16_2_Jul02_SampleSizes
write.xlsx(SCHIG16_2_Jul02_SampleSizes, file = "Output/SCHIG16_2_Jul02_SampleSizes.xlsx")
dput(x = SCHIG16_2_Jul02.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG16_2_Jul02_IDs.txt")



## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG16_2_Jul02", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG16_2_Jul02", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG16_2_Jul02", 
                loci = loci24[-mito.loci24], 
                path = "Genepop/SCHIG16_2_Jul02_23nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/SCHIG16_2_Jul02_23nuclearloci.txt.P")

# Plot Fis values
plot(sort(HWE[, "WC Fis"]), type = "h", lwd = 10, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40"); abline(h = 0, lwd = 5)

# Look at data for any markers out of HWE
HWE[HWE$PValue < 0.05, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Baseline 2012")

## Get baseline objects needed for MSA
Chignik7Populations <- dget(file = "Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "Objects/Groupvec7.txt")
Chignik24BaselineFormat <- dget(file = "Objects/Chignik24BaselineFormat.txt")
Inits <- dget(file = "Objects/Inits.txt")
loci24MSA <- dget(file = "Objects/loci24.txt")
ChignikGroups <- dget(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2013/Objects/ChignikGroups.txt")

## Defining the random seeds as the same as WASSIP mixtures for repeatability.
WASSIPSockeyeSeeds <- dget(file = "V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt")

setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG16_2_Jul02", loci = loci24MSA, IDs = NULL, mixname = "SCHIG16_2_Jul02",
                                          dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "Objects/ChignikMixtureFormat.txt")


## Start Priors are based upon best available information, Birch's logistic fit of the past five years for the 28th: 
Chignik2016Period2Prior <- dget(file = "Objects/Chignik2016Period2Prior.txt")


## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG16_2_Jul02", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2016Period2Prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

## Create output directory
dir.create("BAYES/Output/SCHIG16_2_Jul02")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG16_2_Jul02_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG16_2_Jul02",
                               prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
dput(x = SCHIG16_2_Jul02_Estimates, file="Estimates objects/SCHIG16_2_Jul02_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG16_2_Jul02_Estimates[[1]][[1]][, "GR"]

# View traceplot
par(mfrow = c(2, 1), mar = c(1.1, 4.1, 4.1, 2.1))
plot(SCHIG16_2_Jul02_Estimates$Output$SCHIG16_2_Jul02[seq(from = 1, to = 100000, by = 10), 1], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Black Lake", xlab = "")
abline(v = seq(from = 0, to = 10000, by = 2000))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(SCHIG16_2_Jul02_Estimates$Output$SCHIG16_2_Jul02[seq(from = 1, to = 100000, by = 10), 2], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Chignik Lake", xlab = "Repetitions", cex.lab = 1.2)
abline(v = seq(from = 0, to = 10000, by = 2000))
mtext(text = "Posterior", side = 2, outer = TRUE, line = -1, cex = 1.2)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG16_2_Jul02_Estimates$Stats[[1]],
           file="Estimates tables/SCHIG16_2_Jul02_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="Objects/ChignikInseasonReport.f.txt")

SCHIG16_2_Jul02_SampleSizes


ChignikInseasonReport.f(NewData = SCHIG16_2_Jul02_Estimates, Period = 2, NumSampled = 190, 
                        NumAnalyzed = SCHIG16_2_Jul02_SampleSizes[1, "Genotyped"],
                        Included = SCHIG16_2_Jul02_SampleSizes[1, "Final"], Month = "July", Day = 2)


## Prior for next round
Chignik2016Period3Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG16_2_Jul02_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2016Period3Prior, file = "Objects/Chignik2016Period3Prior.txt")

## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/2016ChignikInseason_2.RData")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 July 07 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016")
## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/2016ChignikInseason_3.RData")
## load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/2016ChignikInseason_3.RData")

#This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")
#source("C:/Users/thdann/Documents/R/Functions.GCL.R")

username = "krshedd"
#username = "thdann"

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = username, password = password)

## Save original LocusControl
loci24 <- LocusControl$locusnames
mito.loci24 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl_3_Jul07.txt")
dput(x = loci24, file = "Objects/loci24.txt")
dput(x = mito.loci24, file = "Objects/mito.loci24.txt")

## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "SCHIG16", username = username, password = password)
rm(username, password)

## Save unaltered .gcl's as back-up:
dput(x = SCHIG16.gcl, file = paste("Raw genotypes/SCHIG16_3_Jul07.gcl.txt", sep = ''))
SCHIG16.gcl$n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2016.1 June 27
unique(SCHIG16.gcl$attributes$CAPTURE_DATE)
SCHIG16_3_Jul07IDs <- AttributesToIDs.GCL(silly = "SCHIG16", attribute = "CAPTURE_DATE", 
                                          matching = unique(SCHIG16.gcl$attributes$CAPTURE_DATE)[3])

SCHIG16_3_Jul07IDs <- list(na.omit(SCHIG16_3_Jul07IDs))
names(SCHIG16_3_Jul07IDs) <- "SCHIG16"

PoolCollections.GCL("SCHIG16", loci = loci24, IDs = SCHIG16_3_Jul07IDs, newname = "SCHIG16_3_Jul07")
SCHIG16_3_Jul07.gcl$n ## 190
table(SCHIG16_3_Jul07.gcl$attributes$CAPTURE_DATE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG16_3_Jul07_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG16_3_Jul07", loci24)
min(Original_SCHIG16_3_Jul07_SampleSizebyLocus) ## Good? Not really, 142.
apply(Original_SCHIG16_3_Jul07_SampleSizebyLocus, 1, min) / SCHIG16_3_Jul07.gcl$n

Original_SCHIG16_3_Jul07_PercentbyLocus <- apply(Original_SCHIG16_3_Jul07_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(Original_SCHIG16_3_Jul07_PercentbyLocus < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_SCHIG16_3_Jul07_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG16_3_Jul07_ColSize <- SCHIG16_3_Jul07.gcl$n

## Remove individuals with >20% missing data
SCHIG16_3_Jul07_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG16_3_Jul07", proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG16_3_Jul07_PostMissLoci <- SCHIG16_3_Jul07.gcl$n

SCHIG16_3_Jul07_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, 
                                      dimnames = list("SCHIG16_3_Jul07", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG16_3_Jul07_SampleSizes[, "Genotyped"] <- Original_SCHIG16_3_Jul07_ColSize
SCHIG16_3_Jul07_SampleSizes[, "Missing"] <- Original_SCHIG16_3_Jul07_ColSize - ColSize_SCHIG16_3_Jul07_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG16_3_Jul07_DuplicateCheck95MinProportion <- 
  CheckDupWithinSilly.GCL(sillyvec = "SCHIG16_3_Jul07", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG16_3_Jul07_RemovedDups <- RemoveDups.GCL(SCHIG16_3_Jul07_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG16_3_Jul07_PostDuplicate <- SCHIG16_3_Jul07.gcl$n

SCHIG16_3_Jul07_SampleSizes[, "Duplicate"] <- ColSize_SCHIG16_3_Jul07_PostMissLoci-ColSize_SCHIG16_3_Jul07_PostDuplicate
SCHIG16_3_Jul07_SampleSizes[, "Final"] <- ColSize_SCHIG16_3_Jul07_PostDuplicate

SCHIG16_3_Jul07_SampleSizes
write.xlsx(SCHIG16_3_Jul07_SampleSizes, file = "Output/SCHIG16_3_Jul07_SampleSizes.xlsx")
dput(x = SCHIG16_3_Jul07.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG16_3_Jul07_IDs.txt")



## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG16_3_Jul07", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG16_3_Jul07", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG16_3_Jul07", 
                loci = loci24[-mito.loci24], 
                path = "Genepop/SCHIG16_3_Jul07_23nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/SCHIG16_3_Jul07_23nuclearloci.txt.P")

# Plot Fis values
plot(sort(HWE[, "WC Fis"]), type = "h", lwd = 10, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40"); abline(h = 0, lwd = 5)

# Look at data for any markers out of HWE
HWE[HWE$PValue < 0.05, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Baseline 2012")

## Get baseline objects needed for MSA
Chignik7Populations <- dget(file = "Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "Objects/Groupvec7.txt")
Chignik24BaselineFormat <- dget(file = "Objects/Chignik24BaselineFormat.txt")
Inits <- dget(file = "Objects/Inits.txt")
loci24MSA <- dget(file = "Objects/loci24.txt")
ChignikGroups <- dget(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2013/Objects/ChignikGroups.txt")

## Defining the random seeds as the same as WASSIP mixtures for repeatability.
WASSIPSockeyeSeeds <- dget(file = "V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt")

setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG16_3_Jul07", loci = loci24MSA, IDs = NULL, mixname = "SCHIG16_3_Jul07",
                                          dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "Objects/ChignikMixtureFormat.txt")


## Start Priors are based upon best available information, Birch's logistic fit of the past five years for the 28th: 
Chignik2016Period3Prior <- dget(file = "Objects/Chignik2016Period3Prior.txt")


## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG16_3_Jul07", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2016Period3Prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

## Create output directory
dir.create("BAYES/Output/SCHIG16_3_Jul07")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG16_3_Jul07_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG16_3_Jul07",
                               prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
dput(x = SCHIG16_3_Jul07_Estimates, file="Estimates objects/SCHIG16_3_Jul07_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG16_3_Jul07_Estimates[[1]][[1]][, "GR"]

# View traceplot
par(mfrow = c(2, 1), mar = c(1.1, 4.1, 4.1, 2.1))
plot(SCHIG16_3_Jul07_Estimates$Output$SCHIG16_3_Jul07[seq(from = 1, to = 100000, by = 10), 1], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Black Lake", xlab = "")
abline(v = seq(from = 0, to = 10000, by = 2000))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(SCHIG16_3_Jul07_Estimates$Output$SCHIG16_3_Jul07[seq(from = 1, to = 100000, by = 10), 2], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Chignik Lake", xlab = "Repetitions", cex.lab = 1.2)
abline(v = seq(from = 0, to = 10000, by = 2000))
mtext(text = "Posterior", side = 2, outer = TRUE, line = -1, cex = 1.2)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG16_3_Jul07_Estimates$Stats[[1]],
           file="Estimates tables/SCHIG16_3_Jul07_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="Objects/ChignikInseasonReport.f.txt")

SCHIG16_3_Jul07_SampleSizes


ChignikInseasonReport.f(NewData = SCHIG16_3_Jul07_Estimates, Period = 3, NumSampled = 190, 
                        NumAnalyzed = SCHIG16_3_Jul07_SampleSizes[1, "Genotyped"],
                        Included = SCHIG16_3_Jul07_SampleSizes[1, "Final"], Month = "July", Day = 7)


## Prior for next round
Chignik2016Period4Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG16_3_Jul07_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2016Period4Prior, file = "Objects/Chignik2016Period4Prior.txt")

## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/2016ChignikInseason_3.RData")



require(beepr)
beep(8)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Hurray! You survived another year of Chignik in-season. Be ready for more big wind next year! ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sex Specific Estimates for Birch Foster ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016")
source("V:\\Analysis\\R files\\Scripts\\PROD\\Functions.GCL.r")
source("H:/R Source Scripts/Functions.GCL_KS.R")
# source("V:/DATA/R_GEN/tempGCL Source Scripts/IndividualAssignmentSummaryKS.GCL.R")


# Get the mixnames for each round of MSA
SCHIG16_mixnames <- list.dirs(path = "BAYES/Output", full.names = FALSE, recursive = FALSE)

#~~~~~~~~~~~~~~~~~~
## Get fish IDs for mixture fish
# dir.create("Final Fish IDs")
# 
# for(i in seq(SCHIG16_mixnames)){
#   load(paste("2016ChignikInseason_", i, ".RData", sep = ''))
#   dput(x = get(paste(SCHIG16_mixnames[i], ".gcl", sep = ''))$attributes$FK_FISH_ID, file = paste("Final Fish IDs/", SCHIG16_mixnames[i], "_IDs.txt", sep = ''))
# }; rm(list = ls(all = TRUE))
#~~~~~~~~~~~~~~~~~~




# Get baseline objects
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017")
Chignik7Populations <- dget(file = "Baseline 2012/Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "Baseline 2012/Objects/Groupvec7.txt")
loci24MSA <- dget(file = "Baseline 2012/Objects/loci24.txt")
ChignikGroups <- dget(file = "Mixtures/2013/Objects/ChignikGroups.txt")
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016")

# Get mixture objects
SCHIG16_IDs <- sapply(SCHIG16_mixnames, function(mix) {dget(file = paste("Final Fish IDs/", mix, "_IDs.txt", sep = ''))} )


# #!#!#!# Note that I did NOT thin the first round, but did thin by 100 the next 6 rounds #!#!#!#
# # Just trying for Round 1 on June 27 (no thinning)
# SCHIG16_1_Jun27_IndvAssign <- IndividualAssignmentSummaryKS.GCL(GroupNames = ChignikGroups, groupvec = Groupvec7, SCHIG16_mixnames = SCHIG16_mixnames[1], BAYESoutputDir = "BAYES/Output", nchains = 5, nreps = 40000, burn = 1/2, thin = 1)
# str(SCHIG16_1_Jun27_IndvAssign)
# hist(SCHIG16_1_Jun27_IndvAssign[[1]][, 1])
# 
# # Doing Round 2 - 6 now (thin = 100)
# SCHIG16_2_6_IndvAssign <- IndividualAssignmentSummaryKS.GCL(GroupNames = ChignikGroups, groupvec = Groupvec7, SCHIG16_mixnames = SCHIG16_mixnames[-1], BAYESoutputDir = "BAYES/Output", nchains = 5, nreps = 40000, burn = 1/2, thin = 100)
# str(SCHIG16_2_6_IndvAssign)
# 
# # Push all 6 mixtures into 1 list object
# SCHIG16_IndvAssign <- c(SCHIG16_1_Jun27_IndvAssign, SCHIG16_2_6_IndvAssign)
# str(SCHIG16_IndvAssign)
# sapply(SCHIG16_IndvAssign, dim)
# sapply(SCHIG16_IDs, length)
# 
# # Rename dimnames so they reflect Fish ID (dna_vial)
# for(i in seq(SCHIG16_mixnames)) {
#   dimnames(SCHIG16_IndvAssign[[i]])[[1]] = SCHIG16_IDs[[i]]
# }
# 
# 
# # Histogram of probablility Black Lake
# invisible(sapply(SCHIG16_mixnames, function(mixture) {hist(SCHIG16_IndvAssign[[mixture]][, "Black Lake"], col = 8, xlab = "Probability Black Lake (individuals)", main = mixture)} ))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get ASL metadata
ASL <- read.table(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/chignikISG2016ASLdata.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
str(ASL)
table(ASL$sex, ASL$sample_date)

ASL <- ASL[, c("dna_vial", "sample_date", "sex")]
str(ASL)

# SCHIG16_IndvAssign_Sex <- lapply(SCHIG16_IndvAssign, function(mix) {
#   cbind(mix, "sex" = ASL[match(rownames(mix), table = ASL$dna_vial), "sex"])
#   } )
# str(SCHIG16_IndvAssign_Sex)
# 
# # Plot male proportion over time
# plot(sapply(SCHIG16_IndvAssign_Sex, function(mix) {table(mix[, "sex"])["1"] / sum(table(mix[, "sex"])[c("1", "2")])} ), type = 'o', pch = 16, cex = 2, lwd = 3, ylim = c(0, 0.6), ylab = "Proportion Male")

SCHIG16_IDs_Sex <- sapply(SCHIG16_IDs, function(mix) {
  cbind("dna_vial" = mix, "sex" = ASL[match(mix, table = ASL$dna_vial), "sex"])
}, simplify = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Attempt Jim's Model

GroupNames = ChignikGroups; groupvec = Groupvec7; mixnames = SCHIG16_mixnames[2:6]; BAYESoutputDir = "BAYES/Output"; nchains = 5; nreps = 40000; burn = 1/2; thin = 100; threshhold=5e-7

# This is actually a more simplified version (no beta dist)
ChignikMSAbySex.f <- function(GroupNames, groupvec, mixnames, BAYESoutputDir, nchains = 5, nreps = 40000, burn = 1/2, thin = 100, threshhold = 5e-7) {
  
  names(mixnames)=mixnames
  
  if(length(nchains)==1){
    
    nchains=rep(nchains,length(mixnames))
    
    names(nchains)=mixnames
    
  }
  
  if(length(nreps)==1){
    
    nreps=rep(nreps,length(mixnames))
    
    names(nreps)=mixnames
    
  }
  
  myoutputdirs=sapply(mixnames,function(mixname){paste(BAYESoutputDir,mixname,sep="/")},simplify=FALSE)
  
  myfilenames=sapply(mixnames,function(mixname){paste(myoutputdirs[mixname],"/",mixnames[mixname],"Chain",seq(nchains[mixname]),"CLS.CLS",sep="")},simplify=FALSE)
  
  output=results=vector("list",length(mixnames))
  
  names(output)=names(results)=mixnames
  
  for(mixname in mixnames){
    
    testoutput=read.table(myfilenames[[mixname]][1],header=FALSE)
    
    mixsampsize=dim(testoutput)[1]/floor(nreps[mixname]/thin)
    
    chains=paste("Chain",seq(nchains[mixname]),sep="")
    
    skip=mixsampsize*floor(nreps[mixname]*burn/thin)
    
    output[[mixname]]=NULL
    
    results[[mixname]]=array(NA,c(mixsampsize,max(groupvec)),list(seq(mixsampsize),GroupNames))
    
    for(chain in chains){
      
      output[[mixname]]=rbind(output[[mixname]],read.table(myfilenames[[mixname]][match(chain,chains)],header=FALSE,skip=skip))
      
    }#chain 
    nits=nrow(output[[mixname]])/mixsampsize
    
    C=ncol(output[[mixname]])
    
    z=apply(output[[mixname]],1,which.max)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Divergence to a new function
    
    # Which individual are we talking about?
    indv.pos <- rep(seq(mixsampsize),nits)
    
    # Which iteration is this?
    its.pos <- rep(seq(nits), each = mixsampsize)
    
    # How many groups
    G <- max(groupvec)
    
    
    # Make z into an array, this should be much easier
    zgroup <- groupvec[z]
    
    zgroup.mat <- array(data = zgroup, dim = c(mixsampsize, nits), dimnames = list(SCHIG16_IDs[[mixname]], seq(nits)))
    
    # t(apply(zgroup.mat, 1, function(ind) {tabulate(ind, G) / nits} ))
    
    indv.males <- which(SCHIG16_IDs_Sex[[mixname]][, "sex"] == 1)
    indv.females <- which(SCHIG16_IDs_Sex[[mixname]][, "sex"] == 2)
    
    n.males <- length(indv.males)
    n.females <- length(indv.females)
    
    n.group <- apply(zgroup.mat, 2, function(ind) {tabulate(ind, G)})
    n.males.group <- apply(zgroup.mat[indv.males, ], 2, function(ind) {tabulate(ind, G)})
    n.females.group <- apply(zgroup.mat[indv.females, ], 2, function(ind) {tabulate(ind, G)})
    
  #   hist((n.males.group / n.males)[1, ], col = "blue", xlim = c(0.7, 1), breaks = seq(0.7, 1, 0.01), ylim = c(0, 200))
  #   hist((n.females.group / n.females)[1, ], col = "red", add = TRUE, breaks = seq(0.7, 1, 0.01))
  #   hist((n.group / mixsampsize)[1, ], col = "grey", breaks = seq(0.7, 1, 0.01))
  #   
  #   
  #   early.posterior <- (n.group / mixsampsize)[1, ]
  #   c(median(early.posterior), quantile(x = early.posterior, probs = c(0.05, 0.95)), sd(early.posterior), mean(early.posterior))
  #   
    
    # Accounts for sampling error
    Output.legit <- list(Male = cbind(apply(n.males.group, 2, function(iter) {rbinom(n = 1, size = sum(iter), prob = iter[1] / sum(iter))}) / n.males,
                                      apply(n.males.group, 2, function(iter) {rbinom(n = 1, size = sum(iter), prob = iter[2] / sum(iter))}) / n.males),
                         Female = cbind(apply(n.females.group, 2, function(iter) {rbinom(n = 1, size = sum(iter), prob = iter[1] / sum(iter))}) / n.females,
                                        apply(n.females.group, 2, function(iter) {rbinom(n = 1, size = sum(iter), prob = iter[2] / sum(iter))}) / n.females),
                         Overall = cbind(apply(n.group, 2, function(iter) {rbinom(n = 1, size = sum(iter), prob = iter[1] / sum(iter))}) / mixsampsize,
                                         apply(n.group, 2, function(iter) {rbinom(n = 1, size = sum(iter), prob = iter[2] / sum(iter))}) / mixsampsize))
  
  
    # Does not account for sampling error, ONLY genetic error
    Output <- list(Male = t(n.males.group / n.males),
                   Female = t(n.females.group / n.females),
                   Overall = t(n.group / mixsampsize))
  
  
    summary <- setNames(vector("list", 3), c("Male", "Female", "Overall"))
  
    for(sex in c("Male", "Female", "Overall")) {
      summary[[sex]] <- array(NA,c(G,6),dimnames=list(GroupNames,c("mean","sd","median","5%","95%","P=0")))
      
      summary[[sex]][GroupNames, 1] <- apply(Output[[sex]],2,mean)
      summary[[sex]][GroupNames, 2] <- apply(Output[[sex]],2,sd)
      summary[[sex]][GroupNames, 3] <- apply(Output[[sex]],2,median)
      summary[[sex]][GroupNames, 4] <- apply(Output[[sex]],2,quantile,probs=0.05)
      summary[[sex]][GroupNames, 5] <- apply(Output[[sex]],2,quantile,probs=0.95)
      summary[[sex]][GroupNames, 6] <- apply(Output[[sex]],2,function(clm){sum(clm<threshhold)/length(clm)})
    }
  
  
    summary.legit <- setNames(vector("list", 3), c("Male", "Female", "Overall"))
    
    for(sex in c("Male", "Female", "Overall")) {
      summary.legit[[sex]] <- array(NA,c(G,6),dimnames=list(GroupNames,c("mean","sd","median","5%","95%","P=0")))
      
      summary.legit[[sex]][GroupNames, 1] <- apply(Output.legit[[sex]],2,mean)
      summary.legit[[sex]][GroupNames, 2] <- apply(Output.legit[[sex]],2,sd)
      summary.legit[[sex]][GroupNames, 3] <- apply(Output.legit[[sex]],2,median)
      summary.legit[[sex]][GroupNames, 4] <- apply(Output.legit[[sex]],2,quantile,probs=0.05)
      summary.legit[[sex]][GroupNames, 5] <- apply(Output.legit[[sex]],2,quantile,probs=0.95)
      summary.legit[[sex]][GroupNames, 6] <- apply(Output.legit[[sex]],2,function(clm){sum(clm<threshhold)/length(clm)})
    }
  
  
    results[[mixname]] = summary.legit
    
  }#mixname 
  
  return(results)
  
}


SCHIG16_2_6_IndvAssign_Sex <- ChignikMSAbySex.f(GroupNames = ChignikGroups, groupvec = Groupvec7, mixnames = SCHIG16_mixnames[2:6], BAYESoutputDir = "BAYES/Output", nchains = 5, nreps = 40000, burn = 1/2, thin = 100, threshhold = 5e-7)
SCHIG16_1_Jun27_IndvAssign_Sex <- ChignikMSAbySex.f(GroupNames = ChignikGroups, groupvec = Groupvec7, mixnames = SCHIG16_mixnames[1], BAYESoutputDir = "BAYES/Output", nchains = 5, nreps = 40000, burn = 1/2, thin = 1, threshhold = 5e-7)

SCHIG16_IndvAssign_Sex_BAYES <- c(SCHIG16_1_Jun27_IndvAssign_Sex, SCHIG16_2_6_IndvAssign_Sex)
str(SCHIG16_IndvAssign_Sex_BAYES)

SCHIG16_sex_date.mat <- table(ASL$sample_date, ASL$sex)

dput(x = SCHIG16_IndvAssign_Sex_BAYES, file = "Objects/SCHIG16_IndvAssign_Sex_BAYES.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots results

ASL <- read.table(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/chignikISG2016ASLdata.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
ASL <- ASL[, c("dna_vial", "sample_date", "sex")]
SCHIG16_sex_date.mat <- table(ASL$sample_date, ASL$sex)
SCHIG16_IndvAssign_Sex_BAYES <- dget(file = "Objects/SCHIG16_IndvAssign_Sex_BAYES.txt")
GeneticEstimates2016 <- dget(file = "Objects/GeneticEstimates2016.txt")


library(devEMF)
emf(file = paste("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2016/SexSpecificEstimates2.emf", sep = ''), width = 7, height = 6, family = "Times")

# par(bg=rgb(red=31,green=73,blue=125,maxColorValue=255))
par(bg = "white")
par(col = "white", col.axis = "white", col.lab = "white")
cex.lab = 1.5
col.axs = "black"
col.fem = "deeppink1"
col.male = "navy"
col.overall = "grey40"

  
  plot(y = sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Overall"]]["Black Lake", "mean"]} ), x = GeneticEstimates2016[, "Day"], type = "o", col = "blue", cex = 2, lwd = 3, pch = 16, ylim = c(0, 1), xlab = "", ylab = "", main = "", bty = 'n', axes = FALSE, xlim = c(min(x = GeneticEstimates2016[, "Day"]) - 1, max(x = GeneticEstimates2016[, "Day"]) + 2), cex.lab = cex.lab)
  axis(side = 2, col = col.axs, col.axis = col.axs, cex.axis = cex.lab)
  axis(side = 1, at = GeneticEstimates2016[, "Day"], labels = rep('', 6), col = col.axs)
  text(x = GeneticEstimates2016[, "Day"], y = -0.15, labels = rownames(SCHIG16_sex_date.mat), srt = 45, xpd = TRUE, col = col.axs)
  mtext(text = "Date", side = 1, line = 4, cex = cex.lab, col = col.axs)
  mtext(text = "Proportion Black Lake", side = 2, line = 3, cex = cex.lab, col = col.axs)

  arrows(x0 = GeneticEstimates2016[, "Day"] + 0.2, x1 = GeneticEstimates2016[, "Day"] + 0.2, y0 = sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Male"]]["Black Lake", "5%"]} ), sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Male"]]["Black Lake", "95%"]} ), angle = 90, col = col.male, code = 3, lwd = 2, length = 0.15)
  arrows(x0 = GeneticEstimates2016[, "Day"] - 0.2, x1 = GeneticEstimates2016[, "Day"] - 0.2, y0 = sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Female"]]["Black Lake", "5%"]} ), sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Female"]]["Black Lake", "95%"]} ), angle = 90, col = col.fem, code = 3, lwd = 2, length = 0.15)
  arrows(x0 = GeneticEstimates2016[, "Day"], x1 = GeneticEstimates2016[, "Day"], y0 = sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Overall"]]["Black Lake", "5%"]} ), sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Overall"]]["Black Lake", "95%"]} ), angle = 90, col = col.overall, code = 3, lwd = 2, length = 0.15)
  
  points(x = GeneticEstimates2016[, "Day"], y = sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Male"]]["Black Lake", "mean"]} ), type = "o", col = col.male, cex = 2, lwd = 3, pch = 16)
  points(x = GeneticEstimates2016[, "Day"], y = sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Female"]]["Black Lake", "mean"]} ), type = "o", col = col.fem, cex = 2, lwd = 3, pch = 16)
  points(x = GeneticEstimates2016[, "Day"], y = sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Overall"]]["Black Lake", "mean"]} ), type = "o", col = col.overall, cex = 2, lwd = 3, pch = 16)
  
  text(x = GeneticEstimates2016[, "Day"]-1, y = rep(0.15, 6), labels = SCHIG16_sex_date.mat[, "1"], col = col.male, cex = SCHIG16_sex_date.mat[, "1"] / rowSums(SCHIG16_sex_date.mat[, c("1", "2")]) * 2, adj = 0)
  text(x = GeneticEstimates2016[, "Day"]-1, y = rep(0.1, 6), labels = SCHIG16_sex_date.mat[, "2"], col = col.fem, cex = SCHIG16_sex_date.mat[, "2"] / rowSums(SCHIG16_sex_date.mat[, c("1", "2")]) * 2, adj = 0.1)
  text(x = GeneticEstimates2016[, "Day"]-1, y = rep(0.05, 6), labels = rowSums(SCHIG16_sex_date.mat[, c("1", "2")]), col = col.overall, cex = 2, adj = 0.1)
  
  
  legend("topright", legend = c("Overall", "Male", "Female"), fill = c(col.overall, col.male, col.fem), bty = 'n', cex = 1.5, text.col = col.axs)
#   legend("topright", legend = c("Overall"), fill = c(col.overall), bty = 'n', cex = 1.5, text.col = col.axs)

dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


round(GeneticEstimates2016, 2)[, c("5%", "95%")]
round(t(sapply(SCHIG16_IndvAssign_Sex_BAYES, function(period) {period[["Overall"]]["Black Lake",]} )), 2)[, c("5%", "95%")]




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Genotype frequencies for 2016 for One_U1004-183 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Original_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG16", loci24)


## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG16_ColSize <- SCHIG16.gcl$n

## Remove individuals with >20% missing data
SCHIG16_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG16", loci = loci24, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG16_PostMissLoci <- SCHIG16.gcl$n

SCHIG16_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, dimnames = list("SCHIG16_", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG16_SampleSizes[, "Genotyped"] <- Original_SCHIG16_ColSize
SCHIG16_SampleSizes[, "Missing"] <- Original_SCHIG16_ColSize - ColSize_SCHIG16_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG16_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SCHIG16", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG16_RemovedDups <- RemoveDups.GCL(SCHIG16_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG16_PostDuplicate <- SCHIG16.gcl$n

SCHIG16_SampleSizes[, "Duplicate"] <- ColSize_SCHIG16_PostMissLoci-ColSize_SCHIG16_PostDuplicate
SCHIG16_SampleSizes[, "Final"] <- ColSize_SCHIG16_PostDuplicate


# Get IDs for each date

chig_dates <- unique(SCHIG16.gcl$attributes$CAPTURE_DATE)
chig_dates <- gsub(pattern = " AKDT", replacement = '', x = chig_dates)

SCHIG16_date_IDs <- NULL
for(i in seq(length(chig_dates))) {
  SCHIG16_date_IDs[[i]] <- AttributesToIDs.GCL(silly = "SCHIG16", attribute = "CAPTURE_DATE", matching = unique(SCHIG16.gcl$attributes$CAPTURE_DATE)[i])
}

names(SCHIG16_date_IDs) <- chig_dates

lapply(SCHIG16_date_IDs, function(ids) {table(SCHIG16.gcl$counts[ids, "One_U1004-183", 1])} )
lapply(SCHIG16_date_IDs, function(ids) {sum(table(SCHIG16.gcl$counts[ids, "One_U1004-183", 1]))} )
genofreqs <- sapply(SCHIG16_date_IDs, function(ids) {round(table(SCHIG16.gcl$counts[ids, "One_U1004-183", 1]) / sum(table(SCHIG16.gcl$counts[ids, "One_U1004-183", 1])), 2)} )

GeneticEstimates2016 <- dget(file = "Objects/GeneticEstimates2016.txt")

par(mar = c(7.1, 5.1, 1.1, 1.1))
plot(y = genofreqs["0", ], x = GeneticEstimates2016[, "Day"], col = "blue", ylim = c(0, 1), type = "o", pch = 16, cex = 2, lwd = 3, axes = FALSE, bty = 'n', xlab = '', ylab = "Genotype Frequency", cex.lab = 2)
mtext(text = "Date", side = 1, line = 6, cex =2)
axis(side = 2, cex.axis = 1.5)
axis(side = 1, at = GeneticEstimates2016[, "Day"], labels = rep('', 6))
text(x = GeneticEstimates2016[, "Day"], y = -0.15, labels = chig_dates <- gsub(pattern = "2016-", replacement = '', x = chig_dates), xpd = TRUE, srt = 45, cex = 2)
points(y = genofreqs["1", ], x = GeneticEstimates2016[, "Day"], col = "green", type = "o", pch = 16, cex = 2, lwd = 3)
points(y = genofreqs["2", ], x = GeneticEstimates2016[, "Day"], col = "red", type = "o", pch = 16, cex = 2, lwd = 3)
legend("topleft", legend = c("aa", "Aa", "AA"), fill = c("blue", "green", "red"), bty = 'n', title = "One_U1004-183", cex = 2, xjust = 0, )


barplot(height = genofreqs, col = c("blue", "green", "red"))
