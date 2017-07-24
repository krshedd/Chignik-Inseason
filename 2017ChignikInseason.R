## 2017 Chignik Escapement mixture analysis - Inseason!
## Kyle Shedd Mon Jun 13 09:43:05 2017

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 June 27 2017 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_1.RData")
## load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_1.RData")

#This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")
#source("C:/Users/thdann/Documents/R/Functions.GCL.R")

username = "krshedd"

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = username, password = password)

## Save original LocusControl
loci24 <- LocusControl$locusnames
mito.loci24 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl_1_Jun25.txt")
dput(x = loci24, file = "Objects/loci24.txt")
dput(x = mito.loci24, file = "Objects/mito.loci24.txt")

## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "SCHIG17", username = username, password = password)
rm(username, password)

## Save unaltered .gcl's as back-up:
dput(x = SCHIG17.gcl, file = paste("Raw genotypes/SCHIG17_1_Jun25.gcl.txt", sep = ''))
SCHIG17.gcl$n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2017.1 June 27
unique(SCHIG17.gcl$attributes$CAPTURE_DATE)
SCHIG17_1_Jun25IDs <- AttributesToIDs.GCL(silly = "SCHIG17", attribute = "CAPTURE_DATE", matching = unique(SCHIG17.gcl$attributes$CAPTURE_DATE)[1:2])

SCHIG17_1_Jun25IDs <- list(na.omit(SCHIG17_1_Jun25IDs))
names(SCHIG17_1_Jun25IDs) <- "SCHIG17"

PoolCollections.GCL("SCHIG17", loci = loci24, IDs = SCHIG17_1_Jun25IDs, newname = "SCHIG17_1_Jun25")
SCHIG17_1_Jun25.gcl$n ## 190


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG17_1_Jun25_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG17_1_Jun25", loci24)
min(Original_SCHIG17_1_Jun25_SampleSizebyLocus)  # Good? Yes, 179.
apply(Original_SCHIG17_1_Jun25_SampleSizebyLocus, 1, min) / apply(Original_SCHIG17_1_Jun25_SampleSizebyLocus, 1, max)

## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG17_1_Jun25_ColSize <- SCHIG17_1_Jun25.gcl$n

## Remove individuals with >20% missing data
SCHIG17_1_Jun25_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG17_1_Jun25", proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG17_1_Jun25_PostMissLoci <- SCHIG17_1_Jun25.gcl$n

SCHIG17_1_Jun25_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, dimnames = list("SCHIG17_1_Jun25", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG17_1_Jun25_SampleSizes[, "Genotyped"] <- Original_SCHIG17_1_Jun25_ColSize
SCHIG17_1_Jun25_SampleSizes[, "Missing"] <- Original_SCHIG17_1_Jun25_ColSize - ColSize_SCHIG17_1_Jun25_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG17_1_Jun25_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SCHIG17_1_Jun25", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG17_1_Jun25_RemovedDups <- RemoveDups.GCL(SCHIG17_1_Jun25_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG17_1_Jun25_PostDuplicate <- SCHIG17_1_Jun25.gcl$n

SCHIG17_1_Jun25_SampleSizes[, "Duplicate"] <- ColSize_SCHIG17_1_Jun25_PostMissLoci-ColSize_SCHIG17_1_Jun25_PostDuplicate
SCHIG17_1_Jun25_SampleSizes[, "Final"] <- ColSize_SCHIG17_1_Jun25_PostDuplicate

write.xlsx(SCHIG17_1_Jun25_SampleSizes, file = "Output/SCHIG17_1_Jun25_SampleSizes.xlsx")
dput(x = SCHIG17_1_Jun25.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG17_1_Jun25_IDs.txt")

## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG17_1_Jun25", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG17_1_Jun25", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG17_1_Jun25", loci = loci24[-mito.loci24], path = "Genepop/SCHIG17_1_Jun25_23nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/SCHIG17_1_Jun25_23nuclearloci.txt.P")

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

setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG17_1_Jun25", loci = loci24MSA, IDs = NULL, mixname = "SCHIG17_1_Jun25",
                                          dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "Objects/ChignikMixtureFormat.txt")


## Start Priors are based upon best available information, mean by day logistic fit of the past five years for the 25th: 
avg.transition.lwr <- dget(file = "Objects/avg.transition.lwr.txt")
avg.transition.upr <- dget(file = "Objects/avg.transition.upr.txt")
start.prior.upr <- mean(x = c(avg.transition.lwr["06/25"], avg.transition.upr["06/25"]))
Chignik2017StartPriorWeights <- c(start.prior.upr, 1-start.prior.upr)
Chignik2017StartPrior <- Prior.GCL(groupvec = Groupvec7, groupweights = Chignik2017StartPriorWeights, minval = 0.01)
dput(x = Chignik2017StartPrior, file = "Objects/Chignik2017StartPrior.txt")

## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG17_1_Jun25", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2017StartPrior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

## Create output directory
dir.create("BAYES/Output/SCHIG17_1_Jun25")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG17_1_Jun25_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG17_1_Jun25",
                                                          prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
dput(x = SCHIG17_1_Jun25_Estimates, file="Estimates objects/SCHIG17_1_Jun25_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG17_1_Jun25_Estimates[[1]][[1]][, "GR"]

# View traceplot
par(mfrow = c(2, 1), mar = c(1.1, 4.1, 4.1, 2.1))
plot(SCHIG17_1_Jun25_Estimates$Output$SCHIG17_1_Jun25[seq(from = 1, to = 100000, by = 10), 1], type = "l", ylim = c(0, 1), ylab = "", main = "Black Lake", xlab = "")
abline(v = seq(from = 0, to = 10000, by = 2000))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(SCHIG17_1_Jun25_Estimates$Output$SCHIG17_1_Jun25[seq(from = 1, to = 100000, by = 10), 2], type = "l", ylim = c(0, 1), ylab = "", main = "Chignik Lake", xlab = "Repetitions", cex.lab = 1.2)
abline(v = seq(from = 0, to = 10000, by = 2000))
mtext(text = "Posterior", side = 2, outer = TRUE, line = -1, cex = 1.2)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG17_1_Jun25_Estimates$Stats[[1]],
           file="Estimates tables/SCHIG17_1_Jun25_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="Objects/ChignikInseasonReport.f.txt")

SCHIG17_1_Jun25_SampleSizes

ChignikInseasonReport.f(NewData = SCHIG17_1_Jun25_Estimates, Period = 1, NumSampled = 190, NumAnalyzed = SCHIG17_1_Jun25_SampleSizes[1, "Genotyped"],
                        Included = SCHIG17_1_Jun25_SampleSizes[1, "Final"], Month = "June", Day = 25)


## Prior for next round
Chignik2017Period2Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG17_1_Jun25_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2017Period2Prior, file = "Objects/Chignik2017Period2Prior.txt")

## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_1.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 July 01 2017 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_2.RData")
## load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_2.RData")

#This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

username = "krshedd"

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = username, password = password)

## Save original LocusControl
loci24 <- LocusControl$locusnames
mito.loci24 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl_2_Jul01.txt")
dput(x = loci24, file = "Objects/loci24.txt")
dput(x = mito.loci24, file = "Objects/mito.loci24.txt")

## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "SCHIG17", username = username, password = password)
rm(username, password)

## Save unaltered .gcl's as back-up:
dput(x = SCHIG17.gcl, file = paste("Raw genotypes/SCHIG17_2_Jul01.gcl.txt", sep = ''))
SCHIG17.gcl$n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2017.2 July 1
unique(SCHIG17.gcl$attributes$CAPTURE_DATE)
SCHIG17_2_Jul01IDs <- AttributesToIDs.GCL(silly = "SCHIG17", attribute = "CAPTURE_DATE", 
                                          matching = unique(SCHIG17.gcl$attributes$CAPTURE_DATE)[3])

SCHIG17_2_Jul01IDs <- list(na.omit(SCHIG17_2_Jul01IDs))
names(SCHIG17_2_Jul01IDs) <- "SCHIG17"

PoolCollections.GCL("SCHIG17", loci = loci24, IDs = SCHIG17_2_Jul01IDs, newname = "SCHIG17_2_Jul01")
SCHIG17_2_Jul01.gcl$n ## 190
table(SCHIG17_2_Jul01.gcl$attributes$CAPTURE_DATE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG17_2_Jul01_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG17_2_Jul01", loci24)
min(Original_SCHIG17_2_Jul01_SampleSizebyLocus) ## Good? Not really, 182.
apply(Original_SCHIG17_2_Jul01_SampleSizebyLocus, 1, min) / SCHIG17_2_Jul01.gcl$n  # 0.96

Original_SCHIG17_2_Jul01_PercentbyLocus <- apply(Original_SCHIG17_2_Jul01_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(Original_SCHIG17_2_Jul01_PercentbyLocus < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_SCHIG17_2_Jul01_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG17_2_Jul01_ColSize <- SCHIG17_2_Jul01.gcl$n

## Remove individuals with >20% missing data
SCHIG17_2_Jul01_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG17_2_Jul01", proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG17_2_Jul01_PostMissLoci <- SCHIG17_2_Jul01.gcl$n

SCHIG17_2_Jul01_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, 
                                      dimnames = list("SCHIG17_2_Jul01", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG17_2_Jul01_SampleSizes[, "Genotyped"] <- Original_SCHIG17_2_Jul01_ColSize
SCHIG17_2_Jul01_SampleSizes[, "Missing"] <- Original_SCHIG17_2_Jul01_ColSize - ColSize_SCHIG17_2_Jul01_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG17_2_Jul01_DuplicateCheck95MinProportion <- 
  CheckDupWithinSilly.GCL(sillyvec = "SCHIG17_2_Jul01", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG17_2_Jul01_RemovedDups <- RemoveDups.GCL(SCHIG17_2_Jul01_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG17_2_Jul01_PostDuplicate <- SCHIG17_2_Jul01.gcl$n

SCHIG17_2_Jul01_SampleSizes[, "Duplicate"] <- ColSize_SCHIG17_2_Jul01_PostMissLoci-ColSize_SCHIG17_2_Jul01_PostDuplicate
SCHIG17_2_Jul01_SampleSizes[, "Final"] <- ColSize_SCHIG17_2_Jul01_PostDuplicate

SCHIG17_2_Jul01_SampleSizes
write.xlsx(SCHIG17_2_Jul01_SampleSizes, file = "Output/SCHIG17_2_Jul01_SampleSizes.xlsx")
dput(x = SCHIG17_2_Jul01.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG17_2_Jul01_IDs.txt")



## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG17_2_Jul01", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG17_2_Jul01", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG17_2_Jul01", 
                loci = loci24[-mito.loci24], 
                path = "Genepop/SCHIG17_2_Jul01_23nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/SCHIG17_2_Jul01_23nuclearloci.txt.P")

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

setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG17_2_Jul01", loci = loci24MSA, IDs = NULL, mixname = "SCHIG17_2_Jul01",
                                          dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "Objects/ChignikMixtureFormat.txt")


## Rolling prior: 
Chignik2017Period2Prior <- dget(file = "Objects/Chignik2017Period2Prior.txt")


## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG17_2_Jul01", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2017Period2Prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

## Create output directory
dir.create("BAYES/Output/SCHIG17_2_Jul01")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG17_2_Jul01_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG17_2_Jul01",
                               prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
dput(x = SCHIG17_2_Jul01_Estimates, file="Estimates objects/SCHIG17_2_Jul01_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG17_2_Jul01_Estimates[[1]][[1]][, "GR"]

# View traceplot
par(mfrow = c(2, 1), mar = c(1.1, 4.1, 4.1, 2.1))
plot(SCHIG17_2_Jul01_Estimates$Output$SCHIG17_2_Jul01[seq(from = 1, to = 100000, by = 10), 1], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Black Lake", xlab = "")
abline(v = seq(from = 0, to = 10000, by = 2000))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(SCHIG17_2_Jul01_Estimates$Output$SCHIG17_2_Jul01[seq(from = 1, to = 100000, by = 10), 2], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Chignik Lake", xlab = "Repetitions", cex.lab = 1.2)
abline(v = seq(from = 0, to = 10000, by = 2000))
mtext(text = "Posterior", side = 2, outer = TRUE, line = -1, cex = 1.2)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG17_2_Jul01_Estimates$Stats[[1]],
           file="Estimates tables/SCHIG17_2_Jul01_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="Objects/ChignikInseasonReport.f.txt")

SCHIG17_2_Jul01_SampleSizes


ChignikInseasonReport.f(NewData = SCHIG17_2_Jul01_Estimates, Period = 2, NumSampled = 190, 
                        NumAnalyzed = SCHIG17_2_Jul01_SampleSizes[1, "Genotyped"],
                        Included = SCHIG17_2_Jul01_SampleSizes[1, "Final"], Month = "July", Day = 1)


## Prior for next round
Chignik2017Period3Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG17_2_Jul01_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2017Period3Prior, file = "Objects/Chignik2017Period3Prior.txt")

## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_2.RData")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 July 07/08 2017 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_3.RData")
## load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_3.RData")

#This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

username = "krshedd"

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = username, password = password)

## Save original LocusControl
loci24 <- LocusControl$locusnames
mito.loci24 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl_3_Jul0708.txt")
dput(x = loci24, file = "Objects/loci24.txt")
dput(x = mito.loci24, file = "Objects/mito.loci24.txt")

## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "SCHIG17", username = username, password = password)
rm(username, password)

## Save unaltered .gcl's as back-up:
dput(x = SCHIG17.gcl, file = paste("Raw genotypes/SCHIG17_3_Jul0708.gcl.txt", sep = ''))
SCHIG17.gcl$n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2017.3 July 7-8
unique(SCHIG17.gcl$attributes$CAPTURE_DATE)
SCHIG17_3_Jul0708IDs <- AttributesToIDs.GCL(silly = "SCHIG17", attribute = "CAPTURE_DATE", 
                                          matching = unique(SCHIG17.gcl$attributes$CAPTURE_DATE)[4:5])

SCHIG17_3_Jul0708IDs <- list(na.omit(SCHIG17_3_Jul0708IDs))
names(SCHIG17_3_Jul0708IDs) <- "SCHIG17"

PoolCollections.GCL("SCHIG17", loci = loci24, IDs = SCHIG17_3_Jul0708IDs, newname = "SCHIG17_3_Jul0708")
SCHIG17_3_Jul0708.gcl$n ## 190
table(SCHIG17_3_Jul0708.gcl$attributes$CAPTURE_DATE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG17_3_Jul0708_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG17_3_Jul0708", loci24)
min(Original_SCHIG17_3_Jul0708_SampleSizebyLocus) ## Good 181.
apply(Original_SCHIG17_3_Jul0708_SampleSizebyLocus, 1, min) / SCHIG17_3_Jul0708.gcl$n  # 0.95

Original_SCHIG17_3_Jul0708_PercentbyLocus <- apply(Original_SCHIG17_3_Jul0708_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(Original_SCHIG17_3_Jul0708_PercentbyLocus < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_SCHIG17_3_Jul0708_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG17_3_Jul0708_ColSize <- SCHIG17_3_Jul0708.gcl$n

## Remove individuals with >20% missing data
SCHIG17_3_Jul0708_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG17_3_Jul0708", proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG17_3_Jul0708_PostMissLoci <- SCHIG17_3_Jul0708.gcl$n

SCHIG17_3_Jul0708_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, 
                                      dimnames = list("SCHIG17_3_Jul0708", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG17_3_Jul0708_SampleSizes[, "Genotyped"] <- Original_SCHIG17_3_Jul0708_ColSize
SCHIG17_3_Jul0708_SampleSizes[, "Missing"] <- Original_SCHIG17_3_Jul0708_ColSize - ColSize_SCHIG17_3_Jul0708_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG17_3_Jul0708_DuplicateCheck95MinProportion <- 
  CheckDupWithinSilly.GCL(sillyvec = "SCHIG17_3_Jul0708", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG17_3_Jul0708_RemovedDups <- RemoveDups.GCL(SCHIG17_3_Jul0708_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG17_3_Jul0708_PostDuplicate <- SCHIG17_3_Jul0708.gcl$n

SCHIG17_3_Jul0708_SampleSizes[, "Duplicate"] <- ColSize_SCHIG17_3_Jul0708_PostMissLoci-ColSize_SCHIG17_3_Jul0708_PostDuplicate
SCHIG17_3_Jul0708_SampleSizes[, "Final"] <- ColSize_SCHIG17_3_Jul0708_PostDuplicate

SCHIG17_3_Jul0708_SampleSizes
write.xlsx(SCHIG17_3_Jul0708_SampleSizes, file = "Output/SCHIG17_3_Jul0708_SampleSizes.xlsx")
dput(x = SCHIG17_3_Jul0708.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG17_3_Jul0708_IDs.txt")



## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG17_3_Jul0708", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG17_3_Jul0708", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG17_3_Jul0708", 
                loci = loci24[-mito.loci24], 
                path = "Genepop/SCHIG17_3_Jul0708_23nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/SCHIG17_3_Jul0708_23nuclearloci.txt.P")

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

setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG17_3_Jul0708", loci = loci24MSA, IDs = NULL, mixname = "SCHIG17_3_Jul0708",
                                          dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "Objects/ChignikMixtureFormat.txt")


## Rolling prior: 
Chignik2017Period3Prior <- dget(file = "Objects/Chignik2017Period3Prior.txt")


## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG17_3_Jul0708", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2017Period3Prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

## Create output directory
dir.create("BAYES/Output/SCHIG17_3_Jul0708")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG17_3_Jul0708_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG17_3_Jul0708",
                               prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
dput(x = SCHIG17_3_Jul0708_Estimates, file="Estimates objects/SCHIG17_3_Jul0708_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG17_3_Jul0708_Estimates[[1]][[1]][, "GR"]

# View traceplot
par(mfrow = c(2, 1), mar = c(1.1, 4.1, 4.1, 2.1))
plot(SCHIG17_3_Jul0708_Estimates$Output$SCHIG17_3_Jul0708[seq(from = 1, to = 100000, by = 10), 1], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Black Lake", xlab = "")
abline(v = seq(from = 0, to = 10000, by = 2000))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(SCHIG17_3_Jul0708_Estimates$Output$SCHIG17_3_Jul0708[seq(from = 1, to = 100000, by = 10), 2], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Chignik Lake", xlab = "Repetitions", cex.lab = 1.2)
abline(v = seq(from = 0, to = 10000, by = 2000))
mtext(text = "Posterior", side = 2, outer = TRUE, line = -1, cex = 1.2)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG17_3_Jul0708_Estimates$Stats[[1]],
           file="Estimates tables/SCHIG17_3_Jul0708_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="Objects/ChignikInseasonReport.f.txt")

SCHIG17_3_Jul0708_SampleSizes

# July 7, n = 107; July 8, n = 83
day <- as.numeric(as.Date(paste("July", 7), format = "%B %d") - as.Date("05/24", format = "%m/%d")) + (83 / 190)
Day <- "7-8"

ChignikInseasonReport.f(NewData = SCHIG17_3_Jul0708_Estimates, Period = 3, NumSampled = 190, 
                        NumAnalyzed = SCHIG17_3_Jul0708_SampleSizes[1, "Genotyped"],
                        Included = SCHIG17_3_Jul0708_SampleSizes[1, "Final"], Month = "July", Day = 7)


## Prior for next round
Chignik2017Period4Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG17_3_Jul0708_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2017Period4Prior, file = "Objects/Chignik2017Period4Prior.txt")

## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_3.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 4 July 13 2017 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_4.RData")
## load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_4.RData")

#This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

username = "krshedd"

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = username, password = password)

## Save original LocusControl
loci24 <- LocusControl$locusnames
mito.loci24 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl_4_Jul13.txt")
dput(x = loci24, file = "Objects/loci24.txt")
dput(x = mito.loci24, file = "Objects/mito.loci24.txt")

## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "SCHIG17", username = username, password = password)
rm(username, password)

## Save unaltered .gcl's as back-up:
dput(x = SCHIG17.gcl, file = paste("Raw genotypes/SCHIG17_4_Jul13.gcl.txt", sep = ''))
SCHIG17.gcl$n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2017.4 July 13
unique(SCHIG17.gcl$attributes$CAPTURE_DATE)
SCHIG17_4_Jul13IDs <- AttributesToIDs.GCL(silly = "SCHIG17", attribute = "CAPTURE_DATE", 
                                            matching = unique(SCHIG17.gcl$attributes$CAPTURE_DATE)[6])

SCHIG17_4_Jul13IDs <- list(na.omit(SCHIG17_4_Jul13IDs))
names(SCHIG17_4_Jul13IDs) <- "SCHIG17"

PoolCollections.GCL("SCHIG17", loci = loci24, IDs = SCHIG17_4_Jul13IDs, newname = "SCHIG17_4_Jul13")
SCHIG17_4_Jul13.gcl$n ## 190
table(SCHIG17_4_Jul13.gcl$attributes$CAPTURE_DATE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG17_4_Jul13_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG17_4_Jul13", loci24)
min(Original_SCHIG17_4_Jul13_SampleSizebyLocus) ## Good 178.
apply(Original_SCHIG17_4_Jul13_SampleSizebyLocus, 1, min) / SCHIG17_4_Jul13.gcl$n  # 0.94

Original_SCHIG17_4_Jul13_PercentbyLocus <- apply(Original_SCHIG17_4_Jul13_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(Original_SCHIG17_4_Jul13_PercentbyLocus < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_SCHIG17_4_Jul13_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG17_4_Jul13_ColSize <- SCHIG17_4_Jul13.gcl$n

## Remove individuals with >20% missing data
SCHIG17_4_Jul13_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG17_4_Jul13", proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG17_4_Jul13_PostMissLoci <- SCHIG17_4_Jul13.gcl$n

SCHIG17_4_Jul13_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, 
                                        dimnames = list("SCHIG17_4_Jul13", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG17_4_Jul13_SampleSizes[, "Genotyped"] <- Original_SCHIG17_4_Jul13_ColSize
SCHIG17_4_Jul13_SampleSizes[, "Missing"] <- Original_SCHIG17_4_Jul13_ColSize - ColSize_SCHIG17_4_Jul13_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG17_4_Jul13_DuplicateCheck95MinProportion <- 
  CheckDupWithinSilly.GCL(sillyvec = "SCHIG17_4_Jul13", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG17_4_Jul13_RemovedDups <- RemoveDups.GCL(SCHIG17_4_Jul13_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG17_4_Jul13_PostDuplicate <- SCHIG17_4_Jul13.gcl$n

SCHIG17_4_Jul13_SampleSizes[, "Duplicate"] <- ColSize_SCHIG17_4_Jul13_PostMissLoci-ColSize_SCHIG17_4_Jul13_PostDuplicate
SCHIG17_4_Jul13_SampleSizes[, "Final"] <- ColSize_SCHIG17_4_Jul13_PostDuplicate

SCHIG17_4_Jul13_SampleSizes
write.xlsx(SCHIG17_4_Jul13_SampleSizes, file = "Output/SCHIG17_4_Jul13_SampleSizes.xlsx")
dput(x = SCHIG17_4_Jul13.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG17_4_Jul13_IDs.txt")



## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG17_4_Jul13", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG17_4_Jul13", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG17_4_Jul13", 
                loci = loci24[-mito.loci24], 
                path = "Genepop/SCHIG17_4_Jul13_23nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/SCHIG17_4_Jul13_23nuclearloci.txt.P")

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

setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 4 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG17_4_Jul13", loci = loci24MSA, IDs = NULL, mixname = "SCHIG17_4_Jul13",
                                          dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "Objects/ChignikMixtureFormat.txt")


## Rolling prior: 
Chignik2017Period4Prior <- dget(file = "Objects/Chignik2017Period4Prior.txt")


## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG17_4_Jul13", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2017Period4Prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

## Create output directory
dir.create("BAYES/Output/SCHIG17_4_Jul13")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG17_4_Jul13_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG17_4_Jul13",
                               prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
dput(x = SCHIG17_4_Jul13_Estimates, file="Estimates objects/SCHIG17_4_Jul13_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG17_4_Jul13_Estimates[[1]][[1]][, "GR"]

# View traceplot
par(mfrow = c(2, 1), mar = c(1.1, 4.1, 4.1, 2.1))
plot(SCHIG17_4_Jul13_Estimates$Output$SCHIG17_4_Jul13[seq(from = 1, to = 100000, by = 10), 1], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Black Lake", xlab = "")
abline(v = seq(from = 0, to = 10000, by = 2000))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(SCHIG17_4_Jul13_Estimates$Output$SCHIG17_4_Jul13[seq(from = 1, to = 100000, by = 10), 2], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Chignik Lake", xlab = "Repetitions", cex.lab = 1.2)
abline(v = seq(from = 0, to = 10000, by = 2000))
mtext(text = "Posterior", side = 2, outer = TRUE, line = -1, cex = 1.2)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG17_4_Jul13_Estimates$Stats[[1]],
           file="Estimates tables/SCHIG17_4_Jul13_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="Objects/ChignikInseasonReport.f.txt")

SCHIG17_4_Jul13_SampleSizes

# July 7, n = 107; July 8, n = 83
# day <- as.numeric(as.Date(paste("July", 7), format = "%B %d") - as.Date("05/24", format = "%m/%d")) + (83 / 190)
# Day <- "7-8"

ChignikInseasonReport.f(NewData = SCHIG17_4_Jul13_Estimates, Period = 4, NumSampled = 190, 
                        NumAnalyzed = SCHIG17_4_Jul13_SampleSizes[1, "Genotyped"],
                        Included = SCHIG17_4_Jul13_SampleSizes[1, "Final"], Month = "July", Day = 13)


## Prior for next round
Chignik2017Period5Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG17_4_Jul13_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2017Period5Prior, file = "Objects/Chignik2017Period5Prior.txt")

## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_4.RData")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 5 July 18 2017 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_5.RData")
## load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_5.RData")

#This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

username = "krshedd"

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = username, password = password)

## Save original LocusControl
loci24 <- LocusControl$locusnames
mito.loci24 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl_5_Jul18.txt")
dput(x = loci24, file = "Objects/loci24.txt")
dput(x = mito.loci24, file = "Objects/mito.loci24.txt")

## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "SCHIG17", username = username, password = password)
rm(username, password)

## Save unaltered .gcl's as back-up:
dput(x = SCHIG17.gcl, file = paste("Raw genotypes/SCHIG17_5_Jul18.gcl.txt", sep = ''))
SCHIG17.gcl$n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2017.5 July 17
unique(SCHIG17.gcl$attributes$CAPTURE_DATE)
SCHIG17_5_Jul18IDs <- AttributesToIDs.GCL(silly = "SCHIG17", attribute = "CAPTURE_DATE", 
                                          matching = unique(SCHIG17.gcl$attributes$CAPTURE_DATE)[7])

SCHIG17_5_Jul18IDs <- list(na.omit(SCHIG17_5_Jul18IDs))
names(SCHIG17_5_Jul18IDs) <- "SCHIG17"

PoolCollections.GCL("SCHIG17", loci = loci24, IDs = SCHIG17_5_Jul18IDs, newname = "SCHIG17_5_Jul18")
SCHIG17_5_Jul18.gcl$n ## 190
table(SCHIG17_5_Jul18.gcl$attributes$CAPTURE_DATE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG17_5_Jul18_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG17_5_Jul18", loci24)
min(Original_SCHIG17_5_Jul18_SampleSizebyLocus) ## Ya 186.
apply(Original_SCHIG17_5_Jul18_SampleSizebyLocus, 1, min) / SCHIG17_5_Jul18.gcl$n  # 0.98

Original_SCHIG17_5_Jul18_PercentbyLocus <- apply(Original_SCHIG17_5_Jul18_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(Original_SCHIG17_5_Jul18_PercentbyLocus < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_SCHIG17_5_Jul18_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG17_5_Jul18_ColSize <- SCHIG17_5_Jul18.gcl$n

## Remove individuals with >20% missing data
SCHIG17_5_Jul18_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG17_5_Jul18", proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG17_5_Jul18_PostMissLoci <- SCHIG17_5_Jul18.gcl$n

SCHIG17_5_Jul18_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, 
                                      dimnames = list("SCHIG17_5_Jul18", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG17_5_Jul18_SampleSizes[, "Genotyped"] <- Original_SCHIG17_5_Jul18_ColSize
SCHIG17_5_Jul18_SampleSizes[, "Missing"] <- Original_SCHIG17_5_Jul18_ColSize - ColSize_SCHIG17_5_Jul18_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG17_5_Jul18_DuplicateCheck95MinProportion <- 
  CheckDupWithinSilly.GCL(sillyvec = "SCHIG17_5_Jul18", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG17_5_Jul18_RemovedDups <- RemoveDups.GCL(SCHIG17_5_Jul18_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG17_5_Jul18_PostDuplicate <- SCHIG17_5_Jul18.gcl$n

SCHIG17_5_Jul18_SampleSizes[, "Duplicate"] <- ColSize_SCHIG17_5_Jul18_PostMissLoci-ColSize_SCHIG17_5_Jul18_PostDuplicate
SCHIG17_5_Jul18_SampleSizes[, "Final"] <- ColSize_SCHIG17_5_Jul18_PostDuplicate

SCHIG17_5_Jul18_SampleSizes
write.xlsx(SCHIG17_5_Jul18_SampleSizes, file = "Output/SCHIG17_5_Jul18_SampleSizes.xlsx")
dput(x = SCHIG17_5_Jul18.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG17_5_Jul18_IDs.txt")



## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG17_5_Jul18", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG17_5_Jul18", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG17_5_Jul18", 
                loci = loci24[-mito.loci24], 
                path = "Genepop/SCHIG17_5_Jul18_23nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/SCHIG17_5_Jul18_23nuclearloci.txt.P")

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

setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 5 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG17_5_Jul18", loci = loci24MSA, IDs = NULL, mixname = "SCHIG17_5_Jul18",
                                          dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "Objects/ChignikMixtureFormat.txt")


## Rolling prior: 
Chignik2017Period5Prior <- dget(file = "Objects/Chignik2017Period5Prior.txt")


## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG17_5_Jul18", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2017Period5Prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

## Create output directory
dir.create("BAYES/Output/SCHIG17_5_Jul18")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG17_5_Jul18_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG17_5_Jul18",
                               prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
dput(x = SCHIG17_5_Jul18_Estimates, file="Estimates objects/SCHIG17_5_Jul18_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG17_5_Jul18_Estimates[[1]][[1]][, "GR"]

# View traceplot
par(mfrow = c(2, 1), mar = c(1.1, 4.1, 4.1, 2.1))
plot(SCHIG17_5_Jul18_Estimates$Output$SCHIG17_5_Jul18[seq(from = 1, to = 100000, by = 10), 1], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Black Lake", xlab = "")
abline(v = seq(from = 0, to = 10000, by = 2000))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(SCHIG17_5_Jul18_Estimates$Output$SCHIG17_5_Jul18[seq(from = 1, to = 100000, by = 10), 2], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Chignik Lake", xlab = "Repetitions", cex.lab = 1.2)
abline(v = seq(from = 0, to = 10000, by = 2000))
mtext(text = "Posterior", side = 2, outer = TRUE, line = -1, cex = 1.2)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG17_5_Jul18_Estimates$Stats[[1]],
           file="Estimates tables/SCHIG17_5_Jul18_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="Objects/ChignikInseasonReport.f.txt")

SCHIG17_5_Jul18_SampleSizes

# July 7, n = 107; July 8, n = 83
# day <- as.numeric(as.Date(paste("July", 7), format = "%B %d") - as.Date("05/24", format = "%m/%d")) + (83 / 190)
# Day <- "7-8"

ChignikInseasonReport.f(NewData = SCHIG17_5_Jul18_Estimates, Period = 5, NumSampled = 190, 
                        NumAnalyzed = SCHIG17_5_Jul18_SampleSizes[1, "Genotyped"],
                        Included = SCHIG17_5_Jul18_SampleSizes[1, "Final"], Month = "July", Day = 18)


## Prior for next round
Chignik2017Period6Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG17_5_Jul18_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2017Period6Prior, file = "Objects/Chignik2017Period6Prior.txt")

## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_5.RData")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 6 July 18 2017 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_6.RData")
## load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_6.RData")

#This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

username = "krshedd"

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = username, password = password)

## Save original LocusControl
loci24 <- LocusControl$locusnames
mito.loci24 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl_6_Jul23.txt")
dput(x = loci24, file = "Objects/loci24.txt")
dput(x = mito.loci24, file = "Objects/mito.loci24.txt")

## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "SCHIG17", username = username, password = password)
rm(username, password)

## Save unaltered .gcl's as back-up:
dput(x = SCHIG17.gcl, file = paste("Raw genotypes/SCHIG17_6_Jul23.gcl.txt", sep = ''))
SCHIG17.gcl$n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2017.5 July 17
unique(SCHIG17.gcl$attributes$CAPTURE_DATE)
SCHIG17_6_Jul23IDs <- AttributesToIDs.GCL(silly = "SCHIG17", attribute = "CAPTURE_DATE", 
                                          matching = unique(SCHIG17.gcl$attributes$CAPTURE_DATE)[8])

SCHIG17_6_Jul23IDs <- list(na.omit(SCHIG17_6_Jul23IDs))
names(SCHIG17_6_Jul23IDs) <- "SCHIG17"

PoolCollections.GCL("SCHIG17", loci = loci24, IDs = SCHIG17_6_Jul23IDs, newname = "SCHIG17_6_Jul23")
SCHIG17_6_Jul23.gcl$n ## 190
table(SCHIG17_6_Jul23.gcl$attributes$CAPTURE_DATE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG17_6_Jul23_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG17_6_Jul23", loci24)
min(Original_SCHIG17_6_Jul23_SampleSizebyLocus) ## Ya 186.
apply(Original_SCHIG17_6_Jul23_SampleSizebyLocus, 1, min) / SCHIG17_6_Jul23.gcl$n  # 0.98

Original_SCHIG17_6_Jul23_PercentbyLocus <- apply(Original_SCHIG17_6_Jul23_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(Original_SCHIG17_6_Jul23_PercentbyLocus < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_SCHIG17_6_Jul23_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG17_6_Jul23_ColSize <- SCHIG17_6_Jul23.gcl$n

## Remove individuals with >20% missing data
SCHIG17_6_Jul23_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG17_6_Jul23", proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG17_6_Jul23_PostMissLoci <- SCHIG17_6_Jul23.gcl$n

SCHIG17_6_Jul23_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, 
                                      dimnames = list("SCHIG17_6_Jul23", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG17_6_Jul23_SampleSizes[, "Genotyped"] <- Original_SCHIG17_6_Jul23_ColSize
SCHIG17_6_Jul23_SampleSizes[, "Missing"] <- Original_SCHIG17_6_Jul23_ColSize - ColSize_SCHIG17_6_Jul23_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG17_6_Jul23_DuplicateCheck95MinProportion <- 
  CheckDupWithinSilly.GCL(sillyvec = "SCHIG17_6_Jul23", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG17_6_Jul23_RemovedDups <- RemoveDups.GCL(SCHIG17_6_Jul23_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG17_6_Jul23_PostDuplicate <- SCHIG17_6_Jul23.gcl$n

SCHIG17_6_Jul23_SampleSizes[, "Duplicate"] <- ColSize_SCHIG17_6_Jul23_PostMissLoci-ColSize_SCHIG17_6_Jul23_PostDuplicate
SCHIG17_6_Jul23_SampleSizes[, "Final"] <- ColSize_SCHIG17_6_Jul23_PostDuplicate

SCHIG17_6_Jul23_SampleSizes
write.xlsx(SCHIG17_6_Jul23_SampleSizes, file = "Output/SCHIG17_6_Jul23_SampleSizes.xlsx")
dput(x = SCHIG17_6_Jul23.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG17_6_Jul23_IDs.txt")



## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG17_6_Jul23", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG17_6_Jul23", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG17_6_Jul23", 
                loci = loci24[-mito.loci24], 
                path = "Genepop/SCHIG17_6_Jul23_23nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/SCHIG17_6_Jul23_23nuclearloci.txt.P")

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

setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 6 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG17_6_Jul23", loci = loci24MSA, IDs = NULL, mixname = "SCHIG17_6_Jul23",
                                          dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "Objects/ChignikMixtureFormat.txt")


## Rolling prior: 
Chignik2017Period6Prior <- dget(file = "Objects/Chignik2017Period6Prior.txt")


## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG17_6_Jul23", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2017Period6Prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

## Create output directory
dir.create("BAYES/Output/SCHIG17_6_Jul23")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG17_6_Jul23_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG17_6_Jul23",
                               prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
dput(x = SCHIG17_6_Jul23_Estimates, file="Estimates objects/SCHIG17_6_Jul23_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG17_6_Jul23_Estimates[[1]][[1]][, "GR"]

# View traceplot
par(mfrow = c(2, 1), mar = c(1.1, 4.1, 4.1, 2.1))
plot(SCHIG17_6_Jul23_Estimates$Output$SCHIG17_6_Jul23[seq(from = 1, to = 100000, by = 10), 1], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Black Lake", xlab = "")
abline(v = seq(from = 0, to = 10000, by = 2000))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(SCHIG17_6_Jul23_Estimates$Output$SCHIG17_6_Jul23[seq(from = 1, to = 100000, by = 10), 2], 
     type = "l", ylim = c(0, 1), ylab = "", main = "Chignik Lake", xlab = "Repetitions", cex.lab = 1.2)
abline(v = seq(from = 0, to = 10000, by = 2000))
mtext(text = "Posterior", side = 2, outer = TRUE, line = -1, cex = 1.2)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG17_6_Jul23_Estimates$Stats[[1]],
           file="Estimates tables/SCHIG17_6_Jul23_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="Objects/ChignikInseasonReport.f.txt")

SCHIG17_6_Jul23_SampleSizes

# July 7, n = 107; July 8, n = 83
# day <- as.numeric(as.Date(paste("July", 7), format = "%B %d") - as.Date("05/24", format = "%m/%d")) + (83 / 190)
# Day <- "7-8"

ChignikInseasonReport.f(NewData = SCHIG17_6_Jul23_Estimates, Period = 6, NumSampled = 190, 
                        NumAnalyzed = SCHIG17_6_Jul23_SampleSizes[1, "Genotyped"],
                        Included = SCHIG17_6_Jul23_SampleSizes[1, "Final"], Month = "July", Day = 23)


## Prior for next round
Chignik2017Period7Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG17_6_Jul23_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2017Period7Prior, file = "Objects/Chignik2017Period7Prior.txt")

## save.image("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/2017ChignikInseason_6.RData")