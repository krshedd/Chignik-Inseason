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
SCHIG17_1_Jun25IDs <- AttributesToIDs.GCL(silly = "SCHIG17", attribute = "CAPTURE_DATE", matching = unique(SCHIG17.gcl$attributes$CAPTURE_DATE)[1])

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
min(Original_SCHIG17_1_Jun25_SampleSizebyLocus) ## Good? Yes, 185.
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


