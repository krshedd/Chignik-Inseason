## 2015 Chignik Escapement mixture analysis - Inseason!
## Kyle Shedd Mon Jun 13 09:43:05 2016

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 June 27 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run by THD on Mon Jun 29 2015

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason")
## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_1.RData")
## load("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_1.RData")

#This sources all of the new GCL functions to this workspace
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")

## ReadLOKI.GCL does not work in RStudio, have to go into base R
## Pull all data for each silly code and create .gcl objects for each
ReadLOKI.GCL(sillyvec = "SCHIG15", markersuite = "Sockeye2013Chignik_24SNPs")
objects(pattern = "\\.gcl") # check objects
str(SCHIG15.gcl)

## Create a locus list from Locus Control object
loci24 <- LocusControl$locusnames

## Save unaltered .gcl's as back-up:
dput(x = SCHIG15.gcl, file = paste("Raw genotypes/SCHIG15_1_Jun27.gcl.txt", sep = ''))

## Save original LocusControl
dput(x = LocusControl, file = "Objects/OriginalLocusControl_1_Jun27.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2015.1 June 27
unique(SCHIG15.gcl$attributes$CAPTURE_DATE)
SCHIG15_1_Jun27IDs <- AttributesToIDs.GCL(silly = "SCHIG15", attribute = "CAPTURE_DATE", matching = unique(SCHIG15.gcl$attributes$CAPTURE_DATE)[1])

SCHIG15_1_Jun27IDs <- list(na.omit(SCHIG15_1_Jun27IDs))
names(SCHIG15_1_Jun27IDs) <- "SCHIG15"

PoolCollections.GCL("SCHIG15", loci = loci24, IDs = SCHIG15_1_Jun27IDs, newname = "SCHIG15_1_Jun27")
SCHIG15_1_Jun27.gcl$n ## 190


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG15_1_Jun27_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG15_1_Jun27", loci24)
min(Original_SCHIG15_1_Jun27_SampleSizebyLocus) ## Good? Yes, 185.
apply(Original_SCHIG15_1_Jun27_SampleSizebyLocus, 1, min) / apply(Original_SCHIG15_1_Jun27_SampleSizebyLocus, 1, max)

## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG15_1_Jun27_ColSize <- SCHIG15_1_Jun27.gcl$n

## Remove individuals with >20% missing data
SCHIG15_1_Jun27_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG15_1_Jun27", loci = loci24, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG15_1_Jun27_PostMissLoci <- SCHIG15_1_Jun27.gcl$n

SCHIG15_1_Jun27_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, dimnames = list("SCHIG15_1_Jun27", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG15_1_Jun27_SampleSizes[, "Genotyped"] <- Original_SCHIG15_1_Jun27_ColSize
SCHIG15_1_Jun27_SampleSizes[, "Missing"] <- Original_SCHIG15_1_Jun27_ColSize - ColSize_SCHIG15_1_Jun27_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG15_1_Jun27_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SCHIG15_1_Jun27", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG15_1_Jun27_RemovedDups <- RemoveDups.GCL(SCHIG15_1_Jun27_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG15_1_Jun27_PostDuplicate <- SCHIG15_1_Jun27.gcl$n

SCHIG15_1_Jun27_SampleSizes[, "Duplicate"] <- ColSize_SCHIG15_1_Jun27_PostMissLoci-ColSize_SCHIG15_1_Jun27_PostDuplicate
SCHIG15_1_Jun27_SampleSizes[, "Final"] <- ColSize_SCHIG15_1_Jun27_PostDuplicate

write.xlsx(SCHIG15_1_Jun27_SampleSizes, file = "Output/SCHIG15_1_Jun27_SampleSizes.xlsx")
dput(x = SCHIG15_1_Jun27.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG15_1_Jun27_IDs.txt")

## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG15_1_Jun27", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG15_1_Jun27", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG15_1_Jun27", loci = loci24[-23], path = "Genepop/SCHIG15_1_Jun27_23nuclearloci.gen", VialNums = TRUE)
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")

# Read in Genepop output .P file
HWE <- ReadGenepopHWEKS.GCL(file = "Genepop/SCHIG15_1_Jun27_23nuclearloci.txt.P")

# Plot Fis values
plot(sort(HWE[, "WC Fis"]), type = "h", lwd = 10, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40"); abline(h = 0, lwd = 5)

# Look at data for any markers out of HWE
HWE[HWE$PValue < 0.05, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get baseline objects needed for MSA
Chignik7Populations <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Groupvec7.txt")
Chignik24BaselineFormat <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik24BaselineFormat.txt")
Inits <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Inits.txt")
loci24MSA <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/loci24.txt")

## Defining the random seeds as the same as WASSIP mixtures for repeatability.
WASSIPSockeyeSeeds <- dget(file = "V:/WORK/WASSIP/Sockeye/Mixture/Objects/WASSIPSockeyeSeeds.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG15_1_Jun27", loci = loci24MSA, IDs = NULL, mixname = "SCHIG15_1_Jun27",
                                          dir = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/BAYES/Mixture", type = "BAYES", PT = FALSE) ## Changed to 2015 directory
dput(ChignikMixtureFormat, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikMixtureFormat.txt") ## Changed to 2015 directory

ChignikGroups <- dget(file = "V:/WORK/Sockeye/Chignik/2013 Chignik Inseason/Objects/ChignikGroups.txt")


## Start Priors are based upon best available information, Birch's logistic fit of the past five years for the 28th: 
# Data from BirchWLS.Estimate values, see "V:\DOC\Power point presentations\Sockeye\Chignik\2015 CRAA\Chignik inseason summary data.R" line 583
Chignik2015StartPriorWeights <- c(0.899, 0.101)
Chignik2015StartPrior <- Prior.GCL(groupvec = Groupvec7, groupweights = Chignik2015StartPriorWeights, minval = 0.01)
dput(x = Chignik2015StartPrior, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015StartPrior.txt") ## Added .txt extension, capitalized Objects

## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG15_1_Jun27", basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2015StartPrior, initmat = Inits, dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,1),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 
## I almost thinned your individual assignments, because the files are ridiculously large, but I thought maybe you didn't for a reason.  I recommend thinning by 100 for future runs. THD 062915

## Create output directory
dir.create("BAYES/Output/SCHIG15_1_Jun27")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG15_1_Jun27_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG15_1_Jun27",
                                                          prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE) ## Same here, will you use the posterior?  They really blow up memory usage.
dput(x = SCHIG15_1_Jun27_Estimates, file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates objects/SCHIG15_1_Jun27_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG15_1_Jun27_Estimates[[1]][[1]][, "GR"]

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG15_1_Jun27_Estimates$Stats[[1]],
           file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates tables/SCHIG15_1_Jun27_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReport.f.txt")
# Old function if the new one doesn't work
# ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReportOLD.f.txt")


SCHIG15_1_Jun27_SampleSizes

ChignikInseasonReport.f(NewData = SCHIG15_1_Jun27_Estimates, Period = 1, NumSampled = 190, NumAnalyzed = SCHIG15_1_Jun27_SampleSizes[1, "Genotyped"],
                        Included = SCHIG15_1_Jun27_SampleSizes[1, "Final"], Month = "June", Day = 27)

## Use Adobe Acrobat Pro on Gen1 to add the footer of page, date, and time
## Create News Release Fishery Update based on the previous round, change the numbers manually in the text
## Post Fishery Update online using the News Release System

## Prior for next round
Chignik2015Period2Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG15_1_Jun27_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2015Period2Prior, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015Period2Prior.txt")

## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_1.RData")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 July 1 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run by KRS on Fri Jul 03 10:00:00 2015
date()

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason")
## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_2.RData")
## load("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_2.RData")

#This sources all of the new GCL functions to this workspace
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")

## ReadLOKI.GCL does not work in RStudio, have to go into base R
## Pull all data for each silly code and create .gcl objects for each
ReadLOKI.GCL(sillyvec = "SCHIG15", markersuite = "Sockeye2013Chignik_24SNPs")
objects(pattern = "\\.gcl") # check objects
str(SCHIG15.gcl)

## Create a locus list from Locus Control object
loci24 <- LocusControl$locusnames

## Save unaltered .gcl's as back-up:
dput(x = SCHIG15.gcl, file = paste("Raw genotypes/SCHIG15_2_Jul01.gcl.txt", sep = ''))

## Save original LocusControl
dput(x = LocusControl, file = "Objects/OriginalLocusControl_2_Jul01.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2015.2 July 1
unique(SCHIG15.gcl$attributes$CAPTURE_DATE)
SCHIG15_2_Jul01IDs <- AttributesToIDs.GCL(silly = "SCHIG15", attribute = "CAPTURE_DATE", matching = unique(SCHIG15.gcl$attributes$CAPTURE_DATE)[2])

SCHIG15_2_Jul01IDs <- list(na.omit(SCHIG15_2_Jul01IDs))
names(SCHIG15_2_Jul01IDs) <- "SCHIG15"

PoolCollections.GCL("SCHIG15", loci = loci24, IDs = SCHIG15_2_Jul01IDs, newname = "SCHIG15_2_Jul01")
SCHIG15_2_Jul01.gcl$n ## 189, one vial was dropped in river and thus not sampled


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG15_2_Jul01_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG15_2_Jul01", loci24)
min(Original_SCHIG15_2_Jul01_SampleSizebyLocus) ## Good? Yes, 184.
apply(Original_SCHIG15_2_Jul01_SampleSizebyLocus, 1, min) / apply(Original_SCHIG15_2_Jul01_SampleSizebyLocus, 1, max) # 0.978

## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG15_2_Jul01_ColSize <- SCHIG15_2_Jul01.gcl$n

## Remove individuals with >20% missing data
SCHIG15_2_Jul01_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG15_2_Jul01", loci = loci24, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG15_2_Jul01_PostMissLoci <- SCHIG15_2_Jul01.gcl$n

SCHIG15_2_Jul01_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, dimnames = list("SCHIG15_2_Jul01", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG15_2_Jul01_SampleSizes[, "Genotyped"] <- Original_SCHIG15_2_Jul01_ColSize
SCHIG15_2_Jul01_SampleSizes[, "Missing"] <- Original_SCHIG15_2_Jul01_ColSize - ColSize_SCHIG15_2_Jul01_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG15_2_Jul01_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SCHIG15_2_Jul01", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG15_2_Jul01_RemovedDups <- RemoveDups.GCL(SCHIG15_2_Jul01_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG15_2_Jul01_PostDuplicate <- SCHIG15_2_Jul01.gcl$n

SCHIG15_2_Jul01_SampleSizes[, "Duplicate"] <- ColSize_SCHIG15_2_Jul01_PostMissLoci-ColSize_SCHIG15_2_Jul01_PostDuplicate
SCHIG15_2_Jul01_SampleSizes[, "Final"] <- ColSize_SCHIG15_2_Jul01_PostDuplicate

write.xlsx(SCHIG15_2_Jul01_SampleSizes, file = "Output/SCHIG15_2_Jul01_SampleSizes.xlsx")
dput(x = SCHIG15_1_Jun27.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG15_1_Jun27_IDs.txt")

## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG15_2_Jul01", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG15_2_Jul01", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG15_2_Jul01", loci = loci24[-23], path = "Genepop/SCHIG15_2_Jul01_23nuclearloci.gen", VialNums = TRUE)
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")

# Read in Genepop output .P file
HWE <- ReadGenepopHWEKS.GCL(file = "Genepop/SCHIG15_2_Jul01_23nuclearloci.txt.P")

# Plot Fis values
plot(sort(HWE[, "WC Fis"]), type = "h", lwd = 10, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40"); abline(h = 0, lwd = 5)

# Look at data for any markers out of HWE
HWE[HWE$PValue < 0.05, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get baseline objects needed for MSA
Chignik7Populations <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Groupvec7.txt")
Chignik24BaselineFormat <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik24BaselineFormat.txt")
Inits <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Inits.txt")
loci24MSA <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/loci24.txt")

## Defining the random seeds as the same as WASSIP mixtures for repeatability.
WASSIPSockeyeSeeds <- dget(file = "V:/WORK/WASSIP/Sockeye/Mixture/Objects/WASSIPSockeyeSeeds.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG15_2_Jul01", loci = loci24MSA, IDs = NULL, mixname = "SCHIG15_2_Jul01",
                                          dir = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikMixtureFormat.txt")

ChignikGroups <- dget(file = "V:/WORK/Sockeye/Chignik/2013 Chignik Inseason/Objects/ChignikGroups.txt")

## Get prior from previous round (sequential priors)
Chignik2015Period2Prior <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015Period2Prior.txt")

## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG15_2_Jul01", basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2015Period2Prior, initmat = Inits, dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") ## Switch is on for outputing CLS file with thinning by 100

## Create output directory
dir.create("BAYES/Output/SCHIG15_2_Jul01")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG15_2_Jul01_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG15_2_Jul01",
                                                          prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
dput(x = SCHIG15_2_Jul01_Estimates, file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates objects/SCHIG15_2_Jul01_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG15_2_Jul01_Estimates$Stats[[1]][, "GR"]
plot(SCHIG15_2_Jul01_Estimates[[2]][[1]][, 1], type = "l", ylim = c(0, 1), ylab = "Proportion Black Lake", xlab = "Interation")
abline(v = seq(from = 0, to = 1e5, by = 2e4), lwd = 2)
text(x = seq(from = 1e4, to = 9e4, by = 2e4), y = 1.1, xpd = TRUE, labels = paste("Chain", 1:5))

# Sneak Peak
SCHIG15_2_Jul01_Estimates$Stats[[1]]
hist(SCHIG15_2_Jul01_Estimates[[2]][[1]][, 1], col = 8, main = "Proportion Black Lake Posterior", xlab = "Proportion")
abline(v = mean(SCHIG15_2_Jul01_Estimates[[2]][[1]][, 1]), lwd = 3)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG15_2_Jul01_Estimates$Stats[[1]],
           file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates tables/SCHIG15_2_Jul01_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReport.f.txt")
# Old function if the new one doesn't work
# ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReportOLD.f.txt")


SCHIG15_2_Jul01_SampleSizes

ChignikInseasonReport.f(NewData = SCHIG15_2_Jul01_Estimates, Period = 2, NumSampled = 189, NumAnalyzed = SCHIG15_2_Jul01_SampleSizes[1, "Genotyped"],
                        Included = SCHIG15_2_Jul01_SampleSizes[1, "Final"], Month = "July", Day = 1)

## Use Adobe Acrobat Pro on Gen1 to add the footer of page, date, and time
## Create News Release Fishery Update based on the previous round, change the numbers manually in the text
## Post Fishery Update online using the News Release System

## Update 2015 Chignik Estimates Summary and send to Birch/Mary-Beth

## Prior for next round
Chignik2015Period3Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG15_2_Jul01_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2015Period3Prior, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015Period3Prior.txt")

## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_2.RData")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 July 5 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run by KRS on Tue Jul 07 15:04:07 2015
date()

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason")
## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_3.RData")
## load("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_3.RData")

#This sources all of the new GCL functions to this workspace
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")

## ReadLOKI.GCL does not work in RStudio, have to go into base R
## Pull all data for each silly code and create .gcl objects for each
ReadLOKI.GCL(sillyvec = "SCHIG15", markersuite = "Sockeye2013Chignik_24SNPs")
objects(pattern = "\\.gcl") # check objects
str(SCHIG15.gcl)

## Create a locus list from Locus Control object
loci24 <- LocusControl$locusnames

## Save unaltered .gcl's as back-up:
dput(x = SCHIG15.gcl, file = paste("Raw genotypes/SCHIG15_3_Jul05.gcl.txt", sep = ''))

## Save original LocusControl
dput(x = LocusControl, file = "Objects/OriginalLocusControl_3_Jul05.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2015.3 July 5
unique(SCHIG15.gcl$attributes$CAPTURE_DATE)
SCHIG15_3_Jul05IDs <- AttributesToIDs.GCL(silly = "SCHIG15", attribute = "CAPTURE_DATE", matching = unique(SCHIG15.gcl$attributes$CAPTURE_DATE)[3])

SCHIG15_3_Jul05IDs <- list(na.omit(SCHIG15_3_Jul05IDs))
names(SCHIG15_3_Jul05IDs) <- "SCHIG15"

PoolCollections.GCL("SCHIG15", loci = loci24, IDs = SCHIG15_3_Jul05IDs, newname = "SCHIG15_3_Jul05")
SCHIG15_3_Jul05.gcl$n ## 189, double axillary in one vial


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG15_3_Jul05_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG15_3_Jul05", loci24)
min(Original_SCHIG15_3_Jul05_SampleSizebyLocus) ## Good? Yes, 175.
apply(Original_SCHIG15_3_Jul05_SampleSizebyLocus, 1, min) / apply(Original_SCHIG15_3_Jul05_SampleSizebyLocus, 1, max) # 0.926

## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG15_3_Jul05_ColSize <- SCHIG15_3_Jul05.gcl$n

## Remove individuals with >20% missing data
SCHIG15_3_Jul05_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG15_3_Jul05", loci = loci24, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG15_3_Jul05_PostMissLoci <- SCHIG15_3_Jul05.gcl$n

SCHIG15_3_Jul05_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, dimnames = list("SCHIG15_3_Jul05", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG15_3_Jul05_SampleSizes[, "Genotyped"] <- Original_SCHIG15_3_Jul05_ColSize
SCHIG15_3_Jul05_SampleSizes[, "Missing"] <- Original_SCHIG15_3_Jul05_ColSize - ColSize_SCHIG15_3_Jul05_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG15_3_Jul05_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SCHIG15_3_Jul05", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG15_3_Jul05_RemovedDups <- RemoveDups.GCL(SCHIG15_3_Jul05_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG15_3_Jul05_PostDuplicate <- SCHIG15_3_Jul05.gcl$n

SCHIG15_3_Jul05_SampleSizes[, "Duplicate"] <- ColSize_SCHIG15_3_Jul05_PostMissLoci-ColSize_SCHIG15_3_Jul05_PostDuplicate
SCHIG15_3_Jul05_SampleSizes[, "Final"] <- ColSize_SCHIG15_3_Jul05_PostDuplicate

write.xlsx(SCHIG15_3_Jul05_SampleSizes, file = "Output/SCHIG15_3_Jul05_SampleSizes.xlsx")
dput(x = SCHIG15_1_Jun27.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG15_1_Jun27_IDs.txt")

## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG15_3_Jul05", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG15_3_Jul05", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG15_3_Jul05", loci = loci24[-23], path = "Genepop/SCHIG15_3_Jul05_23nuclearloci.gen", VialNums = TRUE)
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")

# Read in Genepop output .P file
HWE <- ReadGenepopHWEKS.GCL(file = "Genepop/SCHIG15_3_Jul05_23nuclearloci.txt.P")

# Plot Fis values
plot(sort(HWE[, "WC Fis"]), type = "h", lwd = 10, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40"); abline(h = 0, lwd = 5)

# Look at data for any markers out of HWE
HWE[HWE$PValue < 0.05, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get baseline objects needed for MSA
Chignik7Populations <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Groupvec7.txt")
Chignik24BaselineFormat <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik24BaselineFormat.txt")
Inits <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Inits.txt")
loci24MSA <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/loci24.txt")

## Defining the random seeds as the same as WASSIP mixtures for repeatability.
WASSIPSockeyeSeeds <- dget(file = "V:/WORK/WASSIP/Sockeye/Mixture/Objects/WASSIPSockeyeSeeds.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG15_3_Jul05", loci = loci24MSA, IDs = NULL, mixname = "SCHIG15_3_Jul05",
                                          dir = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikMixtureFormat.txt")

ChignikGroups <- dget(file = "V:/WORK/Sockeye/Chignik/2013 Chignik Inseason/Objects/ChignikGroups.txt")

## Get prior from previous round (sequential priors)
Chignik2015Period3Prior <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015Period3Prior.txt")

## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG15_3_Jul05", basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2015Period3Prior, initmat = Inits, dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") ## Switch is on for outputing CLS file with thinning by 100

## Create output directory
dir.create("BAYES/Output/SCHIG15_3_Jul05")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG15_3_Jul05_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG15_3_Jul05",
                                                          prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
dput(x = SCHIG15_3_Jul05_Estimates, file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates objects/SCHIG15_3_Jul05_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG15_3_Jul05_Estimates$Stats[[1]][, "GR"]

# View trace of 2nd half of all 5 chains
plot(SCHIG15_3_Jul05_Estimates[[2]][[1]][, 1], type = "l", ylim = c(0, 1), ylab = "Proportion Black Lake", xlab = "Interation")
abline(v = seq(from = 0, to = 1e5, by = 2e4), lwd = 2)
text(x = seq(from = 1e4, to = 9e4, by = 2e4), y = 1.1, xpd = TRUE, labels = paste("Chain", 1:5))

# Sneak Peak
SCHIG15_3_Jul05_Estimates$Stats[[1]]
hist(SCHIG15_3_Jul05_Estimates[[2]][[1]][, 1], col = 8, main = "Proportion Black Lake Posterior", xlab = "Proportion")
abline(v = mean(SCHIG15_3_Jul05_Estimates[[2]][[1]][, 1]), lwd = 3)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG15_3_Jul05_Estimates$Stats[[1]],
           file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates tables/SCHIG15_3_Jul05_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReport.f.txt")
# Old function if the new one doesn't work
# ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReportOLD.f.txt")

SCHIG15_3_Jul05_SampleSizes

ChignikInseasonReport.f(NewData = SCHIG15_3_Jul05_Estimates, Period = 3, NumSampled = 190, NumAnalyzed = SCHIG15_3_Jul05_SampleSizes[1, "Genotyped"],
                        Included = SCHIG15_3_Jul05_SampleSizes[1, "Final"], Month = "July", Day = 5)

## Use Adobe Acrobat Pro on Gen1 to add the footer of page, date, and time
## Create News Release Fishery Update based on the previous round, change the numbers manually in the text
## Post Fishery Update online using the News Release System

## Update 2015 Chignik Estimates Summary and send to Birch/Mary-Beth

## Prior for next round
Chignik2015Period4Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG15_3_Jul05_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2015Period4Prior, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015Period4Prior.txt")

## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_3.RData")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 4 July 12 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run by KRS on Tue Jul 14 14:12:44 2015
date()

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason")
## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_4.RData")
## load("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_4.RData")

#This sources all of the new GCL functions to this workspace
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")

## ReadLOKI.GCL does not work in RStudio, have to go into base R
## Pull all data for each silly code and create .gcl objects for each
ReadLOKI.GCL(sillyvec = "SCHIG15", markersuite = "Sockeye2013Chignik_24SNPs")
objects(pattern = "\\.gcl") # check objects
str(SCHIG15.gcl)

## Create a locus list from Locus Control object
loci24 <- LocusControl$locusnames

## Save unaltered .gcl's as back-up:
dput(x = SCHIG15.gcl, file = paste("Raw genotypes/SCHIG15_4_Jul12.gcl.txt", sep = ''))

## Save original LocusControl
dput(x = LocusControl, file = "Objects/OriginalLocusControl_4_Jul12.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2015.4 July 9
unique(SCHIG15.gcl$attributes$CAPTURE_DATE)
SCHIG15_4_Jul12IDs <- AttributesToIDs.GCL(silly = "SCHIG15", attribute = "CAPTURE_DATE", matching = unique(SCHIG15.gcl$attributes$CAPTURE_DATE)[4])

SCHIG15_4_Jul12IDs <- list(na.omit(SCHIG15_4_Jul12IDs))
names(SCHIG15_4_Jul12IDs) <- "SCHIG15"

PoolCollections.GCL("SCHIG15", loci = loci24, IDs = SCHIG15_4_Jul12IDs, newname = "SCHIG15_4_Jul12")
SCHIG15_4_Jul12.gcl$n ## 190


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG15_4_Jul12_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG15_4_Jul12", loci24)
min(Original_SCHIG15_4_Jul12_SampleSizebyLocus) ## Good? Yes, 182.
apply(Original_SCHIG15_4_Jul12_SampleSizebyLocus, 1, min) / apply(Original_SCHIG15_4_Jul12_SampleSizebyLocus, 1, max) # 0.958

## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG15_4_Jul12_ColSize <- SCHIG15_4_Jul12.gcl$n

## Remove individuals with >20% missing data
SCHIG15_4_Jul12_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG15_4_Jul12", loci = loci24, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG15_4_Jul12_PostMissLoci <- SCHIG15_4_Jul12.gcl$n

SCHIG15_4_Jul12_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, dimnames = list("SCHIG15_4_Jul12", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG15_4_Jul12_SampleSizes[, "Genotyped"] <- Original_SCHIG15_4_Jul12_ColSize
SCHIG15_4_Jul12_SampleSizes[, "Missing"] <- Original_SCHIG15_4_Jul12_ColSize - ColSize_SCHIG15_4_Jul12_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG15_4_Jul12_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SCHIG15_4_Jul12", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG15_4_Jul12_RemovedDups <- RemoveDups.GCL(SCHIG15_4_Jul12_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG15_4_Jul12_PostDuplicate <- SCHIG15_4_Jul12.gcl$n

SCHIG15_4_Jul12_SampleSizes[, "Duplicate"] <- ColSize_SCHIG15_4_Jul12_PostMissLoci-ColSize_SCHIG15_4_Jul12_PostDuplicate
SCHIG15_4_Jul12_SampleSizes[, "Final"] <- ColSize_SCHIG15_4_Jul12_PostDuplicate

SCHIG15_4_Jul12_SampleSizes

write.xlsx(SCHIG15_4_Jul12_SampleSizes, file = "Output/SCHIG15_4_Jul12_SampleSizes.xlsx")
dput(x = SCHIG15_1_Jun27.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG15_1_Jun27_IDs.txt")

## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG15_4_Jul12", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG15_4_Jul12", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG15_4_Jul12", loci = loci24[-23], path = "Genepop/SCHIG15_4_Jul12_23nuclearloci.gen", VialNums = TRUE)
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")

# Read in Genepop output .P file
HWE <- ReadGenepopHWEKS.GCL(file = "Genepop/SCHIG15_4_Jul12_23nuclearloci.txt.P")

# Plot Fis values
plot(sort(HWE[, "WC Fis"]), type = "h", lwd = 10, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40"); abline(h = 0, lwd = 5)

# Look at data for any markers out of HWE
HWE[HWE$PValue < 0.05, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get baseline objects needed for MSA
Chignik7Populations <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Groupvec7.txt")
Chignik24BaselineFormat <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik24BaselineFormat.txt")
Inits <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Inits.txt")
loci24MSA <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/loci24.txt")

## Defining the random seeds as the same as WASSIP mixtures for repeatability.
WASSIPSockeyeSeeds <- dget(file = "V:/WORK/WASSIP/Sockeye/Mixture/Objects/WASSIPSockeyeSeeds.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 4 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG15_4_Jul12", loci = loci24MSA, IDs = NULL, mixname = "SCHIG15_4_Jul12",
                                          dir = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikMixtureFormat.txt")

ChignikGroups <- dget(file = "V:/WORK/Sockeye/Chignik/2013 Chignik Inseason/Objects/ChignikGroups.txt")

## Get prior from previous round (sequential priors)
Chignik2015Period4Prior <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015Period4Prior.txt")

## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG15_4_Jul12", basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2015Period4Prior, initmat = Inits, dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") ## Switch is on for outputing CLS file with thinning by 100

## Create output directory
dir.create("BAYES/Output/SCHIG15_4_Jul12")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG15_4_Jul12_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG15_4_Jul12",
                                                          prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
dput(x = SCHIG15_4_Jul12_Estimates, file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates objects/SCHIG15_4_Jul12_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG15_4_Jul12_Estimates$Stats[[1]][, "GR"]

# View trace of 2nd half of all 5 chains
plot(SCHIG15_4_Jul12_Estimates[[2]][[1]][, 1], type = "l", ylim = c(0, 1), ylab = "Proportion Black Lake", xlab = "Interation")
abline(v = seq(from = 0, to = 1e5, by = 2e4), lwd = 2)
text(x = seq(from = 1e4, to = 9e4, by = 2e4), y = 1.1, xpd = TRUE, labels = paste("Chain", 1:5))

# Sneak Peak
SCHIG15_4_Jul12_Estimates$Stats[[1]]
hist(SCHIG15_4_Jul12_Estimates[[2]][[1]][, 1], col = 8, main = "Proportion Black Lake Posterior", xlab = "Proportion")
abline(v = mean(SCHIG15_4_Jul12_Estimates[[2]][[1]][, 1]), lwd = 3)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG15_4_Jul12_Estimates$Stats[[1]],
           file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates tables/SCHIG15_4_Jul12_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReport.f.txt")
# Old function if the new one doesn't work
# ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReportOLD.f.txt")

SCHIG15_4_Jul12_SampleSizes

ChignikInseasonReport.f(NewData = SCHIG15_4_Jul12_Estimates, Period = 4, NumSampled = 190, NumAnalyzed = SCHIG15_4_Jul12_SampleSizes[1, "Genotyped"],
                        Included = SCHIG15_4_Jul12_SampleSizes[1, "Final"], Month = "July", Day = 12)

# Confirm the right fish!!!
SCHIG15_4_Jul12.gcl$attributes$CAPTURE_DATE
SCHIG15_4_Jul12.gcl$attributes$FK_FISH_ID

## Use Adobe Acrobat Pro on Gen1 to add the footer of page, date, and time
## Create News Release Fishery Update based on the previous round, change the numbers manually in the text
## Post Fishery Update online using the News Release System

## Update 2015 Chignik Estimates Summary and send to Birch/Mary-Beth

## Prior for next round
Chignik2015Period5Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG15_4_Jul12_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2015Period5Prior, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015Period5Prior.txt")

## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_4.RData")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 5 July 18 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run by KRS on Mon Jul 20 15:16:37 2015
date()

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason")
## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_5.RData")
## load("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_5.RData")

#This sources all of the new GCL functions to this workspace
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")
source("V:/DATA/R_GEN/Edited GCL Source Scripts/ReadLOKIv2.GCL.r")

## ReadLOKI.GCL does not work in RStudio, have to go into base R
## Pull all data for each silly code and create .gcl objects for each
#ReadLOKI.GCL(sillyvec = "SCHIG15", markersuite = "Sockeye2013Chignik_24SNPs")
ReadLOKIv2.GCL(sillyvec = "SCHIG15", markersuite = "Sockeye2013Chignik_24SNPs")

objects(pattern = "\\.gcl") # check objects
str(SCHIG15.gcl)

## Create a locus list from Locus Control object
loci24 <- LocusControl$locusnames

## Save unaltered .gcl's as back-up:
dput(x = SCHIG15.gcl, file = paste("Raw genotypes/SCHIG15_5_Jul18.gcl.txt", sep = ''))

## Save original LocusControl
dput(x = LocusControl, file = "Objects/OriginalLocusControl_5_Jul18.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2015.5 July 19
unique(SCHIG15.gcl$attributes$CAPTURE_DATE)
SCHIG15_5_Jul18IDs <- AttributesToIDs.GCL(silly = "SCHIG15", attribute = "CAPTURE_DATE", matching = unique(SCHIG15.gcl$attributes$CAPTURE_DATE)[5])

SCHIG15_5_Jul18IDs <- list(na.omit(SCHIG15_5_Jul18IDs))
names(SCHIG15_5_Jul18IDs) <- "SCHIG15"

PoolCollections.GCL("SCHIG15", loci = loci24, IDs = SCHIG15_5_Jul18IDs, newname = "SCHIG15_5_Jul18")
SCHIG15_5_Jul18.gcl$n ## 190


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG15_5_Jul18_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG15_5_Jul18", loci24)
min(Original_SCHIG15_5_Jul18_SampleSizebyLocus) ## Good? Yes, 171.
apply(Original_SCHIG15_5_Jul18_SampleSizebyLocus, 1, min) / apply(Original_SCHIG15_5_Jul18_SampleSizebyLocus, 1, max) # 0.958

## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG15_5_Jul18_ColSize <- SCHIG15_5_Jul18.gcl$n

## Remove individuals with >20% missing data
SCHIG15_5_Jul18_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG15_5_Jul18", loci = loci24, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG15_5_Jul18_PostMissLoci <- SCHIG15_5_Jul18.gcl$n

SCHIG15_5_Jul18_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, dimnames = list("SCHIG15_5_Jul18", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG15_5_Jul18_SampleSizes[, "Genotyped"] <- Original_SCHIG15_5_Jul18_ColSize
SCHIG15_5_Jul18_SampleSizes[, "Missing"] <- Original_SCHIG15_5_Jul18_ColSize - ColSize_SCHIG15_5_Jul18_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG15_5_Jul18_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SCHIG15_5_Jul18", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG15_5_Jul18_RemovedDups <- RemoveDups.GCL(SCHIG15_5_Jul18_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG15_5_Jul18_PostDuplicate <- SCHIG15_5_Jul18.gcl$n

SCHIG15_5_Jul18_SampleSizes[, "Duplicate"] <- ColSize_SCHIG15_5_Jul18_PostMissLoci-ColSize_SCHIG15_5_Jul18_PostDuplicate
SCHIG15_5_Jul18_SampleSizes[, "Final"] <- ColSize_SCHIG15_5_Jul18_PostDuplicate

SCHIG15_5_Jul18_SampleSizes

write.xlsx(SCHIG15_5_Jul18_SampleSizes, file = "Output/SCHIG15_5_Jul18_SampleSizes.xlsx")
dput(x = SCHIG15_1_Jun27.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG15_1_Jun27_IDs.txt")

## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG15_5_Jul18", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG15_5_Jul18", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG15_5_Jul18", loci = loci24[-23], path = "Genepop/SCHIG15_5_Jul18_23nuclearloci.gen", VialNums = TRUE)
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")

# Read in Genepop output .P file
HWE <- ReadGenepopHWEKS.GCL(file = "Genepop/SCHIG15_5_Jul18_23nuclearloci.txt.P")

# Plot Fis values
plot(sort(HWE[, "WC Fis"]), type = "h", lwd = 10, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40"); abline(h = 0, lwd = 5)

# Look at data for any markers out of HWE
HWE[HWE$PValue < 0.05, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get baseline objects needed for MSA
Chignik7Populations <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Groupvec7.txt")
Chignik24BaselineFormat <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik24BaselineFormat.txt")
Inits <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Inits.txt")
loci24MSA <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/loci24.txt")

## Defining the random seeds as the same as WASSIP mixtures for repeatability.
WASSIPSockeyeSeeds <- dget(file = "V:/WORK/WASSIP/Sockeye/Mixture/Objects/WASSIPSockeyeSeeds.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 5 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG15_5_Jul18", loci = loci24MSA, IDs = NULL, mixname = "SCHIG15_5_Jul18",
                                          dir = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikMixtureFormat.txt")

ChignikGroups <- dget(file = "V:/WORK/Sockeye/Chignik/2013 Chignik Inseason/Objects/ChignikGroups.txt")

## Get prior from previous round (sequential priors)
Chignik2015Period5Prior <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015Period5Prior.txt")

## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG15_5_Jul18", basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2015Period5Prior, initmat = Inits, dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") ## Switch is on for outputing CLS file with thinning by 100

## Create output directory
dir.create("BAYES/Output/SCHIG15_5_Jul18")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG15_5_Jul18_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG15_5_Jul18",
                                                          prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
dput(x = SCHIG15_5_Jul18_Estimates, file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates objects/SCHIG15_5_Jul18_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG15_5_Jul18_Estimates$Stats[[1]][, "GR"]

# View trace of 2nd half of all 5 chains
plot(SCHIG15_5_Jul18_Estimates[[2]][[1]][, 1], type = "l", ylim = c(0, 1), ylab = "Proportion Black Lake", xlab = "Interation")
abline(v = seq(from = 0, to = 1e5, by = 2e4), lwd = 2)
text(x = seq(from = 1e4, to = 9e4, by = 2e4), y = 1.1, xpd = TRUE, labels = paste("Chain", 1:5))

# Sneak Peak
SCHIG15_5_Jul18_Estimates$Stats[[1]]
hist(SCHIG15_5_Jul18_Estimates[[2]][[1]][, 1], col = 8, main = "Proportion Black Lake Posterior", xlab = "Proportion")
abline(v = mean(SCHIG15_5_Jul18_Estimates[[2]][[1]][, 1]), lwd = 3)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG15_5_Jul18_Estimates$Stats[[1]],
           file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates tables/SCHIG15_5_Jul18_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReport.f.txt")
# Old function if the new one doesn't work
# ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReportOLD.f.txt")

SCHIG15_5_Jul18_SampleSizes

ChignikInseasonReport.f(NewData = SCHIG15_5_Jul18_Estimates, Period = 5, NumSampled = 190, NumAnalyzed = SCHIG15_5_Jul18_SampleSizes[1, "Genotyped"],
                        Included = SCHIG15_5_Jul18_SampleSizes[1, "Final"], Month = "July", Day = 18)

# Confirm the right fish!!!
SCHIG15_5_Jul18.gcl$attributes$CAPTURE_DATE
SCHIG15_5_Jul18.gcl$attributes$FK_FISH_ID

## Use Adobe Acrobat Pro on Gen1 to add the footer of page, date, and time
## Create News Release Fishery Update based on the previous round, change the numbers manually in the text
## Post Fishery Update online using the News Release System

## Update 2015 Chignik Estimates Summary and send to Birch/Mary-Beth

## Prior for next round
Chignik2015Period6Prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG15_5_Jul18_Estimates$Stats[[1]][, 1], minval = 0.01)
dput(x = Chignik2015Period6Prior, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015Period6Prior.txt")

## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_5.RData")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 6 July 25 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run by KRS on Mon Jul 27 14:38:21 2015
date()

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("Z:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2015")
## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_6.RData")
## load("Z:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2015/2015ChignikInseason_6.RData")

#This sources all of the new GCL functions to this workspace
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")
source("V:/DATA/R_GEN/Edited GCL Source Scripts/ReadLOKIv2.GCL.r")

## ReadLOKI.GCL does not work in RStudio, have to go into base R
## Pull all data for each silly code and create .gcl objects for each
#ReadLOKI.GCL(sillyvec = "SCHIG15", markersuite = "Sockeye2013Chignik_24SNPs")
ReadLOKIv2.GCL(sillyvec = "SCHIG15", markersuite = "Sockeye2013Chignik_24SNPs")

objects(pattern = "\\.gcl") # check objects
str(SCHIG15.gcl)

## Create a locus list from Locus Control object
loci24 <- LocusControl$locusnames

## Save unaltered .gcl's as back-up:
dput(x = SCHIG15.gcl, file = paste("Raw genotypes/SCHIG15_6_Jul25.gcl.txt", sep = ''))

## Save original LocusControl
dput(x = LocusControl, file = "Objects/OriginalLocusControl_6_Jul25.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata-ID associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sample dates defined inseason.

## Chignik 2015.5 July 19
unique(SCHIG15.gcl$attributes$CAPTURE_DATE)
SCHIG15_6_Jul25IDs <- AttributesToIDs.GCL(silly = "SCHIG15", attribute = "CAPTURE_DATE", matching = unique(SCHIG15.gcl$attributes$CAPTURE_DATE)[6])

SCHIG15_6_Jul25IDs <- list(na.omit(SCHIG15_6_Jul25IDs))
names(SCHIG15_6_Jul25IDs) <- "SCHIG15"

PoolCollections.GCL("SCHIG15", loci = loci24, IDs = SCHIG15_6_Jul25IDs, newname = "SCHIG15_6_Jul25")
SCHIG15_6_Jul25.gcl$n ## 190


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('xlsx')

## Get sample size by locus
Original_SCHIG15_6_Jul25_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG15_6_Jul25", loci24)
min(Original_SCHIG15_6_Jul25_SampleSizebyLocus) ## Good? Yes, 181.
apply(Original_SCHIG15_6_Jul25_SampleSizebyLocus, 1, min) / apply(Original_SCHIG15_6_Jul25_SampleSizebyLocus, 1, max) # 0.968

## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG15_6_Jul25_ColSize <- SCHIG15_6_Jul25.gcl$n

## Remove individuals with >20% missing data
SCHIG15_6_Jul25_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG15_6_Jul25", loci = loci24, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG15_6_Jul25_PostMissLoci <- SCHIG15_6_Jul25.gcl$n

SCHIG15_6_Jul25_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, dimnames = list("SCHIG15_6_Jul25", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG15_6_Jul25_SampleSizes[, "Genotyped"] <- Original_SCHIG15_6_Jul25_ColSize
SCHIG15_6_Jul25_SampleSizes[, "Missing"] <- Original_SCHIG15_6_Jul25_ColSize - ColSize_SCHIG15_6_Jul25_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG15_6_Jul25_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SCHIG15_6_Jul25", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG15_6_Jul25_RemovedDups <- RemoveDups.GCL(SCHIG15_6_Jul25_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG15_6_Jul25_PostDuplicate <- SCHIG15_6_Jul25.gcl$n

SCHIG15_6_Jul25_SampleSizes[, "Duplicate"] <- ColSize_SCHIG15_6_Jul25_PostMissLoci-ColSize_SCHIG15_6_Jul25_PostDuplicate
SCHIG15_6_Jul25_SampleSizes[, "Final"] <- ColSize_SCHIG15_6_Jul25_PostDuplicate

SCHIG15_6_Jul25_SampleSizes

write.xlsx(SCHIG15_6_Jul25_SampleSizes, file = "Output/SCHIG15_6_Jul25_SampleSizes.xlsx")
dput(x = SCHIG15_1_Jun27.gcl$attributes$FK_FISH_ID, file = "Final Fish IDs/SCHIG15_1_Jun27_IDs.txt")

## Combine loci
LocusControl
CombineLoci.GCL(sillyvec = "SCHIG15_6_Jul25", markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = "SCHIG15_6_Jul25", markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = "SCHIG15_6_Jul25", loci = loci24[-23], path = "Genepop/SCHIG15_6_Jul25_23nuclearloci.gen", VialNums = TRUE)
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")

# Read in Genepop output .P file
HWE <- ReadGenepopHWEKS.GCL(file = "Genepop/SCHIG15_6_Jul25_23nuclearloci.txt.P")

# Plot Fis values
plot(sort(HWE[, "WC Fis"]), type = "h", lwd = 10, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40"); abline(h = 0, lwd = 5)

# Look at data for any markers out of HWE
HWE[HWE$PValue < 0.05, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get baseline objects needed for MSA
Chignik7Populations <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Groupvec7.txt")
Chignik24BaselineFormat <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Chignik24BaselineFormat.txt")
Inits <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/Inits.txt")
loci24MSA <- dget(file = "V:/WORK/Sockeye/Chignik/2012 Chignik Baseline/Objects/loci24.txt")

## Defining the random seeds as the same as WASSIP mixtures for repeatability.
WASSIPSockeyeSeeds <- dget(file = "V:/WORK/WASSIP/Sockeye/Mixture/Objects/WASSIPSockeyeSeeds.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 6 MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping mixture file:
ChignikMixtureFormat <- CreateMixture.GCL(sillys = "SCHIG15_6_Jul25", loci = loci24MSA, IDs = NULL, mixname = "SCHIG15_6_Jul25",
                                          dir = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(ChignikMixtureFormat, file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikMixtureFormat.txt")

ChignikGroups <- dget(file = "V:/WORK/Sockeye/Chignik/2013 Chignik Inseason/Objects/ChignikGroups.txt")

## Get prior from previous round (sequential priors)
Chignik2015Period6Prior <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/Chignik2015Period6Prior.txt")

## Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci24MSA, mixname = "SCHIG15_6_Jul25", basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = Chignik2015Period6Prior, initmat = Inits, dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") ## Switch is on for outputing CLS file with thinning by 100

## Create output directory
dir.create("BAYES/Output/SCHIG15_6_Jul25")

## Run BAYES, check Raftery-Lewis for each chain, and summarize stats and check for convergence (G-R) in 5th chain

## This is the summarizing and dputting of estimates
SCHIG15_6_Jul25_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG15_6_Jul25",
                                                          prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
dput(x = SCHIG15_6_Jul25_Estimates, file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates objects/SCHIG15_6_Jul25_Estimates.txt")

# Verify that Gelman-Rubin < 1.2
SCHIG15_6_Jul25_Estimates$Stats[[1]][, "GR"]

# View trace of 2nd half of all 5 chains
plot(SCHIG15_6_Jul25_Estimates[[2]][[1]][, 1], type = "l", ylim = c(0, 1), ylab = "Proportion Black Lake", xlab = "Interation")
abline(v = seq(from = 0, to = 1e5, by = 2e4), lwd = 2)
text(x = seq(from = 1e4, to = 9e4, by = 2e4), y = 1.1, xpd = TRUE, labels = paste("Chain", 1:5))

# Sneak Peak
SCHIG15_6_Jul25_Estimates$Stats[[1]]
hist(SCHIG15_6_Jul25_Estimates[[2]][[1]][, 1], col = 8, main = "Proportion Black Lake Posterior", xlab = "Proportion")
abline(v = mean(SCHIG15_6_Jul25_Estimates[[2]][[1]][, 1]), lwd = 3)

## This is the formatting and writing of tables of estimates 
write.xlsx(x = SCHIG15_6_Jul25_Estimates$Stats[[1]],
           file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates tables/SCHIG15_6_Jul25_Estimates Table.xlsx")

## Write Update Report
ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReport.f.txt")
# Old function if the new one doesn't work
# ChignikInseasonReport.f <- dget(file="V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/ChignikInseasonReportOLD.f.txt")

SCHIG15_6_Jul25_SampleSizes

ChignikInseasonReport.f(NewData = SCHIG15_6_Jul25_Estimates, Period = 6, NumSampled = 189, NumAnalyzed = SCHIG15_6_Jul25_SampleSizes[1, "Genotyped"],
                        Included = SCHIG15_6_Jul25_SampleSizes[1, "Final"], Month = "July", Day = 25)

# Confirm the right fish!!!
SCHIG15_6_Jul25.gcl$attributes$CAPTURE_DATE
SCHIG15_6_Jul25.gcl$attributes$FK_FISH_ID

## Use Adobe Acrobat Pro on Gen1 to add the footer of page, date, and time
## Create News Release Fishery Update based on the previous round, change the numbers manually in the text
## Post Fishery Update online using the News Release System

## Update 2015 Chignik Estimates Summary and send to Birch/Mary-Beth

## save.image("V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/2015ChignikInseason_6.RData")

require(beepr)
beep(8)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Hurray! You survived another year of Chignik in-season. Be ready for more big wind next year! ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sex Specific Estimates for Birch Foster ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2015")
source("V:\\Analysis\\R files\\Scripts\\PROD\\Functions.GCL.r")
source("H:/R Source Scripts/Functions.GCL_KS.R")
# source("V:/DATA/R_GEN/tempGCL Source Scripts/IndividualAssignmentSummaryKS.GCL.R")


# Get the mixnames for each round of MSA
SCHIG15_mixnames <- list.dirs(path = "BAYES/Output", full.names = FALSE, recursive = FALSE)

#~~~~~~~~~~~~~~~~~~
## Get fish IDs for mixture fish
# dir.create("Final Fish IDs")
# 
# for(i in seq(SCHIG15_mixnames)){
#   load(paste("2015ChignikInseason_", i, ".RData", sep = ''))
#   dput(x = get(paste(SCHIG15_mixnames[i], ".gcl", sep = ''))$attributes$FK_FISH_ID, file = paste("Final Fish IDs/", SCHIG15_mixnames[i], "_IDs.txt", sep = ''))
# }; rm(list = ls(all = TRUE))
#~~~~~~~~~~~~~~~~~~




# Get baseline objects
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017")
Chignik7Populations <- dget(file = "Baseline 2012/Objects/Chignik7Populations.txt")
Groupvec7 <- dget(file = "Baseline 2012/Objects/Groupvec7.txt")
loci24MSA <- dget(file = "Baseline 2012/Objects/loci24.txt")
ChignikGroups <- dget(file = "Mixtures/2013/Objects/ChignikGroups.txt")
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2015")

# Get mixture objects
SCHIG15_IDs <- sapply(SCHIG15_mixnames, function(mix) {dget(file = paste("Final Fish IDs/", mix, "_IDs.txt", sep = ''))} )


# #!#!#!# Note that I did NOT thin the first round, but did thin by 100 the next 6 rounds #!#!#!#
# # Just trying for Round 1 on June 27 (no thinning)
# SCHIG15_1_Jun27_IndvAssign <- IndividualAssignmentSummaryKS.GCL(GroupNames = ChignikGroups, groupvec = Groupvec7, SCHIG15_mixnames = SCHIG15_mixnames[1], BAYESoutputDir = "BAYES/Output", nchains = 5, nreps = 40000, burn = 1/2, thin = 1)
# str(SCHIG15_1_Jun27_IndvAssign)
# hist(SCHIG15_1_Jun27_IndvAssign[[1]][, 1])
# 
# # Doing Round 2 - 6 now (thin = 100)
# SCHIG15_2_6_IndvAssign <- IndividualAssignmentSummaryKS.GCL(GroupNames = ChignikGroups, groupvec = Groupvec7, SCHIG15_mixnames = SCHIG15_mixnames[-1], BAYESoutputDir = "BAYES/Output", nchains = 5, nreps = 40000, burn = 1/2, thin = 100)
# str(SCHIG15_2_6_IndvAssign)
# 
# # Push all 6 mixtures into 1 list object
# SCHIG15_IndvAssign <- c(SCHIG15_1_Jun27_IndvAssign, SCHIG15_2_6_IndvAssign)
# str(SCHIG15_IndvAssign)
# sapply(SCHIG15_IndvAssign, dim)
# sapply(SCHIG15_IDs, length)
# 
# # Rename dimnames so they reflect Fish ID (dna_vial)
# for(i in seq(SCHIG15_mixnames)) {
#   dimnames(SCHIG15_IndvAssign[[i]])[[1]] = SCHIG15_IDs[[i]]
# }
# 
# 
# # Histogram of probablility Black Lake
# invisible(sapply(SCHIG15_mixnames, function(mixture) {hist(SCHIG15_IndvAssign[[mixture]][, "Black Lake"], col = 8, xlab = "Probability Black Lake (individuals)", main = mixture)} ))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get ASL metadata
ASL <- read.table(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2015/chignikISG2015ASLdata.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
str(ASL)
table(ASL$sex, ASL$sample_date)

ASL <- ASL[, c("dna_vial", "sample_date", "sex")]
str(ASL)

# SCHIG15_IndvAssign_Sex <- lapply(SCHIG15_IndvAssign, function(mix) {
#   cbind(mix, "sex" = ASL[match(rownames(mix), table = ASL$dna_vial), "sex"])
#   } )
# str(SCHIG15_IndvAssign_Sex)
# 
# # Plot male proportion over time
# plot(sapply(SCHIG15_IndvAssign_Sex, function(mix) {table(mix[, "sex"])["1"] / sum(table(mix[, "sex"])[c("1", "2")])} ), type = 'o', pch = 16, cex = 2, lwd = 3, ylim = c(0, 0.6), ylab = "Proportion Male")

SCHIG15_IDs_Sex <- sapply(SCHIG15_IDs, function(mix) {
  cbind("dna_vial" = mix, "sex" = ASL[match(mix, table = ASL$dna_vial), "sex"])
}, simplify = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Attempt Jim's Model

GroupNames = ChignikGroups; groupvec = Groupvec7; mixnames = SCHIG15_mixnames[2:6]; BAYESoutputDir = "BAYES/Output"; nchains = 5; nreps = 40000; burn = 1/2; thin = 100; threshhold=5e-7

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
    
    zgroup.mat <- array(data = zgroup, dim = c(mixsampsize, nits), dimnames = list(SCHIG15_IDs[[mixname]], seq(nits)))
    
    # t(apply(zgroup.mat, 1, function(ind) {tabulate(ind, G) / nits} ))
    
    indv.males <- which(SCHIG15_IDs_Sex[[mixname]][, "sex"] == 1)
    indv.females <- which(SCHIG15_IDs_Sex[[mixname]][, "sex"] == 2)
    
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


SCHIG15_2_6_IndvAssign_Sex <- ChignikMSAbySex.f(GroupNames = ChignikGroups, groupvec = Groupvec7, mixnames = SCHIG15_mixnames[2:6], BAYESoutputDir = "BAYES/Output", nchains = 5, nreps = 40000, burn = 1/2, thin = 100, threshhold = 5e-7)
SCHIG15_1_Jun27_IndvAssign_Sex <- ChignikMSAbySex.f(GroupNames = ChignikGroups, groupvec = Groupvec7, mixnames = SCHIG15_mixnames[1], BAYESoutputDir = "BAYES/Output", nchains = 5, nreps = 40000, burn = 1/2, thin = 1, threshhold = 5e-7)

SCHIG15_IndvAssign_Sex_BAYES <- c(SCHIG15_1_Jun27_IndvAssign_Sex, SCHIG15_2_6_IndvAssign_Sex)
str(SCHIG15_IndvAssign_Sex_BAYES)

SCHIG15_sex_date.mat <- table(ASL$sample_date, ASL$sex)

dput(x = SCHIG15_IndvAssign_Sex_BAYES, file = "Objects/SCHIG15_IndvAssign_Sex_BAYES.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots results

ASL <- read.table(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2015/chignikISG2015ASLdata.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
ASL <- ASL[, c("dna_vial", "sample_date", "sex")]
SCHIG15_sex_date.mat <- table(ASL$sample_date, ASL$sex)
SCHIG15_IndvAssign_Sex_BAYES <- dget(file = "Objects/SCHIG15_IndvAssign_Sex_BAYES.txt")
GeneticEstimates2015 <- dget(file = "Objects/GeneticEstimates2015.txt")


library(devEMF)
emf(file = paste("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2015/SexSpecificEstimates2.emf", sep = ''), width = 7, height = 6, family = "Times")

# par(bg=rgb(red=31,green=73,blue=125,maxColorValue=255))
par(bg = "white")
par(col = "white", col.axis = "white", col.lab = "white")
cex.lab = 1.5
col.axs = "black"
col.fem = "deeppink1"
col.male = "navy"
col.overall = "grey40"

  
  plot(y = sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Overall"]]["Black Lake", "mean"]} ), x = GeneticEstimates2015[, "Day"], type = "o", col = "blue", cex = 2, lwd = 3, pch = 16, ylim = c(0, 1), xlab = "", ylab = "", main = "", bty = 'n', axes = FALSE, xlim = c(min(x = GeneticEstimates2015[, "Day"]) - 1, max(x = GeneticEstimates2015[, "Day"]) + 2), cex.lab = cex.lab)
  axis(side = 2, col = col.axs, col.axis = col.axs, cex.axis = cex.lab)
  axis(side = 1, at = GeneticEstimates2015[, "Day"], labels = rep('', 6), col = col.axs)
  text(x = GeneticEstimates2015[, "Day"], y = -0.15, labels = rownames(SCHIG15_sex_date.mat), srt = 45, xpd = TRUE, col = col.axs)
  mtext(text = "Date", side = 1, line = 4, cex = cex.lab, col = col.axs)
  mtext(text = "Proportion Black Lake", side = 2, line = 3, cex = cex.lab, col = col.axs)

  arrows(x0 = GeneticEstimates2015[, "Day"] + 0.2, x1 = GeneticEstimates2015[, "Day"] + 0.2, y0 = sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Male"]]["Black Lake", "5%"]} ), sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Male"]]["Black Lake", "95%"]} ), angle = 90, col = col.male, code = 3, lwd = 2, length = 0.15)
  arrows(x0 = GeneticEstimates2015[, "Day"] - 0.2, x1 = GeneticEstimates2015[, "Day"] - 0.2, y0 = sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Female"]]["Black Lake", "5%"]} ), sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Female"]]["Black Lake", "95%"]} ), angle = 90, col = col.fem, code = 3, lwd = 2, length = 0.15)
  arrows(x0 = GeneticEstimates2015[, "Day"], x1 = GeneticEstimates2015[, "Day"], y0 = sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Overall"]]["Black Lake", "5%"]} ), sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Overall"]]["Black Lake", "95%"]} ), angle = 90, col = col.overall, code = 3, lwd = 2, length = 0.15)
  
  points(x = GeneticEstimates2015[, "Day"], y = sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Male"]]["Black Lake", "mean"]} ), type = "o", col = col.male, cex = 2, lwd = 3, pch = 16)
  points(x = GeneticEstimates2015[, "Day"], y = sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Female"]]["Black Lake", "mean"]} ), type = "o", col = col.fem, cex = 2, lwd = 3, pch = 16)
  points(x = GeneticEstimates2015[, "Day"], y = sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Overall"]]["Black Lake", "mean"]} ), type = "o", col = col.overall, cex = 2, lwd = 3, pch = 16)
  
  text(x = GeneticEstimates2015[, "Day"]-1, y = rep(0.15, 6), labels = SCHIG15_sex_date.mat[, "1"], col = col.male, cex = SCHIG15_sex_date.mat[, "1"] / rowSums(SCHIG15_sex_date.mat[, c("1", "2")]) * 2, adj = 0)
  text(x = GeneticEstimates2015[, "Day"]-1, y = rep(0.1, 6), labels = SCHIG15_sex_date.mat[, "2"], col = col.fem, cex = SCHIG15_sex_date.mat[, "2"] / rowSums(SCHIG15_sex_date.mat[, c("1", "2")]) * 2, adj = 0.1)
  text(x = GeneticEstimates2015[, "Day"]-1, y = rep(0.05, 6), labels = rowSums(SCHIG15_sex_date.mat[, c("1", "2")]), col = col.overall, cex = 2, adj = 0.1)
  
  
  legend("topright", legend = c("Overall", "Male", "Female"), fill = c(col.overall, col.male, col.fem), bty = 'n', cex = 1.5, text.col = col.axs)
#   legend("topright", legend = c("Overall"), fill = c(col.overall), bty = 'n', cex = 1.5, text.col = col.axs)

dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


round(GeneticEstimates2015, 2)[, c("5%", "95%")]
round(t(sapply(SCHIG15_IndvAssign_Sex_BAYES, function(period) {period[["Overall"]]["Black Lake",]} )), 2)[, c("5%", "95%")]




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Genotype frequencies for 2015 for One_U1004-183 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Original_SampleSizebyLocus <- SampSizeByLocus.GCL("SCHIG15", loci24)


## Get number of individuals per silly before removing missing loci individuals
Original_SCHIG15_ColSize <- SCHIG15.gcl$n

## Remove individuals with >20% missing data
SCHIG15_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG15", loci = loci24, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG15_PostMissLoci <- SCHIG15.gcl$n

SCHIG15_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, dimnames = list("SCHIG15_", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG15_SampleSizes[, "Genotyped"] <- Original_SCHIG15_ColSize
SCHIG15_SampleSizes[, "Missing"] <- Original_SCHIG15_ColSize - ColSize_SCHIG15_PostMissLoci

## Check within collections for duplicate individuals.
SCHIG15_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SCHIG15", loci = loci24, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
SCHIG15_RemovedDups <- RemoveDups.GCL(SCHIG15_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG15_PostDuplicate <- SCHIG15.gcl$n

SCHIG15_SampleSizes[, "Duplicate"] <- ColSize_SCHIG15_PostMissLoci-ColSize_SCHIG15_PostDuplicate
SCHIG15_SampleSizes[, "Final"] <- ColSize_SCHIG15_PostDuplicate


# Get IDs for each date

chig_dates <- unique(SCHIG15.gcl$attributes$CAPTURE_DATE)
chig_dates <- gsub(pattern = " AKDT", replacement = '', x = chig_dates)

SCHIG15_date_IDs <- NULL
for(i in seq(length(chig_dates))) {
  SCHIG15_date_IDs[[i]] <- AttributesToIDs.GCL(silly = "SCHIG15", attribute = "CAPTURE_DATE", matching = unique(SCHIG15.gcl$attributes$CAPTURE_DATE)[i])
}

names(SCHIG15_date_IDs) <- chig_dates

lapply(SCHIG15_date_IDs, function(ids) {table(SCHIG15.gcl$counts[ids, "One_U1004-183", 1])} )
lapply(SCHIG15_date_IDs, function(ids) {sum(table(SCHIG15.gcl$counts[ids, "One_U1004-183", 1]))} )
genofreqs <- sapply(SCHIG15_date_IDs, function(ids) {round(table(SCHIG15.gcl$counts[ids, "One_U1004-183", 1]) / sum(table(SCHIG15.gcl$counts[ids, "One_U1004-183", 1])), 2)} )

GeneticEstimates2015 <- dget(file = "Objects/GeneticEstimates2015.txt")

par(mar = c(7.1, 5.1, 1.1, 1.1))
plot(y = genofreqs["0", ], x = GeneticEstimates2015[, "Day"], col = "blue", ylim = c(0, 1), type = "o", pch = 16, cex = 2, lwd = 3, axes = FALSE, bty = 'n', xlab = '', ylab = "Genotype Frequency", cex.lab = 2)
mtext(text = "Date", side = 1, line = 6, cex =2)
axis(side = 2, cex.axis = 1.5)
axis(side = 1, at = GeneticEstimates2015[, "Day"], labels = rep('', 6))
text(x = GeneticEstimates2015[, "Day"], y = -0.15, labels = chig_dates <- gsub(pattern = "2015-", replacement = '', x = chig_dates), xpd = TRUE, srt = 45, cex = 2)
points(y = genofreqs["1", ], x = GeneticEstimates2015[, "Day"], col = "green", type = "o", pch = 16, cex = 2, lwd = 3)
points(y = genofreqs["2", ], x = GeneticEstimates2015[, "Day"], col = "red", type = "o", pch = 16, cex = 2, lwd = 3)
legend("topleft", legend = c("aa", "Aa", "AA"), fill = c("blue", "green", "red"), bty = 'n', title = "One_U1004-183", cex = 2, xjust = 0, )


barplot(height = genofreqs, col = c("blue", "green", "red"))
