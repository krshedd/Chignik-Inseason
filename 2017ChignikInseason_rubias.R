# Eric Anderson's 'rubias' GSI package in R
# Testing with 2017 Chignik sockeye salmon data, run all 6 mixtures at once!

devtools::install_github("eriqande/rubias")
require(rubias)
?rubias

# chinook is the example baseline ('reference') data.frame
head(chinook[, 1:8])
str(chinook)

# chinook_mix is the example mixture data.frame
head(chinook_mix[, 1:8])
str(chinook_mix)

# Load gcl functions
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Convert Baseline to Rubias ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Baseline 2012/")
load("Chignik2012SockeyeBaseline.RData")

objects(pattern = "\\.gcl")
Chignik7Populations
loci24
table(LocusControl$ploidy[loci24])

Groupvec7
Groups

SBROAD97.SBSPR97.SBOUL97.SFAN97.SALEC97.gcl$scores[1, , ]

require(memisc)
SBROAD97.SBSPR97.SBOUL97.SFAN97.SALEC97.df <- to.data.frame(X = SBROAD97.SBSPR97.SBOUL97.SFAN97.SALEC97.gcl$scores, as.vars = 2)
str(SBROAD97.SBSPR97.SBOUL97.SFAN97.SALEC97.df) # not what we need


create_rubias_baseline <- function(sillyvec, loci, group_names, groupvec) {
  silly_base.lst <- lapply(sillyvec, function(silly) {
    my.gcl <- get(paste0(silly, ".gcl"))
    scores.mat <- t(apply(my.gcl$scores[, loci, ], 1, function(ind) {c(t(ind))} ))
    colnames(scores.mat) <- as.vector(sapply(loci, function(locus) {c(locus, paste(locus, 1, sep = "."))} ))
    scores.df <- data.frame(scores.mat, stringsAsFactors = FALSE)
    scores.df$sample_type <- "reference"
    scores.df$repunit <- group_names[groupvec[sillyvec == silly]]
    scores.df$collection <- silly
    scores.df$indiv <- as.character(my.gcl$attributes$SillySource)
    silly_base.df <- scores.df[, c("sample_type", "repunit", "collection", "indiv", gsub(pattern = "-", replacement = ".", x = colnames(scores.mat)))] } #silly
    )
  return(do.call("rbind", silly_base.lst))
}

chignik_7pops_22loci.rubias_base <- create_rubias_baseline(sillyvec = Chignik7Populations, loci = loci24, group_names = Groups, groupvec = Groupvec7)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Convert Mixture to Rubias ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
