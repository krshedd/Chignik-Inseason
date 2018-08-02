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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Convert Baseline to `rubias` ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Baseline 2012/")
load("Chignik2012SockeyeBaseline.RData")

## Load gcl functions
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

## Check objects
objects(pattern = "\\.gcl")
Chignik7Populations
loci22 <- loci24  # aka loci24MSA in the inseason script
table(LocusControl$ploidy[loci22])

Groupvec7
Groups

SBROAD97.SBSPR97.SBOUL97.SFAN97.SALEC97.gcl$scores[1, , ]

chignik_anderson <- Anderson_etal.GCL(popvec = Chignik7Populations, loci = loci22, groups = Groupvec7, group_names = Groups)  # doesn't work for combined markers
chignik_leaveoneout <- LeaveOneOutDist.GCL(sillyvec = Chignik7Populations, loci = loci22, groupvec = Groupvec7)
chignik_confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = chignik_leaveoneout, groupnames = Groups, groupvec = Groupvec7, sillyvec = Chignik7Populations)
chignik_confusion$GroupByGroup

loci19 <- loci22[-c(4, 6, 8)]
chignik_anderson_loci19 <- Anderson_etal.GCL(popvec = Chignik7Populations, loci = loci19, groups = Groupvec7, group_names = Groups)  # doesn't work for combined markers
chignik_leaveoneout_loci19 <- LeaveOneOutDist.GCL(sillyvec = Chignik7Populations, loci = loci19, groupvec = Groupvec7)
chignik_confusion_loci19 <- ConfusionMatrices.GCL(LeaveOneOutDist = chignik_leaveoneout_loci19, groupnames = Groups, groupvec = Groupvec7, sillyvec = Chignik7Populations)

chignik_7pops_19loci.rubias_base <- create_rubias_baseline(sillyvec = Chignik7Populations, loci = loci19, group_names = Groups, groupvec = Groupvec7)
chignik_7pops_19loci.rubias_base_sa <- self_assign(reference = chignik_7pops_19loci.rubias_base, gen_start_col = 5)
chignik_7pops_19loci.rubias_base_sa %>% 
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood)) %>% 
  group_by(collection, repunit, inferred_repunit) %>% 
  summarise("mean_repu_scaled_like" = mean(repu_scaled_like)) %>% 
  spread(inferred_repunit, mean_repu_scaled_like)
rm(list = setdiff(ls(), c("chignik_7pops_22loci.rubias_base", "chignik_7pops_22loci.rubias_base_sa", "chignik_7pops_19loci.rubias_base", "chignik_7pops_19loci.rubias_base_sa")))
save.image("chignik_sockeye_rubias_example.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create function for rubias baseline
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

chignik_7pops_22loci.rubias_base <- create_rubias_baseline(sillyvec = Chignik7Populations, loci = loci22, group_names = Groups, groupvec = Groupvec7)
str(chignik_7pops_22loci.rubias_base)
dput(x = chignik_7pops_22loci.rubias_base, file = "Objects/chignik_7pops_22loci.rubias_base.txt")
chignik_7pops_22loci.rubias_base <- dget(file = "Objects/chignik_7pops_22loci.rubias_base.txt")

rm(list = setdiff(ls(), c("chignik_7pops_22loci.rubias_base", "Chignik7Populations", "loci22", "Groupvec7", "Groups")))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Convert Mixture to `rubias` ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures/2017/")
# load("2017ChignikInseason_6.RData")

## Load gcl functions
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get genotypes from LOKI
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = "krshedd", password = password)
loci24 <- LocusControl$locusnames
LOKI2R.GCL(sillyvec = "SCHIG17", username = "krshedd", password = password)
rm(password)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Data QC
# Get number of individuals per silly before removing missing loci individuals
Original_SCHIG17_ColSize <- SCHIG17.gcl$n

# Remove individuals with >20% missing data
SCHIG17_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SCHIG17", proportion = 0.8)

# Get number of individuals per silly after removing missing loci individuals
ColSize_SCHIG17_PostMissLoci <- SCHIG17.gcl$n

SCHIG17_SampleSizes <- matrix(data = NA, nrow = 1, ncol = 4, 
                              dimnames = list("SCHIG17", c("Genotyped", "Missing", "Duplicate", "Final")))
SCHIG17_SampleSizes[, "Genotyped"] <- Original_SCHIG17_ColSize
SCHIG17_SampleSizes[, "Missing"] <- Original_SCHIG17_ColSize - ColSize_SCHIG17_PostMissLoci

# Check within collections for duplicate individuals.
SCHIG17_DuplicateCheck95MinProportion <- 
  CheckDupWithinSilly.GCL(sillyvec = "SCHIG17", loci = loci24, quantile = NULL, minproportion = 0.95)

# Remove duplicate individuals
SCHIG17_RemovedDups <- RemoveDups.GCL(SCHIG17_DuplicateCheck95MinProportion)

# Get number of individuals per silly after removing duplicate individuals
ColSize_SCHIG17_PostDuplicate <- SCHIG17.gcl$n

SCHIG17_SampleSizes[, "Duplicate"] <- ColSize_SCHIG17_PostMissLoci-ColSize_SCHIG17_PostDuplicate
SCHIG17_SampleSizes[, "Final"] <- ColSize_SCHIG17_PostDuplicate
SCHIG17_SampleSizes


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Define Strata
PoolCollectionsByDateDF <- function(silly, date.df, loci) {
  sapply(silly, function(mix) {
    mix.dates <- unique(as.Date(get(paste0(mix, ".gcl"))$attributes$CAPTURE_DATE))
    by(data = date.df, INDICES = date.df$Strata, function(x) {
      IDs <- AttributesToIDs.GCL(silly = mix, attribute = "CAPTURE_DATE", matching = mix.dates[mix.dates >= x$Begin & mix.dates <= x$End])
      IDs <- list(na.omit(IDs))
      names(IDs) <- mix
      PoolCollections.GCL(collections = mix, loci = loci, IDs = IDs, newname = paste(mix, as.character(x$Strata), sep = "_"))
      list("First Last Fish" = range(as.numeric(unlist(IDs))), "n" = get(paste0(mix, "_", as.character(x$Strata), ".gcl"))$n)
    } )
  }, simplify = FALSE, USE.NAMES = TRUE)
}

Chignik2017_date.df <- data.frame("Strata" = factor(paste0("Strata", 1:6)),
                                  "Begin" = as.Date(c("2017-06-25", "2017-07-01", "2017-07-07", "2017-07-13", "2017-07-18", "2017-07-23")),
                                  "End" = as.Date(c("2017-06-26", "2017-07-01", "2017-07-08", "2017-07-13", "2017-07-18", "2017-07-23")))
str(Chignik2017_date.df)

PoolCollectionsByDateDF(silly = "SCHIG17", date.df = Chignik2017_date.df, loci = loci24)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Combine loci
CombineLoci.GCL(sillyvec = paste0("SCHIG17_Strata", 1:6), markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = paste0("SCHIG17_Strata", 1:6), markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create function for rubias mixture
# ASSUMES each silly is its own mixture
create_rubias_mixture <- function(sillyvec, loci) {
  silly_mix.lst <- lapply(sillyvec, function(silly) {
    my.gcl <- get(paste0(silly, ".gcl"))
    scores.mat <- t(apply(my.gcl$scores[, loci, ], 1, function(ind) {c(t(ind))} ))
    colnames(scores.mat) <- as.vector(sapply(loci, function(locus) {c(locus, paste(locus, 1, sep = "."))} ))
    scores.df <- data.frame(scores.mat, stringsAsFactors = FALSE)
    scores.df$sample_type <- "mixture"
    scores.df$repunit <- NA
    mode(scores.df$repunit) <- "character"
    scores.df$collection <- silly
    scores.df$indiv <- as.character(my.gcl$attributes$SillySource)
    silly_mix.df <- scores.df[, c("sample_type", "repunit", "collection", "indiv", gsub(pattern = "-", replacement = ".", x = colnames(scores.mat)))] } #silly
  )
  return(do.call("rbind", silly_mix.lst))
}

chignik_2017.rubias_mix <- create_rubias_mixture(sillyvec = paste0("SCHIG17_Strata", 1:6), loci = loci22)
str(chignik_2017.rubias_mix)
dput(x = chignik_2017.rubias_mix, file = "Objects/chignik_2017.rubias_mix.txt")

rm(list = setdiff(ls(), c("chignik_2017.rubias_mix", "chignik_7pops_22loci.rubias_base", "Chignik7Populations", "loci22", "Groupvec7", "Groups", "SCHIG17_Strata1.gcl", "SCHIG17_Strata2.gcl", "SCHIG17_Strata3.gcl", "SCHIG17_Strata4.gcl", "SCHIG17_Strata5.gcl", "SCHIG17_Strata6.gcl", "SCHIG17.gcl", "Chignik2017_date.df")))
# save.image("2017ChignikInseason_rubias.RData")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Run `rubias` ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(rubias)
require(tidyverse)

chignik_2017_mix_est <- infer_mixture(reference = chignik_7pops_22loci.rubias_base, 
                                      mixture = chignik_2017.rubias_mix, 
                                      gen_start_col = 5)
str(chignik_2017_mix_est, max.level = 2)
dput(x = chignik_2017_mix_est, file = "Objects/chignik_2017_mix_est.txt")

rep_mix_ests <- chignik_2017_mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  # adding mixing proportions over collections in the repunit

spread(data = rep_mix_ests, key = repunit, value = repprop)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Baseline testing with `rubias` ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Self assignment
chignik_7pops_22loci.rubias_base_sa <- self_assign(reference = chignik_7pops_22loci.rubias_base, gen_start_col = 5)
str(chignik_7pops_22loci.rubias_base_sa, max.level = 1)

black2black_sa.df <- chignik_7pops_22loci.rubias_base_sa %>% 
  filter(collection == "SBROAD97.SBSPR97.SBOUL97.SFAN97.SALEC97") %>% 
  filter(inferred_collection == "SBROAD97.SBSPR97.SBOUL97.SFAN97.SALEC97")
plot(x = black2black_sa.df$scaled_likelihood,
     y = chignik_leaveoneout[[2]]$SBROAD97.SBSPR97.SBOUL97.SFAN97.SALEC97$SBROAD97.SBSPR97.SBOUL97.SFAN97.SALEC97)

chignik2chignik_sa.df <- chignik_7pops_22loci.rubias_base_sa %>% 
  filter(collection == "SCHIA08.SCHIA97E.SCHIA97M") %>% 
  filter(inferred_collection == "SCHIA08.SCHIA97E.SCHIA97M")
plot(x = chignik2chignik_sa.df$scaled_likelihood,
     y = chignik_leaveoneout[[2]]$SCHIA08.SCHIA97E.SCHIA97M$SCHIA08.SCHIA97E.SCHIA97M)

sa_to_repu <- chignik_7pops_22loci.rubias_base_sa %>% 
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))

str(sa_to_repu)
sa_to_repu %>% 
  group_by(collection, repunit, inferred_repunit) %>% 
  summarise("mean_repu_scaled_like" = mean(repu_scaled_like)) %>% 
  spread(inferred_repunit, mean_repu_scaled_like)

## Leave-one-out
chignik_7pops_22loci.rubias_base_loo <- assess_reference_loo(reference = chignik_7pops_22loci.rubias_base, gen_start_col = 5, reps = 200, mixsize = 200)
chignik_7pops_22loci.rubias_base_loo

tmp_22_loo <- chignik_7pops_22loci.rubias_base_loo %>% 
  group_by(iter, repunit) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))
tmp_22_loo

# Plot
ggplot(tmp_22_loo, aes(x = true_repprop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)

ggplot(tmp_22_loo, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)

tmp_22_loo %>% 
  filter(repunit == "BlackLake") %>% 
  mutate("resid" = repu_n_prop - repprop_posterior_mean) %>% 
  ggplot(aes(x = resid)) + geom_histogram()

## MC
chignik_7pops_22loci.rubias_base_mc <- assess_reference_mc(reference = chignik_7pops_22loci.rubias_base, gen_start_col = 5, reps = 50, mixsize = 200)
chignik_7pops_22loci.rubias_base_mc

tmp_22_mc <- chignik_7pops_22loci.rubias_base_mc %>% 
  group_by(iter, repunit) %>% 
  summarise(true_repprop = sum(omega), repprop_posterior_mean = sum(post_mean), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))
tmp_22_mc

# Plot
ggplot(tmp_22_mc, aes(x = true_repprop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)

ggplot(tmp_22_mc, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 19 diploid SNPs
## Leave-one-out
chignik_7pops_19loci.rubias_base_loo <- assess_reference_loo(reference = chignik_7pops_19loci.rubias_base, gen_start_col = 5, reps = 50, mixsize = 200)
chignik_7pops_19loci.rubias_base_loo

tmp_19_loo <- chignik_7pops_19loci.rubias_base_loo %>% 
  group_by(iter, repunit) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))
tmp_19_loo

# Plot
ggplot(tmp_19_loo, aes(x = true_repprop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)

ggplot(tmp_19_loo, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)


## MC
chignik_7pops_19loci.rubias_base_mc <- assess_reference_mc(reference = chignik_7pops_19loci.rubias_base, gen_start_col = 5, reps = 50, mixsize = 200)
chignik_7pops_19loci.rubias_base_mc

tmp_19_mc <- chignik_7pops_19loci.rubias_base_mc %>% 
  group_by(iter, repunit) %>% 
  summarise(true_repprop = sum(omega), repprop_posterior_mean = sum(post_mean), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))
tmp_19_mc

# Plot
ggplot(tmp_19_mc, aes(x = true_repprop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)

ggplot(tmp_19_mc, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 17 diploid SNPs, no U1004, no Tf_ex
## Leave-one-out
chignik_7pops_17loci.rubias_base <- chignik_7pops_19loci.rubias_base[-c(grep(pattern = "U1004", x = colnames(chignik_7pops_19loci.rubias_base)), 
                                        grep(pattern = "Tf_ex", x = colnames(chignik_7pops_19loci.rubias_base)))]
chignik_7pops_17loci.rubias_base_loo <- assess_reference_loo(reference = chignik_7pops_17loci.rubias_base, gen_start_col = 5, reps = 50, mixsize = 200)
chignik_7pops_17loci.rubias_base_loo

tmp_17_loo <- chignik_7pops_17loci.rubias_base_loo %>% 
  group_by(iter, repunit) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))
tmp_17_loo

# Plot
ggplot(tmp_17_loo, aes(x = true_repprop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)

ggplot(tmp_17_loo, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)


## MC
chignik_7pops_17loci.rubias_base_mc <- assess_reference_mc(reference = chignik_7pops_17loci.rubias_base, gen_start_col = 5, reps = 50, mixsize = 200)
chignik_7pops_17loci.rubias_base_mc

tmp_17_mc <- chignik_7pops_17loci.rubias_base_mc %>% 
  group_by(iter, repunit) %>% 
  summarise(true_repprop = sum(omega), repprop_posterior_mean = sum(post_mean), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))
tmp_17_mc

# Plot
ggplot(tmp_17_mc, aes(x = true_repprop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)

ggplot(tmp_17_mc, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)









repu_scores <- list(all_22_loci = chignik_7pops_22loci.rubias_base_sa,
                    only_19_loci = chignik_7pops_19loci.rubias_base_sa) %>%
  bind_rows(.id = "DataSet") %>%
  group_by(DataSet, indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_lik = sum(scaled_likelihood)) %>%
  ungroup() %>%
  filter(repunit == inferred_repunit) %>%
  spread(DataSet, value = repu_scaled_lik)
# then plot them
ggplot(repu_scores, aes(x = only_19_loci, y = all_22_loci, colour = repunit)) + 
  geom_point() +
  facet_wrap(~ collection) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")
