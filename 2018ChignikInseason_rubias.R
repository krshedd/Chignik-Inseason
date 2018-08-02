#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Chignik Sockeye Postseason with `rubias` 2018 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This inseason analysis serves two purposes:
#   1) Create a rubias version of the 2012 Chignik baseline (7 pops, 24 SNPs),
#   2) Analyze inseason samples postseason against the 2012 Chignik baseline (7 pops, 24 SNPs),
#   3) Provide historical perspecitve (i.e. compare to previous years).

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create baseline and save objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2018/Mixtures/2018/")
load("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2018/Baseline 2012/Chignik2012SockeyeBaseline.RData")
rm(LocusControl)

source("C:/Users/krshedd//R/Functions.GCL.R")

Chignik7Populations
Chignik7CommonNames
Groups
Groupvec7
loci22 <- loci24  # aka loci24MSA in the inseason script
ChignikGroups <- c("Black Lake", "Chignik Lake")
CreateLocusControl.GCL(markersuite = "Sockeye2013Chignik_24SNPs", username = "krshedd", password = password)
loci24 <- LocusControl$locusnames
chignik_colors <- c("red", "blue")
Chignik24BaselineFormat
Inits
WASSIPSockeyeSeeds <- dget(file = "V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save baseline objects
dir.create("Objects")
save_objects(
  objects = c("Chignik7Populations", 
              "Chignik7CommonNames", 
              "Groups", 
              "ChignikGroups", 
              "Groupvec7", 
              "LocusControl", 
              "loci24", 
              "loci22",
              "chignik_colors",
              "Chignik24BaselineFormat",
              "Inits",
              "WASSIPSockeyeSeeds"), 
  path = "Objects/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save baseline genotypes
dir.create("Baseline genotpyes")
save_sillys(sillyvec = Chignik7Populations, path = "Baseline genotpyes/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create rubias baseline
dir.create("rubias")
dir.create("rubias/baseline")

chignik7_loci22.base <- create_rubias_baseline(sillyvec = Chignik7Populations, loci = loci22, group_names = ChignikGroups, groupvec = Groupvec7, path = "rubias/baseline/", baseline_name = "chignik7_loci22")
dput(x = chignik7_loci22.base, file = "rubias/baseline/chignik7_loci22.base.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Load baseline and objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2018/Mixtures/2018/")

library(tidyverse)
library(lubridate)
library(rubias)
library(ggthemes)
library(gridExtra)
library(ggpubr)
# library(conflicted)

source("C:/Users/krshedd//R/Functions.GCL.R")
load_objects(path = "Objects")
load_objects(path = "rubias/baseline")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Test baseline ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Leave one out testing
chignik7_loci22.base_loo <- assess_reference_loo(reference = chignik7_loci22.base, gen_start_col = 5, reps = 100, mixsize = 200)

# Summarize to reporting unit level
loo_loci22_out <- chignik7_loci22.base_loo %>% 
  mutate(repunit_f = factor(x = repunit, levels = ChignikGroups)) %>% 
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))

# Plot LOO results
ggplot(loo_loci22_out, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit_f)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(name = "Reporting Group", values = chignik_colors) +
  facet_wrap(~ repunit_f) +
  xlab("True Reporting Group Proportion") +
  ylab("Posterior Mean Reporting Group Proportion") +
  ggtitle("Chignik Leave-one-out Test Results")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 100% Leave one out testing
arep_100 <- sapply(ChignikGroups, function(x) {tibble(repunit = x, ppn = 1.0)}, simplify = FALSE)
chignik7_loci22.base_loo_100 <- assess_reference_loo(reference = chignik7_loci22.base, gen_start_col = 5, reps = 10, mixsize = 200, alpha_repunit = arep_100)

chignik7_loci22.base_loo_100 %>% 
  mutate(repunit_f = factor(x = repunit, levels = ChignikGroups)) %>% 
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n)) %>% 
  filter(repunit_scenario == repunit_f) %>% 
  ggplot(aes(x = iter, y = repprop_posterior_mean, fill = repunit_f)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0.9) +
  scale_fill_manual(name = "Reporting Group", values = chignik_colors) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  facet_wrap(~ repunit_f) +
  xlab("Iteration") +
  ylab("Posterior Mean Reporting Group Proportion") +
  ggtitle("Chignik 100% Proof Tests")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Re-analyze 2017 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Read in 2017 genotypes
sillys <- paste0("SCHIG17")
LOKI2R.GCL(sillyvec = sillys, username = "krshedd", password = password)
dir.create("Raw genotypes")
save_sillys(sillyvec = sillys, path = "Raw genotypes")
# load_sillys(path = "Raw genotypes", sillyvec = sillys)

sapply(sillys, function(silly) {get(paste0(silly, ".gcl"))$n} )  # 1140/yr
sapply(sillys, function(silly) {table(get(paste0(silly, ".gcl"))$attributes$CAPTURE_DATE, useNA = "always")} )  # 2016 is missing Capture Date

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Stratify mixtures
#~~~~~~~~~~~~~~~~~~
# Function to pool by sample date
PoolCollectionsByDateDF <- function(silly, date.df, loci) {
  sapply(silly, function(mix) {
    mix.dates <- unique(as.Date(get(paste(mix, ".gcl", sep = ''))$attributes$CAPTURE_DATE))
    by(data = date.df, INDICES = date.df$strata, function(x) {
      IDs <- AttributesToIDs.GCL(silly = mix, attribute = "CAPTURE_DATE", matching = mix.dates[mix.dates >= x$begin_date & mix.dates <= x$end_date])
      IDs <- list(na.omit(IDs))
      names(IDs) <- mix
      PoolCollections.GCL(collections = mix, loci = loci, IDs = IDs, newname = paste(mix, as.character(x$strata), sep = "_"))
      list("First Last Fish" = range(as.numeric(unlist(IDs))), "n" = get(paste(mix, "_", as.character(x$strata), ".gcl", sep = ''))$n)
    } )
  }, simplify = FALSE, USE.NAMES = TRUE)
}

date_2017.df <- tibble(strata = 1:6, 
                       begin_date = as.Date(x = c("2017-06-25", "2017-07-01", "2017-07-07", "2017-07-13", "2017-07-18", "2017-07-23")), 
                       end_date = as.Date(x = c("2017-06-26", "2017-07-01", "2017-07-08", "2017-07-13", "2017-07-18", "2017-07-23")))

PoolCollectionsByDateDF(silly = "SCHIG17", date.df = date_2017.df, loci = loci24)

sillys_strata_2017 <- setdiff(
  str_replace(string = objects(pattern = "\\.gcl"), 
              pattern = "\\.gcl", 
              replacement = ''),
  sillys)
save_objects(objects = "sillys_strata_2017", path = "Objects")

dir.create("Raw genotypes/Strata")
save_sillys(sillyvec = sillys_strata_2017, path = "Raw genotypes/Strata/")

sapply(sillys_strata_2017, function(silly) {table(get(paste0(silly, ".gcl"))$attributes$CAPTURE_DATE, useNA = "always")} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Data QC
sillys_strata_2017_n <- matrix(data = NA, nrow = length(sillys_strata_2017), ncol = 4, 
                                    dimnames = list(sillys_strata_2017, c("Genotyped", "Missing", "Duplicate", "Final")))

#~~~~~~~~~~~~~~~~~~
#### Check loci
## Get sample size by locus
original_sillys_strata_2017_n_locus <- SampSizeByLocus.GCL(sillyvec = sillys_strata_2017, loci = loci24)
min(original_sillys_strata_2017_n_locus)  ## 178
round(apply(original_sillys_strata_2017_n_locus, 1, function(locus) {min(locus) / max(locus)}), 2)

original_sillys_strata_2017_percent_locus <- apply(original_sillys_strata_2017_n_locus, 1, function(row) {row / max(row)} )
which(apply(original_sillys_strata_2017_percent_locus, 2, min) < 0.8)  # no re-runs!

# Genotpying percentage by locus and strata
require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(original_sillys_strata_2017_percent_locus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares

#~~~~~~~~~~~~~~~~~~
#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
sillys_strata_2017_n[, "Genotyped"] <- sapply(paste0(sillys_strata_2017, ".gcl"), function(x) get(x)$n)

### Missing
## Remove individuals with >20% missing data
sillys_strata_2017_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = sillys_strata_2017, proportion = 0.8)
save_objects(objects = "sillys_strata_2017_MissLoci", path = "Objects")

## Get number of individuals per silly after removing missing loci individuals
sillys_strata_2017_n[, "Missing"] <- sillys_strata_2017_n[, "Genotyped"] - 
  sapply(paste0(sillys_strata_2017, ".gcl"), function(x) get(x)$n)

### Duplicate
## Check within collections for duplicate individuals.
sillys_strata_2017_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = sillys_strata_2017, loci = loci24, quantile = NULL, minproportion = 0.95)
detach("package:reshape", unload = TRUE)
sillys_strata_2017_DuplicateCheckReportSummary <- sapply(sillys_strata_2017, function(x) sillys_strata_2017_DuplicateCheck95MinProportion[[x]]$report)
sillys_strata_2017_DuplicateCheckReportSummary
save_objects(objects = "sillys_strata_2017_DuplicateCheckReportSummary", path = "Objects")

## Remove duplicate individuals
sillys_strata_2017_RemovedDups <- RemoveDups.GCL(sillys_strata_2017_DuplicateCheck95MinProportion)

### Final
sillys_strata_2017_n[, "Final"] <- sapply(paste0(sillys_strata_2017, ".gcl"), function(x) get(x)$n)
## Get number of individuals per silly after removing duplicate individuals
sillys_strata_2017_n[, "Duplicate"] <- sillys_strata_2017_n[, "Genotyped"] - sillys_strata_2017_n[, "Missing"] - sillys_strata_2017_n[, "Final"]
sillys_strata_2017_n

save_objects(objects = "sillys_strata_2017_n", path = "Objects")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Combine markers
CombineLoci.GCL(sillyvec = sillys_strata_2017, markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = sillys_strata_2017, markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MSA with rubias
#~~~~~~~~~~~~~~~~~~
## Create mixtures
# dir.create("rubias/mixture")
chignik_2017.mix <- create_rubias_mixture(sillyvec = sillys_strata_2017, loci = loci22)
save_objects(objects = "chignik_2017.mix", path = "rubias/mixture")

#~~~~~~~~~~~~~~~~~~
## Run MSA
# dir.create("rubias/output")
chignik_2017.out <- run_rubias_mixture(reference = chignik7_loci22.base,
                                              mixture = chignik_2017.mix, 
                                              group_names = ChignikGroups,
                                              gen_start_col = 5, 
                                              method = "PB", 
                                              reps = 25000, 
                                              burn_in = 5000, 
                                              pb_iter = 100,
                                              sample_int_Pi = 10,
                                              path = "rubias/output")
# str(chignik_2017.out, give.attr = FALSE, max.level = 2)
# save.image("rubias/output/retro_chignik_2017_rubias.RData")

#~~~~~~~~~~~~~~~~~~
## rubais output
chignik_2017.sum <- custom_combine_rubias_output(rubias_output = chignik_2017.out, group_names = ChignikGroups, bias_corr = TRUE)
str(chignik_2017.sum)

chignik_2017.sum %>% 
  select(mixture_collection, repunit, mean) %>% 
  spread(repunit, mean)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Analyze 2018 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Read in 2018 genotypes
sillys <- paste0("SCHIG18")
LOKI2R.GCL(sillyvec = sillys, username = "krshedd", password = password)
save_sillys(sillyvec = sillys, path = "Raw genotypes")
# load_sillys(path = "Raw genotypes", sillyvec = sillys)

sapply(sillys, function(silly) {get(paste0(silly, ".gcl"))$n} )  # 1140/yr
sapply(sillys, function(silly) {table(get(paste0(silly, ".gcl"))$attributes$CAPTURE_DATE, useNA = "always")} )  # 2016 is missing Capture Date


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Read in ASL data
ASL.df <- read_csv(file = "2018Chignik_GeneticsASL2.csv")

ASL.df <- ASL.df %>% 
  unite(euro_age, fw_age, sw_age, sep = ".", remove = FALSE) %>% 
  mutate(euro_age = replace(euro_age, which(euro_age == "NA.NA"), NA)) %>% 
  mutate(euro_age = factor(x = euro_age, levels = sort(unique(euro_age)))) %>% 
  mutate(sample_date = as.Date(ASL.df$sample_date, format = "%m/%d/%Y")) %>% 
  mutate(sex = recode(sex, !!!list("1" = "male", "2" = "female")))

ASL.df %>% 
  group_by(sample_date) %>% 
  summarise(n = n())

# 2018
# att_names_18 <- dimnames(SCHIG18.gcl$attributes)
# SCHIG18.gcl$attributes <- SCHIG18.gcl$attributes %>%
#   left_join(ASL.df, by = c("FK_FISH_ID" = "dna_vial")) %>% 
#   mutate("CAPTURE_DATE" = `sample_date`) %>% 
#   select(att_names)
# dimnames(SCHIG18.gcl$attributes) <- att_names_18

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Stratify mixtures
#~~~~~~~~~~~~~~~~~~
# Function to pool by sample date
PoolCollectionsByDateDF <- function(silly, date.df, loci) {
  sapply(silly, function(mix) {
    mix.dates <- unique(as.Date(get(paste(mix, ".gcl", sep = ''))$attributes$CAPTURE_DATE))
    by(data = date.df, INDICES = date.df$strata, function(x) {
      IDs <- AttributesToIDs.GCL(silly = mix, attribute = "CAPTURE_DATE", matching = mix.dates[mix.dates >= x$begin_date & mix.dates <= x$end_date])
      IDs <- list(na.omit(IDs))
      names(IDs) <- mix
      PoolCollections.GCL(collections = mix, loci = loci, IDs = IDs, newname = paste(mix, as.character(x$strata), sep = "_"))
      list("First Last Fish" = range(as.numeric(unlist(IDs))), "n" = get(paste(mix, "_", as.character(x$strata), ".gcl", sep = ''))$n)
    } )
  }, simplify = FALSE, USE.NAMES = TRUE)
}

date_2018.df <- tibble(strata = 1:6, 
                       begin_date = as.Date(x = c("2018-06-26", "2018-07-02", "2018-07-08", "2018-07-17", "2018-07-22", "2018-07-27")), 
                       end_date = as.Date(x = c("2018-06-27", "2018-07-02", "2018-07-12", "2018-07-17", "2018-07-23", "2018-07-27")))

PoolCollectionsByDateDF(silly = "SCHIG18", date.df = date_2018.df, loci = loci24)

sillys_strata_2018 <- setdiff(
  str_replace(string = objects(pattern = "\\.gcl"), 
              pattern = "\\.gcl", 
              replacement = ''),
  sillys)
save_objects(objects = "sillys_strata_2018", path = "Objects")

# dir.create("Raw genotypes/Strata")
save_sillys(sillyvec = sillys_strata_2018, path = "Raw genotypes/Strata/")

sapply(sillys_strata_2018, function(silly) {table(get(paste0(silly, ".gcl"))$attributes$CAPTURE_DATE, useNA = "always")} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Data QC
sillys_strata_2018_n <- matrix(data = NA, nrow = length(sillys_strata_2018), ncol = 4, 
                               dimnames = list(sillys_strata_2018, c("Genotyped", "Missing", "Duplicate", "Final")))

#~~~~~~~~~~~~~~~~~~
#### Check loci
## Get sample size by locus
original_sillys_strata_2018_n_locus <- SampSizeByLocus.GCL(sillyvec = sillys_strata_2018, loci = loci24)
min(original_sillys_strata_2018_n_locus)  ## 166
round(apply(original_sillys_strata_2018_n_locus, 1, function(locus) {min(locus) / max(locus)}), 2)

original_sillys_strata_2018_percent_locus <- apply(original_sillys_strata_2018_n_locus, 1, function(row) {row / max(row)} )
which(apply(original_sillys_strata_2018_percent_locus, 2, min) < 0.8)  # no re-runs!

# Genotpying percentage by locus and strata
require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(original_sillys_strata_2018_percent_locus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares

#~~~~~~~~~~~~~~~~~~
#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
sillys_strata_2018_n[, "Genotyped"] <- sapply(paste0(sillys_strata_2018, ".gcl"), function(x) get(x)$n)

### Missing
## Remove individuals with >20% missing data
sillys_strata_2018_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = sillys_strata_2018, proportion = 0.8)
save_objects(objects = "sillys_strata_2018_MissLoci", path = "Objects")

## Get number of individuals per silly after removing missing loci individuals
sillys_strata_2018_n[, "Missing"] <- sillys_strata_2018_n[, "Genotyped"] - 
  sapply(paste0(sillys_strata_2018, ".gcl"), function(x) get(x)$n)

### Duplicate
## Check within collections for duplicate individuals.
sillys_strata_2018_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = sillys_strata_2018, loci = loci24, quantile = NULL, minproportion = 0.95)
detach("package:reshape", unload = TRUE)
sillys_strata_2018_DuplicateCheckReportSummary <- sapply(sillys_strata_2018, function(x) sillys_strata_2018_DuplicateCheck95MinProportion[[x]]$report)
sillys_strata_2018_DuplicateCheckReportSummary
save_objects(objects = "sillys_strata_2018_DuplicateCheckReportSummary", path = "Objects")

## Remove duplicate individuals
sillys_strata_2018_RemovedDups <- RemoveDups.GCL(sillys_strata_2018_DuplicateCheck95MinProportion)

### Final
sillys_strata_2018_n[, "Final"] <- sapply(paste0(sillys_strata_2018, ".gcl"), function(x) get(x)$n)
## Get number of individuals per silly after removing duplicate individuals
sillys_strata_2018_n[, "Duplicate"] <- sillys_strata_2018_n[, "Genotyped"] - sillys_strata_2018_n[, "Missing"] - sillys_strata_2018_n[, "Final"]
sillys_strata_2018_n

save_objects(objects = "sillys_strata_2018_n", path = "Objects")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Combine markers
CombineLoci.GCL(sillyvec = sillys_strata_2018, markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = sillys_strata_2018, markerset = c("One_GPDH2", "One_GPDH"), delim=".", update = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MSA with rubias
#~~~~~~~~~~~~~~~~~~
## Create mixtures
# dir.create("rubias/mixture")
chignik_2018.mix <- create_rubias_mixture(sillyvec = sillys_strata_2018, loci = loci22)
save_objects(objects = "chignik_2018.mix", path = "rubias/mixture")

#~~~~~~~~~~~~~~~~~~
## Run MSA
# dir.create("rubias/output")
chignik_2018.out <- run_rubias_mixture(reference = chignik7_loci22.base,
                                       mixture = chignik_2018.mix, 
                                       group_names = ChignikGroups,
                                       gen_start_col = 5, 
                                       method = "PB", 
                                       reps = 25000, 
                                       burn_in = 5000, 
                                       pb_iter = 100,
                                       sample_int_Pi = 10,
                                       path = "rubias/output")
# str(chignik_2018.out, give.attr = FALSE, max.level = 2)
# save.image("rubias/output/retro_chignik_2018_rubias.RData")

#~~~~~~~~~~~~~~~~~~
## rubais output
chignik_2018.sum <- custom_combine_rubias_output(rubias_output = chignik_2018.out, group_names = ChignikGroups, bias_corr = TRUE)
str(chignik_2018.sum)

chignik_2018.sum %>% 
  select(mixture_collection, repunit, mean) %>% 
  spread(repunit, mean)

chignik_2018.sum %>% 
  separate(col = mixture_collection, into = c("silly", "strata"), sep = "_") %>% 
  mutate(strata = as.integer(strata)) %>% 
  left_join(date_2018.df, by = "strata")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### BAYES 2018 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set up file structure
dir.create("BAYES")
sapply(c("Baseline", "Control", "Mixture", "Output"), function(new_dir) {dir.create(paste("BAYES", new_dir, sep = "/"))})
file.copy(from = "../2017/BAYES/ChignikPops24Loci.bse", to = "BAYES/Baseline/")

# Mixture format
ChignikMixtureFormat <- CreateMixture.GCL(sillys = sillys_strata_2018[1], loci = loci22, IDs = NULL, mixname = sillys_strata_2018[1],
                                          dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
save_objects(objects = "ChignikMixtureFormat", path = "Objects/")

# Dump all mixtures
invisible(sapply(sillys_strata_2018, function(silly) {
  CreateMixture.GCL(sillys = silly, loci = loci22, IDs = NULL, mixname = silly,
                    dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
} ))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run 1
# Start prior from model
chignik_2018_start_prior_weight <- (mod_avg$chignik_late_prop[which(mod_avg$date == "2018-06-26")] * 122 + mod_avg$chignik_late_prop[which(mod_avg$date == "2018-06-27")] * 67) / 189
chignik_2018_start_prior_weights <- c(1 - chignik_2018_start_prior_weight, chignik_2018_start_prior_weight)
chignik_1_2018_prior <- Prior.GCL(groupvec = Groupvec7, groupweights = chignik_2018_start_prior_weights, minval = 0.01)
save_objects(objects = "chignik_1_2018_prior", path = "Objects/")

# Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci22, mixname = "SCHIG18_1", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = chignik_1_2018_prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

# Summarize output
dir.create("BAYES/Output/SCHIG18_1")
SCHIG18_1_estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG18_1",
                                                          prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
# dir.create("Estimates objects")
save_objects(objects = "SCHIG18_1_estimates", path = "Estimates objects/")

# Create function to view traceplot
trace_plot.f <- function(estimate_obj) {
  par(mfrow = c(2, 1), mar = c(4.1, 4.1, 2.1, 2.1))
  plot(estimate_obj[[2]][[1]][seq(from = 1, to = 100000, by = 10), 1], type = "l", ylim = c(0, 1), ylab = "", main = "Black Lake", xlab = "")
  abline(v = seq(from = 0, to = 10000, by = 2000))
  par(mar = c(4.1, 4.1, 2.1, 2.1))
  plot(estimate_obj[[2]][[1]][seq(from = 1, to = 100000, by = 10), 2], type = "l", ylim = c(0, 1), ylab = "", main = "Chignik Lake", xlab = "Repetitions", cex.lab = 1.2)
  abline(v = seq(from = 0, to = 10000, by = 2000))
  mtext(text = "Posterior", side = 2, outer = TRUE, line = -1, cex = 1.2)
  
  print(estimate_obj[[1]][[1]][, "GR"])
}

trace_plot.f(SCHIG18_1_estimates)

# Prior for next round
chignik_2_2018_prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG18_1_estimates[[1]][[1]][, 1], minval = 0.01)
save_objects(objects = "chignik_2_2018_prior", path = "Objects/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run 2
# Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci22, mixname = "SCHIG18_2", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = chignik_2_2018_prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

# Summarize output
dir.create("BAYES/Output/SCHIG18_2")
SCHIG18_2_estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG18_2",
                                                    prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
# dir.create("Estimates objects")
save_objects(objects = "SCHIG18_2_estimates", path = "Estimates objects/")

trace_plot.f(SCHIG18_2_estimates)

# Prior for next round
chignik_3_2018_prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG18_2_estimates[[1]][[1]][, 1], minval = 0.01)
save_objects(objects = "chignik_3_2018_prior", path = "Objects/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run 3
# Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci22, mixname = "SCHIG18_3", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = chignik_3_2018_prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

# Summarize output
dir.create("BAYES/Output/SCHIG18_3")
SCHIG18_3_estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG18_3",
                                                    prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
# dir.create("Estimates objects")
save_objects(objects = "SCHIG18_3_estimates", path = "Estimates objects/")

trace_plot.f(SCHIG18_3_estimates)

# Prior for next round
chignik_4_2018_prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG18_3_estimates[[1]][[1]][, 1], minval = 0.01)
save_objects(objects = "chignik_4_2018_prior", path = "Objects/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run 4
# Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci22, mixname = "SCHIG18_4", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = chignik_4_2018_prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

# Summarize output
dir.create("BAYES/Output/SCHIG18_4")
SCHIG18_4_estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG18_4",
                                                    prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
# dir.create("Estimates objects")
save_objects(objects = "SCHIG18_4_estimates", path = "Estimates objects/")

trace_plot.f(SCHIG18_4_estimates)

# Prior for next round
chignik_5_2018_prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG18_4_estimates[[1]][[1]][, 1], minval = 0.01)
save_objects(objects = "chignik_5_2018_prior", path = "Objects/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run 5
# Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci22, mixname = "SCHIG18_5", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = chignik_5_2018_prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

# Summarize output
dir.create("BAYES/Output/SCHIG18_5")
SCHIG18_5_estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG18_5",
                                                    prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
# dir.create("Estimates objects")
save_objects(objects = "SCHIG18_5_estimates", path = "Estimates objects/")

trace_plot.f(SCHIG18_5_estimates)

# Prior for next round
chignik_6_2018_prior <- Prior.GCL(groupvec = Groupvec7, groupweights = SCHIG18_5_estimates[[1]][[1]][, 1], minval = 0.01)
save_objects(objects = "chignik_6_2018_prior", path = "Objects/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run 6
# Dump Control Files
CreateControlFile.GCL(sillyvec = Chignik7Populations, loci = loci22, mixname = "SCHIG18_6", 
                      basename = "ChignikPops24Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = Groupvec7, priorvec = chignik_6_2018_prior, initmat = Inits, 
                      dir = "BAYES/Control", seeds = WASSIPSockeyeSeeds, thin = c(1,1,100),
                      mixfortran = ChignikMixtureFormat, basefortran = Chignik24BaselineFormat, switches = "F T F T T T F") 

# Summarize output
dir.create("BAYES/Output/SCHIG18_6")
SCHIG18_6_estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = "SCHIG18_6",
                                                    prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)  # Yes, I want the Posterior so I can look at trace plot.
# dir.create("Estimates objects")
save_objects(objects = "SCHIG18_6_estimates", path = "Estimates objects/")

trace_plot.f(SCHIG18_6_estimates)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## All run output
SCHIG18_estimates <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = ChignikGroups, maindir = "BAYES/Output", mixvec = sillys_strata_2018,
                                                  prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)
save_objects(objects = "SCHIG18_estimates", path = "Estimates objects/")

chignik_2018.bayes <- bind_rows(sapply(sillys_strata_2018, function(mix) {
  as.tibble(cbind(SCHIG18_estimates[[mix]], "repunit" = ChignikGroups, silly = mix))
}, simplify = FALSE))

chignik_2018_dates.bayes <- chignik_2018.bayes %>% 
  select(silly, repunit, mean, sd, median, `5%`, `95%`, `P=0`, GR) %>% 
  separate(col = silly, into = c("silly", "strata"), sep = "_") %>% 
  mutate(strata = as.integer(strata)) %>% 
  mutate(repunit = factor(x = repunit, levels = ChignikGroups)) %>% 
  mutate(mean = as.double(mean)) %>% 
  mutate(sd = as.double(sd)) %>% 
  mutate(median = as.double(median)) %>% 
  mutate(`5%` = as.double(`5%`)) %>% 
  mutate(`95%` = as.double(`95%`)) %>% 
  mutate(`P=0` = as.double(`P=0`)) %>% 
  mutate(GR = as.double(GR)) %>% 
  left_join(date_avg, by = "strata") %>% 
  mutate(date = as.Date(avg_date, origin = "2017-12-31")) %>% 
  left_join(date_2018.df, by = "strata")

save_objects(objects = "chignik_2018_dates.bayes", path = "Estimates objects/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Average model ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Read in average transition curve for 2018
mod_avg <- read_csv(file = "2018 average GSI model.csv") %>% 
  mutate(date = as.Date(date, format = "%d-%B"))

# Quick plot
mod_avg %>% 
  ggplot(aes(x = date, y = chignik_late_prop)) +
  geom_point()

## Get all attributes to calculate average date for each strata
SCHIG18_attributes <- as.tibble(bind_rows(lapply(sillys_strata_2018, function(silly) {
  cbind(get(paste0(silly, ".gcl"))$attributes, "new_silly" = silly)}
  )))

## Calculate average date per strata
date_avg <- SCHIG18_attributes %>% 
  separate(col = new_silly, into = c("new_silly", "strata"), sep = "_") %>% 
  mutate(strata = as.integer(strata)) %>% 
  mutate(julian = yday(CAPTURE_DATE)) %>% 
  group_by(strata) %>%
  summarize(avg_date = mean(julian))
  
# Join average date with GSI results
chignik_2018_dates.rubias <- chignik_2018.sum %>% 
  separate(col = mixture_collection, into = c("silly", "strata"), sep = "_") %>% 
  mutate(strata = as.integer(strata)) %>% 
  left_join(date_avg, by = "strata") %>% 
  mutate(date = as.Date(avg_date, origin = "2017-12-31")) %>% 
  left_join(date_2018.df, by = "strata")

save_objects(objects = "chignik_2018_dates.rubias", path = "Estimates objects/")

# Plot GSI with avg transition
mod_avg %>% 
  filter(date > "2018-06-15") %>% 
  mutate(black_early_prop = 1 - chignik_late_prop) %>% 
  ggplot(aes(x = date, y = black_early_prop * 100)) +
  geom_line(lwd = 3) +
  geom_point(data = filter(chignik_2018_dates.rubias, repunit == "Black Lake"), aes(x = date, y = mean * 100), colour = "red", cex = 5) +
  geom_errorbar(data = filter(chignik_2018_dates.rubias, repunit == "Black Lake"), aes(x = date, y = mean * 100, ymin = `5%` * 100, ymax = `95%` * 100), colour = "red", lwd = 1.5, width = 3) +
  geom_point(data = filter(chignik_2018_dates.bayes, repunit == "Black Lake"), aes(x = date, y = mean * 100), colour = "blue", cex = 5) +
  geom_errorbar(data = filter(chignik_2018_dates.bayes, repunit == "Black Lake"), aes(x = date, y = mean * 100, ymin = `5%` * 100, ymax = `95%` * 100), colour = "blue", lwd = 1.5, width = 3) +
  xlab("Sample Date") +
  ylab("Black Lake (Early Run) % of Sample") +
  annotate("text", x = as.Date("2018-06-25"), y = c(50, 40), label = c("BAYES", "rubias"), colour = c("blue", "red"), cex = 8) +
  ggtitle("2018 Chignik Sockeye GSI with expected transition")


# Write tables
# dir.create("Estimates tables/")
write_csv(x = chignik_2018_dates.rubias, path = "Estimates tables/Chignik 2018 rubias Estimates.csv")
write_csv(x = chignik_2018_dates.bayes, path = "Estimates tables/Chignik 2018 BAYES Estimates.csv")


save.image("2018ChignikInseason_rubias.RData")
