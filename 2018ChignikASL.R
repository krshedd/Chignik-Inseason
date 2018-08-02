setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2018/Mixtures/2018/")
library(tidyverse)
ASL.df <- read_csv(file = "2018Chignik_GeneticsASL2.csv")

ASL.df <- ASL.df %>% 
  unite(euro_age, fw_age, sw_age, sep = ".", remove = FALSE) %>% 
  mutate(euro_age = replace(euro_age, which(euro_age == "NA.NA"), NA)) %>% 
  mutate(euro_age = factor(x = euro_age, levels = sort(unique(euro_age)))) %>% 
  mutate(sample_date = as.Date(ASL.df$sample_date, format = "%m/%d/%Y")) %>% 
  mutate(sex = recode(sex, !!!list("1" = "male", "2" = "female")))

ages <- sort(unique(ASL.df$euro_age))
ages_fw <- sort(unique(ASL.df$fw_age))
ages_sw <- sort(unique(ASL.df$sw_age))

ASL.df %>% 
  mutate(euro_age = factor(x = euro_age, levels = ages)) %>% 
  group_by(sample_date, euro_age) %>% 
  summarise(age_count = n()) %>% 
  mutate(age_comp = age_count / sum(age_count))

# Fresh water age comp over time
ASL.df %>% 
  mutate(fw_age = factor(x = fw_age, levels = rev(ages_fw))) %>% 
  drop_na(fw_age) %>% 
  group_by(sample_date, fw_age) %>% 
  summarise(age_count = n()) %>% 
  mutate(age_comp = age_count / sum(age_count)) %>% 
  ggplot(aes(x = sample_date, y = age_comp, fill = fw_age)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Sample Date") +
  ylab("Freshwater Age Proportion") +
  ggtitle("Chignik Sockeye - Freshwater Age Composition")

# Salt water age comp over time
ASL.df %>% 
  mutate(sw_age = factor(x = sw_age, levels = rev(ages_sw))) %>% 
  drop_na(sw_age) %>% 
  group_by(sample_date, sw_age) %>% 
  summarise(age_count = n()) %>% 
  mutate(age_comp = age_count / sum(age_count)) %>% 
  ggplot(aes(x = sample_date, y = age_comp, fill = sw_age)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Sample Date") +
  ylab("Saltwater Age Proportion") +
  ggtitle("Chignik Sockeye - Saltwater Age Composition")

# European age comp over time
ASL.df %>% 
  mutate(euro_age = factor(x = euro_age, levels = rev(ages))) %>% 
  drop_na(euro_age) %>% 
  group_by(sample_date, euro_age) %>% 
  summarise(age_count = n()) %>% 
  mutate(age_comp = age_count / sum(age_count)) %>% 
  ggplot(aes(x = sample_date, y = age_comp, fill = euro_age)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Sample Date") +
  ylab("European Age Proportion") +
  ggtitle("Chignik Sockeye - European Age Composition")

# Sex comp over time
ASL.df %>% 
  drop_na(sex) %>% 
  group_by(sample_date, sex) %>% 
  summarise(sex_count = n()) %>% 
  mutate(sex_comp = sex_count / sum(sex_count)) %>% 
  ggplot(aes(x = sample_date, y = sex_comp, fill = sex)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Sample Date") +
  ylab("Sex Proportion") +
  ggtitle("Chignik Sockeye - Sex Composition")
