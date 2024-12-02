# Author: Madison Coots
# Date: October 31, 2023
# ----------------------
# This file pulls all of the requisite NHANES data tables (from 2011-2018)
# and performs a join for the analysis. The following tables are downloaded:
# 1. Demographics (DEMO)
# 2. Reproductive Health (RHQ)
#
# The full data object is written to the path specified in write_path.
# Coots et al. (2023)

library(haven)
library(tidyverse)
library(janitor)
library(stringr)
library(here)
library(rstudioapi)

working_directory <- getwd()
if (basename(working_directory) == "code") {
  # being called from master script
  directory_path <- here(working_directory, 'bc/')
} else {
  # being run on its own
  directory_path <- dirname(getActiveDocumentContext()$path)
}

data_object_write_path <- directory_path %>% str_remove("code/bc") %>% str_c("data/bc/processed/")
table_write_path <- directory_path %>% str_remove("code/bc") %>% str_c("data/bc/raw/")

# 2011-2012
download.file("https://wwwn.cdc.gov/nchs/nhanes/2011-2012/DEMO_G.XPT", demo <- tempfile(), mode="wb")
download.file("https://wwwn.cdc.gov/nchs/nhanes/2011-2012/RHQ_G.XPT", rhq <- tempfile(), mode="wb")

raw_demographics_11_12 <- foreign::read.xport(demo) %>%
  clean_names()
raw_repro_11_12 <- foreign::read.xport(rhq) %>%
  clean_names()


# 2013-2014
download.file("https://wwwn.cdc.gov/nchs/nhanes/2013-2014/DEMO_H.XPT", demo <- tempfile(), mode="wb")
download.file("https://wwwn.cdc.gov/nchs/nhanes/2013-2014/RHQ_H.XPT", rhq <- tempfile(), mode="wb")

raw_demographics_13_14 <- foreign::read.xport(demo) %>%
  clean_names()
raw_repro_13_14 <- foreign::read.xport(rhq) %>%
  clean_names()


# 2015-2016
download.file("https://wwwn.cdc.gov/nchs/nhanes/2015-2016/DEMO_I.XPT", demo <- tempfile(), mode="wb")
download.file("https://wwwn.cdc.gov/nchs/nhanes/2015-2016/RHQ_I.XPT", rhq <- tempfile(), mode="wb")

raw_demographics_15_16 <- foreign::read.xport(demo) %>%
  clean_names()
raw_repro_15_16 <- foreign::read.xport(rhq) %>%
  clean_names()

# 2017-2018
download.file("https://wwwn.cdc.gov/nchs/nhanes/2017-2018/DEMO_J.XPT", demo <- tempfile(), mode="wb")
download.file("https://wwwn.cdc.gov/nchs/nhanes/2017-2018/RHQ_J.XPT", rhq <- tempfile(), mode="wb")

raw_demographics_17_18 <- foreign::read.xport(demo) %>%
  clean_names()
raw_repro_17_18 <- foreign::read.xport(rhq) %>%
  clean_names()

# ------------------- Combine -------------------

# Demographics data
raw_demographics_all <- bind_rows(
  raw_demographics_11_12,
  raw_demographics_13_14,
  raw_demographics_15_16,
  raw_demographics_17_18
)

# Reproductive health data
raw_repro_all <- bind_rows(
  raw_repro_11_12,
  raw_repro_13_14,
  raw_repro_15_16,
  raw_repro_17_18
)

write_csv(raw_demographics_all, file = paste(table_write_path, "raw_demographics_all.csv", sep = ""))
write_csv(raw_repro_all, file = paste(table_write_path, "raw_repro_all.csv", sep = ""))
 
p_rel_bc_white = 0.12
p_rel_bc_black = 0.09
p_rel_bc_asian = 0.08
p_rel_bc_hispanic = 0.09

# rate_white = 133.7/100000
# rate_black = 127.8/100000
# rate_asian = 101.3/100000
# rate_hispanic = 99.2/100000

set.seed(1)
data <- raw_demographics_all %>%
  inner_join(raw_repro_all, by = c("seqn")) %>% # This table only has data for females
  filter(ridageyr >= 35,
         rhq010 != ".",
         rhq010 != 777,
         rhq010 != 999,
         rhd180 != ".",
         rhd180 != 777,
         rhd180 != 999,
         dmdborn4 != 77,
         ridreth3 != 7) %>%
  mutate(
    US_born = case_when(dmdborn4 == 1 ~ TRUE,
                        dmdborn4 == 2 ~ FALSE),
    race_str = case_when(ridreth3 == 1 | ridreth3 == 2 ~ "Hispanic",
                     ridreth3 == 3 ~ "White",
                     ridreth3 == 4 ~ "Black",
                     ridreth3 == 6 ~ "Asian"),
    race_code = case_when(race_str == "White" ~ 1,
                          race_str == "Black" ~ 2,
                          race_str == "Hispanic" & US_born ~ 3,
                          race_str == "Hispanic" & !US_born ~ 5,
                          race_str == "Asian" ~ 11),
    projection_age = ridageyr + 5, # For 5-year risk
    # Adjusting weights for 8 years of data
    wtmec8yr = wtmec2yr/4,
    N_Biop = 99,
    HypPlas = 99) %>%
  # mutate(
  #   rel_1 = case_when(
  #     race_str == "Black" ~ rbinom(n = n(), size = 1, prob = p_bc_black),
  #     race_str == "White" ~ rbinom(n = n(), size = 1, prob = p_bc_white),
  #     race_str == "Hispanic" ~ rbinom(n = n(), size = 1, prob = p_bc_hispanic),
  #     race_str == "Asian" ~ rbinom(n = n(), size = 1, prob = p_bc_asian),
  #     TRUE ~ NA_integer_), # for any race not specified
  #   rel_2 = case_when(
  #     race_str == "Black" ~ rbinom(n = n(), size = 1, prob = p_bc_black),
  #     race_str == "White" ~ rbinom(n = n(), size = 1, prob = p_bc_white),
  #     race_str == "Hispanic" ~ rbinom(n = n(), size = 1, prob = p_bc_hispanic),
  #     race_str == "Asian" ~ rbinom(n = n(), size = 1, prob = p_bc_asian),
  #     TRUE ~ NA_integer_),
  #   N_Rels = rel_1 + rel_2) %>%
  # mutate(
  #   N_Rels = case_when(
  #     race_str == "Black" ~ rpois(n = n(), lambda = n_rel_samples*rate_black),
  #     race_str == "White" ~ rpois(n = n(), lambda = n_rel_samples*rate_white),
  #     race_str == "Hispanic" ~ rpois(n = n(), lambda = n_rel_samples*rate_hispanic),
  #     race_str == "Asian" ~ rpois(n = n(), lambda = n_rel_samples*rate_asian))) %>%
  mutate(
    N_Rels = case_when(
      race_str == "Black" ~ rbinom(n = n(), size = 1, prob = p_rel_bc_black),
      race_str == "White" ~ rbinom(n = n(), size = 1, prob = p_rel_bc_white),
      race_str == "Hispanic" ~ rbinom(n = n(), size = 1, prob = p_rel_bc_hispanic),
      race_str == "Asian" ~ rbinom(n = n(), size = 1, prob = p_rel_bc_asian))) %>%
  select(ID = seqn,
         wtmec8yr,
         T1 = ridageyr,
         T2 = projection_age,
         race_str,
         Race = race_code,
         US_born,
         AgeMen = rhq010,
         Age1st = rhd180,
         N_Biop,
         HypPlas,
         N_Rels)

saveRDS(data, file = paste(data_object_write_path, "bc_data.rds", sep = ""))
