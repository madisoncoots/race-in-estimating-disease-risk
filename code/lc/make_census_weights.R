# Author: Madison Coots
# Date: August 22, 2024
# ======================
# Code for constructing survey weights from U.S. census data for lung cancer 
# analysis data. 
# Coots et al. (2024)
# 
# Data accessed from 

library(tidyverse)
library(janitor)

# =============================================================================
census_data_read_path <- here::here(dirname(rstudioapi::getActiveDocumentContext()$path)) %>% 
  str_remove("code/lc") %>%
  str_c("data/lc/raw/")

lc_data_read_path <- "~/Documents/harvard/research/lung_cancer/data/processed/lc_data.rds"

write_path <- here::here(dirname(rstudioapi::getActiveDocumentContext()$path)) %>% 
  str_remove("code/lc") %>%
  str_c("data/lc/processed/")
# =============================================================================

# Part I: Extract relevant census counts 
# =============================================================================

raw_census_data <- read_csv(str_c(census_data_read_path, 'raw_census_data.csv')) %>%
  clean_names()

raw_lc_data <- readRDS(lc_data_read_path)

cleaned_census_data <- raw_census_data %>%
  filter(month == 6, year == 2024) %>% # Filter to the most recent period in the data (June 2024)
  select(age,
         nhwa_male, nhwa_female,
         nhba_male, nhba_female,
         nhaa_male, nhaa_female,
         nhna_male, nhna_female, # Native Hawaiian and Other Pacific Islander
         h_male, h_female) %>%
  filter(age >= 40, age <= 80, age != 999) %>% # Filter to appropriate age subset for lung cancer model (LCRAT)
  pivot_longer( # Re-shape data into desired format for easy joining
    cols = -age,  # This excludes the age column from the transformation
    names_to = c("race", "sex"),  # Split the column names into 'race' and 'sex'
    names_sep = "_",  # Specify the separator used in your column names
    values_to = "number_of_people"  # The new column name for the values
  ) %>%
  mutate(
    race = case_when(
      str_detect(race, "nhwa") ~ "Non-Hispanic White",
      str_detect(race, "nhba") ~ "Non-Hispanic Black",
      str_detect(race, "nhaa") ~ "Non-Hispanic Asian",
      str_detect(race, "nhna") ~ "Non-Hispanic Pacific Islander",
      str_detect(race, "h") ~ "Hispanic",
      TRUE ~ race
    ),
    sex = case_when(
      sex == "male" ~ "Male",
      sex == "female" ~ "Female",
      TRUE ~ sex
    )
  )

# Part II: Attach weights to data
# =============================================================================
# Some processing to make sure column names are compatible between data sets
cleaned_census_weights <- cleaned_census_data %>%
  mutate(female = if_else(sex == "Female", 1, 0)) %>% 
  select(-sex) %>%
  mutate(race_str = case_when(
    race == "Non-Hispanic Asian" ~ "Asian/PI",
    race == "Non-Hispanic Pacific Islander" ~ "Asian/PI",
    race == "Non-Hispanic Black" ~ "Black",
    race == "Non-Hispanic White" ~ "White",
    race == "Hispanic" ~ "Hispanic",
  )) %>%
  select(-race)

# Calculate proportions in census data
total_population <- sum(cleaned_census_weights$number_of_people)

population_proportions <- cleaned_census_weights %>%
  group_by(age, female, race_str) %>%
  summarise(
    group_population = sum(number_of_people),
    proportion = group_population / total_population,
    .groups = 'drop'
  )

# Calculate study sample proportions
lc_data_proportions <- raw_lc_data %>%
  group_by(age, female, race_str) %>%
  summarise(
    group_count = n(),
    proportion = group_count / nrow(raw_lc_data),
    .groups = 'drop'
  )

# Calculate weights
weights <- left_join(lc_data_proportions, population_proportions, by = c("age", "female", "race_str")) %>%
  mutate(weight = proportion.y / proportion.x) %>%
  select(-proportion.x, -proportion.y, -group_population)

# Join weights back into the study data
lc_data_weights <- raw_lc_data %>%
  left_join(weights, by = c("age", "female", "race_str")) %>%
  select(-group_count) %>%
  select(pid, weight)

saveRDS(lc_data_weights, str_c(write_path, 'lc_data_weights.rds'))

