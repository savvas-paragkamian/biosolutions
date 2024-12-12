#!/usr/bin/env Rscript
## Script name: biosolutions_analysis.R
##
## Purpose of script: 
## Data summary and statistics of biosolutions project experiments.
## There are 3 different plant conditions that are inoculated with 
## ~60 microbial isolates to screen their role association with 
## plant growth.
## The plant model is Arabidopsis thaliana.
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-11-14
##
library(tidyverse)
library(readxl)
#library(vegan)

# load data

pgp_data <- read_delim("../data/sarrislab_biosolutions_project.xlsx - pgp_experiment.tsv",
                       na=c("","NA"),
                       delim="\t")
nacl_data <- read_delim("../data/sarrislab_biosolutions_project.xlsx - nacl_experiment.tsv",
                       na=c("","NA"),
                       delim="\t")
water_deficit_data <- read_delim("../data/sarrislab_biosolutions_project.xlsx - water_deficit_experiment.tsv",
                       na=c("","NA"),
                       delim="\t")

microbes <- read_delim("../data/sarrislab_biosolutions_project.xlsx - microbes.tsv", delim="\t")

plant_batches <- read_delim("../data/sarrislab_biosolutions_project.xlsx - plant_batches.tsv", delim="\t")


# microbe ids from the first report
ids <- c("378", "295", "247", "620", "614", "253", "305", "323", "345", "094", 
         "248", "630", "175", "204", "297", "234", "003", "018", "022", "033", 
         "086", "097", "129", "139", "144", "166", "184", "186", "190", "197", 
         "200", "203", "210", "227", "243", "263", "273", "288", "289", "309", 
         "324", "343", "347", "348", "356", "550", "552", "558", "574", "605", 
         "606", "609", "621", "623", "624", "628", "740", "742", "1020", "1060")

## Combine all data from experiments
all_experiments <- bind_rows(pgp_data,nacl_data,water_deficit_data)
all_experiments$microbe_id <- gsub('"', '', all_experiments$microbe_id)

write_delim(all_experiments,"../results/all_experiments_data.tsv",delim="\t")

print(unique(all_experiments$condition))

microbes_conditions_summary <- all_experiments |>
    distinct(microbe_id,condition) |> 
    group_by(microbe_id) |> 
    summarise(n_conditions=n(),
              conditions_tested=str_c(condition, collapse = ","))

# there are some microbes that have been tested twice 
# in the same condition
microbes_batches_summary <- all_experiments |>
    distinct(microbe_id,condition, batch_id) |> 
    group_by(microbe_id,condition) |> 
    summarise(n_batches=n(), .groups="keep") 

## microbes summary and next steps
# ids from report not in experiments data
ids[!(ids %in% unique(all_experiments$microbe_id))] 

# microbes from the updated list of microbes not in ids from report
microbes_not_in_biosolutions <- microbes[!(microbes$microbe_id %in% ids),]


## Normalisation
## how to compare with the control values?

### calculate mean and sd error values
all_experiments_stats <- all_experiments |>
    group_by(microbe_id,batch_id,condition) |>
    summarise(across(
                     starts_with("var"),
                     list(mean = ~mean(.,na.rm=T),
                          #sd = ~sd(.),
                          sderr = ~sd(.,na.rm=T) / sqrt(n())),
                     .names = "{.col}_{.fn}"),
              var_n_plants=n(),
              .groups="keep") |>
    ungroup()

write_delim(all_experiments_stats,"../results/all_experiments_stats.tsv",delim="\t")

########################### plant batches ##########################
##

batch_condition <- plant_batches |>
    filter(action=="condition") |>
    dplyr::select(batch_id,value) |>
    rename("condition"="value")

batch_microbes <- all_experiments |> 
    distinct(batch_id, microbe_id) |>
    filter(!(microbe_id %in% c("Control", "E. coli"))) |>
    group_by(batch_id) |>
    summarise(microbes=n(),
              microbes_id=str_c(microbe_id, collapse = ","))


plant_batches_w <- plant_batches |> 
    filter(!is.na(date)) |>
    pivot_wider(id_cols=batch_id,
                names_from=action,
                values_from=date) 

batch_steps <- data.frame(action=c("plant_seeds","inoculation","wet_mass_measurement"),
                         steps=c("plant_growth","inoculation","measurements"))

plant_batches_l <- plant_batches |>
    filter(!is.na(date)) |>
    dplyr::select(-c(value,Comments)) |>
    group_by(batch_id) |>
    mutate(date_end=if_else(is.na(lead(date)), date, lead(date))) |>
    left_join(batch_steps) |>
    left_join(batch_condition) |>
    left_join(batch_microbes) |>
    ungroup() |>
    mutate(duration=date_end-date) |>
    filter(duration>0) |>
    mutate(midpoint = date + (date_end - date) / 2) |>
    mutate(endpoint = if_else(steps=="measurements",date_end + 3,NA)) 

write_delim(plant_batches_l,"../results/plant_batches_results.tsv",delim="\t")

################################# Statistics ##############################

# load the data and filter NA values
experiments_fi <- all_experiments |>
    dplyr::distinct(microbe_id,batch_id,var_fresh_weight_mg) |>
    na.omit()

# nest the data based on batch ids, i.e individual experiments
nested_data <- experiments_fi %>%
    ungroup() |>
  group_by(batch_id) %>%
  nest()

# perform the statistics, Anova and Kruskal Wallis, for a single variable
experiment_stats <- nested_data |>
    mutate(
        anova = map(data, ~ aov(var_fresh_weight_mg ~ microbe_id, data = .x)),
        summary_anova = map(anova, broom::tidy),
        kruskal = map(data, ~ kruskal.test(var_fresh_weight_mg ~ microbe_id, data = .x)),
        summary_kruskal = map(kruskal,broom::tidy)
    )

# Trasform the summaries to long format with unnest for 
# easier handling
experiment_stats_long <- experiment_stats |>
    select(batch_id,summary_anova,summary_kruskal) |>
    unnest(summary_anova, names_sep = "_") |>
    unnest(summary_kruskal, names_sep = "_")
