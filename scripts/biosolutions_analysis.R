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
    dplyr::select(microbe_id,batch_id,starts_with("var")) |>
    distinct() |>
    filter(!is.na(microbe_id), batch_id!=15)


# nest the data based on batch ids, i.e individual experiments
nested_data <- experiments_fi |>
    ungroup() |>
    group_by(batch_id) |>
    nest() |>
    mutate(
    # Identify all "var" columns
    vars = map(data, ~ select(.x, starts_with("var")) %>% colnames()))

# perform the statistics, Anova and Kruskal Wallis, for all variables
variables <- select(experiments_fi, starts_with("var")) %>% colnames()

all_stats <- list()
pairwise_stats <- list()

for (i in seq_along(variables)) {
    
    print(reformulate("microbe_id", variables[i]))
# Perform tests the var i
    asasa <- nested_data |>
        mutate(anova=map(
                      data, ~ broom::tidy(aov(reformulate("microbe_id", variables[i]),
                                              data = .x)) %>% mutate(variable = variables[i])),
               kruskal=map(
                      data, ~ broom::tidy(kruskal.test(reformulate("microbe_id", variables[i]),
                                              data = .x)) %>% mutate(variable = variables[i])),
               tukey=map(
                         data, ~ broom::tidy(TukeyHSD(aov(reformulate("microbe_id", variables[i]),
                                                          data=.x))) %>% mutate(variable = variables[i]))
        )
# Unnest results for long-format output
    asasa_long <- asasa |>
        select(batch_id,anova,kruskal) |>
        unnest(anova, names_sep = "_") |>
        unnest(kruskal, names_sep = "_")
    
    tukey_l <- asasa |>
        select(batch_id,tukey) |>
        unnest(tukey, names_sep = "_")
    pairwise_stats[[i]] <- tukey_l

    all_stats[[i]] <- asasa_long
}

all_stats_results <- bind_rows(all_stats) |> distinct()


stats_results_kruskal <- all_stats_results |>
    dplyr::select(batch_id, starts_with("kruskal")) |>
    distinct()

write_delim(stats_results_kruskal,"../results/stats_results_kruskal.tsv",delim="\t")

stats_results_anova <- all_stats_results |>
    dplyr::select(batch_id, starts_with("anova")) |>
    distinct()

write_delim(stats_results_anova,"../results/stats_results_anova.tsv",delim="\t")


####################### Post hoc test against the control ######################

all_pairwise_results <- bind_rows(pairwise_stats) |> distinct()

write_delim(all_pairwise_results,"../results/all_pairwise_results.tsv",delim="\t")

control_pairwise_sig <- all_pairwise_results |>
    filter(grepl("Control",tukey_contrast)) |>
    filter(tukey_adj.p.value < 0.1) |>
    mutate(microbe_id=gsub("^.*-","",tukey_contrast)) |>
    filter(!grepl("Control",microbe_id))

write_delim(control_pairwise_sig,"../results/control_pairwise_sig_microbes.tsv",delim="\t")
