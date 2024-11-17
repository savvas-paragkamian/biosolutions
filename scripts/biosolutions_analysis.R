#!/usr/bin/Rscript
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

pgp_data <- read_delim("../data/sarrislab_in_planta_experiments.xlsx - pgp_experiment.tsv",
                       na=c("","NA"),
                       delim="\t")
nacl_data <- read_delim("../data/sarrislab_in_planta_experiments.xlsx - nacl_experiment.tsv",
                       na=c("","NA"),
                       delim="\t")
water_deficit_data <- read_delim("../data/sarrislab_in_planta_experiments.xlsx - water_deficit_experiment.tsv",
                       na=c("","NA"),
                       delim="\t")

microbes <- read_delim("../data/sarrislab_in_planta_experiments.xlsx - microbes.tsv", delim="\t")

plant_batches <- read_delim("../data/sarrislab_in_planta_experiments.xlsx - plant_batches.tsv", delim="\t")


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

microbes_conditions_summary <- all_experiments |>
    distinct(microbe_id,condition) |> 
    group_by(microbe_id) |> 
    summarise(n_conditions=n(),
              conditions_tested=str_c(condition, collapse = ","))

batch_summary <- all_experiments |> 
    distinct(batch_id, microbe_id,condition)

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


all_experiments_stats <- all_experiments |>
    group_by(microbe_id,batch_id,condition) |>
    summarise(across(
                     starts_with("var"),
                     list(mean = ~mean(.,na.rm=T),
                          #sd = ~sd(.),
                          sderr = ~sd(.,na.rm=T) / sqrt(n())),
                     .names = "{.col}_{.fn}"),
              var_n_plants=n(),
              .groups="keep")


all_experiments_stats_l <- all_experiments_stats |>
    pivot_longer(cols = starts_with("var"),
    names_to = "variable",
    values_to = "value")

#
# Data
#df <- tibble(
#  id = c("A", "A", "B", "B"),
#  condition = c("control", "treated", "control", "treated"),
#  value = c(10, 15, 8, 12)
#)
#
## Calculate percentage change
#df <- df %>%
#  group_by(id) %>%
#  mutate(control_value = value[condition == "control"], # Extract control value
#         percent_change = ((value - control_value) / control_value) * 100) %>%
#  ungroup() %>%
#  select(-control_value) # Optional: Remove intermediate control_value column
#all_experiments_stats_l_n <- all_experiments_stats_l |>
#    group_by()

## Figures

vars <- unique(all_experiments_stats_l$variable)
vars_mean <- vars[grep("mean",vars)]
vars_err <- vars[grep("sderr",vars)]

for (i in seq_along(vars_mean)){
    print(i)

    fig_bar <- ggplot(all_experiments_stats,
                      aes(x = microbe_id,
                          y = !!sym(vars_mean[i]),
                          fill = condition)) + 
                geom_col(width = 0.6,
                         position = position_dodge(width = 0.82)) +
                geom_errorbar(aes(ymin = !!sym(vars_mean[i]) - !!sym(vars_err[i]),
                                  ymax = !!sym(vars_mean[i]) + !!sym(vars_err[i])),
                              position = position_dodge(width = 0.82),
                              width = 0.1) +  # Error bars
                labs(
                     title = "Mean Values with Standard Error",
                     x = "Microbe ID",
                     y = vars_mean[i]) +
                theme_bw() +
                scale_fill_brewer(palette = "Set3") +
                theme(
                      axis.text.x = element_text(angle = 45, hjust = 1,size = 11)
                )
    
    ggsave(paste0("../figures/plant_growth_promotion_mean_",vars_mean[i],".png"),
           plot=fig_bar, 
           height = 20, 
           width = 50,
           dpi = 300, 
           units="cm",
           device="png")
    
}



## Statistics
