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

library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(rstatix)

# load data

pgp_data <- read_delim("../data/sarrislab_biosolutions_project - pgp_experiment.tsv",
                       na=c("","NA"),
                       delim="\t")
nacl_data <- read_delim("../data/sarrislab_biosolutions_project - nacl_experiment.tsv",
                       na=c("","NA"),
                       delim="\t")
water_deficit_data <- read_delim("../data/sarrislab_biosolutions_project - water_deficit_experiment.tsv",
                       na=c("","NA"),
                       delim="\t")

synthetic_community <- read_delim("../data/sarrislab_biosolutions_project - synthetic_community_experiment.tsv",
                       na=c("","NA"),
                       delim="\t")

## add the final file with correct data for microbes
#microbes <- read_delim("../data/sarrislab_biosolutions_project.xlsx - microbes.tsv", delim="\t")

plant_batches <- read_delim("../data/sarrislab_biosolutions_project - plant_batches.tsv", delim="\t")


# microbe ids from the first report
ids <- c("378", "295", "247", "620", "614", "253", "305", "323", "345", "094", 
         "248", "630", "175", "204", "297", "234", "003", "018", "022", "033", 
         "086", "097", "129", "139", "144", "166", "184", "186", "190", "197", 
         "200", "203", "210", "227", "243", "263", "273", "288", "289", "309", 
         "324", "343", "347", "348", "356", "550", "552", "558", "574", "605", 
         "606", "609", "621", "623", "624", "628", "740", "742", "1020", "1060")

## Combine all data from experiments
all_experiments <- bind_rows(pgp_data,nacl_data,water_deficit_data,synthetic_community) |>
    mutate(var_rosette_healthy_leaves= (var_rosette_leaves - var_dry_rosette_leaves)/var_rosette_leaves)

all_experiments$microbe_id <- gsub('"', '', all_experiments$microbe_id)

# filter for microbes that there are sequences
exclude <- c("247", "248", "347","558","606","620","636","740", "742", "1020", "1060")

all_experiments <- all_experiments |>
    filter(!(microbe_id %in% exclude))

write_delim(all_experiments,"../results/all_experiments_data.tsv",delim="\t")

print(unique(all_experiments$condition))

microbes_conditions_summary <- all_experiments |>
    distinct(microbe_id,condition) |> 
    group_by(microbe_id) |> 
    summarise(n_conditions=n(),
              conditions_tested=str_c(condition, collapse = ","))


write_delim(microbes_conditions_summary,"../results/microbes_conditions_summary.tsv",delim="\t")

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
#microbes_not_in_biosolutions <- microbes[!(microbes$microbe_id %in% ids),]


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

all_experiments_mean_l <- all_experiments_stats |>
    dplyr::select(microbe_id,
                  batch_id,
                  condition,
                  contains("mean")) |>
    pivot_longer(
        cols = contains("mean"),  # Select columns with "mean" in their names
        names_to = "variable",    # Create a column for variable names
        values_to = "value")       # Create a column for the value


# Calculate percentage change
all_experiments_norm <- all_experiments_mean_l |>
    group_by(batch_id,variable) |>
    mutate(control_value = value[microbe_id == "Control"], #Extract control value
         percent_change = ((value - control_value) / control_value) * 100,
         difference_value = value - control_value) |>
  ungroup() |>
  select(-control_value) # Optional: Remove intermediate control_value column

all_experiments_percent_change <- all_experiments_norm |>
    select(-c(value,difference_value)) |>  # Remove the original value column
    mutate(percent_change=round(percent_change,2)) |>
    pivot_wider(
                names_from = variable,
                values_from = percent_change
                ) |>
    filter(!(microbe_id %in% c("Control","Control+", "E. coli")))

write_delim(all_experiments_percent_change,"../results/all_experiments_percent_change.tsv",delim="\t")

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
    filter(!is.na(microbe_id))


# nest the data based on batch ids, i.e individual experiments
nested_data <- experiments_fi |>
    ungroup() |>
    group_by(batch_id) |>
    nest() |>
    mutate(
    # Identify all "var" columns
    vars = map(data, ~ select(.x, starts_with("var")) %>% colnames()))


variables <- select(experiments_fi, starts_with("var")) %>% colnames()

# ------------------------ Non parametric------------------------- #

print("Non parametric test")
# initiate the empty lists for the no parametric tests
all_stats_np <- list()
pairwise_stats_np <- list()

# iterate through all variables of the experiments
for (i in seq_along(variables)) {
    print(reformulate("microbe_id", variables[i]))
# Perform tests the var i
    stats_nopara <- nested_data |>
        mutate(
               kruskal=map(
                      data, ~ kruskal_test(reformulate("microbe_id", variables[i]),
                                              data = .x) %>% mutate(variable = variables[i])),
               dunn=map(
                         data, ~ dunn_test(data=.x, formula=reformulate("microbe_id",variables[i]),
                                                          p.adjust.method = "BH", 
                                                          detailed = FALSE) %>% 
                         mutate(variable = variables[i]))
        )
# Unnest results for long-format output
    stats_nopara_l <- stats_nopara |>
        select(batch_id,kruskal) |>
        unnest(kruskal, names_sep = "_")
# unnest the dunn test results
    dunn_l <- stats_nopara |>
        select(batch_id,dunn) |>
        unnest(dunn, names_sep = "_")
# and the results to list, one dataframe for each variable
    pairwise_stats_np[[i]] <- dunn_l
    all_stats_np[[i]] <- stats_nopara_l
}

# -------- combine the results ------- #
all_stats_np_results <- bind_rows(all_stats_np) |> distinct()

stats_results_kruskal <- all_stats_np_results |>
    dplyr::select(batch_id, starts_with("kruskal")) |>
    distinct()

write_delim(stats_results_kruskal,"../results/stats_results_kruskal.tsv",delim="\t")


####################### Post hoc test against the control ######################

# bind together all results from all variables and batches
all_pairwise_dunn_results <- bind_rows(pairwise_stats_np) |> distinct()

write_delim(all_pairwise_dunn_results,"../results/all_pairwise_dunn_results.tsv",delim="\t")

# filter only the pairwise comparisons with the control values
# we need all microbe id in group1 when compared with dunn_group2=="Control"
control_pairwise_dunn_sig <- all_pairwise_dunn_results |>
    filter(dunn_group2=="Control") |>
    filter(dunn_p.adj.signif!="ns") 

write_delim(control_pairwise_dunn_sig,"../results/control_pairwise_dunn_sig.tsv",delim="\t")

#filter the pairwise only if they are significant and higher than the control of
#each batch and variable

all_experiments_dunn_sig <- all_experiments_norm |> 
    filter(if_else(variable=="var_dry_rosette_leaves_mean", percent_change < 0, percent_change > 0)) |>
    mutate(variable=gsub("_mean","",variable)) |>
    inner_join(control_pairwise_dunn_sig,
               by=c("variable"="dunn_variable",
                    "microbe_id"="dunn_group1",
                    "batch_id"="batch_id")
    )

write_delim(all_experiments_dunn_sig,"../results/all_experiments_dunn_sig.tsv",delim="\t")



# ------------------------ parametric------------------------- #
# perform parametric Anova and tukey post hoc for all variables

print("Parametric test")
all_stats <- list()
pairwise_stats <- list()

for (i in seq_along(variables)) {
    print(reformulate("microbe_id", variables[i]))
# Perform tests the var i
    asasa <- nested_data |>
        mutate(
               anova=map(
                      data, ~ broom::tidy(aov(reformulate("microbe_id", variables[i]),
                                              data = .x)) %>% mutate(variable = variables[i])),
               tukey=map(
                         data, ~ broom::tidy(TukeyHSD(aov(reformulate("microbe_id", variables[i]),
                                                          data=.x))) %>% mutate(variable = variables[i]))
        )
# Unnest results for long-format output
    asasa_long <- asasa |>
        select(batch_id,anova) |>
        unnest(anova, names_sep = "_") 
    tukey_l <- asasa |>
        select(batch_id,tukey) |>
        unnest(tukey, names_sep = "_")
    pairwise_stats[[i]] <- tukey_l
    all_stats[[i]] <- asasa_long
}

all_stats_results <- bind_rows(all_stats) |> distinct()

stats_results_anova <- all_stats_results |>
    dplyr::select(batch_id, starts_with("anova")) |>
    distinct()

write_delim(stats_results_anova,"../results/stats_results_anova.tsv",delim="\t")

####################### Post hoc test against the control ######################

all_pairwise_results <- bind_rows(pairwise_stats) |> distinct()

write_delim(all_pairwise_results,"../results/all_pairwise_results.tsv",delim="\t")

control_pairwise_sig <- all_pairwise_results |>
    filter(grepl("Control-",tukey_contrast)) |>
    filter(tukey_adj.p.value < 0.1) |>
    mutate(microbe_id=gsub("^.*-","",tukey_contrast)) |>
    filter(!grepl("Control",microbe_id))

write_delim(control_pairwise_sig,"../results/control_pairwise_sig_microbes.tsv",delim="\t")

#filter the pairwise only if they are significant and higher than the control of
#each batch and variable

all_experiments_tukey_sig <- all_experiments_norm |> 
    filter(if_else(variable=="var_dry_rosette_leaves_mean", percent_change < 0, percent_change > 0)) |>
    mutate(variable=gsub("_mean","",variable)) |>
    inner_join(control_pairwise_sig,
               by=c("variable"="tukey_variable",
                    "microbe_id"="microbe_id",
                    "batch_id"="batch_id")
    )

write_delim(all_experiments_tukey_sig,"../results/all_experiments_tukey_sig.tsv",delim="\t")

