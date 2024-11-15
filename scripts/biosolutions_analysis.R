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

pgp_data <- read_delim("../data/sarrislab_in_planta_experiments.xlsx - template_pgp.tsv", delim="\t")
nacl_data <- read_delim("../data/sarrislab_in_planta_experiments.xlsx - template_nacl.tsv", delim="\t")
water_deficit_data <- read_delim("../data/sarrislab_in_planta_experiments.xlsx - template_water_deficit.tsv", delim="\t")

microbes <- read_delim("../data/sarrislab_in_planta_experiments.xlsx - template_microbes.tsv", delim="\t")

all_experiments$microbe_id <- gsub('"', '', all_experiments$microbe_id)

# ids for report
ids <- c("378", "295", "247", "620", "614", "253", "305", "323", "345", "094", 
         "248", "630", "175", "204", "297", "234", "003", "018", "022", "033", 
         "086", "097", "129", "139", "144", "166", "184", "186", "190", "197", 
         "200", "203", "210", "227", "243", "263", "273", "288", "289", "309", 
         "324", "343", "347", "348", "356", "550", "552", "558", "574", "605", 
         "606", "609", "621", "623", "624", "628", "740", "742", "1020", "1060")

## Combine all data from experiments
all_experiments <- rbind(pgp_data,nacl_data,water_deficit_data)


microbes_conditions_summary <- all_experiments |>
    distinct(microbe_id,condition) |> 
    group_by(microbe_id) |> 
    summarise(n_conditions=n(),
              conditions=str_c(condition, collapse = ","))

batch_summary <- all_experiments |> 
    distinct(batch_id, microbe_id,condition)

## microbes summary and next steps
# ids from report not in experiments data
ids[!(ids %in% unique(all_experiments$microbe_id))] 

# microbes from the updated list of microbes not in ids from report
microbes_not_in_biosolutions <- microbes[!(microbes$microbe_id %in% ids),]

