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


### transform to long format
all_experiments_stats_l <- all_experiments_stats |>
    pivot_longer(cols = starts_with("var"),
    names_to = "variable",
    values_to = "value")


## only controls all data
controls_all <- all_experiments |>
    filter(microbe_id %in% c("Control", "E. coli")) |>
    mutate(controls = paste(microbe_id,batch_id, sep="_"))# ,"E. coli" 

# Data

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
         percent_change = ((value - control_value) / control_value) * 100) |>
  ungroup() |>
  select(-control_value) # Optional: Remove intermediate control_value column

all_experiments_norm_no_control <- all_experiments_norm |>
    select(-value) %>%  # Remove the original value column
    pivot_wider(
                names_from = variable,
                values_from = percent_change
                ) |>
    filter(microbe_id!="Control")


## Figures
#
### for all variables across contitions
vars <- unique(all_experiments_stats_l$variable)
vars_mean <- vars[grep("mean",vars)]
vars_err <- vars[grep("sderr",vars)]

for (i in seq_along(vars_mean)){
    print(i)

    fig_bar_all <- ggplot(all_experiments_stats,
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
                      axis.text.x = element_text(angle = 45,
                                                 hjust = 1,
                                                 size = 11)
                )
    
    ggsave(paste0("../figures/plant_growth_promotion_mean_",vars_mean[i],".png"),
           plot=fig_bar_all, 
           height = 20, 
           width = 50,
           dpi = 300, 
           units="cm",
           device="png")

## percent change

    fig_percent_c <- ggplot(all_experiments_norm_no_control,
                      aes(x = microbe_id,
                          y = !!sym(vars_mean[i]),
                          fill = condition)) + 
                geom_col(width = 0.6,
                         position = position_dodge(width = 0.82)) +
                labs(
                     title = "Percentage Change with mean of Control plants",
                     x = "Microbe ID",
                     y = vars_mean[i]) +
                theme_bw() +
                scale_fill_brewer(palette = "Set3") +
                theme(
                      axis.text.x = element_text(angle = 45,
                                                 hjust = 1,
                                                 size = 11)
                )
    
    ggsave(paste0("../figures/plant_growth_promotion_percent_change_",
                  vars_mean[i],
                  ".png"),
           plot=fig_percent_c, 
           height = 20, 
           width = 50,
           dpi = 300, 
           units="cm",
           device="png")
    
}


### for each batch and each variable

batches <- unique(all_experiments_stats$batch_id)
vars_a <- colnames(all_experiments)[grepl("var", colnames(all_experiments))]

for (b in seq_along(batches)) {

    batch_data <- all_experiments_stats |>
        filter(batch_id==batches[b])

    condition <- unique(batch_data$condition)
    print(batches[b])
    print(condition)
    print("bar plot")
    # bar plot

    for (i in seq_along(vars_mean)){
        print(i)
    
        fig_bar <- ggplot(batch_data,
                          aes(x = microbe_id,
                              y = !!sym(vars_mean[i]))) + 
                    geom_col(width = 0.6,
                             position = "identity")+
                    #position_dodge(width = 0.82)) +
                    geom_errorbar(aes(ymin = !!sym(vars_mean[i]) - !!sym(vars_err[i]),
                                      ymax = !!sym(vars_mean[i]) + !!sym(vars_err[i])),
                                  #position = position_dodge(width = 0.82),
                                  width = 0.1) +  # Error bars
                    labs(
                         title = paste0(condition," batch = ",batches[b]),
                         x = "Microbe ID",
                         y = vars_mean[i]) +
                    theme_bw() +
                    scale_fill_brewer(palette = "Set3") +
                    theme(
                          axis.text.x = element_text(angle = 45,
                                                     hjust = 1,
                                                     size = 11)
                    )
        
        ggsave(paste0("../figures/","batch_",
                      batches[b],
                      "_",condition,"_",
                      vars_mean[i],"_bar.png"),
               plot=fig_bar, 
               height = 20, 
               width = 50,
               dpi = 300, 
               units="cm",
               device="png")
    }
        
    print("box plot")
    # box plot
    batch_data_a <- all_experiments |>
        filter(batch_id==batches[b])


    for (i in seq_along(vars_a)){
        print(i)
    
        fig_bar <- ggplot(batch_data_a,
                          aes(x = microbe_id,
                              y = !!sym(vars_a[i]))) + 
                    geom_boxplot(width = 0.6,
                             position = "identity")+
                    geom_point(
                             position = "identity")+
                    #position_dodge(width = 0.82)) +
                    labs(
                         title = paste0(condition," batch = ",batches[b]),
                         x = "Microbe ID",
                         y = vars_a[i]) +
                    theme_bw() +
                    scale_fill_brewer(palette = "Set3") +
                    theme(
                          axis.text.x = element_text(angle = 45,
                                                     hjust = 1,
                                                     size = 11)
                    )
        
        ggsave(paste0("../figures/","batch_",
                      batches[b],
                      "_",condition,"_",
                      vars_a[i],"_boxplot.png"),
               plot=fig_bar, 
               height = 20, 
               width = 50,
               dpi = 300, 
               units="cm",
               device="png")

    }
}

conditions <- unique(all_experiments$condition)

controls_all_l <- controls_all |>
    pivot_longer(cols = starts_with("var"),
                 names_to = "variable",
                 values_to = "value")

for (i in seq_along(conditions)) {


    controls_all_c <- controls_all_l |>
        filter(condition==conditions[i]) 

    fig_controls <- ggplot(controls_all_c,
                      aes(x = controls,
                          y = value)) + 
                geom_boxplot(width = 0.6,
                         position = "identity")+
                geom_point(
                         position = "identity")+
                #position_dodge(width = 0.82) +
                labs(
                     x = "Microbe ID",
                     y = "value") +
                theme_bw() +
                scale_fill_brewer(palette = "Set3") +
                theme(
                      axis.text.x = element_text(angle = 45,
                                                 hjust = 1,
                                                 size = 11)
                )+
                facet_wrap(~variable, ncol=1, nrow=9,scales="free" )
    
    ggsave(paste0("../figures/","batch_controls_",conditions[i],"_boxplot.png"),
           plot=fig_controls, 
           height = 90, 
           width = 15,
           dpi = 300, 
           units="cm",
           device="png")

}

## plant batches
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


batch_plot <- ggplot(data=plant_batches_l)+
    geom_segment(
             mapping=aes(y=as.factor(batch_id),
                         x=date,
                         xend=date_end,
                         color=steps),
                         linewidth= 10
             )+
    geom_text(aes(x = midpoint,
                  y=as.factor(batch_id),
                  label = duration),
              color = "black") +
    geom_text(aes(x = endpoint,
                  y=as.factor(batch_id),
                  label = paste0(condition,"\nmicrobes=",microbes,sep="")),
              hjust=0,
              size=3,
              color = "black") +
    scale_color_manual(values=c("plant_growth"="forestgreen",
                                "inoculation"="goldenrod1",
                                "measurements"="lightskyblue")) +
    scale_x_date(
        date_breaks = "1 month",                # Breaks every month
        date_labels = "%m, %Y",
        limits=as.Date(c("2023-09-01","2025-02-01"))
        ) +
    labs(
        title = "Gantt Chart of in planta experiments",
        x = "Date",
        y = "Batches",
        color="Step",
        ) +
    theme_bw() +
    theme(
          legend.position = "inside",
          legend.position.inside=c(0.08,0.78),
          panel.grid.minor = element_blank()
    )

ggsave("../figures/batch_dates_barplot.png",
       plot=batch_plot, 
       height =15, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")


## Statistics



