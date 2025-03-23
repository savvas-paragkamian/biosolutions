#!/usr/bin/env Rscript
## Script name: biosolutions_analysis.R
##
## Purpose of script: 
## Visualisation of the results and statistics of biosolutions project experiments.
## There are 3 different plant conditions that are inoculated with 
## ~60 microbial isolates to screen their role association with 
## plant growth.
## The plant model is Arabidopsis thaliana.
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-12-12
##
library(tidyverse)

# load data

plant_batches_l <- read_delim("../results/plant_batches_results.tsv",delim="\t")

all_experiments <- read_delim("../results/all_experiments_data.tsv",delim="\t")

all_experiments_stats <- read_delim("../results/all_experiments_stats.tsv",delim="\t")

## statistics
stats_results_anova <- read_delim("../results/stats_results_anova.tsv",delim="\t")
stats_results_kruskal <- read_delim("../results/stats_results_kruskal.tsv",delim="\t")

control_pairwise_sig <- read_delim("../results/control_pairwise_sig_microbes.tsv",delim="\t")
all_experiments_tukey_sig <- read_delim("../results/all_experiments_tukey_sig.tsv",delim="\t")

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
    select(-value) |>  # Remove the original value column
    pivot_wider(
                names_from = variable,
                values_from = percent_change
                ) |>
    filter(!(microbe_id %in% c("Control","Control+", "E. coli")))


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

    batch_data_a <- all_experiments |>
        filter(batch_id==batches[b])

    condition <- unique(batch_data$condition)
    print(batches[b])
    print(condition)
    print("bar plot")
    # bar plot

    for (i in seq_along(vars_mean)){
        print(i)
        print(vars_mean[i])
        print(vars_a[i])
    
        fig_bar <- ggplot()+
                    geom_col(batch_data,
                             mapping=aes(x = microbe_id,
                              y = !!sym(vars_mean[i])),
                              color="mediumseagreen",
                             fill="white",
                             width = 0.6,
                             position = "identity")+
                    geom_point(data=batch_data_a,
                               mapping=aes(x = microbe_id,
                                           y=!!sym(vars_a[i]),
                                           color="mediumseagreen"),
                               fill="mediumseagreen",
                               shape = 21,
                               alpha = 0.8,
                               show.legend = F,
                               position = position_jitterdodge(0.3))+
                    #position_dodge(width = 0.82)) +
                    geom_errorbar(batch_data,
                                  mapping=aes(x= microbe_id,
                                              ymin = !!sym(vars_mean[i]) - !!sym(vars_err[i]),
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

## plant batches Gantt Chart

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
        limits=as.Date(c("2023-09-01","2025-05-01"))
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


############################# summary statistics ########################

plant_batch_stats_kruskal <- plant_batches_l |>
    distinct(batch_id,condition,microbes,microbes_id) |>
    filter(!is.na(microbes_id)) |>
    right_join(stats_results_kruskal) |>
    rename("variable"="kruskal_variable")

plant_batch_stats_sig_k <- plant_batch_stats_kruskal |>
    filter(kruskal_p.value<0.05) |>
    group_by(condition, variable) |>
    summarise(n_batches=n(),
              batches=str_c(batch_id, collapse = ","),
              .groups="keep") |>
    mutate(statistic="kruskal")

plant_batch_stats_anova <- plant_batches_l |>
    distinct(batch_id,condition,microbes,microbes_id) |>
    filter(!is.na(microbes_id)) |>
    right_join(stats_results_anova) |>
    rename("variable"="anova_variable")

plant_batch_stats_sig <- plant_batch_stats_anova |>
    filter(anova_p.value<0.05) |>
    group_by(condition, variable) |>
    summarise(n_batches=n(),
              batches= str_c(batch_id, collapse = ","),
              .groups="keep") |>
    mutate(statistic="anova") |>
    bind_rows(plant_batch_stats_sig_k) |>
    mutate(midpoint = n_batches / 2)

plant_batch_sig_fig <- ggplot() +
    geom_col(data=plant_batch_stats_sig,
             mapping=aes(y=variable,
                         x=n_batches,
                         fill=statistic),
             width = 0.6,
             position = position_dodge(width = 0.82)) +
    geom_text(data=plant_batch_stats_sig,
              mapping=aes(y = variable,
                          x = midpoint,
                          group = statistic,
                          label = batches),
              color = "black",
              position = position_dodge(width = 0.82)) +
    labs(
         title = "Plant batches with p.value less than 0.05",
         y = "Variables",
         x = "# Batches") +
    theme_bw() +
    scale_fill_manual(values=c("kruskal"="goldenrod3",
                               "anova"="cadetblue")) +
    theme(
          legend.position = "inside",
          legend.position.inside=c(0.96,0.88),
          axis.text.y = element_text(size = 13),
          axis.text.x = element_text(angle = 0,
                                     hjust = 1,
                                     size = 13)
    )+
    facet_wrap(~condition,
               nrow=1,
               ncol=3)

    ggsave("../figures/plant_batches_statistics_sig.png",
           plot=plant_batch_sig_fig, 
           height = 20, 
           width = 50,
           dpi = 300, 
           units="cm",
           device="png")

### significant pairwise microbes with values higher than control
plant_batch_c <- plant_batches_l |>
    distinct(batch_id,condition)

control_pairwise_sig_sum <- all_experiments_tukey_sig |>
    distinct(microbe_id,batch_id, variable, tukey_adj.p.value) |>
    left_join(plant_batch_c) |> 
    mutate(
           significance = case_when(
                                    tukey_adj.p.value < 0.001 ~ "***",
                                    tukey_adj.p.value < 0.01 ~ "**",
                                    tukey_adj.p.value < 0.05 ~ "*",
                                    tukey_adj.p.value < 0.1 ~ ".",
                                    TRUE ~ ""  # For non-significant p-values
           )
    ) |>
    group_by(microbe_id,variable,significance,condition,tukey_adj.p.value) |>
    summarise(n_batches=n(),
              batches=str_c(batch_id, collapse = ","),
              .groups="keep") 


microbe_heatmap <- ggplot()+
      geom_tile(data=control_pairwise_sig_sum,
                aes(y=microbe_id, x=variable,fill=tukey_adj.p.value),
                color="gray90",
                alpha=1,
                show.legend = T)+
      geom_text(data=control_pairwise_sig_sum,
                aes(y=microbe_id, x=variable, label=significance),
                size=4) +
      scale_fill_gradient(low="#CC79A7",
                          high="gray87",
                          breaks=waiver(),
                          n.breaks=4,
                          limits=c(min(control_pairwise_sig_sum$tukey_adj.p.value),0.1),
                          na.value="white",
                          guide = "colorbar")+
      guides(fill = guide_colorbar(ticks = FALSE,
                                   title="tukey_adj\np.value",
                                   label.vjust = 0.8,
                                   title.vjust = 0.8))+
      ylab("") +
      xlab("")+
      theme_bw()+
      theme(
            panel.border=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor=element_blank(),
            axis.text.x = element_text(face="bold",angle = 90, hjust = 0),
            axis.text.y = element_text(face="bold"),
            axis.text = element_text(size=13), 
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.title=element_text(size=9)) +
      facet_wrap(~condition,
               nrow=1,
               ncol=3)

ggsave("../figures/significant_microbes_control.png",
       plot = microbe_heatmap,
       width = 60,
       height = 30,
       units='cm', 
       device = "png",
       dpi = 300)

