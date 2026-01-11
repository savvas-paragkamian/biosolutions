#!/usr/bin/env Rscript

###############################################################################
# script name: phygenetics_figures.R
# developed by: Savvas Paragkamian
# framework: SarrisLab
###############################################################################
# GOAL:
# Aim of this script is to analyse and  visualise the de novo trees from the
# gtdb-tk pipeline.
###############################################################################
# usage:./biosolutions_phygenetics.R
###############################################################################

library(tidyverse)
library(ape)
library(ggtree)
library(tidytree)
library(treeio)
library(ggtreeExtra)
library(ggnewscale)

# load data from gtdb tk results
taxonomy <- read_tsv("../results/gtdbtk_denovo_infer/gtdbtk.bac120.decorated.tree-taxonomy", col_names = c("genome_id", "taxonomy"))
taxonomy_table <- read_tsv("../results/gtdbtk_denovo_infer/gtdbtk.bac120.decorated.tree-table")

## read the de novo tree
tree <- read.tree("../results/gtdbtk_denovo_infer/gtdbtk.bac120.decorated.tree")

## load the batchfile
assembly_taxonomy <- read_delim("../data/sarrislab_biosolutions_assembly_taxonomy.tsv", delim="\t")

# check which microbes are not included. SRL094, SRL243, SRL286, SRL324 are not
# in tree because the have high values of contamination
# 92.86 88.70 88.71 80.73

microbes_not_tree <- assembly_taxonomy[which(!(assembly_taxonomy$sample_alias %in% taxonomy$genome_id)),]

################### Tree with biosolutions assemblies only######################
biosolutions_genomes <- taxonomy %>%
  filter(str_detect(genome_id, "SRL")) |>
  pull(genome_id)


matching_tips <- intersect(tree$tip.label, biosolutions_genomes)
biosolutions_tree <- keep.tip(tree, matching_tips)

biosolutions_taxonomy <- taxonomy |>
    filter(genome_id %in% biosolutions_genomes) |>
    mutate(genus=str_remove(str_extract(taxonomy, "g__[^;]+"),"g__")) |>
    mutate(phylum=str_remove(str_extract(taxonomy, "p__[^;]+"),"p__")) |>
    mutate(class=str_remove(str_extract(taxonomy, "c__[^;]+"),"c__")) |>
    mutate(family=str_remove(str_extract(taxonomy, "f__[^;]+"),"f__")) |>
    mutate(species=str_extract(taxonomy, "s__[^;]+")) |>
    mutate(species= if_else(is.na(species),genome_id,species))

tree_df <- data.frame(label = biosolutions_tree$tip.label)

# Join with taxonomy
tip_data <- left_join(tree_df, biosolutions_taxonomy, by = c("label" = "genome_id")) |>
    left_join(assembly_taxonomy,by=c("species"="sample_alias")) |>
    mutate(id = factor(species, levels = biosolutions_tree$tip.label)) |>
    mutate(genome_size_mb=round(Genome_Size / 1e6,2))

biosolutions_tree_p <- ggtree(biosolutions_tree) %<+% tip_data + 
  geom_tippoint(aes(color = genus), size = 2) +
  geom_tiplab(aes(label = species), size = 1.5) +
  theme(legend.position = "inside",legend.position.inside = c(0.2,0.7))

ggsave(plot=biosolutions_tree_p,
       "../figures/biosolutions_tree_p.pdf",
       device="pdf",
       height = 50,
       width=30,
       units="cm", limitsize = F)


tree_circular <- ggtree(biosolutions_tree,layout="ellipse") %<+% tip_data + 
  geom_tippoint(aes(color = genus),size=6) +
  geom_tiplab(
    aes(label = species),
    size = 2.5,
    align = TRUE,     # extend labels radially outward
    linesize = 0.2,
    offset = 0.02      # push labels slightly outward from tips
  ) +
  #kcoord_cartesian(clip = "off") +  # allows labels outside plot boundary
  theme(legend.position = "bottom",
        plot.margin = margin(3, 3, 3, 3, "cm")
  ) 

ggsave(plot=tree_circular,
       "../figures/biosolutions_tree_circular.png",
       device="png",
       height = 30,
       width = 30,
       units="cm", limitsize = F)

# with data

p <- ggtree(biosolutions_tree, layout = "circular", size = 0.35) %<+% tip_data +
  geom_tiplab2(size = 1.8, offset = 0.01, align = TRUE) +
  theme_tree2() +
  geom_fruit(
    geom = geom_col,
    mapping = aes(y = id, x = genome_size_mb),
    orientation = "y",
    width = 0.12,
    offset = 0.02,
    pwidth = 0.25,
    fill = "grey80"
  ) +
  xlab(NULL)


p <- ggtree(biosolutions_tree, layout = "circular", size = 0.35) %<+% tip_data +
  # 1) labels first, with bigger offset so they sit outside the tree
  geom_tiplab2(
    aes(label = label),
    size = 2.2,
    offset = 0.06,        # increase until you see them
    align = TRUE,
    linesize = 0.2
  ) +
  theme_tree2() +
  # 2) genome size ring as bars, colored by genome size (legend from fill scale)
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = genome_size_mb, fill = genome_size_mb),
    width = 0.12,
    offset = 0.02,
    pwidth = 0.28,
    color = NA
  ) +
  scale_fill_gradient(
    name = "Genome size (Mb)",
    low = "white",
    high = "darkgreen"
  ) +
  # 3) give extra room so labels aren’t clipped
  theme(
    plot.margin = margin(5.5, 60, 5.5, 5.5),  # big right margin for labels
    legend.position = "right"
  ) +
  xlab(NULL)


# colors for families
family_cols <- c(
  # Actinobacteria (cool side)
  "Streptomycetaceae"  = "#2C7BB6",
  "Dermabacteraceae"   = "#4A90C4",
  "Microbacteriaceae"  = "#6BAED6",
  "Micrococcaceae"     = "#9ECAE1",
  # Alphaproteobacteria
  "Beijerinckiaceae"   = "#66C2A5",
  "Rhizobiaceae"       = "#41AE76",
  # Betaproteobacteria
  "Burkholderiaceae"   = "#1B9E77",
  # Gammaproteobacteria
  "Xanthomonadaceae"   = "#A6D854",
  "Enterobacteriaceae" = "#E6AB02",
  "Pseudomonadaceae"   = "#ccb137",
  # Bacillales (warm side, coherent)
  "Brevibacillaceae"   = "#FDAE61",
  "Planococcaceae"     = "#F46D43",
  "Amphibacillaceae"   = "#E6550D",
  "Paenibacillaceae"   = "#D73027",
  # Bacillaceae complex (same hue, different intensity)
  "Bacillaceae"        = "#B2182B",
  "Bacillaceae_B"      = "#CB181D",
  "Bacillaceae_G"      = "#EF3B2C",
  "Bacillaceae_H"      = "#FB6A4A",
  # GTDB placeholder families (neutral warm greys)
  "DSM-1321"           = "#BDBDBD",
  "DSM-18226"          = "#969696"
)

# the master plot with multiple data

p <- ggtree(biosolutions_tree, layout = "rectangular", size = 0.5) %<+% tip_data +
  # 1) show tip labels clearly
  geom_tippoint(aes(color = family), size = 2) +
  geom_tiplab(aes(color = family), size = 3.5, align = TRUE, linesize = 0.5) +
  geom_tiplab(
    aes(label = genus, color = family),
    size   = 3.5,
    offset = 0.15,
    align = TRUE,
    linesize = 0,
    fontface = "italic" # pushes genus to the right of the main label
  ) +
  geom_tiplab(
    aes(label = family, color = family),
    size   = 3.5,
    offset = 0.51,
    align = TRUE,
    linesize = 0,
    fontface = "italic" # pushes genus to the right of the main label
  ) +
  scale_color_manual(
  name   = "Family",
  values = family_cols,
  drop   = FALSE,
  na.value = "grey80",
  guide = "none"
  ) +
  theme_tree2() +
  # 2) genome size as a color-intensity tile band
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = genome_size_mb),
    width  = 0.1,     # vertical height per tip
    offset = 0.5,    # distance from tree/labels
    pwidth = 0.05     # horizontal space for the band
  ) +
  scale_fill_gradient(
    name = "Genome size (Mb)",
    low  = "lightblue1",
    high = "dodgerblue4",
    limits=c(min(tip_data$genome_size_mb),max(tip_data$genome_size_mb))
  ) +
ggnewscale::new_scale_fill() +
  # Track 2: GC content (%)
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = GC_Content),
    width = 0.1, offset = 0.02, pwidth = 0.05
  ) +
  scale_fill_gradient(
    name = "GC (%)",
    low = "ivory2", high = "lightgoldenrod4",
    limits = range(tip_data$GC_Content, na.rm = TRUE)
  ) +
  ggnewscale::new_scale_fill() +
  # Track 3: Total coding sequences
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = Total_Coding_Sequences),
    width = 0.1, offset = 0.02, pwidth = 0.05
  ) +
  scale_fill_gradient(
    name = "CDS (count)",
    low = "darkolivegreen2", high = "darkgreen",
    limits = range(tip_data$Total_Coding_Sequences, na.rm = TRUE)
  ) +
  ggnewscale::new_scale_fill() +
  # Track 4: Completeness (%)
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = Completeness),
    width = 0.1, offset = 0.02, pwidth = 0.05
  ) +
  scale_fill_gradient(
    name = "Completeness (%)",
    low = "white", high = "purple4",
    limits = range(tip_data$Completeness, na.rm = TRUE)
  ) +
  ggnewscale::new_scale_fill() +
  # Track 5: Contamination (%)
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = Contamination),
    width = 0.1, offset = 0.02, pwidth = 0.05
  ) +
  scale_fill_gradient(
    name = "Contamination (%)",
    low = "white", high = "red3",
    limits = range(tip_data$Contamination, na.rm = TRUE)
  ) +
  xlab(NULL) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave(plot=p,
       "../figures/biosolutions_tree_data.png",
       device="png",
       dpi = 600,
       height = 30,
       width = 40,
       units="cm", limitsize = F)
