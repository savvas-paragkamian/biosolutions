#!/usr/bin/env Rscript
## Script name: biosolutions_genomics.R
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
## Date Created: 2025-08-14
##

library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(GSEABase)

# load metadata
sample_metadata <- read_delim("../data/sarrislab_biosolutions_samples_PRJEB103843_metadata.tsv",
                       na=c("","NA"),
                       delim="\t")

sample_locations <- sample_metadata |>
    distinct(latitude, longitude, hostScientificName)

write_delim(sample_locations, "../data/sampling_locations.tsv",delim="\t")

# load the results of stastistical tests

all_experiments_dunn_sig <- read_delim("../results/all_experiments_dunn_sig.tsv",delim="\t")

# load all BAKTA annotations merged into one

bakta_annotations <- read_delim("../results/bakta_merged.tsv", delim="\t", skip=5)

bakta_annotations_l <- bakta_annotations |>
    separate_rows(DbXrefs, sep = ", ")

bakta_go <- bakta_annotations_l |>
    filter(grepl("^GO:",DbXrefs))

bakta_ko <- bakta_annotations_l |>
    filter(grepl("^KEGG:",DbXrefs))

## summaries per isolate
bakta_go_srl <- bakta_go |>
    distinct(assembly,DbXrefs) |>
    group_by(assembly) |>
    summarise(n_gos=n())

bakta_ko_srl <- bakta_ko |>
    distinct(assembly,DbXrefs) |>
    group_by(assembly) |>
    summarise(n_kos=n())

########################### GO slim ######################
go_tbl <- bakta_go |>
    filter(Type == "cds") |>
    transmute(
              assembly,
              gene_id = `Locus Tag`,
              GO = DbXrefs
              ) |>
    distinct()

## load he go slim for gene ontology 
goslim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_generic.obo")

# assign ids
slim_ids <- ids(goslim)   # vector of "GO:...." slim IDs

# prepare the assemblies
assemblies <- unique(go_tbl$assembly)
onts <- c("BP", "MF", "CC")

## initiate the objects
res <- vector("list", length(assemblies) * length(onts))
k <- 1

for (i in seq_along(assemblies)) {

  asm <- assemblies[i]

  go_vec <- unique(go_tbl$GO[go_tbl$assembly == asm])
  go_vec <- go_vec[!is.na(go_vec)]
  go_vec <- go_vec[grepl("^GO:", go_vec)]

  if (length(go_vec) == 0) next

  # get the ids and ontology names for the assembly
  idSrc <- GOCollection(go_vec)

  for (j in seq_along(onts)) {

    ont <- onts[j]

    gs <- goSlim(idSrc, goslim, ontology = ont)

    df <- as.data.frame(gs) |>
      rownames_to_column("GO_slim") |>
      mutate(
        assembly = asm,
        ontology = ont,
        n_input_go = length(go_vec),
        .before = 1
      )

    res[[k]] <- df
    k <- k + 1
  }
}

go_slim_all <- bind_rows(res)

write_delim(go_slim_all,"../results/bakta_go_slim_all.tsv")

# bubble plot

# ---- choose what to plot ----
ont_use <- "BP"          # "BP", "MF", or "CC"
top_terms <- 50          # keep top N GO-slim terms for readability

df <- go_slim_all %>%
  filter(ontology == ont_use, Count > 0) %>%
  mutate(GO_term = paste0(Term, " (", GO_slim, ")"))

# ---- order GO terms by global abundance ----
term_order <- df %>%
  group_by(GO_term) %>%
  summarise(total = sum(Count), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = top_terms) %>%
  pull(GO_term)

df <- df %>% filter(GO_term %in% term_order)

# ---- order assemblies by clustering on profiles ----
mat <- df %>%
  dplyr::select(assembly, GO_term, Percent) %>%
  pivot_wider(names_from = GO_term, values_from = Percent, values_fill = 0) %>%
  arrange(assembly)

m <- as.matrix(mat[,-1])
rownames(m) <- mat$assembly

hc <- hclust(dist(m), method = "average")
asm_order <- hc$labels[hc$order]

df <- df %>%
  mutate(
    assembly = factor(assembly, levels = asm_order),
    GO_term  = factor(GO_term, levels = rev(term_order))
  )

# ---- bubble plot ----
go_slim_bp_p <- ggplot(df, aes(x = assembly, y = GO_term, size = Count, fill = Percent)) +
  geom_point(shape = 21, color = "black", alpha = 0.85) +
  scale_fill_gradient(
    name = "Percent",
    low = "plum2", high = "#654792",
    limits = range(df$Percent, na.rm = TRUE)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size= 13, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size= 13),
    panel.grid.major = element_line(color = "grey90")
  ) +
  labs(
    x = "Isolate (clustered)",
    y = paste0("GO-Slim terms (Top ", top_terms, ", ", ont_use, ")"),
    size = "GO count",
    fill = "Percent"
  )

ggsave("../figures/go_slim_bp_bubble.png",
       plot = go_slim_bp_p,
       width = 45,
       height = 30,
       units='cm', 
       device = "png",
       dpi = 300)

