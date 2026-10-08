# by: Floriane Coulmance: 16/08/2024
# usage:
# Rscript Fig3-4.R 
#___________________________________________________________________


# Clear the work space
rm(list = ls())


# new libraries
source("../scripts/R/helper_functions.R")
library(smartsnp)
library(ggplot2)
library(ggimage)
library(scales)
library(dplyr)
library(stringr)
library(ggtext)
library(ggpubr)
library(patchwork)
library(cowplot)
library(pairwiseAdonis)
library(tidyverse)
library(reshape2)
library(ggnewscale)
library(tibble)
library(png)
library(grid)
library(ggtree)
library(glue)
library(vegan)
library(knitr)
library(plotly)
library(base64enc)
library(scatterplot3d)
library(nlme)
library(lme4)
library(lmerTest)

# ############################
# CONFIG
# ############################


# ============================================================
# Parse command line arguments from Snakemake
# ============================================================
args <- commandArgs(trailingOnly = TRUE)


# ============================================================
# Paths passed from Snakemake
# ============================================================
base_path      <- "/Users/fcoulman/Desktop/hamlet_pheno/3_CHAPTER3/hamlet_pheno/"
figure_path     <- file.path(base_path, "figures")
path_phenotypes <- file.path(base_path, "1_phenotyping/pca")
path_genotypes_all <- file.path(base_path, "2_popgen/byALL")
path_genotypes_loc <- file.path(base_path, "2_popgen/byLOC")
fst_all_file <- file.path(figure_path, "TableS5.tex")
fst_loc_file <- file.path(figure_path, "TableS4.tex")
association_file <- file.path(base_path, "metadata/assortative_mating.csv")
transect_file <- file.path(base_path, "metadata/assortative_counts.csv")
logos_path      <- file.path(base_path, "metadata/logos_hamlet")
spec_colors     <- file.path(base_path, "metadata/species_colors.tsv")
geo_colors      <- file.path(base_path, "metadata/locations_colors.tsv")


# 
# base_path      <- get_arg("--base_path", ".")
# figure_path     <- get_arg("--figure_path", file.path(base_path, "figures"))
# path_phenotypes <- get_arg("--path_phenotypes", file.path(base_path, "1_phenotyping/pca"))
# path_genotypes_all <- file.path(
#   base_path,
#   "2_popgen",
#   "byALL"
# )
# path_genotypes_loc <- file.path(
#   base_path,
#   "2_popgen",
#   "byLOC"
# )
# fst_all_file <- file.path(
#   figure_path,
#   "TableS5.tex"
# )
# fst_loc_file <- file.path(
#   figure_path,
#   "TableS4.tex"
# )
# association_file <- file.path(
#   base_path,
#   "metadata",
#   "assortative_mating.csv"
# )
# transect_file <- file.path(
#   base_path,
#   "metadata",
#   "assortative_counts.csv"
# )
# logos_path      <- get_arg("--logos_path", file.path(base_path, "metadata/logos_hamlet"))
# spec_colors     <- get_arg("--spec_colors", file.path(base_path, "metadata/species_colors.tsv"))
# geo_colors      <- get_arg("--geo_colors", file.path(base_path, "metadata/locations_colors.tsv"))
# 


if (!dir.exists(figure_path)) {
  dir.create(
    figure_path,
    recursive = TRUE
  )
}


# ============================================================
# Print summary for debugging
# ============================================================
cat("---- CONFIG ----\n")
cat("base_path:      ", base_path, "\n")
cat("figure_path:    ", figure_path, "\n")
cat("path_phenotypes:", path_phenotypes, "\n")
cat("path_genotypes_all:", path_genotypes_all, "\n")
cat("path_genotypes_loc:", path_genotypes_loc, "\n")
cat("association_file:", association_file, "\n")
cat("logos_path:     ", logos_path, "\n")
cat("spec_colors:    ", spec_colors, "\n")
cat("geo_colors:     ", geo_colors, "\n")
cat("-----------------\n")


# ============================================================
# GENERAL SETTINGS
# ============================================================
n_perm <- 10000
seed <- 123

species_info <- add_species_logos(spec_colors, logos_path)
# head(species_info)
geo_table <- read.delim(geo_colors, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# head(geo_table)


message("\n========================================")
message("PHENOTYPE ANALYSIS")
message("========================================")

pheno_files <- list.files(
  path = path_phenotypes,
  pattern = "_PCs\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(pheno_files) == 0) {
  stop(
    "No phenotype PCA files found in: ",
    path_phenotypes
  )
}

message(
  "Phenotype files found: ",
  length(pheno_files)
)

print(pheno_files)

pheno_info <- tibble(
  file = pheno_files,
  
  dataset = basename(file) %>%
    str_remove("_PCs\\.csv$")
) %>%
  mutate(
    location_code = str_extract(
      dataset,
      "(?<=lab_)[A-Za-z]+(?=[0-9_])"
    ),
    level = case_when(
      
      str_detect(
        dataset,
        "^lab_571_"
      ) ~ "all",
      
      TRUE ~ "location"
    )
  ) %>%
  filter(
    dataset %in% c(
      "lab_flo29_left_noflash",
      "lab_boc229_left_noflash",
      "lab_bel46_left_noflash",
      "lab_571_left_noflash"
    )
  )


print(pheno_info$level)

pheno_distances <- vector(
  "list",
  nrow(pheno_info)
)

pheno_distances_lda <- vector(
  "list",
  nrow(pheno_info)
)

for (i in seq_len(nrow(pheno_info))) {
  
  message(
    "\nPhenotype dataset ",
    i,
    "/",
    nrow(pheno_info),
    ": ",
    pheno_info$dataset[i]
  )
  
  pheno_distances[[i]] <-
    calculate_pheno_distance(
      pca_file = pheno_info$file[i],
      dataset = pheno_info$dataset[i],
      level = pheno_info$level[i],
      species_info = species_info,
      geo_table = geo_table
    )
  
  pheno_distances_lda[[i]] <-
    calculate_pheno_distance_lda(
      pca_file = pheno_info$file[i],
      dataset = pheno_info$dataset[i],
      level = pheno_info$level[i],
      species_info = species_info,
      geo_table = geo_table
    )
}

# pheno_distances <- bind_rows(
#   pheno_distances
# )

pheno_distances_lda <- bind_rows(
  pheno_distances_lda
)


message("\n========================================")
message("GENOTYPE ANALYSIS")
message("========================================")
# 
# geno_files_all <- list.files(
#   path = path_genotypes_all,
#   pattern = ".agg.ld_pruned_gtmat\\.traw$",
#   recursive = TRUE,
#   full.names = TRUE
# )
# 
# geno_files_loc <- list.files(
#   path = path_genotypes_loc,
#   pattern = "_ld_pruned_gtmat\\.traw$",
#   recursive = TRUE,
#   full.names = TRUE
# )
# 
# geno_files <- unique(
#   c(
#     geno_files_all,
#     geno_files_loc
#   )
# )
# 
# 
# if (length(geno_files) == 0) {
#   stop(
#     "No genotype files found."
#   )
# }
# 
# geno_info <- tibble(
#   file = geno_files,
#   dataset = basename(file) %>%
#     str_remove("_ld_pruned_gtmat\\.traw$")
# ) %>%
#   mutate(
#     level = case_when(
#       str_detect(
#         tolower(file),
#         "all"
#       ) ~ "all",
# 
#       str_detect(
#         tolower(file),
#         "byloc"
#       ) ~ "location",
# 
#       TRUE ~ NA_character_
#     )
#   )
# 
# 
# print(geno_info)
# 
# geno_distances <- vector(
#   "list",
#   nrow(geno_info)
# )
# 
# 
# for (i in seq_len(nrow(geno_info))) {
# 
#   message(
#     "\nGenotype dataset ",
#     i,
#     "/",
#     nrow(geno_info),
#     ": ",
#     geno_info$dataset[i]
#   )
# 
#   geno_distances[[i]] <-
#     calculate_geno_distance(
#       gtmat_file = geno_info$file[i],
#       dataset = geno_info$dataset[i],
#       level = geno_info$level[i],
#       species_info = species_info,
#       geo_table = geo_table
#     )
# }
# 
# geno_distances <- bind_rows(
#   geno_distances
# )


geno_distances_fst <- read_fst_tables(
  s4_file = fst_loc_file,
  s5_file = fst_all_file,
  species_info = species_info,
  geo_table = geo_table
)


message("\n========================================")
message("ASSOCIATION / RI ANALYSIS")
message("========================================")

transect <- read.csv(
  transect_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Location = as.character(Location),
    species = as.character(species),
    count = as.numeric(count)
  )

association <- read.csv(
  association_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Location = as.character(Location),
    species1 = as.character(species1),
    species2 = as.character(species2),
    spawning = as.numeric(spawning)
  )

locations <- sort(
  unique(
    association$Location
  )
)


association_global <- association %>%
  group_by(
    species1,
    species2
  ) %>%
  summarise(
    count = sum(
      spawning,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

global_species_counts <- transect %>%
  group_by(species) %>%
  summarise(
    count = sum(count),
    .groups = "drop"
  ) %>%
  mutate(
    proportion = count / sum(count)
  )

RI_global <- calculate_global_RI(
  association_global,
  global_species_counts,
  species1_col = "species1",
  species2_col = "species2",
  spawning_count_col = "count",
  species_col = "species",
  species_count_col = "count")


asso_RI_global <- RI_global$results %>%
  select(
    species1,
    species2,
    RI
  ) %>%
  rename(
    distance_asso = RI
  ) %>%
  mutate(
    level = "all",
    Location = "all"
  )


asso_RI_location <- vector(
  "list",
  length(locations)
)

for (i in seq_along(locations)) {
  
  loc <- locations[i]
  
  message(
    "\nAssociation location ",
    i,
    "/",
    length(locations),
    ": ",
    loc
  )
  
  
  # --------------------------------------------------------
  # Extract this location
  # --------------------------------------------------------
  
  pairing_table <- association %>%
    filter(
      Location == loc
    ) %>%
    select(
      species1,
      species2,
      spawning
    ) %>%
    rename(
      count = spawning
    )
  
  print(pairing_table)
  
  species_counts <- transect %>%
    filter(
      Location == loc
    ) %>%
    mutate(
      proportion = count / sum(count)
    )

  print(species_counts)
  
  
  # --------------------------------------------------------
  # Run RI permutation
  # --------------------------------------------------------
  
  RI <- calculate_global_RI(
    pairing_table,
    species_counts,
    species1_col = "species1",
    species2_col = "species2",
    spawning_count_col = "count",
    species_col = "species",
    species_count_col = "count")

  print(RI$results)
  
  # --------------------------------------------------------
  # Standardise output
  # --------------------------------------------------------
  
  asso_RI_location[[i]] <-
    RI$results %>%
    select(
      species1,
      species2,
      RI
    ) %>%
    rename(
      distance_asso = RI
    ) %>%
    mutate(
      level = "location",
      Location = loc
    )
}


asso_RI <- bind_rows(
  asso_RI_global,
  bind_rows(
    asso_RI_location
  ) 
) %>%mutate(
  Location = recode(
    Location,
    "Bocas del Toro" = "Panama"
  )
)


message("\n========================================")
message("INTERSECTIONS")
message("========================================")
speciation_hypercube_data <- pheno_distances_lda %>%
  dplyr::inner_join(
    geno_distances_fst,
    by = c(
      "level",
      "Location",
      "species1",
      "species2"
    )
  ) %>%
  dplyr::inner_join(
    asso_RI,
    by = c(
      "level",
      "Location",
      "species1",
      "species2"
    )
  )
print(speciation_hypercube_data)

# correlation_figure <- plot_all_pairwise_correlations(
#   speciation_hypercube_data
# )
# correlation_figure


message("\n========================================")
message("HYPERCUBE")
message("========================================")
hypercube <- plot_speciation_hypercube(
  speciation_hypercube_data %>% filter(level == "location")
)
# print(hypercube)

data_location <- speciation_hypercube_data %>%
  filter(level == "location")

scatter_location <- scatterplot3d(
  data_location[, c(5, 7, 6)],
  pch = 16,
  color = "#D06495",
  angle = 135,
  aspect = c(1, 1, 1),
  xlim = c(0, 1),
  ylim = c(0, 1),
  zlim = c(0, 1),
  xlab = "Phenotypic divergence",
  ylab = "Reproductive isolation",
  zlab = "Genetic divergence (Fst)"
)
grid.echo()
scatter_location <- grid.grab()

p_pheno_geno_location <- plot_speciation_paper(
  data = speciation_hypercube_data,
  species_meta = species_info,
  x = "distance_pheno",
  y = "distance_geno",
  colour = "distance_asso",
  x_label = "Phenotypic divergence",
  y_label = "Genetic divergence (FST)",
  colour_label = "Reproductive isolation",
  panel = "location"
)

p_pheno_geno_location <- plot_pairwise_lmer(
  speciation_hypercube_data,
  data_level = "location",
  x = "distance_pheno",
  y = "distance_geno",
  x_lab = "Phenotypic divergence",
  y_lab = "Genetic differentiation (FST)"
)

p_pheno_asso_location <- plot_pairwise_lmer(
  speciation_hypercube_data,
  data_level = "location",
  x = "distance_pheno",
  y = "distance_asso",
  x_lab = "Phenotypic divergence",
  y_lab = "Reproductive isolation"
)

p_geno_asso_location <- plot_pairwise_lmer(
  speciation_hypercube_data,
  data_level = "location",
  x = "distance_geno",
  y = "distance_asso",
  x_lab = "Genetic differentiation (FST)",
  y_lab = "Reproductive isolation"
)




data_all <- speciation_hypercube_data %>%
  filter(level == "all")

p_pheno_geno_all <- plot_pairwise_lmer(
  speciation_hypercube_data,
  data_level = "all",
  x = "distance_pheno",
  y = "distance_geno",
  x_lab = "Phenotypic divergence",
  y_lab = "Genetic differentiation (FST)"
)


# ############################
# FINAL PLOTS
# ############################
########## FIGURE 3 #################### 
figure3 <- ggarrange(
  scatter_location,
  paper_location,
  nrow = 2,
  labels = c("(a)", "(b)"),
  heights = c(3.5, 5),
  widths = c(2.5, 4)
  )

ggsave(
  filename = file.path(figure_path, "Fig3_speciationLOC.png"),
  plot = figure3,
  width = 16,
  height = 19,
  units = "in",
  dpi = 300,
  type = "cairo-png"
)

htmlwidgets::saveWidget(
  hypercube,
  file.path(figure_path, "Fig3_speciation3D.html"),
  selfcontained = FALSE
)

########## FIGURE S21 ###################
figureS21 <- correlation_figure

ggsave(
  filename = file.path(figure_path, "FigS21_correlations.png"),
  plot = figureS21,
  width = 9,
  height = 12,
  units = "in",
  dpi = 300,
  type = "cairo-png"
)

########## FIGURE S22 #################### 
figureS22 <- figure_all

ggsave(
  filename = file.path(figure_path, "FigS22_speciationALL.png"),
  plot = figureS22,
  width = 12,
  height = 19,
  units = "in",
  dpi = 300,
  type = "cairo-png"
)
