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
base_path      <- get_arg("--base_path", ".")
figure_path     <- get_arg("--figure_path", file.path(base_path, "figures"))
path_phenotypes <- get_arg("--path_phenotypes", file.path(base_path, "1_phenotyping/pca"))
path_genotypes_all <- file.path(
  base_path,
  "2_popgen",
  "byALL"
)
path_genotypes_loc <- file.path(
  base_path,
  "2_popgen",
  "byLOC"
)
fst_all_file <- file.path(
  figure_path,
  "TableS5.tex"
)
fst_loc_file <- file.path(
  figure_path,
  "TableS4.tex"
)
association_file <- file.path(
  base_path,
  "metadata",
  "assortative_mating.csv"
)
logos_path      <- get_arg("--logos_path", file.path(base_path, "metadata/logos_hamlet"))
spec_colors     <- get_arg("--spec_colors", file.path(base_path, "metadata/species_colors.tsv"))
geo_colors      <- get_arg("--geo_colors", file.path(base_path, "metadata/locations_colors.tsv"))

# base_path      <- "/Users/fcoulman/Desktop/hamlet_pheno/3_CHAPTER3/hamlet_pheno/"
# figure_path     <- file.path(base_path, "figures")
# path_phenotypes <- file.path(base_path, "1_phenotyping/pca")
# path_genotypes_all <- file.path(base_path, "2_popgen/byALL")
# path_genotypes_loc <- file.path(base_path, "2_popgen/byLOC")
# fst_all_file <- file.path(figure_path, "TableS5.tex")
# fst_loc_file <- file.path(figure_path, "TableS4.tex")
# association_file <- file.path(base_path, "metadata/assortative_mating.csv")
# logos_path      <- file.path(base_path, "metadata/logos_hamlet")
# spec_colors     <- file.path(base_path, "metadata/species_colors.tsv")
# geo_colors      <- file.path(base_path, "metadata/locations_colors.tsv")

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

pheno_distances <- bind_rows(
  pheno_distances
)

pheno_distances_lda <- bind_rows(
  pheno_distances_lda
)


message("\n========================================")
message("GENOTYPE ANALYSIS")
message("========================================")

geno_files_all <- list.files(
  path = path_genotypes_all,
  pattern = ".agg.ld_pruned_gtmat\\.traw$",
  recursive = TRUE,
  full.names = TRUE
)

geno_files_loc <- list.files(
  path = path_genotypes_loc,
  pattern = "_ld_pruned_gtmat\\.traw$",
  recursive = TRUE,
  full.names = TRUE
)

geno_files <- unique(
  c(
    geno_files_all,
    geno_files_loc
  )
)


if (length(geno_files) == 0) {
  stop(
    "No genotype files found."
  )
}

geno_info <- tibble(
  file = geno_files,
  dataset = basename(file) %>%
    str_remove("_ld_pruned_gtmat\\.traw$")
) %>%
  mutate(
    level = case_when(
      str_detect(
        tolower(file),
        "all"
      ) ~ "all",

      str_detect(
        tolower(file),
        "byloc"
      ) ~ "location",

      TRUE ~ NA_character_
    )
  )


print(geno_info)

geno_distances <- vector(
  "list",
  nrow(geno_info)
)


for (i in seq_len(nrow(geno_info))) {

  message(
    "\nGenotype dataset ",
    i,
    "/",
    nrow(geno_info),
    ": ",
    geno_info$dataset[i]
  )

  geno_distances[[i]] <-
    calculate_geno_distance(
      gtmat_file = geno_info$file[i],
      dataset = geno_info$dataset[i],
      level = geno_info$level[i],
      species_info = species_info,
      geo_table = geo_table
    )
}

geno_distances <- bind_rows(
  geno_distances
)


geno_distances_fst <- read_fst_tables(
  s4_file = fst_loc_file,
  s5_file = fst_all_file,
  species_info = species_info,
  geo_table = geo_table
)


message("\n========================================")
message("ASSOCIATION / RI ANALYSIS")
message("========================================")

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


RI_global <- global_RI_permutation(
  pairing_table = association_global,
  n_perm = n_perm,
  seed = seed
)


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
  
  
  # --------------------------------------------------------
  # Run RI permutation
  # --------------------------------------------------------
  
  RI <- global_RI_permutation(
    pairing_table = pairing_table,
    n_perm = n_perm,
    seed = seed
  )
  
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
# pheno_geno <- inner_join(
#   pheno_distances,
#   geno_distances,
#   by = c("level", "Location", "species1", "species2")
# )

pheno_geno <- inner_join(
  pheno_distances_lda, 
  geno_distances_fst,
  by = c("level", "Location", "species1", "species2")
)

# pheno_asso <- inner_join(
#   pheno_distances,
#   asso_RI,
#   by = c("level", "Location", "species1", "species2")
# )

pheno_asso <- inner_join(
  pheno_distances_lda,
  asso_RI,
  by = c("level", "Location", "species1", "species2")
)

# geno_asso <- inner_join(
#   geno_distances,
#   asso_RI,
#   by = c("level", "Location", "species1", "species2")
# )

geno_asso <- inner_join(
  geno_distances_fst, 
  asso_RI,
  by = c("level", "Location", "species1", "species2")
)

pheno_geno_cor <- plot_pairwise_correlations(
  df = pheno_geno,
  x_col = "distance_pheno",
  y_col = "distance_geno",
  x_lab = "Phenotypic distance",
  y_lab = "Genetic distance",
  title_prefix = "Phenotype vs genotype"
)

pheno_geno_cor$all$plot
pheno_geno_cor$location$plot

pheno_asso_cor <- plot_pairwise_correlations(
  df = pheno_asso,
  x_col = "distance_pheno",
  y_col = "distance_asso",
  x_lab = "Phenotypic distance",
  y_lab = "Reproductive isolation",
  title_prefix = "Phenotype vs reproductive isolation"
)

pheno_asso_cor$all$plot
pheno_asso_cor$location$plot


geno_asso_cor <- plot_pairwise_correlations(
  df = geno_asso,
  x_col = "distance_geno",
  y_col = "distance_asso",
  x_lab = "Genetic distance",
  y_lab = "Reproductive isolation",
  title_prefix = "Genotype vs reproductive isolation"
)

geno_asso_cor$all$plot
geno_asso_cor$location$plot

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

hypercube <- plot_speciation_hypercube(
  speciation_hypercube_data
)
hypercube

hypercube_paper <- plot_speciation_paper(
  speciation_hypercube_data,
  species_info
)
hypercube_paper



speciation_hypercube_data2 <- pheno_distances %>%
  dplyr::inner_join(
    geno_distances,
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

hypercube2 <- plot_speciation_hypercube(
  speciation_hypercube_data
)
hypercube2

hypercube_paper2 <- plot_speciation_paper2(
  speciation_hypercube_data,
  species_info
)
hypercube_paper2

# ############################
# FINAL PLOTS
# ############################

########## FIGURE 3 ###################
figure3 <- plot_grid(
  plotlist = list(pheno_geno_cor$all$plot,
               pheno_geno_cor$location$plot,
               pheno_asso_cor$all$plot,
               pheno_asso_cor$location$plot,
               geno_asso_cor$all$plot,
               geno_asso_cor$location$plot
              ),
  ncol = 2,
  labels = c("(a)", "", "(b)", "", "(c)", ""),
  label_size = 14,
  align = "v"
)

ggsave(
  filename = file.path(figure_path, "Fig3_correlations.png"),
  plot = figure3,
  width = 12,
  height = 19,
  units = "in",
  dpi = 300,
  type = "cairo-png"
)


########## FIGURE 4 ###################
# figure4 <- hypercube
# 
# ggsave(
#   filename = file.path(figure_path, "Fig4_interactiveCUBE.png"),
#   plot = figure4,
#   width = 12,
#   height = 19,
#   units = "in",
#   dpi = 300,
#   type = "cairo-png"
# )

figure4 <- hypercube_paper2

ggsave(
  filename = file.path(figure_path, "Fig4_2dCUBE_v2.png"),
  plot = figure4,
  width = 12,
  height = 19,
  units = "in",
  dpi = 300,
  type = "cairo-png"
)