# by: Floriane Coulmance: 16/08/2024
# usage:
# Rscript FigS22.R 
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
logos_path      <- get_arg("--logos_path", file.path(base_path, "metadata/logos_hamlet"))
spec_colors     <- get_arg("--spec_colors", file.path(base_path, "metadata/species_colors.tsv"))
geo_colors      <- get_arg("--geo_colors", file.path(base_path, "metadata/locations_colors.tsv"))

# ============================================================
# Print summary for debugging
# ============================================================
cat("---- CONFIG ----\n")
cat("base_path:      ", base_path, "\n")
cat("figure_path:    ", figure_path, "\n")
cat("path_phenotypes:", path_phenotypes, "\n")
cat("logos_path:     ", logos_path, "\n")
cat("spec_colors:    ", spec_colors, "\n")
cat("geo_colors:     ", geo_colors, "\n")
cat("-----------------\n")


# ############################
# ANALYSIS
# ############################

species_info <- add_species_logos(spec_colors, logos_path)
# head(species_info)
geo_table <- read.delim(geo_colors, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# head(geo_table)

# Define your locations and PCs of interest
dataset <- list(
    bel = list(name = "lab_bel46_left_noflash", dir = "byLOC", colors = "species"),
    boc = list(name = "lab_boc229_left_noflash", dir = "byLOC", colors = "species"),
    flk = list(name = "lab_flo29_left_noflash", dir = "byLOC", colors = "species"),
    all = list(name = "lab_571_left_noflash", dir = "byALL", colors = "species")
)


results <- list()

for(dat in names(dataset)) {
    dat_info <- dataset[[dat]]
    dat_name <- dat_info$name
    dat_dir <- dat_info$dir
    color <- dat_info$colors

    message("Processing: ", dat)

    #-----------------------------------
    # Read phenotypic PCA and generate average distance matrix
    #-----------------------------------
    pca_file <- file.path(path_phenotypes, paste0(dat_name, "_PCs.csv"))
    pca <- read.csv(pca_file, sep=",")
    
    pc_table <- write_metadata_gxp(pca)
    print(pc_table)

    # Average by species
    pheno_avg <- average_pca_by_species(
        pc_table,
        species_col = "spec",
        pc_prefix = "PC",
        min_n = 5
    )

    print(pheno_avg)

    #-----------------------------------
    # Read genetic PCA and generate average distance matrix
    #-----------------------------------
    gtmat_file <- file.path(base_path, "2_popgen", dat_dir, if (dat %in% "all") "all.agg.ld_pruned_gtmat.traw" else paste0(dat, "_ld_pruned_gtmat.traw"))
    sample_file <- file.path(base_path, if (dat %in% "all") "metadata/geno_names.txt" else paste0("2_popgen/", dat_dir, "/", dat, ".txt"))
 
    pca_res <- pca_analysis(gtmat_file, sample_file, color_by = color)
    pca_eigen <- pca_res$eigen
    print(pca_eigen)

    # Average by species
    geno_avg <- average_pca_by_species(
        pca_eigen,
        species_col = "spec",
        pc_prefix = "PC",
        min_n = 5
    )

    print(geno_avg)

    #-----------------------------------
    # Filter for species pairs present in both pheno and geno datasets
    #-----------------------------------
    shared_species <- intersect(
        pheno_avg$spec,
        geno_avg$spec
    )
    print(shared_species)

    pheno_shared <- pheno_avg %>%
        filter(spec %in% shared_species) %>%
        arrange(spec)
    print(pheno_shared)

    geno_shared <- geno_avg %>%
        filter(spec %in% shared_species) %>%
        arrange(spec)
    print(geno_shared)


    #-----------------------------------
    # Build pairwise distance matrix
    #-----------------------------------
    # Remove species column
    pheno_mat <- pheno_shared %>%
        column_to_rownames("spec")
    print(pheno_mat)

    geno_mat <- geno_shared %>%
        column_to_rownames("spec")
    print(geno_mat)

    # Pairwise Euclidean distances
    pheno_dist <- dist(pheno_mat, method = "euclidean")
    print(pheno_dist)
    geno_dist  <- dist(geno_mat, method = "euclidean")
    print(geno_dist)


    #-----------------------------------
    # Create joint tables for plots
    #-----------------------------------
    # Transform to longer table
    pheno_dist_table <- as.data.frame(as.table(as.matrix(pheno_dist))) %>%
        setNames(c("species1", "species2", "distance_pheno")) %>%
        mutate(
            species1 = as.character(species1),
            species2 = as.character(species2)
        ) %>%
        filter(species1 < species2)
    print(pheno_dist_table)

    geno_dist_table <- as.data.frame(as.table(as.matrix(geno_dist))) %>%
        setNames(c("species1", "species2", "distance_geno")) %>%
        mutate(
            species1 = as.character(species1),
            species2 = as.character(species2)
        ) %>%
        filter(species1 < species2)
    print(geno_dist_table)

    # Join tables
    phenoGeno_t <- pheno_dist_table %>%
        inner_join(geno_dist_table,
                   by = c("species1", "species2"))

    print(phenoGeno_t)

    
    #-----------------------------------
    # Mantel tests and plots
    #-----------------------------------
    pearsonMantel <- mantel(
        pheno_dist,
        geno_dist,
        method = "pearson",
        permutations = 9999
    )

    spearmanMantel <- mantel(
        pheno_dist,
        geno_dist,
        method = "spearman",
        permutations = 9999
    )


    phenoGeno_p <- plot_distance_correlation(
        df = phenoGeno_t,
        x_col = "distance_pheno",
        y_col = "distance_geno",
        x_lab = "Phenotypic distance",
        y_lab = "Genetic distance",
        title = paste0(dat, " phenotype vs genotype distances"),
        pearson_res = pearsonMantel,
        spearman_res = spearmanMantel
    )

    #-----------------------------------
    # Store outputs
    #-----------------------------------
    results[[dat]] <- list(
        pheno_D = pheno_dist,
        geno_D = geno_dist,
        phenoGeno_table = phenoGeno_t,
        phenoGeno_Mantel =  phenoGeno_p
    )
}

# ############################
# FINAL PLOTS
# ############################



#-----------------------------------
# Plot all sympatric species pairs together
#-----------------------------------
selected <- c("bel", "boc", "flk")

phenoD_combined <- bind_rows(
    lapply(selected, function(x) {

        as.data.frame(
            as.table(as.matrix(results[[x]]$pheno_D))
        ) %>%
            setNames(c("species1", "species2", "distance_pheno")) %>%
            mutate(
                species1 = as.character(species1),
                species2 = as.character(species2),
                dataset = x
            ) %>%
            filter(species1 < species2)

    })
)
print(phenoD_combined)


genoD_combined <- bind_rows(
    lapply(selected, function(x) {

        as.data.frame(
            as.table(as.matrix(results[[x]]$geno_D))
        ) %>%
            setNames(c("species1", "species2", "distance_geno")) %>%
            mutate(
                species1 = as.character(species1),
                species2 = as.character(species2),
                dataset = x
            ) %>%
            filter(species1 < species2)

    })
)
print(genoD_combined)


phenoGeno_combined <- bind_rows(
    lapply(selected, function(x) {
        results[[x]]$phenoGeno_table %>%
            mutate(dataset = x)
    })
)
print(phenoGeno_combined)

pearsonMantel <- mantel(
    phenoD_combined,
    genoD_combined,
    method = "pearson",
    permutations = 9999
)

spearmanMantel <- mantel(
    phenoD_combined,
    genoD_combined,
    method = "spearman",
    permutations = 9999
)


phenoGeno_p_combined <- plot_distance_correlation(
    df = phenoGeno_combined,
    x_col = "distance_pheno",
    y_col = "distance_geno",
    x_lab = "Phenotypic distance",
    y_lab = "Genetic distance",
    title = paste0(dat, " phenotype vs genotype distances"),
    pearson_res = pearsonMantel,
    spearman_res = spearmanMantel
)

#-----------------------------------
# Retrieve over all location plots
#-----------------------------------
all_p <- results[["all"]][["phenoGeno_Mantel"]]

#-----------------------------------
# Plot sympatric vs all locations correlations
#-----------------------------------
mantel_p <- plot_grid(
    plotlist = c(all_p, phenoGeno_p_combined),
    ncol = 1,
    labels = c("(a)", "(b)"),
    label_size = 14,
    align = "v"
)

ggsave(
  filename = file.path(figure_path, "FigS22_corr.png"),
  plot = mantel_p,
  width = 6,
  height = 12,
  units = "in",
  dpi = 300,
  type = "cairo-png"
)