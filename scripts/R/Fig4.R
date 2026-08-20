# by: Floriane Coulmance: 19/08/2026
# usage:
# Rscript Fig4.R 
#___________________________________________________________________

# Clear the work space
rm(list = ls())

# new libraries
# source("../scripts/R/helper_functions.R")
source("helper_functions.R")


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
# base_path      <- get_arg("--base_path", ".")
base_path      <- "/Users/fcoulman/Desktop/hamlet_pheno/3_CHAPTER3/hamlet_pheno/"
# figure_path     <- get_arg("--figure_path", file.path(base_path, "figures"))
figure_path     <- file.path(base_path, "figures")
# logos_path      <- get_arg("--logos_path", file.path(base_path, "metadata/logos_hamlet"))
logos_path      <- file.path(base_path, "metadata/logos_hamlet")
# spec_colors     <- get_arg("--spec_colors", file.path(base_path, "metadata/species_colors.tsv"))
spec_colors     <- file.path(base_path, "metadata/species_colors.tsv")
# geo_colors      <- get_arg("--geo_colors", file.path(base_path, "metadata/locations_colors.tsv"))
geo_colors      <- file.path(base_path, "metadata/locations_colors.tsv")


# ============================================================
# Print summary for debugging
# ============================================================
cat("---- CONFIG ----\n")
cat("base_path:      ", base_path, "\n")
cat("figure_path:    ", figure_path, "\n")
cat("logos_path:     ", logos_path, "\n")
cat("spec_colors:    ", spec_colors, "\n")
cat("geo_colors:     ", geo_colors, "\n")
cat("-----------------\n")


# ############################
# ANALYSIS
# ############################

# 1. ASSORTATIVE MATING PERMUTATIONS
pairing_table <- read.delim(file.path(base_path, "metadata", "pairing_counts.txt"),
                            sep=" ", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)

RI <- global_RI_permutation(pairing_table, n_perm = 10000, seed = 123)
print(RI$results)

# 2. PHENOTYPIC DIVERGENCE

# 3. DXY/FST


# ############################
# FINAL PLOTS
# ############################

########## FIGURE 4 ##########


