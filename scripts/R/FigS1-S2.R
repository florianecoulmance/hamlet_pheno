# by: Floriane Coulmance: 07/11/2025
# usage:
# Rscript FigS1-S2.R 
#___________________________________________________________________

# Clear the work space
rm(list = ls())


# new libraries
source("helper_functions.R")
library(dplyr)
library(data.table)
library(tidyverse)
library(grid)
library(janitor)
library(ggplot2)
library(forcats)
library(readr)
library(ggtext)
library(purrr)
library(reshape2)
library(tidyr)
library(ggmap)
library(scatterpie)
library(ggpubr)
library(png)
library(ggimage)


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
logos_path      <- get_arg("--logos_path", file.path(base_path, "metadata/logos_hamlet"))
spec_colors     <- get_arg("--spec_colors", file.path(base_path, "metadata/species_colors.tsv"))
geo_colors      <- get_arg("--geo_colors", file.path(base_path, "metadata/locations_colors.tsv"))
pheno_meta      <- get_arg("--pheno_meta", file.path(base_path, "metadata/pheno_metadata.csv"))
geno_names      <- get_arg("--geno_names", file.path(base_path, "metadata/geno_names.txt"))

# ============================================================
# Print summary for debugging
# ============================================================
cat("---- CONFIG ----\n")
cat("base_path:      ", base_path, "\n")
cat("figure_path:    ", figure_path, "\n")
cat("logos_path:     ", logos_path, "\n")
cat("spec_colors:    ", spec_colors, "\n")
cat("geo_colors:     ", geo_colors, "\n")
cat("pheno_meta:     ", pheno_meta, "\n")
cat("geno_names:     ", geno_names, "\n")
cat("-----------------\n")


# ############################
# ANALYSIS
# ############################

species_info <- add_species_logos(spec_colors, logos_path)
head(species_info)
geo_table <- read.delim(geo_colors, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
head(geo_table)


# Read metadata and extract geo/spec
pheno_dat <- read.csv(pheno_meta, header = TRUE) %>%
  mutate(
    geo = case_when(
      grepl("PL17_0", sample_label) ~ substr(sample_label, 12, 14),
      grepl("PL17_1", sample_label) ~ substr(sample_label, 12, 14),
      grepl("PL22", sample_label)   ~ substr(sample_label, 13, 15),
      grepl("PL17_65", sample_label) ~ substr(sample_label, 11, 13),
      grepl("PL17_84", sample_label) ~ substr(sample_label, 11, 13)
    ),
    spec = case_when(
      grepl("PL17_0", sample_label) ~ substr(sample_label, 9, 11),
      grepl("PL17_160", sample_label) ~ substr(sample_label, 12, 14),
      grepl("PL17_1", sample_label) ~ substr(sample_label, 9, 11),
      grepl("PL22", sample_label)   ~ substr(sample_label, 10, 12),
      grepl("PL17_65", sample_label) ~ substr(sample_label, 8, 10),
      grepl("PL17_84", sample_label) ~ substr(sample_label, 8, 10)
    )
  )

geno_dat <- read_lines(geno_names) %>%
  as_tibble() %>% 
  rename(Sample = value) %>%
  mutate(
    SampleID = str_sub(Sample, end = -7),
    spec = str_sub(Sample, -6, -4),
    geo = str_sub(Sample, -3, -1)
  ) %>%
  filter(!(spec %in% c("tor", "tab", "tig")))


# List of datasets and options
datasets <- list(
  pheno = list(data = pheno_dat, label_col = "Locations", radius_factor = 0.9),
  geno  = list(data = geno_dat, label_col = "Locations", radius_factor = 0.9)
)

# Initialize results list
results <- list()

# Loop through datasets and generate plots
for(name in names(datasets)) {

    #-----------------------------------
    # Get datasets and parameters
    #-----------------------------------
    dat <- datasets[[name]]$data
    label_col <- datasets[[name]]$label_col
    radius_factor <- datasets[[name]]$radius_factor

    #-----------------------------------
    # Plot tally
    #-----------------------------------
    tally <- plot_species_geo_overview(dat, species_info, geo_table)
    
    #-----------------------------------
    # Shape data for mapping
    #-----------------------------------
    dcast_dat <- make_summary_by_location(dat, geo_table, radius_factor)
    print(dcast_dat)
    #-----------------------------------
    # Plot sample map
    #-----------------------------------
    map <- plot_location_pies(
    data_table = dcast_dat,
    species_info = species_info,
    radius_factor = radius_factor
    )

    #-----------------------------------
    # Store outputs
    #-----------------------------------
    results[[name]] <- list(
    p_tally = tally,
    p_map = map
  )
}

# ############################
# FINAL PLOTS
# ############################

########## FIGURE S1 ###################
pheno_t <- results[["pheno"]][["p_tally"]]
pheno_m <- results[["pheno"]][["p_map"]]

figureS1 <- ggarrange(pheno_t, pheno_m, ncol = 1, nrow = 2, labels=c("(a)","(b)"))

ggsave(filename = file.path(figure_path, "FigS1_pSampling.png"),
  plot = figureS1,
  width = 18.5,    # A4 width in inches
  height = 24,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S2 ###################
geno_t <- results[["geno"]][["p_tally"]]
geno_m <- results[["geno"]][["p_map"]]

figureS2 <- ggarrange(geno_t, geno_m, ncol = 1, nrow = 2, labels=c("(a)","(b)"))
ggsave(filename = file.path(figure_path, "FigS2_gSampling.png"),
  plot = figureS2,
  width = 18.5,    # A4 width in inches
  height = 24,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)