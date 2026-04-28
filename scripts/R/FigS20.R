# by: Floriane Coulmance: 20/04/2026
# usage:
# Rscript FigS20.R
#___________________________________________________________________

# Clear the work space
rm(list = ls())


# new libraries
source("../scripts/R/helper_functions.R")
library(tidyverse)
library(stringr)
library(patchwork)
library(ggtext)
library(paletteer)
library(prismatic)
library(ggthemes)

args <- commandArgs(trailingOnly = TRUE)

base_dir <- args[1]
fig_dir  <- args[2]

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# LOAD METADATA
# ---------------------------------------------------------

species_meta <- read_tsv(file.path(base_dir, "metadata/species_colors.tsv"))
print(species_meta)

# ---------------------------------------------------------
# LOAD RESULTS
# ---------------------------------------------------------

paths <- get_nh_result_paths(base_dir)
print(paths)

data <- map_dfr(paths, read_pofz)
print(data)


# ---------------------------------------------------------
# SPLIT BY LOCATION
# ---------------------------------------------------------

plots <- data %>%
  split(.$loc) %>%
  map(~ plot_location(.x, species_meta))

# ---------------------------------------------------------
# FINAL FIGURE
# ---------------------------------------------------------

combine_plots(plots, file.path(fig_dir, "FigS20_newhybrid.png"))