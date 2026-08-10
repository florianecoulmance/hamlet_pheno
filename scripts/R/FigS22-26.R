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
library(cowplot)
library(png)
library(dplyr)
library(ggplot2)
library(ggnewscale)


args <- commandArgs(trailingOnly = TRUE)

base_dir <- args[1]
fig_dir  <- args[2]

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# LOAD METADATA
# ---------------------------------------------------------
logos_path      <- get_arg("--logos_path", file.path(base_dir, "metadata/logos_hamlet"))
spec_colors     <- get_arg("--spec_colors", file.path(base_dir, "metadata/species_colors.tsv"))
cat("logos_path:     ", logos_path, "\n")
cat("spec_colors:    ", spec_colors, "\n")

species_info <- add_species_logos(spec_colors, logos_path)
head(species_info)

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
  map(~ {
    if (nrow(.x) == 0) return(NULL)
    plot_location(.x, species_info)
  }) %>%
  compact()   # removes NULLs

# ---------------------------------------------------------
# FINAL FIGURE
# ---------------------------------------------------------

# combine_plots(plots, file.path(fig_dir, "FigS20_newhybrid.png"))

# purrr::iwalk(plots, function(p, loc) {
#   ggsave(
#     filename = file.path(fig_dir, paste0("newHybrids_", loc, ".png")),
#     plot = p,
#     height = 20,
#     width = 20,
#     dpi = 600
#   )
# })

purrr::iwalk(seq_along(plots), function(i, ...) {

  loc <- names(plots)[i]

  ggsave(
    filename = file.path(
      fig_dir,
      paste0("FigS", 21 + i, "_newHybrids_", loc, ".png")
    ),
    plot = plots[[i]],
    height = 20,
    width = 20,
    dpi = 600
  )
})