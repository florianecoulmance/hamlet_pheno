
# by: Floriane Coulmance: 04/05/2026
# usage:
# Rscript FigS21.R
#___________________________________________________________________

# Clear the work space
rm(list = ls())


# new libraries
source("../scripts/R/helper_functions.R")
library(data.table)
library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

base_dir <- args[1]
fig_dir  <- args[2]

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


ld_datasets <- c("all", "bel", "boc", "hon", "pri", "pue", "uni", "nig", "abe")


# -----------------------------
# loop over datasets
# -----------------------------
plots <- lapply(ld_datasets, build_plot)

# ensure 3x3 layout (pad if needed)
length(plots) <- 9

# -----------------------------
# plot grid (600 dpi PNG)
# -----------------------------
png(
  filename = file.path(fig_dir, "FigS21_LDboxplots.png"),
  width = 12,
  height = 12,
  units = "in",
  res = 600
)

grid.arrange(
  grobs = plots,
  ncol = 3,
  nrow = 3
)

dev.off()