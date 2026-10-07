# by: Floriane Coulmance: 16/08/2024
# usage:
# Rscript Fig1-2_S3-S9.R 
#___________________________________________________________________

# Clear the work space
rm(list = ls())


# new libraries
source("../scripts/R/helper_functions.R")
library(ggplot2)
library(ggimage)
library(scales)
# library(dplyr)
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
# library(MASS)
# library(tidyr)
library(purrr)




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
path_phenotypes <- get_arg("--path_phenotypes", file.path(base_path, "1_phenotyping/pca"))
heatmap_path    <- get_arg("--heatmap_path", file.path(base_path, "1_phenotyping/heatmaps"))
figure_path     <- get_arg("--figure_path", file.path(base_path, "figures"))
logos_path      <- get_arg("--logos_path", file.path(base_path, "metadata/logos_hamlet"))
spec_colors     <- get_arg("--spec_colors", file.path(base_path, "metadata/species_colors.tsv"))
geo_colors      <- get_arg("--geo_colors", file.path(base_path, "metadata/locations_colors.tsv"))

# ============================================================
# Print summary for debugging
# ============================================================
cat("---- CONFIG ----\n")
cat("base_path:      ", base_path, "\n")
cat("path_phenotypes:", path_phenotypes, "\n")
cat("heatmap_path:   ", heatmap_path, "\n")
cat("figure_path:    ", figure_path, "\n")
cat("logos_path:     ", logos_path, "\n")
cat("spec_colors:    ", spec_colors, "\n")
cat("geo_colors:     ", geo_colors, "\n")
cat("-----------------\n")


# ############################
# ANALYSIS
# ############################

species_info <- add_species_logos(spec_colors, logos_path)
head(species_info)
geo_table <- read.delim(geo_colors, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
head(geo_table)

# Define your locations and PCs of interest
dataset <- list(
  ver = list(name = "lab_ver60_left_noflash", pcs = c("PC1", "PC2", "PC3", "PC4"), colors = "species"),
  tob = list(name = "lab_tob69_left_noflash", pcs = c("PC1", "PC3", "PC2", "PC4"), colors = "species"),
  flo = list(name = "lab_flo29_left_noflash", pcs = c("PC1", "PC4", "PC2", "PC3"), colors = "species"),
  bel = list(name = "lab_bel46_left_noflash", pcs = c("PC1", "PC2", "PC3", "PC4"), colors = "species"),
  boc = list(name = "lab_boc229_left_noflash", pcs = c("PC1", "PC3", "PC2", "PC4"), colors = "species"),
  uvi = list(name = "lab_uvi136_left_noflash", pcs = c("PC1", "PC4", "PC2", "PC3"), colors = "species"),
  all = list(name = "lab_571_left_noflash", pcs = c("PC1","PC4", "PC2", "PC3"), colors = "species"),
  pue = list(name = "lab_pue187_left_noflash", pcs = c("PC1","PC2", "PC3", "PC4"), colors = "location"),
  nig = list(name = "lab_nig111_left_noflash", pcs = c("PC1","PC2", "PC3", "PC4"), colors = "location"),
  uni = list(name = "lab_uni74_left_noflash", pcs = c("PC1","PC2", "PC3", "PC4"), colors = "location"),
  chl = list(name = "lab_chl34_left_noflash", pcs = c("PC1","PC2", "PC3", "PC4"), colors = "location"),
  abe = list(name = "lab_abe30_left_noflash", pcs = c("PC1","PC2", "PC3", "PC4"), colors = "location"),
  ind = list(name = "lab_ind28_left_noflash", pcs = c("PC1","PC2", "PC3", "PC4"), colors = "location")
)

# Create a list to store plots per location
results <- list()

for(dat in names(dataset)) {
  
  dat_info <- dataset[[dat]]
  dat_name <- dat_info$name
  pcs      <- dat_info$pcs
  color <- dat_info$colors
  
  message("Processing: ", dat_name)
  
  #-----------------------------------
  # Read PCA + variance + metadata
  #-----------------------------------
  pca_file <- file.path(path_phenotypes, paste0(dat_name, "_PCs.csv"))
  var_file <- file.path(path_phenotypes, paste0(dat_name, "_var.csv"))
  print(pca_file)
  print(var_file)

  pca <- read.csv(pca_file, sep=",")
  var <- read.csv(var_file, sep=",")
  
  pc_table <- write_metadata_gxp(pca)
  
  #-----------------------------------
  # PCA plot
  #-----------------------------------
  p_pca <- pca_plot(pc_table, pcs[1], pcs[2], species_info, geo_table, var, color_by = color, extract_legend = FALSE)
  s_pca <- pca_plot(pc_table, pcs[3], pcs[4], species_info, geo_table, var, color_by = color, extract_legend = FALSE)

  #-----------------------------------
  # VAR plot
  #-----------------------------------
  p_var <- plot_variance(var, dat_name)
  
  #-----------------------------------
  # PERMANOVA (filter <5 inds per species inside perm_f)
  #-----------------------------------
  p_perm <- perm_f(pc_table, species_info, geo_table, color_by = color)
  
  #-----------------------------------
  # Hierarchical clustering
  #-----------------------------------
  p_hclust <- hierClustering(path_phenotypes, paste0(dat_name, "_PCs.csv"), species_info, geo_table, color_by = color)
  
  #-----------------------------------
  # HEATMAPS fish body for PC combination
  #-----------------------------------
  p_heat <- heat_plots(heatmap_path, dat_name, c(pcs[1], pcs[2]), species_info, geo_table, color_by = color)


  #-----------------------------------
  # Discriminant analysis
  #-----------------------------------
  da <- discriminant_analysis(
    pc_table = pc_table,
    pcs = paste0("PC", 1:15),
    group_col = if (color == "species") "spec" else "geo"
  )
  print(da)

  # LDA plot
  p_lda <- lda_plot(
    da = da,
    species_info = species_info,
    geo_info = geo_table,
    color_by = color,
    dataset = dat_name,
    extract_legend = FALSE
  )


  #-----------------------------------
  # Store outputs
  #-----------------------------------
  results[[dat]] <- list(
    data = pc_table,
    variance = var,
    pca = p_pca,
    sup_pca = s_pca,
    variance_plot = p_var,
    permanova = p_perm,
    hclust = p_hclust,
    heatmap = p_heat,
    discriminant = da,
    lda = p_lda
  )
}




# ############################
# FINAL PLOTS
# ############################
# Create common legend
leg <- legend_plot(species_info)
leg_g <- legend_geo(geo_table)


########## FIGURE 1 ###################
# PCA plots for all, mexico, belize and usvi
resFig1 <- results[names(results) %in% c("all", "ver", "bel", "uvi")]
keep_names <- names(resFig1)
print(keep_names)

# Extract PCA and PERMANOVA plots
# Fig1_pcas <- lapply(resFig1, `[[`, "pca")
# Fig1_permanova <- lapply(resFig1, `[[`, "permanova")

right_pca <- plot_grid(
  results[["ver"]][["pca"]],
  results[["bel"]][["pca"]],
  results[["uvi"]][["pca"]],
  ncol = 1,
  align = "v",
  scale = 0.95
)

# Two columns: PCA on left, PERMANOVA on right
figure1_top <- ggarrange(
  results[["all"]][["pca"]],
  right_pca,
  ncol = 2,
  widths = c(2, 1),
  align = "h"
)

right_permanova <- plot_grid(
  results[["ver"]][["permanova"]],
  results[["bel"]][["permanova"]],
  results[["uvi"]][["permanova"]],
  ncol = 1,
  align = "v",
  scale = 0.95
)

figure1_bottom <- plot_grid(
  results[["all"]][["permanova"]],
  right_permanova,
  ncol = 2,
  rel_widths = c(2, 1),
  align = "h"
)

# plot_grid(
#   plot_grid(plotlist = Fig1_pcas, ncol = 1, labels = c("(a)", "(b)", "(c)", "(d)")),
#   plot_grid(plotlist = Fig1_permanova, ncol = 1),
#   ncol = 2,
#   rel_widths = c(1, 1)
# )

# Add legend
figure1 <- ggarrange(
  figure1_top,
  NULL,
  leg,
  NULL,
  figure1_bottom,
  NULL,
  nrow = 6,
  ncol = 1,
  heights = c(6, 0.3, 0.75, 0.1, 6, 0.05),
  labels = c("(a)", "", "", "", "(b)", "")
)

ggsave(
  filename = file.path(figure_path, "Fig1_pPCA.png"),
  plot = figure1,
  width = 18, 
  height = 29, 
  units = "in",
  dpi = 150,
  type = "cairo-png"
)



########## FIGURE S3 ###################
# Variance of Principal Components for combined phenotypic space
figureS3 <- results[["all"]][["variance_plot"]]
ggsave(
  filename = file.path(figure_path, "FigS3_pAllVAR.png"),
  plot = figureS3,
  width = 8.27, 
  height = 5.22, 
  units = "in",      # inches
  dpi = 150,         # moderate dpi to reduce file size but keep quality
  type = "cairo-png" # better compression and anti-aliasing
)


########## FIGURE S4 ###################
figureS4 <- ggarrange(
  results[["all"]][["sup_pca"]],
  NULL,
  leg,
  NULL,
  # results[["all"]][["heatmap"]],
  nrow = 4,
  heights = c(8, 0.3, 1, 0.05)
)

ggsave(
  filename = file.path(figure_path, "FigS4_pPCA_all.png"),
  plot = figureS4,
  width = 12,    # A4 width in inches
  height = 16,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)


########## FIGURE S5 ###################
# Additional combined phenotypic space PCA and other locations
resFigS5 <- results[names(results) %in% c("tob", "flo", "boc")]
keep_names <- names(resFigS5)
print(keep_names)

FigS5_pcas <- lapply(resFigS5, `[[`, "pca") # extract per location pcas
FigS5_permanova <- lapply(resFigS5, `[[`, "permanova")

figureS5_top <- plot_grid(
  plot_grid(plotlist = FigS5_pcas, ncol = 1, labels = c("(a)", "(b)", "(c)")),
  plot_grid(plotlist = FigS5_permanova, ncol = 1),
  ncol = 2,
  rel_widths = c(1, 1)
)

# Add legend
figureS5 <- ggarrange(
  figureS5_top,
  NULL,
  leg,
  NULL,
  nrow = 4,
  heights = c(8, 0.3, 0.75, 0.05)
)

ggsave(
  filename = file.path(figure_path, "FigS5_pPCA_loc.png"),
  plot = figureS5,
  width = 12, 
  height = 20, 
  units = "in",
  dpi = 150,
  type = "cairo-png"
)


########## FIGURE S6 ###################
# Other PCs combination for phenotypes per location
# Remove the overall entry before extracting plots
results_no_overall <- results[!names(results) %in% c("all", "pue", "nig", "uni", "chl", "abe", "ind")]
keep_names <- names(results_no_overall)
print(keep_names)

all_sup <- lapply(results_no_overall, `[[`, "sup_pca") # extract per location supplementary pcas
sup_grid <- plot_grid(plotlist = all_sup, ncol = 2, rel_widths = c(1, 1), scale = 0.95, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")) # bundle location pcas in one plot
# Combine supplementary PCA grid with legend at the bottom
figureS6 <- ggarrange(
  sup_grid,
  NULL,
  leg,
  NULL,
  nrow = 4,
  heights = c(9, 0.3, 1, 0.05)
  )

# Save Figure S5 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS6_pPCA_sup.png"),
  plot = figureS6,
  width = 12,    # A4 width in inches
  height = 17,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)



# ########## FIGURE S8 ###################
# # Combined phenotypic space: PERMANOVA + hierarchical clustering + heatmaps
# perm <- results[["all"]][["permanova"]]
# hier <- results[["all"]][["hclust"]]
# heat <- results[["all"]][["heatmap"]]

# # # Bottom row: hier + heat
# bottom_row <- ggarrange(hier, heat, ncol = 2, labels=c("(b)","(c)"), font.label=list(color="black",size=20))

# # Combine top (perm) with bottom row
# figureS8 <- ggarrange(perm, bottom_row, nrow = 2, ncol = 1)
# figureS8 <- ggarrange(
#   NULL,
#   figureS8,
#   NULL,
#   leg,
#   NULL,
#   nrow = 5,
#   heights = c(0.2, 9, 0.2, 1, 0.05) 
# )


# # Save as PNG (A4 size)
# ggsave(
#   filename = file.path(figure_path, "FigS8_pAll.png"),
#   plot = figureS8,
#   width = 12,    # A4 width in inches
#   height = 16,  # A4 height in inches
#   units = "in",
#   dpi = 150,
#   type = "cairo-png"
# )


# ########## FIGURE S9 ###################
# # PERMANOVA heatmaps for each location
# all_perm <- lapply(results_no_overall, `[[`, "permanova") # extract per location pcas
# perm_grid <- plot_grid(plotlist = all_perm, ncol = 2, rel_widths = c(1, 1), scale=0.95) # bundle location permanovas in one plot
# # Combine permanova grid with legend at the bottom
# figureS9 <- ggarrange(
#   perm_grid
# )

# # Save Figure S7 as A4 PNG, optimized for small file size
# ggsave(filename = file.path(figure_path, "FigS9_pLocPERM.png"),
#        plot = figureS9,
#        width = 14.2,    # A4 width in inches
#        height = 17,  # A4 height in inches
#        units = "in",
#        dpi = 150,       # good quality but light (~1 MB)
#        type = "cairo-png" # smoother text rendering, smaller file
# )

# ########## FIGURE S10 ###################
# # Hierarchical clustering plots for all locations with legend
# all_hier <- lapply(results_no_overall, `[[`, "hclust") # extract per location pcas
# hier_grid <- plot_grid(plotlist = all_hier, ncol = 2, rel_widths = c(1, 1), scale=0.95) # bundle location pcas in one plot
# # Combine hierarchical clustering grid with legend at the bottom
# figureS10 <- ggarrange(
#   hier_grid,
#   NULL,
#   leg,
#   NULL,
#   nrow = 4,
#   heights = c(8, 0.2, 1, 0.05)
#   )

# # Save Figure S8 as A4 PNG, optimized for small file size
# ggsave(filename = file.path(figure_path, "FigS10_pLocHCLUST.png"),
#        plot = figureS10,
#        width = 12,    # A4 width in inches
#        height = 14,  # A4 height in inches
#        units = "in",
#        dpi = 150,       # good quality but light (~1 MB)
#        type = "cairo-png" # smoother text rendering, smaller file
# )


########## FIGURE S7 ###################
# # Heatmap PC images for each location
# all_heat <- lapply(results_no_overall, `[[`, "heatmap") # extract per location pcas
# heat_grid <- plot_grid(plotlist = all_heat, ncol = 2, rel_widths = c(1, 1), scale=0.95) # bundle location pcas in one plot
# # Combine heatmaps grid with legend at the bottom
# figureS7 <- ggarrange(
#   heat_grid
#   )

# # Save Figure S9 as A4 PNG, optimized for small file size
# ggsave(filename = file.path(figure_path, "FigS7_pLocHEAT.png"),
#        plot = figureS7,
#        width = 14.2,    # A4 width in inches
#        height = 17,  # A4 height in inches
#        units = "in",
#        dpi = 150,       # good quality but light (~1 MB)
#        type = "cairo-png" # smoother text rendering, smaller file
# )


########## FIGURE S7 ###################
# Per species phenotypic space: PCA + heatmaps + hierarchical clustering + PERMANOVA
pca_pue <- results[["pue"]][["pca"]]
heat_pue <- results[["pue"]][["heatmap"]]
perm_pue <- results[["pue"]][["permanova"]]
hier_pue <- results[["pue"]][["hclust"]]
pue <- plot_grid(pca_pue, perm_pue, ncol = 2, rel_widths = c(1, 1), align = "h", axis = "tb")


pca_nig <- results[["nig"]][["pca"]]
heat_nig <- results[["nig"]][["heatmap"]]
perm_nig <- results[["nig"]][["permanova"]]
hier_nig <- results[["nig"]][["hclust"]]
nig <- plot_grid(pca_nig, perm_nig, ncol = 2, rel_widths = c(1, 1), align = "h", axis = "tb")

pca_uni <- results[["uni"]][["pca"]]
heat_uni <- results[["uni"]][["heatmap"]]
perm_uni <- results[["uni"]][["permanova"]]
hier_uni <- results[["uni"]][["hclust"]]
uni <- plot_grid(pca_uni, perm_uni, ncol = 2, rel_widths = c(1, 1), align = "h", axis = "tb")

pca_chl <- results[["chl"]][["pca"]]
heat_chl <- results[["chl"]][["heatmap"]]
perm_chl <- results[["chl"]][["permanova"]]
hier_chl <- results[["chl"]][["hclust"]]
chl <- plot_grid(pca_chl, perm_chl, ncol = 2, rel_widths = c(1, 1), align = "h", axis = "tb")

pca_abe <- results[["abe"]][["pca"]]
heat_abe <- results[["abe"]][["heatmap"]]
perm_abe <- results[["abe"]][["permanova"]]
hier_abe <- results[["abe"]][["hclust"]]
abe <- plot_grid(pca_abe, perm_abe, ncol = 2, rel_widths = c(1, 1), align = "h", axis = "tb")

pca_ind <- results[["ind"]][["pca"]]
heat_ind <- results[["ind"]][["heatmap"]]
perm_ind <- results[["ind"]][["permanova"]]
hier_ind <- results[["ind"]][["hclust"]]
ind <- plot_grid(pca_ind, perm_ind, ncol = 2, rel_widths = c(1, 1), align = "h", axis = "tb")

# Arrange all six species
figureS7_top <- plot_grid(
  pue,
  nig,
  uni,
  chl,
  abe,
  ind,
  nrow = 6,
  rel_widths = c(1, 1, 1, 1, 1, 1),
  rel_heights = c(1, 1, 1, 1, 1, 1)
)

# Add legend
figureS7 <- ggarrange(
  figureS7_top,
  leg_g,
  nrow = 2,
  heights = c(10, 1)
)


# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS7_pSpe.png"),
  plot = figureS7,
  width = 14,    # A4 width in inches
  height = 32,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)