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
  boc = list(name = "lab_boc229_left_noflash", pcs = c("PC1", "PC3", "PC2", "PC4"), colors = "species"),
  uvi = list(name = "lab_uvi136_left_noflash", pcs = c("PC1", "PC4", "PC2", "PC3"), colors = "species"),
  tob = list(name = "lab_tob69_left_noflash", pcs = c("PC1", "PC3", "PC2", "PC4"), colors = "species"),
  ver = list(name = "lab_ver60_left_noflash", pcs = c("PC1", "PC2", "PC3", "PC4"), colors = "species"),
  bel = list(name = "lab_bel46_left_noflash", pcs = c("PC1", "PC2", "PC3", "PC4"), colors = "species"),
  flo = list(name = "lab_flo29_left_noflash", pcs = c("PC1", "PC4", "PC2", "PC3"), colors = "species"),
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
    heatmap = p_heat
  )
}




# ############################
# FINAL PLOTS
# ############################

# Remove the overall entry before extracting plots
results_no_overall <- results[!names(results) %in% c("all", "pue", "nig", "uni", "chl", "abe", "ind")]
keep_names <- names(results_no_overall)
print(keep_names)

# Create common legend
leg <- legend_plot(species_info)

# Combine all data for legend
# df_all <- do.call(rbind, lapply(results[!names(results) %in% c("all", "pue", "nig", "uni", "chl", "abe", "ind")], function(x) x$data))
# var_all <- do.call(rbind, lapply(results[!names(results) %in% c("all", "pue", "nig", "uni", "chl", "abe", "ind")], function(x) x$variance))

# plot <- pca_plot(results[["all"]][["data"]], "PC1", "PC2", species_info, geo_table, results[["all"]][["variance"]], color_by = "species", extract_legend = TRUE)
# print("Plot for all pheno: ")
# print(class(plot))
# leg_combined <- get_legend(plot)
# print("Legend for all : ")
# print(leg_combined)
# print(class(leg_combined))

########## FIGURE 1 ###################
# PCA plots for all locations with legend
all_pcas <- lapply(results_no_overall, `[[`, "pca") # extract per location pcas
pca_grid <- plot_grid(plotlist = all_pcas, ncol = 2, labels=c("(a)","(b)","(c)","(d)","(e)","(f)")) # bundle location pcas in one plot
# Combine PCA grid with legend at the bottom
figure1 <- ggarrange(
  pca_grid,
  NULL,
  leg,
  nrow = 3,
  heights = c(10, 0.3, 1)#, widths = c (3,1)
  )

# Save Figure 1 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "Fig1_pLocPCA.png"),
  plot = figure1,
  width = 18.5,    # A4 width in inches
  height = 25,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)


########## FIGURE 2 ###################
# Combined phenotypic space with legend
figure2 <- ggarrange(
  results[["all"]][["pca"]],
  NULL,
  leg,
  nrow = 3,
  heights = c(10, 0.3, 1)
  )
  
ggsave(
  filename = file.path(figure_path, "Fig2_pAllPCA.png"),
  plot = figure2,
  width = 18, 
  height = 19, 
  units = "in",      # inches
  dpi = 150,         # moderate dpi to reduce file size but keep quality
  type = "cairo-png" # better compression and anti-aliasing
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
# Other PCs combination for phenotypes per location
all_sup <- lapply(results_no_overall, `[[`, "sup_pca") # extract per location supplementary pcas
sup_grid <- plot_grid(plotlist = all_sup, ncol = 2, labels=c("(a)","(b)","(c)","(d)","(e)","(f)")) # bundle location pcas in one plot
# Combine supplementary PCA grid with legend at the bottom
figureS4 <- ggarrange(
  sup_grid,
  NULL,
  leg,
  nrow = 3,
  heights = c(10, 0.3, 1)
  )

# Save Figure S4 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS4_pLocSUP.png"),
  plot = figureS4,
  width = 18.5,    # A4 width in inches
  height = 25,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S5 ###################
# PERMANOVA heatmaps for each location
all_perm <- lapply(results_no_overall, `[[`, "permanova") # extract per location pcas
perm_grid <- plot_grid(plotlist = all_perm, ncol = 2, labels=c("(a)","(b)","(c)","(d)","(e)","(f)")) # bundle location permanovas in one plot
# Combine permanova grid with legend at the bottom
figureS5 <- ggarrange(
  perm_grid
)

# Save Figure S5 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS5_pLocPERM.png"),
       plot = figureS5,
       width = 14.2,    # A4 width in inches
       height = 17,  # A4 height in inches
       units = "in",
       dpi = 150,       # good quality but light (~1 MB)
       type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S6 ###################
# Hierarchical clustering plots for all locations with legend
all_hier <- lapply(results_no_overall, `[[`, "hclust") # extract per location pcas
hier_grid <- plot_grid(plotlist = all_hier, ncol = 2, labels=c("(a)","(b)","(c)","(d)","(e)","(f)")) # bundle location pcas in one plot
# Combine hierarchical clustering grid with legend at the bottom
figureS6 <- ggarrange(
  hier_grid,
  NULL,
  leg,
  nrow = 3,
  heights = c(10, 0.3, 1)
  )

# Save Figure S6 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS6_pLocHCLUST.png"),
       plot = figureS6,
       width = 14.2,    # A4 width in inches
       height = 20,  # A4 height in inches
       units = "in",
       dpi = 150,       # good quality but light (~1 MB)
       type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S7 ###################
# Heatmap PC images for each location
all_heat <- lapply(results_no_overall, `[[`, "heatmap") # extract per location pcas
heat_grid <- plot_grid(plotlist = all_heat, ncol = 2, labels=c("(a)","(b)","(c)","(d)","(e)","(f)")) # bundle location pcas in one plot
# Combine heatmaps grid with legend at the bottom
figureS7 <- ggarrange(
  heat_grid
  )

# Save Figure S7 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS7_pLocHEAT.png"),
       plot = figureS7,
       width = 14.2,    # A4 width in inches
       height = 17,  # A4 height in inches
       units = "in",
       dpi = 150,       # good quality but light (~1 MB)
       type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S8 ###################
# Combined phenotypic space: supplementary PCA + PERMANOVA + hierarchical clustering + heatmaps
perm <- results[["all"]][["permanova"]]
hier <- results[["all"]][["hclust"]]
heat <- results[["all"]][["heatmap"]]

# # Bottom row: hier + heat
bottom_row <- ggarrange(hier, heat, ncol = 2, labels=c("(b)","(c)"))

# Combine top (perm) with bottom row
figureS8 <- ggarrange(perm, bottom_row, nrow = 2, ncol = 1, labels=c("(a)",""))
figureS8 <- ggarrange(
  figureS8,
  NULL,
  leg,
  nrow = 3,
  heights = c(10, 0.3, 1) 
)


# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS8_pAll.png"),
  plot = figureS8,
  width = 14.2,    # A4 width in inches
  height = 19,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)


########## FIGURE S9 ###################
# Additional combined phenotypic space PCA
figureS9 <- ggarrange(
  results[["all"]][["sup_pca"]],
  NULL,
  leg,
  nrow = 3,
  heights = c(10, 0.3, 1)
)

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS9_pAllSUP.png"),
  plot = figureS9,
  width = 18,    # A4 width in inches
  height = 19,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)

########## FIGURE S10 ###################
# Per species phenotypic space: PCA + heatmaps + hierarchical clustering + PERMANOVA
pca_pue <- results[["pue"]][["pca"]]
leg_geo <- get_legend(results[["pue"]][["pca"]])
heat_pue <- results[["pue"]][["heatmap"]]
perm_pue <- results[["pue"]][["permanova"]]
hier_pue <- results[["pue"]][["hclust"]]
pue <- plot_grid(pca_pue, heat_pue, perm_pue, hier_pue, ncol = 4, rel_widths = c(1, 1, 1, 1))

pca_nig <- results[["nig"]][["pca"]]
heat_nig <- results[["nig"]][["heatmap"]]
perm_nig <- results[["nig"]][["permanova"]]
hier_nig <- results[["nig"]][["hclust"]]
nig <- plot_grid(pca_nig, heat_nig, perm_nig, hier_nig, ncol = 4, rel_widths = c(1, 1, 1, 1))

pca_uni <- results[["uni"]][["pca"]]
heat_uni <- results[["uni"]][["heatmap"]]
perm_uni <- results[["uni"]][["permanova"]]
hier_uni <- results[["uni"]][["hclust"]]
uni <- plot_grid(pca_uni, heat_uni, perm_uni, hier_uni, ncol = 4, rel_widths = c(1, 1, 1, 1))

pca_chl <- results[["chl"]][["pca"]]
heat_chl <- results[["chl"]][["heatmap"]]
perm_chl <- results[["chl"]][["permanova"]]
hier_chl <- results[["chl"]][["hclust"]]
chl <- plot_grid(pca_uni, heat_uni, perm_uni, hier_uni, ncol = 4, rel_widths = c(1, 1, 1, 1))

pca_abe <- results[["abe"]][["pca"]]
heat_abe <- results[["abe"]][["heatmap"]]
perm_abe <- results[["abe"]][["permanova"]]
hier_abe <- results[["abe"]][["hclust"]]
abe <- plot_grid(pca_abe, heat_abe, perm_abe, hier_abe, ncol = 4, rel_widths = c(1, 1, 1, 1))

pca_ind <- results[["ind"]][["pca"]]
heat_ind <- results[["ind"]][["heatmap"]]
perm_ind <- results[["ind"]][["permanova"]]
hier_ind <- results[["ind"]][["hclust"]]
ind <- plot_grid(pca_ind, heat_ind, perm_ind, hier_ind, ncol = 4, rel_widths = c(1, 1, 1, 1))

figureS10 <- ggarrange(pue, nig, uni, chl, abe, ind, nrow = 6, labels=c("(a)","(b)","(c)","(d)","(e)","(f)"))
figureS10 <- ggarrange(
  figureS10,
  leg_geo,
  nrow = 2,
  heights = c(10, 0.5)
  )

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS10_pSpe.png"),
  plot = figureS10,
  width = 24,    # A4 width in inches
  height = 40,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)