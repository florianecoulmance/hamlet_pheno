# by: Floriane Coulmance: 03/11/2025
# usage:
# Rscript Fig3-4_S11-S19.R 
#___________________________________________________________________

# Clear the work space
rm(list = ls())


# new libraries
source("../scripts/R/helper_functions.R")
library(smartsnp)
library(ggplot2)
library(stringi)
library(ggtext)
library(dplyr)
library(ggimage)
library(scales)
library(stringr)
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
library(viridis)
library(scico)
library(data.table)
library(hierfstat)
library(SNPRelate)
library(genoscapeRtools)
library(gridExtra)
library(purrr)
library(forcats)
library(tidytext)

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

species_info <- add_species_logos(spec_colors, logos_path)
# head(species_info)
geo_table <- read.delim(geo_colors, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# head(geo_table)

# Define your locations and PCs of interest
dataset <- list(
    arc = list(dir = "byLOC", colors = "species"),
    bar = list(dir = "byLOC", colors = "species"),
    bel = list(dir = "byLOC", colors = "species"),
    boc = list(dir = "byLOC", colors = "species"),
    flk = list(dir = "byLOC", colors = "species"),
    gun = list(dir = "byLOC", colors = "species"),
    hon = list(dir = "byLOC", colors = "species"),
    liz = list(dir = "byLOC", colors = "species"),
    pri = list(dir = "byLOC", colors = "species"),
    qui = list(dir = "byLOC", colors = "species"),
    all_s = list(dir = "byALL", colors = "species"),
    all_l = list(dir = "byALL", colors = "location"),
    abe = list(dir = "bySPC", colors = "location"),
    aff = list(dir = "bySPC", colors = "location"),
    atl = list(dir = "bySPC", colors = "location"),
    chl = list(dir = "bySPC", colors = "location"),
    gem = list(dir = "bySPC", colors = "location"),
    gum = list(dir = "bySPC", colors = "location"),
    gut = list(dir = "bySPC", colors = "location"),
    ind = list(dir = "bySPC", colors = "location"),
    nig = list(dir = "bySPC", colors = "location"),
    pue = list(dir = "bySPC", colors = "location"),
    ran = list(dir = "bySPC", colors = "location"),
    tan = list(dir = "bySPC", colors = "location"),
    uni = list(dir = "bySPC", colors = "location")
)

# Create a list to store plots per location
results <- list()

for(dat in names(dataset)) {
    dat_info <- dataset[[dat]]
    dat_dir <- dat_info$dir
    # print(dat_dir)
    color <- dat_info$colors
    # print(color)  
    message("Processing: ", dat)
  
    #-----------------------------------
    # Read GTMAT file + Sample file + PERMANOVA & PERMDISP result table
    #-----------------------------------
    gtmat_file <- file.path(base_path, "2_popgen", dat_dir, if (dat %in% c("all_s", "all_l")) "all.agg.ld_pruned_gtmat.traw" else paste0(dat, "_ld_pruned_gtmat.traw"))
    sample_file <- file.path(base_path, if (dat %in% c("all_s", "all_l")) "metadata/geno_names.txt" else paste0("2_popgen/", dat_dir, "/", dat, ".txt"))
    # perm_file <- list.files(file.path(base_path, "2_popgen", dat_dir, "permanova_results"), pattern = paste0(dat, ".*\\.csv$"), full.names = TRUE)
    perm_file <- if (dat == "all_s") file.path(base_path, "2_popgen", dat_dir, "permanova_results/all.lm.pairwise.csv") else if (dat == "all_l") file.path(base_path, "2_popgen", dat_dir, "permanova_results/all.sm.pairwise.csv") else list.files(file.path(base_path, "2_popgen", dat_dir, "permanova_results"), pattern = paste0(dat, ".*\\.csv$"), full.names = TRUE)

    # print(gtmat_file)
    # print(sample_file)
    # print(perm_file)

    #-----------------------------------
    # PCA
    #-----------------------------------
    pca_res <- pca_analysis(gtmat_file, sample_file, color_by = color)
    pca_eigen <- pca_res$eigen
    pca_var   <- pca_res$var
    # print(pca_res)
    # print(pca_eigen)
    # print(pca_var)

    #-----------------------------------
    # PCA plot
    #-----------------------------------
    if (dat %in% c("all_s", "all_l")) {
        # For the "all" dataset, generate three PCA plots (PC1-2, PC3-4, PC5-6)
        p_pca1 <- pca_plot_all(pca_eigen, "PC1", "PC2", species_info, pca_var) %>% annotate_figure(., top = text_grob("(a)", color = "black", face = "bold", size = 20, x = unit(0, "lines"), vjust=0, hjust=0))
        p_pca2 <- pca_plot_all(pca_eigen, "PC3", "PC4", species_info, pca_var) %>% annotate_figure(., top = text_grob("(b)", color = "black", face = "bold", size = 20, x = unit(0, "lines"), vjust=0, hjust=0))
        p_pca3 <- pca_plot_all(pca_eigen, "PC5", "PC6", species_info, pca_var) %>% annotate_figure(., top = text_grob("(c)", color = "black", face = "bold", size = 20, x = unit(0, "lines"), vjust=0, hjust=0))
    } else {
        # For all other datasets, generate only PC1-2 and PC3-4
        p_pca1 <- pca_plot(pca_eigen, "PC1", "PC2", species_info, geo_table, pca_var, color_by = color)
        p_pca2 <- pca_plot(pca_eigen, "PC3", "PC4", species_info, geo_table, pca_var, color_by = color)
        p_pca3 <- NULL
    }

    #-----------------------------------
    # VAR plot
    #-----------------------------------
    pca_var_df <- data.frame(PC=as.numeric(rownames(pca_var)), Variance = pca_var$X0)
    # print(pca_var_df)
    p_var <- plot_variance(pca_var_df, dat)

    # -----------------------------------
    # PERMANOVA + PERMDISP (filter <5 inds per species inside perm_f)
    #-----------------------------------
    # p_perm <- plot_permanova_permdisp(perm_file, species_info, geo_table, color_by = color, params_legend = if(dat %in% c("bel", "uni")) c(0.4, 0.8) else if(dat %in% c("all_s", "all_l")) c(0.3, 0.7) else "none")


    # -----------------------------------
    # FST (filter <3 inds per species)
    #-----------------------------------
    # p_fst <- fst_analysis(gtmat_file, color_by = color, species_info, geo_table, dat)

    #-----------------------------------
    # Store outputs
    #-----------------------------------
    results[[dat]] <- list(
      pca_f = p_pca1,
      pca_s = p_pca2,
      pca_t = p_pca3,
      variance_plot = p_var #,
      # permanova = p_perm
    )
}


# # ############################
# # FINAL PLOTS
# # ############################
# # Set datasets for plots
results_locations <- results[names(results) %in% c("hon", "bel", "boc", "pri")]
keep_names <- names(results_locations)
# print(keep_names)

# # Create common legend
leg <- legend_plot(species_info, gen = TRUE)
leg_g <- legend_geo(geo_table, gen = TRUE)

########## FIGURE 2 ###################
all_pcas <- lapply(results_locations, `[[`, "pca_f") # extract per location pcas
pca_grid <- plot_grid(plotlist = all_pcas, ncol = 4, rel_widths = c(1, 1, 1, 1), scale = 0.95) # bundle location pcas in one plot

# 1. READ PIXY OUTPUT
fst <- read.table(
  file.path(base_path, "2_popgen", "byALL", "all.flt_fst.min3.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
head(fst)

# dxy <- read.table(
#   file.path(base_path, "2_popgen", "byALL", "all.flt_dxy.min3.txt"),
#   header = TRUE,
#   sep = "\t",
#   stringsAsFactors = FALSE
# )
# head(dxy)

# # 2. AGGREGATE FST AND DXY BY WINDOW
# fst_window <- aggregate_fst_by_window(fst)
# head(fst_window)
# dxy_window <- aggregate_dxy_by_window(dxy)
# head(dxy_window)


# # 3. CREATE GENOMIC COORDINATES
# fst_genome <- add_genome_position(fst_window)
# head(fst_genome)
# fst_window <- fst_genome$data
# head(fst_window)


# chromosome_info <- fst_genome$chromosomes
# head(chromosome_info)

# dxy_genome <- add_genome_position(dxy_window)
# head(dxy_genome)
# dxy_window <- dxy_genome$data
# head(dxy_window)

fst_plot <- fst %>%
  select(pop1, pop2, avg_wc_fst)
head(fst_plot)
# dxy_plot <- dxy %>%
#   select(pop1, pop2, avg_dxy)
# head(dxy_plot)


# # 4. PLOT A — FST
# pA <- plot_genome_fst(
#   fst_window,
#   chromosome_info
# )

# # 5. PLOT B — DXY
# pB <- plot_genome_dxy(
#   dxy_window,
#   chromosome_info
# )

print(plot_pairwise_metric)
typeof(plot_pairwise_metric)
class(plot_pairwise_metric)

location_colors <- setNames(
  geo_table$Color,
  geo_table$geo
)

print("BEFORE pC")
# 6. PLOT C — FST boxplots
pC <- plot_pairwise_metric(
  fst_plot,
  metric = "avg_wc_fst",
  xlab = "Pairwise FST",
  location_colors = location_colors
)
print("AFTER pC")

# # 7. PLOT D — DXY boxplots
# print("BEFORE pD")
# pD <- plot_pairwise_metric(
#   dxy_plot,
#   metric = "avg_dxy",
#   xlab = "Pairwise DXY",
#   location_colors = location_colors
# )
# print("AFTER pD")

# Combine PCA grid with legend at the bottom
figure2 <- ggarrange(
  pca_grid,
  NULL,
  pC,
  NULL,
  leg,
  NULL,
  nrow = 6, 
  heights = c(8, 0.3, 8, 0.3, 1, 0.05)
) # adjust if legend is too big/small

ggsave(
  filename = file.path(figure_path, "Fig2_pairFST.png"),
  plot = pC,
  width = 8,
  height = 14,
  units = "in",
  dpi = 150,
  type = "cairo-png"
)


# ########## FIGURE S14 ###################
figureS14 <- plot_fst_categories(
  df = fst
)

ggsave(
  filename = file.path(figure_path, "FigS14_gFSTviolin.png"),
  plot = figureS14,
  width = 6.5,
  height = 6.5
)


# ########## FIGURE S15 ###################
# # Combined genetic space with legend
# # Access each PCA plot for the "all" dataset
# pca1 <- results[["all_s"]][["pca_f"]]
# pca2 <- results[["all_s"]][["pca_s"]]
# pca3 <- results[["all_s"]][["pca_t"]]

# plot <- ggarrange(
#   pca1,
#   pca2,
#   pca3,
#   nrow = 3,
#   ncol = 1,
#   common.legend=T,
#   legend = "right"
#   )

# figureS15 <- ggarrange(
#   NULL,
#   plot,
#   nrow = 2,
#   ncol = 1,
#   heights = c(0.02, 15)
#   )

# ggsave(
#   filename = file.path(figure_path, "FigS15_gAllPCA.png"),
#   plot = figureS15,
#   width = 6, 
#   height = 16, 
#   units = "in",      # inches
#   dpi = 150,         # moderate dpi to reduce file size but keep quality
#   type = "cairo-png" # better compression and anti-aliasing
# )


########## FIGURE S16 ###################
# Variance of Principal Components for combined genetic space
figureS16 <- results[["all_s"]][["variance_plot"]]
ggsave(
  filename = file.path(figure_path, "FigS16_gAllVAR.png"),
  plot = figureS16,
  width = 8.27, 
  height = 5.22, 
  units = "in",      # inches
  dpi = 150,         # moderate dpi to reduce file size but keep quality
  type = "cairo-png" # better compression and anti-aliasing
)


# ########## FIGURE S17 ###################
# # PCA plots for all locations with legend
# all_pcas <- lapply(results_locations, `[[`, "pca_f") # extract per location pcas
# pca_grid <- plot_grid(plotlist = all_pcas, ncol = 2, rel_widths = c(1, 1), scale = 0.95) # bundle location pcas in one plot
# # Combine PCA grid with legend at the bottom
# figureS17 <- ggarrange(
#   pca_grid,
#   NULL,
#   leg,
#   NULL,
#   nrow = 4, 
#   heights = c(8, 0.3, 1, 0.05)
# ) # adjust if legend is too big/small

# # Save Figure S13 as A4 PNG, optimized for small file size
# ggsave(filename = file.path(figure_path, "FigS17_gLocPCA.png"),
#   plot = figureS17,
#   width = 12,    # A4 width in inches
#   height = 14,  # A4 height in inches
#   units = "in",
#   dpi = 150,       # good quality but light (~1 MB)
#   type = "cairo-png" # smoother text rendering, smaller file
# )


# ########## FIGURE S18 ###################
# # Other PCs combination for genotypes per location
# all_sup <- lapply(results_locations, `[[`, "pca_s") # extract per location pcas
# sup_grid <- plot_grid(plotlist = all_sup, ncol = 2, rel_widths = c(1, 1), scale = 0.95) # bundle location pcas in one plot
# # Combine PCA grid with legend at the bottom
# figureS18 <- ggarrange(
#   sup_grid,
#   NULL,
#   leg,
#   NULL,
#   nrow = 4,
#   heights = c(8, 0.3, 1, 0.05)
# ) # adjust if legend is too big/small

# # Save Figure S12 as A4 PNG, optimized for small file size
# ggsave(filename = file.path(figure_path, "FigS18_gLocSUP.png"),
#   plot = figureS18,
#   width = 12,    # A4 width in inches
#   height = 14,  # A4 height in inches
#   units = "in",
#   dpi = 150,       # good quality but light (~1 MB)
#   type = "cairo-png" # smoother text rendering, smaller file
# )


# ########## FIGURE S19 ###################
# # Combined genotypic space: PERMANOVA
# figureS19 <- results[["all_s"]][["permanova"]]

# # Save as PNG (A4 size)
# ggsave(
#   filename = file.path(figure_path, "FigS19_gAllPERM.png"),
#   plot = figureS19,
#   width = 7.5,    # A4 width in inches
#   height = 7.5,  # A4 height in inches
#   units = "in",
#   dpi = 150,
#   type = "cairo-png"
# )


# ########## FIGURE S20 ###################
# # PERMANOVA heatmaps for each location
# all_perm <- lapply(results_locations, `[[`, "permanova") # extract per location pcas
# figureS20 <- plot_grid(plotlist = all_perm, ncol = 2, rel_widths = c(1, 1), scale = 0.95)# bundle location pcas in one plot

# # Save Figure S16 as A4 PNG, optimized for small file size
# ggsave(filename = file.path(figure_path, "FigS20_gLocPERM.png"),
#        plot = figureS20,
#        width = 10,    # A4 width in inches
#        height = 10,  # A4 height in inches
#        units = "in",
#        dpi = 150,       # good quality but light (~1 MB)
#        type = "cairo-png" # smoother text rendering, smaller file
# )


# ########## FIGURE S21 ###################
# # Per species genotypic space: PCA + PERMANOVA + PERMDISP
# # Remove the overall entry before extracting plots
# results_spc <- results[names(results) %in% c("pue", "nig", "uni")]
# keep_spc <- names(results_spc)
# # print(keep_spc)

# pca_pue_f <- results[["pue"]][["pca_f"]]
# pca_pue_s <- results[["pue"]][["pca_s"]] %>% annotate_figure(., top=NULL)
# perm_pue <- results[["pue"]][["permanova"]]
# pue <- plot_grid(pca_pue_f, pca_pue_s, perm_pue, ncol = 3, rel_widths = c(1, 1, 1), scale=0.95)


# pca_nig_f <- results[["nig"]][["pca_f"]]
# pca_nig_s <- results[["nig"]][["pca_s"]]
# perm_nig <- results[["nig"]][["permanova"]]
# nig <- plot_grid(pca_nig_f, pca_nig_s, perm_nig, ncol = 3, rel_widths = c(1, 1, 1), scale=0.95)

# pca_uni_f <- results[["uni"]][["pca_f"]]
# pca_uni_s <- results[["uni"]][["pca_s"]]
# perm_uni <- results[["uni"]][["permanova"]]
# uni <- plot_grid(pca_uni_f, pca_uni_s, perm_uni, ncol = 3, rel_widths = c(1, 1, 1), scale=0.95)

# figureS21 <- plot_grid(
#   pue,
#   NULL,
#   nig,
#   NULL,
#   uni,
#   NULL,
#   leg_g,
#   NULL,
#   ncol = 1,
#   nrow = 8,
#   rel_heights = c(6, 0.2, 6, 0.2, 6, 0.2, 2, 0.05),
#   rel_widths = c(20, 20, 20, 20, 20, 20, 16, 20)
#   )

# # Save as PNG (A4 size)
# ggsave(
#   filename = file.path(figure_path, "FigS21_gSpe.png"),
#   plot = figureS21,
#   width = 15,    # A4 width in inches
#   height = 18,  # A4 height in inches
#   units = "in",
#   dpi = 150,
#   type = "cairo-png"
# )