# by: Floriane Coulmance: 03/11/2025
# usage:
# Rscript Fig3-4_S11-S15.R 
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
FST_RESULTS <- list()

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
    gtmat_file <- file.path(base_path, "2_popgen", dat_dir, if (dat %in% c("all_s", "all_l")) "thinned_all_ld_pruned_gtmat.traw" else paste0(dat, "_ld_pruned_gtmat.traw"))
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
        p_pca1 <- pca_plot_all(pca_eigen, "PC1", "PC2", species_info, pca_var)
        p_pca2 <- pca_plot_all(pca_eigen, "PC3", "PC4", species_info, pca_var)
        p_pca3 <- pca_plot_all(pca_eigen, "PC5", "PC6", species_info, pca_var)
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
    p_perm <- plot_permanova_permdisp(perm_file, species_info, geo_table, color_by = color, params_legend = if(dat %in% c("bel", "nig")) c(0.3, 0.7) else if(dat %in% c("all_s", "all_l")) c(0.2, 0.7) else "none")


    # -----------------------------------
    # FST (filter <3 inds per species)
    #-----------------------------------
    p_fst <- fst_analysis(gtmat_file, color_by = color, species_info, geo_table, dat)

    #-----------------------------------
    # Store outputs
    #-----------------------------------
    results[[dat]] <- list(
      pca_f = p_pca1,
      pca_s = p_pca2,
      pca_t = p_pca3,
      variance_plot = p_var,
      permanova = p_perm,
      fst = p_fst
    )
}

print(names(FST_RESULTS))


# ############################
# FINAL PLOTS
# ############################

# Set datasets for plots
results_locations <- results[names(results) %in% c("hon", "bel", "boc", "pri")]
keep_names <- names(results_locations)
# print(keep_names)

# Create common legend
leg <- legend_plot(species_info, gen = TRUE)
leg_g <- legend_geo(geo_table, gen = TRUE)

########## FIGURE 3 ###################
# PCA plots for all locations with legend
all_pcas <- lapply(results_locations, `[[`, "pca_f") # extract per location pcas
pca_grid <- plot_grid(plotlist = all_pcas, ncol = 2, labels=c("(a)","(b)","(c)","(d)"), label_size=30) # bundle location pcas in one plot
# Combine PCA grid with legend at the bottom
figure3 <- ggarrange(
  pca_grid,
  NULL,
  leg,
  NULL,
  nrow = 4, 
  heights = c(8, 0.3, 1, 0.05)
) # adjust if legend is too big/small

# Save Figure 3 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "Fig3_gLocPCA.png"),
  plot = figure3,
  width = 12,    # A4 width in inches
  height = 14,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)


########## FIGURE 4 ###################
# Combined genetic space with legend
# Access each PCA plot for the "all" dataset
pca1 <- results[["all_s"]][["pca_f"]]
pca2 <- results[["all_s"]][["pca_s"]]
pca3 <- results[["all_s"]][["pca_t"]]

figure4 <- ggarrange(
  pca1,
  pca2,
  pca3,
  nrow = 3,
  ncol = 1,
  labels = c("(a)", "(b)", "(c)"),
  font.label=list(color="black",size=30),
  common.legend=T,
  legend = "right"
  )

ggsave(
  filename = file.path(figure_path, "Fig4_gAllPCA.png"),
  plot = figure4,
  width = 8, 
  height = 16, 
  units = "in",      # inches
  dpi = 150,         # moderate dpi to reduce file size but keep quality
  type = "cairo-png" # better compression and anti-aliasing
)

########## FIGURE S11 ###################
# Variance of Principal Components for combined genetic space
figureS11 <- results[["all_s"]][["variance_plot"]]
ggsave(
  filename = file.path(figure_path, "FigS11_gAllVAR.png"),
  plot = figureS11,
  width = 8.27, 
  height = 5.22, 
  units = "in",      # inches
  dpi = 150,         # moderate dpi to reduce file size but keep quality
  type = "cairo-png" # better compression and anti-aliasing
)

########## FIGURE S12 ###################
# Other PCs combination for genotypes per location
all_sup <- lapply(results_locations, `[[`, "pca_s") # extract per location pcas
sup_grid <- plot_grid(plotlist = all_sup, ncol = 2, labels = c("(a)", "(b)", "(c)", "(d)"), label_size=30) # bundle location pcas in one plot
# Combine PCA grid with legend at the bottom
figureS12 <- ggarrange(
  sup_grid,
  NULL,
  leg,
  NULL,
  nrow = 4,
  heights = c(8, 0.3, 1, 0.05)
) # adjust if legend is too big/small

# Save Figure S12 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS12_gLocSUP.png"),
  plot = figureS12,
  width = 12,    # A4 width in inches
  height = 14,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S13 ###################
# PERMANOVA heatmaps for each location
all_perm <- lapply(results_locations, `[[`, "permanova") # extract per location pcas
figureS13 <- plot_grid(plotlist = all_perm, ncol = 2, labels = c("(a)", "(b)", "(c)", "(d)"), label_size=30) # bundle location pcas in one plot

# Save Figure S13 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS13_gLocPERM.png"),
       plot = figureS13,
       width = 10,    # A4 width in inches
       height = 10,  # A4 height in inches
       units = "in",
       dpi = 150,       # good quality but light (~1 MB)
       type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S14 ###################
# Combined genotypic space: PERMANOVA
figureS14 <- results[["all_s"]][["permanova"]]

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS14_gAllPERM.png"),
  plot = figureS14,
  width = 12,    # A4 width in inches
  height = 12,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)

########## FIGURE S15 ###################
# Per species genotypic space: PCA + PERMANOVA + PERMDISP
# Remove the overall entry before extracting plots
results_spc <- results[names(results) %in% c("pue", "nig", "uni", "pri")]
keep_spc <- names(results_spc)
# print(keep_spc)

pca_pue_f <- results[["pue"]][["pca_f"]]
pca_pue_s <- results[["pue"]][["pca_s"]]
perm_pue <- results[["pue"]][["permanova"]]
pue <- plot_grid(pca_pue_f, pca_pue_s, perm_pue, ncol = 3, rel_widths = c(1, 1, 1))


pca_nig_f <- results[["nig"]][["pca_f"]]
pca_nig_s <- results[["nig"]][["pca_s"]]
perm_nig <- results[["nig"]][["permanova"]]
nig <- plot_grid(pca_nig_f, pca_nig_s, perm_nig, ncol = 3, rel_widths = c(1, 1, 1))

pca_uni_f <- results[["uni"]][["pca_f"]]
pca_uni_s <- results[["uni"]][["pca_s"]]
perm_uni <- results[["uni"]][["permanova"]]
uni <- plot_grid(pca_uni_f, pca_uni_s, perm_uni, ncol = 3, rel_widths = c(1, 1, 1))

figureS15 <- plot_grid(
  pue,
  nig,
  uni,
  NULL,
  leg_g,
  ncol = 1,
  nrow = 5,
  rel_heights = c(5, 5, 5, 0.5, 1),
  rel_widths = c(20, 20, 20, 20, 16),
  labels = c("(a)", "(b)", "(c)", ""),
  label_size=30
  )

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS15_gSpe.png"),
  plot = figureS15,
  width = 20,    # A4 width in inches
  height = 16,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)


########## FIGURE S16 ###################
# FST boxplots
fst_species <- FST_RESULTS[names(FST_RESULTS) %in% "all_s"]
fst_symp <- FST_RESULTS[names(FST_RESULTS) %in% c("hon", "bel", "boc", "pri", "arc", "bar", "flk", "gun", "qui")]
fst_locations <- FST_RESULTS[names(FST_RESULTS) %in% c("atl", "pue", "nig", "uni", "abe", "aff", "chl", "gem", "gum", "ind", "tan")]

# for (nm in names(fst_species)) {
#   print(nm)
#   print(fst_species[[nm]])
# }

fst_species_df <- map_df(
  names(fst_species),
  ~ fst_species[[.x]] %>% mutate(dataset = .x)
)
print(fst_species_df)

fst_symp_df <- map_df(
  names(fst_symp),
  ~ fst_symp[[.x]] %>% mutate(dataset = .x)
)
print(fst_symp_df)

fst_locations_df <- map_df(
  names(fst_locations),
  ~ fst_locations[[.x]] %>% mutate(dataset = .x)
)
print(fst_locations_df)

# Combine everything into one long df
fst_all_df <- bind_rows(
  fst_species_df,
  fst_symp_df,
  fst_locations_df
)

fst_all_df$dataset <- factor(
  fst_all_df$dataset,
  levels = c(names(fst_species), names(fst_symp), names(fst_locations))
)

fst_all_df <- fst_all_df %>%
  mutate(group =
    case_when(
      dataset %in% names(fst_species) ~ "between species",
      dataset %in% names(fst_symp) ~ "between sympatric species",
      dataset %in% names(fst_locations) ~ "between locations per species"
    )
  )

print(fst_all_df)

figureS16 <- ggplot(fst_all_df, aes(x = dataset, y = Fst, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "transparent") +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1) +
  theme_bw() +
  labs(
    x = "",
    y = "Pairwise FST"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS16_gFSTviolin.png"),
  plot = figureS16,
  width = 10,    # A4 width in inches
  height = 7,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)


########## FIGURE S17 ###################
# Per location FST
figureS17 <- ggarrange(
  results[["bel"]][["fst"]],
  results[["boc"]][["fst"]],
  results[["hon"]][["fst"]],
  results[["pri"]][["fst"]],
  results[["arc"]][["fst"]],
  results[["bar"]][["fst"]],
  results[["flk"]][["fst"]],
  results[["gun"]][["fst"]],
  results[["qui"]][["fst"]],
  ncol = 3,
  nrow = 3,
  labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"),
  font.label=list(color="black",size=30)
)

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS17_gLocFST.png"),
  plot = figureS17,
  width = 14,    # A4 width in inches
  height = 12,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)


########## FIGURE S18 ###################
# Combined genotypic space: FST

figureS18 <- results[["all_s"]][["fst"]]

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS18_gAllsFST.png"),
  plot = figureS18,
  width = 15,    # A4 width in inches
  height = 15,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)


########## FIGURE S19 ###################
# Per species: FST
figureS19 <- ggarrange(
  results[["pue"]][["fst"]],
  results[["nig"]][["fst"]],
  results[["uni"]][["fst"]],
  results[["abe"]][["fst"]],
  results[["aff"]][["fst"]],
  results[["chl"]][["fst"]],
  results[["gem"]][["fst"]],
  results[["gum"]][["fst"]],
  results[["ind"]][["fst"]],
  results[["tan"]][["fst"]],
  ncol = 2,
  nrow = 5,
  labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)"),
  font.label=list(color="black",size=30)
)

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS19_gSpeFST.png"),
  plot = figureS19,
  width = 14,    # A4 width in inches
  height = 19.5,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)


# ########## FIGURE S19 ###################
# # Combined genotypic space: FST

# figureS19 <- results[["all_l"]][["fst"]]

# # Save as PNG (A4 size)
# ggsave(
#   filename = file.path(figure_path, "FigS19_gAlllFST.png"),
#   plot = figureS19,
#   width = 15,    # A4 width in inches
#   height = 15,  # A4 height in inches
#   units = "in",
#   dpi = 150,
#   type = "cairo-png"
# )