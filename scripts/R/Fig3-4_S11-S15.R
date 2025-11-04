# by: Floriane Coulmance: 03/11/2025
# usage:
# Rscript Fig3-4_S11-S15.R 
#___________________________________________________________________

# Clear the work space
rm(list = ls())


# new libraries
source("Fig1-2_S3-S10.R")
library(smartsnp)
library(ggplot2)
library(stringi)
library(ggtext)
library(dplyr)
library(ggimage)


# ############################
# FUNCTIONS
# ############################

# ============================================================
# Function: pca_analysis
# Purpose : Run SmartSNP PCA on genotype data and extract key PCA components.
# Input   :
#   - gtfile     : Path to the PLINK .traw genotype file.
#   - samplefile : Path to the text file containing sample names (one per line).
#   - color_by   : Defines grouping mode:
#                  "species" → use last 3 characters of sample names.
#                  "location" → use characters -5 to -3 of sample names.
# Output  :
#   - A list containing:
#       $eigen : Data frame of PCA sample coordinates.
#       $var   : Numeric vector of eigenvalues (variance explained).
# Notes   :
#   - The function automatically determines grouping mode from color_by.
#   - If there is only one unique group, the function exits with a warning.
# ============================================================
pca_analysis <- function(gtfile, samplefile, color_by) {
    # Determine mode based on color_by
    speciesMODE <- color_by == "species"
    
    # Read sample names and define groups
    allsamples <- readLines(samplefile)
    # Extract spec and geo depending on mode
    if (speciesMODE) {
        spec <- substr(allsamples, nchar(allsamples) - 2, nchar(allsamples))
        geo  <- substr(allsamples, nchar(allsamples) - 5, nchar(allsamples) - 3)
    } else {
        geo  <- substr(allsamples, nchar(allsamples) - 2, nchar(allsamples))
        spec <- substr(allsamples, nchar(allsamples) - 5, nchar(allsamples) - 3)
    }
    
    # Check there is enough variation
    if (length(unique(group_id)) <= 1) {
        message("Not enough different samples")
        quit(save = "no", status = 1)
    }
    
    # Run SmartPCA
    sm.pca <- smart_pca(
        snp_data = gtfile,
        sample_group = if (speciesMODE) spec else geo,
        missing_value = NA
    )
    
    # Extract PCA outputs
    pca_coords <- sm.pca$pca.sample_coordinates
    var_explained <- sm.pca$pca.eigenvalues
    
    # Ensure the PCA coordinates rows are in the same order as the original samples
    pca_coords <- pca_coords[match(allsamples, rownames(pca_coords)), ]
    
    # Then add spec, geo, and sample
    pca_coords$spec   <- spec
    pca_coords$geo    <- geo
    pca_coords$sample <- allsamples
    
    # Return both
    return(list(eigen = pca_coords, var = var_explained))
}


# ============================================================
# Function: pca_plot_all
# Purpose : Plot PCA for the entire dataset ("all") with species logos,
#           colored by species or location, and shapes for focal groups (e.g., "boc").
# Input   :
#   - pca_data     : Data frame of PCA coordinates (must contain PC1, PC2, sample, spec, geo)
#   - pc_first     : Name of first PC to plot ("PC1")
#   - pc_second    : Name of second PC to plot ("PC2")
#   - species_info : Data frame with species metadata (columns: spec, Color, link, Species)
#   - variance     : Data frame of variance explained per PC (column: variation)
#   - color_by     : Either "species" or "location"; controls coloring
# Output  :
#   - ggplot2 object of PCA scatter plot with logos and customized colors/shapes
# ============================================================
pca_plot_all <- function(pca_data, pc_first, pc_second, species_info, variance) {
    
    pca_data$boc <- ifelse(pca_data$geo == "boc", "Panama", "all other locations")

    # Always species mode
    group_col <- "spec"
    info_table <- species_info
    color_map <- setNames(info_table$Color, info_table$spec)
    label_map <- setNames(
        paste0("<img src='", info_table$link, "' width='120' /><br>*H. ", info_table$Species, "*"),
        info_table$spec
    )


    # Build PCA plot
    p <- ggplot(pca_data, aes(x = .data[[pc_first]], y = .data[[pc_second]])) +
        geom_point(aes(color = .data[[group_col]], shape = boc), size = 5, alpha = 0.8) +
        scale_color_manual(values = color_map, labels = label_map) +
        scale_shape_manual(values = c(19, 8)) + # filled circle for Panama, star/other for rest
        theme_minimal() +
        theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.text = element_markdown(size = 12),
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)
        ) +
        labs(
        x = paste0(pc_first, ", variance = ", format(round(variance$variation[1], 1), nsmall = 1), " %"),
        y = paste0(pc_second, ", variance = ", format(round(variance$variation[2], 1), nsmall = 1), " %")
        ) +
        guides(
        color = guide_legend(ncol = 1, byrow = TRUE),
        shape = guide_legend(ncol = 1)
        )

    return(p)

}


# ============================================================
# Function: plot_permanova_permdisp
# Purpose : Create a symmetric PERMANOVA & PERMDISP plot for all pairwise species comparisons,
#           excluding comparisons where either group has fewer than 4 samples.
# Input   :
#   - pair_file      : Path to the pairwise CSV containing columns spc1, spc2, n_spc1, n_spc2, permanova_corr_pval, permadisp_corr_pval.
#   - params_legend  : Position for legend ("right", "bottom", etc.).
# Output  :
#   - A ggplot object with pairwise categories colored and shaped by significance.
# ============================================================
plot_permanova_permdisp <- function(pair_file, params_legend = "none") {
    # Read CSV and select only needed columns
    pair_table <- read.table(file = pair_file, sep = ",", header = TRUE) %>%
        select(spc1, spc2, n_spc1, n_spc2, permanova_corr_pval, permadisp_corr_pval) %>%
        filter(n_spc1 > 4 & n_spc2 > 4) %>% # Exclude rows where either group has fewer than 5 individuals
        select(-n_spc1, -n_spc2) 
    
    if(nrow(pair_table) == 0) return(NULL)  # nothing to plot

    # Make symmetric
    df_sym <- rbind(
        pair_table,
        data.frame(
            spc1 = pair_table$spc2,
            spc2 = pair_table$spc1,
            permanova_corr_pval = pair_table$permanova_corr_pval,
            permadisp_corr_pval = pair_table$permadisp_corr_pval
        )
    ) %>%
        mutate(
            permanova_cat = case_when(
                permanova_corr_pval <= 0.01 ~ "<= 0.01",
                permanova_corr_pval <= 0.05 ~ "<= 0.05",
                permanova_corr_pval > 0.05  ~ "not significant"
            ),
            permadisp_cat = case_when(
                permadisp_corr_pval <= 0.01 ~ "<= 0.01",
                permadisp_corr_pval <= 0.05 ~ "<= 0.05",
                permadisp_corr_pval > 0.05  ~ "not significant"
            )
        )

    # Classify pairwise relationships
    df_sig <- df_sym %>%
        mutate(
            category = case_when(
                permanova_corr_pval <= 0.05 & permadisp_corr_pval <= 0.05 ~ "diff_withinvar",
                permanova_corr_pval <= 0.05 & permadisp_corr_pval > 0.05  ~ "struct",
                permanova_corr_pval > 0.05  & permadisp_corr_pval <= 0.05 ~ "sim_withinvar",
                permanova_corr_pval > 0.05  & permadisp_corr_pval > 0.05  ~ "no"
            )
        )

    # Reshape to matrix for plotting
    sig_pair <- dcast(as.data.table(df_sig), spc1 ~ spc2, value.var = "category")
    sig_pair[upper.tri(as.matrix(sig_pair[,-1]))] <- NA
    sig_pair_melt <- melt(sig_pair, id.vars = "spc1")

    # Build plot
    p_allLoc <- ggplot(sig_pair_melt) +
        geom_point(aes(x = spc1, y = variable, color = value, shape = value), size = 10) +
        scale_shape_manual(
        name = "permanova & permadisp\np-values",
        breaks = c("struct", "diff_withinvar", "no"),
        labels = c(
            expression(atop(""<="0.05" ~ "&" ~ "">="0.05", "structure evidence")),
            expression(atop(""<="0.05" ~ "&" ~ " "<="0.05", "within group variation")),
            expression(atop("">"0.05" ~ "&" ~ "">"0.05", "no structure"))
        ),
        values = c(19,21,18)
        ) +
        scale_colour_manual(
        name = "permanova & permadisp\np-values",
        breaks = c("struct", "diff_withinvar", "no"),
        labels = c(
            expression(atop(""<="0.05" ~ "&" ~ "">="0.05", "structure evidence")),
            expression(atop(""<="0.05" ~ "&" ~ " "<="0.05", "within group variation")),
            expression(atop("">"0.05" ~ "&" ~ "">"0.05", "no structure"))
        ),
        values = c("#FFB055", "#FE7309", "#D0210D")
        ) +
        scale_y_discrete(position = "right") +
        labs(x = "", y = "") +
        theme_minimal() +
        theme(
        legend.position = params_legend,
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.text = element_text(size = 12, margin = margin(0, 50, 0, 0)),
        legend.title = element_text(size = 15),
        legend.key.height = unit(1.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1
        )

    return(p_allLoc)
}

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
head(species_info)
geo_table <- read.delim(geo_colors, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
head(geo_table)

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
    all = list(dir = "byALL", colors = "species"),
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
    uni = list(dir = "bySPC", colors = "location"),
)

# Create a list to store plots per location
results <- list()

for(dat in names(dataset)) {
    dat_info <- dataset[[dat]]
    dat_dir <- dat_info$dir
    color <- dat_info$colors
  
    message("Processing: ", dat_info)
  
    #-----------------------------------
    # Read GTMAT file + Sample file + PERMANOVA & PERMDISP result table
    #-----------------------------------
    gtmat_file <- file.path(base_path, "2_popgen", dat_dir, if (dat_info == "all") "thinned_all_ld_pruned_gtmat.traw" else paste0(dat_info, "_ld_pruned_gtmat.traw"))
    sample_file <- file.path(base_path, if (dat_info == "all") "metadata/geno_names.txt" else paste0("2_popgen/", dat_dir, "/", dat_info, ".txt"))
    perm_file <- list.files(file.path(base_path, "2_popgen", dat_dir, "permanova_results"), pattern = paste0(dat_info, ".*\\.csv$"), full.names = TRUE)
    print(gtmat_file)
    print(sample_file)
    print(perm_file)

    #-----------------------------------
    # PCA
    #-----------------------------------
    pca_res <- pca_analysis(gtmat_file, sample_file, color_by = color)
    pca_eigen <- pca_res$eigen
    pca_var   <- pca_res$var

    #-----------------------------------
    # PCA plot
    #-----------------------------------
    if (dat == "all") {
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
    p_var <- plot_variance(var, dat)

    # -----------------------------------
    # PERMANOVA + PERMDISP (filter <5 inds per species inside perm_f)
    #-----------------------------------
    p_perm <- plot_permanova_permdisp(perm_file, params_legend = if(dat == "boc") c(0.3, 0.7) else if(dat == "all") c(0.2, 0.7) else "none")

    #-----------------------------------
    # Store outputs
    #-----------------------------------
    results[[dat]] <- list(
      pca_f = p_pca1,
      pca_s = p_pca2,
      pca_t = p_pca3,
      variance_plot = p_var,
      permanova = p_perm
    )
}


# ############################
# FINAL PLOTS
# ############################

# Set datasets for plots
results_locations <- results[names(results) %in% c("hon", "bel", "boc", "pri")]
keep_names <- names(results_locations)
print(keep_names)

########## FIGURE 3 ###################
# PCA plots for all locations with legend
all_pcas <- lapply(results_locations, `[[`, "pca_f") # extract per location pcas
pca_grid <- plot_grid(plotlist = all_pcas, ncol = 2) # bundle location pcas in one plot
# Combine PCA grid with legend at the bottom
figure3 <- ggarrange(pca_grid, labels=c('(a)','(b)','(c)','(d)'), common.legend=T, legend = "bottom") # adjust if legend is too big/small

# Save Figure 3 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "Fig3_gLocPCA.png"),
  plot = figure3,
  width = 18.5,    # A4 width in inches
  height = 25.5,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)


########## FIGURE 4 ###################
# Combined genetic space with legend
# Access each PCA plot for the "all" dataset
pca1 <- results[["all"]][["pca_f"]]
pca2 <- results[["all"]][["pca_s"]]
pca3 <- results[["all"]][["pca_t"]]

figure4 <- ggarrange(pca1, pca2, pca3, ncol = 1, nrow = 3, labels=c('(a)','(b)','(c)'), common.legend=T, legend = "right")

ggsave(
  filename = file.path(figure_path, "Fig4_gAllPCA.png"),
  plot = figure4,
  width = 13, 
  height = 26, 
  units = "in",      # inches
  dpi = 150,         # moderate dpi to reduce file size but keep quality
  type = "cairo-png" # better compression and anti-aliasing
)

########## FIGURE S11 ###################
# Variance of Principal Components for combined genetic space
figureS11 <- results[["all"]][["variance_plot"]]
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
sup_grid <- plot_grid(plotlist = all_sup, ncol = 2) # bundle location pcas in one plot
# Combine PCA grid with legend at the bottom
figureS12 <- ggarrange(sup_grid, labels=c('(a)','(b)','(c)','(d)'), common.legend=T, legend = "bottom") # adjust if legend is too big/small

# Save Figure S12 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS12_gLocSUP.png"),
  plot = figureS12,
  width = 18.5,    # A4 width in inches
  height = 25.5,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S13 ###################
# PERMANOVA heatmaps for each location
all_perm <- lapply(results_no_overall, `[[`, "permanova") # extract per location pcas
figureS13 <- plot_grid(plotlist = all_perm, ncol = 2) # bundle location pcas in one plot

# Save Figure S13 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS13_gLocPERM.png"),
       plot = figureS13,
       width = 14.2,    # A4 width in inches
       height = 17,  # A4 height in inches
       units = "in",
       dpi = 150,       # good quality but light (~1 MB)
       type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S14 ###################
# Combined genotypic space: PERMANOVA
figureS14 <- results[["all"]][["permanova"]]

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS14_gAllPERM.png"),
  plot = figureS14,
  width = 14.2,    # A4 width in inches
  height = 17,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)

########## FIGURE S15 ###################
# Per species genotypic space: PCA + PERMANOVA + PERMDISP
# Remove the overall entry before extracting plots
results_spc <- results[names(results) %in% c("pue", "nig", "uni", "pri")]
keep_spc <- names(results_spc)
print(keep_spc)

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

figureS15 <- plot_grid(pue, nig, uni, nrow = 3, rel_heights = c(1, 1, 1))

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS15_gSpe.png"),
  plot = figureS15,
  width = 18.5,    # A4 width in inches
  height = 25.5,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)