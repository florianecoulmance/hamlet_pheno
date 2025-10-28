# by: Floriane Coulmance: 16/08/2024
# usage:
# Rscript Fig1-2_S3-S9.R 
#___________________________________________________________________

# Clear the work space
rm(list = ls())

# Load needed library
# library(magick)
# library(vegan)
# library(factoextra)
# library(reshape2)
# library(pairwiseAdonis)
# library(ggrepel)
# library(stringr)
# library(tidyverse)
# library(ggimage)
# library(ggtext)
# library(ggplot2)
# library(scales)
# library(ggnewscale)
# library(ggpubr)
# library(cluster)    # clustering algorithms
# library(dendextend)
# library(hypoimg)
# library(ggtree)
# library(ggtreeExtra)
# library(grid)
# library(png)
# library(ape)
# library(dplyr)



# new libraries
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

# ############################
# FUNCTIONS
# ############################

# ============================================================
# Helper function to extract named arguments
# ============================================================
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && length(args) > idx) {
    return(args[idx + 1])
  } else {
    return(default)
  }
}


# ============================================================
# Function to add logo links to species color table
# ============================================================
add_species_logos <- function(color_file, logos_p) {
  # Read species-color table
  species_info <- read.delim(color_file, stringsAsFactors = FALSE)
  
  # Ensure 'Species' is lowercase and trimmed
  species_info$Species <- tolower(trimws(species_info$Species))
  
  # Add column with full path to logos based on species abbreviation
  species_info$link <- file.path(logos_path, paste0("H_", species_info$Species, ".l.cairo.png"))
  
  # Return species color table with links to logo
  return(species_info)
}


# ============================================================
# Function to get legend ffrom a plot
# ============================================================
# get_legend <- function(my_plot) {
#   tmp <- ggplotGrob(my_plot + theme(legend.position = "bottom")) # make sure it has a legend
#   leg <- gtable::gtable_filter(tmp, "guide-box")
#   if (length(leg) == 0) return(NULL)
#   leg
# }


# ============================================================
# Function: write_metadata_gxp
# Purpose : Prepare metadata table for downstream analyses
# Input   :
#   - PCs:   Data frame containing PCA results (with "images" column)
#   - path_meta: Path to metadata (not used here since no writing)
#   - effect, m_type, dat: extra arguments kept for compatibility
# Output  :
#   - Returns a cleaned metadata table with columns:
#       sample, geo, spec, pop, im
# ============================================================
write_metadata_gxp <- function(PCs) {
  
  # Extract sample name from image filename (remove everything after first "-")
  PCs$sample <- gsub("\\-.*", "", PCs$images)
  
  # Derive metadata columns from sample names
  PCs$geo  <- stringr::str_sub(PCs$sample, -3, -1)   # last 3 chars = location
  PCs$spec <- stringr::str_sub(PCs$sample, -6, -4)   # chars -6 to -4 = species code
  PCs$pop  <- stringr::str_sub(PCs$sample, -6, -1)   # last 6 chars = population
  PCs$im   <- gsub('.{4}$', '', PCs$images)          # remove last 4 chars from image name
  
  # Drop unused column if present
  PCs$X <- NULL
  
  # Return metadata table
  return(PCs)
  
}


# ============================================================
# Function: pca_plot
# Purpose : Create PCA scatter plot with species or location colors,
#           including centroids, ellipses, and optional logos.
# Input   :
#   - pca_data       : Data frame of PCA results; must contain "spec" and "geo" columns.
#   - pc_first       : Name of the first principal component to plot (e.g., "PC1").
#   - pc_second      : Name of the second principal component to plot (e.g., "PC2").
#   - species_info   : Data frame with species metadata (columns: spec, Color, link, Species).
#   - geo_info       : Data frame with location metadata (columns: geo, Color, Locations).
#   - var            : Data frame containing explained variance per PC (column X0).
#   - color_by       : Either "species" (default) or "location"; controls grouping and coloring.
#   - extract_legend : If TRUE, return only the legend (default: FALSE).
#   - legend_rows    : Number of rows for the legend (default: 2).
# Output  :
#   - If extract_legend = FALSE: a ggplot2 PCA scatter plot (annotated with title).
#   - If extract_legend = TRUE : a legend grob for combined plotting.
# ============================================================
pca_plot <- function(pca_data, pc_first, pc_second, species_info, geo_info, var, color_by = "species", extract_legend = FALSE, legend_rows = 2) {

  # -----------------------------
  # 1. Choose color/grouping mode
  # -----------------------------
  color_by <- match.arg(color_by, c("species", "location"))
  
  if (color_by == "species") {
    group_col <- "spec"
    info_table <- species_info
    color_map <- setNames(info_table$Color, info_table$spec)
    label_map <- setNames(
      paste0("<img src='", info_table$link, "' width='100' /><br>*", info_table$Species, "*"),
      info_table$spec)
  } else {
    group_col <- "geo"
    info_table <- geo_info
    color_map <- setNames(info_table$Color, info_table$geo)
    label_map <- setNames(info_table$Locations, info_table$geo)
  }
  
  # print(head(pca_data))
  # -----------------------------
  # 2. Merge PCA data with info
  # -----------------------------
  plot_data <- pca_data %>%
    left_join(info_table, by = group_col)
  
  # -----------------------------
  # 3. Calculate centroids
  # -----------------------------
  centroids <- plot_data %>%
    group_by(.data[[group_col]]) %>%
    summarise(x = mean(.data[[pc_first]], na.rm = TRUE),
              y = mean(.data[[pc_second]], na.rm = TRUE),
              Color = first(Color),
              link = if (color_by == "species") first(link) else NA_character_,
              .groups = "drop")
  print(centroids$link)
  
  # PCA scatter plot with centroids and ellipses
  p <- ggplot(plot_data, aes(x = .data[[pc_first]], y = .data[[pc_second]], color = .data[[group_col]])) +
    geom_point(size = 5, alpha = 0.5) +
    stat_ellipse(aes(color = .data[[group_col]]), linetype = 5, lwd = 1) +
    # geom_point(data = centroids, aes(x = x, y = y, color = .data[[group_col]]), size = 15, alpha = 1) +
    # geom_image(data = centroids, aes(x = x, y = y, image = link), vjust=1, hjust=0, size = 0.15, asp = 1.1, alpha=1) +
    scale_color_manual(values = color_map, labels = label_map) +
    theme_minimal() +
    theme(
      legend.position = ifelse(extract_legend, "bottom", "none"),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10)
    ) +
    scale_x_continuous(position = "bottom",labels = unit_format(unit = "k", scale = 1e-3)) +
    scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3)) +
    labs(
      x = paste0(pc_first,", variance =  ", format(round(var$X0[as.numeric(str_sub(pc_first, 3, -1))] * 100, 1), nsmall = 1), " %"),
      y = paste0(pc_second,", variance = ", format(round(var$X0[as.numeric(str_sub(pc_second, 3, -1))] * 100, 1), nsmall = 1), " %")
    )
  
  # -----------------------------
  # 5. Add species logos (if applicable)
  # -----------------------------
  print(names(centroids))
  if (color_by == "species" && "link" %in% names(centroids)) {
    p <- p +
      geom_image(
        data = centroids,
        aes(x = x, y = y, image = link),
        inherit.aes = FALSE,
        vjust = 1, hjust = 0, size = 0.1, asp = 1.1, alpha = 1
      )
  }

  # -----------------------------
  # 6. Add title (per location or per species)
  # -----------------------------
  # ---- Get title depending on color_by mode ----
  title_val <- NULL
  if (color_by == "species") {
    # if coloring by species, title is the location name
    geo_val <- unique(plot_data$geo)
    if (length(geo_val) == 1) {
      title_val <- geo_info$Locations[geo_info$geo == geo_val]
    }
  } else if (color_by == "location") {
    # if coloring by location, title is the species name
    spec_val <- unique(plot_data$spec)
    if (length(spec_val) == 1) {
      title_val <- species_info$Species[species_info$spec == spec_val]
    }
  }

  # -----------------------------
  # 7. Return plot or legend
  # -----------------------------
  if (extract_legend) {
    p_legend <- p +
      theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = legend_rows))
    legend <- get_legend(p_legend)
    return(legend)
  } else {
    # ---- Annotate with location title ----
    p_annot <- annotate_figure(
      p,
      top = text_grob(title_val, color = "black", face = "bold", size = 15,
                      x = unit(5.5, "pt"), hjust = 0)
    )
    
    return(p_annot)
  }

}


# ============================================================
# Function: plot_variance
# Purpose : Plot the proportion of variance explained by each PC
# Input   :
#   - var       : data frame containing PC numbers and variance values
#   - loc_name  : character string used for the plot title (location name)
# Output  :
#   - ggplot object showing explained variance per PC
# ============================================================
plot_variance <- function(var, loc_name) {
  # Ensure proper column names and indexing
  if (ncol(var) == 2) {
    colnames(var) <- c("PC", "Variance")
  } else if (ncol(var) == 1) {
    var$PC <- seq_len(nrow(var))
    colnames(var) <- c("Variance", "PC")
  }
  
  # Create explained variance plot
  p <- ggplot(var, aes(x = as.factor(PC + 1), y = Variance)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = sprintf("%.1f%%", Variance * 100)),
              vjust = -0.5, size = 3) +
    labs(x = "Principal Component", 
         y = "Explained variance (%)",
         title = paste0("Explained variance - ", loc_name)) +
    theme_minimal(base_size = 12)
  
  return(p)
}


# ============================================================
# Function: perm_f
# Purpose : Run pairwise PERMANOVA on first 15 PCs and visualize
#           R² and F.Model values as a combined heatmap.
# Input   :
#   - pc_table     : PCA table (must include 'spec', 'geo', and 'sample' columns)
#   - species_info : Data frame with species metadata (spec, Species)
#   - geo_map      : Data frame with location metadata (geo, Locations)
#   - color_by     : Either "species" or "location"; controls grouping.
# Output  :
#   - Annotated ggplot2 heatmap of PERMANOVA R² (upper) and F.Model (lower)
# ============================================================
perm_f <- function(pc_table, species_col, geo_map, color_by = "species") {
  
  # ---- Filter groups with <5 individuals ----
  if(color_by == "species") {
    keep <- names(table(pc_table$spec)[table(pc_table$spec) >= 5])
    pc_table <- pc_table %>% filter(spec %in% keep)
    pc_table$group <- factor(pc_table$spec)
  } else if(color_by == "location") {
    keep <- names(table(pc_table$geo)[table(pc_table$geo) >= 5])
    pc_table <- pc_table %>% filter(geo %in% keep)
    pc_table$group <- factor(pc_table$geo)
  }
  
  # ---- Prepare data ----
  PC_table_ordered <- pc_table[order(pc_table$group), ]
  rownames(PC_table_ordered) <- PC_table_ordered$sample
  
  # ---- Compute Euclidean distance ----
  dist_t <- dist(PC_table_ordered[,1:15], method = "euclidean")
  
  # ---- Pairwise PERMANOVA ----
  pairwise_result <- pairwise.adonis(dist_t, PC_table_ordered$group, perm = 10000)
  
  # ---- Reshape R² (upper triangle) ----
  vis_R <- pairwise_result %>%
    separate(pairs, into = c("group1","group2"), sep = " vs ") %>%
    select(group1, group2, R2)
  
  visH_R <- dcast(vis_R, group1 ~ group2, value.var = "R2") %>%
    complete(group1 = levels(pc_table$group), fill = list()) %>%
    select(group1, intersect(levels(pc_table$group), colnames(.)))
  melt_R <- melt(visH_R, id.vars = "group1")
  melt_R$value <- as.numeric(melt_R$value)
  
  # ---- Reshape F.Model (lower triangle) ----
  vis_F <- pairwise_result %>%
    separate(pairs, into = c("group1","group2"), sep = " vs ") %>%
    select(group1, group2, F.Model)
  
  visH_F <- dcast(vis_F, group1 ~ group2, value.var = "F.Model") %>%
    complete(group1 = levels(pc_table$group), fill = list()) %>%
    select(group1, intersect(levels(pc_table$group), colnames(.))) %>%
    arrange(desc(group1))
  melt_F <- melt(visH_F, id.vars = "group1")
  melt_F$value <- as.numeric(melt_F$value)
  
  # ---- Plot heatmap ----
  p <- ggplot() +
    geom_tile(aes(x = melt_R$group1, y = melt_R$variable, fill = melt_R$value), color = "transparent") +
    scale_fill_gradient(low = "#ffead1", high = "#fdae53", na.value = "transparent", name = "R²") +
    new_scale_fill() +
    geom_tile(aes(x = melt_F$variable, y = melt_F$group1, fill = melt_F$value), color = "transparent") +
    scale_fill_gradient(low = "#ffedec", high = "#ff5c52", na.value = "transparent", name = "F") +
    labs(x = "", y = "", fill = "F") +
    scale_x_discrete(position = "bottom") +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      legend.direction = "vertical",
      legend.box = "horizontal",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.key.height = unit(1, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, size = 12),
      axis.text.y = element_text(size = 12),
      axis.ticks = element_blank(),
      aspect.ratio = 1
    )
  
  # ---- Get title ----
  if(color_by == "species"){
    geo_val <- unique(pc_table$geo)
    title_val <- if(length(geo_val) == 1) geo_map$Locations[geo_map$geo == geo_val] else ""
  } else if(color_by == "location"){
    species_val <- unique(pc_table$spec)
    if (length(species_val) == 1 && species_val %in% species_col$spec){
      title_val <- species_col$Species[species_col$spec == species_val] 
      }
  }
  
  # ---- Annotate with location title ----
  p_annot <- annotate_figure(
    p,
    top = text_grob(title_val, color = "black", face = "bold", size = 30,
                    x = unit(5.5, "pt"), hjust = -0.4)
  )
  
  return(p_annot)
  
}


# ============================================================
# Function: hierClustering
# Purpose : Perform hierarchical clustering on the first 15 scaled PCs,
#           coloring tips by either species or location.
# Input   :
#   - data_path   : Directory containing PCA files
#   - pca_file    : PCA CSV file (must include an "images" column)
#   - species_col : TSV with species color info (spec, Color, Species)
#   - geo_map     : TSV with location color info (geo, Color, Locations)
#   - color_by    : "species" (default) or "location"
#   - extract_legend : If TRUE, returns legend instead of tree
#   - legend_rows : Number of rows in the legend (default: 2)
# Output  :
#   - Annotated ggtree plot (or legend if extract_legend = TRUE)
# Notes   :
#   - Uses the first 15 PCs, scaled
#   - Distance: Euclidean; Clustering: Ward.D
# ============================================================
hierClustering <- function(data_path, pca_file, species_col, geo_map, color_by = "species", extract_legend = FALSE, legend_rows = 2) 
{

  # -----------------------------
  # 1. Load PCA data
  # -----------------------------
  pca_t <- read.csv(file.path(data_path, pca_file), header = TRUE) %>%
    mutate(
      # Fix specific naming error
      images = ifelse(images == "PL17_160pueflo-l1-s4-f4-c2-d1.png",
                      "PL17_160floflo-l1-s4-f4-c2-d1.png", images),
      label = gsub("-.*", "", images)
    ) %>%
    select(label, everything(), -X, -images) %>%
    column_to_rownames(var = "label") %>%
    select(1:15)   # always 15 PCs
  
  # Scale data
  pca_scaled <- scale(pca_t)
  
  # -----------------------------
  # 2. Hierarchical clustering
  # -----------------------------
  d <- dist(pca_scaled, method = "euclidean")
  hc <- hclust(d, method = "ward.D")
  
  # Extract tree data
  tree <- ggtree(hc)$data %>%
    mutate(
      spec = if_else(isTip, str_sub(label, -6, -4), "ungrouped"),
      geo  = if_else(isTip, str_sub(label, -3, -1), "ungrouped")
    )
  
  # Adjust positions for aesthetics
  tree$x[tree$isTip] <- tree$x[tree$isTip] * 1.1
  tree$branch.length <- scale(tree$branch.length)
  
  # -----------------------------
  # 3. Select color mode
  # -----------------------------
  color_by <- match.arg(color_by, c("species", "location"))
  
  if (color_by == "species") {
    color_table <- species_col
    color_table$spec <- tolower(trimws(color_table$spec))
    color_map <- setNames(color_table$Color, color_table$spec)
    group_col <- "spec"
    legend_name <- "Species"
  } else {
    color_table <- geo_map
    color_table$geo <- tolower(trimws(color_table$geo))
    color_map <- setNames(color_table$Color, color_table$geo)
    group_col <- "geo"
    legend_name <- "Location"
  }

  # -----------------------------
  # 4. Join color info
  # -----------------------------
  scol <- tree %>%
    left_join(color_table, by = setNames(group_col, group_col)) %>%
    replace_na(list(Color = "gray20")) %>%
    pull(Color)
  
  # -----------------------------
  # 5. Plot tree
  # -----------------------------
  t <- ggtree(tree, layout = "fan", size = 0.5, color = scol) +
    geom_tippoint(aes(color = .data[[group_col]]), size = 3, alpha = 0.5) +
    scale_color_manual(values = color_map, name = legend_name) +
    theme(
      legend.position = ifelse(extract_legend, "right", "none"),
      plot.title = element_text(size = 14, color = "gray20", face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray20")
    )
  
  # -----------------------------
  # 6. Get title
  # -----------------------------
  title_val <- ""

  if (color_by == "species") {
    # If clustering by species (i.e. showing all species from one location)
    geo_val <- unique(tree$loc)
    if (length(geo_val) == 1 && geo_val %in% geo_map$geo) {
      title_val <- geo_map$Locations[geo_map$geo == geo_val]
    }
    
  } else if (color_by == "location") {
    # If clustering by location (i.e. showing all locations for one species)
    spec_val <- unique(tree$spec)
    if (length(spec_val) == 1 && spec_val %in% species_col$spec) {
      title_val <- species_col$Species[species_col$spec == spec_val]
    }
  }
  
  # -----------------------------
  # 7. Return legend or annotated tree
  # -----------------------------
  if (extract_legend) {
    legend <- cowplot::get_legend(
      t + guides(color = guide_legend(ncol = legend_rows))
    )
    return(legend)
  } else {      # ---- Annotate with location title ----
    t_annot <- annotate_figure(
      t,
      top = text_grob(title_val, color = "black", face = "bold", size = 30,
                      x = unit(5.5, "pt"), hjust = -0.4)
    )
    return(t_annot)
  }
}


# ============================================================
# Function: heat_plots
# Purpose : Generate annotated heatmaps of selected PCs
#           (e.g. PC1, PC3) with legend and location/species title.
# Input   :
#   - im_p      : Directory containing PC image files.
#   - name      : File prefix used in image filenames 
#                 (must include a 3-letter geo or species code 
#                  before the first number, e.g. "lab_boc229...").
#   - pcs       : Numeric vector of PC numbers to plot (e.g. c(1, 3)).
#   - spec_map  : Data frame with columns: Species, spec, Color.
#   - geo_map   : Data frame with columns: Locations, geo, Color.
#   - color_by  : Either "species" (default; between species within a location)
#                 or "location" (within species across locations).
# Output  :
#   - A ggarrange (cowplot) object with PC heatmaps, shared legend, and
#     a title corresponding to the selected mode (location or species).
# Notes   :
#   - Assumes PC images exist as PNG files named:
#       <name>_PC<pc>_originalrescaled.png
#   - Cropping indices are hard-coded for 1000×684 image size.
#   - Title logic matches PCA and hierarchical clustering functions.
# ============================================================
heat_plots <- function(im_p, name, pcs, spec_map, geo_map, color_by = "species") {
  
  # -----------------------------
  # Helper: Load and crop PC image
  # -----------------------------
  load_pc_img <- function(pc) {
    file <- file.path(im_p, paste0(name, "_", pc, "_originalrescaled.png"))
    if (!file.exists(file)) stop(paste("Missing image file:", file))  
    img <- readPNG(file)
    grob <- rasterGrob(img[1:500,100:1000,], interpolate = TRUE)
    grob <- annotate_figure(grob, top = text_grob(
      paste0("PC", pc),
      color = "black", face = "bold", size = 20,
      x = unit(5.5, "pt"), hjust = -2.5, vjust = 4
    ))
    return(grob)
  }
  
  # main PCs
  plots <- lapply(pcs, load_pc_img)
  
  # legend (from bottom of last PC img)
  file_last <- file.path(im_p, paste0(name, "_", tail(pcs,1), "_originalrescaled.png"))
  if (!file.exists(file_last)) stop(paste("Missing legend source image:", file_last))
  
  img_last <- readPNG(file_last)
  legend <- rasterGrob(img_last[600:684,,], interpolate = TRUE)
  
  # arrange
  combined <- ggarrange(plotlist = c(plots, list(legend)), 
                        nrow = length(pcs)+1, 
                        heights = c(rep(1, length(pcs)), 0.3),
                        align = "v")
  
  # extract geo code from filename (3 letters before first number)
  abbrev <- sub(".*?([a-z]{3})(?=[0-9]).*", "\\1", name, perl = TRUE)
  

  # -----------------------------
  # 6. Get title
  # -----------------------------
  if (color_by == "species") {
    # Between species within a location → use location name
    if (abbrev %in% geo_map$geo) {
      title_val <- geo_map$Locations[match(abbrev, geo_map$geo)]
    } else {
      title_val <- ""
    }
    
  } else if (color_by == "location") {
    # Within species across locations → use full species name
    if (abbrev %in% spec_map$spec) {
      title_val <- spec_map$Species[match(abbrev, spec_map$spec)]
    } else {
      title_val <- ""
    }
    
  } else {
    stop("Invalid value for color_by. Use 'species' or 'location'.")
  }
  
  # add title
  combined <- annotate_figure(combined, top = text_grob(
    title_val, color = "black", face = "bold", size = 30,
    x = unit(5.5, "pt"), hjust = -0.4
  ))
  
  return(combined)

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

# base_path <- "/Users/fcoulman/Desktop/hamlet_pheno/3_CHAPTER3/hamlet_pheno/"
# path_phenotypes <- file.path(base_path, "alignment_folder/571_left_noflash/pca/")
# heatmap_path <- file.path(base_path, "alignment_folder/571_left_noflash/images/LAB_greenBG/")
# figure_path <- file.path(base_path, "figures/")
# logos_path <- file.path(base_path, "metadata/logos_hamlet/")
# spec_colors <- file.path(base_path, "metadata/species_colors.tsv")
# geo_colors <- file.path(base_path, "metadata/locations_colors.tsv")



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

# Combine all data for legend
df_all <- do.call(rbind, lapply(results[!names(results) %in% c("all", "pue", "nig", "uni", "chl", "abe", "ind")], function(x) x$data))
var_all <- do.call(rbind, lapply(results[!names(results) %in% c("all", "pue", "nig", "uni", "chl", "abe", "ind")], function(x) x$variance))

# Make a dummy plot to extract legend
# p_dummy <- pca_plot(df_all, "PC1", "PC2", species_info, geo_table, var_all, color_by = "species", extract_legend = TRUE)
# print(class(p_dummy))
# Extract combined legend
combined_legend <- pca_plot(df_all, "PC1", "PC2", species_info, geo_table, var_all, color_by = "species", extract_legend = TRUE)
print(class(combined_legend))

########## FIGURE 1 ###################
# PCA plots for all locations with legend
all_pcas <- lapply(results_no_overall, `[[`, "pca") # extract per location pcas
pca_grid <- plot_grid(plotlist = all_pcas, ncol = 2) # bundle location pcas in one plot
# Combine PCA grid with legend at the bottom
figure1 <- plot_grid(pca_grid, combined_legend, ncol = 1, rel_heights = c(1, 0.15)) # adjust if legend is too big/small

# Save Figure 1 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "Fig1_pLocPCA.png"),
  plot = figure1,
  width = 8.27,    # A4 width in inches
  height = 11.69,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)


########## FIGURE 2 ###################
# Combined phenotypic space with legend
figure2 <- results[["all"]][["pca"]]
print(dim(results[["all"]][["pca"]]))
print(head(results[["all"]][["pca"]]))
ggsave(
  filename = file.path(figure_path, "Fig2_pAllPCA.png"),
  plot = figure2,
  width = 8.27, 
  height = 5.22, 
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
sup_grid <- plot_grid(plotlist = all_sup, ncol = 2) # bundle location pcas in one plot
# Combine supplementary PCA grid with legend at the bottom
figureS4 <- plot_grid(sup_grid, combined_legend, ncol = 1, rel_heights = c(1, 0.15)) # adjust if legend is too big/small

# Save Figure S4 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS4_pLocSUP.png"),
  plot = figureS4,
  width = 8.27,    # A4 width in inches
  height = 11.69,  # A4 height in inches
  units = "in",
  dpi = 150,       # good quality but light (~1 MB)
  type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S5 ###################
# PERMANOVA heatmaps for each location
all_perm <- lapply(results_no_overall, `[[`, "permanova") # extract per location pcas
figureS5 <- plot_grid(plotlist = all_perm, ncol = 2) # bundle location pcas in one plot

# Save Figure S3 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS5_pLocPERM.png"),
       plot = figureS5,
       width = 8.27,    # A4 width in inches
       height = 11.69,  # A4 height in inches
       units = "in",
       dpi = 150,       # good quality but light (~1 MB)
       type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S6 ###################
# Hierarchical clustering plots for all locations with legend
all_hier <- lapply(results_no_overall, `[[`, "hclust") # extract per location pcas
hier_grid <- plot_grid(plotlist = all_hier, ncol = 2) # bundle location pcas in one plot
# Combine PCA grid with legend at the bottom
figureS6 <- plot_grid(hier_grid, combined_legend, ncol = 1, rel_heights = c(1, 0.15)) # adjust if legend is too big/small

# Save Figure S4 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS6_pLocHCLUST.png"),
       plot = figureS6,
       width = 8.27,    # A4 width in inches
       height = 11.69,  # A4 height in inches
       units = "in",
       dpi = 150,       # good quality but light (~1 MB)
       type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S7 ###################
# Heatmap PC images for each location
all_heat <- lapply(results_no_overall, `[[`, "heatmap") # extract per location pcas
figureS7 <- plot_grid(plotlist = all_heat, ncol = 2) # bundle location pcas in one plot

# Save Figure S5 as A4 PNG, optimized for small file size
ggsave(filename = file.path(figure_path, "FigS7_pLocHEAT.png"),
       plot = figureS7,
       width = 8.27,    # A4 width in inches
       height = 11.69,  # A4 height in inches
       units = "in",
       dpi = 150,       # good quality but light (~1 MB)
       type = "cairo-png" # smoother text rendering, smaller file
)

########## FIGURE S8 ###################
# Combined phenotypic space: supplementary PCA + PERMANOVA + hierarchical clustering + heatmaps
sup <- results[["all"]][["sup_pca"]]
perm <- results[["all"]][["permanova"]]
hier <- results[["all"]][["hclust"]]
heat <- results[["all"]][["heatmap"]]

# # Bottom row: hier + heat
# bottom_row <- plot_grid(hier, heat, ncol = 2, rel_widths = c(1,1))

# Combine top (perm) with bottom row
figureS8 <- plot_grid(sup, perm, hier, heat, nrow = 2, ncol = 2, rel_widths = c(1,1), rel_heights = c(1, 1))

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS8_pAll.png"),
  plot = figureS8,
  width = 8.27,    # A4 width in inches
  height = 11.69,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)

########## FIGURE S9 ###################
# Per species phenotypic space: PCA + heatmaps + hierarchical clustering + PERMANOVA
pca_pue <- results[["pue"]][["pca"]]
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

figureS9 <- plot_grid(pue, nig, uni, chl, abe, ind, nrow = 6, rel_heights = c(1, 1, 1, 1, 1, 1))

# Save as PNG (A4 size)
ggsave(
  filename = file.path(figure_path, "FigS9_pSpe.png"),
  plot = figureS9,
  width = 8.27,    # A4 width in inches
  height = 11.69,  # A4 height in inches
  units = "in",
  dpi = 150,
  type = "cairo-png"
)