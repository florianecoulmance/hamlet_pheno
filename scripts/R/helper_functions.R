# by: Floriane Coulmance: 05/11/2025
# usage:
# source("../scripts/R/helper_functions.R") 
# This file contains functions only.
#___________________________________________________________________

# Clear the work space
rm(list = ls())


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
# Function to extract legend from a plot
# ============================================================
exc_legend <- function(pca_file, species_col, num_row = 2) {
  # create a list of species logo to integrate to plots
  new_spec <- species_col %>%
    mutate(
      spec = spec,
      html = glue(
        "<img src='{link}' width='110' /><br>*H. {Species}*"
      )
    )
  
  logos_spec <- setNames(new_spec$html, new_spec$spec)
  # Inspect the result
  print(logos_spec)

  spec_colors <- setNames(species_col$Color, species_col$spec)
  print(spec_colors)

  # PC1 vs. PC2
  p1 <- ggplot(pca_file,aes(x=PC1,y=PC2,color=spec)) + geom_point(size = 15, alpha = 0.5) 
  p1 <- p1 + scale_color_manual(values=spec_colors, labels = logos_spec) +
    scale_x_continuous(position = "bottom") +
    theme(legend.position="bottom",legend.title=element_blank(),legend.box = "horizontal",legend.text =  element_markdown(size = 15),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.9),
          text = element_text(size=30), legend.key=element_blank()) +
    guides(color = guide_legend(nrow = num_row))
  
  
  legend <- get_legend(p1)

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
  print(extract_legend)
  # -----------------------------
  # 1. Choose color/grouping mode
  # -----------------------------
  color_by <- match.arg(color_by, c("species", "location"))
  
  if (color_by == "species") {
    group_col <- "spec"
    info_table <- species_info
    color_map <- setNames(info_table$Color, info_table$spec)
    label_map <- setNames(
      paste0("<img src='", info_table$link, "' width='120' /><br>*H. ", info_table$Species, "*"),
      info_table$spec)
  } else {
    group_col <- "geo"
    info_table <- geo_info
    color_map <- setNames(info_table$Color, info_table$geo)
    label_map <- setNames(info_table$Locations, info_table$geo)
  }
  print(color_map)
  print(label_map)

  
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

    # -----------------------------
  # 5. Add species logos if applicable
  # -----------------------------
  if (color_by == "species" && "link" %in% names(centroids)) {
    add_logos <- geom_image(
        data = centroids,
        aes(x = x, y = y, image = link),
        inherit.aes=FALSE,
        vjust = 1, hjust = 0, size = 0.12, asp = 1.1, alpha = 1
      )
  } else {
    add_logos <- NULL
  }
  print(add_logos)

  # -----------------------------
  # 4. Build the actual PCA plot
  # -----------------------------
  p <- ggplot(plot_data, aes(x = .data[[pc_first]], y = .data[[pc_second]], color = .data[[group_col]])) +
    geom_point(size = 5, alpha = 0.5) +
    stat_ellipse(aes(color = .data[[group_col]]), linetype = 5, lwd = 1) +
    geom_point(data = centroids, aes(x = x, y = y, color = .data[[group_col]]), size = 15, alpha = 0) +
    # geom_image(data = centroids, aes(x = x, y = y, image = link), vjust=1, hjust=0, size = 0.15, asp = 1.1, alpha=1) +
    scale_color_manual(values = color_map, labels = label_map) +
    theme_minimal() +
    theme(
      legend.position = "bottom", #if (extract_legend) "bottom" else "none",
      legend.box = "horizontal",
      legend.text = element_markdown(size = 10),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 20)
    ) +
    scale_x_continuous(position = "bottom",labels = unit_format(unit = "k", scale = 1e-3)) +
    scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3)) +
    labs(
      x = paste0(pc_first,", variance =  ", format(round(var$X0[as.numeric(str_sub(pc_first, 3, -1))] * 100, 1), nsmall = 1), " %"),
      y = paste0(pc_second,", variance = ", format(round(var$X0[as.numeric(str_sub(pc_second, 3, -1))] * 100, 1), nsmall = 1), " %")
    ) +
    guides(color = guide_legend(nrow = legend_rows))

  # print(class(p))
  # if(extract_legend){
  #   leg <- cowplot::get_legend(p)
  #   print(class(leg))
  # }

  
  if(!is.null(add_logos)) p <- p + add_logos
  print(class(p))


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
  # 7. Return plot
  # -----------------------------
  p_annot <- annotate_figure(
    p,
    top = text_grob(title_val, color = "black", face = "bold", size = 15,
                      x = unit(5.5, "pt"), hjust = 0)
  )
    
  return(p_annot)

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
  print(var)
  # Ensure var is a data frame
  var <- as.data.frame(var)
  print(var)
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
hierClustering <- function(data_path, pca_file, species_col, geo_map, color_by = "species", extract_legend = FALSE, legend_rows = 2) {

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
      pc,
      color = "black", face = "bold", size = 15,
      x = unit(2, "pt"), hjust = -2.5, vjust = 4
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
    print(allsamples)

    geo  <- substr(allsamples, nchar(allsamples) - 2, nchar(allsamples))
    spec <- substr(allsamples, nchar(allsamples) - 5, nchar(allsamples) - 3)

    # Extract spec and geo depending on mode
    if (speciesMODE) {
        if(length(unique(spec)) <= 1) {
            message("Not enough different samples")
            quit(save = "no", status = 1)
        }
    } else {
        if(length(unique(geo)) <= 1) {
            message("Not enough different samples")
            quit(save = "no", status = 1)
        }
    }
    
    # Run SmartPCA
    sm.pca <- smart_pca(
        snp_data = gtfile,
        sample_group = if (speciesMODE) spec else geo,
        missing_value = NA,
        pc_axes = 6
    )
    
    # Extract PCA outputs
    pca_coords <- sm.pca$pca.sample_coordinates
    var_explained <- sm.pca$pca.eigenvalues
    print(head(pca_coords))
    print(head(var_explained))

    # Ensure the PCA coordinates rows are in the same order as the original samples
    #pca_coords <- pca_coords[match(allsamples, rownames(pca_coords)), ]
    #print(head(pca_coords))

    # Then add spec, geo, and sample
    pca_coords$spec   <- spec
    pca_coords$geo    <- geo
    pca_coords$sample <- allsamples
    print(head(pca_coords))    

    # Extract the "variance explained" row as numeric
    var_explained <- as.numeric(var_explained["variance explained", ])
    print(var_explained)

    # If currently in percentages, convert to fractions
    var_explained <- var_explained / 100  # remove if already in 0-1 range

    # Create the table in the same format as read.csv would
    table <- data.frame(X0 = var_explained)
    print(table)
    # Add row numbers starting from 0 if needed
    rownames(table) <- 0:(length(var_explained)-1)
    print(table)

    # Access the values like
    (print(table$X0))

    # Return both
    return(list(eigen = pca_coords, var = table))
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
        x = paste0(pc_first, ", variance = ", format(round(variance$X0[as.numeric(str_sub(pc_first, 3, -1))] * 100, 1), nsmall = 1), " %"),
        y = paste0(pc_second, ", variance = ", format(round(variance$X0[as.numeric(str_sub(pc_second, 3, -1))] * 100, 1), nsmall = 1), " %")
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

    print(head(df_sym))
    print(head(df_sig))
    
    # Reshape to matrix for plotting
    sig_pair <- dcast(as.data.table(df_sig), spc1 ~ spc2, value.var = "category")
    print(sig_pair)
    
    # Get the names of the columns to modify (all except the first)
    cols_to_modify <- setdiff(names(sig_pair), "spc1")

    # Convert all columns except the first to matrix
    sig_pair_mat <- as.matrix(sig_pair[, ..cols_to_modify])
    

    # Set upper triangle to NA
    sig_pair_mat[upper.tri(sig_pair_mat)] <- NA

    # Assign back to the data.table
    sig_pair[, (cols_to_modify) := as.data.table(sig_pair_mat)]

    #sig_pair[upper.tri(as.matrix(sig_pair[,-1]))] <- NA
    print(sig_pair)
    sig_pair_melt <- melt(sig_pair, id.vars = "spc1")
    print(sig_pair_melt)

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
