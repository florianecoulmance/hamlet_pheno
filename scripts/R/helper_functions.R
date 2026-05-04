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
  # print(logos_spec)

  spec_colors <- setNames(species_col$Color, species_col$spec)
  # print(spec_colors)

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
# Function: legend_plot
# Purpose : Create a custom horizontal legend associating each species
#           with its color, logo, and italicized label.
# Input   :
#   - info_table: Data frame containing species information with columns:
#       - spec: species code
#       - Species: species name (without "H.")
#       - Color: color assigned to the species
#       - link: path or URL to the species logo image
# Output  :
#   - Returns a ggplot object displaying the legend horizontally, where:
#       • each species is represented by a colored dot above its logo
#       • species names are displayed below the logo in italics
# Notes   :
#   - Only selected species (aff, eco, esp, lib, ran) are displayed
#   - Layout parameters (dot/logo/text positions, n_cols) can be adjusted
# ============================================================
legend_plot <- function(info_table, gen = FALSE) {

  if (gen) {
    selected_specs <- c("atl", "cas", "eco", "esp", "flo", "gem", "lib")
  } else {
    selected_specs <- c("aff", "eco", "esp", "lib", "ran")
  }

  legend_info <- info_table %>%
    dplyr::filter(!(spec %in% selected_specs))

  color_map <- setNames(legend_info$Color, legend_info$spec)
  # label_map <- setNames(
  #   paste0("<img src='", legend_info$link, "' width='90' /><br>*H. ", legend_info$Species, "*"),
  #   legend_info$spec)

  n_cols <- 8   # number of columns per row (adjust as needed)
  n_cols <- ceiling(nrow(legend_info)/2)
  n_rows <- 2
  # print(n_rows)

  # x_spacing <- 1.2  # smaller = more compact columns

  legend_df <- legend_info %>%
    mutate(idx = row_number(),
          row = n_rows - ((idx - 1) %/% n_cols + 1),   # vertical position
          col = ((idx - 1) %% n_cols) * 0.4,         # horizontal position
          y_dot = row,                     # adjust vertical dot position
          y_logo = row,                          # logo y
          y_text = row - 0.3)                    # text y
  
  legend_plot <- ggplot(legend_df) +
  # colored dot
  geom_point(aes(x = col - 0.2, y = y_dot, color = spec), size = 10, show.legend = FALSE) +
  scale_color_manual(values = color_map) +
  # logo
  geom_image(aes(x = col, y = y_logo, image = link), size = 0.6, asp = 1.1) +
  # species name
  geom_text(aes(x = col - 0.1, y = y_text, label = paste0("H. ", Species)), size = 4.5, vjust = 1, fontface = "italic") +
  theme_void() +
  theme(plot.margin = margin(0.3,0,5,0)) +
  coord_cartesian(clip = "off")

  return(legend_plot)

}


# ============================================================
# Function: legend_geo
# Purpose : Create a custom horizontal legend associating each location
#           with its color, and label.
# Input   :
#   - info_table: Data frame containing location information with columns:
#       - geo: location code
#       - Locations: location name
#       - Color: color assigned to the species
# Output  :
#   - Returns a ggplot object displaying the legend horizontally, where:
#       • each location is represented by a colored dot
#       • location names are displayed next to the dot in italics
# Notes   :
#   - Only selected locations are displayed
#   - Layout parameters (dot/text positions, n_cols) can be adjusted
# ============================================================
legend_geo <- function(info_table, gen = FALSE) {

  if (gen) {
    selected_locs <- c("hon", "boc", "bel", "gun", "qui", "pri", "bar", "arc", "flk", "are", "ala")
    n_cols <- 11   # number of columns per row (adjust as needed)
  } else {
    selected_locs <- c("boc", "uvi", "bel", "flo", "tob")
    n_cols <- 5   # number of columns per row (adjust as needed)
  }

  legend_info <- info_table %>%
    dplyr::filter(geo %in% selected_locs)

  color_map <- setNames(legend_info$Color, legend_info$geo)
  # label_map <- setNames(
  #   paste0("<img src='", legend_info$link, "' width='90' /><br>*H. ", legend_info$Species, "*"),
  #   legend_info$spec)

  # n_cols <- ceiling(nrow(legend_info)/2)
  n_rows <- 1
  # print(n_rows)

  # x_spacing <- 1.2  # smaller = more compact columns

  legend_df <- legend_info %>%
    mutate(idx = row_number(),
          row = n_rows - ((idx - 1) %/% n_cols + 1),   # vertical position
          col = ((idx - 1) %% n_cols) * 0.5,         # horizontal position
          y_dot = row,                     # adjust vertical dot position
          y_logo = row,                          # logo y
          y_text = row)                    # text y
  
  legend_plot <- ggplot(legend_df) +
  # colored dot
  geom_point(aes(x = col - 0.15, y = y_dot, color = geo), size = 10, show.legend = FALSE) +
  scale_color_manual(values = color_map) +
  # logo
  # geom_image(aes(x = col, y = y_logo, image = link), size = 0.8, asp = 1.1) +
  # species name
  geom_text(aes(x = col, y = y_text, label = Locations), size = 4, hjust = 1, vjust = 2, fontface = "bold") +
  theme_void() +
  theme(plot.margin = margin(0.3,0,3,0)) +
  coord_cartesian(clip = "off")

  return(legend_plot)

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
  # print(extract_legend)
  # -----------------------------
  # 1. Choose color/grouping mode
  # -----------------------------
  color_by <- match.arg(color_by, c("species", "location"))
  
  if (color_by == "species") {
    group_col <- "spec"
    info_table <- species_info
    color_map <- setNames(info_table$Color, info_table$spec)
    label_map <- setNames(
      paste0("<img src='", info_table$link, "' width='90' /><br>*H. ", info_table$Species, "*"),
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
  # pca_data <- pca_data %>%
  #   mutate(
  #     # Fix specific naming error
  #     spec = ifelse(images == "PL17_160pueflo-l1-s4-f4-c2-d1.png",
  #                     "flo", spec))
  # print(pca_data$images)

  plot_data <- pca_data %>%
    left_join(info_table, by = group_col)
  
  # print(plot_data)

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
  # print(add_logos)

  # -----------------------------
  # 4. Build the actual PCA plot
  # -----------------------------
  p <- ggplot(plot_data, aes(x = .data[[pc_first]], y = .data[[pc_second]], color = .data[[group_col]])) +
    geom_point(size = 4, alpha = 1) +
    stat_ellipse(aes(color = .data[[group_col]]), linetype = 5, lwd = 1.5) +
    geom_point(data = centroids, aes(x = x, y = y, color = .data[[group_col]]), size = 15, alpha = 0) +
    # geom_image(data = centroids, aes(x = x, y = y, image = link), vjust=1, hjust=0, size = 0.15, asp = 1.1, alpha=1) +
    scale_color_manual(values = color_map, labels = label_map) +
    # theme_minimal() +
    theme(
      legend.position = if (extract_legend) "bottom" else "none",
      legend.box = "horizontal",
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 18),
      plot.margin = margin(3,3,0,0)
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

  # if(!is.null(add_logos)) p <- p + add_logos
  
  # print(class(p))


  # -----------------------------
  # 6. Add title (per location or per species)
  # -----------------------------
  # ---- Get title depending on color_by mode ----
  title_val <- NULL
  if (color_by == "species") {
    # if coloring by species, title is the location name
    geo_val <- unique(plot_data$geo)
    print(geo_val)
    print(length(geo_val))
    if (length(geo_val) == 1) {
      title_val <- geo_info$Locations[geo_info$geo == geo_val]
      if (title_val=="Panama") {
        title_val <- "(b) Panama"
      } else if (title_val=="USVI") {
        title_val <- "(f) USVI"
      } else if (title_val=="Belize") {
        title_val <- "(a) Belize"
      } else if (title_val=="Florida Keys") {
        title_val <- "(c) Florida Keys"
      } else if (title_val=="Tobago") {
        title_val <- "(b) Tobago"
      } else if (title_val=="Mexico") {
        title_val <- "(a) Mexico"
      } else if (title_val=="Honduras") { 
        title_val <- "(c) Honduras"
      } else if (title_val=="Puerto Rico") { 
        title_val <- "(d) Puerto Rico"
      } else {
        title_val <- ""
      }
    } else {
      title_val <- ""
    }
  } else if (color_by == "location") {
    # if coloring by location, title is the species name
    spec_val <- unique(plot_data$spec)
    # print(spec_val)
    if (length(spec_val) == 1) {
      title_val <- paste0("H. ", species_info$Species[species_info$spec == spec_val])
      if (title_val=="H. puella") {
        title_val <- "(a) H. puella"
      } else if (title_val=="H. nigricans") {
        title_val <- "(b) H. nigricans"
      } else if (title_val=="H. unicolor") {
        title_val <- "(c) H. unicolor"
      } else if (title_val=="H. chlorurus") {
        title_val <- "(d) H. chlorurus"
      } else if (title_val=="H. aberrans") {
        title_val <- "(e) H. aberrans"
      } else if (title_val=="H. indigo") {
        title_val <- "(f) H. indigo"
      } else {
        title_val <- ""
      }
    }
  } else {
    title_val <- ""
  }
  print(title_val)

  # -----------------------------
  # 7. Return plot
  # -----------------------------
  if (color_by == "species") {
    p_annot <- annotate_figure(
      p,
      top = text_grob(title_val, color = "black", face = "bold", size = 20, x = unit(0, "lines"), vjust=0, hjust=0)#, fig.lab.pos = "top.left"
    )
  } else if (color_by == "location" && pc_first=="PC1") {
      p_annot <- annotate_figure(
      p,
      top = text_grob(title_val, color = "black", face = "bold.italic", size = 20, x = unit(0, "lines"), vjust=0, hjust=0)#, fig.lab.pos = "top.left"
    )
  } else {
    p_annot <- annotate_figure(
      p,
      top = text_grob("   ",color = "black", face = "bold.italic", size = 20, x = unit(0, "lines"), vjust=0, hjust=0)
    )
  }
    
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
  # Ensure var is a data frame
  var <- as.data.frame(var)
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
         y = "Explained variance (%)"#,
        #  title = paste0("Explained variance - ", loc_name)
         ) +
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

  pc_table <- pc_table %>%
  mutate(
    # Fix specific naming error
    spec = ifelse(images == "PL17_160pueflo-l1-s4-f4-c2-d1.png",
                    "flo", spec)
  )
  
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
  print(pairwise_result)

  # ---- Reshape R² (upper triangle) ----
  vis_R <- pairwise_result %>%
    separate(pairs, into = c("group1","group2"), sep = " vs ") %>%
    select(group1, group2, R2)

  # Mirror the table so it’s symmetric
  vis_R_sym <- bind_rows(
    vis_R,
    vis_R %>% rename(group1 = group2, group2 = group1)
  )

  # Make sure all combinations exist
  all_R_groups <- unique(c(vis_R$group1, vis_R$group2))
  vis_R_complete <- vis_R_sym %>%
    complete(group1 = all_R_groups, group2 = all_R_groups, fill = list(R2 = NA))

  visH_R <- dcast(vis_R_complete, group1 ~ group2, value.var = "R2") %>%
    complete(group1 = levels(pc_table$group), fill = list()) %>%
    select(group1, intersect(levels(pc_table$group), colnames(.)))
  
  # ---- Reshape F.Model (lower triangle) ----
  vis_F <- pairwise_result %>%
    separate(pairs, into = c("group1","group2"), sep = " vs ") %>%
    select(group1, group2, F.Model)

  # Mirror the table so it’s symmetric
  vis_F_sym <- bind_rows(
    vis_F,
    vis_F %>% rename(group1 = group2, group2 = group1)
  )

  # Make sure all combinations exist
  all_F_groups <- unique(c(vis_F$group1, vis_F$group2))
  vis_F_complete <- vis_F_sym %>%
    complete(group1 = all_F_groups, group2 = all_F_groups, fill = list(F.Model = NA))

  visH_F <- dcast(vis_F_complete, group1 ~ group2, value.var = "F.Model") %>%
    complete(group1 = levels(pc_table$group), fill = list()) %>%
    select(group1, intersect(levels(pc_table$group), colnames(.))) #%>%
    # arrange(desc(group1))
  
  # ---- Reshape p-values (overlay text) ----
  vis_p <- pairwise_result %>%
    separate(pairs, into = c("group1","group2"), sep = " vs ") %>%
    select(group1, group2, p.adjusted)
  
  # Mirror the table so it’s symmetric
  vis_p_sym <- bind_rows(
    vis_p,
    vis_p %>% rename(group1 = group2, group2 = group1)
  )

  # Make sure all combinations exist
  all_p_groups <- unique(c(vis_p$group1, vis_p$group2))
  vis_p_complete <- vis_p_sym %>%
    complete(group1 = all_p_groups, group2 = all_p_groups, fill = list(p.adjusted = NA))

  visH_p <- dcast(vis_p_complete, group1 ~ group2, value.var = "p.adjusted") %>%
    complete(group1 = levels(pc_table$group), fill = list()) %>%
    select(group1, intersect(levels(pc_table$group), colnames(.))) #%>%
    # arrange(desc(group1))

  # ensure all matrices have same row/column order
  visH_R <- visH_R[match(visH_F$group1, visH_R$group1), c("group1", visH_F$group1)]
  visH_F <- visH_F[match(visH_R$group1, visH_F$group1), c("group1", visH_R$group1)]
  visH_p <- visH_p[match(visH_R$group1, visH_p$group1), c("group1", visH_R$group1)]
  
  # ---- Blank lower triangle (keep only upper) ----
  mat_R <- as.matrix(visH_R[,-1])
  rownames(mat_R) <- visH_R$group1
  mat_R[lower.tri(mat_R, diag = TRUE)] <- NA
  visH_R[,-1] <- mat_R

  melt_R <- melt(visH_R, id.vars = "group1")
  melt_R$value <- as.numeric(melt_R$value)

  mat_F <- as.matrix(visH_F[,-1])
  rownames(mat_F) <- visH_F$group1
  mat_F[upper.tri(mat_F, diag = TRUE)] <- NA
  visH_F[,-1] <- mat_F

  melt_F <- melt(visH_F, id.vars = "group1")
  melt_F$value <- as.numeric(melt_F$value)

  mat_p <- as.matrix(visH_p[,-1])
  rownames(mat_p) <- visH_p$group1
  mat_p[upper.tri(mat_p, diag = TRUE)] <- NA
  visH_p[,-1] <- mat_p

  melt_p <- melt(visH_p, id.vars = "group1")
  melt_p$value <- as.numeric(melt_p$value)
  
  # ---- Plot heatmap ----
  p <- ggplot() +
    geom_tile(data = melt_R, aes(x = group1, y = variable, fill = value), color = "transparent") +
    scale_fill_gradient(low = "#ffead1", high = "#fdae53", na.value = "transparent", name = "R²") +
    new_scale_fill() +
    geom_tile(data = melt_F, aes(x = group1, y = variable, fill = value), color = "transparent") +
    scale_fill_gradient(low = "#ffedec", high = "#ff5c52", na.value = "transparent", name = "F") +
    geom_text(data = melt_p, aes(x = group1, y = variable,
          label = ifelse(is.na(value), "", ifelse(value < 0.001, "***", ifelse(value < 0.01, "**", ifelse(value < 0.05, "*", "ns"))))),
      size = 10, color = "black", fontface="bold") +
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
    ) +
    guides(color = guide_legend(ncol = 1))
  
  # ---- Get title ----
  if(color_by == "species"){
    geo_val <- unique(pc_table$geo)
    title_val <- if(length(geo_val) == 1) geo_map$Locations[geo_map$geo == geo_val] else ""
    if (title_val=="Panama") {
        title_val <- "(e) Panama"
      } else if (title_val=="USVI") {
        title_val <- "(f) USVI"
      } else if (title_val=="Belize") {
        title_val <- "(d) Belize"
      } else if (title_val=="Florida Keys") {
        title_val <- "(c) Florida Keys"
      } else if (title_val=="Tobago") {
        title_val <- "(b) Tobago"
      } else if (title_val=="Mexico") {
        title_val <- "(a) Mexico"
      } else {
        title_val <- "(a)"
      }
  } else if (color_by == "location"){
    species_val <- unique(pc_table$spec)
    title_val <- "" #if (length(species_val) == 1) paste0("H. ", species_col$Species[species_col$spec == species_val]) else ""
  }
  
  # ---- Annotate with location title ----
  p_annot <- annotate_figure(
    p,
    top = text_grob(title_val, color = "black", face = "bold", size = 20, x = unit(0, "lines"), vjust=0, hjust=0)#, fig.lab.pos = "top.left"
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
  # tree$x[tree$isTip] <- tree$x[tree$isTip] * 1.1
  # tree$branch.length <- scale(tree$branch.length)
  
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
  t <- ggtree(hc, layout = "fan", size = 0.5, color = scol) %<+% tree +
    geom_tippoint(aes(color = .data[[group_col]]), size = 3, alpha = 0.5) +
    scale_color_manual(values = color_map, name = legend_name) +
    theme(
      legend.position = ifelse(extract_legend, "right", "none"),
      plot.title = element_text(size = 14, color = "gray20", face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray20")
    )

  # Extend terminal branches (multiply x for tip nodes)
  t$data$x <- ifelse(t$data$isTip, t$data$x * 1.2, t$data$x)
  
  # -----------------------------
  # 6. Get title
  # -----------------------------
  title_val <- NULL

  if (!(str_detect(pca_file, "lab_571_"))) {
    if(color_by == "species"){
      geo_val <- unique(tree$geo)[1]
      title_val <- geo_map$Locations[geo_map$geo == geo_val]
      if (title_val=="Panama") {
        title_val <- "(e) Panama"
      } else if (title_val=="USVI") {
        title_val <- "(f) USVI"
      } else if (title_val=="Belize") {
        title_val <- "(d) Belize"
      } else if (title_val=="Florida Keys") {
        title_val <- "(c) Florida Keys"
      } else if (title_val=="Tobago") {
        title_val <- "(b) Tobago"
      } else if (title_val=="Mexico") {
        title_val <- "(a) Mexico"
      } else {
        title_val <- ""
      }
    } else if (color_by == "location"){
      species_val <- unique(tree$spec)[1]
      title_val <- "" #if (length(species_val) == 1) paste0("H. ", species_col$Species[species_col$spec == species_val]) else ""
    }
  } else {
    title_val <- NULL
  }

  # print(title_val)
  
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
      top = text_grob(title_val, color = "black", face = "bold", size = 20, x = unit(0, "lines"), vjust=0, hjust=0)#, fig.lab.pos = "top.left"
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
    grob <- rasterGrob(img[1:500,100:1100,], interpolate = TRUE)
    grob <- annotate_figure(grob, top = text_grob(
      pc,
      color = "black", face = "bold", size = 20,
      x = unit(2, "pt"), hjust = -1, vjust = 1
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
      if (title_val=="Panama") {
        title_val <- "(e) Panama"
      } else if (title_val=="USVI") {
        title_val <- "(f) USVI"
      } else if (title_val=="Belize") {
        title_val <- "(d) Belize"
      } else if (title_val=="Florida Keys") {
        title_val <- "(c) Florida Keys"
      } else if (title_val=="Tobago") {
        title_val <- "(b) Tobago"
      } else if (title_val=="Mexico") {
        title_val <- "(a) Mexico"
      } else {
        title_val <- ""
      }
    } else {
      title_val <- ""
    }
    
  } else if (color_by == "location") {
    # Within species across locations → use full species name
    if (abbrev %in% spec_map$spec) {
      title_val <- "" #paste0("H. ", spec_map$Species[match(abbrev, spec_map$spec)])
    } else {
      title_val <- ""
    }
    
  } else {
    stop("Invalid value for color_by. Use 'species' or 'location'.")
  }
  
  # add title
  combined <- annotate_figure(
    combined,
    top = text_grob(title_val, color = "black", face = "bold", size = 20, x = unit(0, "lines"), vjust=0, hjust=0)#, fig.lab.pos = "top.left"
  )
  
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
    # print(allsamples)

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
    # print(head(pca_coords))
    # print(head(var_explained))

    # Ensure the PCA coordinates rows are in the same order as the original samples
    #pca_coords <- pca_coords[match(allsamples, rownames(pca_coords)), ]
    #print(head(pca_coords))

    # Then add spec, geo, and sample
    pca_coords$spec   <- spec
    pca_coords$geo    <- geo
    pca_coords$sample <- allsamples
    # print(head(pca_coords))    

    # Extract the "variance explained" row as numeric
    var_explained <- as.numeric(var_explained["variance explained", ])
    # print(var_explained)

    # If currently in percentages, convert to fractions
    var_explained <- var_explained / 100  # remove if already in 0-1 range

    # Create the table in the same format as read.csv would
    table <- data.frame(X0 = var_explained)
    # print(table)
    # Add row numbers starting from 0 if needed
    rownames(table) <- 0:(length(var_explained)-1)
    # print(table)

    # Access the values like
    # (print(table$X0))

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
    
    pca_data$boc <- factor(ifelse(pca_data$geo == "boc", "Panama", "other"),
                           levels = c("Panama", "other"))

    # Always species mode
    group_col <- "spec"
    info_table <- species_info
    color_map <- setNames(info_table$Color, info_table$spec)
    label_map <- setNames(
        paste0("<img src='", info_table$link, "' width='70' /><br>*H. ", info_table$Species, "*"),
        info_table$spec
    )


    # Build PCA plot
    p <- ggplot(pca_data, aes(x = .data[[pc_first]], y = .data[[pc_second]])) +
        geom_point(aes(color = .data[[group_col]], shape = boc), size = 5, alpha = 0.8) +
        scale_color_manual(
          values = color_map,
          labels = label_map,
          limits = names(color_map),
          breaks = names(color_map),
          drop = FALSE
          ) +
        scale_shape_manual(
          values = c(Panama = 8, other = 19),
          limits = c("Panama", "other"),
          breaks = c("Panama", "other"),
          drop = FALSE
          ) + # filled circle for Panama, star/other for rest
        theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.text = if (pc_first %in% c("PC1", "PC3", "PC5")) element_blank() else element_markdown(size = 12),
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.margin = margin(0,0,3,0)
        ) +
        labs(
        x = paste0(pc_first, ", variance = ", format(round(variance$X0[as.numeric(str_sub(pc_first, 3, -1))] * 100, 1), nsmall = 1), " %"),
        y = paste0(pc_second, ", variance = ", format(round(variance$X0[as.numeric(str_sub(pc_second, 3, -1))] * 100, 1), nsmall = 1), " %")
        ) +
        guides(
        color = guide_legend(ncol = 1, byrow = TRUE, override.aes = if (pc_first %in% c("PC1", "PC3", "PC5")) list(alpha = 0, size = 1) else list()),
        shape = guide_legend(ncol = 1, override.aes = if (pc_first %in% c("PC1", "PC3", "PC5")) list(alpha = 0, size = 1) else list())
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
plot_permanova_permdisp <- function(pair_file, species_col, geo_map, color_by = "species", params_legend = "none") {
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

    # print(head(df_sym))
    # print(head(df_sig))
    
    # Reshape to matrix for plotting
    sig_pair <- dcast(as.data.table(df_sig), spc1 ~ spc2, value.var = "category")
    # print(sig_pair)
    
    # Get the names of the columns to modify (all except the first)
    cols_to_modify <- setdiff(names(sig_pair), "spc1")

    # Convert all columns except the first to matrix
    sig_pair_mat <- as.matrix(sig_pair[, ..cols_to_modify])
    

    # Set upper triangle to NA
    sig_pair_mat[upper.tri(sig_pair_mat)] <- NA

    # Assign back to the data.table
    sig_pair[, (cols_to_modify) := as.data.table(sig_pair_mat)]

    #sig_pair[upper.tri(as.matrix(sig_pair[,-1]))] <- NA
    # print(sig_pair)
    sig_pair_melt <- melt(sig_pair, id.vars = "spc1")
    # print(sig_pair_melt)

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
        legend.text = element_text(size = 13, margin = margin(0, 53, 0, 0)),
        legend.title = element_text(size = 16),
        legend.key.height = unit(1.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0,5,0,0),
        axis.text.x = element_text(angle = 0, size = 16),
        axis.text.y = element_text(size = 16),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1
        )

    # ---- Get title ----
    print(pair_file)
    if(color_by == "species"){
      geo_val <- substr(pair_file, nchar(pair_file) - 19 + 1, nchar(pair_file) - 17 + 1)
      print(geo_val)
      title_val <- if(length(geo_val) == 1) geo_map$Locations[geo_map$geo == geo_val] else ""
      if (length(title_val) == 1 && title_val=="Panama") {
        title_val <- "(b) Panama"
      } else if (length(title_val) == 1 && title_val=="Belize") {
        title_val <- "(a) Belize"
      } else if (length(title_val) == 1 && title_val=="Honduras") {
        title_val <- "(c) Honduras"
      } else if (length(title_val) == 1 && title_val=="Puerto Rico") {
        title_val <- "(d) Puerto Rico"
      } else {
        title_val <- ""
      }
    } else if(color_by == "location"){
      species_val <- substr(pair_file, nchar(pair_file) - 19 + 1, nchar(pair_file) - 17 + 1)
      print(species_val)
      title_val <- "" #if (length(species_val) == 1) paste0("H. ", species_col$Species[species_col$spec == species_val]) else ""
    }
    print(title_val)

    # ---- Annotate with location title ----
    p_annot <- annotate_figure(
      p_allLoc,
      top = text_grob(title_val, color = "black", face = "bold", size = 20, x = unit(0, "lines"), vjust=0, hjust=0)#, fig.lab.pos = "top.left"
    )

    return(p_annot)
}


# ============================================================
# Function: fst_analysis
# Purpose : Compute pairwise Weir & Cockerham FST values from a 
#           gtraw genotype file, filtering out groups with <3 samples.
# Input   :
#   - gtfile      : Path to gtraw file.
#   - color_by    : "species" or "location" (defines populations).
#   - species_col : Species lookup table.
#   - geo_map     : Location lookup table.
# Output  :
#   - Pairwise FST matrix (optionally reshaped for plotting).
# ============================================================
fst_analysis <- function(gtfile, color_by, species_col, geo_map, label) {
  
  # ---- Load genotype data ----
  gtraw <- read.table(gtfile, header = TRUE, check.names = FALSE)
  geno <- t(gtraw[, 7:ncol(gtraw)])   # genotypes start at column 7
  samples <- colnames(gtraw)[7:ncol(gtraw)]

  # ---- Extract species & location ----
  samples_clean <- ifelse(
    grepl("PL17", samples),
    sapply(strsplit(samples, "_"), function(x) paste(x[1:2], collapse = "_")),
    sub("_.*", "", samples)
  )
  # print(samples_clean)
  spec <- substr(samples_clean, nchar(samples_clean) - 5, nchar(samples_clean) - 3)
  # print(spec)
  geo <- substr(samples_clean, nchar(samples_clean) - 2, nchar(samples_clean))
  # print(geo)
  # print(head(data.frame(samples_clean, spec, geo)))

  pop <- if (color_by == "species") spec else geo
  # print(pop)

  # ---- Filter out populations with <3 individuals ----
  pop_counts <- table(pop)
  keep_pops <- names(pop_counts[pop_counts >= 3])
  # print(keep_pops)
  if (length(keep_pops) < 2) {
    message("Not enough populations with ≥3 individuals to compute pairwise FST.")
    return(NULL)
  }

  keep_idx <- pop %in% keep_pops
  pop <- pop[keep_idx]
  geno <- geno[keep_idx, , drop = FALSE]
  samples_clean <- samples_clean[keep_idx]
  # print(samples_clean)

  # Assign rownames
  rownames(geno) <- samples_clean

  # ---- Build sample_groups tibble for pairwise_fst ----
  sample_groups <- tibble(
    sample = samples_clean,
    group  = pop
  )

  # print(sample_groups)

  # ---- Compute pairwise FST ----
  fst_mat <- pairwise_fst(geno, sample_groups)
  # print(fst_mat$Fst)
  # print(fst_mxat)
  print(paste("label inside function =", label))
  FST_RESULTS[[label]] <<- fst_mat$Fst
  print(names(FST_RESULTS))

  # ---- Reshape for plotting ----
  fst_dt <- as.data.table(fst_mat$Fst)
  # print(fst_dt)
  # Remove diagonal
  fst_dt[pop1 == pop2, Fst := NA]
  # print(fst_dt)

  # Build plot
  p <- ggplot() +
    geom_tile(data = fst_dt, aes(x = pop1, y = pop2, fill = Fst), color = "transparent") +
    scale_fill_gradient(low = "#F3D6F3", high = "#A964B7", na.value = "transparent", name = "FST") +
    geom_text(data = fst_dt,
              aes(x = pop1, y = pop2, label = round(Fst, 3)),
              size = 3, color = "black", fontface="bold") +
    labs(x = "", y = "", fill = "FST") +
    scale_x_discrete(position = "top") +
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
    ) +
    guides(color = guide_legend(ncol = 1))

  # ---- Get title ----
  if(color_by == "species"){
    # print(geo)
    geo_val <- unique(geo)
    title_val <- if(length(geo_val) == 1) geo_map$Locations[geo_map$geo == geo_val] else ""
    if (title_val=="Panama") {
        title_val <- "(b) Panama"
      } else if (title_val=="Cayos Arcas") {
        title_val <- "(e) Cayos Arcas"
      } else if (title_val=="Belize") {
        title_val <- "(a) Belize"
      } else if (title_val=="Florida Keys") {
        title_val <- "(g) Florida Keys"
      } else if (title_val=="Barbados") {
        title_val <- "(f) Barbados"
      } else if (title_val=="Guna Yala") {
        title_val <- "(h) Guna Yala"
      } else if (title_val=="Honduras") { 
        title_val <- "(c) Honduras"
      } else if (title_val=="Puerto Rico") { 
        title_val <- "(d) Puerto Rico"
      } else if (title_val=="Quintana Roo") {
        title_val <- "(i) Quintana Roo"
      } else {
        title_val <- ""
      }
  } else if(color_by == "location"){
    # print(spec)
    species_val <- unique(spec)
    title_val <- if (length(species_val) == 1) paste0("H. ", species_col$Species[species_col$spec == species_val]) else ""
    if (title_val=="H. puella") {
        title_val <- "(a) H. puella"
      } else if (title_val=="H. nigricans") {
        title_val <- "(b) H. nigricans"
      } else if (title_val=="H. unicolor") {
        title_val <- "(c) H. unicolor"
      } else if (title_val=="H. chlorurus") {
        title_val <- "(f) H. chlorurus"
      } else if (title_val=="H. aberrans") {
        title_val <- "(d) H. aberrans"
      } else if (title_val=="H. indigo") {
        title_val <- "(i) H. indigo"
      } else if (title_val=="H. affinis") {
        title_val <- "(e) H. affinis"
      } else if (title_val=="H. gemma") {
        title_val <- "(g) H. gemma"
      } else if (title_val=="H. gummigutta") {
        title_val <- "(h) H. gummigutta"
      } else if (title_val=="H. sp1") {
        title_val <- "(j) H. sp1"
      } else {
        title_val <- ""
      }
  }
  print(title_val)

  # ---- Annotate with location title ----
  if (color_by == "species") {
    p_annot <- annotate_figure(
      p,
      top = text_grob(title_val, color = "black", face = "bold", size = 20, x = unit(0, "lines"), vjust=0, hjust=0)#, fig.lab.pos = "top.left"
    )
  } else if (color_by == "location") {
      p_annot <- annotate_figure(
      p,
      top = text_grob(title_val, color = "black", face = "bold.italic", size = 20, x = unit(0, "lines"), vjust=0, hjust=0)#, fig.lab.pos = "top.left"
    )
  }

  return(p_annot)

}



# ============================================================
# Function: make_summary_by_location
# Purpose : Summarize sample counts per location, compute scaled values 
#           (e.g., for pie plot radii), and merge coordinate and color metadata.
# Input   :
#   - sample_data   : Data frame containing sample-level information with columns 'geo' and 'spec'.
#   - location_data : Data frame containing location metadata with columns 
#                     'geo', 'Locations', 'coord_N.x', 'coord_W.x', and optionally 'Color'.
#   - scale_factor  : Numeric scaling factor used to compute the radius (default = 0.9).
# Output  :
#   - A data frame containing one row per location, counts of each category ('spec'),
#     total count, scaled radius values, and merged coordinates/colors.
# ============================================================
make_summary_by_location <- function(sample_data, location_data, scale_factor = 0.9) {
  # Step 1: join to add Locations info
  sample_data <- sample_data %>%
    left_join(location_data %>% select(geo, Locations), by = "geo")
  
  # Step 2: summarize category counts per Location
  summary_dat <- dcast(
    sample_data,
    Locations ~ spec,
    fill = 0,
    value.var = "spec",
    fun.aggregate = length
  ) %>%
    mutate(
      total = rowSums(across(-Locations)),
      radius = log10(total) * scale_factor
    )
  
  # Step 3: add coordinates and color (if available)
  summary_dat <- summary_dat %>%
    left_join(
      location_data %>% select(Locations, coord_N.x, coord_W.x),
      by = "Locations"
    )
  
  return(summary_dat)
}


# ============================================================
# Function: custom_geom_scatterpie_legend
# Purpose : Create a custom legend for scatterpie plots, allowing finer control
#           over pie radius scaling and label text size for geographic pie charts.
# Input   :
#   - radius   : Numeric vector of pie radii to include in the legend.
#   - x, y     : Numeric coordinates specifying the position of the legend.
#   - n        : Number of legend items to display (default = 5).
#   - labeller : Optional function to convert radius values into label text.
#   - textsize : Numeric value controlling the size of the legend text labels (default = 6).
# Output  :
#   - A list of ggplot2 layers (arcs, segments, text) forming a custom legend.
# Notes   :
#   - Adapted from an example on the RStudio Community forum:
#     https://community.rstudio.com/t/pie-chart-world-map-for-genetics/97192/2
# ============================================================
custom_geom_scatterpie_legend <- function(radius, x, y, n = 5, labeller, textsize = 6) {
  
  # If more radii than desired, sample evenly spaced unique values within range
  if (length(radius) > n) {
    radius <- unique(sapply(
      seq(min(radius), max(radius), length.out = n),
      scatterpie:::round_digit
    ))
  }
  
  label <- FALSE  # default: no custom label function
  
  # If a custom labeller function is provided, validate and activate labeling
  if (!missing(labeller)) {
    if (!inherits(labeller, "function")) {
      stop("labeller should be a function for converting radius")
    }
    label <- TRUE
  }
  
  # Build a dataframe for legend components
  dd <- data.frame(
    r = radius,
    start = 0,
    end = 2 * pi,
    x = x,
    y = y + radius - max(radius),
    maxr = max(radius)
  )
  
  # Apply label transformation if provided
  dd$label <- if (label) labeller(dd$r) else dd$r
  
  # Return list of ggplot2 layers (arc bars, guide lines, and text labels)
  list(
    ggforce:::geom_arc_bar(
      aes_(x0 = ~x, y0 = ~y, r0 = ~r, r = ~r, start = ~start, end = ~end),
      data = dd,
      inherit.aes = FALSE
    ),
    geom_segment(
      aes_(x = ~x, xend = ~x + maxr * 1.5, y = ~y + r, yend = ~y + r),
      data = dd,
      inherit.aes = FALSE
    ),
    geom_text(
      aes_(x = ~x + maxr * 1.6, y = ~y + r, label = ~label),
      data = dd,
      hjust = "left",
      inherit.aes = FALSE,
      size = textsize
    )
  )
}


# ============================================================
# Function: plot_location_pies
# Purpose : Plot scatterpie charts of sample distributions on a world map,
#           with colors/logos by species and customizable site labels.
# Input   :
#   - data_table     : Data frame with one row per site, species counts in columns,
#                      columns for Site, coord_N.x, coord_W.x, radius
#   - species_info   : Table mapping species codes to colors and logos
#   - radius_factor  : Scaling factor used when computing radius (default 0.9)
#   - label_column   : Column name in data_table to use for pie chart labels
#   - map_limits     : Vector of c(xmin, xmax, ymin, ymax) to set map limits
#   - show_legend    : Logical, whether to display the legend
# Output  :
#   - ggplot object of the world map with pie charts
# ============================================================
plot_location_pies <- function(data_table, species_info, radius_factor = 0.9,
                               map_limits = c(-100.15, -55.12, 6.00, 31.97),
                               show_legend = FALSE) {
  
  # Load world map
  world <- map_data('world')

  info_table <- species_info
  color_map <- setNames(info_table$Color, info_table$spec)
  label_map <- setNames(
      paste0("<img src='", info_table$link, "' width='120' /><br>*H. ", info_table$Species, "*"),
      info_table$spec
  )
  
  # Identify species columns (numeric counts)
  species_cols <- colnames(data_table)[!(colnames(data_table) %in% c("Locations", "coord_N.x", "coord_W.x", "radius", "total"))]
  
  # # Compute radius if not already present
  # if(!"radius" %in% colnames(data_table)) {
  #   data_table <- data_table %>%
  #     mutate(radius = log10(rowSums(across(all_of(species_cols)))) * radius_factor)
  # }
  
  # Base map
  p <- ggplot(data = world, aes(long, lat)) +
    geom_map(map = world, aes(map_id = region), fill = "grey", color = "#F6F6F6", alpha = 0.8) +
    
    # Scatterpie
    geom_scatterpie(aes(x = coord_W.x, y = coord_N.x, group = Locations, r = radius),
                    data = data_table, cols = species_cols, color = NA, alpha = 1) +
    
    # Optional legend for pie sizes
    custom_geom_scatterpie_legend(data_table$radius, x = -64, y = 26, n = 4,
                                   labeller = function(x) round(10^(x / radius_factor), digits = 0)) +
    
    # Site labels (flexible column)
    geom_text(aes(x = coord_W.x, y = coord_N.x, label = Locations,
                  hjust = ifelse(.data$Locations %in% c("Quintana Roo", "Alacranes Reef", "San Andrés", "Guna Yala"), NA, 1.2),
                  vjust = ifelse(.data$Locations %in% c("Cayo Arenas"), NA, -2)),
              data = data_table, color = "grey20", size = 6, fontface = "italic", position = position_dodge(width = 1)) +
    
    # Map annotations
    annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", fontface = "italic", color = "grey", size = 4) +
    annotate(geom = "text", x = -77, y = 15, label = "Caribbean Sea", fontface = "italic", color = "grey", size = 4) +
    annotate(geom = "text", x = -73, y = 30, label = "Atlantic", fontface = "italic", color = "grey", size = 4) +
    annotate(geom = "text", x = -64, y = 30, label = "Sample size", fontface = "bold", color = "black", size = 6) +

    # Map limits
    coord_sf(xlim = map_limits[1:2], ylim = map_limits[3:4], expand = FALSE) +
    
    # Species colors + logos
    scale_fill_manual(values = color_map,
                      labels = label_map) +
    
    # Theme
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      legend.text = element_markdown(size = 8),
      legend.position = if(show_legend) "bottom" else "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  return(p)
}


# ============================================================
# Function: plot_species_geo_overview
# Purpose : Create a tally plot of species counts per location/geo,
#           displaying species logos and names, ordered by region and species groups.
# Input   :
#   - dat           : A data frame containing at least columns 'spec' (species) and 'geo' (location codes).
#   - species_info  : A data frame containing species metadata with columns:
#                     'spec' (species code), 'Species' (full name), 'Color' (hex color), and 'link' (logo image path).
# Output  :
#   - A ggplot object showing counts of each species per location with logos and color coding.
# ============================================================
plot_species_geo_overview <- function(dat, species_info, geo_table) {
  
  # Use your species_info to create color and label maps
  color_map <- setNames(species_info$Color, species_info$spec)
  label_map <- setNames(
    paste0("<img src='", species_info$link, "' width='30' /><br>*H. ", species_info$Species, "*"),
    species_info$spec
  )



  # Define groups
  ingroup <- c("abe", "aff", "atl", "cas", "chl", "eco", "esp", "flo", "gem", "gum", 
              "gut", "ind", "lib", "may", "nig", "pro", "pue", "ran", "tan", "uni")
  outgroup <- c("Total")

  gulf <- c("liz", "tam", "ala", "arc", "are", "flk", "ver")
  west_carib <- c("bel", "boc", "gun", "hon", "san", "qui")
  east_carib <- c("bar", "hai", "pri", "uvi", "tob")
  
  # Summarize counts
  counts <- dat %>%
    dplyr::count(spec, geo) %>%
    replace(is.na(.), 0) %>%
    arrange(spec) %>%
    mutate(geo = ifelse(geo == "flo", "flk", geo),
           geo = ifelse(geo == "por", "pri", geo))
  
  # Add totals by species automatically
  test1 <- counts %>%
    split(.$spec) %>%
    purrr::map_df(~ janitor::adorn_totals(., where = "row")) %>%
    # group_by(spec) %>%
    mutate(geo = ifelse(spec == "Total", "Total", geo),
           spec = ifelse(spec == "Total", NA, spec)) %>%
    fill(spec, .direction = "down") #%>%
    # ungroup()
    # print(test1)
  
  # Add totals by location automatically
  test2 <- counts %>%
    split(.$geo) %>%
    purrr::map_df(~ janitor::adorn_totals(., where = "row")) %>%
    # group_by(geo) %>%
    mutate(geo = ifelse(geo == "-", NA, geo)) %>%
    fill(geo, .direction = "down") #%>%
    # ungroup()
  # print(test2)

  t <- unique(rbind(test1, test2))
  # print(t)

  sum <- textGrob(sum(counts$n), gp = gpar(col="black", fontsize = 22, fontface = "bold"))

  # Order species within groups
  order_sp <- t %>%
    mutate(group_sp = case_when(
      spec %in% ingroup ~ "ingroup",
      spec %in% outgroup ~ "outgroup"
    )) %>%
    dplyr::count(spec, group_sp) %>%
    group_by(group_sp) %>%
    arrange(desc(n), .by_group = TRUE) %>%
    pull(spec)

  # print(order_sp)
  
  # Order locations within regions
  order_loc1 <- t %>%
    mutate(Region = case_when(
      geo %in% gulf ~ "Gulf of Mexico",
      geo %in% west_carib ~ "Western Caribbean",
      geo %in% east_carib ~ "Eastern Caribbean",
      geo %in% "Total" ~ " "
    ))
    
  order_loc <- order_loc1 %>%
    dplyr::count(geo, Region) %>%
    group_by(Region) %>%
    arrange(match(Region, c("Western Caribbean", "Eastern Caribbean", "Gulf of Mexico", " ")),
            desc(n)) %>%
    pull(geo)

  # Compute vertical line positions dynamically
  region_boundaries <- order_loc1 %>%
    filter(Region != " ") %>%
    distinct(geo, Region) %>%
    mutate(geo_num = as.numeric(factor(geo, levels = order_loc))) %>%
    group_by(Region) %>%
    summarise(max_x = max(geo_num)) %>%
    ungroup() %>%
    # mutate(vline_x = max_x + 0.5) %>%
    pull(max_x)

  # print(max(region_boundaries))

  region_labels <- order_loc1 %>%
    filter(Region != " ") %>%
    distinct(geo, Region) %>%
    mutate(geo_num = as.numeric(factor(geo, levels = order_loc))) %>%
    group_by(Region) %>%
    summarise(mid_x = mean(geo_num)) %>%
    ungroup()
  # print(region_labels)

  # Keep only the locations actually in your dataset
  geo_table_filtered <- geo_table %>%
    filter(geo %in% t$geo) %>%
    arrange(match(geo, order_loc))  # optional: match the plotting order

  # Create a named vector for labels
  x_labels <- setNames(geo_table_filtered$Locations, geo_table_filtered$geo)
  print(x_labels)
  
  # Create plot
  p <- t %>%
    mutate(
      geo = fct_relevel(geo, order_loc),
      spec = fct_relevel(spec, rev(order_sp))
    ) %>%
    ggplot(aes(x = geo, y = spec)) +
    geom_count(aes(size = n, color = spec), show.legend = FALSE) +
    geom_text(aes(label = n), size = 6, nudge_x = 0.4, color = "gray50") +
    geom_vline(xintercept = region_boundaries + 0.6, # remove last
             col = "gray80", linetype = "dashed") +
    geom_hline(yintercept = 1.5, col = "gray80", linetype = "dashed") +
    scale_color_manual(values = color_map) +
    scale_size_area(max_size = 20) +
    # scale_y_discrete(labels = label_map) +
    scale_x_discrete(labels = function(geo) {
                                ifelse(
                                  is.na(x_labels[geo]), 
                                  "Total",   # keep NA labels blank
                                  ifelse(
                                    seq_along(x_labels[geo]) %% 2 == 0,
                                    paste0("\n", x_labels[geo]),
                                    paste0(x_labels[geo], "\n")
                                  )
                                )
                              }
                    ) +  # <-- this maps geo codes to full location names
    coord_cartesian(clip = "off") +
    labs(title = NULL, x = NULL, y = NULL) +
    # annotate region names dynamically
    geom_text(data = region_labels,
            aes(x = mid_x, y = length(unique(t$spec)) + 0.8, label = Region),
            color = "gray20", size = 7) +
    annotation_custom(sum, xmin = max(region_boundaries) + 1, xmax = max(region_boundaries) + 1, ymin = 1, ymax = 1) +
    # annotate(geom = "text", 
    #          x = c(1.5, 4, 6.8), 
    #          y = c(16.5, 16.5, 16.5), 
    #          label = c("Western Caribbean", "Eastern Caribbean", "Gulf of Mexico"),
    #          color = "gray20") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(face = "bold", size = 15),
      axis.text.y = element_blank(),
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.25, l = 1.5, unit = "cm")
    )

  # Filter only the species actually present in t
  species_info_subset <- species_info %>%
  filter(spec %in% unique(t$spec))

  p <- p +
  geom_image(
    data = species_info_subset,
    aes(x = -0.95, y = spec, image = link),
    inherit.aes = FALSE,
    size = 0.08, by = "width"
  ) +
  geom_text(
    data = species_info_subset,
    aes(x = -0.85, y = spec, label = paste0("H. ", Species)),
    hjust = 0, vjust = 2, size = 5, fontface = "italic"
  )

  return(p)                                    

}


# ---------------------------------------------------------
# GET RESULT PATHS (NEW STRUCTURE)
# ---------------------------------------------------------
get_nh_result_paths <- function(base_dir) {
  list.dirs(base_dir, recursive = TRUE, full.names = TRUE) %>%
    keep(~ grepl("NH.Results/.+_Results$", .x))
}

# ---------------------------------------------------------
# EXTRACT PofZ
# ---------------------------------------------------------
read_pofz <- function(path) {

  pofz_file <- list.files(path, pattern = "PofZ.txt", full.names = TRUE)
  print(pofz_file)
  indiv_file <- list.files(path, pattern = "individuals.txt", full.names = TRUE)
  print(indiv_file)

  if (length(pofz_file) == 0 | length(indiv_file) == 0) {
    return(NULL)
  }

  # extract pair + location
  pair_dir <- dirname(dirname(dirname(path)))
  print(pair_dir)

  runname <- basename(pair_dir)
  print(runname) 



  loc <- basename(dirname(dirname(dirname(dirname(path)))))
  print(loc)

  pops <- str_split(runname, "_")[[1]]
  print(pops)


  df <- vroom::vroom(pofz_file,
                     delim = "\t",
                     skip = 1,
                     col_names = c("indNR", "IndivName",
                                   "P1", "P2", "F1", "F2",
                                   "P1_bc", "P2_bc"))

  inds <- readLines(indiv_file)

  df$IndivName <- inds
  # if (length(inds) != nrow(df)) {
  #   warning("Mismatch in individuals vs PofZ rows: ", indiv_file)
  #   return(NULL)
  #   }

  df <- df %>%
    pivot_longer(cols = -c(indNR, IndivName),
                 names_to = "class",
                 values_to = "prob") %>%
    mutate(
      run = runname,
      pop1 = pops[1],
      pop2 = pops[2],
      loc = loc
    )

  print(df)
  return(df)

}

# ---------------------------------------------------------
# COLOR PALETTE (IDENTICAL TO OLD FIGURE)
# ---------------------------------------------------------
get_hybrid_colors <- function() {
  clr <- paletteer_c("ggthemes::Red-Green-Gold Diverging", 3) %>%
    c(., clr_lighten(.)) %>%
    color()

  clr[c(1,4,2,5,6,3)] %>%
    set_names(c("P1","P1_bc","F1","F2","P2_bc","P2"))
}

# ---------------------------------------------------------
# PLOT ONE LOCATION
# ---------------------------------------------------------
plot_location <- function(df, species_meta) {

  # detect hybrids
  hybrids <- df %>%
    filter(class %in% c("F1", "F2"), prob > 0.99) %>%
    pull(IndivName) %>%
    unique()

  df <- df %>%
      mutate(
        spec = stringr::str_extract(IndivName, "[a-z]{3}"),
      ) %>%
      arrange(spec, IndivName)
  
  df$ind_label <- factor(df$IndivName, levels = unique(df$IndivName))
  
  print(df)

  color_map <- setNames(species_meta$Color, species_meta$spec)
  species_labels <- setNames(paste0("*H. ", species_meta$Species, "*"), species_meta$spec)

  # labels with logos (like PCA)
  label_map <- setNames(
    paste0(
      "<img src='", species_meta$link, "' width='60'/><br>*H. ",
      species_meta$Species, "*"
    ),
    species_meta$spec
  )

  df <- df %>%
    mutate(
      run_label = paste0(
        label_map[pop1],
        " - ",
        label_map[pop2]
      )
    )

  colors <- get_hybrid_colors()

  df$class <- factor(df$class,
                   levels = c("P1", "P1_bc", "F1", "F2", "P2_bc", "P2"))

  loc_names <- c(
    bel = "Belize",
    boc = "Panama",
    hon = "Honduras",
    pri = "Puerto Rico",
    flk = "Florida Keys"
  )

  loc_title <- loc_names[unique(df$loc)[1]]

  # build strip data
  strip_df <- df %>%
    distinct(ind_label, spec)

  # plot
  ggplot(df, aes(x = ind_label, y = prob, fill = class)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = colors,
      breaks = c("P1", "P1_bc", "F1", "F2", "P2_bc", "P2"),
      name = "Ancestry",
      drop = FALSE
    ) +
    new_scale_fill() +
    geom_tile(
      data = strip_df,
      aes(x = ind_label, y = -0.05, fill = spec),
      height = 0.05,
      inherit.aes = FALSE
    ) +
    
    scale_fill_manual(values = color_map, labels = species_labels, name = "Species", drop = FALSE) +

    scale_x_discrete(
      labels = function(x) {
        ifelse(
          x %in% hybrids,
          paste0("<b>", x, "</b>"),
          x
        )
      }
    ) +

     # y axis fixed
    scale_y_continuous(
      breaks = c(0, 1),
      limits = c(-0.1, 1),
      expand = c(0, 0)
    ) +
    facet_grid(run_label ~ .) +
    # theme_minimal() +
    theme(
      legend.position = "left",
      strip.text.y = ggtext::element_markdown(angle = 0, size = 7),
      legend.text = element_markdown(size = 7),
      axis.text.x = ggtext::element_markdown(angle = 90, size = 7),
      axis.title.x = element_blank(),
      axis.title.y = element_text(vjust = 4),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(title = loc_title)
}

# ---------------------------------------------------------
# FINAL FIGURE
# ---------------------------------------------------------
combine_plots <- function(plot_list, output_path) {

  legend <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = "bottom"))

  p <- cowplot::plot_grid(
    cowplot::plot_grid(plotlist = plot_list, ncol = 1),
    legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
  )

  ggsave(
  filename = output_path,
  plot = p,
  height = 16,
  width = 10,
  dpi = 600,
  bg = "white"
  )
}

# -----------------------------
# Read pairwise LD file
# -----------------------------
read_ld_file <- function(file, label) {
  dt <- fread(file)
  
  r2 <- dt[[ncol(dt)]]
  
  data.frame(
    r2 = r2,
    dataset = label
  )
}

# -----------------------------
# Read global LD matrix (upper triangle)
# -----------------------------
read_global_ld <- function(file, label = "global") {
  mat <- as.matrix(fread(file))
  
  r2 <- mat[upper.tri(mat)]
  
  data.frame(
    r2 = r2,
    dataset = label
  )
}

# -----------------------------
# Make boxplot for one dataset
# -----------------------------
plot_ld_box <- function(df, title = NULL) {
  
  # define comparisons (all pairwise)
  comparisons <- list(
    c("global", "LG04_LG12_1"),
    c("global", "LG04_LG12_2"),
    c("global", "LG12_1_LG12_2")
  )
  
  ggplot(df, aes(x = dataset, y = r2, group = dataset)) +
    geom_boxplot(outlier.size = 0.3) +
    
    # stat_compare_means(
    #   comparisons = comparisons,
    #   method = "wilcox.test",
    #   label = "p.signif",
      
    #   # keep everything inside 0–0.1
    #   label.y = c(0.02, 0.22, 0.24),
    #   step.increase = 0
    # )  +
    
    coord_cartesian(ylim = c(0, 0.1)) +
    
    # scale_y_continuous(expand = expansion(mult = c(-0.15, 0.15))) +
    
    theme_classic() +    
    labs(
      x = NULL,
      y = expression(r^2),
      title = title
    )
}


# -----------------------------
# Build one dataset plot
# -----------------------------
build_plot <- function(ds) {

  message("Processing:    ", ds)
  
  files <- list(
    global = file.path(base_dir, "/2_popgen/ld/", paste0(ds, "_global.ld")),
    LG04_LG12_1 = file.path(base_dir, "/2_popgen/ld/", paste0(ds, ".LG04_LG12_1.ld")),
    LG04_LG12_2 = file.path(base_dir, "/2_popgen/ld/", paste0(ds, ".LG04_LG12_2.ld")),
    LG12_1_LG12_2 = file.path(base_dir, "/2_popgen/ld/", paste0(ds, ".LG12_1_LG12_2.ld"))
  )
  
  # use YOUR existing functions
  df_global <- read_global_ld(files$global, "global")
  df_1      <- read_ld_file(files$LG04_LG12_1, "LG04_LG12_1")
  df_2      <- read_ld_file(files$LG04_LG12_2, "LG04_LG12_2")
  df_3      <- read_ld_file(files$LG12_1_LG12_2, "LG12_1_LG12_2")
  
  df_all <- rbind(df_global, df_1, df_2, df_3)
  
  plot_ld_box(df_all, title = ds)
}