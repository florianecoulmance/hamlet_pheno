# by: Floriane Coulmance: 16/08/2024
# usage:
# Rscript Fig3_S21.R 
#___________________________________________________________________


# Clear the work space
rm(list = ls())


# new libraries
source("../scripts/R/helper_functions.R")
library(smartsnp)
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
library(vegan)
library(knitr)
library(plotly)
library(base64enc)
library(scatterplot3d)
library(nlme)
library(lme4)
library(lmerTest)
library(openxlsx)

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
base_path      <- "/Users/fcoulman/Desktop/hamlet_pheno/3_CHAPTER3/hamlet_pheno/"
figure_path     <- file.path(base_path, "figures")
path_phenotypes <- file.path(base_path, "1_phenotyping/pca")
path_genotypes_all <- file.path(base_path, "2_popgen/byALL")
path_genotypes_loc <- file.path(base_path, "2_popgen/byLOC")
fst_all_file <- file.path(figure_path, "TableS5.tex")
fst_loc_file <- file.path(figure_path, "TableS4.tex")
association_file <- file.path(base_path, "metadata/assortative_mating.csv")
transect_file <- file.path(base_path, "metadata/assortative_counts.csv")
logos_path      <- file.path(base_path, "metadata/logos_hamlet")
spec_colors     <- file.path(base_path, "metadata/species_colors.tsv")
geo_colors      <- file.path(base_path, "metadata/locations_colors.tsv")


if (!dir.exists(figure_path)) {
  dir.create(
    figure_path,
    recursive = TRUE
  )
}


# ============================================================
# Print summary for debugging
# ============================================================
cat("---- CONFIG ----\n")
cat("base_path:      ", base_path, "\n")
cat("figure_path:    ", figure_path, "\n")
cat("path_phenotypes:", path_phenotypes, "\n")
cat("path_genotypes_all:", path_genotypes_all, "\n")
cat("path_genotypes_loc:", path_genotypes_loc, "\n")
cat("association_file:", association_file, "\n")
cat("logos_path:     ", logos_path, "\n")
cat("spec_colors:    ", spec_colors, "\n")
cat("geo_colors:     ", geo_colors, "\n")
cat("-----------------\n")


# ============================================================
# GENERAL SETTINGS
# ============================================================
n_perm <- 10000
seed <- 123

species_info <- add_species_logos(spec_colors, logos_path)
# head(species_info)
geo_table <- read.delim(geo_colors, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# head(geo_table)


message("\n========================================")
message("PHENOTYPE ANALYSIS")
message("========================================")

pheno_files <- list.files(
  path = path_phenotypes,
  pattern = "_PCs\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(pheno_files) == 0) {
  stop(
    "No phenotype PCA files found in: ",
    path_phenotypes
  )
}

message(
  "Phenotype files found: ",
  length(pheno_files)
)

print(pheno_files)

pheno_info <- tibble(
  file = pheno_files,
  
  dataset = basename(file) %>%
    str_remove("_PCs\\.csv$")
) %>%
  mutate(
    location_code = str_extract(
      dataset,
      "(?<=lab_)[A-Za-z]+(?=[0-9_])"
    ),
    level = case_when(
      
      str_detect(
        dataset,
        "^lab_571_"
      ) ~ "all",
      
      TRUE ~ "location"
    )
  ) %>%
  filter(
    dataset %in% c(
      "lab_flo29_left_noflash",
      "lab_boc229_left_noflash",
      "lab_bel46_left_noflash",
      "lab_571_left_noflash"
    )
  )


print(pheno_info$level)

pheno_distances <- vector(
  "list",
  nrow(pheno_info)
)

pheno_distances_lda <- vector(
  "list",
  nrow(pheno_info)
)

for (i in seq_len(nrow(pheno_info))) {
  
  message(
    "\nPhenotype dataset ",
    i,
    "/",
    nrow(pheno_info),
    ": ",
    pheno_info$dataset[i]
  )
  
  pheno_distances[[i]] <-
    calculate_pheno_distance(
      pca_file = pheno_info$file[i],
      dataset = pheno_info$dataset[i],
      level = pheno_info$level[i],
      species_info = species_info,
      geo_table = geo_table
    )
  
  pheno_distances_lda[[i]] <-
    calculate_pheno_distance_lda(
      pca_file = pheno_info$file[i],
      dataset = pheno_info$dataset[i],
      level = pheno_info$level[i],
      species_info = species_info,
      geo_table = geo_table
    )
}

# pheno_distances <- bind_rows(
#   pheno_distances
# )

pheno_distances_lda <- bind_rows(
  pheno_distances_lda
)


message("\n========================================")
message("GENOTYPE ANALYSIS")
message("========================================")
geno_distances_fst <- read_fst_tables(
  s4_file = fst_loc_file,
  s5_file = fst_all_file,
  species_info = species_info,
  geo_table = geo_table
)


message("\n========================================")
message("ASSOCIATION / RI ANALYSIS")
message("========================================")

transect <- read.csv(
  transect_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Location = as.character(Location),
    species = as.character(species),
    count = as.numeric(count)
  )

association <- read.csv(
  association_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Location = as.character(Location),
    species1 = as.character(species1),
    species2 = as.character(species2),
    spawning = as.numeric(spawning)
  )

locations <- sort(
  unique(
    association$Location
  )
)


association_global <- association %>%
  group_by(
    species1,
    species2
  ) %>%
  summarise(
    count = sum(
      spawning,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

global_species_counts <- transect %>%
  group_by(species) %>%
  summarise(
    count = sum(count),
    .groups = "drop"
  ) %>%
  mutate(
    proportion = count / sum(count)
  )

RI_global <- calculate_global_RI(
  association_global,
  global_species_counts,
  species1_col = "species1",
  species2_col = "species2",
  spawning_count_col = "count",
  species_col = "species",
  species_count_col = "count")


asso_RI_global <- RI_global$results %>%
  select(
    species1,
    species2,
    observed,
    expected,
    RI
  ) %>%
  rename(
    distance_asso = RI
  ) %>%
  mutate(
    level = "all",
    Location = "all"
  )

asso_prop_global <- RI_global$species_proportions %>%
  mutate(
    level = "all",
    Location = "all"
  )
  
asso_RI_location <- vector(
  "list",
  length(locations)
)

asso_prop_location <- vector(
  "list",
  length(locations)
)

for (i in seq_along(locations)) {
  
  loc <- locations[i]
  
  message(
    "\nAssociation location ",
    i,
    "/",
    length(locations),
    ": ",
    loc
  )
  
  
  # --------------------------------------------------------
  # Extract this location
  # --------------------------------------------------------
  
  pairing_table <- association %>%
    filter(
      Location == loc
    ) %>%
    select(
      species1,
      species2,
      spawning
    ) %>%
    rename(
      count = spawning
    )
  
  print(pairing_table)
  
  species_counts <- transect %>%
    filter(
      Location == loc
    ) %>%
    mutate(
      proportion = count / sum(count)
    )

  print(species_counts)
  
  
  # --------------------------------------------------------
  # Run RI permutation
  # --------------------------------------------------------
  
  RI <- calculate_global_RI(
    pairing_table,
    species_counts,
    species1_col = "species1",
    species2_col = "species2",
    spawning_count_col = "count",
    species_col = "species",
    species_count_col = "count")

  print(RI$results)
  
  # --------------------------------------------------------
  # Standardise output
  # --------------------------------------------------------
  
  asso_RI_location[[i]] <-
    RI$results %>%
    select(
      species1,
      species2,
      observed,
      expected,
      RI
    ) %>%
    rename(
      distance_asso = RI
    ) %>%
    mutate(
      level = "location",
      Location = loc
    )
  
  asso_prop_location[[i]] <- 
    RI$species_proportions %>%
    mutate(
      level = "location",
      Location = loc
    )
}


asso_RI <- bind_rows(
  asso_RI_global,
  bind_rows(
    asso_RI_location
  ) 
) %>% mutate(
  Location = dplyr::recode(
    Location,
    "Bocas del Toro" = "Panama"
  )
)

asso_RI %>% filter(species1 == "aberrans", species2 == "puella")

# Within-species RI
RI_within <- asso_RI %>%
  group_by(level, Location) %>%
  filter(species1 == species2) %>%
  select(species = species1, RI_within = distance_asso)

# Between-species RI
RI_between <- asso_RI %>%
  group_by(level, Location) %>%
  filter(species1 != species2) %>%
  select(species1, species2, RI_between = distance_asso)

# Calculate mean RI for each inter-species pair
RI_inter_species <- RI_between %>%
  left_join(
    RI_within %>% rename(species1 = species, RI_species1 = RI_within),
    by = c("species1", "level", "Location")
  ) %>%
  left_join(
    RI_within %>% rename(species2 = species, RI_species2 = RI_within),
    by = c("species2", "level", "Location")
  ) %>%
  mutate(
    RI_mean = (RI_species1 + RI_species2 + RI_between) / 3
  )

RI_inter_species <- RI_inter_species %>%
  filter(!(Location %in% c("all", "Jamaica - Guna Yala"))) %>%
  rename(RI = RI_mean,
         component_species1 = RI_species1,
         component_species2 = RI_species2,
         component_inter = RI_between) %>%
  ungroup() %>%
  select(Location, species1, species2, RI, component_species1, component_species2, component_inter)

# write.xlsx(RI_inter_species, file.path(figure_path, "TableS8.xlsx"))
# writeLines(
#   knitr::kable(
#     RI_inter_species,
#     format = "latex",
#     booktabs = TRUE,
#     row.names = FALSE,
#     escape = TRUE
#   ),
#   con = file.path(
#     figure_path,
#     "TableS8.tex"
#   )
# )

asso_RI <- asso_RI %>% filter(!(Location %in% c("all", "Jamaica - Guna Yala"))) %>%
  rename(observed_pairing = observed,
         expected_pairing = expected) %>%
  ungroup() %>%
  select(Location, species1, species2, expected_pairing, observed_pairing)

# write.xlsx(asso_RI, file.path(figure_path, "TableS9.xlsx"))
# writeLines(
#   knitr::kable(
#     asso_RI,
#     format = "latex",
#     booktabs = TRUE,
#     row.names = FALSE,
#     escape = TRUE
#   ),
#   con = file.path(
#     figure_path,
#     "TableS9.tex"
#   )
# )


RI_inter_species <- RI_inter_species %>%
  mutate(level = "location") %>%
  select(level, Location, species1, species2, RI, component_species1, component_species2, component_inter)

species_pairs <- tibble::tribble(
  ~species1,    ~species2,
  "aberrans",   "indigo",
  "indigo",     "puella",
  "indigo",     "unicolor",
  "nigricans",  "unicolor",
  "aberrans",   "unicolor",
  "indigo",     "nigricans",
  "nigricans",  "puella",
  "puella",     "unicolor",
  "aberrans",   "nigricans",
  "aberrans",   "puella"
)

RI_all <- RI_inter_species %>%
  filter(level == "location") %>%
  filter(is.finite(RI)) %>%
  inner_join(species_pairs, by = c("species1", "species2")) %>%
  group_by(species1, species2) %>%
  summarise(
    RI = mean(RI),
    .groups = "drop"
  ) %>%
  mutate(
    level = "all",
    Location = "all"
  ) %>%
  select(level, Location, species1, species2, RI)

RI <- bind_rows(
  RI_inter_species,
  bind_rows(
    RI_all
  ) 
)

message("\n========================================")
message("INTERSECTIONS")
message("========================================")
speciation_hypercube_data <- pheno_distances_lda %>%
  dplyr::inner_join(
    geno_distances_fst,
    by = c(
      "level",
      "Location",
      "species1",
      "species2"
    )
  ) %>%
  dplyr::inner_join(
    RI,
    by = c(
      "level",
      "Location",
      "species1",
      "species2"
    )
  ) %>%
  mutate(distance_asso = RI)
print(speciation_hypercube_data)

pheno_geno_data_all <- pheno_distances_lda %>%
  dplyr::inner_join(
    geno_distances_fst,
    by = c(
      "level",
      "Location",
      "species1",
      "species2"
    )
  ) %>%
  filter(level == "all")

pheno_geno_data_location <- pheno_distances_lda %>%
  dplyr::inner_join(
    geno_distances_fst,
    by = c(
      "level",
      "Location",
      "species1",
      "species2"
    )
  ) %>%
  filter(level == "location")

pheno_asso_data_all <- pheno_distances_lda %>%
  dplyr::inner_join(
    RI,
    by = c(
      "level",
      "Location",
      "species1",
      "species2"
    )
  ) %>%
  filter(level == "all") %>%
  mutate(distance_asso = RI)

pheno_asso_data_location <- pheno_distances_lda %>%
  dplyr::inner_join(
    RI,
    by = c(
      "level",
      "Location",
      "species1",
      "species2"
    )
  ) %>%
  filter(level == "location") %>%
  mutate(distance_asso = RI)

geno_asso_data_all <- geno_distances_fst %>%
  dplyr::inner_join(
    RI,
    by = c(
      "level",
      "Location",
      "species1",
      "species2"
    )
  ) %>%
  filter(level == "all") %>%
  mutate(distance_asso = RI)

geno_asso_data_location <- geno_distances_fst %>%
  dplyr::inner_join(
    RI,
    by = c(
      "level",
      "Location",
      "species1",
      "species2"
    )
  ) %>%
  filter(level == "location",
         RI != "-Inf") %>%
  mutate(distance_asso = RI)


message("\n========================================")
message("HYPERCUBE")
message("========================================")
# hypercube <- plot_speciation_hypercube(
#   speciation_hypercube_data %>% filter(level == "location")
# )
# print(hypercube)

data_location <- speciation_hypercube_data %>%
  filter(level == "location")

data_all <- speciation_hypercube_data %>%
  filter(level == "all")


scatter_location <- scatterplot3d(
  data_location[, c(5, 6, 7)],
  pch = 16,
  color = "#D06495",
  angle = 40,
  cex.symbols = 2,
  type = "h",
  aspect = c(1, 1, 1),
  xlim = c(0, 1),
  ylim = c(0, 1),
  zlim = c(0, 1),
  xlab = "Phenotypic divergence",
  ylab = "",
  zlab = "RI",
  lab = c(3, 3, 3),
  x.ticklabs = c(0, 0.5, 1),
  y.ticklabs = c(0, 0.5, 1),
  z.ticklabs = c(0, 0.5, 1),
  font.lab = 2
)
# title(main = "Between sympatric species", font.main = 2, cex.main = 1.3)
grid.echo()

grid.gedit("graphics-plot-1-lines-4",
           gp = gpar(lwd = 0))

grid.gedit("graphics-plot-1-lines-6",
           gp = gpar(lwd = 0))

grid.gedit("graphics-plot-1-lines-10",
           gp = gpar(lwd = 0.2, lty = 2))
grid.gedit("graphics-plot-1-lines-11",
           gp = gpar(lwd = 0.2, lty = 2))
grid.gedit("graphics-plot-1-lines-12",
           gp = gpar(lwd = 0.2, lty = 2))
grid.gedit(
  "graphics-plot-1-lines-8",
  gp = gpar(lwd = 0.2, lty = 2)
)

grid.gedit("graphics-plot-1-lines-2", gp = gpar(lwd = 3))
grid.gedit("graphics-plot-1-lines-3", gp = gpar(lwd = 3))
g <- grid.get("graphics-plot-1-lines-1")
g$gp$lwd[1] <- 3
grid.set("graphics-plot-1-lines-1", g)

g_p <- grid.get("graphics-plot-1-points-1")
g_p$gp$fontsize <- 18
g_p$gp$alpha[1] <- 0.75
grid.set("graphics-plot-1-points-1", g_p)

grid.text(
  "Genetic differentiation",
  x = unit(0.76, "npc"),
  y = unit(0.23, "npc"),
  rot = 33,
  gp = gpar(fontface = "bold")
)
scatter_location <- grid.grab()
png(
  file.path(figure_path, "scatter_location.png"),
  width = 2200,
  height = 2000,
  res = 300
)
grid.draw(scatter_location)
dev.off()
dev.off()

scatter_all <- scatterplot3d(
  data_all[, c(5, 6, 7)],
  pch = 15,
  color = "#D09F64",
  angle = 40,
  cex.symbols = 2,
  type = "h",
  aspect = c(1, 1, 1),
  xlim = c(0, 1),
  ylim = c(0, 1),
  zlim = c(0, 1),
  xlab = "Phenotypic divergence",
  ylab = "",
  zlab = "RI",
  lab = c(3, 3, 3),
  x.ticklabs = c(0, 0.5, 1),
  y.ticklabs = c(0, 0.5, 1),
  z.ticklabs = c(0, 0.5, 1),
  font.lab = 2
)
# title(main = "Between species", font.main = 2, cex.main = 1.3)
grid.echo()

grid.gedit("graphics-plot-1-lines-4",
           gp = gpar(lwd = 0))

grid.gedit("graphics-plot-1-lines-6",
           gp = gpar(lwd = 0))

grid.gedit("graphics-plot-1-lines-10",
           gp = gpar(lwd = 0.2, lty = 2))
grid.gedit("graphics-plot-1-lines-11",
           gp = gpar(lwd = 0.2, lty = 2))
grid.gedit("graphics-plot-1-lines-12",
           gp = gpar(lwd = 0.2, lty = 2))
grid.gedit(
  "graphics-plot-1-lines-8",
  gp = gpar(lwd = 0.2, lty = 2)
)

grid.gedit("graphics-plot-1-lines-2", gp = gpar(lwd = 3))
grid.gedit("graphics-plot-1-lines-3", gp = gpar(lwd = 3))
g <- grid.get("graphics-plot-1-lines-1")
g$gp$lwd[1] <- 3
grid.set("graphics-plot-1-lines-1", g)

g_p <- grid.get("graphics-plot-1-points-1")
g_p$gp$fontsize <- 18
g_p$gp$alpha[1] <- 0.75
grid.set("graphics-plot-1-points-1", g_p)
grid.text(
  "Genetic differentiation",
  x = unit(0.76, "npc"),
  y = unit(0.23, "npc"),
  rot = 33,
  gp = gpar(fontface = "bold")
)
scatter_all <- grid.grab()
png(
  file.path(figure_path, "scatter_all.png"),
  width = 2200,
  height = 2000,
  res = 300
)
grid.draw(scatter_all)
dev.off()
dev.off()


p_pheno_geno_location <- plot_speciation_paper(
  data = speciation_hypercube_data,
  species_meta = species_info,
  x = "distance_pheno",
  y = "distance_geno",
  colour = "distance_asso",
  x_label = "Phenotypic divergence",
  y_label = expression(Genetic~differentiation~(F[ST])),
  colour_label = expression(RI),
  label_species_pair = TRUE,
  label_location = TRUE,
  legend_bottom_right = FALSE,
  panel = "location"
)
p_pheno_geno_location

# p_pheno_asso_location <- plot_speciation_paper(
#   speciation_hypercube_data,
#   species_meta = species_info,
#   x = "distance_pheno",
#   y = "distance_asso",
#   colour = NULL,
#   x_label = "Phenotypic divergence",
#   y_label = expression(RI),
#   colour_label = NULL,
#   label_species_pair = FALSE,
#   label_location = FALSE,
#   legend_bottom_right = FALSE,
#   panel = "location"
# )
# p_pheno_asso_location
# 
# p_geno_asso_location <- plot_speciation_paper(
#   speciation_hypercube_data,
#   species_meta = species_info,
#   x = "distance_geno",
#   y = "distance_asso",
#   colour = NULL,
#   x_label = expression(Genetic~divergence~(F[ST])),
#   y_label = expression(RI),
#   colour_label = NULL,
#   label_species_pair = FALSE,
#   label_location = FALSE,
#   legend_bottom_right = TRUE,
#   panel = "location"
# )
# p_geno_asso_location



p_pheno_geno_all <- plot_speciation_paper(
  speciation_hypercube_data,
  species_meta = species_info,
  x = "distance_pheno",
  y = "distance_geno",
  colour = "distance_asso",
  x_label = "Phenotypic divergence",
  y_label = expression(Genetic~differentiation~(F[ST])),
  colour_label = expression(RI),
  label_species_pair = TRUE,
  label_location = FALSE,
  legend_bottom_right = FALSE,
  panel = "all"
)
p_pheno_geno_all

# p_pheno_asso_all <- plot_speciation_paper(
#   speciation_hypercube_data,
#   species_meta = species_info,
#   x = "distance_pheno",
#   y = "distance_asso",
#   colour = NULL,
#   x_label = "Phenotypic divergence",
#   y_label = expression(RI),
#   colour_label = NULL,
#   label_species_pair = FALSE,
#   label_location = FALSE,
#   legend_bottom_right = FALSE,
#   panel = "all"
# )
# p_pheno_asso_all
# 
# p_geno_asso_all <- plot_speciation_paper(
#   speciation_hypercube_data,
#   species_meta = species_info,
#   x = "distance_geno",
#   y = "distance_asso",
#   colour = NULL,
#   x_label = expression(Genetic~divergence~(F[ST])),
#   y_label = expression(RI),
#   colour_label = NULL,
#   label_species_pair = FALSE,
#   label_location = FALSE,
#   legend_bottom_right = TRUE,
#   panel = "all"
# )
# p_geno_asso_all


# ############################
# FINAL PLOTS
# ############################
########## FIGURE 3 #################### 
# figure3_bottom <- plot_grid(
#   p_pheno_asso_location,
#   p_geno_asso_location,
#   ncol = 2,
#   labels = c("(b)", "(c)"),
#   scale = 0.95
# )

figure3 <- ggarrange(
  p_pheno_geno_location,
  NULL,
  p_pheno_geno_all,
  nrow = 3,
  heights = c(10, 0.2, 10),
  labels = c("(a)", "", "(b)")
)

ggsave(
  filename = file.path(figure_path, "Fig3_cube.png"),
  plot = figure3,
  width = 16,
  height = 23,
  units = "in",
  dpi = 300,
  type = "cairo-png"
)

# htmlwidgets::saveWidget(
#   hypercube,
#   file.path(figure_path, "Fig3_speciation3D.html"),
#   selfcontained = FALSE
# )

# ########## FIGURE S21 ###################
# figureS21_bottom <- plot_grid(
#   p_pheno_asso_all,
#   p_geno_asso_all,
#   ncol = 2,
#   labels = c("(b)", "(c)"),
#   scale = 0.95
# )
# 
# figureS21 <- ggarrange(
#   p_pheno_geno_all,
#   NULL,
#   figureS21_bottom,
#   nrow = 3,
#   heights = c(10, 0.2, 5),
#   labels = c("(a)", "", "")
# )
# 
# ggsave(
#   filename = file.path(figure_path, "FigS21_speciationALL.png"),
#   plot = figureS21,
#   width = 17,
#   height = 18,
#   units = "in",
#   dpi = 300,
#   type = "cairo-png"
# )


########## FIGURE S21 ###################
p_pheno_geno_all_2 <- plot_speciation_paper(
  pheno_geno_data_all,
  species_meta = species_info,
  x = "distance_pheno",
  y = "distance_geno",
  colour = NULL,
  x_label = "Phenotypic divergence",
  y_label = expression(Genetic~differentiation~(F[ST])),
  colour_label = NULL,
  label_species_pair = FALSE,
  label_location = FALSE,
  legend_bottom_right = FALSE,
  panel = "all"
)
p_pheno_geno_all_2

p_pheno_asso_all_2 <- plot_speciation_paper(
  pheno_asso_data_all,
  species_meta = species_info,
  x = "distance_pheno",
  y = "distance_asso",
  colour = NULL,
  x_label = "Phenotypic divergence",
  y_label = expression(RI),
  colour_label = NULL,
  label_species_pair = FALSE,
  label_location = FALSE,
  legend_bottom_right = FALSE,
  panel = "all"
)
p_pheno_asso_all_2

p_geno_asso_all_2 <- plot_speciation_paper(
  geno_asso_data_all,
  species_meta = species_info,
  x = "distance_geno",
  y = "distance_asso",
  colour = NULL,
  x_label = expression(Genetic~differentiation~(F[ST])),
  y_label = expression(RI),
  colour_label = NULL,
  label_species_pair = FALSE,
  label_location = FALSE,
  legend_bottom_right = TRUE,
  panel = "all"
)
p_geno_asso_all_2

all <- plot_grid(
  p_pheno_geno_all_2,
  p_pheno_asso_all_2,
  p_geno_asso_all_2,
  nrow = 3,
  labels = c("(a)", "(b)", "(c)")
)
all

p_pheno_geno_loc_2 <- plot_speciation_paper(
  pheno_geno_data_location,
  species_meta = species_info,
  x = "distance_pheno",
  y = "distance_geno",
  colour = NULL,
  x_label = "Phenotypic divergence",
  y_label = expression(Genetic~differentiation~(F[ST])),
  colour_label = NULL,
  label_species_pair = FALSE,
  label_location = FALSE,
  legend_bottom_right = FALSE,
  panel = "location"
)
p_pheno_geno_loc_2

p_pheno_asso_loc_2 <- plot_speciation_paper(
  pheno_asso_data_location,
  species_meta = species_info,
  x = "distance_pheno",
  y = "distance_asso",
  colour = NULL,
  x_label = "Phenotypic divergence",
  y_label = expression(RI),
  colour_label = NULL,
  label_species_pair = FALSE,
  label_location = FALSE,
  legend_bottom_right = FALSE,
  panel = "location"
)
p_pheno_asso_loc_2

p_geno_asso_loc_2 <- plot_speciation_paper(
  geno_asso_data_location,
  species_meta = species_info,
  x = "distance_geno",
  y = "distance_asso",
  colour = NULL,
  x_label = expression(Genetic~differentiation~(F[ST])),
  y_label = expression(RI),
  colour_label = NULL,
  label_species_pair = FALSE,
  label_location = FALSE,
  legend_bottom_right = TRUE,
  panel = "location"
)
p_geno_asso_loc_2

loc <- plot_grid(
  p_pheno_geno_loc_2,
  p_pheno_asso_loc_2,
  p_geno_asso_loc_2,
  nrow = 3
)
loc

figureS21 <- ggarrange(
  all,
  loc,
  ncol = 2,
  widths = c(1, 1)
)

ggsave(
  filename = file.path(figure_path, "FigS21_correlations.png"),
  plot = figureS21,
  width = 17,
  height = 22,
  units = "in",
  dpi = 300,
  type = "cairo-png"
)

