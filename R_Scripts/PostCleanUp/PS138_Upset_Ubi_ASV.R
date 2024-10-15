
#############################################
## ASV pres-abs 
#############################################
setwd("~/PS138")
load("~/PS138/ampvis_subsets.Rdata")


ASV.hel <- ampvis_hel_nofobs$abund
ENV.hel <-ampvis_hel_nofobs$metadata
TAX.hel <- ampvis_hel_nofobs$tax
library(ampvis2)
library("ape")
library("phyloseq")
library(ggplot2)
library(dplyr)
library(DECIPHER)
library(phangorn)
library(dada2)
library(ggtree)
library(phyloseq)
library(dendextend)
library(phyloseq)
library(SRS)
library(tidyr)
library(microbiome)
library(scico)
library(UpSetR)
library(tibble)
ice_depths <- c( "MP", "Ice-TS","Ice-BS","UIW" )
euphotic_zone_depths <- c( "2", "10","25", "chl-max", "50", "100", "200")
pelagic_water_depths <- c("500", "1000", "1500", "2000", "3000")
bottom_water_depths <- c("bottom", "20m from bottom", "5m from bottom")


core1 <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name = "Abundance") %>%
  dplyr::rename(sample_title = variable) %>%
  left_join(ENV.hel, by = c("sample_title" = "CTD #")) %>%
  filter(depth_category %in% ice_depths) %>%
  group_by(asv) %>% 
  summarize(Abundance = mean(Abundance)) %>% 
  filter(Abundance >= 0.005) %>%
  distinct(asv, .keep_all = TRUE) %>%
  column_to_rownames("asv")

core2 <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name = "Abundance") %>%
  dplyr::rename(sample_title = variable) %>%
  left_join(ENV.hel, by = c("sample_title" = "CTD #")) %>%
  filter(depth_category %in% euphotic_zone_depths) %>%
  group_by(asv) %>% 
  summarize(Abundance = mean(Abundance)) %>% 
  filter(Abundance >= 0.005) %>%
  distinct(asv, .keep_all = TRUE) %>%
  column_to_rownames("asv")

core3 <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name = "Abundance") %>%
  dplyr::rename(sample_title = variable) %>%
  left_join(ENV.hel, by = c("sample_title" = "CTD #")) %>%
  filter(depth_category %in% pelagic_water_depths) %>%
  group_by(asv) %>% 
  summarize(Abundance = mean(Abundance)) %>% 
  filter(Abundance >= 0.005) %>%
  distinct(asv, .keep_all = TRUE) %>%
  column_to_rownames("asv")

core4 <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name = "Abundance") %>%
  dplyr::rename(sample_title = variable) %>%
  left_join(ENV.hel, by = c("sample_title" = "CTD #")) %>%
  filter(depth_category %in% bottom_water_depths) %>%
  group_by(asv) %>% 
  summarize(Abundance = mean(Abundance)) %>% 
  filter(Abundance >= 0.005) %>%
  distinct(asv, .keep_all = TRUE) %>%
  column_to_rownames("asv")



# Combine in list
core <- list()
core[["ice"]] <- as.character(row.names(core1))
core[["euphotic zone"]] <- as.character(row.names(core2))
core[["pelagic water"]] <- as.character(row.names(core3))
core[["bottom water"]] <- as.character(row.names(core4))

# Function for presence-absence matrix
presAbs <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

overlaps <- presAbs(core) %>%
  rownames_to_column("asv")

set_order <- c("bottom water", "pelagic water", "euphotic zone", "ice")
core_colors <- c("bottom water" = "brown4", 
                 "pelagic water" = "navyblue", 
                 "euphotic zone" = "lawngreen", 
                 "ice" = "turquoise")

upset(
  overlaps,
  sets = set_order,
  keep.order = TRUE,
  number.angles = 0, 
  nsets=4,
  main.bar.color = "gray22",
  sets.bar.color = core_colors,
  matrix.color = "gray22",
  set_size.show = TRUE,
  point.size = 2.44, line.size = 0.8, text.scale = 1.2,
  mainbar.y.label = "Shared ASVs",
  order.by = c("freq"))
#decreasing = c(FALSE))

all_depths_asvs <- overlaps %>%
  filter(ice == 1 & `euphotic zone` == 1 & `pelagic water` == 1 & `bottom water` == 1) %>%
  pull(asv)

print(paste("Number of ASVs present in all depth classes:", length(all_depths_asvs)))
print("ASVs present in all depth classes:")
print(all_depths_asvs)


# Select the ASVs you want to plot
selected_asvs <- all_depths_asvs

# Extract abundance data
abundance_data <- ampvis_hel_nofobs$abund %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ASV") %>%
  filter(ASV %in% selected_asvs) %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "Abundance")
library(dplyr)
# Extract metadata
metadata <- ampvis_hel_nofobs$metadata %>%
  dplyr::select(sample_id = `CTD #`, depth_category)

# Prepare the data
plot_data <- abundance_data %>%
  left_join(metadata, by = "sample_id") %>%
  mutate(depth_category = factor(depth_category, levels = c(
    "MP", "Ice-TS", "Ice-BS", "UIW", "2", "10", "chl-max", "25", "50", 
    "100", "200", "500", "1000", "1500", "2000", "3000", 
    "20m from bottom", "5m from bottom", "bottom"
  )))

# Create the plot
p <- ggplot(plot_data, aes(x = depth_category, y = Abundance, fill = ASV)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = 21, outlier.size = 1.5) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.major.x = element_blank()) +
  labs(x = "Depth Category", 
       y = "Abundance (log scale)", 
       title = "Selected ASV Abundances Across Depth Categories",
       fill = "ASV") +
  scale_fill_brewer(palette = "Set1")

# Print the plot
print(p)

plot_data_mean <- plot_data %>%
  group_by(ASV, depth_category) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# Create the line plot
p2 <- ggplot(plot_data_mean, aes(x = depth_category, y = mean_abundance, color = ASV, group = ASV)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.major.x = element_blank()) +
  labs(x = "Depth Category", 
       y = "Mean Abundance (log scale)", 
       title = "Selected ASV Abundances Across Depth Categories",
       color = "ASV") +
  scale_color_brewer(palette = "Set1")

# Print the plot
print(p2)

# Get taxonomic information
tax_info <- ampvis_hel_nofobs$tax %>%
  filter(OTU %in% all_depths_asvs) %>%
  dplyr::select(ASV = OTU, Phylum, Genus)

# Extract abundance data
abundance_data <- ampvis_hel_nofobs$abund %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ASV") %>%
  filter(ASV %in% all_depths_asvs) %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "Abundance")

# Extract metadata
metadata <- ampvis_hel_nofobs$metadata %>%
  dplyr::select(sample_id = `CTD #`, depth_category)

metadata$depth_category[metadata$depth_category == "25"] <- "chl-max"
# Prepare the data
plot_data <- abundance_data %>%
  left_join(metadata, by = "sample_id") %>%
  left_join(tax_info, by = "ASV") %>%
  mutate(depth_category = factor(depth_category, levels = rev(c(
    "MP", "Ice-TS", "Ice-BS", "UIW", "2", "10", "chl-max", "50", 
    "100", "200", "500", "1000", "1500", "2000", "3000", 
    "20m from bottom", "5m from bottom", "bottom"
  ))))

# Calculate mean abundance for each ASV in each depth category
plot_data_mean <- plot_data %>%
  group_by(ASV, Phylum, Genus, depth_category) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# Define a color palette with enough colors to repeat
color_palette <- c(
  "red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "yellow",
  "darkred", "darkblue", "darkgreen", "purple4", "darkorange", "brown4", "pink4", "darkcyan", "darkmagenta", "yellow4"
)

# Create the faceted line plot with repeating colors
p3 <- ggplot(plot_data_mean, aes(x = depth_category, y = mean_abundance, color = ASV, group = ASV)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~ Phylum, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "lightgrey")) +
  labs(x = "Depth Category", 
       y = "Mean Abundance (log scale)", 
       title = "ASV Abundances Across Depth Categories by Phylum",
       color = "ASV") +
  scale_color_manual(values = rep(color_palette, length.out = length(unique(plot_data_mean$ASV))))+
  coord_flip()

# Print the plot
print(p3)
#####################

# First, filter out Bacteroidota and Proteobacteria
plot_data_filtered <- plot_data_mean %>%
  filter(!(Phylum %in% c("Bacteroidota", "Proteobacteria")))

# Create a color palette
n_colors <- length(unique(plot_data_filtered$ASV))
color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(n_colors)

# Create the plot
p3.5 <- ggplot(plot_data_filtered, aes(x = depth_category, y = mean_abundance, color = Genus, group = ASV)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~ Phylum, scales = "fixed",ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "lightgrey")) +
  labs(x = "Depth Category", 
       y = "Mean Abundance (log scale)", 
       title = "ASV Abundances Across Depth Categories by Phylum (excluding Bacteroidota and Proteobacteria)",
       color = "ASV") +
  scale_color_manual(values = color_palette) +
  coord_flip()

# Print the plot
print(p3.5)

# Save the plot
ggsave("ASV_abundances_by_phylum_excluding_Bacteroidota_Proteobacteria.png", p3, width = 16, height = 12, dpi = 300)
# If there are too many ASVs for a single legend, you can use the Genus instead
p4 <- ggplot(plot_data_mean, aes(x = depth_category, y = mean_abundance, color = Genus, group = ASV)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~ Phylum, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "lightgrey")) +
  labs(x = "Depth Category", 
       y = "Mean Abundance (log scale)", 
       title = "ASV Abundances Across Depth Categories by Phylum",
       color = "Genus") +
  scale_color_manual(values = rep(color_palette, length.out = length(unique(plot_data_mean$Genus))))

# Print the plot
print(p4)
################
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# Filter the taxonomic information for the selected phyla and ASVs
selected_phyla <- c("Proteobacteria", "Bacteroidota")
tax_info_filtered <- ampvis_hel_nofobs$tax %>%
  filter(OTU %in% all_depths_asvs, Phylum %in% selected_phyla) %>%
  dplyr::select(ASV = OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species)

# Extract abundance data for the selected ASVs
abundance_data_filtered <- ampvis_hel_nofobs$abund %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ASV") %>%
  filter(ASV %in% all_depths_asvs) %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "Abundance")



# Prepare the data for plotting
plot_data <- abundance_data_filtered %>%
  left_join(metadata, by = "sample_id") %>%
  left_join(tax_info_filtered, by = "ASV") %>%
  mutate(depth_category = factor(depth_category, levels = rev(c(
    "MP", "Ice-TS", "Ice-BS", "UIW", "2", "10", "chl-max", "50", 
    "100", "200", "500", "1000", "1500", "2000", "3000", 
    "20m from bottom", "5m from bottom", "bottom"
  ))))

# Calculate mean abundance for each ASV in each depth category
plot_data_mean <- plot_data %>%
  group_by(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, depth_category) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

#new
# Filter the data for Proteobacteria
proteobacteria_data <- plot_data_mean %>%
  filter(Phylum == "Proteobacteria")

# Similarly, you might want to create bacteroidota_data
bacteroidota_data <- plot_data_mean %>%
  filter(Phylum == "Bacteroidota")
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggnewscale)


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

create_phylum_plot <- function(data, phylum_name) {
  # Ensure the data is for the specified phylum
  if (unique(data$Phylum) != phylum_name) {
    stop("The provided data does not match the specified phylum.")
  }
  
  # Get unique classes and genera
  classes <- unique(data$Class)
  genera_by_class <- split(unique(data$Genus), data$Class)
  
  # Create color palettes for each class
  color_palettes <- list(
    brewer.pal(9, "Set1"),
    brewer.pal(8, "Set2"),
    brewer.pal(12, "Set3"),
    brewer.pal(8, "Dark2"),
    brewer.pal(8, "Paired")
  )
  
  # Assign colors to genera within each class
  color_scale <- lapply(seq_along(classes), function(i) {
    class_genera <- genera_by_class[[i]]
    n_genera <- length(class_genera)
    palette <- colorRampPalette(color_palettes[[i]])(n_genera)
    setNames(palette, class_genera)
  })
  names(color_scale) <- classes
  
  # Combine all color scales
  all_colors <- unlist(color_scale)
  
  # Create the plot
  p <- ggplot(data, aes(x = depth_category, y = mean_abundance, group = ASV)) +
    geom_line(aes(color = Genus), linewidth = 1) +
    geom_point(aes(color = Genus), size = 2) +
    scale_y_log10() +
    facet_wrap(Class ~ Order, scales = "free_y", nrow = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "right",
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "lightgrey"),
          strip.placement = "outside") +
    labs(x = "Depth Category", 
         y = "Mean Abundance (log scale)", 
         title = paste(phylum_name, "ASV Abundances Across Depth Categories"),
         color = "Genus") +
    scale_color_manual(values = all_colors) +
    coord_flip()
  
  # Create separate legends for each class
  legend_list <- lapply(classes, function(class) {
    guide_legend(title = paste("Genus (", class, ")"),
                 override.aes = list(color = color_scale[[class]]),
                 ncol = 1,
                 byrow = TRUE)
  })
  names(legend_list) <- paste0("color_", seq_along(classes))
  
  # Apply the split legends
  p <- p + guides(color = legend_list)
  
  return(p)
}
p_proteobacteria <- create_phylum_plot(proteobacteria_data, "Proteobacteria")
print(p_proteobacteria)
ggsave("heatmaps/proteobacteria_abundance_by_class.png", p_proteobacteria, width = 20, height = 12, dpi = 300)

# Create and save the plot for Bacteroidota
p_bacteroidota <- create_phylum_plot(bacteroidota_data, "Bacteroidota")
print(p_bacteroidota)
ggsave("heatmaps/bacteroidota_abundance_by_class.png", p_bacteroidota, width = 20, height= 12, dpi = 300)



##########
#####

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(grid)

create_proteobacteria_plot <- function(data, class_name) {
  # Filter data for the specified class
  class_data <- data %>% filter(Class == class_name)
  
  # Create a color palette for genera
  n_genera <- length(unique(class_data$Genus))
  color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(n_genera)
  
  # Assign colors to genera
  color_scale <- setNames(color_palette, unique(class_data$Genus))
  
  # Ensure depth_category is a factor with correct levels
  depth_order <- rev(c("MP", "Ice-TS", "Ice-BS", "UIW", "2", "10", "25", "chl-max", "50", "100", "200",
                       "500", "1000", "1500", "2000", "3000", "20m from bottom", "5m from bottom", "bottom"))
  class_data$depth_category <- factor(class_data$depth_category, levels = depth_order)
  
  
  # Create the plot
  p <- ggplot(class_data, aes(x = depth_category, y = mean_abundance, group = ASV)) +
    geom_line(aes(color = Genus), linewidth = 1) +
    geom_point(aes(color = Genus), size = 2) +
    scale_y_log10() +
    facet_wrap(~ Order, scales = "fixed", ncol = 5) +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "lightgrey"),
      strip.placement = "outside",
    ) +
    labs(x = "Depth Category", 
         y = "Mean Abundance ", 
         title = paste(class_name, "ASV Abundances Across Depth Categories"),
         color = "Genus") +
    scale_color_manual(values = color_scale) +
    coord_flip()
  
  
  
  
  return(p)
}

# Create plots for Alpha and Gammaproteobacteria
p_alpha <- create_proteobacteria_plot(proteobacteria_data, "Alphaproteobacteria")
p_gamma <- create_proteobacteria_plot(proteobacteria_data, "Gammaproteobacteria")

# Combine the plots
combined_plot <- grid.arrange(p_alpha, p_gamma, ncol = 1,
                              top = textGrob("Proteobacteria ASV Abundances Across Depth Categories",
                                             gp = gpar(fontsize = 16, font = 2)))

# Save the combined plot
ggsave("heatmaps/proteobacteria_abundance_by_class.png", combined_plot, width = 24, height = 16, dpi = 300)
# Save the combined plot
print(combined_plot)

# Save the combined plot
##################
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(grid)

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

create_bacteroidota_plot <- function(data, class_name) {
  # Filter data for the specified class and exclude Flavobacteriaceae
  class_data <- data %>% 
    filter(Class == class_name, Family != "Flavobacteriaceae")
  
  # Create a color palette for genera
  n_genera <- length(unique(class_data$Genus))
  color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(n_genera)
  
  # Assign colors to genera
  color_scale <- setNames(color_palette, unique(class_data$Genus))
  
  # Ensure depth_category is a factor with correct levels
  depth_order <- rev(c("MP", "Ice-TS", "Ice-BS", "UIW", "2", "10", "25", "chl-max", "50", "100", "200",
                       "500", "1000", "1500", "2000", "3000", "20m from bottom", "5m from bottom", "bottom"))
  class_data$depth_category <- factor(class_data$depth_category, levels = depth_order)
  
  # Create the plot
  p <- ggplot(class_data, aes(x = depth_category, y = mean_abundance, group = ASV)) +
    geom_line(aes(color = Genus), linewidth = 1) +
    geom_point(aes(color = Genus), size = 2) +
    scale_y_log10() +
    facet_wrap(~ Family, scales = "fixed", ncol = 3) +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "lightgrey"),
      strip.placement = "outside",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.ticks.y = element_blank(),
      plot.margin = unit(c(1, 1, 1, 4), "lines")
    ) +
    labs(x = "Depth Category", 
         y = "Mean Abundance", 
         title = paste(class_name, "ASV Abundances Across Depth Categories (excluding Flavobacteriaceae)"),
         color = "Genus") +
    scale_color_manual(values = color_scale) +
    coord_flip()
  
  return(p)
}

# Create plot for Bacteroidia class, excluding Flavobacteriaceae
p_bacteroidia <- create_bacteroidota_plot(bacteroidota_data, "Bacteroidia")

# Print the plot
print(p_bacteroidia)

# Save the plot
ggsave("heatmaps/bacteroidia_abundance_by_family_excluding_flavobacteriaceae.png", p_bacteroidia, width = 16, height = 12, dpi = 300)
# Assuming you have a bacteroidota_data dataframe similar to proteobacteria_data
# Create plots for Bacteroidota classes
p_bacteroidia <- create_bacteroidota_plot(bacteroidota_data, "Bacteroidia")
print(p_bacteroidia)
p_flavobacteriia <- create_bacteroidota_plot(bacteroidota_data, "Flavobacteriia")
print(p_flavobacteriia)
# Combine the plots
combined_plot_bacteroidota <- grid.arrange(p_bacteroidia, p_flavobacteriia, ncol = 1,
                                           top = textGrob("Bacteroidota ASV Abundances Across Depth Categories",
                                                          gp = gpar(fontsize = 16, font = 2)))

# Save the combined plot
ggsave("heatmaps/bacteroidota_abundance_by_class.png", combined_plot_bacteroidota, width = 24, height = 16, dpi = 300)

# Print the combined plot
print(combined_plot_bacteroidota)

##########
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

create_flavobacteriaceae_plot <- function(data) {
  # Filter data for Flavobacteriaceae family
  family_data <- data %>% 
    filter(Family == "Flavobacteriaceae")
  
  # Create a color palette for ASVs
  n_asvs <- length(unique(family_data$ASV))
  color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(n_asvs)
  
  # Ensure depth_category is a factor with correct levels
  depth_order <- rev(c("MP", "Ice-TS", "Ice-BS", "UIW", "2", "10", "25", "chl-max", "50", "100", "200",
                       "500", "1000", "1500", "2000", "3000", "20m from bottom", "5m from bottom", "bottom"))
  family_data$depth_category <- factor(family_data$depth_category, levels = depth_order)
  
  # Create the plot
  p <- ggplot(family_data, aes(x = depth_category, y = mean_abundance, group = ASV)) +
    geom_line(aes(color = ASV), linewidth = 1) +
    geom_point(aes(color = ASV), size = 2) +
    scale_y_log10() +
    facet_wrap(~ Genus, scales = "fixed", ncol = 5) +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "lightgrey"),
      strip.placement = "outside",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.ticks.y = element_blank(),
      plot.margin = unit(c(1, 1, 1, 4), "lines")
    ) +
    labs(x = "Depth Category", 
         y = "Mean Abundance", 
         title = "Flavobacteriaceae ASV Abundances Across Depth Categories",
         color = "ASV") +
    scale_color_manual(values = color_palette) +
    coord_flip()
  
  return(p)
}

# Create plot for Flavobacteriaceae
p_flavobacteriaceae <- create_flavobacteriaceae_plot(bacteroidota_data)
print(p_flavobacteriaceae)
# Save the plot
ggsave("heatmaps/flavobacteriaceae_abundance_by_genus.png", p_flavobacteriaceae, width = 20, height = 16, dpi = 300)

# Print the plot
print(p_flavobacteriaceae)

#library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(grid)

# Function to create plot for Bacteroidia excluding Flavobacteriaceae
create_bacteroidota_plot <- function(data, class_name) {
  # Filter data for the specified class and exclude Flavobacteriaceae
  class_data <- data %>% 
    filter(Class == class_name, Family != "Flavobacteriaceae")
  
  # Create a color palette for genera
  n_genera <- length(unique(class_data$Genus))
  color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(n_genera)
  
  # Assign colors to genera
  color_scale <- setNames(color_palette, unique(class_data$Genus))
  
  # Ensure depth_category is a factor with correct levels
  depth_order <- rev(c("MP", "Ice-TS", "Ice-BS", "UIW", "2", "10", "25", "chl-max", "50", "100", "200",
                       "500", "1000", "1500", "2000", "3000", "20m from bottom", "5m from bottom", "bottom"))
  class_data$depth_category <- factor(class_data$depth_category, levels = depth_order)
  
  # Create the plot
  p <- ggplot(class_data, aes(x = depth_category, y = mean_abundance, group = ASV)) +
    geom_line(aes(color = Genus), linewidth = 1) +
    geom_point(aes(color = Genus), size = 2) +
    scale_y_log10() +
    facet_wrap(~ Family, scales = "fixed", ncol = 3) +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "lightgrey"),
      strip.placement = "outside",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.ticks.y = element_blank(),
      plot.margin = unit(c(1, 1, 1, 4), "lines"),
      axis.title.x = element_blank()
    ) +
    labs(x = "Depth Category", 
         title = paste(class_name, "ASV Abundances Across Depth Categories (excluding Flavobacteriaceae)"),
         color = "Genus") +
    scale_color_manual(values = color_scale) +
    coord_flip()
  
  return(p)
}

# Function to create plot for Flavobacteriaceae faceted by genus
create_flavobacteriaceae_plot <- function(data) {
  # Filter data for Flavobacteriaceae family
  family_data <- data %>% 
    filter(Family == "Flavobacteriaceae")
  
  # Define genera to be grouped together
  grouped_genera <- c("Flavobacteriaceae uc", "Formosa", "NS2b marine group", 
                      "NS4 marine group", "NS5 marine group", "Ulvibacter")
  
  # Create a new column for faceting
  family_data <- family_data %>%
    mutate(facet_group = ifelse(Genus %in% grouped_genera, "Grouped Genera", Genus))
  
  # Create a color palette for ASVs
  n_asvs <- length(unique(family_data$ASV))
  color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(n_asvs)
  
  # Ensure depth_category is a factor with correct levels
  depth_order <- rev(c("MP", "Ice-TS", "Ice-BS", "UIW", "2", "10", "25", "chl-max", "50", "100", "200",
                       "500", "1000", "1500", "2000", "3000", "20m from bottom", "5m from bottom", "bottom"))
  family_data$depth_category <- factor(family_data$depth_category, levels = depth_order)
  
  # Create the plot
  p <- ggplot(family_data, aes(x = depth_category, y = mean_abundance, group = ASV)) +
    geom_line(aes(color = Genus), linewidth = 1) +
    geom_point(aes(color = Genus), size = 2) +
    scale_y_log10() +
    facet_wrap(~ facet_group, scales = "fixed", ncol = 3) +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "lightgrey"),
      strip.placement = "outside",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.ticks.y = element_blank(),
      plot.margin = unit(c(1, 1, 1, 4), "lines")
    ) +
    labs(x = "Depth Category", 
         y = "Mean Abundance", 
         title = "Flavobacteriaceae ASV Abundances Across Depth Categories",
         color = "ASV") +
    scale_color_manual(values = color_palette) +
    coord_flip()
  
  return(p)
}

# Create plot for Bacteroidia class, excluding Flavobacteriaceae
p_bacteroidia <- create_bacteroidota_plot(bacteroidota_data, "Bacteroidia")

# Create plot for Flavobacteriaceae
p_flavobacteriaceae <- create_flavobacteriaceae_plot(bacteroidota_data)

combined_plot <- grid.arrange(
  p_bacteroidia, p_flavobacteriaceae, 
  ncol = 1, 
  heights = c(1.3, 2.2),  # This makes the Flavobacteriaceae plot twice as tall
  top = textGrob("Bacteroidota ASV Abundances Across Depth Categories",
                 gp = gpar(fontsize = 16, font = 2))
)
# Save the combined plot
ggsave("heatmaps/bacteroidota_combined_abundance.png", combined_plot, width = 24, height = 16, dpi = 300)

# Print the combined plot
print(combined_plot)
